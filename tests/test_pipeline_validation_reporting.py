"""Pipeline-only common validation, audit publication and fail-before-run tests."""

from argparse import Namespace
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from postgwas.config import load_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.core.input_validation import (
    InputValidationSession,
    current_validation_session,
    record_file_validation,
    validate_once,
)
from postgwas.core.preflight import PipelinePreflightEvidence
from postgwas.pipeline import cli
from postgwas.pipeline.validation_reporting import (
    save_pipeline_validation_report,
    validation_report_path,
)
from preflight_support import pipeline_input_vcf_evidence


def configuration_for(tmp_path):
    return load_configuration(cli_overrides={
        "run.output_directory": str(tmp_path / "out"),
        "run.dataset_id": "STUDY",
    })


def test_independent_preflight_errors_are_aggregated_and_inferences_not_published(monkeypatch):
    initial = pipeline_input_vcf_evidence()
    visited = []

    def validator(name):
        def check(args, *, preflight_evidence):
            visited.append(name)
            if name != "third":
                args.partial_inference = name
                raise ValueError("Invalid resource for %s" % name)
            assert not hasattr(args, "partial_inference")
            return PipelinePreflightEvidence(name, initial["input_vcf"])
        return check

    monkeypatch.setattr(cli.REGISTRY, "get", lambda name: SimpleNamespace(preflight=name))
    monkeypatch.setattr(cli, "resolve_reference", validator)
    monkeypatch.setattr("postgwas.pipeline.resource_validation.record_pipeline_resource_validation", lambda *_: None)
    with InputValidationSession() as session:
        with pytest.raises(cli.PipelinePreflightError) as caught:
            cli._run_pipeline_preflights(Namespace(), ("first", "second", "third"), initial_evidence=initial)
    assert visited == ["first", "second", "third"]
    assert "Invalid resource for first" in str(caught.value)
    assert "Invalid resource for second" in str(caught.value)
    assert "third" in caught.value.evidence
    assert [record.consumers for record in session.records] == [("first",), ("second",)]


def test_entry_failure_saves_failed_audit_and_never_starts_analysis(tmp_path, monkeypatch):
    config = configuration_for(tmp_path)
    args = Namespace(vcf=str(tmp_path / "missing.vcf.gz"), output_directory=config.run.output_directory)
    plan = SimpleNamespace(active_modules=("formatter", "magma"))
    monkeypatch.setattr(cli, "_validate_pipeline_entry_vcf", lambda *_: (_ for _ in ()).throw(ValueError("Malformed GWAS-VCF")))
    monkeypatch.setattr(cli, "execute_pipeline", lambda *_args, **_kwargs: pytest.fail("Analysis started"))
    with pytest.raises(cli.PipelinePreflightError, match="Malformed GWAS-VCF"):
        cli._prepare_and_execute_pipeline(args, plan, config)
    report = yaml.safe_load(validation_report_path(args, config).read_text())
    assert report["status"] == "failed"
    assert report["phase"] == "pipeline_startup"
    assert {consumer for item in report["files"] if item["status"] == "blocked"
            for consumer in item["consumers"]} == {"formatter", "magma"}
    assert current_validation_session() is None


def test_successful_preflight_session_is_reused_during_execution(tmp_path, monkeypatch):
    config = configuration_for(tmp_path)
    args = Namespace(vcf=str(tmp_path / "study.vcf.gz"), output_directory=config.run.output_directory)
    plan = SimpleNamespace(active_modules=("formatter",))
    resource = tmp_path / "resource.txt"
    resource.write_text("resource\n")
    initial = pipeline_input_vcf_evidence()
    calls = []

    def check_resource():
        calls.append(1)
        record_file_validation(resource, "Example resource", checks=("example content check",), metrics={"rows": 1})
        return "validated"

    def preflight(_args, _modules, *, initial_evidence):
        with current_validation_session().scope("formatter"):
            assert validate_once((resource,), {"validator": "test-resource"}, check_resource) == "validated"
        return {**initial_evidence, "formatter": PipelinePreflightEvidence(
            "formatter", initial["input_vcf"], deferred_checks=("Validate generated table.",),
        )}

    def execute(*_args, **_kwargs):
        assert current_validation_session() is not None
        assert validate_once((resource,), {"validator": "test-resource"}, check_resource) == "validated"
        _kwargs["finalize_validation"]()
        return "executed"

    monkeypatch.setattr(cli, "_validate_pipeline_entry_vcf", lambda *_: initial["input_vcf"])
    monkeypatch.setattr(cli, "_run_pipeline_preflights", preflight)
    monkeypatch.setattr(cli, "resolve_compute_args", lambda *_: None)
    monkeypatch.setattr(cli, "execute_pipeline", execute)
    assert cli._prepare_and_execute_pipeline(args, plan, config) == "executed"
    report = yaml.safe_load(validation_report_path(args, config).read_text())
    assert report["status"] == "passed"
    assert report["phase"] == "pipeline_execution"
    assert report["stages"][0]["deferred_checks"] == ["Validate generated table."]
    assert calls == [1]
    assert current_validation_session() is None


def test_report_contains_every_file_when_screen_is_limited(tmp_path, monkeypatch):
    config = configuration_for(tmp_path)
    config.logging.file_validation.max_screen_files = 1
    output = tmp_path / "audit.yaml"
    captured = []
    monkeypatch.setattr("postgwas.core.validation_reporting.print_screen_block", captured.append)
    with InputValidationSession() as session:
        with session.scope("magma"):
            record_file_validation(tmp_path / "one.bim", "BIM", checks=("ID convention",), metrics={"rsids": 2})
            record_file_validation(tmp_path / "two.bim", "BIM", status="failed", message="Mixed IDs")
        save_pipeline_validation_report(session, {}, ("magma",), output, error=ValueError("Mixed IDs"))
        from postgwas.core.validation_reporting import FileValidationDisplay
        FileValidationDisplay(session, config).flush(report_path=output)
    report = yaml.safe_load(output.read_text())
    assert len(report["files"]) == 2
    assert "two.bim" in captured[0] and "Mixed IDs" in captured[0]
    assert "Additional files" in captured[0]
    assert "one.bim" not in captured[0]


def test_audit_cannot_overwrite_a_validated_resource(tmp_path):
    resource = tmp_path / "resource.txt"
    resource.write_text("original\n")
    with InputValidationSession() as session:
        record_file_validation(resource, "Resource", checks=("availability",))
        with pytest.raises(ValueError, match="overwrite a validated input"):
            save_pipeline_validation_report(session, {}, (), resource)
    assert resource.read_text() == "original\n"


def test_report_write_failure_stops_before_analysis(tmp_path, monkeypatch):
    config = configuration_for(tmp_path)
    args = Namespace(vcf=str(tmp_path / "study.vcf.gz"), output_directory=config.run.output_directory)
    plan = SimpleNamespace(active_modules=())
    initial = pipeline_input_vcf_evidence()
    monkeypatch.setattr(cli, "_validate_pipeline_entry_vcf", lambda *_: initial["input_vcf"])
    monkeypatch.setattr(cli, "_run_pipeline_preflights", lambda *_args, initial_evidence: initial_evidence)
    monkeypatch.setattr(cli, "save_pipeline_validation_report", lambda *_args, **_kwargs: (_ for _ in ()).throw(OSError("disk full")))
    monkeypatch.setattr(cli, "execute_pipeline", lambda *_args, **_kwargs: pytest.fail("Analysis started"))
    with pytest.raises(cli.PipelinePreflightError, match="disk full"):
        cli._prepare_and_execute_pipeline(args, plan, config)


def test_audit_cannot_overwrite_an_unvisited_resource(tmp_path):
    resource = tmp_path / "unvisited.tsv"
    resource.write_text("SNP\tP\nrs1\t0.5\n")
    with InputValidationSession() as session:
        with pytest.raises(ValueError, match="unrecognised validation audit"):
            save_pipeline_validation_report(session, {}, (), resource)
    assert resource.read_text() == "SNP\tP\nrs1\t0.5\n"


def test_execution_failure_preserves_runtime_validation_evidence(tmp_path, monkeypatch):
    config = configuration_for(tmp_path)
    args = Namespace(vcf=str(tmp_path / "study.vcf.gz"), output_directory=config.run.output_directory)
    plan = SimpleNamespace(active_modules=("formatter",))
    initial = pipeline_input_vcf_evidence()
    monkeypatch.setattr(cli, "_validate_pipeline_entry_vcf", lambda *_: initial["input_vcf"])
    monkeypatch.setattr(cli, "_run_pipeline_preflights", lambda *_args, initial_evidence: initial_evidence)
    monkeypatch.setattr(cli, "resolve_compute_args", lambda *_: None)

    def execute(*_args, **_kwargs):
        with current_validation_session().scope("formatter"):
            record_file_validation(tmp_path / "generated.tsv", "Generated input", status="failed", message="Invalid allele")
        raise ValueError("Invalid generated input")

    monkeypatch.setattr(cli, "execute_pipeline", execute)
    with pytest.raises(ValueError, match="Invalid generated input"):
        cli._prepare_and_execute_pipeline(args, plan, config)
    report = yaml.safe_load(validation_report_path(args, config).read_text())
    assert report["phase"] == "pipeline_execution"
    assert report["status"] == "failed"
    assert any(item["role"] == "Generated input" and item["status"] == "failed" for item in report["files"])


@pytest.mark.parametrize("resume_checkpoint", (False, True))
def test_final_audit_failure_cannot_report_pipeline_completion(tmp_path, monkeypatch, resume_checkpoint):
    from io import StringIO

    from rich.console import Console

    from postgwas.core.errors import ModuleExecutionError
    from postgwas.pipeline import executor
    from postgwas.pipeline.planner import PipelinePlan

    config = configuration_for(tmp_path)
    args = Namespace(output_directory=config.run.output_directory, resume=True, overwrite=False)
    plan = PipelinePlan(("formatter",), ("formatter",), ("formatter",))
    stream = StringIO()
    monkeypatch.setattr(executor, "console", Console(file=stream, color_system=None, width=120))
    monkeypatch.setattr(executor.REGISTRY, "require_pipeline_enabled", lambda _: SimpleNamespace(
        runner="example", description="Format example", pipeline_output_name=None,
    ))

    def runner(_args, context):
        output = tmp_path / "out" / "01_formatter" / "result.tsv"
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("SNP\nrs1\n")
        context["formatter"] = {"output": str(output)}
        return context["formatter"]

    monkeypatch.setattr(executor, "resolve_reference", lambda _: runner)

    if resume_checkpoint:
        executor.execute_pipeline(args, plan, config)
        stream.seek(0)
        stream.truncate()

    def fail_audit():
        raise OSError("audit disk full")

    with pytest.raises(ModuleExecutionError, match="audit disk full"):
        executor.execute_pipeline(args, plan, config, finalize_validation=fail_audit)
    assert "Completed 1/1" not in stream.getvalue()
    assert "All tasks completed successfully" not in stream.getvalue()
    assert "Pipeline Failed at step: formatter" in stream.getvalue()


@pytest.mark.parametrize("path", ("../audit.yaml", "/tmp/audit.yaml", "."))
def test_validation_report_configuration_refuses_unsafe_paths(path):
    with pytest.raises(ConfigurationError, match="safe relative file path"):
        load_configuration(cli_overrides={"pipeline.validation.report_file": path})


def test_report_filename_is_yaml_configurable(tmp_path):
    config = configuration_for(tmp_path)
    config.pipeline.validation.report_file = Path("audit/custom.yaml")
    args = Namespace(output_directory=config.run.output_directory)
    assert validation_report_path(args, config) == config.run.output_directory / "audit/custom.yaml"


@pytest.mark.parametrize("cli_output", (None, "explicit"))
def test_yaml_common_outputs_resolve_before_preflight(tmp_path, monkeypatch, cli_output):
    config = configuration_for(tmp_path)
    explicit = tmp_path / "explicit" if cli_output else None
    args = Namespace(vcf=str(tmp_path / "missing.vcf.gz"),
                     output_directory=explicit, dataset_id=None)
    plan = SimpleNamespace(active_modules=("formatter",))

    def fail_entry(observed, configuration):
        assert observed.output_directory == (explicit or config.run.output_directory)
        assert observed.dataset_id == config.run.dataset_id
        raise ValueError("Expected invalid input")

    monkeypatch.setattr(cli, "_validate_pipeline_entry_vcf", fail_entry)
    with pytest.raises(cli.PipelinePreflightError, match="Expected invalid input"):
        cli._prepare_and_execute_pipeline(args, plan, config)
    report = yaml.safe_load(validation_report_path(args, config).read_text())
    assert report["status"] == "failed"


@pytest.mark.parametrize("indexed", (True, False))
def test_public_pipeline_hidden_screen_validation_audit(tmp_path, indexed):
    import os
    import shutil
    import subprocess
    import sys

    bcftools = shutil.which("bcftools")
    if bcftools is None:
        pytest.skip("Native bcftools is required for indexed VCF integration")
    repository = Path(__file__).resolve().parents[1]
    fixture = repository / "tests" / "fixtures" / "formatting" / "harmonised.vcf"
    vcf = tmp_path / "study.vcf.gz"
    subprocess.run([bcftools, "view", "-Oz", "-o", str(vcf), str(fixture)], check=True, capture_output=True)
    if indexed:
        subprocess.run([bcftools, "index", str(vcf)], check=True, capture_output=True)
    config = configuration_for(tmp_path)
    environment = {**os.environ, "PYTHONPATH": str(repository / "src"), "PYTHONDONTWRITEBYTECODE": "1"}
    completed = subprocess.run([
        sys.executable, "-m", "postgwas", "pipeline", "--modules", "formatter",
        "--format", "magma", "--vcf", str(vcf), "--dataset-id", "STUDY",
        "--output-directory", str(config.run.output_directory), "--hide-screen",
    ], cwd=repository, env=environment, capture_output=True, text=True, check=False)
    assert completed.stdout == "" and completed.stderr == ""
    transcript = (config.run.output_directory / config.logging.screen_log_file).read_text()
    assert completed.returncode == (0 if indexed else 2), transcript
    report = yaml.safe_load((config.run.output_directory / config.pipeline.validation.report_file).read_text())
    assert report["status"] == ("passed" if indexed else "failed")
    assert "File validation" in transcript
    assert "Pipeline input validation" in transcript
    if indexed:
        assert report["phase"] == "pipeline_execution"
        assert any(item["metrics"].get("variants_from_index") == 4 for item in report["files"])
        assert "All tasks completed successfully" in transcript
    else:
        assert report["phase"] == "pipeline_startup"
        assert "Starting Execution Chain" not in transcript
        assert not (config.run.output_directory / "01_formatter").exists()


@pytest.mark.parametrize("formats_override", (None, "gcta_gene"))
def test_public_pipeline_uses_yaml_output_and_formatter_formats(tmp_path, formats_override):
    import os
    import shutil
    import subprocess
    import sys

    bcftools = shutil.which("bcftools")
    if bcftools is None:
        pytest.skip("Native bcftools required for pipeline integration")
    repository = Path(__file__).resolve().parents[1]
    fixture = repository / "tests/fixtures/formatting/harmonised.vcf"
    vcf = tmp_path / "study.vcf.gz"
    subprocess.run([bcftools, "view", "-Oz", "-o", str(vcf), str(fixture)], check=True)
    subprocess.run([bcftools, "index", str(vcf)], check=True)
    output = tmp_path / "yaml_outputs"
    run_config = tmp_path / "run.yaml"
    run_config.write_text(yaml.safe_dump({
        "run": {"dataset_id": "STUDY", "output_directory": str(output)},
        "modules": {"formatting": {"formats": ["magma"]}},
    }))
    command = [sys.executable, "-B", "-m", "postgwas", "pipeline",
               "--modules", "formatter", "--vcf", str(vcf),
               "--run-config", str(run_config), "--hide-screen"]
    if formats_override:
        command.extend(["--format", formats_override])
    completed = subprocess.run(command, cwd=repository, capture_output=True, text=True,
                               env={**os.environ, "PYTHONPATH": str(repository / "src"),
                                    "PYTHONDONTWRITEBYTECODE": "1"})
    transcript = (output / "run_metadata/screen.log").read_text()
    assert completed.returncode == 0, transcript + completed.stdout + completed.stderr
    report = yaml.safe_load((output / "run_metadata/input_validation.yaml").read_text())
    assert report["status"] == "passed"
    assert "All tasks completed successfully" in transcript
    if formats_override:
        assert list(output.rglob("STUDY_gcta.ma"))
        assert not list(output.rglob("STUDY_magma_p_values.tsv"))
    else:
        assert list(output.rglob("STUDY_magma_p_values.tsv"))
