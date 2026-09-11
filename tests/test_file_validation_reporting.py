"""Evidence-only file cards for every family, direct commands and pipelines."""

from functools import partial
from io import StringIO
from pathlib import Path
import gzip
from contextlib import nullcontext
import os
import shutil
import subprocess
import sys

import pytest
from rich.console import Console
import yaml

from postgwas.config import load_configuration
from postgwas.core.direct_execution import run_direct_with_checkpoint
from postgwas.core.errors import ConfigurationError
from postgwas.core.input_validation import (
    FileValidationRecord, InputValidationSession, current_validation_session,
    record_file_validation, validate_once,
)
from postgwas.core.ui import StageProgress, run_with_progress
from postgwas.core.validation_reporting import (
    FileValidationDisplay, group_file_validations, render_file_validation,
    run_with_file_validation,
    register_file_availability_bundle, register_file_validation_bundle,
    validation_audit_path,
    write_validation_audit,
)
from postgwas.pipeline.registry import REGISTRY


def record(path, *, role="Reference", checks=("schema",), status="passed", metrics=None, message=""):
    return FileValidationRecord(str(path) if path is not None else None, role, checks, status,
                                metrics or {}, message, ("example",))


def test_duplicate_checks_and_indexes_render_once_without_file_io(tmp_path, monkeypatch):
    config = load_configuration()
    vcf = tmp_path / "study.vcf.gz"
    records = [
        record(vcf, role="input VCF", checks=("regular file", "non-empty file")),
        record(vcf, role="Reading VCF header", checks=("Reading VCF header",)),
        record(vcf, role="Input GWAS-VCF", checks=("indexed VCF header", "required tag declarations"),
               metrics={"variants_from_index": 12519186, "declared_genome_build": "GRCh37",
                        "contigs": [*map(str, range(1, 23)), "X"],
                        "sample_count": 1, "sample": "STUDY"},
               message="Verbose implementation details belong in the audit."),
        record(str(vcf) + ".tbi", role="VCF index", checks=("index available",)),
        record(str(vcf) + ".tbi", checks=("regular file", "non-empty file")),
    ]
    monkeypatch.setattr(Path, "open", lambda *_a, **_k: pytest.fail("Renderer read a file"))
    monkeypatch.setattr(Path, "stat", lambda *_a, **_k: pytest.fail("Renderer inspected a file"))
    text = render_file_validation(records, config)
    assert text.count("study.vcf.gz") == 1
    assert "tbi: PASSED" in text
    assert "12,519,186" in text and "1–22, X" in text
    assert "Variants (from index)" in text and "GRCh37" in text
    assert "VCF sample columns" in text and "VCF sample ID" in text
    assert "Verbose implementation" not in text and "Used by" not in text
    assert "1 — passed: 1" in text
    assert len(records) == 5  # Presentation has not discarded evidence.


@pytest.mark.parametrize("role,metric,value,label", [
    ("PLINK BIM", "identifier_type", "rsid", "Variant ID type"),
    ("PLINK FAM", "samples", 200, "Samples"),
    ("PLINK BED", "variants", 1000, "Variants"),
    ("BED4 annotation reference", "rows", 1703, "Rows"),
    ("Gene coordinates", "genes", 10, "Genes"),
    ("SNP sets", "sets", 3, "Sets"),
    ("Gene sets", "gene_sets", 7, "Gene sets"),
    ("Gene covariates", "properties", 4, "Properties"),
    ("Feature matrix", "columns", 20, "Columns"),
    ("Binary kernel", "dtype", "float64", "Numeric type"),
    ("H5AD expression atlas", "cells", 300, "Cells"),
    ("Reference manifest", "declared_fields", 6, "Declared fields"),
    ("LDSC cell-type manifest", "cell_types", 12, "Cell types"),
    ("Feature names", "names", 30, "Names"),
])
def test_all_file_families_use_one_layout(tmp_path, role, metric, value, label):
    text = render_file_validation([record(tmp_path / "resource", role=role,
                                           metrics={metric: value})], load_configuration())
    assert role + " — CHECKS PASSED" in text
    assert label in text and (f"{value:,}" if isinstance(value, int) else value) in text
    # Count the actual filename field, not the broader section heading.
    assert text.count(": resource") == 1


def test_availability_is_not_promoted_and_empty_check_record_is_not_certified(tmp_path):
    text = render_file_validation([
        record(tmp_path / "unknown.bin", checks=("regular nonempty file",)),
    ], load_configuration())
    assert "availability only" in text and "CHECKS PASSED" not in text


def test_declared_availability_bundle_is_summarized_without_changing_audit(
    tmp_path, capsys,
):
    paths = tuple(tmp_path / ("%d.l2.ldscore.gz" % number) for number in range(1, 4))
    destination = tmp_path / "audit.yaml"
    with InputValidationSession() as session:
        display = FileValidationDisplay(session, load_configuration())
        for path in paths:
            record_file_validation(
                path,
                "LDSC chromosome-split reference",
                checks=("regular file", "non-empty file"),
            )
        register_file_availability_bundle(
            paths,
            "LDSC chromosome reference bundle",
            (
                ("count", "chromosomes", ["1", "2", "3"]),
                ("success", "ldsc_reference_ld_files", "3 / 3"),
                ("count", "ldsc_unique_bundle_files", 3),
            ),
        )
        display.flush()
        write_validation_audit(
            {"report_kind": "test", "status": "passed"},
            session.records,
            destination,
        )

    screen = capsys.readouterr().out
    assert "Files 1–3" in screen
    assert "LDSC chromosome reference resources — AVAILABLE — availability only" in screen
    assert "Chromosomes" in screen and "1–3" in screen
    assert "Reference LD-score files" in screen and "3 / 3" in screen
    assert "Additional files" not in screen
    assert all(path.name not in screen for path in paths)
    audit = yaml.safe_load(destination.read_text())
    assert {row["path"] for row in audit["files"]} == {str(path) for path in paths}


def test_summarized_bundle_member_failure_is_shown_individually(tmp_path, capsys):
    paths = tuple(tmp_path / ("%d.l2.ldscore.gz" % number) for number in range(1, 3))
    with InputValidationSession() as session:
        display = FileValidationDisplay(session, load_configuration())
        for path in paths:
            record_file_validation(path, "Reference", checks=("regular file",))
        register_file_availability_bundle(
            paths,
            "LDSC chromosome reference bundle",
            (("count", "ldsc_unique_bundle_files", len(paths)),),
        )
        display.flush()
        capsys.readouterr()
        record_file_validation(
            paths[0], "LDSC resource", status="failed",
            message="Reference changed after the availability check.",
        )
        display.flush()

    screen = capsys.readouterr().out
    assert "Additional findings for file 1" in screen
    assert "FAILED" in screen
    assert "Reference changed after the availability check." in screen


@pytest.mark.parametrize(
    "covered_checks,covered_metric_keys",
    (
        (
            ("companion dimensions",),
            ("rows", "columns", "dtype"),
        ),
        (
            ("companion dimensions", "all values finite and numeric"),
            ("rows", "columns"),
        ),
    ),
)
def test_validated_bundle_requires_every_hidden_check_and_metric_to_be_covered(
    tmp_path, capsys, covered_checks, covered_metric_keys,
):
    paths = tuple(tmp_path / ("chunk%d.npy" % number) for number in range(2))
    with InputValidationSession() as session:
        display = FileValidationDisplay(session, load_configuration())
        for path in paths:
            record_file_validation(
                path,
                "Feature matrix",
                checks=("companion dimensions", "all values finite and numeric"),
                metrics={"rows": 10, "columns": 2, "dtype": "float32"},
            )
        register_file_validation_bundle(
            paths,
            "PoPS feature matrix bundle",
            (("success", "pops_matrix_validation", "all covered checks passed"),),
            covered_checks=covered_checks,
            covered_metric_keys=covered_metric_keys,
        )
        display.flush()

    screen = capsys.readouterr().out
    assert "PoPS feature-matrix resources — CHECKS PASSED" not in screen
    assert all(path.name in screen for path in paths)


@pytest.mark.parametrize("status", ["failed", "blocked", "warning", "deferred"])
@pytest.mark.parametrize("suffix", [".tbi", ".csi", ".fai", ".gzi"])
def test_companion_or_cross_file_problem_overrides_pass(tmp_path, status, suffix):
    vcf = tmp_path / "study.vcf.gz"
    panel = tmp_path / "panel.bim"
    records = [record(vcf), record(str(vcf) + suffix, status=status, message="Index problem"),
               record(panel), record(None, status="failed", metrics={"paths": [str(panel)]},
                                     message="Incompatible identifiers")]
    groups, _ = group_file_validations(records)
    assert {group.path: group.status for group in groups} == {str(vcf): status, str(panel): "failed"}
    text = render_file_validation(records, load_configuration())
    assert "Index problem" in text and "Incompatible identifiers" in text
    assert "CHECKS PASSED" not in text


def test_default_limit_compacts_passes_and_never_hides_problems(tmp_path):
    config = load_configuration()
    records = [record(tmp_path / ("panel%d" % i)) for i in range(25)]
    assert config.logging.file_validation.max_screen_files == 20
    compact = render_file_validation(records, config)
    assert compact.count("— CHECKS PASSED") == 20
    assert "Additional files" in compact
    records.extend(record(tmp_path / status, status=status, message=status + " details")
                   for status in ("failed", "warning", "blocked", "deferred"))
    config.logging.file_validation.max_screen_files = 1
    text = render_file_validation(records, config)
    for status in ("failed", "warning", "blocked", "deferred"):
        assert status + " details" in text
    assert "Additional files" in text


def test_same_basename_and_orphan_index_are_not_conflated(tmp_path):
    paths = [tmp_path / "one" / "same.bim", tmp_path / "two" / "same.bim", tmp_path / "orphan.csi"]
    records = [record(path) for path in paths]
    text = render_file_validation(records, load_configuration())
    unwrapped = "".join(text.split())
    assert str(paths[0]) in unwrapped and str(paths[1]) in unwrapped
    assert len(group_file_validations(records)[0]) == 3


def test_record_only_context_preserves_direct_execution_and_cache_boundaries(tmp_path):
    source = tmp_path / "input.txt"
    source.write_text("first")
    calls = []

    def inspect():
        calls.append(source.read_text())
        record_file_validation(source, "Example", checks=("content",))

    with InputValidationSession() as pipeline:
        with InputValidationSession(reuse_checks=False) as direct:
            assert current_validation_session() is None
            validate_once((source,), "test", inspect)
            source.write_text("changed")
            validate_once((source,), "test", inspect)
        assert current_validation_session() is pipeline
        assert not pipeline.records
    assert calls == ["first", "changed"]
    assert len(direct.records) == 1


def test_direct_stages_display_once_and_keep_complete_audit(tmp_path, capsys):
    config = load_configuration()
    source = tmp_path / "input.txt"
    source.write_text("fixture")
    calls = []

    def operation():
        progress = StageProgress("Example", enabled=True)
        for number in (1, 2):
            with progress.step(number, 2, "Existing stage"):
                calls.append(number)
                record_file_validation(source, "Example", checks=("schema",), metrics={"rows": 1},
                                       message="Detailed successful explanation.")
            if number == 1:
                assert "File validation" in capsys.readouterr().out

    run_with_file_validation(operation, command="example", configuration=config, output_directory=tmp_path / "out")
    report = yaml.safe_load((tmp_path / "out/run_metadata/example_input_validation.yaml").read_text())
    assert report["status"] == "passed" and calls == [1, 2]
    assert len(report["files"]) == 1
    assert report["files"][0]["message"] == "Detailed successful explanation."
    assert "input.txt" not in capsys.readouterr().out
    assert current_validation_session() is None


@pytest.mark.parametrize("exit_mode", ["exception", "return", "system_exit"])
def test_direct_failures_keep_original_failure_and_file_evidence(tmp_path, exit_mode):
    config = load_configuration()

    def operation():
        def check():
            raise ValueError("Malformed position")
        try:
            validate_once((tmp_path / "broken.tsv",), "test", check)
        except ValueError:
            if exit_mode == "return":
                return 2
            if exit_mode == "system_exit":
                raise SystemExit(2)
            raise

    run = partial(run_with_file_validation, operation, command="example", configuration=config,
                  output_directory=tmp_path / "out")
    if exit_mode == "return":
        assert run() == 2
    else:
        with pytest.raises(ValueError if exit_mode == "exception" else SystemExit):
            run()
    report = yaml.safe_load((tmp_path / "out/run_metadata/example_input_validation.yaml").read_text())
    assert report["status"] == "failed"
    assert any(row["message"] == "Malformed position" for row in report["files"])
    assert current_validation_session() is None


def test_direct_checkpoint_reuse_does_not_rewrite_audit_or_repeat_checks(tmp_path):
    output = tmp_path / "out"
    config = load_configuration(cli_overrides={"run.output_directory": str(output)})
    source = tmp_path / "source.txt"
    source.write_text("fixture")
    calls = []

    def operation():
        calls.append(1)
        record_file_validation(source, "Example", checks=("schema",))
        output.mkdir(exist_ok=True)
        (output / "result.txt").write_text("validated result")
        return 0

    recorded = partial(run_with_file_validation, operation, command="manhattan", configuration=config,
                       output_directory=output)
    run = partial(run_direct_with_checkpoint, recorded, module_name="manhattan", public_command="manhattan",
                  arguments=(), configuration=config, output_directory=output,
                  console=Console(file=StringIO()))
    assert run() == 0
    audit = output / "run_metadata/manhattan_input_validation.yaml"
    before = audit.stat().st_mtime_ns, audit.read_bytes()
    assert run() == 0 and calls == [1]
    assert (audit.stat().st_mtime_ns, audit.read_bytes()) == before


def test_unknown_file_and_cross_file_input_are_protected_from_audit_overwrite(tmp_path):
    path = tmp_path / "resource.tsv"
    path.write_text("SNP\tP\nrs1\t0.5\n")
    document = {"report_kind": "example"}
    with pytest.raises(ValueError, match="unrecognised"):
        write_validation_audit(document, [], path)
    with pytest.raises(ValueError, match="validated input"):
        write_validation_audit(document, [record(None, metrics={"paths": [str(path)]})], path)
    assert path.read_text() == "SNP\tP\nrs1\t0.5\n"


def test_symlinked_audit_parent_is_refused_even_when_it_stays_inside_output(tmp_path):
    (tmp_path / "real").mkdir()
    (tmp_path / "link").symlink_to(tmp_path / "real", target_is_directory=True)
    with pytest.raises(ValueError, match="symlink"):
        validation_audit_path(tmp_path, "link/audit.yaml")
    assert not (tmp_path / "real/audit.yaml").exists()


def test_audit_write_failure_stops_before_successful_top_level_progress(tmp_path, monkeypatch):
    config = load_configuration()
    stream = StringIO()

    def operation():
        record_file_validation(tmp_path / "input", "Example", checks=("schema",))

    def fail_write(*_args):
        raise OSError("audit disk full")

    monkeypatch.setattr("postgwas.core.validation_reporting.write_validation_audit", fail_write)
    with pytest.raises(OSError, match="audit disk full"):
        run_with_progress(partial(run_with_file_validation, operation, command="example", configuration=config,
                                  output_directory=tmp_path), label="Run", title="Analysis", enabled=True,
                          console=Console(file=stream))
    assert "Completed 1/1" not in stream.getvalue()
    assert "Failed 1/1" in stream.getvalue()


@pytest.mark.parametrize("command", [name for name, spec in REGISTRY.commands().items()
                                     if spec.direct_checkpoint != "not_applicable"])
def test_every_scientific_direct_dispatcher_includes_audit_publication_in_progress(tmp_path, monkeypatch, command):
    from types import SimpleNamespace
    from postgwas import __main__ as dispatcher

    stream = StringIO()
    console = Console(file=stream, color_system=None)
    config = load_configuration()

    def analysis():
        progress = StageProgress("Native analysis", enabled=True, console=console)
        with progress.step(1, 2, "Check resource"):
            record_file_validation(tmp_path / "input", "Example", checks=("schema",))
        with progress.step(2, 2, "Finish analysis"):
            pass

    def fail_write(*_args):
        raise OSError("audit disk full")

    monkeypatch.setattr(dispatcher, "console", console)
    monkeypatch.setattr(dispatcher, "resolve_reference", lambda _reference: analysis)
    monkeypatch.setattr(dispatcher, "resolve_screen_settings", lambda *_args: SimpleNamespace(
        configuration=config, output_directory=tmp_path,
    ))
    monkeypatch.setattr(dispatcher, "record_screen", lambda _settings: nullcontext())
    monkeypatch.setattr("postgwas.core.direct_execution.run_direct_with_checkpoint",
                        lambda operation, **_kwargs: operation())
    monkeypatch.setattr("postgwas.core.validation_reporting.write_validation_audit", fail_write)
    monkeypatch.setattr(sys, "argv", ["postgwas", command, "--output-directory", str(tmp_path)])
    assert dispatcher.main() == 1
    assert "Completed 2/2" in stream.getvalue()  # The native analysis really finished.
    assert "Completed 1/1" not in stream.getvalue()  # The complete command did not.
    assert "Failed 1/1" in stream.getvalue() and "audit disk full" in stream.getvalue()


@pytest.mark.parametrize("value", ["not-an-integer", "-1"])
def test_invalid_index_counts_cannot_leave_a_passed_file_card(tmp_path, monkeypatch, value):
    from postgwas.core import vcf

    source = tmp_path / "input.vcf.gz"
    monkeypatch.setattr(vcf, "_read_vcf_validation_command", lambda *_args, **_kwargs: value)
    with InputValidationSession(reuse_checks=False) as session:
        with pytest.raises(vcf.VcfQueryError):
            vcf.count_indexed_vcf_records(source, "unused")
    assert group_file_validations(session.records)[0][0].status == "failed"


def test_direct_native_checkpoint_does_not_own_an_unfinished_audit(tmp_path):
    output = tmp_path / "out"
    config = load_configuration(cli_overrides={"run.output_directory": str(output)})
    calls = []

    def analysis():
        calls.append(1)
        record_file_validation(tmp_path / "input", "Example", checks=("schema",))
        output.mkdir(exist_ok=True)
        (output / "result.txt").write_text("validated result")

    native = partial(run_direct_with_checkpoint, analysis, module_name="manhattan", public_command="manhattan",
                     arguments=(), configuration=config, output_directory=output,
                     console=Console(file=StringIO()))
    run = partial(run_with_file_validation, native, command="manhattan", configuration=config,
                  output_directory=output)
    run()
    audit = output / "run_metadata/manhattan_input_validation.yaml"
    before = audit.stat().st_mtime_ns, audit.read_bytes()
    checkpoint = yaml.safe_load((output / "run_metadata/checkpoints/manhattan_direct.yaml").read_text())
    assert "run_metadata/manhattan_input_validation.yaml" not in checkpoint["outputs"]
    run()
    assert calls == [1]
    assert (audit.stat().st_mtime_ns, audit.read_bytes()) == before


@pytest.mark.parametrize("key,value", [
    ("max_screen_files", 0), ("max_checks_per_file", 0), ("max_list_items", 0),
    ("direct_report_file", "../audit.yaml"), ("direct_report_file", "{unknown}.yaml"),
])
def test_display_configuration_is_schema_validated(key, value):
    with pytest.raises(ConfigurationError):
        load_configuration(cli_overrides={"logging.file_validation." + key: value})


@pytest.mark.parametrize("pipeline", [False, True])
@pytest.mark.parametrize("valid_bed", [False, True])
def test_public_ld_annotation_file_cards_and_unchanged_scientific_output(tmp_path, pipeline, valid_bed):
    bcftools = shutil.which("bcftools")
    if bcftools is None:
        pytest.skip("Native bcftools is required")
    root = Path(__file__).resolve().parents[1]
    vcf = tmp_path / "study.vcf.gz"
    subprocess.run([bcftools, "view", "-Oz", "-o", str(vcf),
                    str(root / "tests/fixtures/formatting/harmonised.vcf")], check=True, capture_output=True)
    subprocess.run([bcftools, "index", str(vcf)], check=True, capture_output=True)
    beds = tmp_path / "beds"
    beds.mkdir()
    with gzip.open(beds / "GRCh37_EUR_ldetect.bed.gz", "wt") as handle:
        handle.write("1\t0\t1000\tblock1\n2\t0\t1000\tblock2\n" if valid_bed else "1\tbad\t1000\tblock1\n")
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump({"modules": {"ld_annotation": {"populations": ["EUR"]}}}))
    output = tmp_path / "out"
    command = ["pipeline", "--modules", "annot_ldblock"] if pipeline else ["annot_ldblock"]
    completed = subprocess.run([
        sys.executable, "-m", "postgwas", *command, "--vcf", str(vcf), "--ld-region-dir", str(beds),
        "--dataset-id", "STUDY", "--output-directory", str(output), "--run-config", str(config), "--hide-screen",
    ], cwd=root, env={**os.environ, "PYTHONPATH": str(root / "src"), "PYTHONDONTWRITEBYTECODE": "1"},
        text=True, capture_output=True)
    screen = (output / "run_metadata/screen.log").read_text()
    assert completed.stdout == "" and completed.stderr == ""
    assert (completed.returncode == 0) is valid_bed, screen
    assert "File validation" in screen
    if valid_bed:
        assert screen.count(": study.vcf.gz") == 1
        input_section = screen.split(": study.vcf.gz", 1)[1].split("Checks", 1)[0]
        assert input_section.count(": GRCh37") == 1
        assert input_section.count("Variants (from index)") == 1
    assert "BED4 annotation reference — " + ("CHECKS PASSED" if valid_bed else "FAILED") in screen
    assert "LD block BED file missing or empty : PASSED" not in screen
    report_name = "input_validation.yaml" if pipeline else "annot_ldblock_input_validation.yaml"
    audit = yaml.safe_load((output / "run_metadata" / report_name).read_text())
    assert audit["status"] == ("passed" if valid_bed else "failed")
    assert any(row["metrics"].get("variants_from_index") == 4 for row in audit["files"])
    delivered = output / ("01_annot_ldblock" if pipeline else "") / "STUDY_ldblock.vcf.gz"
    if valid_bed:
        query = subprocess.run([bcftools, "query", "-f", "%CHROM\t%POS\t%INFO/EUR_LDblock\n", str(delivered)],
                               text=True, capture_output=True, check=True)
        assert query.stdout.splitlines() == ["1\t100\tblock1", "1\t200\tblock1", "2\t300\tblock2", "2\t400\tblock2"]
        scientific_query = "%CHROM\t%POS\t%REF\t%ALT[\t%ES\t%SE\t%LP\t%AF]\n"
        values = [subprocess.run([bcftools, "query", "-f", scientific_query, str(path)],
                                 text=True, capture_output=True, check=True).stdout
                  for path in (vcf, delivered)]
        assert values[0] == values[1]
    else:
        assert not delivered.exists()
