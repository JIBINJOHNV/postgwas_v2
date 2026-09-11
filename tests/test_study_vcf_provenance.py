"""Study-origin validation is shared, header-only and not a reference contract."""

from argparse import Namespace
import json
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace

import pytest
import yaml

from postgwas.config import load_configuration
from postgwas.core import vcf
from postgwas.core.errors import FormattingError
from postgwas.core.input_validation import InputValidationSession
from postgwas.core.preflight import PreflightLogBuffer
from postgwas.pipeline import cli
from postgwas.pipeline.validation_reporting import validation_report_path


@pytest.fixture
def contract():
    return load_configuration().modules.formatting.input_contract


def _header(contract, values=None):
    values = {"version": "0.9.0", "dataset_id": "SOURCE_STUDY",
              "status": "created and structurally validated", **(values or {})}
    return "##fileformat=VCFv4.2\n##genome_build=GRCh37\n" + "".join(
        "##%s=%s\n" % (name, json.dumps(values[role]))
        for role, name in contract.provenance_headers.model_dump().items()
        if values[role] is not None
    )


def test_declared_provenance_accepts_different_software_version_and_source_id(contract, tmp_path):
    header = _header(contract)
    source = tmp_path / "study.vcf"
    source.write_text(header)
    logger = PreflightLogBuffer()
    with InputValidationSession() as session:
        evidence = vcf.validate_postgwas_vcf_provenance(
            header, contract.provenance_headers.model_dump(), vcf_path=source, logger=logger,
        )
    assert evidence == {"postgwas_version": "0.9.0", "postgwas_dataset_id": "SOURCE_STUDY",
                        "postgwas_status": "created and structurally validated"}
    assert source.read_text() == header
    assert session.records[0].metrics == evidence
    assert session.records[0].status == "passed"
    assert logger.events


@pytest.mark.parametrize("role", ["version", "dataset_id", "status"])
@pytest.mark.parametrize("value,label", [(None, "missing"), ("", "empty"), ("   ", "empty"),
                                        ("<ID=wrong>", "invalid text"), ("bad\nvalue", "invalid text")])
def test_invalid_origin_is_actionable_and_logged(contract, tmp_path, role, value, label):
    source = tmp_path / "foreign.vcf"
    logger = PreflightLogBuffer()
    with InputValidationSession() as session:
        with pytest.raises(FormattingError) as caught:
            vcf.validate_postgwas_vcf_provenance(
                _header(contract, {role: value}), contract.provenance_headers.model_dump(),
                vcf_path=source, logger=logger, error_type=FormattingError,
            )
    message = str(caught.value)
    assert str(source) in message and label in message
    assert getattr(contract.provenance_headers, role) in message
    assert "Re-run PostGWAS harmonisation" in message
    assert session.records[0].status == "failed"
    assert logger.events


@pytest.mark.parametrize("suffix,expected", [
    ('##{version}="second"\n', "more than once"),
    ('##{version}="unterminated\n', "invalid quoted value"),
])
def test_duplicate_and_malformed_origin_cannot_pass(contract, suffix, expected):
    names = contract.provenance_headers.model_dump()
    header = _header(contract, {"version": None}) + suffix.format(**names)
    if "second" in suffix:
        header = _header(contract) + suffix.format(**names)
    with pytest.raises(vcf.VcfQueryError, match=expected):
        vcf.validate_postgwas_vcf_provenance(header, names)


def test_all_missing_headers_reported_together_and_custom_names_honoured(contract):
    names = {role: "custom_" + name for role, name in contract.provenance_headers.model_dump().items()}
    with pytest.raises(vcf.VcfQueryError) as caught:
        vcf.validate_postgwas_vcf_provenance(_header(contract), names)
    assert all(name in str(caught.value) for name in names.values())
    header = "\n".join("##%s=custom-value" % name for name in names.values())
    assert set(vcf.validate_postgwas_vcf_provenance(header, names).values()) == {"custom-value"}


def _mock_indexed_input(tmp_path, monkeypatch, header):
    source = tmp_path / "study.vcf.gz"
    source.write_bytes(b"unchanged fixture")
    Path(str(source) + ".tbi").write_bytes(b"index")
    calls = []

    def run(command, *args, **kwargs):
        operation = tuple(command[1:3])
        calls.append(operation)
        return {("view", "--header-only"): header, ("index", "-n"): "2\n",
                ("query", "-l"): "SOURCE_STUDY\n"}[operation]

    monkeypatch.setattr(vcf, "run_checked_command", run)
    return source, calls


def test_origin_checked_before_count_or_extraction(contract, tmp_path, monkeypatch):
    source, calls = _mock_indexed_input(tmp_path, monkeypatch, _header(contract, {"version": None}))
    with pytest.raises(vcf.VcfQueryError, match="missing PostGWAS"):
        vcf.validate_indexed_vcf(
            source, "OUTPUT_LABEL", "bcftools",
            genome_build_header="##" + contract.genome_build_metadata,
            supported_genome_builds=contract.supported_genome_builds,
            provenance_headers=contract.provenance_headers.model_dump(),
        )
    assert calls == [("view", "--header-only")]


def test_reference_contract_and_cached_header_do_not_bypass_new_study_gate(contract, tmp_path, monkeypatch):
    source, calls = _mock_indexed_input(tmp_path, monkeypatch, "##genome_build=GRCh37\n")
    kwargs = dict(genome_build_header="##" + contract.genome_build_metadata,
                  supported_genome_builds=contract.supported_genome_builds)
    # Generic reference/internal validation still works without PostGWAS markers.
    cached = vcf.validate_indexed_vcf(source, "SOURCE_STUDY", "bcftools", **kwargs)
    assert cached.variant_count == 2
    with pytest.raises(vcf.VcfQueryError, match="missing PostGWAS"):
        vcf.validate_indexed_vcf(source, "SOURCE_STUDY", "bcftools", cached=cached,
                                provenance_headers=contract.provenance_headers.model_dump(), **kwargs)
    assert len(calls) == 3  # No second header read when the same file is reused.


def test_module_fields_remain_independent_and_header_is_shared(contract, tmp_path, monkeypatch):
    source, calls = _mock_indexed_input(tmp_path, monkeypatch, _header(contract))
    kwargs = dict(genome_build_header="##" + contract.genome_build_metadata,
                  supported_genome_builds=contract.supported_genome_builds,
                  provenance_headers=contract.provenance_headers.model_dump())
    with InputValidationSession() as session:
        for consumer in ("manhattan", "ld_clump"):
            with session.scope(consumer):
                assert vcf.validate_indexed_vcf(source, "SOURCE_STUDY", "bcftools", **kwargs).variant_count == 2
        with pytest.raises(vcf.VcfQueryError, match="FORMAT/LP"):
            vcf.validate_indexed_vcf(source, "SOURCE_STUDY", "bcftools", required_fields=["FORMAT/LP"], **kwargs)
    assert calls == [("view", "--header-only"), ("index", "-n"), ("query", "-l")]


@pytest.mark.parametrize("target", ["formatter", "magma", "pops", "flames"])
def test_simple_pipeline_missing_origin_fails_before_any_module(target, tmp_path, monkeypatch):
    config = load_configuration(cli_overrides={"run.output_directory": str(tmp_path / "out"),
                                              "run.dataset_id": "OUTPUT_LABEL"})
    header = (Path(__file__).parent / "fixtures/formatting/harmonised.vcf").read_text()
    name = config.modules.formatting.input_contract.provenance_headers.version
    header = "\n".join(line for line in header.splitlines() if not line.startswith("##" + name + "="))
    source, calls = _mock_indexed_input(tmp_path, monkeypatch, header)
    monkeypatch.setattr(cli, "resolve_executable", lambda *args, **kwargs: sys.executable)
    monkeypatch.setattr(cli, "_run_pipeline_preflights", lambda *args, **kwargs: pytest.fail("Module preflight ran"))
    monkeypatch.setattr(cli, "execute_pipeline", lambda *args, **kwargs: pytest.fail("Analysis ran"))
    args = Namespace(vcf=str(source), output_directory=config.run.output_directory)
    plan = SimpleNamespace(active_modules=(target,))
    with pytest.raises(cli.PipelinePreflightError, match="missing PostGWAS") as caught:
        cli._prepare_and_execute_pipeline(args, plan, config)
    assert str(source) in str(caught.value)
    report = yaml.safe_load(validation_report_path(args, config).read_text())
    assert report["status"] == "failed"
    assert any(row["status"] == "failed" and "PostGWAS" in row.get("message", "") for row in report["files"])
    assert source.read_bytes() == b"unchanged fixture"


@pytest.mark.parametrize("missing_origin", [False, True])
def test_native_bcftools_origin_gate_preserves_vcf(tmp_path, missing_origin):
    config = load_configuration().modules.formatting
    executable = shutil.which("bcftools")
    if executable is None:
        pytest.skip("Native bcftools is not installed")
    text = (Path(__file__).parent / "fixtures/formatting/harmonised.vcf").read_text()
    if missing_origin:
        name = config.input_contract.provenance_headers.version
        text = "\n".join(line for line in text.splitlines() if not line.startswith("##" + name + "=")) + "\n"
    plain = tmp_path / "study.vcf"
    plain.write_text(text)
    compressed = tmp_path / "study.vcf.gz"
    subprocess.run([executable, "view", "-Oz", "-o", str(compressed), str(plain)], check=True, capture_output=True)
    subprocess.run([executable, "index", "--tbi", str(compressed)], check=True, capture_output=True)
    before = compressed.read_bytes()
    with InputValidationSession() as session:
        if missing_origin:
            with pytest.raises(FormattingError, match="missing PostGWAS"):
                vcf.validate_harmonised_indexed_vcf(compressed, "STUDY", executable, config)
        else:
            result = vcf.validate_harmonised_indexed_vcf(compressed, "STUDY", executable, config)
            assert result["harmonised"]["postgwas_dataset_id"] == "STUDY"
            assert result["indexed"].variant_count == len([line for line in text.splitlines() if not line.startswith("#")])
    records = [record for record in session.records if record.role == "GWAS-VCF provenance"]
    assert len(records) == 1
    assert records[0].status == ("failed" if missing_origin else "passed")
    assert compressed.read_bytes() == before
