"""Public study-VCF origin gates must precede filtering and QC analysis."""

from argparse import Namespace
from pathlib import Path
import re
from unittest.mock import Mock

import pytest

from postgwas.config import load_configuration
from postgwas.core.paths import configured_output_path
from postgwas.modules.filtering import service as filtering_service
from postgwas.modules.filtering import sumstat_filter
from postgwas.modules.qc_summary import service as qc_service


def _header() -> str:
    fixture = Path(__file__).parent / "fixtures" / "qc_summary" / "mini_merged.vcf"
    lines = []
    for line in fixture.read_text(encoding="utf-8").splitlines():
        lines.append(line)
        if line.startswith("#CHROM\t"):
            break
    return "\n".join(lines) + "\n"


def _direct_boundary(tmp_path, monkeypatch, command, header):
    vcf = tmp_path / "study.vcf"
    vcf.write_text(header, encoding="utf-8")
    args = Namespace(
        vcf=vcf,
        dataset_id="study",
        output_directory=tmp_path / "output",
        threads=1,
    )
    if command == "filtering":
        configuration = filtering_service.resolve_filtering_configuration(args)
        monkeypatch.setattr(
            filtering_service, "_resolve_filtering_executable", lambda value, _: value,
        )
        reader = Mock(return_value=header)
        monkeypatch.setattr(sumstat_filter, "read_vcf_header", reader)
        next_operation = Mock(side_effect=RuntimeError("after provenance gate"))
        monkeypatch.setattr(sumstat_filter, "run_cmd", next_operation)
        log = configured_output_path(
            args.output_directory,
            configuration.modules.filtering.output_layout.preflight_log,
            dataset_id=args.dataset_id,
        )
        runner = filtering_service.run_sumstat_filter_direct
    else:
        configuration = qc_service.resolve_qc_summary_configuration(args)
        monkeypatch.setattr(qc_service, "resolve_executable", lambda value, *a, **k: value)
        reader = Mock(return_value=header)
        monkeypatch.setattr(qc_service, "read_vcf_header", reader)
        next_operation = Mock(side_effect=RuntimeError("after provenance gate"))
        monkeypatch.setattr(qc_service, "run_checked_command", next_operation)
        log = qc_service._output_paths(
            args.output_directory, args.dataset_id, "GRCh37",
            configuration.modules.qc_summary,
        )["log"]
        runner = qc_service.run_qc_summary_direct
    configuration.logging.show_progress = False
    return args, configuration, reader, next_operation, log, runner


@pytest.mark.parametrize("command", ["filtering", "qc"])
@pytest.mark.parametrize("key", ["version", "dataset_id", "status"])
@pytest.mark.parametrize("invalid", ["missing", "empty", "duplicate", "malformed"])
def test_direct_origin_failure_precedes_analysis_and_preserves_outputs(
    tmp_path, monkeypatch, command, key, invalid,
):
    names = load_configuration().modules.formatting.input_contract.provenance_headers
    name = names.model_dump()[key]
    header = _header()
    original = next(line for line in header.splitlines() if line.startswith("##%s=" % name))
    replacement = {
        "missing": "",
        "empty": '##%s="   "' % name,
        "duplicate": original + "\n" + original,
        "malformed": '##%s="unterminated' % name,
    }[invalid]
    header = header.replace(original, replacement)
    args, configuration, reader, next_operation, log, runner = _direct_boundary(
        tmp_path, monkeypatch, command, header,
    )
    existing = args.output_directory / "existing_scientific_output.txt"
    existing.parent.mkdir(parents=True)
    existing.write_text("previous completed result", encoding="utf-8")
    assessment = Mock(side_effect=AssertionError("QC analysis must not start"))
    audit = Mock(side_effect=AssertionError("filtering audit must not start"))
    monkeypatch.setattr(qc_service, "run_qc_assessment", assessment)
    monkeypatch.setattr(sumstat_filter, "_collect_filter_reason_statistics", audit)

    with pytest.raises((ValueError, qc_service.VcfAssessmentError), match=name) as error:
        runner(args, configuration=configuration)

    reader.assert_called_once()
    next_operation.assert_not_called()
    assessment.assert_not_called()
    audit.assert_not_called()
    assert existing.read_text(encoding="utf-8") == "previous completed result"
    assert set(path for path in args.output_directory.rglob("*") if path.is_file()) == {
        existing, log,
    }
    text = log.read_text(encoding="utf-8")
    assert "FAILED" in text
    assert name in text
    assert str(args.vcf) in text
    # Canonical QC logs repeat their timestamp/level prefix on wrapped lines.
    unwrapped = " ".join(re.sub(r"(?m)^\[[^\n]+?\] RUN\s+", "", text).split())
    for guidance in (
        "Re-run PostGWAS harmonisation",
        "do not add provenance headers manually",
    ):
        assert guidance in str(error.value)
        assert guidance in unwrapped


@pytest.mark.parametrize("command", ["filtering", "qc"])
def test_valid_origin_uses_existing_header_and_reaches_next_operation(
    tmp_path, monkeypatch, command,
):
    args, configuration, reader, next_operation, log, runner = _direct_boundary(
        tmp_path, monkeypatch, command, _header(),
    )
    if command == "filtering":
        # Stop at the first count request; the real counter catches ordinary
        # subprocess failures to try its configured unindexed-input fallback.
        next_operation.side_effect = KeyboardInterrupt("after provenance gate")
    with pytest.raises((RuntimeError, KeyboardInterrupt), match="after provenance gate"):
        runner(args, configuration=configuration)

    reader.assert_called_once()
    next_operation.assert_called_once()
    assert "test-version" in log.read_text(encoding="utf-8")


def test_internal_qc_assessment_does_not_apply_public_origin_gate(tmp_path, monkeypatch):
    configuration = load_configuration()
    delegate = Mock(return_value={"internal": "assessment"})
    gate = Mock(side_effect=AssertionError("public provenance gate reached"))
    monkeypatch.setattr(qc_service, "run_vcf_qc_assessment", delegate)
    monkeypatch.setattr(qc_service, "validate_postgwas_vcf_provenance", gate)

    result = qc_service.run_qc_assessment(
        vcf_path=tmp_path / "internal.vcf",
        output_directory=tmp_path,
        dataset_id="study",
        external_af_name=configuration.modules.qc_summary.reference_af_column,
        configuration=configuration.modules.qc_summary,
        bcftools_bin=configuration.resources.executables.bcftools,
        genome_build_header=(
            configuration.modules.harmonisation.vcf_processing.genome_build_header
        ),
        supported_genome_builds=tuple(configuration.resources.genomes),
        threads=1,
    )

    assert result == {"internal": "assessment"}
    delegate.assert_called_once()
    assert delegate.call_args.kwargs["provenance_headers"] is None
    gate.assert_not_called()
