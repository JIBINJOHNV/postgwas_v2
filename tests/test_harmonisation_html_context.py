"""HTML evidence enrichment reuses existing inputs and QC, without analysis."""

from argparse import Namespace
from copy import deepcopy
import json

import pytest

from postgwas.config import load_configuration
from postgwas.modules.harmonisation import cli
from postgwas.modules.harmonisation.sample_sheet import HarmonisationSampleSheetRow


def _row(dataset_id="study"):
    return HarmonisationSampleSheetRow(
        config_version=2,
        dataset_id=dataset_id,
        input_file="/input/study.tsv.gz",
        chromosome_column="#chrom",
        position_column="position",
        effect_allele_column="ALT",
        other_allele_column="REF",
        effect_allele_frequency_column="alt_allele_freq",
        effect_column="effect",
        p_value_column="P",
        control_count=1000,
        imputation_info_column="INFO",
    )


def _configuration(tmp_path, **overrides):
    return load_configuration(cli_overrides={
        "resources.root": str(tmp_path),
        "run.output_directory": str(tmp_path / "output"),
        "execution.retries": 0,
        "logging.show_screen": False,
        **overrides,
    })


@pytest.mark.parametrize("requested", [False, True])
def test_html_context_preserves_exact_mappings_and_resolved_choices(tmp_path, requested):
    config = _configuration(tmp_path, **{
        "modules.harmonisation.concordance_validation.enabled": requested,
    })
    row = _row()
    original = row.model_dump(mode="json")

    context = cli._html_report_context(config, row, sample_sheet="studies.csv")

    assert context["input_mapping"] == original
    assert context["input_mapping"]["chromosome_column"] == "#chrom"
    assert context["input_mapping"]["z_score_column"] is None
    assert context["concordance_requested"] is requested
    assert context["configuration"]["harmonisation"] == (
        config.modules.harmonisation.model_dump(mode="json")
    )
    assert context["sample_sheet"] == "studies.csv"
    assert row.model_dump(mode="json") == original
    json.dumps(context)


def test_summary_reuses_one_qc_read_without_changing_manifest_or_csv(tmp_path, monkeypatch):
    assessment = {
        "raw": {"num_snps": 10, "num_records": 11},
        "qc_passed": {"num_snps": 8, "num_records": 9},
        "rules": [{"name": "custom recorded rule", "failed": 2}],
    }
    manifest = {"qc_assessment": {"raw_variants": 11}}
    original = deepcopy(manifest)
    reads = []

    def read_assessment(*args, **kwargs):
        reads.append((args, kwargs))
        return assessment

    monkeypatch.setattr(cli, "_read_qc_assessment_report", read_assessment)
    context = {}
    record = cli._run_summary_record(
        _row(), status="FAILED", manifest_path=tmp_path / "manifest.json",
        manifest=manifest, strand_reference_panel="panel",
        strand_reference_population="population", comparison_af_panel="panel",
        comparison_af_population="population", report_context=context,
    )

    assert len(reads) == 1
    assert context["qc_assessment"] is assessment
    assert record["final_vcf_snps"] == 10
    assert record["qc_passed_snps"] == 8
    assert manifest == original
    assert set(record) == set(cli._RUN_SUMMARY_FIELDS)
    assert "report_context" not in record


def test_preflight_failure_keeps_input_context_without_analysis(tmp_path, monkeypatch):
    config = _configuration(tmp_path)
    reports = []
    monkeypatch.setattr(cli, "_write_html_reports", lambda **kwargs: reports.append(kwargs))

    cli._write_preflight_failure_run_summary(
        config, [_row()], "required plugin unavailable", sample_sheet="studies.csv",
    )

    context = reports[0]["manifests"]["study"]["report_context"]
    assert context["input_mapping"]["chromosome_column"] == "#chrom"
    assert context["sample_sheet"] == "studies.csv"
    assert "qc_assessment" not in context
    assert reports[0]["records"]["study"]["status"] == "PREFLIGHT_FAILED"


def test_initial_and_failed_datasets_retain_their_own_report_context(tmp_path, monkeypatch):
    config = _configuration(tmp_path)
    reports = []
    monkeypatch.setattr(
        cli, "_write_html_reports", lambda **kwargs: reports.append(deepcopy(kwargs)),
    )

    def fail_engine(**kwargs):
        raise cli.PipelineError("synthetic failure before scientific work")

    monkeypatch.setattr(cli, "run_harmonisation_pipeline", fail_engine)
    with pytest.raises(cli.PipelineError, match="All harmonisation datasets failed"):
        cli._run_validated_rows(Namespace(), config, [_row("one"), _row("two")])

    for dataset_id in ("one", "two"):
        assert reports[0]["records"][dataset_id]["status"] == "NOT_RUN"
        initial_context = reports[0]["manifests"][dataset_id]["report_context"]
        final_context = reports[-1]["manifests"][dataset_id]["report_context"]
        assert initial_context["input_mapping"]["dataset_id"] == dataset_id
        assert final_context["input_mapping"] == initial_context["input_mapping"]
        assert final_context["concordance_requested"] is False
        assert final_context["qc_assessment"] == {}
        assert reports[-1]["records"][dataset_id]["status"] == "FAILED"
