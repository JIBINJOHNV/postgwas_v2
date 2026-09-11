"""Pipeline startup summaries replace only repeated static inventory display."""

from argparse import Namespace
from contextlib import nullcontext
import sys

import pytest

from postgwas.config import load_configuration
from postgwas.core.contracts import RunContext
from postgwas.core.input_validation import InputValidationSession
from postgwas.core.preflight import PipelinePreflightEvidence
from postgwas.modules.formatting.reporting import (
    render_formatter_html_report,
    render_formatter_screen_summary,
)
from postgwas.modules.magma.analysis import MagmaReferencePreflight
from postgwas.modules.magma.service import MagmaPipelineResources
from postgwas.pipeline import runners


@pytest.mark.parametrize("shared", [False, True])
def test_static_resource_display_retains_warnings_and_progress(
    tmp_path, monkeypatch, shared,
):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_text("fixture; VCF evidence supplied below\n")
    configuration = load_configuration()
    resources = MagmaPipelineResources(
        configuration=configuration,
        reference=MagmaReferencePreflight(
            ld_reference_prefix=str(tmp_path / "panel"),
            executable=sys.executable,
            version="1.10",
            gene_set_plans={"positional": {"status": "not_requested"}},
            analysis_scope={},
        ),
        include_gene_sets=False,
        primary_gene_location_metrics={
            "path": str(tmp_path / "genes.loc"),
            "genes": 3,
            "alternate_unique_ids": 2,
            "ambiguous_alternate_ids": 1,
        },
        file_identities=(),
        log_events=(),
    )
    context = RunContext(validations={
        "magma": PipelinePreflightEvidence("magma", {}, resources),
    })
    started = []
    completed = {}
    args = Namespace(
        modules=["magmacovar"], vcf=str(vcf), format=None,
        output_directory=str(tmp_path / "output"), run_config=None,
        _pipeline_stage_progress=Namespace(
            start=started.append,
            complete=lambda number, **values: completed.update({number: values}),
        ),
        _pipeline_progress_plan={
            "kind": "magmacovar",
            "magma_stage_numbers": {"vcf": 1, "ld_reference": 2, "gene_location": 3},
        },
    )
    calls = []

    def current_vcf(*_):
        calls.append("current_vcf")
        return {
            "indexed": Namespace(variant_count=100),
            "harmonised": {"genome_build": "GRCh37", "postgwas_dataset_id": "STUDY"},
        }

    def identifiers(args, *_):
        calls.append("identifier_contract")
        args.variant_id_observations = {
            "magma": {"variant_id_type": "unique", "variants": 250},
        }

    def formatter(*_):
        calls.append("formatter")
        return {"magma": {"rows_out": 80}}

    monkeypatch.setattr(runners, "_validate_current_pipeline_vcf", current_vcf)
    monkeypatch.setattr(runners, "configure_reference_variant_identifiers", identifiers)
    monkeypatch.setattr(
        "postgwas.modules.formatting.service.run_formatter_direct", formatter,
    )
    with InputValidationSession() if shared else nullcontext():
        assert runners.run_formatter_runner(args, context) == {"magma": {"rows_out": 80}}

    assert calls == ["current_vcf", "identifier_contract", "formatter"]
    assert started == [1, 2, 3]
    assert list(completed) == [1, 2, 3]
    fields = {
        number: {field[1]: field[2] for field in result["outcome_fields"] if len(field) == 3}
        for number, result in completed.items()
    }
    # The current VCF may have changed after filtering or imputation.
    assert fields[1]["Total variants"] == 100
    assert fields[1]["Genome build inferred from VCF"] == "GRCh37"
    assert "not independently verifiable" in fields[2]["Build and population provenance"]
    assert fields[3]["Duplicated alternate IDs"] == 1
    assert ("Total BIM variants" in fields[2]) is not shared
    # Gene content/policy summaries are not yet all present in the shared audit.
    assert fields[3]["Reference genes"] == 3


def test_formatter_pipeline_omits_repeat_scan_but_keeps_final_provenance(tmp_path):
    configuration = load_configuration()
    observation = {
        "gcta_gene": {
            "variant_id_type": "rsid", "variants": 250,
            "bim_file": str(tmp_path / "reference.bim"), "consumers": ["GCTA gene"],
        },
    }
    results = {"gcta_gene": {
        "rows_in": 100, "rows_out": 80, "rows_excluded": 20,
        "identifier_uniqueness_validated": True,
        "identifier_duplicate_policy": "error", "variant_id_type": "rsid",
    }}
    values = dict(
        dataset_id="STUDY", vcf=tmp_path / "study.vcf.gz", study_design=None,
        results=results, log_path=tmp_path / "formatter.log",
        html_report_path=tmp_path / "formatter.html",
        case_count_column="N_CASE", control_count_column="N_CONTROL",
        label_width=52, variant_id_observations=observation,
        output_destinations={"gcta_gene": str(tmp_path / "study.ma")},
        input_evidence={"genome_build": "GRCh37", "postgwas_dataset_id": "STUDY"},
        minimum_p_value=configuration.modules.formatting.minimum_p_value,
    )
    direct = render_formatter_screen_summary(**values)
    with InputValidationSession():
        pipeline = render_formatter_screen_summary(**values)
        html = render_formatter_html_report(
            dataset_id="STUDY", input_vcf=values["vcf"], output_directory=tmp_path,
            report_path=values["html_report_path"], study_design=None, results=results,
            schema_reports={"gcta_gene": {}},
            output_destinations=values["output_destinations"], log_path=values["log_path"],
            resolved_config_path=tmp_path / "config.yaml",
            completion_manifest_path=tmp_path / "manifest.json",
            case_count_column="N_CASE", control_count_column="N_CONTROL",
            variant_id_observations=observation, resolved_bcftools=sys.executable,
            module=configuration.modules.formatting,
        )
    assert "detected after scanning 250 BIM variants" in direct
    assert "detected after scanning" not in pipeline
    for report in (direct, pipeline):
        assert "BIM identifier convention" in report
        assert "reference.bim" in report
        assert "not assessed by formatter" in report
        assert "Final identifier validation" in report
        assert "80 retained identifiers are unique" in report
        assert "study.ma" in report
        assert "GRCh37" in report
    assert "Reference variants scanned" in html
    assert "250" in html
    assert "Final uniqueness" in html


def test_startup_resource_presentation_does_not_discard_cautions():
    fields = [
        ("analysis", "Reference"), ("count", "Rows", 10),
        ("warning", "Provenance", "declared only"),
        ("error", "Failure", "incompatible"),
        ("loss", "Excluded variants", 2),
    ]
    assert runners._startup_resource_outcome_fields(fields) is fields
    with InputValidationSession():
        assert runners._startup_resource_outcome_fields(fields) == fields[2:]
    assert len(fields) == 5


@pytest.mark.parametrize("display_fields", [None, [], [("warning", "Provenance", "declared")]])
@pytest.mark.parametrize("deferred", [False, True])
def test_gcta_compact_screen_keeps_full_canonical_stage_evidence(display_fields, deferred):
    from postgwas.modules.gcta_gene.stages import complete_pipeline_stage

    displayed = []
    logged = []
    fields = [("count", "Variants", 250), ("warning", "Provenance", "declared")]
    args = Namespace(
        _pipeline_progress_plan={
            "kind": "gcta_gene", "stage_numbers": {"reference": 1},
            "stages": ("Reference",),
        },
        _pipeline_stage_progress=Namespace(
            complete=lambda number, **values: displayed.append((number, values)),
        ),
    )
    if deferred:
        args._pipeline_progress_completion = lambda: None
    logger = Namespace(record=lambda *event, **values: logged.append((event, values)))
    assert complete_pipeline_stage(
        args, "reference", outcome_fields=fields,
        screen_outcome_fields=display_fields, logger=logger,
    ) == 1
    expected_display = fields if display_fields is None else display_fields
    if deferred:
        assert not displayed
        assert args._pipeline_stage_completion["outcome_fields"] == expected_display
    else:
        assert displayed[0][1]["outcome_fields"] == expected_display
    assert logged[0][1]["details"] == {"Variants": 250, "Provenance": "declared"}
    assert logged[1][1]["status"] == (
        "VALIDATED_PENDING_PIPELINE_CHECKPOINT" if deferred else "COMPLETED"
    )
