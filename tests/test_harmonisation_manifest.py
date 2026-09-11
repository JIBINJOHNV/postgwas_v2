"""Run-manifest terminal-state contracts for dataset failures."""

import inspect
import json
import time
from datetime import datetime
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from postgwas.config import load_configuration
from postgwas.core.paths import configured_output_path
from postgwas.modules.harmonisation.cli import _engine_defaults
from postgwas.modules.harmonisation.sample_sheet import (
    load_harmonisation_sample_sheet,
    to_harmonisation_input,
)
from postgwas.modules.harmonisation.service import (
    PipelineError,
    _combined_dataset_status,
    _finalise_failed_run_manifest,
    _run_post_merge_stages,
    run_harmonisation_pipeline,
)


FIXTURE_DIRECTORY = Path(__file__).parent / "data" / "harmonisation"
SAMPLE_SHEET = FIXTURE_DIRECTORY / "manifest_v2.csv"


class _Logger:
    def __init__(self):
        self.errors = []

    def error(self, message):
        self.errors.append(str(message))


def _pipeline_arguments(tmp_path):
    config = load_configuration()
    defaults = _engine_defaults(
        config,
        resolved_executables={
            "bash": "bash",
            "bcftools": "bcftools",
            "python": "python",
            "tabix": "tabix",
        },
    )
    resource_directory = tmp_path / "resources"
    resource_directory.mkdir()
    (resource_directory / "fixture.txt").write_text("fixture\n", encoding="utf-8")
    row = load_harmonisation_sample_sheet(SAMPLE_SHEET)[0]
    engine_input = to_harmonisation_input(
        row,
        resource_directory=resource_directory,
        output_directory=tmp_path,
        output_layout=defaults["output_layout"],
    )
    harmonisation_directory = Path(engine_input["output_folder"])
    configured_output_path(
        harmonisation_directory,
        defaults["output_layout"]["gwas2vcf_summary"],
        dataset_id=row.dataset_id,
    ).parent.mkdir(parents=True, exist_ok=True)
    return engine_input, defaults


def _chromosome_result():
    return {
        "__dataset__": {
            "status": "OK",
            "completed": ["1"],
            "chromosomes": {},
            "study_decisions": {
                "effect_type": "beta",
                "effect_type_source": "sample_sheet",
                "se_scale": None,
                "se_scale_source": "not_applicable_to_beta",
                "pvalue_type": "raw",
                "pvalue_type_source": "sample_sheet",
                "frequency_type": "effect_allele_frequency",
                "eaf_is_maf_source": "study_level_statistic",
                "strand": "forward",
            },
            "resource_preflight": {
                "require_default_eaf": False,
                "resources_by_chromosome": {"1": {}},
            },
            "chromosome_summaries": {"1": {"rows_out": 1}},
            "reconciliation": {"rows_exported": 1},
        },
        "total_variant_infile": 1,
        "total_variant_read": 1,
        "total_variant_remaining_for_harmonisation": 1,
    }, "GRCh37"


def _qc_frame():
    return pd.DataFrame({
        "section": ["effect", "effect"],
        "metric": ["missing_eaf_count", "variants_removed_invalid_beta_se"],
        "1": [0, 0],
    })


def _assessment():
    return {
        "definition": "test assessment",
        "genome_build": "GRCh37",
        "raw": {
            "num_records": 1,
            "effective_sample_size_reference_quantile_value": 1000.0,
            "effective_sample_size_minimum_threshold": 666.6666667,
            "effective_sample_size_below_minimum_threshold": 0,
            "effective_sample_size_below_minimum_threshold_fraction": 0.0,
        },
        "qc_passed": {
            "num_records": 1,
            "effective_sample_size_below_minimum_threshold": 0,
            "effective_sample_size_below_minimum_threshold_fraction": 0.0,
        },
        "sample_size_reference_quantile": 0.9,
        "sample_size_minimum_fraction_of_reference": 2 / 3,
        "excluded_total": 0,
        "retained_fraction": 1.0,
        "reports": {
            "summary": "summary.tsv",
            "rules": "rules.tsv",
            "json": "assessment.json",
        },
    }


def _cleanup_result():
    return {
        "summary": "summary.tsv",
        "moved": [],
        "removed": [],
        "adapter_input_rows_by_chromosome": {"1": 1},
        "total_adapter_input_rows": 1,
    }


def test_required_merge_partial_downgrades_successful_chromosome_status():
    assert _combined_dataset_status("OK", "PARTIAL") == "PARTIAL"
    assert _combined_dataset_status("PARTIAL", "OK") == "PARTIAL"
    assert _combined_dataset_status("OK", "OK") == "OK"


def test_post_merge_stage_order_and_manifest_finalisation_remain_explicit():
    post_merge_source = inspect.getsource(_run_post_merge_stages)
    ordered_operations = [
        "concat_vcfs_by_build(",
        "run_population_frequency_qc(",
        "finalise_harmonisation_outputs(",
        "run_qc_assessment(",
        "harmonisation_qc_summary_lines(",
        "_resolved_harmonisation_outputs(",
        "_write_combined_log(",
    ]
    positions = [post_merge_source.index(name) for name in ordered_operations]
    assert positions == sorted(positions)
    assert "pd.read_csv(" not in post_merge_source
    assert ".to_csv(" not in post_merge_source

    pipeline_source = inspect.getsource(run_harmonisation_pipeline)
    assert pipeline_source.index("_save_qc_results(") < pipeline_source.index(
        "_run_post_merge_stages("
    )
    assert pipeline_source.index("_run_post_merge_stages(") < (
        pipeline_source.index('manifest["finished"]')
    )
    assert pipeline_source.index('manifest["finished"]') < (
        pipeline_source.rindex("_write_run_manifest(")
    )


def test_population_frequency_inversion_stops_before_output_cleanup(tmp_path):
    engine_input, defaults = _pipeline_arguments(tmp_path)
    merged = {
        "grch37": "input.vcf.gz",
        "grch38": "target.vcf.gz",
        "gwas2vcf": "raw.vcf.gz",
        "merge_failures": [],
        "required_merge_failures": [],
        "optional_merge_failures": [],
        "merge_status": "OK",
    }
    inversion = {
        "status": "frequency_inversion_suspected",
        "decision_reason": "study AF matches 1 - EUR AF",
        "report": str(tmp_path / "frequency_qc.json"),
        "warnings": ["frequency inversion suspected"],
        "external_file_checks": [],
    }

    with patch(
        "postgwas.modules.harmonisation.service.validate_config",
        return_value=(True, []),
    ), patch(
        "postgwas.modules.harmonisation.service.validate_header",
        return_value=(True, []),
    ), patch(
        "postgwas.modules.harmonisation.service.validate_path",
        return_value=lambda _path: None,
    ), patch(
        "postgwas.modules.harmonisation.service.inspect_summary_statistics_file",
        return_value=([], 1),
    ), patch(
        "postgwas.modules.harmonisation.service.harmonise_chromosomes",
        return_value=_chromosome_result(),
    ), patch(
        "postgwas.modules.harmonisation.service.qc_results_to_dataframe",
        return_value=_qc_frame(),
    ), patch(
        "postgwas.modules.harmonisation.service.concat_vcfs_by_build",
        return_value=merged,
    ), patch(
        "postgwas.modules.harmonisation.service.run_population_frequency_qc",
        return_value=inversion,
    ), patch(
        "postgwas.modules.harmonisation.service.finalise_harmonisation_outputs",
    ) as cleanup, patch(
        "postgwas.modules.harmonisation.service._announce",
    ):
        with pytest.raises(
            PipelineError,
            match="likely non-effect-allele frequency column",
        ):
            run_harmonisation_pipeline(engine_input, defaults, threads=1)

    cleanup.assert_not_called()
    manifest_path = configured_output_path(
        engine_input["output_folder"],
        defaults["output_layout"]["run_manifest"],
        dataset_id=engine_input["gwas_outputname"],
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["status"] == "FAILED"
    assert "likely non-effect-allele frequency column" in manifest[
        "failure"
    ]["message"]


@pytest.mark.parametrize(
    "failure_stage",
    ["merge", "cleanup", "qc_assessment", "final_output_validation"],
)
def test_post_chromosome_failures_leave_a_terminal_manifest(tmp_path, failure_stage):
    engine_input, defaults = _pipeline_arguments(tmp_path)
    failure = RuntimeError("%s failed unexpectedly" % failure_stage)

    merge_effect = failure if failure_stage == "merge" else None
    cleanup_effect = failure if failure_stage == "cleanup" else None
    assessment_effect = failure if failure_stage == "qc_assessment" else None
    output_effect = failure if failure_stage == "final_output_validation" else None

    with patch(
        "postgwas.modules.harmonisation.service.validate_config",
        return_value=(True, []),
    ), patch(
        "postgwas.modules.harmonisation.service.validate_header",
        return_value=(True, []),
    ), patch(
        "postgwas.modules.harmonisation.service.validate_path",
        return_value=lambda _path: None,
    ), patch(
        "postgwas.modules.harmonisation.service.inspect_summary_statistics_file",
        return_value=([], 1),
    ), patch(
        "postgwas.modules.harmonisation.service.harmonise_chromosomes",
        return_value=_chromosome_result(),
    ), patch(
        "postgwas.modules.harmonisation.service.qc_results_to_dataframe",
        return_value=_qc_frame(),
    ), patch(
        "postgwas.modules.harmonisation.service.concat_vcfs_by_build",
        return_value={"merge_failures": []},
        side_effect=merge_effect,
    ), patch(
        "postgwas.modules.harmonisation.service.finalise_harmonisation_outputs",
        return_value=_cleanup_result(),
        side_effect=cleanup_effect,
    ), patch(
        "postgwas.modules.harmonisation.service.run_qc_assessment",
        return_value=_assessment(),
        side_effect=assessment_effect,
    ), patch(
        "postgwas.modules.harmonisation.service.harmonisation_qc_summary_lines",
        return_value=[],
    ), patch(
        "postgwas.modules.harmonisation.service._resolved_harmonisation_outputs",
        return_value={"GRCh37": "input.vcf.gz", "GRCh38": "target.vcf.gz"},
        side_effect=output_effect,
    ), patch(
        "postgwas.modules.harmonisation.service._announce",
    ):
        with pytest.raises(RuntimeError) as caught:
            run_harmonisation_pipeline(engine_input, defaults, threads=1)

    assert caught.value is failure
    output_directory = Path(engine_input["output_folder"])
    manifest_path = configured_output_path(
        output_directory,
        defaults["output_layout"]["run_manifest"],
        dataset_id=engine_input["gwas_outputname"],
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["status"] == "FAILED"
    assert manifest["failure"] == {
        "exception_type": "RuntimeError",
        "message": "%s failed unexpectedly" % failure_stage,
    }
    assert manifest["finished"]
    assert manifest["elapsed_seconds"] >= 0
    dataset_log = configured_output_path(
        output_directory,
        defaults["output_layout"]["dataset_log"],
        dataset_id=engine_input["gwas_outputname"],
    ).read_text(encoding="utf-8")
    assert "%s failed unexpectedly" % failure_stage in dataset_log
    assert "status=FAILED." in dataset_log


@pytest.mark.parametrize("keep_intermediate", [False, True])
def test_success_manifest_records_gwas2vcf_cleanup_policy(
    tmp_path, keep_intermediate,
):
    engine_input, defaults = _pipeline_arguments(tmp_path)
    defaults["population_frequency_qc"]["enabled"] = False
    defaults["policies"].setdefault("vcf", {})[
        "keep_gwas2vcf_intermediate"
    ] = keep_intermediate
    raw_gwas2vcf = Path(
        engine_input["output_folder"]
    ) / "study_gwas2vcf_GRCh37_merged.vcf.gz"
    raw_gwas2vcf.write_bytes(b"validated raw VCF\n")
    Path(str(raw_gwas2vcf) + ".tbi").write_bytes(b"validated index\n")
    merged = {
        "grch37": "input.vcf.gz",
        "grch38": "target.vcf.gz",
        "gwas2vcf": str(raw_gwas2vcf),
        "merge_failures": [],
        "required_merge_failures": [],
        "optional_merge_failures": [],
        "merge_status": "OK",
    }
    captured_merge = {}

    def merge_with_provenance(**kwargs):
        captured_merge.update(kwargs)
        return merged

    with patch(
        "postgwas.modules.harmonisation.service.validate_config",
        return_value=(True, []),
    ), patch(
        "postgwas.modules.harmonisation.service.validate_header",
        return_value=(True, []),
    ), patch(
        "postgwas.modules.harmonisation.service.validate_path",
        return_value=lambda _path: None,
    ), patch(
        "postgwas.modules.harmonisation.service.inspect_summary_statistics_file",
        return_value=([], 1),
    ), patch(
        "postgwas.modules.harmonisation.service.harmonise_chromosomes",
        return_value=_chromosome_result(),
    ), patch(
        "postgwas.modules.harmonisation.service.qc_results_to_dataframe",
        return_value=_qc_frame(),
    ), patch(
        "postgwas.modules.harmonisation.service.concat_vcfs_by_build",
        side_effect=merge_with_provenance,
    ), patch(
        "postgwas.modules.harmonisation.service.finalise_harmonisation_outputs",
        return_value=_cleanup_result(),
    ), patch(
        "postgwas.modules.harmonisation.service.run_qc_assessment",
        return_value=_assessment(),
    ) as qc_assessment, patch(
        "postgwas.modules.harmonisation.service.harmonisation_qc_summary_lines",
        return_value=[],
    ), patch(
        "postgwas.modules.harmonisation.service.harmonisation_qc_takeaway_lines",
        return_value=[],
    ), patch(
        "postgwas.modules.harmonisation.service._resolved_harmonisation_outputs",
        return_value={
            "GRCh37": "input.vcf.gz",
            "GRCh38": "target.vcf.gz",
            "gwas2vcf": str(raw_gwas2vcf),
        },
    ), patch(
        "postgwas.modules.harmonisation.service._announce",
    ):
        result = run_harmonisation_pipeline(engine_input, defaults, threads=1)

    manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
    assert manifest["status"] == "OK"
    completed_at = datetime.fromisoformat(manifest["completed_at"])
    assert completed_at.utcoffset() is not None
    assert manifest["vcf_provenance"]["output_directory"] == str(
        tmp_path.resolve()
    )
    assert manifest["vcf_provenance"]["dataset_output_directory"] == str(
        Path(engine_input["output_folder"]).resolve()
    )
    assert captured_merge["provenance"] == manifest["vcf_provenance"]
    qc_arguments = qc_assessment.call_args.kwargs
    assert qc_arguments["vcf_path"] == configured_output_path(
        engine_input["output_folder"],
        defaults["output_layout"]["merged_build_vcf"],
        dataset_id=engine_input["gwas_outputname"],
        build="GRCh37",
    )
    assert "genome_build" not in qc_arguments
    assert qc_arguments["genome_build_header"] == (
        defaults["vcf_processing"]["genome_build_header"]
    )
    assert qc_arguments["supported_genome_builds"] == tuple(
        defaults["vcf_processing"]["target_builds"]
    )
    assert manifest["qc_genome_build"] == "GRCh37"
    assert ("gwas2vcf" in result) is keep_intermediate
    assert raw_gwas2vcf.exists() is keep_intermediate
    assert Path(str(raw_gwas2vcf) + ".tbi").exists() is keep_intermediate
    assert manifest["gwas2vcf_intermediate"]["keep_requested"] is keep_intermediate
    assert manifest["gwas2vcf_intermediate"]["retained"] is keep_intermediate
    assert manifest["gwas2vcf_intermediate"]["retention_reason"] == (
        "explicit_keep_policy"
        if keep_intermediate
        else "successful_run_default_cleanup"
    )
    low_neff = manifest["qc_assessment"]["low_neff"]
    assert low_neff["minimum_fraction_of_reference"] == pytest.approx(2 / 3)
    assert low_neff | {"minimum_fraction_of_reference": None} == {
        "reference_quantile": 0.9,
        "reference_value": 1000.0,
        "minimum_fraction_of_reference": None,
        "minimum_threshold": 666.6666667,
        "raw_below_threshold": 0,
        "raw_below_threshold_fraction": 0.0,
        "qc_passed_below_threshold": 0,
        "qc_passed_below_threshold_fraction": 0.0,
    }


def test_terminal_manifest_preserves_partial_chromosome_status(tmp_path):
    path = tmp_path / "manifest.json"
    manifest = {"status": "PARTIAL", "failed_chromosomes": ["2"]}
    error = PipelineError("chromosome failed", chromosomes=["2"], status="PARTIAL")

    written = _finalise_failed_run_manifest(
        path, manifest, error, time.time(), _Logger(), status="FAILED",
    )

    assert written is True
    saved = json.loads(path.read_text(encoding="utf-8"))
    assert saved["status"] == "PARTIAL"
    assert saved["failed_chromosomes"] == ["2"]
    assert saved["failure"]["exception_type"] == "PipelineError"


def test_keyboard_interrupt_is_recorded_as_interrupted(tmp_path):
    path = tmp_path / "manifest.json"
    manifest = {"status": "running"}

    written = _finalise_failed_run_manifest(
        path, manifest, KeyboardInterrupt(), time.time(), _Logger(),
        status="INTERRUPTED",
    )

    assert written is True
    saved = json.loads(path.read_text(encoding="utf-8"))
    assert saved["status"] == "INTERRUPTED"
    assert saved["failure"] == {
        "exception_type": "KeyboardInterrupt",
        "message": "Run interrupted by the user.",
    }


def test_manifest_write_failure_does_not_replace_original_error(tmp_path):
    logger = _Logger()
    manifest = {"status": "running"}

    with patch(
        "postgwas.modules.harmonisation.service._write_run_manifest",
        side_effect=OSError("manifest is not writable"),
    ):
        written = _finalise_failed_run_manifest(
            tmp_path / "manifest.json",
            manifest,
            RuntimeError("original analysis failure"),
            time.time(),
            logger,
        )

    assert written is False
    assert manifest["status"] == "FAILED"
    assert "original analysis failure" in logger.errors[0]
    assert "manifest is not writable" in logger.errors[0]
    assert "original exception will be preserved" in logger.errors[0]
