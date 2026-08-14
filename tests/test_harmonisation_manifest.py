"""Run-manifest terminal-state contracts for dataset failures."""

import json
import time
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
        "raw": {"num_records": 1},
        "qc_passed": {"num_records": 1},
        "excluded_total": 0,
        "retained_fraction": 1.0,
        "reports": {
            "summary": "summary.tsv",
            "rules": "rules.tsv",
            "json": "assessment.json",
        },
    }


def test_required_merge_partial_downgrades_successful_chromosome_status():
    assert _combined_dataset_status("OK", "PARTIAL") == "PARTIAL"
    assert _combined_dataset_status("PARTIAL", "OK") == "PARTIAL"
    assert _combined_dataset_status("OK", "OK") == "OK"


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
        return_value={"moved": [], "removed": []},
        side_effect=cleanup_effect,
    ), patch(
        "postgwas.modules.harmonisation.service.pd.read_csv",
        return_value=pd.DataFrame({
            "chromosome": ["1"], "key": ["snp_id_col"], "num_rows": [1],
        }),
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
