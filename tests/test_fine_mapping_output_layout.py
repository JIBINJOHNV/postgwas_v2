"""Regression tests for the categorized fine-mapping output layout."""

from __future__ import annotations

import pytest

from postgwas.config import load_configuration, load_run_configuration_for_module
from postgwas.core.errors import ConfigurationError
from postgwas.modules.fine_mapping.engines.susie.adapter import merge_workers
from postgwas.modules.fine_mapping.output_layout import (
    cleanup_successful_workers,
    resolve_output_paths,
)


def _layout():
    return load_configuration().modules.fine_mapping.output_layout.model_dump()


def test_canonical_layout_separates_published_qc_handoff_and_intermediate_paths(
    tmp_path,
):
    paths = resolve_output_paths(tmp_path, _layout(), "study")

    assert paths["combined_results_directory"] == (
        tmp_path / "results" / "combined_results"
    )
    assert paths["susie_qc_file"] == (
        tmp_path / "quality_control" / "study_susie_locus_qc.tsv"
    )
    assert paths["downstream_flames_directory"] == (
        tmp_path / "downstream_inputs" / "flames"
    )
    assert paths["pipeline_log_file"] == (
        tmp_path / "run_metadata" / "pipeline_summary.log"
    )
    assert paths["workers_directory"] == (
        tmp_path / "intermediate_files" / "susie_workers"
    )
    assert paths["fitted_models_directory"] == (
        tmp_path / "intermediate_files" / "fitted_models"
    )


def test_layout_rejects_a_path_that_escapes_the_run_directory(tmp_path):
    config_file = tmp_path / "fine_mapping.yaml"
    config_file.write_text(
        "output_layout:\n  workers_directory: ../workers\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="must be a non-empty path"):
        load_run_configuration_for_module("fine_mapping", config_file)


def test_worker_artifacts_publish_to_categories_then_duplicate_workers_are_removed(
    tmp_path,
):
    paths = resolve_output_paths(tmp_path, _layout(), "study")
    worker = paths["workers_directory"] / "worker_1"
    (worker / "plots").mkdir(parents=True)
    (worker / "rds_files").mkdir()
    (worker / "plots" / "locus.png").write_bytes(b"plot")
    (worker / "rds_files" / "locus.rds").write_bytes(b"model")
    (worker / "study_run_susie.log").write_text("worker log\n", encoding="utf-8")

    merge_workers(
        worker_dirs=[worker],
        output_paths=paths,
        sample_id="study",
        recovery_audit_filename="susie_recovery_audit.tsv",
        index_filename="indexfile.txt",
        annotation_prefix="annotated_",
    )

    assert (paths["diagnostic_plots_directory"] / "locus.png").is_file()
    assert (paths["fitted_models_directory"] / "locus.rds").is_file()
    assert not (tmp_path / "plots").exists()
    assert not paths["primary_credible_sets_directory"].exists()
    assert not paths["locus_logs_directory"].exists()
    assert not paths["ld_diagnostics_directory"].exists()
    assert not paths["primary_flames_annotations_directory"].exists()
    assert cleanup_successful_workers(
        paths["workers_directory"], enabled=True
    )
    assert not paths["workers_directory"].exists()


def test_worker_retention_policy_can_keep_intermediate_workers(tmp_path):
    workers = tmp_path / "intermediate_files" / "susie_workers"
    workers.mkdir(parents=True)

    assert not cleanup_successful_workers(workers, enabled=False)
    assert workers.is_dir()
