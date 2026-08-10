"""Regression tests for the PostGWAS boundary around upstream PoPS v0.2."""

from __future__ import annotations

import argparse
import logging
from argparse import Namespace
from pathlib import Path

import numpy as np
import pytest
import yaml

from postgwas.config import load_configuration, load_module_configuration
from postgwas.modules.pops.cli import build_parser
from postgwas.modules.pops.errors import PopsError
from postgwas.modules.pops.service import (
    _resolved_configuration,
    _require_pops_runtime,
    preflight_pops_pipeline,
    run_pops_direct,
)
from postgwas.pipeline.registry import REGISTRY


def _pops_resources(tmp_path: Path, *, custom_target: bool = False) -> dict[str, Path]:
    genes = ["ENSG%03d" % index for index in range(1, 11)]
    annotation = tmp_path / "genes.tsv"
    annotation.write_text(
        "ENSGID\tCHR\tTSS\n"
        + "".join(
            "%s\t%s\t%d\n" % (gene, 1 + index % 2, 1000 + index)
            for index, gene in enumerate(genes)
        ),
        encoding="utf-8",
    )
    feature_prefix = tmp_path / "features"
    Path(str(feature_prefix) + ".rows.txt").write_text(
        "\n".join(genes) + "\n", encoding="utf-8",
    )
    Path(str(feature_prefix) + ".cols.0.txt").write_text(
        "feature_a\nfeature_b\n", encoding="utf-8",
    )
    np.save(str(feature_prefix) + ".mat.0.npy", np.arange(20, dtype=float).reshape(10, 2))

    resources = {"annotation": annotation, "feature_prefix": feature_prefix}
    if custom_target:
        target = tmp_path / "target.tsv"
        target.write_text(
            "ENSGID\tScore\n"
            + "".join("%s\t%.2f\n" % (gene, index / 10) for index, gene in enumerate(genes)),
            encoding="utf-8",
        )
        resources["target"] = target
    else:
        magma_prefix = tmp_path / "magma"
        Path(str(magma_prefix) + ".genes.out").write_text(
            "GENE\tZSTAT\n"
            + "".join("%s\t%.2f\n" % (gene, index / 10) for index, gene in enumerate(genes)),
            encoding="utf-8",
        )
        Path(str(magma_prefix) + ".genes.raw").write_text(
            "non-empty upstream MAGMA raw placeholder\n", encoding="utf-8",
        )
        resources["magma_prefix"] = magma_prefix
    return resources


def _arguments(tmp_path: Path, resources: dict[str, Path], **overrides) -> Namespace:
    values = {
        "dataset_id": "STUDY",
        "output_directory": str(tmp_path / "output"),
        "genome_build": "GRCh37",
        "feature_matrix_prefix": str(resources["feature_prefix"]),
        "feature_matrix_chunks": 1,
        "pops_gene_location_file": str(resources["annotation"]),
    }
    if "magma_prefix" in resources:
        values["magma_association_prefix"] = str(resources["magma_prefix"])
    if "target" in resources:
        values["target_score_file"] = str(resources["target"])
    values.update(overrides)
    return Namespace(**values)


def _write_upstream_outputs(config: dict, progress_callback=None) -> None:
    _emit_upstream_progress(progress_callback)
    prefix = config["out_prefix"]
    Path(prefix + ".preds").write_text(
        "ENSGID\tPoPS_Score\ttraining_gene\n"
        "ENSG001\t1\tTrue\nENSG002\t2\tTrue\n",
        encoding="utf-8",
    )
    Path(prefix + ".coefs").write_text("parameter\tbeta\nfeature_a\t1\n", encoding="utf-8")
    Path(prefix + ".marginals").write_text(
        "feature\tpval\tselected\nfeature_a\t0.01\tTrue\n",
        encoding="utf-8",
    )
    Path(prefix + ".log").write_text("completed\n", encoding="utf-8")


def _mock_upstream(monkeypatch, function) -> None:
    from postgwas.modules.pops.pops import get_pops_args

    monkeypatch.setattr(
        "postgwas.modules.pops.service._load_upstream_entrypoints",
        lambda: (get_pops_args, function),
    )


def _emit_upstream_progress(progress_callback) -> None:
    """Exercise the same five scientific boundaries emitted by upstream PoPS."""
    if progress_callback is None:
        return
    progress_callback("target_loading", {
        "target_source": "MAGMA gene Z statistics",
        "target_genes": 10,
        "covariates": 6,
        "error_covariance": "used",
    })
    progress_callback("covariate_adjustment", {
        "status": "applied",
        "genes": 10,
        "covariates": 6,
        "hla_policy": "excluded",
    })
    progress_callback("feature_selection", {
        "strategy": "marginal association p-value filtering",
        "genes": 10,
        "tested_features": 2,
        "selected_features": 1,
    })
    progress_callback("model_fitting", {
        "method": "ridge",
        "training_genes": 10,
        "model_features": 1,
        "selected_cv_alpha": 1.0,
    })
    progress_callback("gene_scoring", {
        "genes_scored": 2,
        "model_features": 1,
        "outputs_written": 3,
    })


def test_cli_has_no_independent_defaults():
    parser = build_parser()
    destinations = {
        "magma_association_prefix", "feature_matrix_prefix",
        "feature_matrix_chunks", "pops_gene_location_file", "genome_build",
        "use_magma_covariates", "use_magma_error_covariance",
        "feature_selection_p_cutoff", "method", "save_matrix_files", "verbose",
        "run_config", "resume", "overwrite", "dataset_id", "output_directory",
        "threads", "memory_gb", "seed",
    }
    for action in parser._actions:
        if action.dest in destinations:
            assert action.default == argparse.SUPPRESS
    assert vars(parser.parse_args([])) == {}


def test_missing_scikit_learn_has_actionable_error(monkeypatch):
    monkeypatch.setattr(
        "postgwas.modules.pops.service.find_spec",
        lambda package: None,
    )
    with pytest.raises(PopsError, match="requires scikit-learn"):
        _require_pops_runtime()


def test_configuration_uses_upstream_method_names_and_validates_ranges(tmp_path):
    assert load_configuration().modules.pops.method == "ridge"
    invalid_method = tmp_path / "invalid_method.yaml"
    invalid_method.write_text("method: linear\n", encoding="utf-8")
    with pytest.raises(Exception, match="method"):
        load_module_configuration("pops", invalid_method)
    invalid_count = tmp_path / "invalid_count.yaml"
    invalid_count.write_text("feature_matrix_chunks: 0\n", encoding="utf-8")
    with pytest.raises(Exception, match="feature_matrix_chunks"):
        load_module_configuration("pops", invalid_count)
    invalid_reporting = tmp_path / "invalid_reporting.yaml"
    invalid_reporting.write_text(
        "reporting:\n  top_gene_count: 0\n", encoding="utf-8",
    )
    with pytest.raises(Exception, match="top_gene_count"):
        load_module_configuration("pops", invalid_reporting)

    conflicting_targets = tmp_path / "conflicting_targets.yaml"
    conflicting_targets.write_text(
        "magma_association_prefix: magma/study\n"
        "target_score_file: targets.tsv\n",
        encoding="utf-8",
    )
    with pytest.raises(Exception, match="mutually exclusive"):
        load_module_configuration("pops", conflicting_targets)


def test_configuration_precedence_and_compute_seed_are_resolved_once(tmp_path):
    config = tmp_path / "pops.yaml"
    config.write_text(
        "genome_build: GRCh37\nmethod: lasso\nfeature_matrix_chunks: 3\n",
        encoding="utf-8",
    )
    resolved = _resolved_configuration(
        Namespace(run_config=str(config), method="linreg", seed=17)
    )
    assert resolved.modules.pops.method == "linreg"
    assert resolved.modules.pops.feature_matrix_chunks == 3
    assert resolved.execution.random_seed == 17


def test_pipeline_registry_runs_pops_resource_preflight(tmp_path):
    assert REGISTRY.get("pops").preflight.endswith(":preflight_pops_pipeline")
    resources = _pops_resources(tmp_path)
    args = _arguments(tmp_path, resources)
    del args.magma_association_prefix
    preflight_pops_pipeline(args)
    Path(str(resources["feature_prefix"]) + ".mat.0.npy").unlink()
    with pytest.raises(PopsError, match="feature matrix chunk"):
        preflight_pops_pipeline(args)


def test_pipeline_preflight_rejects_magma_build_mismatch(tmp_path):
    resources = _pops_resources(tmp_path)
    args = _arguments(tmp_path, resources)
    del args.magma_association_prefix
    del args.genome_build
    args.pops_genome_build = "GRCh38"
    with pytest.raises(PopsError, match="does not match pipeline MAGMA"):
        preflight_pops_pipeline(args)


def test_standalone_uses_canonical_seed_and_does_not_require_context(tmp_path, monkeypatch):
    resources = _pops_resources(tmp_path)
    args = _arguments(tmp_path, resources)
    observed = {}

    def fake_pops_main(config, progress_callback=None):
        observed.update(config)
        _emit_upstream_progress(progress_callback)
        _write_upstream_outputs(config)

    _mock_upstream(monkeypatch, fake_pops_main)
    result = run_pops_direct(args)
    assert result["status"] == "success"
    assert Path(result["pops_file"]).is_file()
    assert observed["random_seed"] == load_configuration().execution.random_seed
    assert Path(result["completion_manifest"]).is_file()
    assert result["summary"]["genes_scored"] == 2
    assert result["summary"]["selected_features"] == 1
    assert result["summary"]["top_genes"][0]["gene_id"] == "ENSG002"
    assert yaml.safe_load(Path(result["completion_manifest"]).read_text())["status"] == "COMPLETED"


def test_terminal_summary_explains_ranking_and_incomplete_target_coverage(
    tmp_path, monkeypatch, capsys,
):
    resources = _pops_resources(tmp_path)
    magma_output = Path(str(resources["magma_prefix"]) + ".genes.out")
    magma_output.write_text(
        "GENE\tZSTAT\nENSG001\t1\nENSG002\t2\n", encoding="utf-8",
    )
    args = _arguments(tmp_path, resources)
    run_config = tmp_path / "pops_warning.yaml"
    run_config.write_text("minimum_gene_count: 2\n", encoding="utf-8")
    args.run_config = str(run_config)
    _mock_upstream(monkeypatch, _write_upstream_outputs)

    result = run_pops_direct(args)

    screen = capsys.readouterr().out
    assert "PoPS analysis progress" in screen
    assert "Completed 7/7 · Validate and publish PoPS results" in screen
    assert "Target genes loaded" in screen
    assert "Features selected" in screen
    assert "Training genes" in screen
    assert "PoPS gene-prioritisation summary" in screen
    assert "Genes assigned PoPS scores" in screen
    assert "None. PoPS scores are relative rankings, not p-values." in screen
    assert "COMPLETED WITH SCIENTIFIC WARNINGS" in screen
    assert "Target scores cover 2 of 10 compatible genes" in screen
    assert result["summary"]["warnings"]


def test_real_upstream_reports_ordered_stages_without_internal_log_noise(
    tmp_path, capsys,
):
    resources = _pops_resources(tmp_path, custom_target=True)
    args = _arguments(tmp_path, resources)

    result = run_pops_direct(args)

    screen = capsys.readouterr().out
    completed = [
        screen.index("Completed %d/7" % step)
        for step in range(1, 8)
    ]
    assert completed == sorted(completed)
    assert "Target source" in screen
    assert "target-score" in screen
    assert "Selection strategy" in screen
    assert "Prediction model" in screen
    assert "Config dict =" not in screen
    assert "Computing marginal association table" not in screen
    assert Path(result["pops_file"]).is_file()
    upstream_log = next(
        Path(path) for path in result["published_files"] if path.endswith(".log")
    )
    assert "Config dict =" in upstream_log.read_text(encoding="utf-8")


def test_context_and_root_logging_are_restored_after_upstream_run(tmp_path, monkeypatch):
    resources = _pops_resources(tmp_path)
    args = _arguments(tmp_path, resources)
    root = logging.getLogger()
    original_handlers = list(root.handlers)
    original_level = root.level

    def fake_pops_main(config, progress_callback=None):
        root.handlers.clear()
        root.addHandler(logging.StreamHandler())
        root.setLevel(logging.DEBUG)
        _emit_upstream_progress(progress_callback)
        _write_upstream_outputs(config)

    _mock_upstream(monkeypatch, fake_pops_main)
    context = {}
    result = run_pops_direct(args, context)
    assert context["pops_output"] == result["pops_file"]
    assert root.handlers == original_handlers
    assert root.level == original_level


def test_numpy_target_covariance_is_validated_and_converted_for_upstream(tmp_path, monkeypatch):
    resources = _pops_resources(tmp_path, custom_target=True)
    covariance = tmp_path / "covariance.npy"
    np.save(covariance, np.eye(10))
    args = _arguments(
        tmp_path,
        resources,
        target_error_covariance_file=str(covariance),
    )
    observed = {}

    def fake_pops_main(config, progress_callback=None):
        observed.update(config)
        _emit_upstream_progress(progress_callback)
        _write_upstream_outputs(config)

    _mock_upstream(monkeypatch, fake_pops_main)
    run_pops_direct(args)
    assert observed["y_error_cov_path"].endswith(".npz")


def test_failed_upstream_run_never_publishes_partial_results(tmp_path, monkeypatch):
    resources = _pops_resources(tmp_path)
    args = _arguments(tmp_path, resources)

    def incomplete_pops_main(config, progress_callback=None):
        _emit_upstream_progress(progress_callback)
        Path(config["out_prefix"] + ".preds").write_text("partial\n", encoding="utf-8")

    _mock_upstream(monkeypatch, incomplete_pops_main)
    with pytest.raises(PopsError, match="did not create required"):
        run_pops_direct(args)
    output = Path(args.output_directory)
    assert not (output / "STUDY_pops.preds").exists()
    assert not (output / "run_metadata" / "STUDY_pops_completion.yaml").exists()
