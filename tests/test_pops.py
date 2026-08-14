"""Regression tests for the PostGWAS boundary around upstream PoPS v0.2."""

from __future__ import annotations

import argparse
import logging
from argparse import Namespace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml
from rich.cells import cell_len

from postgwas.config import load_configuration, load_module_configuration
from postgwas.modules.pops.cli import build_parser
from postgwas.modules.pops.errors import PopsError
from postgwas.modules.pops.service import (
    _resolved_configuration,
    _require_pops_runtime,
    preflight_pops_pipeline,
    run_pops_direct,
    validate_pops_configuration,
)
from postgwas.pipeline.registry import REGISTRY


def _pops_resources(tmp_path: Path, *, custom_target: bool = False) -> dict[str, Path]:
    genes = ["ENSG%03d" % index for index in range(1, 11)]
    annotation = tmp_path / "genes.tsv"
    annotation.write_text(
        "ENSGID\tCHR\tTSS\n"
        + "".join(
            "%s\t%s\t%d\n" % (gene, 1 if index < 5 else 2, 1000 + index)
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
            "# VERSION = test\n# COVAR = NSAMP MAC\n"
            + "".join(
                "%s %s 1 2 1 1 100 1 0\n"
                % (gene, 1 if index < 5 else 2)
                for index, gene in enumerate(genes)
            ),
            encoding="utf-8",
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
        "gene_universe_policy",
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


def test_missing_required_inputs_are_reported_together_before_runtime(
    tmp_path, monkeypatch,
):
    monkeypatch.setattr(
        "postgwas.modules.pops.service._require_pops_runtime",
        lambda: pytest.fail("runtime validation must follow required arguments"),
    )
    args = Namespace(
        dataset_id="STUDY",
        output_directory=str(tmp_path / "output"),
    )

    with pytest.raises(Exception) as captured:
        run_pops_direct(args)

    message = str(captured.value)
    assert "Required argument not provided: --genome-build." in message
    assert "Required argument not provided: --feature-matrix-prefix." in message
    assert "Required argument not provided: --pops-gene-location-file." in message
    assert "--magma-association-prefix or --target-score-file" in message
    service_log = tmp_path / "output" / "logs" / "STUDY_pops_service.log"
    assert "Required argument not provided" in service_log.read_text(encoding="utf-8")
    assert not (tmp_path / "output" / "STUDY_pops.preds").exists()


def test_configuration_uses_upstream_method_names_and_validates_ranges(tmp_path):
    assert load_configuration().modules.pops.method == "ridge"
    assert load_configuration().modules.pops.gene_universe_policy == "strict"
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

    unsupported_custom_intersection = tmp_path / "custom_intersection.yaml"
    unsupported_custom_intersection.write_text(
        "gene_universe_policy: intersect\n"
        "target_score_file: targets.tsv\n",
        encoding="utf-8",
    )
    with pytest.raises(Exception, match="supported only for MAGMA targets"):
        load_module_configuration("pops", unsupported_custom_intersection)


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


def test_preflight_rejects_magma_gene_universe_mismatch_with_diagnosis(tmp_path):
    resources = _pops_resources(tmp_path)
    output_path = Path(str(resources["magma_prefix"]) + ".genes.out")
    raw_path = Path(str(resources["magma_prefix"]) + ".genes.raw")
    output_path.write_text(
        output_path.read_text(encoding="utf-8") + "ENSG999\t1.0\n",
        encoding="utf-8",
    )
    raw_path.write_text(
        raw_path.read_text(encoding="utf-8")
        + "ENSG999 2 1 2 1 1 100 1 0\n",
        encoding="utf-8",
    )

    with pytest.raises(PopsError) as captured:
        run_pops_direct(_arguments(tmp_path, resources))

    message = str(captured.value)
    assert "PoPS input files are scientifically incompatible" in message
    assert "MAGMA genes=11" in message
    assert "PoPS annotation genes=10" in message
    assert "PoPS feature-row genes=10" in message
    assert "shared genes=10" in message
    assert "ENSG999" in message
    assert str(output_path) in message
    assert not (Path(tmp_path) / "output" / "STUDY_pops.preds").exists()


def test_intersection_policy_preserves_covariance_and_publishes_detailed_audit(
    tmp_path, monkeypatch, capsys,
):
    from postgwas.modules.pops.pops import munge_magma_covariance_metadata

    resources = _pops_resources(tmp_path)
    original_genes = [
        "ENSG001", "ENSG002", "ENSG003", "ENSG004", "ENSG999",
        "ENSG005", "ENSG006", "ENSG007", "ENSG008", "ENSG009", "ENSG010",
    ]
    magma_prefix = resources["magma_prefix"]
    Path(str(magma_prefix) + ".genes.out").write_text(
        "GENE\tZSTAT\n"
        + "".join(
            "%s\t%.6f\n" % (gene, index / 10)
            for index, gene in enumerate(original_genes)
        ),
        encoding="utf-8",
    )
    raw_lines = ["# VERSION = test", "# COVAR = NSAMP MAC"]
    chromosome_row = {"1": 0, "2": 0}
    for gene in original_genes:
        chromosome = "1" if gene in {
            "ENSG001", "ENSG002", "ENSG003", "ENSG004", "ENSG005", "ENSG999",
        } else "2"
        row_index = chromosome_row[chromosome]
        correlations = [
            0.001 * (row_index + column_index + 1)
            for column_index in range(row_index)
        ]
        raw_lines.append(" ".join(
            [gene, chromosome, "1", "2", "2", "2", "100", "10", "0"]
            + ["%.17g" % value for value in correlations]
        ))
        chromosome_row[chromosome] += 1
    original_raw = Path(str(magma_prefix) + ".genes.raw")
    original_raw.write_text("\n".join(raw_lines) + "\n", encoding="utf-8")
    original_raw_text = original_raw.read_text(encoding="utf-8")
    observed = {}

    def fake_pops_main(config, progress_callback=None):
        derived_prefix = Path(config["magma_prefix"])
        observed["derived_prefix"] = str(derived_prefix)
        observed["compatible_ids"] = pd.read_csv(
            str(derived_prefix) + ".genes.out", sep=r"\s+",
        )["GENE"].tolist()
        source_sigmas, _ = munge_magma_covariance_metadata(str(original_raw))
        derived_sigmas, _ = munge_magma_covariance_metadata(
            str(derived_prefix) + ".genes.raw"
        )
        retained_by_block = ([0, 1, 2, 3, 5], [0, 1, 2, 3, 4])
        observed["expected_covariance"] = [
            sigma[np.ix_(indices, indices)]
            for sigma, indices in zip(source_sigmas, retained_by_block)
        ]
        observed["derived_covariance"] = derived_sigmas
        _write_upstream_outputs(config, progress_callback=progress_callback)

    _mock_upstream(monkeypatch, fake_pops_main)
    result = run_pops_direct(_arguments(
        tmp_path, resources, gene_universe_policy="intersect",
    ))

    assert observed["compatible_ids"] == [
        gene for gene in original_genes if gene != "ENSG999"
    ]
    for derived, expected in zip(
        observed["derived_covariance"], observed["expected_covariance"],
    ):
        np.testing.assert_allclose(derived, expected)
    assert original_raw.read_text(encoding="utf-8") == original_raw_text
    published = {Path(path).name: Path(path) for path in result["published_files"]}
    compatible_out = published["STUDY_pops.compatible.genes.out"]
    compatible_raw = published["STUDY_pops.compatible.genes.raw"]
    excluded_out = published["STUDY_pops.excluded.genes.out"]
    excluded_raw = published["STUDY_pops.excluded.genes.raw"]
    audit_path = published["STUDY_pops.gene_compatibility.tsv"]
    report_path = published["STUDY_pops.gene_compatibility.yaml"]
    assert all(path.is_file() for path in (
        compatible_out, compatible_raw, excluded_out, excluded_raw,
        audit_path, report_path,
    ))
    assert pd.read_csv(excluded_out, sep=r"\s+")["GENE"].tolist() == ["ENSG999"]
    assert [
        line.split()[0]
        for line in excluded_raw.read_text(encoding="utf-8").splitlines()[2:]
    ] == ["ENSG999"]
    audit = pd.read_csv(audit_path, sep="\t")
    excluded_audit = audit.loc[audit["gene_id"] == "ENSG999"].iloc[0]
    assert len(audit) == 11
    assert not bool(excluded_audit["retained_for_pops"])
    assert excluded_audit["decision"] == "missing_annotation_and_features"
    report = yaml.safe_load(report_path.read_text(encoding="utf-8"))
    assert report["original_target_genes"] == 11
    assert report["retained_target_genes"] == 10
    assert report["excluded_target_genes"] == 1
    assert report["absent_from_both"] == 1
    assert report["chromosomes"]["1"] == {
        "original": 6, "retained": 5, "excluded": 1,
    }
    assert report["chromosomes"]["2"] == {
        "original": 5, "retained": 5, "excluded": 0,
    }
    assert report["exclusion_reason_counts"] == {
        "retained": 10,
        "missing_annotation": 0,
        "missing_features": 0,
        "missing_annotation_and_features": 1,
    }
    assert report["original_files_unchanged"] is True
    assert result["summary"]["gene_compatibility"]["excluded_target_genes"] == 1
    screen = capsys.readouterr().out
    assert "MAGMA–PoPS gene compatibility" in screen
    assert "Original target genes" in screen
    assert screen.index("Original target genes") < screen.index("Started 2/7")
    assert "Original MAGMA target genes" in screen
    assert "Excluded target genes" in screen
    aligned_labels = (
        "Dataset", "Genes assigned PoPS scores", "Original MAGMA target genes",
        "1. ENSG002", "PoPS significance cutoff", "Complete PoPS results",
    )
    summary_screen = screen[screen.index("PoPS gene-prioritisation summary"):]
    colon_columns = {
        cell_len(line.split(":", 1)[0])
        for line in summary_screen.splitlines()
        if any(label in line for label in aligned_labels) and ":" in line
    }
    assert colon_columns == {57}
    log_text = (Path(tmp_path) / "output" / "logs" / "STUDY_pops_service.log").read_text(
        encoding="utf-8"
    )
    assert "pops_gene_universe_intersection" in log_text
    assert "pops_gene_compatibility_chromosome" in log_text
    assert "excluded_target_genes=1" in log_text


def test_preflight_rejects_magma_raw_gene_order_mismatch(tmp_path):
    resources = _pops_resources(tmp_path)
    raw_path = Path(str(resources["magma_prefix"]) + ".genes.raw")
    lines = raw_path.read_text(encoding="utf-8").splitlines()
    lines[2], lines[3] = lines[3], lines[2]
    raw_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(PopsError, match=r"gene order differs at data row 1"):
        run_pops_direct(_arguments(tmp_path, resources))


def test_preflight_rejects_magma_annotation_chromosome_mismatch(tmp_path):
    resources = _pops_resources(tmp_path)
    annotation_path = resources["annotation"]
    annotation_path.write_text(
        annotation_path.read_text(encoding="utf-8").replace(
            "ENSG001\t1\t1000", "ENSG001\t2\t1000",
        ),
        encoding="utf-8",
    )

    with pytest.raises(PopsError) as captured:
        run_pops_direct(_arguments(
            tmp_path, resources, gene_universe_policy="intersect",
        ))

    message = str(captured.value)
    assert "disagree on chromosome for 1 shared genes" in message
    assert "ENSG001 (MAGMA=1, PoPS=2)" in message
    assert "intersection cannot repair chromosome disagreement" in message


def test_preflight_rejects_annotation_feature_gene_universe_mismatch(tmp_path):
    resources = _pops_resources(tmp_path, custom_target=True)
    rows_path = Path(str(resources["feature_prefix"]) + ".rows.txt")
    rows = rows_path.read_text(encoding="utf-8").replace("ENSG010\n", "ENSG999\n")
    rows_path.write_text(rows, encoding="utf-8")

    with pytest.raises(PopsError) as captured:
        run_pops_direct(_arguments(tmp_path, resources))

    message = str(captured.value)
    assert "feature rows contain genes absent from the gene annotation" in message
    assert "absent from annotation=1" in message
    assert "ENSG999" in message


def test_preflight_allows_annotation_genes_without_feature_rows(tmp_path):
    resources = _pops_resources(tmp_path)
    annotation_path = resources["annotation"]
    annotation_path.write_text(
        annotation_path.read_text(encoding="utf-8") + "ENSG999\t1\t9999\n",
        encoding="utf-8",
    )

    _, features, annotation, outcome = validate_pops_configuration(
        _arguments(tmp_path, resources)
    )

    assert features["row_count"] == 10
    assert annotation["gene_count"] == 11
    assert outcome["gene_count"] == 10


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
    Path(str(resources["magma_prefix"]) + ".genes.raw").write_text(
        "# VERSION = test\n# COVAR = NSAMP MAC\n"
        "ENSG001 1 1 2 1 1 100 1 0\n"
        "ENSG002 1 1 2 1 1 100 1 0\n",
        encoding="utf-8",
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
