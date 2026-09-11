"""Regression tests for the PostGWAS boundary around upstream PoPS v0.2."""

from __future__ import annotations

import argparse
import logging
from argparse import Namespace
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml
from rich.cells import cell_len
from rich.console import Console

from postgwas.config import load_configuration, load_module_configuration
from postgwas.core.preflight import PipelinePreflightEvidence
from postgwas.core.input_validation import InputValidationSession
from postgwas.core.ui import PipelineStageController
from postgwas.core.validation_reporting import FileValidationDisplay
from postgwas.modules.magma.reporting import (
    MAGMA_GENE_ONLY_STAGE_KEYS,
    MAGMA_STAGE_TITLES,
)
from postgwas.modules.pops.cli import build_parser
from postgwas.modules.pops.errors import PopsError
from postgwas.modules.pops.service import (
    _resolved_configuration,
    _require_pops_runtime,
    _validate_feature_resources,
    preflight_pops_pipeline,
    run_pops_direct,
    validate_pops_configuration,
)
from postgwas.modules.pops.stages import POPS_STAGES, pops_pipeline_progress_plan
from postgwas.pipeline.planner import build_pipeline_plan
from postgwas.pipeline.registry import REGISTRY
from preflight_support import pipeline_input_vcf_evidence
from postgwas.pipeline.runners import run_pops_runner


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


def test_pops_feature_companions_have_one_validated_bundle_card(
    tmp_path, capsys,
):
    resources = _pops_resources(tmp_path)
    module = load_module_configuration("pops").model_copy(update={
        "feature_matrix_prefix": str(resources["feature_prefix"]),
        "feature_matrix_chunks": 1,
    })
    with InputValidationSession() as session:
        display = FileValidationDisplay(session, load_configuration())
        features = _validate_feature_resources(module)
        display.flush()

    screen = capsys.readouterr().out
    assert "PoPS feature-matrix resources — CHECKS PASSED" in screen
    assert "Matrix validation" in screen
    assert "features.cols.0.txt" not in screen
    assert "features.mat.0.npy" not in screen
    assert features["feature_count"] == 2


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
    feature_ids = np.atleast_1d(np.loadtxt(
        config["feature_mat_prefix"] + ".rows.txt", dtype=str,
    )).reshape(-1).tolist()
    if config["magma_prefix"] is not None:
        target = pd.read_csv(config["magma_prefix"] + ".genes.out", sep=r"\s+")
        target_scores = dict(zip(target["GENE"].astype(str), target["ZSTAT"]))
        has_covariates = config["use_magma_covariates"]
    else:
        target = pd.read_csv(config["y_path"], sep="\t")
        target_scores = dict(zip(target["ENSGID"].astype(str), target["Score"]))
        has_covariates = config["y_covariates_path"] is not None
    prediction = pd.DataFrame({
        "ENSGID": feature_ids,
        "PoPS_Score": [
            2.0 if gene == "ENSG002" else 1.0 if gene == "ENSG001" else -index
            for index, gene in enumerate(feature_ids, 1)
        ],
        "Y": [target_scores.get(gene) for gene in feature_ids],
        "feature_selection_gene": [gene in target_scores for gene in feature_ids],
        "training_gene": [gene in target_scores for gene in feature_ids],
    })
    if has_covariates:
        prediction["Y_proj"] = prediction["Y"]
        prediction["project_out_covariates_gene"] = prediction["Y"].notna()
    prediction.to_csv(prefix + ".preds", sep="\t", index=False)
    Path(prefix + ".coefs").write_text("parameter\tbeta\nfeature_a\t1\n", encoding="utf-8")
    Path(prefix + ".marginals").write_text(
        "feature\tpval\tselected\nfeature_a\t0.01\tTrue\n",
        encoding="utf-8",
    )
    Path(prefix + ".log").write_text("completed\n", encoding="utf-8")


def _write_annotated_magma_results(
    path: Path,
    genes: list[str],
    scores: dict[str, float],
    chromosomes: dict[str, str],
) -> None:
    rows = []
    for index, gene in enumerate(genes, 1):
        pvalue = min(1.0, 0.001 * index)
        rows.append({
            "GENE": gene,
            "CHR": chromosomes[gene],
            "START": 1000 + index,
            "STOP": 1100 + index,
            "NSNPS": 2,
            "NPARAM": 1,
            "N": 100,
            "ZSTAT": scores[gene],
            "P": pvalue,
            "GENE_REFERENCE_CHR": chromosomes[gene],
            "GENE_REFERENCE_START": 1000 + index,
            "GENE_REFERENCE_END": 1100 + index,
            "GENE_REFERENCE_STRAND": "+",
            "GENE_SYMBOL": "SYMBOL_%03d" % index,
            "P_bonferroni_corr": min(1.0, pvalue * len(genes)),
            "P_fdr_bh_corr": min(1.0, pvalue * 2),
        })
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


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
        "genes_scored": 10,
        "model_features": 1,
        "outputs_written": 3,
    })


def test_cli_has_no_independent_defaults():
    parser = build_parser()
    destinations = {
        "magma_association_prefix", "magma_annotated_results_file",
        "feature_matrix_prefix",
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
    assert (
        load_configuration().modules.pops.output_layout.integrated_results_suffix
        == ".integrated_gene_results.tsv"
    )
    annotated_without_magma = tmp_path / "annotated_without_magma.yaml"
    annotated_without_magma.write_text(
        "magma_annotated_results_file: magma_annotated.tsv\n",
        encoding="utf-8",
    )
    with pytest.raises(Exception, match="requires magma_association_prefix"):
        load_module_configuration("pops", annotated_without_magma)
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


def test_pipeline_plan_preserves_magma_validation_before_pops_stages():
    plan = pops_pipeline_progress_plan(argparse.Namespace())
    execution = build_pipeline_plan(["pops"])

    assert plan is not None
    assert plan["kind"] == "pops"
    assert plan["label"] == "PoPS pipeline execution progress"
    upstream_stages = tuple(
        MAGMA_STAGE_TITLES[key] for key in MAGMA_GENE_ONLY_STAGE_KEYS
    )
    assert plan["stages"] == (*upstream_stages, *POPS_STAGES)
    assert len(plan["stages"]) == 19
    assert all("pathway" not in stage.lower() for stage in plan["stages"])
    assert plan["modules"] == {
        "formatter": (1, 4),
        "magma": (4, 8),
        "pops": (9, 19),
    }
    assert plan["stage_numbers"] == {
        "target_inputs": 9,
        "gene_annotation": 10,
        "feature_matrix": 11,
        "feature_controls": 12,
        "gene_compatibility": 13,
        "target_loading": 14,
        "covariate_adjustment": 15,
        "feature_selection": 16,
        "model_fitting": 17,
        "gene_scoring": 18,
        "results": 19,
    }
    assert REGISTRY.get("pops").pipeline_progress_factory == (
        "postgwas.modules.pops.stages:pops_pipeline_progress_plan"
    )
    assert execution.steps == ("formatter", "magma", "pops")

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
    preflight_pops_pipeline(
        args, preflight_evidence=pipeline_input_vcf_evidence(),
    )
    Path(str(resources["feature_prefix"]) + ".mat.0.npy").unlink()
    with pytest.raises(PopsError, match="feature matrix chunk"):
        preflight_pops_pipeline(
            args, preflight_evidence=pipeline_input_vcf_evidence(),
        )


def test_feature_matrix_validation_rejects_nonfinite_values(tmp_path, capsys):
    resources = _pops_resources(tmp_path)
    matrix_path = Path(str(resources["feature_prefix"]) + ".mat.0.npy")
    matrix = np.load(matrix_path)
    matrix[0, 0] = np.nan
    np.save(matrix_path, matrix)

    with pytest.raises(PopsError, match="contains non-finite values"):
        run_pops_direct(_arguments(tmp_path, resources))
    compact_screen = " ".join(capsys.readouterr().out.split())
    assert "Failed 3/11 · Validate the PoPS feature-matrix files" in compact_screen
    assert "All 11 stages completed" not in compact_screen
    canonical_log = (
        tmp_path / "output" / "logs" / "STUDY_pops_service.log"
    ).read_text(encoding="utf-8")
    assert "Feature matrix chunk 0 contains non-finite values" in canonical_log


def test_feature_controls_cannot_be_silently_removed_by_feature_subset(tmp_path):
    resources = _pops_resources(tmp_path)
    subset = tmp_path / "subset.txt"
    subset.write_text("feature_a\n", encoding="utf-8")
    controls = tmp_path / "controls.txt"
    controls.write_text("feature_b\n", encoding="utf-8")

    with pytest.raises(PopsError, match="silently discarded"):
        validate_pops_configuration(_arguments(
            tmp_path,
            resources,
            feature_subset_file=str(subset),
            control_features_file=str(controls),
        ))


def test_pipeline_preflight_defers_vcf_build_comparison_until_execution(tmp_path):
    resources = _pops_resources(tmp_path)
    args = _arguments(tmp_path, resources)
    del args.magma_association_prefix
    del args.genome_build
    args.pops_genome_build = "GRCh38"

    validated = preflight_pops_pipeline(
        args, preflight_evidence=pipeline_input_vcf_evidence(),
    )

    assert isinstance(validated, PipelinePreflightEvidence)
    assert validated.module == "pops"
    configuration, features, annotation, controls = validated.resources
    assert configuration.modules.pops.genome_build == "GRCh38"
    assert features["row_count"] == annotation["gene_count"] == 10
    assert controls is features["control_resources"]
    assert validated.deferred_checks == (
        "Validate the pipeline-generated MAGMA gene-association result.",
    )
    assert not hasattr(args, "_pops_resource_preflight")


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
    assert "MAGMA gene Z-scores=11" in message
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
    assert "Original genes with input scores" in screen
    assert screen.index("Original genes with input scores") < screen.index(
        "Started 6/11"
    )
    assert "Original MAGMA genes with Z-scores" in screen
    assert "Excluded MAGMA genes with Z-scores" in screen
    aligned_labels = (
        "Dataset", "Genes receiving finite PoPS scores",
        "Original MAGMA genes with Z-scores", "1. ENSG002",
        "PoPS significance cutoff", "Official PoPS predictions",
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


def test_integrated_results_merge_full_gene_union_and_annotated_magma(
    tmp_path, monkeypatch,
):
    resources = _pops_resources(tmp_path)
    genes = ["ENSG%03d" % index for index in range(1, 10)] + ["ENSG999"]
    scores = {gene: index / 10 for index, gene in enumerate(genes)}
    chromosomes = {
        gene: "1" if gene in {
            "ENSG001", "ENSG002", "ENSG003", "ENSG004", "ENSG005",
        }
        else "2"
        for gene in genes
    }
    prefix = resources["magma_prefix"]
    Path(str(prefix) + ".genes.out").write_text(
        "GENE\tZSTAT\n" + "".join(
            "%s\t%.6f\n" % (gene, scores[gene]) for gene in genes
        ),
        encoding="utf-8",
    )
    Path(str(prefix) + ".genes.raw").write_text(
        "# VERSION = test\n# COVAR = NSAMP MAC\n" + "".join(
            "%s %s 1 2 1 1 100 1 0\n" % (gene, chromosomes[gene])
            for gene in genes
        ),
        encoding="utf-8",
    )
    annotated = tmp_path / "magma_genes_annotated.tsv"
    _write_annotated_magma_results(annotated, genes, scores, chromosomes)
    source_text = annotated.read_text(encoding="utf-8")

    def fake_pops_main(config, progress_callback=None):
        _write_upstream_outputs(config, progress_callback=progress_callback)
        predictions = pd.read_csv(config["out_prefix"] + ".preds", sep="\t")
        predictions.loc[
            predictions["ENSGID"].eq("ENSG009"), "training_gene"
        ] = False
        predictions.to_csv(config["out_prefix"] + ".preds", sep="\t", index=False)

    _mock_upstream(monkeypatch, fake_pops_main)
    run_config = tmp_path / "pops_integrated.yaml"
    run_config.write_text("minimum_gene_count: 2\n", encoding="utf-8")
    result = run_pops_direct(_arguments(
        tmp_path,
        resources,
        gene_universe_policy="intersect",
        magma_annotated_results_file=str(annotated),
        run_config=str(run_config),
    ))

    integrated_path = Path(result["integrated_results_file"])
    html_path = Path(result["integrated_report_file"])
    integrated_table = pd.read_csv(
        integrated_path, sep="\t", dtype=str, keep_default_na=False,
    )
    assert integrated_table.columns[:17].tolist() == [
        "gene_id",
        "gene_symbol",
        "pops_chromosome",
        "pops_tss",
        "pops_score",
        "pops_rank",
        "gene_analysis_status",
        "magma_zstat",
        "magma_pvalue",
        "magma_bonferroni_pvalue",
        "magma_fdr_pvalue",
        "magma_chromosome",
        "magma_start",
        "magma_end",
        "magma_snp_count",
        "magma_parameter_count",
        "magma_sample_size",
    ]
    integrated = integrated_table.set_index("gene_id")
    assert len(integrated) == 11
    assert integrated.index[-1] == "ENSG999"
    assert integrated.loc["ENSG002", "pops_rank"] == "1"
    assert integrated.loc["ENSG001", "gene_symbol"] == "SYMBOL_001"
    assert integrated.loc["ENSG001", "magma_pvalue"] == "0.001"
    assert integrated.loc["ENSG001", "gene_analysis_status"] == (
        "scored_and_used_for_model_fitting"
    )
    assert integrated.loc["ENSG009", "gene_analysis_status"] == (
        "scored_with_target_not_used_for_model_fitting"
    )
    assert integrated.loc["ENSG010", "gene_analysis_status"] == (
        "scored_without_input_target_score"
    )
    assert integrated.loc["ENSG010", "input_target_score"] == "NA"
    assert integrated.loc["ENSG999", "gene_analysis_status"] == (
        "input_target_excluded_from_pops"
    )
    assert integrated.loc["ENSG999", "pops_score"] == "NA"
    assert integrated.loc["ENSG999", "magma_zstat"] == "0.9"
    assert result["summary"]["integrated_gene_count"] == 11
    assert result["summary"]["genes_scored_without_target"] == 1
    assert result["summary"]["genes_scored_with_target_not_fitted"] == 1
    assert result["summary"]["input_target_genes_excluded"] == 1
    html = html_path.read_text(encoding="utf-8")
    assert "Integrated PoPS gene results" in html
    assert "Download TSV" in html
    assert "scored_without_input_target_score" in html
    assert "SYMBOL_001" in html
    assert annotated.read_text(encoding="utf-8") == source_text


def test_annotated_magma_mismatch_fails_before_upstream_execution(
    tmp_path, monkeypatch,
):
    resources = _pops_resources(tmp_path)
    genes = ["ENSG%03d" % index for index in range(1, 11)]
    scores = {gene: index / 10 for index, gene in enumerate(genes)}
    chromosomes = {
        gene: "1" if index < 5 else "2"
        for index, gene in enumerate(genes)
    }
    annotated = tmp_path / "magma_genes_annotated.tsv"
    _write_annotated_magma_results(annotated, genes, scores, chromosomes)
    table = pd.read_csv(annotated, sep="\t")
    table.loc[0, "ZSTAT"] = 99
    table.to_csv(annotated, sep="\t", index=False)
    called = False

    def fake_pops_main(config, progress_callback=None):
        nonlocal called
        called = True

    _mock_upstream(monkeypatch, fake_pops_main)
    with pytest.raises(PopsError, match="do not match"):
        run_pops_direct(_arguments(
            tmp_path,
            resources,
            magma_annotated_results_file=str(annotated),
        ))
    assert not called
    assert not (tmp_path / "output" / "STUDY_pops.integrated_gene_results.tsv").exists()


def test_pipeline_passes_validated_annotated_magma_artifact_to_pops(
    tmp_path, monkeypatch,
):
    annotated = tmp_path / "upstream_magma_annotated.tsv"
    prefix = tmp_path / "upstream_magma"
    observed = {}

    def fake_run_pops_direct(args, ctx):
        observed["prefix"] = args.magma_association_prefix
        observed["annotated"] = args.magma_annotated_results_file
        observed["output_directory"] = args.output_directory
        return {"status": "success"}

    monkeypatch.setattr(
        "postgwas.modules.pops.service.run_pops_direct", fake_run_pops_direct,
    )
    args = Namespace(
        output_directory=str(tmp_path / "pipeline"),
        _step_num=3,
    )
    ctx = {"magma": {
        "primary_mapping": "positional",
        "mapping_analyses": {
            "positional": {"result_statistic_type": "calibrated_gene_p_value"},
        },
        "magma_genes_prefix": str(prefix),
        "magma_genes_annotated": str(annotated),
    }}

    assert run_pops_runner(args, ctx) == {"status": "success"}
    assert observed["prefix"] == str(prefix)
    assert observed["annotated"] == str(annotated)
    assert observed["output_directory"].endswith("03_pops")
    assert args.output_directory == str(tmp_path / "pipeline")


def test_preflight_rejects_magma_raw_gene_order_mismatch(tmp_path):
    resources = _pops_resources(tmp_path)
    raw_path = Path(str(resources["magma_prefix"]) + ".genes.raw")
    lines = raw_path.read_text(encoding="utf-8").splitlines()
    lines[2], lines[3] = lines[3], lines[2]
    raw_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(PopsError, match=r"gene order differs at data row 1"):
        run_pops_direct(_arguments(tmp_path, resources))


def test_preflight_rejects_invalid_magma_technical_covariate_metadata(tmp_path):
    resources = _pops_resources(tmp_path)
    raw_path = Path(str(resources["magma_prefix"]) + ".genes.raw")
    raw_path.write_text(
        raw_path.read_text(encoding="utf-8").replace(
            "ENSG001 1 1 2 1 1 100 1 0",
            "ENSG001 1 1 2 0 1 100 1 0",
            1,
        ),
        encoding="utf-8",
    )

    with pytest.raises(PopsError, match="non-positive or non-finite NSNPS"):
        validate_pops_configuration(_arguments(tmp_path, resources))


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
    assert result["summary"]["genes_scored"] == 10
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
    compact_screen = " ".join(screen.split())
    assert "PoPS analysis progress" in screen
    assert "Completed 11/11 · Validate and publish PoPS results" in screen
    assert "MAGMA gene Z-score inputs" in screen
    assert "Scores used to fit PoPS" in screen
    assert "MAGMA gene-level association Z-scores (ZSTAT)" in compact_screen
    assert "PoPS gene-location annotation" in screen
    assert "PoPS feature matrices" in screen
    assert "PoPS feature controls" in screen
    assert "PoPS gene-identifier compatibility" in screen
    assert "MAGMA gene Z-scores loaded" in screen
    assert "Features selected" in screen
    assert "Genes used to fit the model" in screen
    assert "PoPS gene-prioritisation summary" in screen
    assert "Genes receiving finite PoPS scores" in screen
    assert "None. PoPS scores are relative rankings, not p-values." in screen
    assert "COMPLETED WITH SCIENTIFIC WARNINGS" in screen
    assert "MAGMA gene Z-scores cover 2 of 10 compatible genes" in screen
    assert result["summary"]["warnings"]


def test_real_upstream_reports_ordered_stages_without_internal_log_noise(
    tmp_path, capsys,
):
    resources = _pops_resources(tmp_path, custom_target=True)
    args = _arguments(tmp_path, resources)

    result = run_pops_direct(args)

    screen = capsys.readouterr().out
    completed = [
        screen.index("Completed %d/11" % step)
        for step in range(1, 12)
    ]
    assert completed == sorted(completed)
    assert "Scores used to fit PoPS" in screen
    assert "Custom gene scores loaded" in screen
    assert "MAGMA gene Z-scores loaded" not in screen
    assert "Selection strategy" in screen
    assert "Prediction model" in screen
    assert "Config dict =" not in screen
    assert "Computing marginal association table" not in screen
    assert Path(result["pops_file"]).is_file()
    upstream_log = next(
        Path(path) for path in result["published_files"] if path.endswith(".log")
    )
    assert "Config dict =" in upstream_log.read_text(encoding="utf-8")


def test_detailed_pipeline_progress_continues_from_magma_through_pops_validation(
    tmp_path, monkeypatch,
):
    resources = _pops_resources(tmp_path)
    args = _arguments(tmp_path, resources)
    plan = pops_pipeline_progress_plan(args)
    assert plan is not None
    stream = StringIO()
    controller = PipelineStageController(
        plan["label"],
        plan["stages"],
        console=Console(
            file=stream,
            force_terminal=True,
            color_system=None,
            width=120,
        ),
        outcome_label_width=38,
    )
    for number in range(1, len(MAGMA_GENE_ONLY_STAGE_KEYS) + 1):
        controller.start(number)
        controller.complete(number, outcome="validated upstream MAGMA stage")
    args._pipeline_progress_plan = plan
    args._pipeline_stage_progress = controller
    _mock_upstream(monkeypatch, _write_upstream_outputs)

    result = run_pops_direct(args)

    assert Path(result["pops_file"]).is_file()
    assert controller.completed == 18
    assert controller.current == 19
    completion = args._pipeline_stage_completion
    controller.complete(
        19,
        outcome=completion.get("outcome"),
        outcome_fields=completion.get("outcome_fields"),
    )
    controller.close()
    compact_screen = " ".join(stream.getvalue().split())
    assert "PoPS pipeline execution progress" in compact_screen
    assert "Completed 19/19 · Validate and publish PoPS results" in compact_screen
    assert "All 19 stages completed" in compact_screen
    assert "pathway" not in compact_screen.lower()
    assert "MAGMA gene-statistic file" in compact_screen
    assert "Transcription-start-site column" in compact_screen
    assert "Matrix chunks" in compact_screen
    assert "Control-features file" in compact_screen
    assert "Genes with scores shared by all inputs" in compact_screen
    canonical_log = (
        Path(args.output_directory) / "logs" / "STUDY_pops_service.log"
    ).read_text(encoding="utf-8")
    assert canonical_log.count("stage_outcome") == len(POPS_STAGES)
    assert "target_score_file=" in canonical_log
    assert "gene_annotation_file=" in canonical_log
    assert "pops_feature_matrix_chunk" in canonical_log
    assert "matrix_feature_compatibility=passed" in canonical_log


def test_detailed_pipeline_rejects_pops_build_that_differs_from_validated_vcf(
    tmp_path,
):
    resources = _pops_resources(tmp_path)
    args = _arguments(tmp_path, resources, pops_genome_build="GRCh38")
    del args.genome_build
    plan = pops_pipeline_progress_plan(args)
    assert plan is not None
    stream = StringIO()
    controller = PipelineStageController(
        plan["label"],
        plan["stages"],
        console=Console(
            file=stream,
            force_terminal=True,
            color_system=None,
            width=120,
        ),
    )
    for number in range(1, len(MAGMA_GENE_ONLY_STAGE_KEYS) + 1):
        controller.start(number)
        controller.complete(number)
    args._pipeline_progress_plan = plan
    args._pipeline_stage_progress = controller
    args._pipeline_vcf_header_evidence = {"genome_build": "GRCh37"}

    with pytest.raises(PopsError, match="does not match pipeline MAGMA genome build"):
        run_pops_direct(args)
    controller.fail_active()
    controller.close()

    compact_screen = " ".join(stream.getvalue().split())
    assert (
        "Failed 13/19 · Compare input-score, annotation, and feature gene identifiers"
        in compact_screen
    )
    assert "All 19 stages completed" not in compact_screen


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
