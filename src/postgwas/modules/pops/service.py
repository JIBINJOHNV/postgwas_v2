"""Validated PostGWAS service boundary for the upstream PoPS v0.2 program."""

from __future__ import annotations

import argparse
from importlib.util import find_spec
from itertools import groupby
import logging
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd
from scipy import sparse

from postgwas.config import (
    load_configuration,
    load_run_configuration_for_module,
    write_resolved_configuration,
)
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.core.matrix_validation import read_unique_names, validate_feature_matrix_bundle
from postgwas.core.paths import (
    configured_output_path,
    require_nonempty_file,
    validate_filename_component,
)
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.core.preflight import (
    PipelinePreflightEvidence,
    pipeline_preflight_evidence,
    require_pipeline_input_vcf,
)
from postgwas.core.required_arguments import (
    RequiredAlternative,
    RequiredArgument,
    require_resolved_arguments,
)
from postgwas.core.gene_annotation_validation import validate_gene_tss_annotation
from postgwas.core.io.tables import read_pandas_table, require_table_columns
from postgwas.core.io.reports import write_delimited_report, write_yaml_report
from postgwas.core.ui.progress import StageProgress
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.core.validation_reporting import register_file_validation_bundle
from postgwas.core.values import format_percentage
from postgwas.modules.pops.errors import PopsError
from postgwas.modules.pops.reporting import write_integrated_gene_report
from postgwas.modules.pops.stages import (
    POPS_STAGES,
    pops_pipeline_stage_numbers,
)


_POPS_UPSTREAM_PROGRESS_ORDER = (
    "target_loading",
    "covariate_adjustment",
    "feature_selection",
    "model_fitting",
    "gene_scoring",
)


def _score_wording(source: str) -> dict[str, str]:
    """Return source-specific terms without changing PoPS's generic interface."""
    if source == "MAGMA":
        return {
            "heading": "MAGMA gene Z-score inputs",
            "source": "MAGMA gene-level association Z-scores (ZSTAT)",
            "plural": "MAGMA gene Z-scores",
            "coverage": "MAGMA Z-scores",
            "column": "MAGMA gene Z-score column",
            "loaded": "MAGMA gene Z-scores loaded",
        }
    return {
        "heading": "Custom gene-score inputs",
        "source": "custom gene scores",
        "plural": "custom gene scores",
        "coverage": "custom gene scores",
        "column": "Custom gene-score column",
        "loaded": "Custom gene scores loaded",
    }


def _target_validation_fields(module, outcome: dict) -> list[tuple[str, str, object]]:
    wording = _score_wording(outcome["label"])
    if outcome["label"] == "MAGMA":
        fields = [
            ("analysis", wording["heading"]),
            ("genetic", "Scores used to fit PoPS", wording["source"]),
            ("info", "MAGMA gene-statistic file", outcome["output_path"].name),
            ("info", "MAGMA covariance file", outcome["raw_path"].name),
            ("count", wording["plural"], outcome["gene_count"]),
            ("genetic", wording["column"], outcome["score_column"]),
            (
                "genetic",
                "Chromosomes in covariance metadata",
                ", ".join(outcome["raw_chromosomes"]),
            ),
            (
                "success",
                "Gene identity and order",
                "identical in .genes.out and .genes.raw",
            ),
            (
                "success",
                "MAGMA technical-covariate metadata",
                "finite and strictly positive",
            ),
            (
                "success",
                "MAGMA covariance blocks",
                "%d finite symmetric blocks" % outcome["covariance_blocks"],
            ),
            ("success", "Gene-score input validation", "passed"),
        ]
        if outcome.get("annotated") is None:
            fields.append((
                "info", "Annotated MAGMA enrichment", "not supplied",
            ))
        else:
            fields.extend((
                (
                    "info", "Annotated MAGMA gene-results file",
                    outcome["annotated"]["path"].name,
                ),
                (
                    "success", "Annotated MAGMA consistency",
                    "gene IDs, chromosomes, and Z statistics validated",
                ),
            ))
        return fields
    fields = [
        ("analysis", wording["heading"]),
        ("genetic", "Scores used to fit PoPS", wording["source"]),
        ("info", "Custom gene-score file", outcome["path"].name),
        ("count", wording["plural"], outcome["gene_count"]),
        ("genetic", wording["column"], outcome["score_column"]),
    ]
    if outcome["covariates_path"] is None:
        fields.append(("info", "Custom gene-score covariates", "not requested"))
    else:
        fields.extend((
            (
                "info", "Custom gene-score covariates file",
                outcome["covariates_path"].name,
            ),
            ("count", "Custom gene-score covariates", outcome["covariate_count"]),
            ("success", "Gene-score/covariate gene IDs", "identical"),
        ))
    covariance = outcome.get("covariance")
    if covariance is None:
        fields.append(("info", "Custom gene-error covariance", "not requested"))
    else:
        fields.extend((
            (
                "info", "Custom gene-error covariance file",
                covariance["path"].name,
            ),
            (
                "count",
                "Custom gene-error covariance dimensions",
                "%d × %d" % (covariance["dimension"], covariance["dimension"]),
            ),
            ("analysis", "Covariance representation", covariance["representation"]),
            (
                "success",
                "Covariance validation",
                "finite, symmetric, and positive definite",
            ),
        ))
    fields.append(("success", "Gene-score input validation", "passed"))
    return fields


def _annotation_validation_fields(
    module, annotation: dict,
) -> list[tuple[str, str, object]]:
    schema = module.input_schema
    return [
        ("analysis", "PoPS gene-location annotation"),
        ("info", "Reference file", annotation["path"].name),
        ("genetic", "Declared genome build", module.genome_build.value),
        ("count", "Annotated genes", annotation["gene_count"]),
        (
            "genetic",
            "Chromosomes represented",
            ", ".join(annotation["chromosomes"]),
        ),
        ("genetic", "Gene-ID column", schema.gene_annotation_id_column),
        ("genetic", "Chromosome column", schema.gene_annotation_chromosome_column),
        (
            "genetic",
            "Transcription-start-site column",
            schema.gene_annotation_tss_column,
        ),
        (
            "count",
            "Transcription-start-site range",
            "%g–%g" % (annotation["tss_minimum"], annotation["tss_maximum"]),
        ),
        (
            "success" if annotation["name_column_present"] else "info",
            "Optional gene-name column",
            (
                schema.gene_annotation_name_column
                if annotation["name_column_present"] else "not supplied"
            ),
        ),
        ("success", "Gene-location validation", "passed"),
    ]


def _feature_validation_fields(
    module, features: dict,
) -> list[tuple[str, str, object]]:
    fields = [
        ("analysis", "PoPS feature matrices"),
        ("info", "Feature prefix", Path(module.feature_matrix_prefix).name),
        ("info", "Feature-row file", features["rows_path"].name),
        ("count", "Feature-matrix genes", features["row_count"]),
        ("count", "Matrix chunks", len(features["chunks"])),
        ("count", "Total features", features["feature_count"]),
        (
            "success",
            "Column companion files",
            "%s · %d/%d present and non-empty"
            % (
                module.input_schema.feature_columns_pattern,
                len(features["chunks"]),
                module.feature_matrix_chunks,
            ),
        ),
        (
            "success",
            "Matrix companion files",
            "%s · %d/%d present and non-empty"
            % (
                module.input_schema.feature_matrix_pattern,
                len(features["chunks"]),
                module.feature_matrix_chunks,
            ),
        ),
    ]
    for metrics in features["chunks"]:
        fields.append((
            "analysis",
            "Chunk %d dimensions" % metrics["chunk"],
            "%s genes × %s features · %s"
            % (metrics["rows"], metrics["columns"], metrics["dtype"]),
        ))
    fields.extend((
        ("success", "Matrix dimensions", "match row and column companion files"),
        ("success", "Matrix numeric-value validation", "finite numeric values"),
        ("success", "Feature-matrix validation", "passed"),
    ))
    return fields


def _feature_control_validation_fields(
    controls: dict,
) -> list[tuple[str, str, object]]:
    fields = [("analysis", "PoPS feature controls")]
    for key, label in (
        ("feature_subset", "Feature-subset file"),
        ("controls", "Control-features file"),
    ):
        resource = controls[key]
        if resource is None:
            fields.append(("info", label, "not requested"))
        else:
            fields.extend((
                ("info", label, resource["path"].name),
                (
                    "count",
                    "%s entries" % label.removesuffix(" file"),
                    resource["count"],
                ),
                (
                    "success",
                    "%s compatibility" % label.removesuffix(" file"),
                    "all names found in feature matrices",
                ),
            ))
    retained = controls["controls_retained_by_subset"]
    if retained is not None:
        fields.append((
            "success",
            "Controls retained by feature subset",
            "%d/%d" % (retained, controls["controls"]["count"]),
        ))
    fields.append(("success", "Feature-control validation", "passed"))
    return fields


def _gene_compatibility_fields(
    module, features: dict, annotation: dict, outcome: dict,
) -> list[tuple[str, str, object]]:
    compatibility = outcome["compatibility"]
    retained = compatibility["retained_target_genes"]
    original = compatibility["original_target_genes"]
    excluded = compatibility["excluded_target_genes"]
    fields = [
        ("analysis", "PoPS gene-identifier compatibility"),
        ("count", "Original genes with input scores", original),
        ("genetic", "Gene-location identifiers", annotation["gene_count"]),
        ("genetic", "Feature-row identifiers", features["row_count"]),
        (
            "success",
            "Genes with scores shared by all inputs",
            "%s/%s (%s)"
            % (retained, original, format_percentage(retained, original)),
        ),
        (
            "warning" if compatibility["absent_from_annotation"] else "success",
            "Scored genes absent from annotation",
            compatibility["absent_from_annotation"],
        ),
        (
            "warning" if compatibility["absent_from_features"] else "success",
            "Scored genes absent from feature rows",
            compatibility["absent_from_features"],
        ),
        (
            "warning" if excluded else "success",
            "Scored genes excluded",
            excluded,
        ),
        ("info", "Required compatible genes", module.minimum_gene_count),
    ]
    if outcome.get("pipeline_genome_build") is not None:
        fields.append((
            "success",
            "GWAS/MAGMA and PoPS genome build",
            "%s matched" % outcome["pipeline_genome_build"],
        ))
    fields.extend((
        ("decision", "Gene-universe policy", module.gene_universe_policy),
        (
            "success",
            "Compatibility decision",
            (
                "use all genes with input scores"
                if not excluded else "derive aligned MAGMA inputs"
            ),
        ),
    ))
    return fields


def _progress_outcome_fields(
    stage: str, metrics: dict, *, score_source: str,
) -> list[tuple[str, str, object]]:
    """Translate upstream stage metrics into stable, scientific screen fields."""
    if stage == "target_loading":
        wording = _score_wording(score_source)
        return [
            ("genetic", "Scores used to fit PoPS", wording["source"]),
            ("count", wording["loaded"], metrics["target_genes"]),
            ("count", "Gene-score covariates", metrics["covariates"]),
            ("analysis", "Gene-error covariance", metrics["error_covariance"]),
        ]
    if stage == "covariate_adjustment":
        return [
            ("analysis", "Covariate adjustment", metrics["status"]),
            ("count", "Genes eligible for adjustment", metrics["genes"]),
            ("count", "Covariates projected out", metrics["covariates"]),
            ("genetic", "HLA-region policy", metrics["hla_policy"]),
        ]
    if stage == "feature_selection":
        return [
            ("analysis", "Selection strategy", metrics["strategy"]),
            ("count", "Genes used for feature tests", metrics["genes"]),
            ("count", "Features tested", metrics["tested_features"]),
            ("success", "Features selected", metrics["selected_features"]),
        ]
    if stage == "model_fitting":
        fields = [
            ("analysis", "Prediction model", metrics["method"]),
            ("count", "Genes used to fit the model", metrics["training_genes"]),
            ("count", "Model features", metrics["model_features"]),
        ]
        if metrics.get("selected_cv_alpha") is not None:
            fields.append(
                ("analysis", "Selected regularisation alpha", metrics["selected_cv_alpha"])
            )
        return fields
    if stage == "gene_scoring":
        return [
            (
                "count", "Feature-row genes receiving PoPS scores",
                metrics["genes_scored"],
            ),
            ("count", "Features contributing to scores", metrics["model_features"]),
            ("success", "Output tables written", metrics["outputs_written"]),
        ]
    raise PopsError("Unknown PoPS progress stage: %s" % stage)


def _require_pops_runtime() -> None:
    """Fail early when the published PoPS engine dependency is unavailable."""
    if find_spec("sklearn") is None:
        raise PopsError(
            "PoPS requires scikit-learn. Update the Conda environment from "
            "environment.yml or install the PostGWAS analysis dependencies with "
            "pip install 'postgwas[analysis]'."
        )


def _load_upstream_entrypoints():
    """Import the published PoPS engine only after dependency preflight."""
    _require_pops_runtime()
    try:
        from postgwas.modules.pops.pops import get_pops_args, pops_main
    except ModuleNotFoundError as exc:
        raise PopsError(
            "Cannot import the PoPS runtime because Python package %r is missing."
            % exc.name
        ) from exc
    return get_pops_args, pops_main


def _resolved_configuration(args: argparse.Namespace):
    module_overrides = explicit_overrides(args, {
        "genome_build": "genome_build",
        "pops_genome_build": "genome_build",
        "magma_association_prefix": "magma_association_prefix",
        "magma_annotated_results_file": "magma_annotated_results_file",
        "feature_matrix_prefix": "feature_matrix_prefix",
        "feature_matrix_chunks": "feature_matrix_chunks",
        "gene_universe_policy": "gene_universe_policy",
        "pops_gene_location_file": "gene_location_file",
        "control_features_file": "control_features_file",
        "target_score_file": "target_score_file",
        "target_covariates_file": "target_covariates_file",
        "target_error_covariance_file": "target_error_covariance_file",
        "use_magma_covariates": "use_magma_covariates",
        "use_magma_error_covariance": "use_magma_error_covariance",
        "covariate_projection_chromosomes": "covariate_projection_chromosomes",
        "remove_hla_during_covariate_projection": (
            "remove_hla_during_covariate_projection"
        ),
        "feature_subset_file": "feature_subset_file",
        "feature_selection_chromosomes": "feature_selection_chromosomes",
        "feature_selection_p_cutoff": "feature_selection_p_cutoff",
        "maximum_selected_features": "maximum_selected_features",
        "forward_selected_features": "forward_selected_features",
        "remove_hla_during_feature_selection": (
            "remove_hla_during_feature_selection"
        ),
        "training_chromosomes": "training_chromosomes",
        "remove_hla_during_training": "remove_hla_during_training",
        "method": "method",
        "save_matrix_files": "save_matrix_files",
        "verbose": "verbose",
    })
    global_overrides = explicit_overrides(args, {
        "dataset_id": "run.dataset_id",
        "output_directory": "run.output_directory",
        "threads": "execution.threads",
        "memory_gb": "execution.memory_gb",
        "seed": "execution.random_seed",
        "resume": "run.resume",
        "overwrite": "run.overwrite",
    })
    return load_run_configuration_for_module(
        "pops",
        getattr(args, "run_config", None),
        module_overrides=module_overrides,
        global_overrides=global_overrides,
    )


def _required_file(value: str | None, label: str) -> Path:
    return require_nonempty_file(value, label, error_type=PopsError)


def _prefix_file(prefix: str | None, suffix: str, label: str) -> Path:
    if prefix is None:
        raise PopsError("%s prefix is required." % label)
    return _required_file(str(Path(prefix).expanduser()) + suffix, label)


def _validate_feature_resources(module) -> dict:
    """Reuse a complete, unchanged matrix-bundle check within this pipeline only."""
    schema = module.input_schema
    prefix = module.feature_matrix_prefix
    rows_path = _prefix_file(prefix, schema.feature_rows_suffix, "PoPS feature rows")
    chunks = tuple(
        (
            _prefix_file(
                prefix, schema.feature_columns_pattern.format(chunk=chunk),
                "PoPS feature columns chunk %d" % chunk,
            ),
            _prefix_file(
                prefix, schema.feature_matrix_pattern.format(chunk=chunk),
                "PoPS feature matrix chunk %d" % chunk,
            ),
        )
        for chunk in range(module.feature_matrix_chunks)
    )
    features = validate_feature_matrix_bundle(
        rows_path, chunks, error_type=PopsError,
    )
    bundle_paths = tuple(
        path
        for chunk in features["chunks"]
        for path in (chunk["columns_path"], chunk["matrix_path"])
    )
    register_file_validation_bundle(
        bundle_paths,
        "PoPS feature matrix bundle",
        (
            ("count", "pops_feature_chunks", len(features["chunks"])),
            (
                "success",
                "pops_feature_column_files",
                "%d / %d"
                % (len(features["chunks"]), module.feature_matrix_chunks),
            ),
            (
                "success",
                "pops_feature_matrix_files",
                "%d / %d"
                % (len(features["chunks"]), module.feature_matrix_chunks),
            ),
            ("count", "genes", features["row_count"]),
            ("count", "features", features["feature_count"]),
            (
                "count",
                "pops_matrix_dtypes",
                sorted({chunk["dtype"] for chunk in features["chunks"]}),
            ),
            ("count", "pops_unique_bundle_files", len(bundle_paths)),
            (
                "success",
                "pops_matrix_validation",
                "unique feature names; companion dimensions matched; all "
                "matrix values finite and numeric",
            ),
        ),
        covered_checks=(
            "nonempty unique names",
            "companion dimensions",
            "all values finite and numeric",
            "unique feature names across all chunks",
            "matrix column-count agreement",
            "matrix shape agrees with companion files",
            "all values: finite numeric data",
        ),
        covered_metric_keys=(
            "names", "chunk", "features", "rows", "columns", "dtype",
        ),
    )
    return features


def _validate_feature_control_resources(module, features: dict) -> dict:
    """Validate optional feature subset and control lists against the matrix."""
    all_columns = set(features["feature_names"])
    resources = {}
    for key, value, label in (
        (
            "feature_subset",
            module.feature_subset_file,
            "PoPS feature-subset file",
        ),
        (
            "controls",
            module.control_features_file,
            "PoPS control-features file",
        ),
    ):
        if value is None:
            resources[key] = None
            continue
        path = _required_file(value, label)
        names = read_unique_names(path, label, error_type=PopsError)
        unknown = sorted(set(names.tolist()) - all_columns)
        if unknown:
            raise PopsError(
                "%s contains features absent from the matrix: %s"
                % (label, ", ".join(unknown[:10]))
            )
        resources[key] = {
            "path": path,
            "count": len(names),
            "names": tuple(names.tolist()),
        }
    subset = resources["feature_subset"]
    controls = resources["controls"]
    if subset is not None and controls is not None:
        excluded_controls = sorted(set(controls["names"]) - set(subset["names"]))
        if excluded_controls:
            raise PopsError(
                "PoPS control-features file contains features excluded by the "
                "feature-subset file: %s. Upstream PoPS applies the subset first, "
                "so these controls would otherwise be silently discarded."
                % ", ".join(excluded_controls[:10])
            )
        resources["controls_retained_by_subset"] = controls["count"]
    else:
        resources["controls_retained_by_subset"] = None
    return resources


def _validate_gene_annotation(module) -> dict:
    schema = module.input_schema
    annotation = validate_gene_tss_annotation(
        module.gene_location_file, delimiter=schema.table_delimiter_pattern,
        id_column=schema.gene_annotation_id_column,
        chromosome_column=schema.gene_annotation_chromosome_column,
        tss_column=schema.gene_annotation_tss_column,
        name_column=schema.gene_annotation_name_column,
        require_names=False, require_nonempty_identifiers=False,
        label="PoPS gene annotation", error_type=PopsError,
    )
    configured_chromosomes = {
        value
        for values in (
            module.covariate_projection_chromosomes,
            module.feature_selection_chromosomes,
            module.training_chromosomes,
        )
        if values is not None
        for value in values
    }
    unknown = sorted(configured_chromosomes - set(annotation["chromosomes"]))
    if unknown:
        raise PopsError(
            "Configured PoPS chromosomes are absent from the gene annotation: %s"
            % ", ".join(unknown)
        )
    return annotation


def _clean_table_value(value):
    """Convert pandas/NumPy scalar values to JSON- and CSV-safe Python values."""
    if pd.isna(value):
        return None
    if isinstance(value, np.generic):
        return value.item()
    return value


def _validate_annotated_magma_results(
    module,
    result_identifiers: list[str],
    result_scores: pd.Series,
    gene_chromosomes: dict[str, str],
) -> dict | None:
    """Validate an optional annotated MAGMA table against authoritative results."""
    if module.magma_annotated_results_file is None:
        return None
    schema = module.input_schema.magma_annotated
    path = _required_file(
        module.magma_annotated_results_file,
        "annotated MAGMA gene-results file",
    )
    table = read_pandas_table(path, schema.delimiter, "annotated MAGMA gene results", error_type=PopsError)
    configured_columns = [
        value
        for name, value in schema.model_dump().items()
        if name.endswith("_column")
    ]
    require_table_columns(table, configured_columns, "annotated MAGMA gene results", error_type=PopsError)
    identifiers = table[schema.gene_id_column].astype(str)
    if identifiers.duplicated().any():
        raise PopsError(
            "Annotated MAGMA gene results contain duplicate gene identifiers."
        )
    expected = set(result_identifiers)
    observed = set(identifiers)
    if observed != expected:
        raise PopsError(
            "Annotated MAGMA gene results must contain exactly the genes in "
            ".genes.out (annotated=%d, .genes.out=%d, absent from annotated=%d, "
            "absent from .genes.out=%d). Files: %s; %s"
            % (
                len(observed), len(expected), len(expected - observed),
                len(observed - expected), path,
                str(Path(module.magma_association_prefix).expanduser())
                + module.input_schema.magma_genes_out_suffix,
            )
        )
    indexed = table.assign(_gene_id=identifiers).set_index("_gene_id", drop=True)
    zstats = pd.to_numeric(indexed[schema.zstat_column], errors="coerce")
    if zstats.isna().any() or not np.isfinite(zstats.to_numpy()).all():
        raise PopsError("Annotated MAGMA Z statistics must all be finite.")
    expected_scores = pd.Series(
        result_scores.to_numpy(), index=result_identifiers, dtype=float,
    ).loc[zstats.index]
    if not np.allclose(
        zstats.to_numpy(), expected_scores.to_numpy(), rtol=1e-10, atol=1e-12,
    ):
        raise PopsError(
            "Annotated MAGMA Z statistics do not match the authoritative "
            ".genes.out values."
        )
    for column, label in (
        (schema.pvalue_column, "MAGMA p-values"),
        (schema.bonferroni_pvalue_column, "Bonferroni-adjusted MAGMA p-values"),
        (schema.fdr_pvalue_column, "FDR-adjusted MAGMA p-values"),
    ):
        values = pd.to_numeric(indexed[column], errors="coerce")
        if (
            values.isna().any()
            or not np.isfinite(values.to_numpy()).all()
            or ((values < 0) | (values > 1)).any()
        ):
            raise PopsError("%s must be finite and between zero and one." % label)
    numeric = {
        column: pd.to_numeric(indexed[column], errors="coerce")
        for column in (
            schema.start_column,
            schema.end_column,
            schema.snp_count_column,
            schema.parameter_count_column,
            schema.sample_size_column,
            schema.reference_start_column,
            schema.reference_end_column,
        )
    }
    if any(
        values.isna().any() or not np.isfinite(values.to_numpy()).all()
        for values in numeric.values()
    ):
        raise PopsError(
            "Annotated MAGMA positions, counts, and sample sizes must be finite "
            "numeric values."
        )
    for column in (
        schema.snp_count_column,
        schema.parameter_count_column,
        schema.sample_size_column,
    ):
        if (numeric[column] <= 0).any():
            raise PopsError(
                "Annotated MAGMA gene counts and sample sizes must be positive."
            )
    for start_column, end_column, label in (
        (schema.start_column, schema.end_column, "MAGMA result"),
        (
            schema.reference_start_column,
            schema.reference_end_column,
            "MAGMA reference",
        ),
    ):
        if (
            (numeric[start_column] < 0).any()
            or (numeric[end_column] < numeric[start_column]).any()
        ):
            raise PopsError(
                "%s gene intervals must have non-negative starts and ends not "
                "less than starts." % label
            )
    chromosomes = indexed[schema.chromosome_column].astype(str)
    mismatches = [
        gene
        for gene, chromosome in chromosomes.items()
        if chromosome != gene_chromosomes[gene]
    ]
    if mismatches:
        raise PopsError(
            "Annotated MAGMA chromosome values disagree with .genes.raw for %d "
            "genes. Example genes: %s"
            % (
                len(mismatches),
                _example_ids(mismatches, module.reporting.top_gene_count),
            )
        )
    records = {
        gene: {
            column: _clean_table_value(value)
            for column, value in row.items()
            if column in configured_columns
        }
        for gene, row in indexed.iterrows()
    }
    return {"path": path, "records": records, "columns": configured_columns}


def _validate_magma(module) -> dict:
    schema = module.input_schema
    output_path = _prefix_file(
        module.magma_association_prefix,
        schema.magma_genes_out_suffix,
        "MAGMA gene results",
    )
    raw_path = _prefix_file(
        module.magma_association_prefix,
        schema.magma_genes_raw_suffix,
        "MAGMA raw gene results",
    )
    table = read_pandas_table(output_path, schema.table_delimiter_pattern, "MAGMA gene results", error_type=PopsError)
    require_table_columns(
        table,
        [schema.magma_gene_id_column, schema.magma_score_column],
        "MAGMA gene results", error_type=PopsError)
    identifiers = table[schema.magma_gene_id_column].astype(str)
    scores = pd.to_numeric(table[schema.magma_score_column], errors="coerce")
    if identifiers.duplicated().any():
        raise PopsError("MAGMA gene results contain duplicate gene identifiers.")
    if scores.isna().any() or not np.isfinite(scores.to_numpy()).all():
        raise PopsError("MAGMA gene Z statistics must all be finite.")
    result_identifiers = identifiers.tolist()
    _, raw_records = _read_magma_raw_records(raw_path)
    for row_number, record in enumerate(raw_records, 1):
        for field_index, field_name in (
            (4, "NSNPS"),
            (5, "NPARAM"),
            (7, "MAC"),
        ):
            try:
                value = float(record[field_index])
            except ValueError as exc:
                raise PopsError(
                    "MAGMA raw gene results row %d has a non-numeric %s value: %s"
                    % (row_number, field_name, raw_path)
                ) from exc
            if not np.isfinite(value) or value <= 0:
                raise PopsError(
                    "MAGMA raw gene results row %d has a non-positive or "
                    "non-finite %s value: %s" % (row_number, field_name, raw_path)
                )
    raw_identifiers = [record[0] for record in raw_records]
    if raw_identifiers != result_identifiers:
        result_set = set(result_identifiers)
        raw_set = set(raw_identifiers)
        example_count = module.reporting.top_gene_count
        absent_from_raw = sorted(result_set - raw_set)
        absent_from_results = sorted(raw_set - result_set)
        if absent_from_raw or absent_from_results:
            raise PopsError(
                "MAGMA .genes.out and .genes.raw contain different gene sets "
                "(.genes.out=%d, .genes.raw=%d, absent from .genes.raw=%d, "
                "absent from .genes.out=%d). Example genes absent from .genes.raw: "
                "%s. Example genes absent from .genes.out: %s. Files: %s; %s"
                % (
                    len(result_identifiers), len(raw_identifiers),
                    len(absent_from_raw), len(absent_from_results),
                    _example_ids(absent_from_raw, example_count),
                    _example_ids(absent_from_results, example_count),
                    output_path, raw_path,
                )
            )
        mismatch = next(
            index
            for index, (result_id, raw_id) in enumerate(
                zip(result_identifiers, raw_identifiers), 1,
            )
            if result_id != raw_id
        )
        raise PopsError(
            "MAGMA .genes.out and .genes.raw gene order differs at data row %d "
            "(.genes.out=%s, .genes.raw=%s). PoPS covariance metadata requires "
            "identical gene order. Files: %s; %s"
            % (
                mismatch, result_identifiers[mismatch - 1],
                raw_identifiers[mismatch - 1], output_path, raw_path,
            )
        )
    covariance_blocks = _magma_covariance_blocks(raw_path, result_identifiers)
    gene_chromosomes = {record[0]: record[1] for record in raw_records}
    annotated = _validate_annotated_magma_results(
        module, result_identifiers, scores, gene_chromosomes,
    )
    result_records = {
        gene: {
            column: _clean_table_value(value)
            for column, value in row.items()
        }
        for gene, row in table.assign(_gene_id=identifiers).set_index(
            "_gene_id", drop=True,
        ).iterrows()
    }
    return {
        "genes": set(identifiers),
        "gene_ids": result_identifiers,
        "gene_count": len(identifiers),
        "score_column": schema.magma_score_column,
        "raw_chromosomes": sorted(set(record[1] for record in raw_records)),
        "covariance_blocks": len(covariance_blocks),
        "label": "MAGMA",
        "output_path": output_path,
        "raw_path": raw_path,
        "gene_chromosomes": gene_chromosomes,
        "target_scores": dict(zip(result_identifiers, scores)),
        "target_rows": {
            gene: row_number
            for row_number, gene in enumerate(result_identifiers, 1)
        },
        "result_records": result_records,
        "annotated": annotated,
    }


def _read_magma_raw_records(path: Path) -> tuple[list[str], list[list[str]]]:
    """Read the two MAGMA headers and validated raw records."""
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise PopsError("Cannot read MAGMA raw gene results %s: %s" % (path, exc)) from exc
    if len(lines) < 3:
        raise PopsError("MAGMA raw gene results contain no gene records: %s" % path)
    headers = lines[:2]
    records = []
    for row_number, line in enumerate(lines[2:], 1):
        if not line.strip():
            raise PopsError(
                "MAGMA raw gene results contain a blank data row at row %d: %s"
                % (row_number, path)
            )
        fields = line.split()
        if len(fields) < 9:
            raise PopsError(
                "MAGMA raw gene results data row %d has %d fields; at least 9 "
                "metadata fields are required: %s"
                % (row_number, len(fields), path)
            )
        records.append(fields)
    identifiers = [record[0] for record in records]
    if len(identifiers) != len(set(identifiers)):
        raise PopsError(
            "MAGMA raw gene results contain duplicate gene identifiers: %s" % path
        )
    chromosomes = [record[1] for record in records]
    chromosome_blocks = [chromosome for chromosome, _ in groupby(chromosomes)]
    if len(chromosome_blocks) != len(set(chromosome_blocks)):
        raise PopsError(
            "MAGMA raw chromosomes must occur in contiguous blocks: %s" % path
        )
    return headers, records


def _read_magma_raw_gene_ids(path: Path) -> list[str]:
    """Read MAGMA raw gene IDs without materialising its covariance matrices."""
    _, records = _read_magma_raw_records(path)
    return [record[0] for record in records]


def _example_ids(identifiers, limit: int) -> str:
    values = list(identifiers)[:limit]
    return ", ".join(values) if values else "none"


def _validate_target(module) -> dict:
    schema = module.input_schema
    path = _required_file(module.target_score_file, "PoPS custom gene-score file")
    table = read_pandas_table(
        path, schema.target_table_delimiter, "PoPS custom gene scores", error_type=PopsError)
    require_table_columns(
        table,
        [schema.target_gene_id_column, schema.target_score_column],
        "PoPS custom gene scores", error_type=PopsError)
    identifiers = table[schema.target_gene_id_column].astype(str)
    scores = pd.to_numeric(table[schema.target_score_column], errors="coerce")
    if identifiers.duplicated().any():
        raise PopsError("PoPS custom gene scores contain duplicate gene identifiers.")
    if scores.isna().any() or not np.isfinite(scores.to_numpy()).all():
        raise PopsError("PoPS custom gene scores must all be finite.")
    covariate_path = None
    covariate_count = 0
    if module.target_covariates_file is not None:
        covariate_path = _required_file(
            module.target_covariates_file,
            "PoPS custom gene-score covariates file",
        )
        covariates = read_pandas_table(
            covariate_path,
            schema.target_table_delimiter,
            "PoPS custom gene-score covariates", error_type=PopsError)
        require_table_columns(
            covariates,
            [schema.target_gene_id_column],
            "PoPS custom gene-score covariates", error_type=PopsError)
        covariate_ids = covariates[schema.target_gene_id_column].astype(str)
        if covariate_ids.duplicated().any() or set(covariate_ids) != set(identifiers):
            raise PopsError(
                "Custom gene-score covariates must contain each scored gene "
                "exactly once."
            )
        numeric = covariates.drop(columns=[schema.target_gene_id_column]).apply(
            pd.to_numeric, errors="coerce",
        )
        if (
            numeric.shape[1] == 0
            or numeric.isna().any().any()
            or not np.isfinite(numeric.to_numpy()).all()
        ):
            raise PopsError(
                "Custom gene-score covariates must contain finite numeric columns."
            )
        covariate_count = numeric.shape[1]
    return {
        "genes": set(identifiers),
        "gene_ids": identifiers.tolist(),
        "gene_count": len(identifiers),
        "label": "Custom target",
        "path": path,
        "score_column": schema.target_score_column,
        "covariates_path": covariate_path,
        "covariate_count": covariate_count,
        "target_scores": dict(zip(identifiers, scores)),
        "target_rows": {
            gene: row_number
            for row_number, gene in enumerate(identifiers, 1)
        },
        "result_records": {},
        "annotated": None,
    }


def _gene_compatibility_decision(
    gene: str,
    annotation_genes: set[str],
    feature_genes: set[str],
    retained_genes: set[str],
) -> str:
    """Return the shared retained/excluded decision for one input-score gene."""
    if gene in retained_genes:
        return "retained"
    if gene not in annotation_genes and gene not in feature_genes:
        return "missing_annotation_and_features"
    if gene not in annotation_genes:
        return "missing_annotation"
    return "missing_features"


def _validate_gene_universes(module, features: dict, annotation: dict, outcome) -> dict:
    """Require input-score, annotation, and feature IDs to satisfy indexing."""
    annotation_genes = annotation["genes"]
    feature_genes = set(features["rows"])
    feature_only = sorted(feature_genes - annotation_genes)
    example_count = module.reporting.top_gene_count
    if feature_only:
        raise PopsError(
            "PoPS feature rows contain genes absent from the gene annotation "
            "(annotation=%d, feature rows=%d, absent from annotation=%d). "
            "Example genes absent from annotation: %s. Every feature-row gene "
            "can receive a PoPS prediction and therefore must be annotated. "
            "Files: annotation=%s; feature rows=%s"
            % (
                len(annotation_genes), len(feature_genes), len(feature_only),
                _example_ids(feature_only, example_count), annotation["path"],
                features["rows_path"],
            )
        )
    if len(feature_genes) < module.minimum_gene_count:
        raise PopsError(
            "The PoPS feature matrix contains only %d genes; the configured "
            "minimum is %d."
            % (len(feature_genes), module.minimum_gene_count)
        )
    if outcome is None:
        return {
            "reference_genes": len(feature_genes),
            "shared_target_genes": None,
        }

    missing_annotation = sorted(outcome["genes"] - annotation_genes)
    missing_features = sorted(outcome["genes"] - feature_genes)
    shared = outcome["genes"] & annotation_genes & feature_genes
    if outcome["label"] == "MAGMA":
        chromosome_mismatches = sorted(
            gene
            for gene in shared
            if outcome["gene_chromosomes"][gene]
            != annotation["gene_chromosomes"][gene]
        )
        if chromosome_mismatches:
            examples = ", ".join(
                "%s (MAGMA=%s, PoPS=%s)" % (
                    gene,
                    outcome["gene_chromosomes"][gene],
                    annotation["gene_chromosomes"][gene],
                )
                for gene in chromosome_mismatches[:example_count]
            )
            raise PopsError(
                "MAGMA raw metadata and the PoPS annotation disagree on "
                "chromosome for %d shared genes. Examples: %s. Gene-universe "
                "intersection cannot repair chromosome disagreement. Files: "
                "MAGMA raw=%s; annotation=%s"
                % (
                    len(chromosome_mismatches), examples, outcome["raw_path"],
                    annotation["path"],
                )
            )
    missing_annotation_set = set(missing_annotation)
    missing_features_set = set(missing_features)
    absent_from_both = missing_annotation_set & missing_features_set
    compatibility = {
        "policy": module.gene_universe_policy,
        "original_target_genes": outcome["gene_count"],
        "retained_target_genes": len(shared),
        "excluded_target_genes": outcome["gene_count"] - len(shared),
        "retained_percent": 100 * len(shared) / outcome["gene_count"],
        "excluded_percent": (
            100 * (outcome["gene_count"] - len(shared)) / outcome["gene_count"]
        ),
        "absent_from_annotation": len(missing_annotation),
        "absent_from_features": len(missing_features),
        "absent_from_both": len(absent_from_both),
        "absent_only_from_annotation": len(
            missing_annotation_set - missing_features_set
        ),
        "absent_only_from_features": len(
            missing_features_set - missing_annotation_set
        ),
        "retained_gene_ids": [
            gene for gene in outcome["gene_ids"] if gene in shared
        ],
        "excluded_gene_ids": [
            gene for gene in outcome["gene_ids"] if gene not in shared
        ],
        "missing_annotation_gene_ids": missing_annotation,
        "missing_feature_gene_ids": missing_features,
    }
    outcome["compatibility"] = compatibility
    incompatible = bool(missing_annotation or missing_features)
    unsupported_intersection = (
        module.gene_universe_policy == "intersect"
        and outcome["label"] != "MAGMA"
    )
    if len(shared) < module.minimum_gene_count or unsupported_intersection or (
        incompatible and module.gene_universe_policy == "strict"
    ):
        outcome_path = outcome.get("output_path", outcome.get("path"))
        score_description = _score_wording(outcome["label"])["plural"]
        policy_guidance = (
            "Custom gene-score intersection is not supported; provide gene-score, "
            "covariate, and covariance files already aligned to the PoPS gene "
            "universe."
            if unsupported_intersection
            else "Use matching resources or explicitly select "
            "--gene-universe-policy intersect to derive aligned MAGMA inputs."
        )
        raise PopsError(
            "PoPS input files are scientifically incompatible. %s=%d; "
            "PoPS annotation genes=%d; PoPS feature-row genes=%d; shared genes=%d; "
            "%s absent from annotation=%d; %s absent from feature "
            "rows=%d. Example genes absent from annotation: %s. Example genes "
            "absent from feature rows: %s. This usually indicates different gene "
            "annotation releases. %s Files: outcome=%s; annotation=%s; feature "
            "rows=%s"
            % (
                score_description, outcome["gene_count"], len(annotation_genes),
                len(feature_genes), len(shared), score_description,
                len(missing_annotation), score_description, len(missing_features),
                _example_ids(missing_annotation, example_count),
                _example_ids(missing_features, example_count), policy_guidance,
                outcome_path, annotation["path"], features["rows_path"],
            )
        )
    return {
        "reference_genes": len(feature_genes),
        "shared_target_genes": len(shared),
    }


def _load_magma_covariance_parser():
    """Load the published MAGMA raw parser after dependency preflight."""
    _require_pops_runtime()
    try:
        from postgwas.modules.pops.pops import munge_magma_covariance_metadata
    except ModuleNotFoundError as exc:
        raise PopsError(
            "Cannot import the PoPS MAGMA parser because Python package %r is "
            "missing." % exc.name
        ) from exc
    return munge_magma_covariance_metadata


def _magma_covariance_blocks(path: Path, expected_gene_ids: list[str]):
    """Parse and validate MAGMA chromosome covariance blocks using upstream code."""
    parser = _load_magma_covariance_parser()
    try:
        sigmas, metadata = parser(str(path))
    except (AssertionError, IndexError, OSError, TypeError, ValueError) as exc:
        raise PopsError(
            "Cannot reconstruct MAGMA covariance metadata from %s: %s"
            % (path, exc)
        ) from exc
    parsed_ids = metadata["GENE"].astype(str).tolist()
    if parsed_ids != expected_gene_ids:
        raise PopsError(
            "Published PoPS MAGMA parsing changed gene identity or order in %s."
            % path
        )
    for block_number, sigma in enumerate(sigmas, 1):
        if (
            sigma.ndim != 2
            or sigma.shape[0] != sigma.shape[1]
            or not np.isfinite(sigma).all()
            or not np.allclose(sigma, sigma.T)
        ):
            raise PopsError(
                "MAGMA covariance block %d is not a finite symmetric matrix: %s"
                % (block_number, path)
            )
    return sigmas


def _raw_record_blocks(records: list[list[str]]) -> list[list[int]]:
    """Return original row indices grouped by contiguous chromosome blocks."""
    return [
        [index for index, _ in rows]
        for _, rows in groupby(
            enumerate(records), key=lambda item: item[1][1],
        )
    ]


def _write_magma_raw_subset(
    headers: list[str],
    records: list[list[str]],
    sigmas,
    selected_genes: set[str],
    path: Path,
) -> None:
    """Write a covariance-preserving principal subset in MAGMA raw format."""
    blocks = _raw_record_blocks(records)
    if len(blocks) != len(sigmas):
        raise PopsError(
            "MAGMA raw chromosome blocks (%d) do not match covariance blocks (%d)."
            % (len(blocks), len(sigmas))
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with path.open("w", encoding="utf-8") as handle:
            handle.write("\n".join(headers) + "\n")
            for block_indices, sigma in zip(blocks, sigmas):
                if sigma.shape != (len(block_indices), len(block_indices)):
                    raise PopsError(
                        "MAGMA covariance block shape %s does not match its %d "
                        "raw rows." % (sigma.shape, len(block_indices))
                    )
                selected_local = [
                    local_index
                    for local_index, record_index in enumerate(block_indices)
                    if records[record_index][0] in selected_genes
                ]
                if not selected_local:
                    continue
                subset = sigma[np.ix_(selected_local, selected_local)]
                for subset_index, local_index in enumerate(selected_local):
                    record = records[block_indices[local_index]]
                    fields = list(record[:9])
                    correlations = subset[subset_index, :subset_index]
                    nonzero = np.flatnonzero(correlations)
                    if nonzero.size:
                        fields.extend(
                            "%.17g" % value
                            for value in correlations[int(nonzero[0]):]
                        )
                    handle.write(" ".join(fields) + "\n")
    except OSError as exc:
        raise PopsError("Cannot write derived MAGMA raw file %s: %s" % (path, exc)) from exc


def _derived_magma_prefix(module, genes_out: Path, genes_raw: Path) -> str:
    """Validate paired configured suffixes and return their shared prefix."""
    schema = module.input_schema
    output_text = str(genes_out)
    if not output_text.endswith(schema.magma_genes_out_suffix):
        raise PopsError(
            "Configured derived MAGMA output does not end with %s: %s"
            % (schema.magma_genes_out_suffix, genes_out)
        )
    prefix = output_text[:-len(schema.magma_genes_out_suffix)]
    if str(genes_raw) != prefix + schema.magma_genes_raw_suffix:
        raise PopsError(
            "Configured derived MAGMA .genes.out and .genes.raw paths do not "
            "share one prefix: %s; %s" % (genes_out, genes_raw)
        )
    return prefix


def _prepare_intersected_magma(
    module,
    annotation: dict,
    features: dict,
    outcome: dict,
    paths: dict[str, Path],
    published_paths: dict[str, Path],
) -> tuple[str, dict]:
    """Create aligned retained/excluded MAGMA pairs and complete provenance."""
    compatibility = outcome["compatibility"]
    retained = set(compatibility["retained_gene_ids"])
    excluded = set(compatibility["excluded_gene_ids"])
    output_table = read_pandas_table(
        outcome["output_path"],
        module.input_schema.table_delimiter_pattern,
        "MAGMA gene results", error_type=PopsError)
    gene_column = module.input_schema.magma_gene_id_column
    output_ids = output_table[gene_column].astype(str)
    if output_ids.tolist() != outcome["gene_ids"]:
        raise PopsError("MAGMA gene-result order changed after input preflight.")

    retained_out = paths["compatible_genes_out"]
    retained_raw = paths["compatible_genes_raw"]
    excluded_out = paths["excluded_genes_out"]
    excluded_raw = paths["excluded_genes_raw"]
    output_table.loc[output_ids.isin(retained)].to_csv(
        retained_out, sep="\t", index=False,
    )
    output_table.loc[output_ids.isin(excluded)].to_csv(
        excluded_out, sep="\t", index=False,
    )

    headers, records = _read_magma_raw_records(outcome["raw_path"])
    raw_ids = [record[0] for record in records]
    sigmas = _magma_covariance_blocks(outcome["raw_path"], raw_ids)
    _write_magma_raw_subset(headers, records, sigmas, retained, retained_raw)
    _write_magma_raw_subset(headers, records, sigmas, excluded, excluded_raw)
    retained_prefix = _derived_magma_prefix(module, retained_out, retained_raw)
    _derived_magma_prefix(module, excluded_out, excluded_raw)

    retained_ids = _read_magma_raw_gene_ids(retained_raw)
    if retained_ids != compatibility["retained_gene_ids"]:
        raise PopsError(
            "Derived compatible MAGMA raw file failed gene-order validation: %s"
            % retained_raw
        )
    if excluded and _read_magma_raw_gene_ids(excluded_raw) != compatibility[
        "excluded_gene_ids"
    ]:
        raise PopsError(
            "Derived excluded MAGMA raw file failed gene-order validation: %s"
            % excluded_raw
        )

    annotation_genes = annotation["genes"]
    feature_genes = set(features["rows"])
    audit_rows = []
    chromosome_statistics = {}
    for row_number, record in enumerate(records, 1):
        gene, chromosome = record[0], record[1]
        in_annotation = gene in annotation_genes
        in_features = gene in feature_genes
        is_retained = gene in retained
        reason = _gene_compatibility_decision(
            gene, annotation_genes, feature_genes, retained,
        )
        audit_rows.append({
            "magma_row": row_number,
            "gene_id": gene,
            "chromosome": chromosome,
            "present_in_pops_annotation": in_annotation,
            "present_in_feature_rows": in_features,
            "retained_for_pops": is_retained,
            "decision": reason,
        })
        counts = chromosome_statistics.setdefault(
            chromosome, {"original": 0, "retained": 0, "excluded": 0},
        )
        counts["original"] += 1
        counts["retained" if is_retained else "excluded"] += 1
    pd.DataFrame(audit_rows).to_csv(
        paths["gene_compatibility_table"], sep="\t", index=False,
    )

    report = {
        key: value
        for key, value in compatibility.items()
        if not key.endswith("_gene_ids")
    }
    report.update({
        "decision": (
            "MAGMA genes with Z-scores were restricted to genes present in both "
            "the PoPS annotation and feature rows."
        ),
        "scientific_effect": (
            "Feature selection and model fitting use the retained MAGMA genes "
            "with Z-scores; "
            "the resulting PoPS scores may differ from a run using a matched "
            "MAGMA gene annotation."
        ),
        "covariance_handling": (
            "Retained and excluded .genes.raw files contain covariance-preserving "
            "principal submatrices reconstructed from the original MAGMA blocks."
        ),
        "original_files_unchanged": True,
        "excluded_gene_examples": compatibility["excluded_gene_ids"][
            :module.reporting.top_gene_count
        ],
        "exclusion_reason_counts": {
            reason: sum(row["decision"] == reason for row in audit_rows)
            for reason in (
                "retained", "missing_annotation", "missing_features",
                "missing_annotation_and_features",
            )
        },
        "chromosomes": chromosome_statistics,
        "files": {
            "original_genes_out": str(outcome["output_path"]),
            "original_genes_raw": str(outcome["raw_path"]),
            "compatible_genes_out": str(published_paths["compatible_genes_out"]),
            "compatible_genes_raw": str(published_paths["compatible_genes_raw"]),
            "excluded_genes_out": str(published_paths["excluded_genes_out"]),
            "excluded_genes_raw": str(published_paths["excluded_genes_raw"]),
            "gene_audit_table": str(published_paths["gene_compatibility_table"]),
            "compatibility_report": str(
                published_paths["gene_compatibility_report"]
            ),
        },
    })
    write_yaml_report(report, paths["gene_compatibility_report"])
    return retained_prefix, report


def _load_target_covariance(path: Path, expected_size: int) -> tuple[np.ndarray, bool]:
    converted = False
    try:
        covariance = sparse.load_npz(path).toarray()
    except (AttributeError, OSError, TypeError, ValueError):
        try:
            covariance = np.load(path, allow_pickle=False)
            converted = True
        except (OSError, TypeError, ValueError) as exc:
            raise PopsError(
                "Cannot read custom gene-error covariance %s: %s" % (path, exc)
            ) from exc
    if covariance.shape != (expected_size, expected_size):
        raise PopsError(
            "Custom gene-error covariance has shape %s; expected (%d, %d)."
            % (covariance.shape, expected_size, expected_size)
        )
    if not np.isfinite(covariance).all() or not np.allclose(covariance, covariance.T):
        raise PopsError("Custom gene-error covariance must be finite and symmetric.")
    try:
        np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise PopsError(
            "Custom gene-error covariance must be positive definite."
        ) from exc
    return covariance, converted


def _validate_resolved_pops_configuration(
    args: argparse.Namespace, configuration, *, pipeline: bool = False,
):
    """Preflight resources for one already resolved PoPS configuration."""
    module = configuration.modules.pops
    target_alternatives = () if pipeline else (
        RequiredAlternative((
            RequiredArgument(
                "--magma-association-prefix",
                "modules.pops.magma_association_prefix",
                module.magma_association_prefix,
            ),
            RequiredArgument(
                "--target-score-file",
                "modules.pops.target_score_file",
                module.target_score_file,
            ),
        )),
    )
    require_resolved_arguments(
        (
            RequiredArgument(
                "--genome-build", "modules.pops.genome_build", module.genome_build,
            ),
            RequiredArgument(
                "--feature-matrix-prefix",
                "modules.pops.feature_matrix_prefix",
                module.feature_matrix_prefix,
            ),
            RequiredArgument(
                "--pops-gene-location-file",
                "modules.pops.gene_location_file",
                module.gene_location_file,
            ),
        ),
        alternatives=target_alternatives,
    )
    _require_pops_runtime()
    features = _validate_feature_resources(module)
    annotation = _validate_gene_annotation(module)
    controls = _validate_feature_control_resources(module, features)
    features["control_resources"] = controls
    if pipeline:
        outcome = None
    elif module.magma_association_prefix is not None:
        outcome = _validate_magma(module)
    elif module.target_score_file is not None:
        outcome = _validate_target(module)
    else:
        raise PopsError(
            "Provide --magma-association-prefix PREFIX or --target-score-file "
            "PATH (a custom gene-score table)."
        )
    _validate_gene_universes(module, features, annotation, outcome)
    if outcome is not None:
        if module.target_error_covariance_file is not None:
            covariance_path = _required_file(
                module.target_error_covariance_file,
                "PoPS custom gene-error covariance",
            )
            _load_target_covariance(covariance_path, outcome["gene_count"])
    return configuration, features, annotation, outcome


def validate_pops_configuration(args: argparse.Namespace, *, pipeline: bool = False):
    """Resolve and preflight PoPS resources without running upstream PoPS."""
    return _validate_resolved_pops_configuration(
        args, _resolved_configuration(args), pipeline=pipeline,
    )


def preflight_pops_pipeline(
    args: argparse.Namespace,
    *,
    preflight_evidence=None,
) -> PipelinePreflightEvidence:
    """Fail before upstream pipeline steps if configured PoPS resources are invalid."""
    require_pipeline_input_vcf(preflight_evidence)
    configuration, features, annotation, _ = validate_pops_configuration(
        args, pipeline=True,
    )
    return pipeline_preflight_evidence(
        "pops",
        preflight_evidence,
        resources=(configuration, features, annotation, features["control_resources"]),
        deferred_checks=(
            "Validate the pipeline-generated MAGMA gene-association result.",
        ),
    )


def _flag(arguments: list[str], condition: bool, enabled: str, disabled: str) -> None:
    arguments.append(enabled if condition else disabled)


def _build_upstream_arguments(
    module,
    seed: int,
    output_prefix: Path,
    covariance_path=None,
    magma_prefix: str | None = None,
):
    arguments = [
        "--gene_annot_path", module.gene_location_file,
        "--feature_mat_prefix", module.feature_matrix_prefix,
        "--num_feature_chunks", str(module.feature_matrix_chunks),
        "--out_prefix", str(output_prefix),
    ]
    effective_magma_prefix = magma_prefix or module.magma_association_prefix
    if effective_magma_prefix is not None:
        arguments += ["--magma_prefix", effective_magma_prefix]
    if module.control_features_file is not None:
        arguments += ["--control_features_path", module.control_features_file]
    _flag(
        arguments,
        module.use_magma_covariates,
        "--use_magma_covariates",
        "--ignore_magma_covariates",
    )
    _flag(
        arguments,
        module.use_magma_error_covariance,
        "--use_magma_error_cov",
        "--ignore_magma_error_cov",
    )
    if module.target_score_file is not None:
        arguments += ["--y_path", module.target_score_file]
    if module.target_covariates_file is not None:
        arguments += ["--y_covariates_path", module.target_covariates_file]
    if covariance_path is not None:
        arguments += ["--y_error_cov_path", str(covariance_path)]
    elif module.target_error_covariance_file is not None:
        arguments += ["--y_error_cov_path", module.target_error_covariance_file]
    if module.covariate_projection_chromosomes is not None:
        arguments += [
            "--project_out_covariates_chromosomes",
            *module.covariate_projection_chromosomes,
        ]
    _flag(
        arguments,
        module.remove_hla_during_covariate_projection,
        "--project_out_covariates_remove_hla",
        "--project_out_covariates_keep_hla",
    )
    if module.feature_subset_file is not None:
        arguments += ["--subset_features_path", module.feature_subset_file]
    if module.feature_selection_chromosomes is not None:
        arguments += [
            "--feature_selection_chromosomes",
            *module.feature_selection_chromosomes,
        ]
    arguments += [
        "--feature_selection_p_cutoff",
        str(module.feature_selection_p_cutoff),
    ]
    if module.maximum_selected_features is not None:
        arguments += [
            "--feature_selection_max_num",
            str(module.maximum_selected_features),
        ]
    if module.forward_selected_features is not None:
        arguments += [
            "--feature_selection_fss_num_features",
            str(module.forward_selected_features),
        ]
    _flag(
        arguments,
        module.remove_hla_during_feature_selection,
        "--feature_selection_remove_hla",
        "--feature_selection_keep_hla",
    )
    if module.training_chromosomes is not None:
        arguments += ["--training_chromosomes", *module.training_chromosomes]
    _flag(
        arguments,
        module.remove_hla_during_training,
        "--training_remove_hla",
        "--training_keep_hla",
    )
    arguments += ["--method", module.method, "--random_seed", str(seed)]
    _flag(
        arguments,
        module.save_matrix_files,
        "--save_matrix_files",
        "--no_save_matrix_files",
    )
    _flag(arguments, module.verbose, "--verbose", "--no_verbose")
    return arguments


def _output_paths(output: Path, dataset: str, module) -> tuple[Path, dict[str, Path]]:
    prefix = configured_output_path(
        output,
        module.output_layout.output_prefix,
        error_type=PopsError,
        dataset_id=dataset,
    )
    layout = module.output_layout
    suffixes = {
        "pops_file": layout.predictions_suffix,
        "coefficients_file": layout.coefficients_suffix,
        "marginals_file": layout.marginals_suffix,
        "upstream_log_file": layout.upstream_log_suffix,
        "integrated_results_file": layout.integrated_results_suffix,
        "integrated_report_file": layout.integrated_report_suffix,
    }
    if module.save_matrix_files:
        suffixes["training_data_file"] = layout.training_data_suffix
        suffixes["matrix_data_file"] = layout.matrix_data_suffix
    if module.gene_universe_policy == "intersect":
        suffixes.update({
            "compatible_genes_out": layout.compatible_genes_out_suffix,
            "compatible_genes_raw": layout.compatible_genes_raw_suffix,
            "excluded_genes_out": layout.excluded_genes_out_suffix,
            "excluded_genes_raw": layout.excluded_genes_raw_suffix,
            "gene_compatibility_table": layout.gene_compatibility_table_suffix,
            "gene_compatibility_report": layout.gene_compatibility_report_suffix,
        })
    return prefix, {name: Path(str(prefix) + suffix) for name, suffix in suffixes.items()}


def _complete_outputs(paths: dict[str, Path]) -> bool:
    return all(path.is_file() and path.stat().st_size > 0 for path in paths.values())


def _result(output: Path, prefix: Path, paths: dict[str, Path]) -> dict:
    return {
        "status": "success",
        "pops_file": str(paths["pops_file"]),
        "integrated_results_file": str(paths["integrated_results_file"]),
        "integrated_report_file": str(paths["integrated_report_file"]),
        "output_dir": str(output),
        "out_prefix": str(prefix),
        "published_files": [str(path) for path in paths.values()],
    }


def _boolean_series(values: pd.Series, label: str) -> pd.Series:
    """Validate and normalise an upstream Boolean output column."""
    normalized = values.astype(str).str.strip().str.lower()
    invalid = values.isna() | ~normalized.isin({"true", "false"})
    if invalid.any():
        examples = sorted(set(values.loc[invalid].astype(str)))[:5]
        raise PopsError(
            "%s must contain only true or false values. Invalid examples: %s"
            % (label, ", ".join(examples))
        )
    return normalized.eq("true")


def _validated_predictions(paths, module, annotation, outcome, features) -> dict:
    """Validate the official PoPS prediction table and its model-use flags."""
    schema = module.input_schema
    predictions = read_pandas_table(
        paths["pops_file"], schema.target_table_delimiter, "PoPS predictions", error_type=PopsError)
    gene_column = schema.gene_annotation_id_column
    score_column = schema.prediction_score_column
    required = [
        gene_column,
        score_column,
        schema.prediction_target_column,
        schema.prediction_feature_selection_column,
        schema.prediction_training_column,
    ]
    require_table_columns(predictions, required, "PoPS predictions", error_type=PopsError)
    gene_ids = predictions[gene_column].astype(str)
    scores = pd.to_numeric(predictions[score_column], errors="coerce")
    if gene_ids.duplicated().any():
        raise PopsError("PoPS predictions contain duplicate gene identifiers.")
    if scores.isna().any() or not np.isfinite(scores.to_numpy()).all():
        raise PopsError("PoPS prediction scores must all be finite.")
    feature_ids = [str(gene) for gene in features["rows"]]
    if gene_ids.tolist() != feature_ids:
        missing_predictions = set(feature_ids) - set(gene_ids)
        unexpected_predictions = set(gene_ids) - set(feature_ids)
        if missing_predictions or unexpected_predictions:
            raise PopsError(
                "PoPS predictions must contain every feature-row gene exactly "
                "once (missing predictions=%d, unexpected predictions=%d)."
                % (len(missing_predictions), len(unexpected_predictions))
            )
        raise PopsError(
            "PoPS prediction gene order differs from the feature-row order."
        )
    missing_annotation = set(gene_ids) - annotation["genes"]
    if missing_annotation:
        raise PopsError(
            "Every PoPS prediction gene must occur in the gene annotation "
            "(missing annotation=%d)." % len(missing_annotation)
        )

    targets = pd.to_numeric(
        predictions[schema.prediction_target_column], errors="coerce",
    )
    invalid_targets = (
        predictions[schema.prediction_target_column].notna() & targets.isna()
    )
    if invalid_targets.any() or not np.isfinite(targets.dropna().to_numpy()).all():
        raise PopsError(
            "Non-missing PoPS input gene scores must be finite numeric values."
        )
    retained_genes = set(outcome["compatibility"]["retained_gene_ids"])
    target_gene_ids = set(gene_ids.loc[targets.notna()])
    if target_gene_ids != retained_genes:
        raise PopsError(
            "Non-missing PoPS input gene scores do not match the validated "
            "retained scored genes (prediction scores=%d, retained genes=%d)."
            % (len(target_gene_ids), len(retained_genes))
        )
    if retained_genes:
        target_by_gene = pd.Series(targets.to_numpy(), index=gene_ids)
        observed = target_by_gene.loc[list(retained_genes)].astype(float)
        expected = pd.Series(outcome["target_scores"]).loc[observed.index].astype(float)
        if not np.allclose(
            observed.to_numpy(), expected.to_numpy(), rtol=1e-10, atol=1e-12,
        ):
            raise PopsError(
                "PoPS prediction input-score values do not match the validated "
                "input gene scores."
            )

    covariate_projection = pd.Series(False, index=predictions.index)
    if schema.prediction_covariate_projection_column in predictions.columns:
        covariate_projection = _boolean_series(
            predictions[schema.prediction_covariate_projection_column],
            "PoPS covariate-projection gene flags",
        )
    flags = {
        "covariate_projection": covariate_projection,
        "feature_selection": _boolean_series(
            predictions[schema.prediction_feature_selection_column],
            "PoPS feature-selection gene flags",
        ),
        "model_fitting": _boolean_series(
            predictions[schema.prediction_training_column],
            "PoPS model-fitting gene flags",
        ),
    }
    for label, values in flags.items():
        invalid_genes = set(gene_ids.loc[values]) - retained_genes
        if invalid_genes:
            raise PopsError(
                "PoPS %s flags include %d genes without retained input gene scores."
                % (label.replace("_", "-"), len(invalid_genes))
            )

    projected = None
    has_projected = schema.prediction_projected_target_column in predictions.columns
    has_projection_flag = (
        schema.prediction_covariate_projection_column in predictions.columns
    )
    if has_projected != has_projection_flag:
        raise PopsError(
            "PoPS predictions must contain both covariate-adjusted gene-score values "
            "and covariate-projection gene flags, or neither."
        )
    if has_projected:
        projected = pd.to_numeric(
            predictions[schema.prediction_projected_target_column], errors="coerce",
        )
        invalid_projected = (
            predictions[schema.prediction_projected_target_column].notna()
            & projected.isna()
        )
        if (
            invalid_projected.any()
            or not np.isfinite(projected.dropna().to_numpy()).all()
        ):
            raise PopsError(
                "Non-missing covariate-adjusted PoPS gene scores must be finite."
            )
        projected_gene_ids = set(gene_ids.loc[projected.notna()])
        if projected_gene_ids != retained_genes:
            raise PopsError(
                "Non-missing covariate-adjusted PoPS gene scores do not match "
                "the validated retained genes with input scores."
            )
    return {
        "table": predictions,
        "gene_ids": gene_ids,
        "scores": scores,
        "targets": targets,
        "projected_targets": projected,
        "flags": flags,
    }


def _magma_values(module, outcome, gene: str) -> dict[str, object]:
    schema = module.input_schema.magma_annotated
    if outcome["label"] != "MAGMA" or gene not in outcome["genes"]:
        record = {}
    elif outcome.get("annotated") is not None:
        record = outcome["annotated"]["records"].get(gene, {})
    else:
        record = outcome["result_records"].get(gene, {})
    return {
        "gene_symbol": record.get(schema.gene_symbol_column),
        "chromosome": record.get(schema.chromosome_column),
        "start": record.get(schema.start_column),
        "end": record.get(schema.end_column),
        "snp_count": record.get(schema.snp_count_column),
        "parameter_count": record.get(schema.parameter_count_column),
        "sample_size": record.get(schema.sample_size_column),
        "zstat": record.get(schema.zstat_column),
        "pvalue": record.get(schema.pvalue_column),
        "reference_chromosome": record.get(schema.reference_chromosome_column),
        "reference_start": record.get(schema.reference_start_column),
        "reference_end": record.get(schema.reference_end_column),
        "reference_strand": record.get(schema.reference_strand_column),
        "bonferroni_pvalue": record.get(schema.bonferroni_pvalue_column),
        "fdr_pvalue": record.get(schema.fdr_pvalue_column),
    }


def _integrated_result_data(paths, module, annotation, outcome, features) -> dict:
    """Build one full-union record per scored or input-target gene."""
    validated = _validated_predictions(paths, module, annotation, outcome, features)
    predictions = validated["table"]
    gene_ids = validated["gene_ids"]
    scores = validated["scores"]
    targets = validated["targets"]
    projected = validated["projected_targets"]
    flags = validated["flags"]
    score_order = np.argsort(-scores.to_numpy(), kind="stable")
    ranked_prediction_indices = [int(index) for index in score_order]
    rank_by_gene = {
        str(gene_ids.iloc[index]): rank
        for rank, index in enumerate(ranked_prediction_indices, 1)
    }
    prediction_index = {
        str(gene): index for index, gene in enumerate(gene_ids)
    }
    prediction_genes = set(prediction_index)
    ordered_genes = [str(gene_ids.iloc[index]) for index in ranked_prediction_indices]
    ordered_genes.extend(
        gene for gene in outcome["gene_ids"] if gene not in prediction_genes
    )
    feature_genes = set(str(gene) for gene in features["rows"])
    retained_genes = set(outcome["compatibility"]["retained_gene_ids"])
    columns = module.integrated_results.columns
    statuses = module.integrated_results.status_labels
    records = []
    status_counts = {
        statuses.scored_and_fitted: 0,
        statuses.scored_with_target_not_fitted: 0,
        statuses.scored_without_target: 0,
        statuses.target_excluded_from_pops: 0,
    }
    for gene in ordered_genes:
        index = prediction_index.get(gene)
        scored = index is not None
        target_available = gene in outcome["genes"]
        fitted = bool(flags["model_fitting"].iloc[index]) if scored else None
        if scored and fitted:
            status = statuses.scored_and_fitted
        elif scored and target_available:
            status = statuses.scored_with_target_not_fitted
        elif scored:
            status = statuses.scored_without_target
        else:
            status = statuses.target_excluded_from_pops
        status_counts[status] += 1
        magma = _magma_values(module, outcome, gene)
        gene_symbol = magma["gene_symbol"] or annotation["gene_names"].get(gene)
        row = {
            columns.gene_id: gene,
            columns.gene_symbol: gene_symbol,
            columns.pops_chromosome: annotation["gene_chromosomes"].get(gene),
            columns.pops_tss: _clean_table_value(annotation["gene_tss"].get(gene)),
            columns.target_source: outcome["label"] if target_available else None,
            columns.target_row: outcome["target_rows"].get(gene),
            columns.input_target_score_available: target_available,
            columns.input_target_score: _clean_table_value(
                outcome["target_scores"].get(gene)
            ),
            columns.present_in_pops_annotation: gene in annotation["genes"],
            columns.present_in_feature_rows: gene in feature_genes,
            columns.retained_as_pops_target: gene in retained_genes,
            columns.compatibility_decision: (
                _gene_compatibility_decision(
                    gene,
                    annotation["genes"],
                    feature_genes,
                    retained_genes,
                )
                if target_available else None
            ),
            columns.pops_score_available: scored,
            columns.pops_score: _clean_table_value(scores.iloc[index]) if scored else None,
            columns.pops_rank: rank_by_gene.get(gene),
            columns.pops_target_score_available: (
                bool(pd.notna(targets.iloc[index])) if scored else False
            ),
            columns.pops_target_score: (
                _clean_table_value(targets.iloc[index]) if scored else None
            ),
            columns.pops_adjusted_target_score: (
                _clean_table_value(projected.iloc[index])
                if scored and projected is not None else None
            ),
            columns.used_for_covariate_projection: (
                bool(flags["covariate_projection"].iloc[index]) if scored else None
            ),
            columns.used_for_feature_selection: (
                bool(flags["feature_selection"].iloc[index]) if scored else None
            ),
            columns.used_for_model_fitting: fitted,
            columns.gene_analysis_status: status,
            columns.magma_result_available: (
                outcome["label"] == "MAGMA" and target_available
            ),
            columns.magma_chromosome: magma["chromosome"],
            columns.magma_start: magma["start"],
            columns.magma_end: magma["end"],
            columns.magma_snp_count: magma["snp_count"],
            columns.magma_parameter_count: magma["parameter_count"],
            columns.magma_sample_size: magma["sample_size"],
            columns.magma_zstat: magma["zstat"],
            columns.magma_pvalue: magma["pvalue"],
            columns.magma_reference_chromosome: magma["reference_chromosome"],
            columns.magma_reference_start: magma["reference_start"],
            columns.magma_reference_end: magma["reference_end"],
            columns.magma_reference_strand: magma["reference_strand"],
            columns.magma_bonferroni_pvalue: magma["bonferroni_pvalue"],
            columns.magma_fdr_pvalue: magma["fdr_pvalue"],
        }
        records.append(row)
    return {
        "predictions": validated,
        "records": records,
        "columns": list(columns.model_dump().values()),
        "status_counts": status_counts,
    }


def _write_integrated_results(
    paths, module, annotation, outcome, features, *, dataset: str,
) -> dict:
    integrated = _integrated_result_data(
        paths, module, annotation, outcome, features,
    )
    settings = module.integrated_results
    write_delimited_report(
        integrated["records"],
        paths["integrated_results_file"],
        fieldnames=integrated["columns"],
        delimiter=settings.delimiter,
        null_value=settings.null_value,
    )
    display_counts = {
        "Scored and used for model fitting": integrated["status_counts"][
            settings.status_labels.scored_and_fitted
        ],
        "Scored with input gene score, not used for fitting": integrated[
            "status_counts"
        ][settings.status_labels.scored_with_target_not_fitted],
        "Scored without an input gene score": integrated["status_counts"][
            settings.status_labels.scored_without_target
        ],
        "Input-score gene excluded before PoPS": integrated["status_counts"][
            settings.status_labels.target_excluded_from_pops
        ],
    }
    write_integrated_gene_report(
        integrated["records"],
        integrated["columns"],
        {"status_counts": display_counts},
        dataset_id=dataset,
        tsv_path=paths["integrated_results_file"],
        report_path=paths["integrated_report_file"],
        page_size=settings.html_page_size,
        null_value=settings.null_value,
    )
    return integrated


def _validate_integrated_results(
    paths, module, annotation, outcome, features,
) -> dict:
    expected = _integrated_result_data(
        paths, module, annotation, outcome, features,
    )
    settings = module.integrated_results
    try:
        observed = pd.read_csv(
            paths["integrated_results_file"],
            sep=settings.delimiter,
            dtype=str,
            keep_default_na=False,
        )
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise PopsError(
            "Cannot read integrated PoPS gene results %s: %s"
            % (paths["integrated_results_file"], exc)
        ) from exc
    require_table_columns(observed, expected["columns"], "integrated PoPS gene results", error_type=PopsError)
    if observed.columns.tolist() != expected["columns"]:
        raise PopsError(
            "Integrated PoPS gene-result columns or column order changed."
        )
    for column in expected["columns"]:
        expected_values = [
            settings.null_value if record[column] is None else str(record[column])
            for record in expected["records"]
        ]
        if observed[column].tolist() != expected_values:
            raise PopsError(
                "Integrated PoPS gene results do not match validated source "
                "outputs in column %s." % column
            )
    if (
        not paths["integrated_report_file"].is_file()
        or paths["integrated_report_file"].stat().st_size <= 0
    ):
        raise PopsError("Integrated PoPS HTML report is missing or empty.")
    return expected


def _summarise_outputs(
    paths, module, annotation, outcome, features, *, published_paths=None,
    integrated=None,
) -> dict:
    """Validate published scientific results and build the terminal summary."""
    schema = module.input_schema
    integrated = integrated or _validate_integrated_results(
        paths, module, annotation, outcome, features,
    )
    validated = integrated["predictions"]
    predictions = validated["table"]
    gene_ids = validated["gene_ids"]
    scores = validated["scores"]
    ranked = predictions.assign(score_value=scores, gene_value=gene_ids).nlargest(
        module.reporting.top_gene_count, "score_value",
    )
    top_genes = [
        {
            "rank": rank,
            "gene_id": row.gene_value,
            "gene_name": annotation["gene_names"].get(
                row.gene_value, row.gene_value,
            ),
            "score": float(row.score_value),
        }
        for rank, row in enumerate(ranked.itertuples(index=False), 1)
    ]
    training_count = int(validated["flags"]["model_fitting"].sum())
    selected_features = None
    marginals = read_pandas_table(
        paths["marginals_file"], schema.target_table_delimiter, "PoPS marginals", error_type=PopsError)
    if schema.marginal_selected_column in marginals.columns:
        selected_features = int(_boolean_series(
            marginals[schema.marginal_selected_column],
            "PoPS marginal feature-selection flags",
        ).sum())

    reference_genes = set(str(gene) for gene in features["rows"])
    target_genes = set(outcome["compatibility"]["retained_gene_ids"])
    missing_chromosomes = sorted(
        set(annotation["gene_chromosomes"].values())
        - {
            annotation["gene_chromosomes"][gene]
            for gene in target_genes
            if gene in annotation["gene_chromosomes"]
        }
    )
    score_wording = _score_wording(outcome["label"])
    warnings = []
    if len(target_genes) < len(reference_genes):
        warnings.append(
            "%s cover %d of %d compatible genes (%.1f%%); rankings "
            "were trained on incomplete genome-wide gene outcomes."
            % (
                score_wording["plural"][0].upper() + score_wording["plural"][1:],
                len(target_genes), len(reference_genes),
                100 * len(target_genes) / len(reference_genes),
            )
        )
    if missing_chromosomes:
        warnings.append(
            "No genes with %s were available on chromosome(s): %s."
            % (score_wording["plural"], ", ".join(missing_chromosomes))
        )
    compatibility = {
        key: value
        for key, value in outcome.get("compatibility", {}).items()
        if not key.endswith("_gene_ids")
    }
    if compatibility.get("excluded_target_genes", 0):
        warnings.append(
            "Gene-universe intersection excluded %d of %d MAGMA genes with "
            "Z-scores "
            "(%.1f%% retained). Feature selection and fitting used the retained "
            "MAGMA Z-score gene universe; review the compatibility audit before "
            "interpreting rankings."
            % (
                compatibility["excluded_target_genes"],
                compatibility["original_target_genes"],
                compatibility["retained_percent"],
            )
        )
    if module.gene_universe_policy == "intersect":
        reported = published_paths or paths
        compatibility["audit_table"] = str(reported["gene_compatibility_table"])
        compatibility["report"] = str(reported["gene_compatibility_report"])
        compatibility["compatible_genes_out"] = str(reported["compatible_genes_out"])
        compatibility["compatible_genes_raw"] = str(reported["compatible_genes_raw"])
        compatibility["excluded_genes_out"] = str(reported["excluded_genes_out"])
        compatibility["excluded_genes_raw"] = str(reported["excluded_genes_raw"])
    statuses = module.integrated_results.status_labels
    status_counts = integrated["status_counts"]
    reported = published_paths or paths
    return {
        "genes_scored": len(predictions),
        "score_source": outcome["label"],
        "genes_scored_with_target": int(validated["targets"].notna().sum()),
        "genes_scored_without_target": status_counts[statuses.scored_without_target],
        "genes_scored_with_target_not_fitted": status_counts[
            statuses.scored_with_target_not_fitted
        ],
        "input_target_genes": outcome["gene_count"],
        "input_target_genes_excluded": status_counts[
            statuses.target_excluded_from_pops
        ],
        "target_genes": len(target_genes),
        "compatible_genes": len(reference_genes),
        "training_genes": training_count,
        "selected_features": selected_features,
        "top_genes": top_genes,
        "gene_compatibility": compatibility,
        "integrated_gene_count": len(integrated["records"]),
        "integrated_results_file": str(reported["integrated_results_file"]),
        "integrated_report_file": str(reported["integrated_report_file"]),
        "warnings": warnings,
    }


def _render_summary(
    summary: dict, dataset: str, module, paths, log_path: Path, label_width: int,
) -> str:
    width = label_width
    outer_width = width + 4
    status = "COMPLETED WITH SCIENTIFIC WARNINGS" if summary["warnings"] else "COMPLETED"
    score_wording = _score_wording(summary["score_source"])
    lines = [
        "",
        screen_line("analysis", "PoPS gene-prioritisation summary", indent=2),
        screen_field("info", "Dataset", dataset, indent=6, label_width=outer_width),
        screen_field(
            "warning" if summary["warnings"] else "success",
            "Analysis status", status, indent=6, label_width=outer_width,
        ),
        "",
        screen_line("genetic", "PoPS scoring coverage", indent=6),
        screen_field(
            "count", "Genes receiving finite PoPS scores",
            f'{summary["genes_scored"]:,}', indent=10, label_width=width,
        ),
        screen_field(
            "count", "Scored genes with %s" % score_wording["coverage"],
            "%s of %s" % (
                f'{summary["genes_scored_with_target"]:,}',
                f'{summary["genes_scored"]:,}',
            ),
            indent=10, label_width=width,
        ),
        screen_field(
            "warning" if summary["genes_scored_without_target"] else "success",
            "Scored genes without %s" % score_wording["coverage"],
            f'{summary["genes_scored_without_target"]:,}',
            indent=10, label_width=width,
        ),
        "",
        screen_line("analysis", "PoPS model estimation", indent=6),
        screen_field(
            "count", "Genes used to fit the PoPS model",
            f'{summary["training_genes"]:,}', indent=10, label_width=width,
        ),
        screen_field(
            "count",
            "Genes with %s not used for fitting" % score_wording["coverage"],
            f'{summary["genes_scored_with_target_not_fitted"]:,}',
            indent=10, label_width=width,
        ),
    ]
    if summary["selected_features"] is not None:
        lines.append(screen_field(
            "count", "Features selected", f'{summary["selected_features"]:,}',
            indent=10, label_width=width,
        ))
    compatibility = summary.get("gene_compatibility", {})
    if compatibility.get("policy") == "intersect":
        lines.extend([
            "",
            screen_line("warning", "MAGMA–PoPS gene compatibility", indent=6),
            screen_field(
                "count", "Original MAGMA genes with Z-scores",
                f'{compatibility["original_target_genes"]:,}',
                indent=10, label_width=width,
            ),
            screen_field(
                "success", "Retained MAGMA genes with Z-scores",
                "%s (%.1f%%)" % (
                    f'{compatibility["retained_target_genes"]:,}',
                    compatibility["retained_percent"],
                ),
                indent=10, label_width=width,
            ),
            screen_field(
                "warning", "Excluded MAGMA genes with Z-scores",
                "%s (%.1f%%)" % (
                    f'{compatibility["excluded_target_genes"]:,}',
                    compatibility["excluded_percent"],
                ),
                indent=10, label_width=width,
            ),
            screen_field(
                "warning", "Absent from annotation",
                f'{compatibility["absent_from_annotation"]:,}',
                indent=10, label_width=width,
            ),
            screen_field(
                "warning", "Absent from feature rows",
                f'{compatibility["absent_from_features"]:,}',
                indent=10, label_width=width,
            ),
            screen_field(
                "warning", "Absent from both resources",
                f'{compatibility["absent_from_both"]:,}',
                indent=10, label_width=width,
            ),
            screen_field(
                "info", "Gene-level audit", compatibility["audit_table"],
                indent=10, label_width=width, break_long_values=True,
            ),
            screen_field(
                "info", "Compatibility report", compatibility["report"],
                indent=10, label_width=width, break_long_values=True,
            ),
        ])
    lines.extend(["", screen_line(
        "decision", "Top prioritized genes", indent=6,
    )])
    precision = module.reporting.score_decimal_places
    for gene in summary["top_genes"]:
        label = "%d. %s" % (gene["rank"], gene["gene_name"])
        value = "%s · PoPS score %.*f" % (
            gene["gene_id"], precision, gene["score"],
        )
        lines.append(screen_field(
            "genetic", label, value, indent=10, label_width=width,
        ))
    lines.extend([
        "",
        screen_line("decision", "How to interpret", indent=6),
        screen_field(
            "info", "PoPS significance cutoff",
            "None. PoPS scores are relative rankings, not p-values.",
            indent=10, label_width=width,
        ),
        screen_field(
            "info", "Relevant-gene count",
            "Not statistically defined; report a prespecified top-ranked set and "
            "support it with independent genetic or functional evidence.",
            indent=10, label_width=width,
        ),
        screen_field(
            "info", "Causal interpretation",
            "A high PoPS score prioritizes a gene; it does not establish causality.",
            indent=10, label_width=width,
        ),
    ])
    if summary["warnings"]:
        lines.extend(["", screen_line("warning", "Scientific warnings", indent=6)])
        lines.extend(
            screen_field(
                "warning", "Warning %d" % number, warning,
                indent=10, label_width=width,
            )
            for number, warning in enumerate(summary["warnings"], 1)
        )
    lines.extend([
        "",
        screen_field(
            "success", "Integrated gene-results table",
            paths["integrated_results_file"], indent=6,
            label_width=outer_width, break_long_values=True,
        ),
        screen_field(
            "success", "Interactive gene-results report",
            paths["integrated_report_file"], indent=6,
            label_width=outer_width, break_long_values=True,
        ),
        screen_field(
            "info", "Official PoPS predictions", paths["pops_file"],
            indent=6, label_width=outer_width, break_long_values=True,
        ),
        screen_field(
            "info", "Full log", log_path, indent=6,
            label_width=outer_width, break_long_values=True,
        ),
        "",
    ])
    return "\n".join(lines)


def _restore_root_logging(handlers: list[logging.Handler], level: int) -> None:
    root = logging.getLogger()
    for handler in list(root.handlers):
        if handler not in handlers:
            root.removeHandler(handler)
            handler.close()
    root.handlers[:] = handlers
    root.setLevel(level)


def run_pops_direct(args: argparse.Namespace, ctx=None):
    """Validate inputs, run upstream PoPS in isolation, and publish complete outputs."""
    try:
        configuration = _resolved_configuration(args)
        module = configuration.modules.pops
        output = Path(configuration.run.output_directory).expanduser().resolve()
        dataset = validate_filename_component(
            configuration.run.dataset_id,
            "dataset_id",
            error_type=PopsError,
        )
        output.mkdir(parents=True, exist_ok=True)
        prefix, final_paths = _output_paths(output, dataset, module)
        completion_manifest = configured_output_path(
            output,
            module.output_layout.completion_manifest,
            error_type=PopsError,
            dataset_id=dataset,
        )
        log_path = configured_output_path(
            output,
            module.output_layout.service_log_file,
            error_type=PopsError,
            dataset_id=dataset,
        )
    except BaseException as exc:
        fallback = load_configuration()
        output = Path(
            getattr(args, "output_directory", None) or fallback.run.output_directory
        ).expanduser().resolve()
        raw_dataset = getattr(args, "dataset_id", None) or fallback.run.dataset_id
        try:
            dataset = validate_filename_component(raw_dataset, "dataset_id")
        except ValueError:
            dataset = fallback.run.dataset_id
        log_path = configured_output_path(
            output,
            fallback.modules.pops.output_layout.service_log_file,
            error_type=PopsError,
            dataset_id=dataset,
        )
        write_log_record(
            log_path,
            "ERROR",
            "PoPS configuration or input validation failed: %s: %s"
            % (type(exc).__name__, exc),
            sample_id=dataset,
            file_level=fallback.logging.file_level,
            screen_level=fallback.logging.console_level,
        )
        raise

    pipeline_progress = getattr(args, "_pipeline_stage_progress", None)
    pipeline_stage_numbers = pops_pipeline_stage_numbers(args)
    logger = PipelineLogger(
        dataset,
        "run",
        str(log_path.parent),
        level=configuration.logging.file_level,
        screen_level=configuration.logging.console_level,
        log_path=str(log_path),
        stage_progress=StageProgress(
            "PoPS analysis progress",
            enabled=pipeline_progress is None or pipeline_stage_numbers is None,
            outcome_label_width=configuration.logging.terminal_label_width,
        ),
    )

    def start_pipeline_stage(key: str) -> None:
        if pipeline_progress is not None and pipeline_stage_numbers is not None:
            pipeline_progress.start(int(pipeline_stage_numbers[key]))

    def complete_pipeline_stage(key: str, fields, *, defer: bool = False) -> None:
        if pipeline_progress is None or pipeline_stage_numbers is None:
            return
        if defer:
            args._pipeline_stage_completion = {"outcome_fields": fields}
        else:
            pipeline_progress.complete(
                int(pipeline_stage_numbers[key]), outcome_fields=fields,
            )

    active_upstream_manager = None
    active_upstream_step = None
    upstream_stage_index = 0

    def open_upstream_stage() -> None:
        nonlocal active_upstream_manager, active_upstream_step
        key = _POPS_UPSTREAM_PROGRESS_ORDER[upstream_stage_index]
        local_number = 6 + upstream_stage_index
        start_pipeline_stage(key)
        active_upstream_manager = logger.step(
            local_number,
            len(POPS_STAGES),
            POPS_STAGES[local_number - 1],
            "pops_%s" % key,
        )
        active_upstream_step = active_upstream_manager.__enter__()

    def record_upstream_progress(stage: str, metrics: dict) -> None:
        nonlocal active_upstream_manager, active_upstream_step, upstream_stage_index
        expected_stage = _POPS_UPSTREAM_PROGRESS_ORDER[upstream_stage_index]
        if stage != expected_stage:
            raise PopsError(
                "PoPS progress stages are out of order: expected %s, received %s"
                % (expected_stage, stage)
            )
        if active_upstream_manager is None or active_upstream_step is None:
            raise PopsError("PoPS reported a progress stage before it was started.")
        fields = _progress_outcome_fields(
            stage, metrics, score_source=outcome["label"],
        )
        active_upstream_step.observed("pops_%s" % stage, **metrics)
        active_upstream_step.outcome(
            "PoPS %s completed with validated metrics." % stage.replace("_", " "),
            fields=fields,
            **metrics,
        )
        active_upstream_manager.__exit__(None, None, None)
        active_upstream_manager = None
        active_upstream_step = None
        complete_pipeline_stage(stage, fields)
        upstream_stage_index += 1
        if upstream_stage_index < len(_POPS_UPSTREAM_PROGRESS_ORDER):
            open_upstream_stage()

    try:
        target_alternatives = (
            RequiredAlternative((
                RequiredArgument(
                    "--magma-association-prefix",
                    "modules.pops.magma_association_prefix",
                    module.magma_association_prefix,
                ),
                RequiredArgument(
                    "--target-score-file",
                    "modules.pops.target_score_file",
                    module.target_score_file,
                ),
            )),
        )
        start_pipeline_stage("target_inputs")
        with logger.step(
            1,
            len(POPS_STAGES),
            POPS_STAGES[0],
            "validate_pops_target_inputs",
        ) as step:
            require_resolved_arguments(
                (
                    RequiredArgument(
                        "--genome-build",
                        "modules.pops.genome_build",
                        module.genome_build,
                    ),
                    RequiredArgument(
                        "--feature-matrix-prefix",
                        "modules.pops.feature_matrix_prefix",
                        module.feature_matrix_prefix,
                    ),
                    RequiredArgument(
                        "--pops-gene-location-file",
                        "modules.pops.gene_location_file",
                        module.gene_location_file,
                    ),
                ),
                alternatives=target_alternatives,
            )
            _require_pops_runtime()
            if module.magma_association_prefix is not None:
                outcome = _validate_magma(module)
            elif module.target_score_file is not None:
                outcome = _validate_target(module)
            else:
                raise PopsError(
                    "Provide --magma-association-prefix PREFIX or "
                    "--target-score-file PATH (a custom gene-score table)."
                )
            if module.target_error_covariance_file is not None:
                covariance_path = _required_file(
                    module.target_error_covariance_file,
                    "PoPS custom gene-error covariance",
                )
                _, converted = _load_target_covariance(
                    covariance_path, outcome["gene_count"],
                )
                outcome["covariance"] = {
                    "path": covariance_path,
                    "dimension": outcome["gene_count"],
                    "representation": (
                        "NumPy dense array; converted for PoPS"
                        if converted else "SciPy sparse NPZ"
                    ),
                }
            target_fields = _target_validation_fields(module, outcome)
            step.outcome(
                "Gene-score inputs and optional companions are valid.",
                fields=target_fields,
                target_source=outcome["label"],
                target_genes=outcome["gene_count"],
                target_score_file=str(
                    outcome.get("output_path", outcome.get("path"))
                ),
                magma_annotated_results_file=(
                    str(outcome["annotated"]["path"])
                    if outcome.get("annotated") is not None else None
                ),
                magma_covariance_file=(
                    str(outcome["raw_path"])
                    if outcome.get("raw_path") is not None else None
                ),
                target_covariates_file=(
                    str(outcome["covariates_path"])
                    if outcome.get("covariates_path") is not None else None
                ),
                target_error_covariance_file=(
                    str(outcome["covariance"]["path"])
                    if outcome.get("covariance") is not None else None
                ),
                target_error_covariance_dimension=(
                    outcome["covariance"]["dimension"]
                    if outcome.get("covariance") is not None else None
                ),
                structural_validation="passed",
            )
        complete_pipeline_stage("target_inputs", target_fields)

        start_pipeline_stage("gene_annotation")
        with logger.step(
            2,
            len(POPS_STAGES),
            POPS_STAGES[1],
            "validate_pops_gene_annotation",
        ) as step:
            annotation = _validate_gene_annotation(module)
            annotation_fields = _annotation_validation_fields(module, annotation)
            step.outcome(
                "PoPS gene-location annotation is structurally valid.",
                fields=annotation_fields,
                gene_annotation_file=str(annotation["path"]),
                genome_build=module.genome_build.value,
                annotated_genes=annotation["gene_count"],
                chromosomes=annotation["chromosomes"],
                gene_id_column=module.input_schema.gene_annotation_id_column,
                chromosome_column=(
                    module.input_schema.gene_annotation_chromosome_column
                ),
                transcription_start_site_column=(
                    module.input_schema.gene_annotation_tss_column
                ),
                tss_minimum=annotation["tss_minimum"],
                tss_maximum=annotation["tss_maximum"],
                structural_validation="passed",
            )
        complete_pipeline_stage("gene_annotation", annotation_fields)

        start_pipeline_stage("feature_matrix")
        with logger.step(
            3,
            len(POPS_STAGES),
            POPS_STAGES[2],
            "validate_pops_feature_matrices",
        ) as step:
            features = _validate_feature_resources(module)
            feature_fields = _feature_validation_fields(module, features)
            for metrics in features["chunks"]:
                step.observed(
                    "pops_feature_matrix_chunk",
                    chunk=metrics["chunk"],
                    columns_file=str(metrics["columns_path"]),
                    matrix_file=str(metrics["matrix_path"]),
                    rows=metrics["rows"],
                    columns=metrics["columns"],
                    dtype=metrics["dtype"],
                    finite_numeric_values=True,
                )
            step.outcome(
                "All PoPS feature-matrix companions, dimensions, and values are valid.",
                fields=feature_fields,
                feature_genes=features["row_count"],
                feature_count=features["feature_count"],
                chunks=len(features["chunks"]),
                feature_rows_file=str(features["rows_path"]),
                structural_validation="passed",
                finite_numeric_values=True,
            )
        complete_pipeline_stage("feature_matrix", feature_fields)

        start_pipeline_stage("feature_controls")
        with logger.step(
            4,
            len(POPS_STAGES),
            POPS_STAGES[3],
            "validate_pops_feature_controls",
        ) as step:
            controls = _validate_feature_control_resources(module, features)
            features["control_resources"] = controls
            control_fields = _feature_control_validation_fields(controls)
            step.outcome(
                "Optional PoPS feature-control resources are valid.",
                fields=control_fields,
                feature_subset_count=(
                    controls["feature_subset"]["count"]
                    if controls["feature_subset"] is not None else 0
                ),
                control_feature_count=(
                    controls["controls"]["count"]
                    if controls["controls"] is not None else 0
                ),
                feature_subset_file=(
                    str(controls["feature_subset"]["path"])
                    if controls["feature_subset"] is not None else None
                ),
                control_features_file=(
                    str(controls["controls"]["path"])
                    if controls["controls"] is not None else None
                ),
                matrix_feature_compatibility="passed",
                controls_retained_by_feature_subset=(
                    controls["controls_retained_by_subset"]
                ),
            )
        complete_pipeline_stage("feature_controls", control_fields)

        start_pipeline_stage("gene_compatibility")
        with logger.step(
            5,
            len(POPS_STAGES),
            POPS_STAGES[4],
            "validate_pops_gene_universes",
        ) as step:
            if pipeline_progress is not None and pipeline_stage_numbers is not None:
                vcf_evidence = getattr(args, "_pipeline_vcf_header_evidence", {})
                magma_build = (
                    vcf_evidence.get("genome_build")
                    or getattr(args, "genome_build", None)
                    or configuration.modules.magma.genome_build
                )
                magma_build_value = getattr(magma_build, "value", str(magma_build))
                if module.genome_build.value != magma_build_value:
                    raise PopsError(
                        "PoPS genome build %s does not match pipeline MAGMA "
                        "genome build %s."
                        % (module.genome_build.value, magma_build_value)
                    )
                outcome["pipeline_genome_build"] = magma_build_value
            _validate_gene_universes(module, features, annotation, outcome)
            compatibility_fields = _gene_compatibility_fields(
                module, features, annotation, outcome,
            )
            step.outcome(
                "Input-score, annotation, and feature gene identifiers are "
                "compatible.",
                fields=compatibility_fields,
                **{
                    key: value
                    for key, value in outcome["compatibility"].items()
                    if not key.endswith("_gene_ids")
                },
                annotation_genes=annotation["gene_count"],
                feature_genes=features["row_count"],
                declared_genome_build=module.genome_build.value,
                chromosome_consistency="passed",
            )
        complete_pipeline_stage("gene_compatibility", compatibility_fields)

        logger.record(
            "PARAM", "pops_run",
            genome_build=module.genome_build.value,
            method=module.method,
            random_seed=configuration.execution.random_seed,
            feature_matrix_chunks=module.feature_matrix_chunks,
            feature_selection_p_cutoff=module.feature_selection_p_cutoff,
            minimum_gene_count=module.minimum_gene_count,
            gene_universe_policy=module.gene_universe_policy,
        )
        logger.record(
            "OBSERVED", "pops_inputs",
            annotation_genes=annotation["gene_count"],
            outcome_genes=outcome["gene_count"],
            feature_genes=features["row_count"],
            features=features["feature_count"],
        )
        compatibility = outcome["compatibility"]
        logger.record(
            "OBSERVED", "pops_gene_compatibility",
            policy=compatibility["policy"],
            original_target_genes=compatibility["original_target_genes"],
            retained_target_genes=compatibility["retained_target_genes"],
            excluded_target_genes=compatibility["excluded_target_genes"],
            retained_percent=compatibility["retained_percent"],
            absent_from_annotation=compatibility["absent_from_annotation"],
            absent_from_features=compatibility["absent_from_features"],
            absent_from_both=compatibility["absent_from_both"],
            absent_only_from_annotation=(
                compatibility["absent_only_from_annotation"]
            ),
            absent_only_from_features=compatibility["absent_only_from_features"],
        )
        resolved_path = configured_output_path(
            output,
            module.output_layout.resolved_config_file,
            error_type=PopsError,
            dataset_id=dataset,
        )
        write_resolved_configuration(configuration, resolved_path, modules="pops")
        logger.record("OUTPUT", "resolved_configuration", path=str(resolved_path))

        if (
            _complete_outputs(final_paths)
            and completion_manifest.is_file()
            and completion_manifest.stat().st_size > 0
            and configuration.run.resume
            and not configuration.run.overwrite
        ):
            result = _result(output, prefix, final_paths)
            result["completion_manifest"] = str(completion_manifest)
            summary = _summarise_outputs(
                final_paths, module, annotation, outcome, features,
            )
            result["summary"] = summary
            logger.record("SKIP", "pops_run", reason="validated_complete_outputs")
            for index, key in enumerate(_POPS_UPSTREAM_PROGRESS_ORDER, 6):
                start_pipeline_stage(key)
                with logger.step(
                    index,
                    len(POPS_STAGES),
                    POPS_STAGES[index - 1],
                    "resume_validated_pops_stage",
                ) as step:
                    resumed_fields = [
                        ("success", "Resume decision", "reused validated complete outputs"),
                    ]
                    step.outcome(
                        "Reused completed PoPS analysis stage.",
                        fields=resumed_fields,
                        resume_mode="validated_complete_outputs",
                    )
                complete_pipeline_stage(key, resumed_fields)
            start_pipeline_stage("results")
            with logger.step(
                11,
                len(POPS_STAGES),
                POPS_STAGES[10],
                "validate_resumed_pops_results",
            ) as step:
                result_fields = [
                    ("success", "Resume decision", "reused validated complete outputs"),
                    ("count", "Genes receiving PoPS scores", summary["genes_scored"]),
                    ("count", "Integrated gene rows", summary["integrated_gene_count"]),
                    ("count", "Published result files", len(final_paths)),
                ]
                step.outcome(
                    "Existing PoPS outputs were validated and reused.",
                    fields=result_fields,
                    genes_scored=summary["genes_scored"],
                )
            complete_pipeline_stage("results", result_fields, defer=True)
            print(_render_summary(
                summary, dataset, module, final_paths, log_path,
                configuration.logging.terminal_label_width,
            ))
            if ctx is not None:
                ctx["pops_output"] = result["pops_file"]
                ctx["pops_integrated_results"] = result["integrated_results_file"]
                ctx["pops_integrated_report"] = result["integrated_report_file"]
            return result
        existing = [path for path in final_paths.values() if path.exists()]
        if existing and not configuration.run.overwrite:
            raise PopsError(
                "Existing PoPS outputs are incomplete or resume is disabled. "
                "Use --overwrite after reviewing: %s"
                % ", ".join(str(path) for path in existing)
            )

        staging_root = configured_output_path(
            output,
            module.output_layout.staging_directory,
            error_type=PopsError,
            dataset_id=dataset,
        )
        staging_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="run_", dir=staging_root) as directory:
            staging_prefix = Path(directory) / prefix.name
            staged_paths = {
                name: Path(str(staging_prefix) + path.name.removeprefix(prefix.name))
                for name, path in final_paths.items()
            }
            effective_magma_prefix = None
            if module.gene_universe_policy == "intersect":
                effective_magma_prefix, compatibility_report = (
                    _prepare_intersected_magma(
                        module, annotation, features, outcome, staged_paths,
                        final_paths,
                    )
                )
                logger.record(
                    "ACTION", "pops_gene_universe_intersection",
                    **{
                        key: value
                        for key, value in compatibility_report.items()
                        if key not in {"chromosomes", "files"}
                    },
                )
                for chromosome, counts in compatibility_report["chromosomes"].items():
                    logger.record(
                        "OBSERVED", "pops_gene_compatibility_chromosome",
                        chromosome=chromosome, **counts,
                    )
                logger.record(
                    "OUTPUT", "pops_gene_compatibility_files",
                    **compatibility_report["files"],
                )
                if compatibility_report["excluded_target_genes"]:
                    logger.warning(
                        "PoPS gene-universe intersection retained %d of %d MAGMA "
                        "genes with Z-scores (%.1f%%) and excluded %d; originals "
                        "were "
                        "not modified."
                        % (
                            compatibility_report["retained_target_genes"],
                            compatibility_report["original_target_genes"],
                            compatibility_report["retained_percent"],
                            compatibility_report["excluded_target_genes"],
                        )
                    )
            covariance_path = None
            if module.target_error_covariance_file is not None:
                source = _required_file(
                    module.target_error_covariance_file,
                    "PoPS custom gene-error covariance",
                )
                covariance, converted = _load_target_covariance(source, outcome["gene_count"])
                if converted:
                    covariance_path = Path(directory) / "target_covariance.npz"
                    sparse.save_npz(covariance_path, sparse.csr_matrix(covariance))
                    logger.record(
                        "ACTION", "target_covariance",
                        source=str(source), representation="scipy_npz_for_upstream_compatibility",
                    )
            upstream_arguments = _build_upstream_arguments(
                module,
                configuration.execution.random_seed,
                staging_prefix,
                covariance_path=covariance_path,
                magma_prefix=effective_magma_prefix,
            )
            upstream_handlers = list(logging.getLogger().handlers)
            upstream_level = logging.getLogger().level
            try:
                get_pops_args, pops_main = _load_upstream_entrypoints()
                open_upstream_stage()
                pops_main(
                    vars(get_pops_args(upstream_arguments)),
                    progress_callback=record_upstream_progress,
                )
                if upstream_stage_index != len(_POPS_UPSTREAM_PROGRESS_ORDER):
                    raise PopsError(
                        "Upstream PoPS did not report every required progress stage."
                    )
            except SystemExit as exc:
                raise PopsError("Upstream PoPS exited before completion.") from exc
            finally:
                _restore_root_logging(upstream_handlers, upstream_level)

            start_pipeline_stage("results")
            with logger.step(
                11,
                len(POPS_STAGES),
                POPS_STAGES[10],
                "validate_and_publish_pops_results",
            ) as step:
                derived_names = {
                    "integrated_results_file", "integrated_report_file",
                }
                missing = [
                    path for name, path in staged_paths.items()
                    if name not in derived_names
                    if not path.is_file() or path.stat().st_size <= 0
                ]
                if missing:
                    raise PopsError(
                        "PoPS staging did not create required non-empty outputs: %s"
                        % ", ".join(str(path) for path in missing)
                    )
                integrated = _write_integrated_results(
                    staged_paths,
                    module,
                    annotation,
                    outcome,
                    features,
                    dataset=dataset,
                )
                if not _complete_outputs(staged_paths):
                    raise PopsError(
                        "PoPS staging did not create every required non-empty output."
                    )
                summary = _summarise_outputs(
                    staged_paths, module, annotation, outcome, features,
                    published_paths=final_paths,
                    integrated=integrated,
                )
                all_known_suffixes = {
                    module.output_layout.predictions_suffix,
                    module.output_layout.coefficients_suffix,
                    module.output_layout.marginals_suffix,
                    module.output_layout.upstream_log_suffix,
                    module.output_layout.training_data_suffix,
                    module.output_layout.matrix_data_suffix,
                    module.output_layout.compatible_genes_out_suffix,
                    module.output_layout.compatible_genes_raw_suffix,
                    module.output_layout.excluded_genes_out_suffix,
                    module.output_layout.excluded_genes_raw_suffix,
                    module.output_layout.gene_compatibility_table_suffix,
                    module.output_layout.gene_compatibility_report_suffix,
                    module.output_layout.integrated_results_suffix,
                    module.output_layout.integrated_report_suffix,
                }
                if configuration.run.overwrite:
                    completion_manifest.unlink(missing_ok=True)
                    for suffix in all_known_suffixes:
                        Path(str(prefix) + suffix).unlink(missing_ok=True)
                for name, source in staged_paths.items():
                    destination = final_paths[name]
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    source.replace(destination)
                result = _result(output, prefix, final_paths)
                result["summary"] = summary
                write_yaml_report(
                    {
                        "status": "COMPLETED",
                        "dataset_id": dataset,
                        "genome_build": module.genome_build.value,
                        "method": module.method,
                        "summary": summary,
                        "files": result["published_files"],
                    },
                    completion_manifest,
                )
                result["completion_manifest"] = str(completion_manifest)
                logger.record(
                    "OUTPUT", "pops_outputs", files=result["published_files"],
                )
                logger.record(
                    "OUTPUT", "pops_integrated_gene_results",
                    table=result["integrated_results_file"],
                    html=result["integrated_report_file"],
                    full_union_genes=summary["integrated_gene_count"],
                    genes_scored=summary["genes_scored"],
                    genes_scored_without_input_target=(
                        summary["genes_scored_without_target"]
                    ),
                    input_target_genes_excluded=(
                        summary["input_target_genes_excluded"]
                    ),
                )
                logger.record("OBSERVED", "pops_summary", **{
                    key: value
                    for key, value in summary.items()
                    if key != "top_genes"
                })
                for warning in summary["warnings"]:
                    logger.warning(warning)
                result_fields = [
                    ("success", "Output validation", "all required files are non-empty"),
                    ("count", "Genes receiving PoPS scores", summary["genes_scored"]),
                    ("count", "Integrated gene rows", summary["integrated_gene_count"]),
                    (
                        "count",
                        "Selected features",
                        summary["selected_features"]
                        if summary["selected_features"] is not None
                        else "not recorded",
                    ),
                    (
                        "warning" if summary["warnings"] else "success",
                        "Scientific warnings",
                        len(summary["warnings"]),
                    ),
                    ("count", "Published result files", len(final_paths)),
                ]
                step.outcome(
                    "PoPS results were validated and published.",
                    fields=result_fields,
                    genes_scored=summary["genes_scored"],
                    published_files=len(final_paths),
                )
            complete_pipeline_stage("results", result_fields, defer=True)

        logger.record("STATUS", "pops_run", status="COMPLETED")
        if ctx is not None:
            ctx["pops_output"] = result["pops_file"]
            ctx["pops_integrated_results"] = result["integrated_results_file"]
            ctx["pops_integrated_report"] = result["integrated_report_file"]
        print(_render_summary(
            summary, dataset, module, final_paths, log_path,
            configuration.logging.terminal_label_width,
        ))
        return result
    except BaseException as exc:
        if active_upstream_manager is not None:
            active_upstream_manager.__exit__(type(exc), exc, exc.__traceback__)
        if not logger.summary()["failed"]:
            logger.error("PoPS analysis failed: %s: %s" % (type(exc).__name__, exc))
        raise
    finally:
        logger.close()


__all__ = [
    "preflight_pops_pipeline",
    "run_pops_direct",
    "validate_pops_configuration",
]
