"""Validated PostGWAS service boundary for the upstream PoPS v0.2 program."""

from __future__ import annotations

import argparse
from importlib.util import find_spec
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
from postgwas.core.paths import configured_output_path, validate_filename_component
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.core.io.reports import write_yaml_report
from postgwas.core.ui.progress import StageProgress
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.modules.pops.errors import PopsError


_POPS_PROGRESS_STAGES = (
    "Validate PoPS inputs and reference resources",
    "Load target gene scores",
    "Adjust target scores for configured covariates",
    "Test and select predictive features",
    "Fit the PoPS prediction model",
    "Calculate genome-wide PoPS scores",
    "Validate and publish PoPS results",
)
_POPS_UPSTREAM_PROGRESS_ORDER = (
    "target_loading",
    "covariate_adjustment",
    "feature_selection",
    "model_fitting",
    "gene_scoring",
)


def _progress_outcome_fields(stage: str, metrics: dict) -> list[tuple[str, str, object]]:
    """Translate upstream stage metrics into stable, scientific screen fields."""
    if stage == "target_loading":
        return [
            ("genetic", "Target source", metrics["target_source"]),
            ("count", "Target genes loaded", metrics["target_genes"]),
            ("count", "Target covariates", metrics["covariates"]),
            ("analysis", "Error covariance", metrics["error_covariance"]),
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
            ("count", "Training genes", metrics["training_genes"]),
            ("count", "Model features", metrics["model_features"]),
        ]
        if metrics.get("selected_cv_alpha") is not None:
            fields.append(
                ("analysis", "Selected regularisation alpha", metrics["selected_cv_alpha"])
            )
        return fields
    if stage == "gene_scoring":
        return [
            ("count", "Genes assigned PoPS scores", metrics["genes_scored"]),
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
        "feature_matrix_prefix": "feature_matrix_prefix",
        "feature_matrix_chunks": "feature_matrix_chunks",
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
    if value is None:
        raise PopsError("%s is required." % label)
    path = Path(value).expanduser().resolve()
    if not path.is_file() or path.stat().st_size <= 0:
        raise PopsError("%s does not exist or is empty: %s" % (label, path))
    return path


def _prefix_file(prefix: str | None, suffix: str, label: str) -> Path:
    if prefix is None:
        raise PopsError("%s prefix is required." % label)
    return _required_file(str(Path(prefix).expanduser()) + suffix, label)


def _read_names(path: Path, label: str) -> np.ndarray:
    values = np.atleast_1d(np.loadtxt(path, dtype=str)).reshape(-1)
    if values.size == 0 or any(not str(value).strip() for value in values):
        raise PopsError("%s contains no usable names: %s" % (label, path))
    if len(values) != len(set(values.tolist())):
        raise PopsError("%s contains duplicate names: %s" % (label, path))
    return values


def _read_table(path: Path, delimiter: str, label: str) -> pd.DataFrame:
    try:
        table = pd.read_csv(path, sep=delimiter)
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise PopsError("Cannot read %s %s: %s" % (label, path, exc)) from exc
    if table.empty:
        raise PopsError("%s contains no data rows: %s" % (label, path))
    return table


def _require_columns(table: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise PopsError(
            "%s is missing required columns: %s"
            % (label, ", ".join(missing))
        )


def _validate_feature_resources(module) -> dict:
    schema = module.input_schema
    prefix = module.feature_matrix_prefix
    rows_path = _prefix_file(prefix, schema.feature_rows_suffix, "PoPS feature rows")
    rows = _read_names(rows_path, "PoPS feature rows")
    all_columns: list[str] = []
    chunk_metrics = []
    for chunk in range(module.feature_matrix_chunks):
        columns_path = _prefix_file(
            prefix,
            schema.feature_columns_pattern.format(chunk=chunk),
            "PoPS feature columns chunk %d" % chunk,
        )
        matrix_path = _prefix_file(
            prefix,
            schema.feature_matrix_pattern.format(chunk=chunk),
            "PoPS feature matrix chunk %d" % chunk,
        )
        columns = _read_names(columns_path, "PoPS feature columns chunk %d" % chunk)
        try:
            matrix = np.load(matrix_path, mmap_mode="r", allow_pickle=False)
        except (OSError, ValueError) as exc:
            raise PopsError("Cannot read feature matrix %s: %s" % (matrix_path, exc)) from exc
        if matrix.ndim != 2 or matrix.shape != (len(rows), len(columns)):
            raise PopsError(
                "Feature matrix chunk %d has shape %s; expected (%d, %d)."
                % (chunk, matrix.shape, len(rows), len(columns))
            )
        all_columns.extend(columns.tolist())
        chunk_metrics.append({"chunk": chunk, "rows": matrix.shape[0], "columns": matrix.shape[1]})
    if len(all_columns) != len(set(all_columns)):
        raise PopsError("Feature names must be unique across all matrix chunks.")
    for value, label in (
        (module.feature_subset_file, "PoPS feature-subset file"),
        (module.control_features_file, "PoPS control-features file"),
    ):
        if value is not None:
            names = _read_names(_required_file(value, label), label)
            unknown = sorted(set(names.tolist()) - set(all_columns))
            if unknown:
                raise PopsError(
                    "%s contains features absent from the matrix: %s"
                    % (label, ", ".join(unknown[:10]))
                )
    return {
        "rows": rows,
        "row_count": len(rows),
        "feature_count": len(all_columns),
        "chunks": chunk_metrics,
    }


def _validate_gene_annotation(module) -> dict:
    schema = module.input_schema
    path = _required_file(module.gene_location_file, "PoPS gene-location file")
    table = _read_table(path, schema.table_delimiter_pattern, "PoPS gene annotation")
    required = [
        schema.gene_annotation_id_column,
        schema.gene_annotation_chromosome_column,
        schema.gene_annotation_tss_column,
    ]
    _require_columns(table, required, "PoPS gene annotation")
    identifiers = table[schema.gene_annotation_id_column].astype(str)
    if identifiers.duplicated().any():
        raise PopsError("PoPS gene annotation contains duplicate gene identifiers.")
    chromosome = table[schema.gene_annotation_chromosome_column].astype(str)
    if chromosome.str.strip().eq("").any():
        raise PopsError("PoPS gene annotation contains empty chromosome labels.")
    tss = pd.to_numeric(table[schema.gene_annotation_tss_column], errors="coerce")
    if tss.isna().any() or not np.isfinite(tss.to_numpy()).all() or (tss < 0).any():
        raise PopsError("PoPS gene annotation TSS values must be finite and non-negative.")
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
    unknown = sorted(configured_chromosomes - set(chromosome))
    if unknown:
        raise PopsError(
            "Configured PoPS chromosomes are absent from the gene annotation: %s"
            % ", ".join(unknown)
        )
    names = identifiers
    if schema.gene_annotation_name_column in table.columns:
        configured_names = table[schema.gene_annotation_name_column].astype(str)
        names = configured_names.where(configured_names.str.strip().ne(""), identifiers)
    return {
        "path": path,
        "genes": set(identifiers),
        "gene_count": len(identifiers),
        "chromosomes": sorted(set(chromosome)),
        "gene_names": dict(zip(identifiers, names)),
        "gene_chromosomes": dict(zip(identifiers, chromosome)),
    }


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
    table = _read_table(output_path, schema.table_delimiter_pattern, "MAGMA gene results")
    _require_columns(
        table,
        [schema.magma_gene_id_column, schema.magma_score_column],
        "MAGMA gene results",
    )
    identifiers = table[schema.magma_gene_id_column].astype(str)
    scores = pd.to_numeric(table[schema.magma_score_column], errors="coerce")
    if identifiers.duplicated().any():
        raise PopsError("MAGMA gene results contain duplicate gene identifiers.")
    if scores.isna().any() or not np.isfinite(scores.to_numpy()).all():
        raise PopsError("MAGMA gene Z statistics must all be finite.")
    return {
        "genes": set(identifiers),
        "gene_count": len(identifiers),
        "output_path": output_path,
        "raw_path": raw_path,
    }


def _validate_target(module) -> dict:
    schema = module.input_schema
    path = _required_file(module.target_score_file, "PoPS target-score file")
    table = _read_table(
        path, schema.target_table_delimiter, "PoPS target scores",
    )
    _require_columns(
        table,
        [schema.target_gene_id_column, schema.target_score_column],
        "PoPS target scores",
    )
    identifiers = table[schema.target_gene_id_column].astype(str)
    scores = pd.to_numeric(table[schema.target_score_column], errors="coerce")
    if identifiers.duplicated().any():
        raise PopsError("PoPS target scores contain duplicate gene identifiers.")
    if scores.isna().any() or not np.isfinite(scores.to_numpy()).all():
        raise PopsError("PoPS target scores must all be finite.")
    if module.target_covariates_file is not None:
        covariate_path = _required_file(
            module.target_covariates_file, "PoPS target-covariates file",
        )
        covariates = _read_table(
            covariate_path,
            schema.target_table_delimiter,
            "PoPS target covariates",
        )
        _require_columns(
            covariates,
            [schema.target_gene_id_column],
            "PoPS target covariates",
        )
        covariate_ids = covariates[schema.target_gene_id_column].astype(str)
        if covariate_ids.duplicated().any() or set(covariate_ids) != set(identifiers):
            raise PopsError(
                "Target covariates must contain each target-score gene exactly once."
            )
        numeric = covariates.drop(columns=[schema.target_gene_id_column]).apply(
            pd.to_numeric, errors="coerce",
        )
        if (
            numeric.shape[1] == 0
            or numeric.isna().any().any()
            or not np.isfinite(numeric.to_numpy()).all()
        ):
            raise PopsError("Target covariates must contain finite numeric columns.")
    return {"genes": set(identifiers), "gene_count": len(identifiers), "path": path}


def _load_target_covariance(path: Path, expected_size: int) -> tuple[np.ndarray, bool]:
    converted = False
    try:
        covariance = sparse.load_npz(path).toarray()
    except (AttributeError, OSError, TypeError, ValueError):
        try:
            covariance = np.load(path, allow_pickle=False)
            converted = True
        except (OSError, TypeError, ValueError) as exc:
            raise PopsError("Cannot read target covariance %s: %s" % (path, exc)) from exc
    if covariance.shape != (expected_size, expected_size):
        raise PopsError(
            "Target covariance has shape %s; expected (%d, %d)."
            % (covariance.shape, expected_size, expected_size)
        )
    if not np.isfinite(covariance).all() or not np.allclose(covariance, covariance.T):
        raise PopsError("Target covariance must be finite and symmetric.")
    try:
        np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise PopsError("Target covariance must be positive definite.") from exc
    return covariance, converted


def _resource_signature(module) -> tuple:
    return (
        module.feature_matrix_prefix,
        module.feature_matrix_chunks,
        module.gene_location_file,
        module.feature_subset_file,
        module.control_features_file,
        tuple(module.covariate_projection_chromosomes or ()),
        tuple(module.feature_selection_chromosomes or ()),
        tuple(module.training_chromosomes or ()),
    )


def _validate_resolved_pops_configuration(
    args: argparse.Namespace, configuration, *, pipeline: bool = False,
):
    """Preflight resources for one already resolved PoPS configuration."""
    _require_pops_runtime()
    module = configuration.modules.pops
    if module.genome_build is None:
        raise PopsError(
            "Declare the shared PoPS genome build with --genome-build or "
            "modules.pops.genome_build; PostGWAS will not infer it."
        )
    signature = _resource_signature(module)
    cached = getattr(args, "_pops_resource_preflight", None)
    if cached is not None and cached[0] == signature:
        features, annotation = cached[1], cached[2]
    else:
        features = _validate_feature_resources(module)
        annotation = _validate_gene_annotation(module)
    reference_common = annotation["genes"] & set(features["rows"])
    if len(reference_common) < module.minimum_gene_count:
        raise PopsError(
            "Only %d genes are shared by the PoPS annotation and feature rows; "
            "the configured minimum is %d."
            % (len(reference_common), module.minimum_gene_count)
        )
    if pipeline:
        magma_build = configuration.modules.magma.genome_build
        if module.genome_build != magma_build:
            raise PopsError(
                "PoPS genome build %s does not match pipeline MAGMA genome build %s."
                % (module.genome_build.value, magma_build.value)
            )
        outcome = None
    elif module.magma_association_prefix is not None:
        outcome = _validate_magma(module)
    elif module.target_score_file is not None:
        outcome = _validate_target(module)
    else:
        raise PopsError(
            "Provide a MAGMA association prefix or a custom target-score file."
        )
    if outcome is not None:
        common = outcome["genes"] & annotation["genes"] & set(features["rows"])
        if len(common) < module.minimum_gene_count:
            raise PopsError(
                "Only %d genes are shared by outcomes, annotation and features; "
                "the configured minimum is %d."
                % (len(common), module.minimum_gene_count)
            )
        missing_annotation = outcome["genes"] - annotation["genes"]
        missing_features = outcome["genes"] - set(features["rows"])
        if missing_annotation or missing_features:
            raise PopsError(
                "Every outcome gene must be present in both the PoPS annotation "
                "and feature rows (missing annotation=%d, missing features=%d)."
                % (len(missing_annotation), len(missing_features))
            )
        if module.target_error_covariance_file is not None:
            covariance_path = _required_file(
                module.target_error_covariance_file, "PoPS target covariance",
            )
            _load_target_covariance(covariance_path, outcome["gene_count"])
    return configuration, features, annotation, outcome


def validate_pops_configuration(args: argparse.Namespace, *, pipeline: bool = False):
    """Resolve and preflight PoPS resources without running upstream PoPS."""
    return _validate_resolved_pops_configuration(
        args, _resolved_configuration(args), pipeline=pipeline,
    )


def preflight_pops_pipeline(args: argparse.Namespace) -> None:
    """Fail before upstream pipeline steps if configured PoPS resources are invalid."""
    if hasattr(args, "_pops_resource_preflight"):
        del args._pops_resource_preflight
    configuration, features, annotation, _ = validate_pops_configuration(
        args, pipeline=True,
    )
    args._pops_resource_preflight = (
        _resource_signature(configuration.modules.pops), features, annotation,
    )


def _flag(arguments: list[str], condition: bool, enabled: str, disabled: str) -> None:
    arguments.append(enabled if condition else disabled)


def _build_upstream_arguments(module, seed: int, output_prefix: Path, covariance_path=None):
    arguments = [
        "--gene_annot_path", module.gene_location_file,
        "--feature_mat_prefix", module.feature_matrix_prefix,
        "--num_feature_chunks", str(module.feature_matrix_chunks),
        "--out_prefix", str(output_prefix),
    ]
    if module.magma_association_prefix is not None:
        arguments += ["--magma_prefix", module.magma_association_prefix]
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
    }
    if module.save_matrix_files:
        suffixes["training_data_file"] = layout.training_data_suffix
        suffixes["matrix_data_file"] = layout.matrix_data_suffix
    return prefix, {name: Path(str(prefix) + suffix) for name, suffix in suffixes.items()}


def _complete_outputs(paths: dict[str, Path]) -> bool:
    return all(path.is_file() and path.stat().st_size > 0 for path in paths.values())


def _result(output: Path, prefix: Path, paths: dict[str, Path]) -> dict:
    return {
        "status": "success",
        "pops_file": str(paths["pops_file"]),
        "output_dir": str(output),
        "out_prefix": str(prefix),
        "published_files": [str(path) for path in paths.values()],
    }


def _true_count(values: pd.Series) -> int:
    """Count upstream boolean values without treating the string 'False' as true."""
    return int(values.astype(str).str.strip().str.lower().eq("true").sum())


def _summarise_outputs(paths, module, annotation, outcome, features) -> dict:
    """Validate published scientific results and build the terminal summary."""
    schema = module.input_schema
    predictions = _read_table(
        paths["pops_file"], schema.target_table_delimiter, "PoPS predictions",
    )
    gene_column = schema.gene_annotation_id_column
    score_column = schema.prediction_score_column
    _require_columns(predictions, [gene_column, score_column], "PoPS predictions")
    gene_ids = predictions[gene_column].astype(str)
    scores = pd.to_numeric(predictions[score_column], errors="coerce")
    if gene_ids.duplicated().any():
        raise PopsError("PoPS predictions contain duplicate gene identifiers.")
    if scores.isna().any() or not np.isfinite(scores.to_numpy()).all():
        raise PopsError("PoPS prediction scores must all be finite.")
    missing_annotation = set(gene_ids) - annotation["genes"]
    missing_features = set(gene_ids) - set(features["rows"])
    if missing_annotation or missing_features:
        raise PopsError(
            "Every PoPS prediction gene must occur in both annotation and features "
            "(missing annotation=%d, missing features=%d)."
            % (len(missing_annotation), len(missing_features))
        )

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
    training_count = outcome["gene_count"]
    if schema.prediction_training_column in predictions.columns:
        training_count = _true_count(predictions[schema.prediction_training_column])
    selected_features = None
    marginals = _read_table(
        paths["marginals_file"], schema.target_table_delimiter, "PoPS marginals",
    )
    if schema.marginal_selected_column in marginals.columns:
        selected_features = _true_count(marginals[schema.marginal_selected_column])

    reference_genes = annotation["genes"] & set(features["rows"])
    target_genes = outcome["genes"] & reference_genes
    missing_chromosomes = sorted(
        set(annotation["gene_chromosomes"].values())
        - {
            annotation["gene_chromosomes"][gene]
            for gene in target_genes
            if gene in annotation["gene_chromosomes"]
        }
    )
    warnings = []
    if len(target_genes) < len(reference_genes):
        warnings.append(
            "Target scores cover %d of %d compatible genes (%.1f%%); rankings "
            "were trained on incomplete genome-wide gene outcomes."
            % (
                len(target_genes), len(reference_genes),
                100 * len(target_genes) / len(reference_genes),
            )
        )
    if missing_chromosomes:
        warnings.append(
            "No target-scored genes were available on chromosome(s): %s."
            % ", ".join(missing_chromosomes)
        )
    return {
        "genes_scored": len(predictions),
        "target_genes": len(target_genes),
        "compatible_genes": len(reference_genes),
        "training_genes": training_count,
        "selected_features": selected_features,
        "top_genes": top_genes,
        "warnings": warnings,
    }


def _render_summary(
    summary: dict, dataset: str, module, paths, log_path: Path, label_width: int,
) -> str:
    width = label_width
    status = "COMPLETED WITH SCIENTIFIC WARNINGS" if summary["warnings"] else "COMPLETED"
    lines = [
        "",
        screen_line("analysis", "PoPS gene-prioritisation summary", indent=2),
        screen_field("info", "Dataset", dataset, indent=6, label_width=width),
        screen_field(
            "warning" if summary["warnings"] else "success",
            "Analysis status", status, indent=6, label_width=width,
        ),
        "",
        screen_line("genetic", "Scientific findings", indent=6),
        screen_field(
            "count", "Genes assigned PoPS scores",
            f'{summary["genes_scored"]:,}', indent=10, label_width=width,
        ),
        screen_field(
            "count", "Genes with target scores",
            "%s of %s" % (
                f'{summary["target_genes"]:,}', f'{summary["compatible_genes"]:,}',
            ),
            indent=10, label_width=width,
        ),
        screen_field(
            "count", "Genes used for training", f'{summary["training_genes"]:,}',
            indent=10, label_width=width,
        ),
    ]
    if summary["selected_features"] is not None:
        lines.append(screen_field(
            "count", "Features selected", f'{summary["selected_features"]:,}',
            indent=10, label_width=width,
        ))
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
            "success", "Complete PoPS results", paths["pops_file"],
            indent=6, label_width=width,
        ),
        screen_field("info", "Full log", log_path, indent=6, label_width=width),
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
    progress = None
    active_progress_step = 0
    try:
        configuration = _resolved_configuration(args)
        progress = StageProgress(
            "PoPS analysis progress",
            enabled=configuration.logging.show_progress,
            outcome_label_width=configuration.logging.terminal_label_width,
        )
        active_progress_step = 1
        progress.start_step(
            active_progress_step,
            len(_POPS_PROGRESS_STAGES),
            _POPS_PROGRESS_STAGES[active_progress_step - 1],
        )
        configuration, features, annotation, outcome = (
            _validate_resolved_pops_configuration(args, configuration)
        )
        compatible_target_genes = len(
            outcome["genes"] & annotation["genes"] & set(features["rows"])
        )
        progress.complete_step(
            active_progress_step,
            len(_POPS_PROGRESS_STAGES),
            _POPS_PROGRESS_STAGES[active_progress_step - 1],
            outcome_fields=[
                ("genetic", "Declared genome build", configuration.modules.pops.genome_build.value),
                ("count", "Annotation genes", annotation["gene_count"]),
                ("count", "Feature-matrix genes", features["row_count"]),
                ("count", "Available features", features["feature_count"]),
                ("success", "Compatible target genes", compatible_target_genes),
            ],
        )
    except BaseException as exc:
        if progress is not None and active_progress_step:
            progress.fail_step(
                active_progress_step,
                len(_POPS_PROGRESS_STAGES),
                _POPS_PROGRESS_STAGES[active_progress_step - 1],
            )
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
    module = configuration.modules.pops
    output = Path(configuration.run.output_directory).expanduser().resolve()
    dataset = configuration.run.dataset_id
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
    logger = PipelineLogger(
        dataset,
        "run",
        str(log_path.parent),
        level=configuration.logging.file_level,
        screen_level=configuration.logging.console_level,
        log_path=str(log_path),
    )

    def record_upstream_progress(stage: str, metrics: dict) -> None:
        nonlocal active_progress_step
        expected_stage = _POPS_UPSTREAM_PROGRESS_ORDER[active_progress_step - 2]
        if stage != expected_stage:
            raise PopsError(
                "PoPS progress stages are out of order: expected %s, received %s"
                % (expected_stage, stage)
            )
        logger.record("OBSERVED", "pops_%s" % stage, **metrics)
        progress.complete_step(
            active_progress_step,
            len(_POPS_PROGRESS_STAGES),
            _POPS_PROGRESS_STAGES[active_progress_step - 1],
            outcome_fields=_progress_outcome_fields(stage, metrics),
        )
        active_progress_step += 1
        progress.start_step(
            active_progress_step,
            len(_POPS_PROGRESS_STAGES),
            _POPS_PROGRESS_STAGES[active_progress_step - 1],
        )

    try:
        logger.record(
            "PARAM", "pops_run",
            genome_build=module.genome_build.value,
            method=module.method,
            random_seed=configuration.execution.random_seed,
            feature_matrix_chunks=module.feature_matrix_chunks,
            feature_selection_p_cutoff=module.feature_selection_p_cutoff,
            minimum_gene_count=module.minimum_gene_count,
        )
        logger.record(
            "OBSERVED", "pops_inputs",
            annotation_genes=annotation["gene_count"],
            outcome_genes=outcome["gene_count"],
            feature_genes=features["row_count"],
            features=features["feature_count"],
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
            active_progress_step = len(_POPS_PROGRESS_STAGES)
            progress.start_step(
                active_progress_step,
                len(_POPS_PROGRESS_STAGES),
                _POPS_PROGRESS_STAGES[active_progress_step - 1],
            )
            progress.complete_step(
                active_progress_step,
                len(_POPS_PROGRESS_STAGES),
                _POPS_PROGRESS_STAGES[active_progress_step - 1],
                outcome_fields=[
                    ("success", "Resume decision", "reused validated complete outputs"),
                    ("count", "Genes assigned PoPS scores", summary["genes_scored"]),
                    ("count", "Published result files", len(final_paths)),
                ],
            )
            print(_render_summary(
                summary, dataset, module, final_paths, log_path,
                configuration.logging.terminal_label_width,
            ))
            if ctx is not None:
                ctx["pops_output"] = result["pops_file"]
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
            covariance_path = None
            if module.target_error_covariance_file is not None:
                source = _required_file(
                    module.target_error_covariance_file, "PoPS target covariance",
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
            )
            upstream_handlers = list(logging.getLogger().handlers)
            upstream_level = logging.getLogger().level
            try:
                get_pops_args, pops_main = _load_upstream_entrypoints()
                active_progress_step = 2
                progress.start_step(
                    active_progress_step,
                    len(_POPS_PROGRESS_STAGES),
                    _POPS_PROGRESS_STAGES[active_progress_step - 1],
                )
                pops_main(
                    vars(get_pops_args(upstream_arguments)),
                    progress_callback=record_upstream_progress,
                )
            except SystemExit as exc:
                raise PopsError("Upstream PoPS exited before completion.") from exc
            finally:
                _restore_root_logging(upstream_handlers, upstream_level)

            staged_paths = {
                name: Path(str(staging_prefix) + path.name.removeprefix(prefix.name))
                for name, path in final_paths.items()
            }
            missing = [
                path for path in staged_paths.values()
                if not path.is_file() or path.stat().st_size <= 0
            ]
            if missing:
                raise PopsError(
                    "Upstream PoPS did not create required non-empty outputs: %s"
                    % ", ".join(str(path) for path in missing)
                )
            summary = _summarise_outputs(
                staged_paths, module, annotation, outcome, features,
            )
            all_known_suffixes = {
                module.output_layout.predictions_suffix,
                module.output_layout.coefficients_suffix,
                module.output_layout.marginals_suffix,
                module.output_layout.upstream_log_suffix,
                module.output_layout.training_data_suffix,
                module.output_layout.matrix_data_suffix,
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
        logger.record("OUTPUT", "pops_outputs", files=result["published_files"])
        logger.record("OBSERVED", "pops_summary", **{
            key: value for key, value in summary.items() if key != "top_genes"
        })
        for warning in summary["warnings"]:
            logger.warning(warning)
        logger.record("STATUS", "pops_run", status="COMPLETED")
        if ctx is not None:
            ctx["pops_output"] = result["pops_file"]
        progress.complete_step(
            active_progress_step,
            len(_POPS_PROGRESS_STAGES),
            _POPS_PROGRESS_STAGES[active_progress_step - 1],
            outcome_fields=[
                (
                    "success", "Output validation",
                    "all required files are non-empty",
                ),
                (
                    "count", "Genes assigned PoPS scores",
                    summary["genes_scored"],
                ),
                (
                    "count", "Selected features",
                    summary["selected_features"]
                    if summary["selected_features"] is not None
                    else "not recorded",
                ),
                (
                    "warning" if summary["warnings"] else "success",
                    "Scientific warnings", len(summary["warnings"]),
                ),
                ("count", "Published result files", len(final_paths)),
            ],
        )
        active_progress_step = 0
        print(_render_summary(
            summary, dataset, module, final_paths, log_path,
            configuration.logging.terminal_label_width,
        ))
        return result
    except BaseException as exc:
        if progress is not None and active_progress_step:
            progress.fail_step(
                active_progress_step,
                len(_POPS_PROGRESS_STAGES),
                _POPS_PROGRESS_STAGES[active_progress_step - 1],
            )
        if not logger.summary()["failed"]:
            logger.error("PoPS analysis failed: %s: %s" % (type(exc).__name__, exc))
        raise
    finally:
        if progress is not None:
            progress.close()
        logger.close()


__all__ = [
    "preflight_pops_pipeline",
    "run_pops_direct",
    "validate_pops_configuration",
]
