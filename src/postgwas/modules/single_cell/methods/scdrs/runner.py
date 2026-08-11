"""Validated execution of native scDRS score and downstream analyses."""

from __future__ import annotations

import csv
from dataclasses import dataclass, replace
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd
from scipy import sparse

from postgwas.core.completion import (
    configuration_digest,
    validate_completion_manifest,
    write_completion_manifest,
)
from postgwas.core.io.reports import write_delimited_report, write_yaml_report
from postgwas.core.paths import (
    configured_output_path,
    remove_owned_directory,
    require_nonempty_file,
    resolve_executable,
    validate_filename_component,
)
from postgwas.core.processes import run_checked_command
from postgwas.modules.single_cell.errors import SingleCellError


# scDRS 1.0.3 exposes no CLI seed option and calls score_cell without a seed;
# the supported upstream implementation therefore uses score_cell's fixed default.
_SCDRS_INTERNAL_RANDOM_SEED = 0


@dataclass(frozen=True)
class ScdrsPreflight:
    """Resolved inputs and scientific summaries needed for one scDRS run."""

    h5ad_file: Path
    gene_set_file: Path | None
    magma_gene_results_file: Path | None
    gene_identifier_map_file: Path | None
    covariate_file: Path | None
    executable: str
    software_version: str
    traits: tuple[str, ...]
    h5ad_summary: Mapping[str, Any]
    gene_set_summary: Mapping[str, Any]
    covariate_summary: Mapping[str, Any] | None
    gene_universe: tuple[str, ...]
    cell_names: tuple[str, ...]


@dataclass(frozen=True)
class ScdrsExecution:
    """Native scDRS outputs and execution metrics."""

    outputs: Mapping[str, Path]
    traits: tuple[str, ...]
    metrics: Mapping[str, Any]
    resumed: bool


def _load_anndata(path: Path):
    try:
        import anndata
    except ImportError as exc:
        raise SingleCellError(
            "scDRS H5AD validation requires the optional single-cell "
            "dependencies; install PostGWAS with the 'single-cell' extra"
        ) from exc
    try:
        return anndata.read_h5ad(path, backed="r")
    except Exception as exc:
        raise SingleCellError(
            "Cannot open the scDRS H5AD input %s: %s" % (path, exc)
        ) from exc


def _matrix_values(chunk):
    if sparse.issparse(chunk):
        return np.asarray(chunk.data)
    return np.asarray(chunk)


def _expressed_counts(chunk) -> tuple[np.ndarray, np.ndarray]:
    if sparse.issparse(chunk):
        expressed = chunk > 0
        return (
            np.asarray(expressed.sum(axis=1)).ravel(),
            np.asarray(expressed.sum(axis=0)).ravel(),
        )
    values = np.asarray(chunk)
    return (
        np.count_nonzero(values > 0, axis=1),
        np.count_nonzero(values > 0, axis=0),
    )


def _group_count_summary(obs, names: list[str]) -> dict[str, list[dict]]:
    summary = {}
    for name in names:
        counts = obs[name].value_counts(dropna=False)
        summary[name] = [
            {
                "value": None if pd.isna(value) else str(value),
                "cells": int(count),
            }
            for value, count in counts.items()
        ]
    return summary


def _continuous_annotation_summary(
    obs,
    names: list[str],
    *,
    context: str,
) -> dict[str, dict]:
    summary = {}
    for name in names:
        try:
            values = np.asarray(obs[name], dtype=float)
        except (TypeError, ValueError) as exc:
            raise SingleCellError(
                "scDRS correlation annotation %r must be numeric" % name
            ) from exc
        if np.isinf(values).any():
            raise SingleCellError(
                "scDRS correlation annotation %r contains infinite values" % name
            )
        observed = values[np.isfinite(values)]
        if observed.size < 2 or np.ptp(observed) <= 0:
            raise SingleCellError(
                "scDRS correlation annotation %r requires at least two "
                "finite, nonconstant values %s" % (name, context)
            )
        summary[name] = {
            "observed_cells": int(observed.size),
            "minimum": float(np.min(observed)),
            "maximum": float(np.max(observed)),
        }
    return summary


def _requested_annotation_names(method) -> tuple[list[str], list[str]]:
    categorical = list(method.downstream.group_analysis)
    if method.adjust_proportion_column is not None:
        categorical.append(method.adjust_proportion_column)
    continuous = list(method.downstream.correlation_analysis)
    return categorical, continuous


def _validate_annotation_values(adata, method) -> dict[str, Any]:
    categorical, continuous = _requested_annotation_names(method)
    requested = categorical + continuous
    missing_columns = [name for name in requested if name not in adata.obs.columns]
    if missing_columns:
        raise SingleCellError(
            "scDRS annotations are absent from adata.obs: %s"
            % ", ".join(missing_columns)
        )

    missing_counts: dict[str, int] = {}
    for name in requested:
        missing = int(adata.obs[name].isna().sum())
        missing_counts[name] = missing
        if missing and name == method.adjust_proportion_column:
            raise SingleCellError(
                "scDRS proportion-adjustment annotation %r contains %d missing "
                "values; the upstream cell-weight lookup requires every cell "
                "to have a group" % (name, missing)
            )
        if missing and not method.validation.allow_missing_annotation_values:
            raise SingleCellError(
                "scDRS annotation %r contains %d missing values; either curate "
                "the H5AD file or explicitly permit missing annotation values"
                % (name, missing)
            )
    return {
        "missing_values": missing_counts,
        "input_group_counts": _group_count_summary(adata.obs, categorical),
        "input_continuous_summaries": _continuous_annotation_summary(
            adata.obs,
            continuous,
            context="in the input H5AD",
        ),
    }


def validate_scdrs_h5ad(path: str | Path, method) -> dict[str, Any]:
    """Validate AnnData structure and the exact expression matrix scDRS reads."""
    h5ad_file = require_nonempty_file(
        path, "scDRS H5AD file", error_type=SingleCellError,
    )
    adata = _load_anndata(h5ad_file)
    try:
        if int(adata.n_obs) <= 0 or int(adata.n_vars) <= 0:
            raise SingleCellError("scDRS H5AD must contain cells and genes")
        if not adata.obs_names.is_unique:
            raise SingleCellError("scDRS H5AD cell identifiers must be unique")
        if not adata.var_names.is_unique:
            raise SingleCellError("scDRS H5AD gene identifiers must be unique")
        if adata.X is None:
            raise SingleCellError(
                "scDRS reads adata.X, but this H5AD file has no X matrix"
            )

        annotation_summary = _validate_annotation_values(adata, method)
        gene_cell_counts = np.zeros(int(adata.n_vars), dtype=np.int64)
        passing_cell_mask = np.zeros(int(adata.n_obs), dtype=bool)
        passing_cells = 0
        maximum_fractional_part = 0.0
        chunk_rows = method.validation.matrix_chunk_rows
        for start in range(0, int(adata.n_obs), chunk_rows):
            stop = min(start + chunk_rows, int(adata.n_obs))
            chunk = adata.X[start:stop]
            values = _matrix_values(chunk)
            if values.size and not np.isfinite(values).all():
                raise SingleCellError(
                    "scDRS H5AD adata.X contains NaN or infinite expression values"
                )
            if values.size and np.any(values < 0):
                raise SingleCellError(
                    "scDRS H5AD adata.X contains negative expression values"
                )
            if method.matrix_state == "raw_counts" and values.size:
                deviation = float(np.max(np.abs(values - np.rint(values))))
                maximum_fractional_part = max(
                    maximum_fractional_part, deviation,
                )
                if deviation > method.validation.raw_count_integer_tolerance:
                    raise SingleCellError(
                        "scDRS matrix_state is raw_counts, but adata.X contains "
                        "non-integer values beyond the configured tolerance"
                    )

            cell_gene_counts, _ = _expressed_counts(chunk)
            if method.filter_data:
                passing = cell_gene_counts >= method.minimum_genes_per_cell
            else:
                passing = np.ones(stop - start, dtype=bool)
            passing_cell_mask[start:stop] = passing
            passing_cells += int(np.sum(passing))
            if np.any(passing):
                _, filtered_gene_counts = _expressed_counts(chunk[passing])
                gene_cell_counts += filtered_gene_counts.astype(np.int64)

        if passing_cells <= 0:
            raise SingleCellError(
                "No H5AD cells pass the configured scDRS minimum_genes_per_cell"
            )
        if method.filter_data:
            passing_genes = gene_cell_counts >= method.minimum_cells_per_gene
        else:
            passing_genes = np.ones(int(adata.n_vars), dtype=bool)
        filtered_genes = int(np.sum(passing_genes))
        if filtered_genes <= 0:
            raise SingleCellError(
                "No H5AD genes pass the configured scDRS minimum_cells_per_gene"
            )
        categorical, continuous = _requested_annotation_names(method)
        analysis_obs = adata.obs.iloc[np.flatnonzero(passing_cell_mask)]
        if method.adjust_proportion_column is not None:
            # scDRS 1.0.3 enforces fewer than one group per ten analyzed cells
            # before it computes inverse group-size weights.
            group_count = int(
                analysis_obs[method.adjust_proportion_column].nunique()
            )
            if group_count >= 0.1 * passing_cells:
                raise SingleCellError(
                    "scDRS proportion adjustment has %d groups among %d "
                    "post-filter cells; the supported upstream algorithm "
                    "requires fewer than one group per ten cells"
                    % (group_count, passing_cells)
                )
        annotation_summary["analysis_group_counts"] = _group_count_summary(
            analysis_obs, categorical,
        )
        annotation_summary["analysis_continuous_summaries"] = (
            _continuous_annotation_summary(
                analysis_obs,
                continuous,
                context="after configured scDRS cell filtering",
            )
        )
        gene_universe = tuple(
            str(value)
            for value in adata.var_names[np.asarray(passing_genes, dtype=bool)]
        )
        cell_names = tuple(str(value) for value in adata.obs_names)
        return {
            "cells": int(adata.n_obs),
            "genes": int(adata.n_vars),
            "cells_after_configured_filter": passing_cells,
            "genes_after_configured_filter": filtered_genes,
            "gene_universe": gene_universe,
            "cell_names": cell_names,
            "matrix_source": "X",
            "matrix_state": method.matrix_state,
            "maximum_raw_count_fractional_part": maximum_fractional_part,
            "annotation_validation": annotation_summary,
        }
    finally:
        file_manager = getattr(adata, "file", None)
        if file_manager is not None:
            file_manager.close()


def _parse_gene_token(token: str, format_config) -> tuple[str, float | None]:
    cleaned = token.strip()
    if not cleaned:
        raise SingleCellError("scDRS gene sets must not contain empty gene entries")
    separator = format_config.weight_separator
    if separator not in cleaned:
        return cleaned, None
    gene, raw_weight = cleaned.rsplit(separator, 1)
    gene = gene.strip()
    if not gene or not raw_weight.strip():
        raise SingleCellError("Invalid scDRS gene-weight entry: %r" % cleaned)
    try:
        weight = float(raw_weight)
    except ValueError as exc:
        raise SingleCellError(
            "Invalid scDRS gene weight in entry %r" % cleaned
        ) from exc
    if not math.isfinite(weight):
        raise SingleCellError(
            "scDRS gene weights must be finite: %r" % cleaned
        )
    return gene, weight


def validate_scdrs_gene_sets(
    path: str | Path,
    method,
    *,
    gene_universe: tuple[str, ...],
) -> dict[str, Any]:
    """Validate native .gs syntax and effective overlap after H5AD filtering."""
    gene_set_file = require_nonempty_file(
        path, "scDRS .gs gene-set file", error_type=SingleCellError,
    )
    format_config = method.gene_set_format
    universe = set(gene_universe)
    records: list[dict[str, Any]] = []
    observed_traits: set[str] = set()
    try:
        handle = gene_set_file.open("r", encoding="utf-8", newline="")
    except OSError as exc:
        raise SingleCellError(
            "Cannot read scDRS gene-set file %s: %s" % (gene_set_file, exc)
        ) from exc
    with handle:
        reader = csv.DictReader(handle, delimiter=format_config.delimiter)
        expected_header = [
            format_config.trait_column, format_config.gene_set_column,
        ]
        if reader.fieldnames != expected_header:
            raise SingleCellError(
                "scDRS .gs header must be exactly %s"
                % format_config.delimiter.join(expected_header)
            )
        for line_number, row in enumerate(reader, start=2):
            trait = validate_filename_component(
                row.get(format_config.trait_column, ""),
                "scDRS trait on line %d" % line_number,
                error_type=SingleCellError,
            )
            if "@" in trait:
                raise SingleCellError(
                    "scDRS trait %r contains '@', which is reserved by the "
                    "upstream multi-score filename pattern" % trait
                )
            if trait in observed_traits:
                raise SingleCellError("Duplicate scDRS trait: %s" % trait)
            observed_traits.add(trait)
            raw_gene_set = (row.get(format_config.gene_set_column) or "").strip()
            tokens = raw_gene_set.split(format_config.gene_separator)
            parsed = [_parse_gene_token(token, format_config) for token in tokens]
            genes = [gene for gene, _ in parsed]
            weighted_genes = sum(weight is not None for _, weight in parsed)
            if weighted_genes not in {0, len(parsed)}:
                raise SingleCellError(
                    "scDRS trait %s mixes weighted and unweighted genes; all "
                    "entries must use the same native .gs representation" % trait
                )
            if len(genes) != len(set(genes)):
                raise SingleCellError(
                    "scDRS trait %s contains duplicate genes" % trait
                )
            effective = set(genes).intersection(universe)
            effective_count = len(effective)
            effective_fraction = effective_count / len(universe)
            validation = method.validation
            if effective_count < validation.minimum_effective_genes:
                raise SingleCellError(
                    "scDRS trait %s has %d genes in the filtered H5AD universe; "
                    "the configured minimum is %d"
                    % (
                        trait, effective_count,
                        validation.minimum_effective_genes,
                    )
                )
            if effective_fraction >= validation.maximum_effective_gene_fraction:
                raise SingleCellError(
                    "scDRS trait %s covers %.6g of the filtered H5AD gene "
                    "universe; it must be below the configured maximum %.6g"
                    % (
                        trait, effective_fraction,
                        validation.maximum_effective_gene_fraction,
                    )
                )
            records.append({
                "trait": trait,
                "input_genes": len(genes),
                "weighted_genes": weighted_genes,
                "effective_genes": effective_count,
                "effective_gene_fraction": effective_fraction,
            })
    if not records:
        raise SingleCellError("scDRS .gs file must contain at least one trait")
    return {"traits": records, "trait_count": len(records)}


def validate_scdrs_covariates(
    path: str | Path,
    method,
    *,
    cell_names: tuple[str, ...],
) -> dict[str, Any]:
    """Validate numeric scDRS covariates against AnnData cell identifiers."""
    covariate_file = require_nonempty_file(
        path, "scDRS covariate file", error_type=SingleCellError,
    )
    configured = method.covariate_format
    cells: set[str] = set()
    row_count = 0
    maximum_constant_deviation = 0.0
    try:
        handle = covariate_file.open("r", encoding="utf-8", newline="")
    except OSError as exc:
        raise SingleCellError(
            "Cannot read scDRS covariate file %s: %s" % (covariate_file, exc)
        ) from exc
    with handle:
        reader = csv.reader(handle, delimiter=configured.delimiter)
        header = next(reader, None)
        if header is None or len(header) < 2 or any(not value.strip() for value in header):
            raise SingleCellError(
                "scDRS covariate file requires a cell-ID column and at least "
                "one named numeric covariate"
            )
        if len(header) != len(set(header)):
            raise SingleCellError("scDRS covariate column names must be unique")
        constant_index = None
        if configured.constant_column is not None:
            if configured.constant_column not in header[1:]:
                raise SingleCellError(
                    "scDRS covariates require the configured constant column %r"
                    % configured.constant_column
                )
            constant_index = header.index(configured.constant_column)
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise SingleCellError(
                    "scDRS covariate line %d has %d fields; expected %d"
                    % (line_number, len(row), len(header))
                )
            cell = row[0].strip()
            if not cell:
                raise SingleCellError(
                    "scDRS covariate line %d has an empty cell identifier"
                    % line_number
                )
            if cell in cells:
                raise SingleCellError(
                    "Duplicate scDRS covariate cell identifier: %s" % cell
                )
            cells.add(cell)
            for column, raw_value in zip(header[1:], row[1:]):
                try:
                    value = float(raw_value)
                except ValueError as exc:
                    raise SingleCellError(
                        "scDRS covariate %r is non-numeric on line %d"
                        % (column, line_number)
                    ) from exc
                if not math.isfinite(value):
                    raise SingleCellError(
                        "scDRS covariate %r is not finite on line %d"
                        % (column, line_number)
                    )
            if constant_index is not None:
                constant_value = float(row[constant_index])
                deviation = abs(constant_value - 1.0)
                maximum_constant_deviation = max(
                    maximum_constant_deviation, deviation,
                )
                if deviation > configured.constant_tolerance:
                    raise SingleCellError(
                        "scDRS covariate constant column %r differs from one "
                        "beyond the configured tolerance on line %d"
                        % (configured.constant_column, line_number)
                    )
            row_count += 1
    h5ad_cells = set(cell_names)
    missing = h5ad_cells - cells
    unexpected = cells - h5ad_cells
    if configured.require_exact_cell_ids and (missing or unexpected):
        raise SingleCellError(
            "scDRS covariate cell identifiers must exactly match adata.obs_names "
            "(missing=%d, unexpected=%d)" % (len(missing), len(unexpected))
        )
    overlap = len(cells.intersection(h5ad_cells))
    if overlap <= 0:
        raise SingleCellError(
            "scDRS covariates do not overlap adata.obs_names"
        )
    overlap_fraction = overlap / len(h5ad_cells)
    if overlap_fraction <= configured.minimum_cell_overlap_fraction:
        raise SingleCellError(
            "scDRS covariates overlap %.6g of H5AD cells; upstream scDRS "
            "requires overlap above the configured fraction %.6g"
            % (
                overlap_fraction,
                configured.minimum_cell_overlap_fraction,
            )
        )
    return {
        "rows": row_count,
        "covariates": header[1:],
        "overlapping_cells": overlap,
        "overlap_fraction": overlap_fraction,
        "missing_h5ad_cells": len(missing),
        "unexpected_covariate_cells": len(unexpected),
        "constant_column": configured.constant_column,
        "maximum_constant_deviation": maximum_constant_deviation,
    }


def preflight_scdrs(
    method,
    executable_value: str | Path,
    *,
    dataset_id: str,
    pipeline_pending_magma: bool = False,
) -> ScdrsPreflight:
    """Resolve exact atlas inputs and either validate or defer the GWAS gene set."""
    h5ad_file = require_nonempty_file(
        method.input.h5ad_file, "scDRS H5AD file", error_type=SingleCellError,
    )
    executable = resolve_executable(
        executable_value, "scDRS executable", error_type=SingleCellError,
    )
    version_output = run_checked_command(
        [executable, *method.version_probe_arguments],
        "scDRS version probe",
        error_type=SingleCellError,
    )
    version_lines = [
        line.strip() for line in version_output.splitlines() if line.strip()
    ]
    if not version_lines:
        raise SingleCellError("scDRS version probe returned no version")
    software_version = version_lines[-1]
    if software_version not in method.supported_versions:
        raise SingleCellError(
            "Unsupported scDRS version %r; supported versions are %s"
            % (software_version, ", ".join(method.supported_versions))
        )
    h5ad_summary = validate_scdrs_h5ad(h5ad_file, method)
    gene_set_file = None
    magma_gene_results_file = None
    gene_identifier_map_file = None
    construction = method.magma_gene_set
    if construction.source == "file":
        gene_set_file = require_nonempty_file(
            method.input.gene_set_file,
            "scDRS .gs gene-set file",
            error_type=SingleCellError,
        )
        gene_set_summary = validate_scdrs_gene_sets(
            gene_set_file,
            method,
            gene_universe=h5ad_summary["gene_universe"],
        )
        traits = tuple(
            record["trait"] for record in gene_set_summary["traits"]
        )
    else:
        trait = validate_filename_component(
            construction.trait_pattern.format(dataset_id=dataset_id),
            "generated scDRS trait",
            error_type=SingleCellError,
        )
        if "@" in trait:
            raise SingleCellError(
                "Generated scDRS trait contains the reserved '@' character"
            )
        traits = (trait,)
        if not pipeline_pending_magma:
            magma_gene_results_file = require_nonempty_file(
                method.input.magma_gene_results_file,
                "MAGMA gene result used to construct scDRS gene sets",
                error_type=SingleCellError,
            )
        if construction.mapping_mode == "crosswalk":
            gene_identifier_map_file = require_nonempty_file(
                method.input.gene_identifier_map_file,
                "scDRS gene-identifier crosswalk",
                error_type=SingleCellError,
            )
        gene_set_summary = {
            "source": "magma",
            "status": (
                "awaiting_pipeline_magma"
                if pipeline_pending_magma else "ready_for_construction"
            ),
            "trait_count": 1,
            "traits": [{"trait": trait}],
        }
    covariate_file = None
    covariate_summary = None
    if method.input.covariate_file is not None:
        covariate_file = require_nonempty_file(
            method.input.covariate_file,
            "scDRS covariate file",
            error_type=SingleCellError,
        )
        covariate_summary = validate_scdrs_covariates(
            covariate_file,
            method,
            cell_names=h5ad_summary["cell_names"],
        )
    return ScdrsPreflight(
        h5ad_file=h5ad_file,
        gene_set_file=gene_set_file,
        magma_gene_results_file=magma_gene_results_file,
        gene_identifier_map_file=gene_identifier_map_file,
        covariate_file=covariate_file,
        executable=executable,
        software_version=software_version,
        traits=traits,
        h5ad_summary={
            key: value for key, value in h5ad_summary.items()
            if key not in {"gene_universe", "cell_names"}
        },
        gene_set_summary=gene_set_summary,
        covariate_summary=covariate_summary,
        gene_universe=h5ad_summary["gene_universe"],
        cell_names=h5ad_summary["cell_names"],
    )


def _normalize_gene_identifier(value: object) -> str | None:
    if value is None or pd.isna(value):
        return None
    if isinstance(value, (float, np.floating)) and float(value).is_integer():
        return str(int(value))
    cleaned = str(value).strip()
    return cleaned or None


def _read_magma_gene_statistics(path: Path, construction) -> pd.DataFrame:
    try:
        frame = pd.read_csv(
            path,
            sep=construction.table_delimiter_pattern,
            comment="#",
            engine="python",
            dtype={construction.gene_id_column: str},
        )
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise SingleCellError(
            "Cannot read MAGMA gene statistics %s: %s" % (path, exc)
        ) from exc
    required = {construction.gene_id_column, construction.z_score_column}
    missing = sorted(required - set(frame.columns))
    if missing or frame.empty:
        raise SingleCellError(
            "MAGMA gene statistics are empty or missing columns: %s"
            % ", ".join(missing or sorted(required))
        )
    genes = frame[construction.gene_id_column].map(_normalize_gene_identifier)
    if genes.isna().any():
        raise SingleCellError(
            "MAGMA gene statistics contain missing or empty gene identifiers"
        )
    if genes.duplicated().any():
        raise SingleCellError(
            "MAGMA gene statistics contain duplicate normalized gene identifiers"
        )
    z_scores = pd.to_numeric(frame[construction.z_score_column], errors="coerce")
    if z_scores.isna().any() or not np.isfinite(z_scores.to_numpy(dtype=float)).all():
        raise SingleCellError(
            "MAGMA gene statistics contain missing or non-finite Z scores"
        )
    return pd.DataFrame({"source_gene": genes, "z_score": z_scores.astype(float)})


def _read_identifier_crosswalk(path: Path, construction) -> tuple[dict[str, str], dict]:
    try:
        frame = pd.read_csv(
            path,
            sep=construction.mapping_delimiter,
            dtype=str,
            keep_default_na=False,
        )
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise SingleCellError(
            "Cannot read scDRS gene-identifier crosswalk %s: %s" % (path, exc)
        ) from exc
    required = {
        construction.mapping_source_column,
        construction.mapping_target_column,
    }
    missing = sorted(required - set(frame.columns))
    if missing or frame.empty:
        raise SingleCellError(
            "scDRS gene-identifier crosswalk is empty or missing columns: %s"
            % ", ".join(missing or sorted(required))
        )
    source = frame[construction.mapping_source_column].map(
        _normalize_gene_identifier
    )
    target = frame[construction.mapping_target_column].map(
        _normalize_gene_identifier
    )
    if source.isna().any() or target.isna().any():
        raise SingleCellError(
            "scDRS gene-identifier crosswalk contains empty source or target IDs"
        )
    pairs = pd.DataFrame({"source": source, "target": target})
    duplicate_sources = pairs["source"].duplicated(keep=False)
    duplicate_targets = pairs["target"].duplicated(keep=False)
    if duplicate_sources.any() or duplicate_targets.any():
        raise SingleCellError(
            "scDRS gene-identifier crosswalk is not one-to-one "
            "(repeated_source_rows=%d, repeated_target_rows=%d); ambiguous "
            "identifier mappings require manual curation"
            % (int(duplicate_sources.sum()), int(duplicate_targets.sum()))
        )
    mapping = dict(zip(pairs["source"], pairs["target"]))
    return mapping, {
        "crosswalk_rows": len(pairs),
        "repeated_source_rows": 0,
        "repeated_target_rows": 0,
    }


def prepare_scdrs_gene_set(
    preflight: ScdrsPreflight,
    method,
    *,
    statistics_file: Path,
    gene_set_file: Path,
    mapping_report: Path,
    logger,
) -> ScdrsPreflight:
    """Convert validated MAGMA Z statistics to one auditable native scDRS .gs."""
    construction = method.magma_gene_set
    if construction.source != "magma":
        return preflight
    if preflight.magma_gene_results_file is None:
        raise SingleCellError(
            "MAGMA gene results are required to construct the pipeline scDRS gene set"
        )
    statistics = _read_magma_gene_statistics(
        preflight.magma_gene_results_file, construction,
    )
    mapping_metrics = {
        "mapping_mode": construction.mapping_mode,
        "source_identifier_type": construction.source_identifier_type,
        "target_identifier_type": construction.target_identifier_type,
        "input_magma_genes": len(statistics),
    }
    if construction.mapping_mode == "crosswalk":
        if preflight.gene_identifier_map_file is None:
            raise SingleCellError(
                "A gene-identifier crosswalk is required for scDRS mapping_mode "
                "crosswalk"
            )
        mapping, crosswalk_metrics = _read_identifier_crosswalk(
            preflight.gene_identifier_map_file, construction,
        )
        mapping_metrics.update(crosswalk_metrics)
        statistics["target_gene"] = statistics["source_gene"].map(mapping)
        unmapped = int(statistics["target_gene"].isna().sum())
        mapping_metrics["unmapped_magma_genes"] = unmapped
        if unmapped and construction.unmapped_policy == "error":
            raise SingleCellError(
                "%d MAGMA genes are absent from the configured identifier "
                "crosswalk and unmapped_policy is error" % unmapped
            )
        statistics = statistics.loc[
            statistics["target_gene"].notna()
        ].copy()
    else:
        statistics["target_gene"] = statistics["source_gene"]
        mapping_metrics.update({
            "crosswalk_rows": None,
            "unmapped_magma_genes": 0,
        })
    if statistics["target_gene"].duplicated().any():
        raise SingleCellError(
            "MAGMA-to-H5AD identifier conversion produced duplicate target genes"
        )
    if len(statistics) < construction.minimum_genes:
        raise SingleCellError(
            "Only %d MAGMA genes remain after identifier mapping; scDRS "
            "munge-gs requires at least the configured minimum %d"
            % (len(statistics), construction.minimum_genes)
        )
    trait = preflight.traits[0]
    if trait == construction.statistics_gene_column:
        raise SingleCellError(
            "Generated scDRS trait must differ from statistics_gene_column"
        )
    records = [
        {
            construction.statistics_gene_column: row.target_gene,
            trait: float(row.z_score),
        }
        for row in statistics.itertuples(index=False)
    ]
    write_delimited_report(
        records,
        statistics_file,
        fieldnames=[construction.statistics_gene_column, trait],
        delimiter=method.gene_set_format.delimiter,
        null_value=construction.statistics_null_value,
    )
    run_checked_command(
        [
            preflight.executable,
            "munge-gs",
            "--out-file", str(gene_set_file),
            "--zscore-file", str(statistics_file),
            "--weight", construction.weight,
            "--n-min", str(construction.minimum_genes),
            "--n-max", str(construction.maximum_genes),
        ],
        "Construct scDRS gene set from MAGMA Z scores",
        logger=logger,
        error_type=SingleCellError,
        expected_outputs=[gene_set_file],
    )
    gene_set_summary = validate_scdrs_gene_sets(
        gene_set_file,
        method,
        gene_universe=preflight.gene_universe,
    )
    mapping_metrics.update({
        "mapped_magma_genes": len(statistics),
        "mapped_genes_in_h5ad": int(
            statistics["target_gene"].isin(set(preflight.gene_universe)).sum()
        ),
        "unmapped_policy": construction.unmapped_policy,
        "ambiguous_policy": construction.ambiguous_policy,
        "selection_weight": construction.weight,
        "selection_minimum_genes": construction.minimum_genes,
        "selection_maximum_genes": construction.maximum_genes,
        "generated_gene_set": gene_set_summary,
    })
    write_yaml_report(
        {
            "schema_version": 1,
            "status": "COMPLETED",
            "source_gene_results": str(preflight.magma_gene_results_file),
            "crosswalk": (
                None if preflight.gene_identifier_map_file is None
                else str(preflight.gene_identifier_map_file)
            ),
            **mapping_metrics,
        },
        mapping_report,
    )
    logger.record("RESULT", "scdrs_gene_mapping", **mapping_metrics)
    return replace(
        preflight,
        gene_set_file=gene_set_file,
        gene_set_summary=gene_set_summary,
    )


def expected_scdrs_outputs(
    directory: Path,
    traits: tuple[str, ...],
    method,
    output_layout,
) -> dict[str, Path]:
    """Resolve every native output implied by the configured scDRS analyses."""
    outputs: dict[str, Path] = {}
    for trait_index, trait in enumerate(traits, start=1):
        values = {"trait": trait}
        outputs["score_%d" % trait_index] = configured_output_path(
            directory,
            output_layout.scdrs_score_file_pattern,
            error_type=SingleCellError,
            **values,
        )
        outputs["full_score_%d" % trait_index] = configured_output_path(
            directory,
            output_layout.scdrs_full_score_file_pattern,
            error_type=SingleCellError,
            **values,
        )
        for annotation_index, annotation in enumerate(
            method.downstream.group_analysis, start=1,
        ):
            validate_filename_component(
                annotation,
                "scDRS group annotation",
                error_type=SingleCellError,
            )
            outputs[
                "group_%d_%d" % (trait_index, annotation_index)
            ] = configured_output_path(
                directory,
                output_layout.scdrs_group_file_pattern,
                error_type=SingleCellError,
                trait=trait,
                annotation=annotation,
            )
        if method.downstream.correlation_analysis:
            outputs["correlation_%d" % trait_index] = configured_output_path(
                directory,
                output_layout.scdrs_correlation_file_pattern,
                error_type=SingleCellError,
                **values,
            )
        if method.downstream.gene_analysis:
            outputs["gene_%d" % trait_index] = configured_output_path(
                directory,
                output_layout.scdrs_gene_file_pattern,
                error_type=SingleCellError,
                **values,
            )
    if len(outputs) != len(set(outputs.values())):
        raise SingleCellError("Configured scDRS output patterns collide")
    return outputs


def _bool_argument(value: bool) -> str:
    return "True" if value else "False"


def _compute_command(preflight: ScdrsPreflight, method, output: Path) -> list[str]:
    command = [
        preflight.executable,
        "compute-score",
        "--h5ad-file", str(preflight.h5ad_file),
        "--h5ad-species", method.h5ad_species,
        "--gs-file", str(preflight.gene_set_file),
        "--gs-species", method.gene_set_species,
        "--out-folder", str(output),
        "--weight-opt", method.weight_option,
        "--flag-filter-data", _bool_argument(method.filter_data),
        "--flag-raw-count", _bool_argument(method.matrix_state == "raw_counts"),
        "--n-ctrl", str(method.control_gene_sets),
        "--min-genes", str(method.minimum_genes_per_cell),
        "--min-cells", str(method.minimum_cells_per_gene),
        "--flag-return-ctrl-raw-score",
        _bool_argument(method.return_control_raw_score),
        "--flag-return-ctrl-norm-score",
        _bool_argument(method.return_control_normalized_score),
    ]
    if preflight.covariate_file is not None:
        command.extend(["--cov-file", str(preflight.covariate_file)])
    if method.adjust_proportion_column is not None:
        command.extend(["--adj-prop", method.adjust_proportion_column])
    return command


def _downstream_command(
    preflight: ScdrsPreflight,
    method,
    output: Path,
    output_layout,
) -> list[str]:
    full_scores = configured_output_path(
        output,
        output_layout.scdrs_full_score_file_pattern,
        error_type=SingleCellError,
        trait="@",
    )
    command = [
        preflight.executable,
        "perform-downstream",
        "--h5ad-file", str(preflight.h5ad_file),
        "--score-file", str(full_scores),
        "--out-folder", str(output),
        "--flag-filter-data", _bool_argument(method.filter_data),
        "--flag-raw-count", _bool_argument(method.matrix_state == "raw_counts"),
        "--min-genes", str(method.minimum_genes_per_cell),
        "--min-cells", str(method.minimum_cells_per_gene),
        "--knn-n-neighbors", str(method.downstream.knn_neighbors),
        "--knn-n-pcs", str(method.downstream.knn_principal_components),
    ]
    if method.downstream.group_analysis:
        command.extend([
            "--group-analysis", ",".join(method.downstream.group_analysis),
        ])
    if method.downstream.correlation_analysis:
        command.extend([
            "--corr-analysis",
            ",".join(method.downstream.correlation_analysis),
        ])
    if method.downstream.gene_analysis:
        command.append("--gene-analysis")
    return command


def run_scdrs(
    *,
    preflight: ScdrsPreflight,
    configuration,
    paths: Mapping[str, Path],
    dataset: str,
    logger,
) -> ScdrsExecution:
    """Run scDRS atomically, or resume only checksum-matched native outputs."""
    method = configuration.modules.single_cell.scdrs
    layout = configuration.modules.single_cell.output_layout
    engine = paths["scdrs_engine_directory"]
    staging = paths["scdrs_staging_directory"]
    qc_report = paths["scdrs_qc_report"]
    manifest_path = paths["scdrs_completion_manifest"]
    if engine not in qc_report.parents:
        raise SingleCellError(
            "modules.single_cell.output_layout.scdrs_qc_report must be inside "
            "scdrs_engine_directory"
        )
    generated_paths = {
        paths["scdrs_generated_gene_statistics"],
        paths["scdrs_generated_gene_set"],
        paths["scdrs_gene_mapping_report"],
    }
    if any(engine not in path.parents for path in generated_paths):
        raise SingleCellError(
            "Configured generated scDRS gene-set artifacts must be inside "
            "scdrs_engine_directory"
        )
    outputs = expected_scdrs_outputs(
        engine, preflight.traits, method, layout,
    )
    outputs["qc_report"] = qc_report
    construction = method.magma_gene_set
    if construction.source == "magma":
        outputs.update({
            "generated_gene_statistics": paths[
                "scdrs_generated_gene_statistics"
            ],
            "generated_gene_set": paths["scdrs_generated_gene_set"],
            "gene_mapping_report": paths["scdrs_gene_mapping_report"],
        })
        if len(outputs) != len(set(outputs.values())):
            raise SingleCellError("Configured scDRS output paths collide")
        completion_inputs = {
            "h5ad": preflight.h5ad_file,
            "magma_gene_results": preflight.magma_gene_results_file,
        }
        if preflight.gene_identifier_map_file is not None:
            completion_inputs["gene_identifier_map"] = (
                preflight.gene_identifier_map_file
            )
    else:
        completion_inputs = {
            "h5ad": preflight.h5ad_file,
            "gene_sets": preflight.gene_set_file,
        }
    if preflight.covariate_file is not None:
        completion_inputs["covariates"] = preflight.covariate_file
    digest = configuration_digest({
        "scdrs": method.model_dump(mode="json"),
        "magma_gene_set_source": (
            configuration.modules.magma.model_dump(mode="json")
            if construction.source == "magma" else None
        ),
        "output_layout": {
            key: value for key, value in layout.model_dump().items()
            if key.startswith("scdrs_")
        },
        "executable": preflight.executable,
    })
    if (
        configuration.run.resume
        and not configuration.run.overwrite
        and manifest_path.is_file()
    ):
        manifest = validate_completion_manifest(
            manifest_path,
            dataset_id=dataset,
            module="single_cell.scdrs",
            genome_build="not_applicable_gene_level_input",
            configuration_sha256=digest,
            inputs=completion_inputs,
            outputs=outputs,
            error_type=SingleCellError,
        )
        return ScdrsExecution(
            outputs=outputs,
            traits=preflight.traits,
            metrics=dict(manifest.get("metrics", {})),
            resumed=True,
        )
    if (engine.exists() or manifest_path.exists()) and not configuration.run.overwrite:
        raise SingleCellError(
            "Existing or incomplete scDRS output was found; use --resume for "
            "a matching complete run or --overwrite to replace it"
        )
    if staging.exists():
        if not configuration.run.overwrite:
            raise SingleCellError(
                "An isolated incomplete scDRS run exists at %s; review it or "
                "use --overwrite" % staging
            )
        remove_owned_directory(
            staging,
            configuration.run.output_directory,
            "incomplete scDRS staging directory",
            error_type=SingleCellError,
        )
    if configuration.run.overwrite:
        remove_owned_directory(
            engine,
            configuration.run.output_directory,
            "scDRS engine directory",
            error_type=SingleCellError,
        )
        manifest_path.unlink(missing_ok=True)

    staging.mkdir(parents=True, exist_ok=False)
    staged_outputs = expected_scdrs_outputs(
        staging, preflight.traits, method, layout,
    )
    execution_preflight = preflight
    if construction.source == "magma":
        staged_statistics = staging / paths[
            "scdrs_generated_gene_statistics"
        ].relative_to(engine)
        staged_gene_set = staging / paths[
            "scdrs_generated_gene_set"
        ].relative_to(engine)
        staged_mapping_report = staging / paths[
            "scdrs_gene_mapping_report"
        ].relative_to(engine)
        execution_preflight = prepare_scdrs_gene_set(
            preflight,
            method,
            statistics_file=staged_statistics,
            gene_set_file=staged_gene_set,
            mapping_report=staged_mapping_report,
            logger=logger,
        )
        staged_outputs.update({
            "generated_gene_statistics": staged_statistics,
            "generated_gene_set": staged_gene_set,
            "gene_mapping_report": staged_mapping_report,
        })
    run_checked_command(
        _compute_command(execution_preflight, method, staging),
        "scDRS cell-score computation",
        logger=logger,
        error_type=SingleCellError,
        expected_outputs=[
            path for name, path in staged_outputs.items()
            if name.startswith(("score_", "full_score_"))
        ],
    )
    downstream_requested = bool(
        method.downstream.group_analysis
        or method.downstream.correlation_analysis
        or method.downstream.gene_analysis
    )
    if downstream_requested:
        run_checked_command(
            _downstream_command(execution_preflight, method, staging, layout),
            "scDRS downstream analyses",
            logger=logger,
            error_type=SingleCellError,
            expected_outputs=[
                path for name, path in staged_outputs.items()
                if not name.startswith(("score_", "full_score_"))
            ],
        )

    staged_qc = staging / qc_report.relative_to(engine)
    metrics = {
        "traits": len(execution_preflight.traits),
        "trait_ids": list(execution_preflight.traits),
        "cells": preflight.h5ad_summary["cells"],
        "genes": preflight.h5ad_summary["genes"],
        "cells_after_configured_filter": (
            preflight.h5ad_summary["cells_after_configured_filter"]
        ),
        "genes_after_configured_filter": (
            preflight.h5ad_summary["genes_after_configured_filter"]
        ),
        "control_gene_sets": method.control_gene_sets,
        "software_version": preflight.software_version,
        "matrix_source": "X",
        "matrix_state": method.matrix_state,
        "native_output_files": len(
            [
                name for name in staged_outputs
                if name not in {
                    "generated_gene_statistics", "generated_gene_set",
                    "gene_mapping_report",
                }
            ]
        ),
        "gene_set_source": construction.source,
        "upstream_internal_random_seed": _SCDRS_INTERNAL_RANDOM_SEED,
        "postgwas_execution_seed_applied": False,
    }
    magma_provenance = None
    if construction.source == "magma":
        magma = configuration.modules.magma
        primary_mapping = magma.mapping.primary
        definition = magma.mapping.definitions[primary_mapping]
        magma_provenance = {
            "genome_build": magma.genome_build.value,
            "population": magma.population.value,
            "gene_window_upstream_kb": magma.gene_window_upstream_kb,
            "gene_window_downstream_kb": magma.gene_window_downstream_kb,
            "gene_model": magma.gene_model,
            "primary_mapping": primary_mapping,
            "mapping_method": definition.method,
            "gene_id_type": definition.gene_id_type,
            "annotation_source": definition.source_name,
            "annotation_version": definition.source_version,
            "result_statistic_type": definition.result_statistic_type,
        }
    write_yaml_report(
        {
            "schema_version": 1,
            "status": "COMPLETED",
            "dataset_id": dataset,
            "method": "scdrs",
            "software_version": preflight.software_version,
            "inputs": {
                "h5ad": str(preflight.h5ad_file),
                "gene_sets": (
                    str(paths["scdrs_generated_gene_set"])
                    if construction.source == "magma"
                    else str(execution_preflight.gene_set_file)
                ),
                "magma_gene_results": (
                    None if preflight.magma_gene_results_file is None
                    else str(preflight.magma_gene_results_file)
                ),
                "gene_identifier_map": (
                    None if preflight.gene_identifier_map_file is None
                    else str(preflight.gene_identifier_map_file)
                ),
                "covariates": (
                    None if preflight.covariate_file is None
                    else str(preflight.covariate_file)
                ),
            },
            "h5ad_validation": dict(preflight.h5ad_summary),
            "gene_set_validation": dict(execution_preflight.gene_set_summary),
            "covariate_validation": (
                None if preflight.covariate_summary is None
                else dict(preflight.covariate_summary)
            ),
            "magma_provenance": magma_provenance,
            "scientific_settings": {
                "h5ad_species": method.h5ad_species,
                "gene_set_species": method.gene_set_species,
                "gene_set_source": construction.source,
                "matrix_source": "X",
                "matrix_state": method.matrix_state,
                "filter_data": method.filter_data,
                "minimum_genes_per_cell": method.minimum_genes_per_cell,
                "minimum_cells_per_gene": method.minimum_cells_per_gene,
                "control_gene_sets": method.control_gene_sets,
                "upstream_internal_random_seed": (
                    _SCDRS_INTERNAL_RANDOM_SEED
                ),
                "postgwas_execution_seed_applied": False,
                "weight_option": method.weight_option,
                "group_analysis": list(method.downstream.group_analysis),
                "correlation_analysis": list(
                    method.downstream.correlation_analysis
                ),
                "gene_analysis": method.downstream.gene_analysis,
            },
            "metrics": metrics,
        },
        staged_qc,
    )
    engine.parent.mkdir(parents=True, exist_ok=True)
    staging.replace(engine)
    write_completion_manifest(
        manifest_path,
        dataset_id=dataset,
        module="single_cell.scdrs",
        genome_build="not_applicable_gene_level_input",
        configuration_sha256=digest,
        inputs=completion_inputs,
        outputs=outputs,
        metrics=metrics,
        error_type=SingleCellError,
    )
    return ScdrsExecution(
        outputs=outputs,
        traits=execution_preflight.traits,
        metrics=metrics,
        resumed=False,
    )


__all__ = [
    "ScdrsExecution",
    "ScdrsPreflight",
    "expected_scdrs_outputs",
    "preflight_scdrs",
    "run_scdrs",
    "validate_scdrs_covariates",
    "validate_scdrs_gene_sets",
    "validate_scdrs_h5ad",
]
