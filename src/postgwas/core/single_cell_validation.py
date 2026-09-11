"""Read-only single-cell file contracts, independent of analysis modules.

Expression arrays are scanned in backed chunks and never retained in the session
cache. Only identifier sets and compact summaries are reusable. Policies are
provided by the consumer's validated configuration; this file does not infer
species, expression scale, gene identifiers, or cross-file compatibility.

The supported scDRS X-matrix/intercept/proportion contracts are preserved from
the pinned upstream implementation:
https://github.com/martinjzhang/scDRS/tree/1518b56c73eb23c6ca20999d439aaa31a6b7624b
"""

from __future__ import annotations

import csv
from copy import deepcopy
from dataclasses import asdict, dataclass
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

from postgwas.core.input_validation import record_file_validation, validate_once
from postgwas.core.paths import require_nonempty_file, validate_filename_component


@dataclass(frozen=True)
class H5adValidationPolicy:
    """Only settings that affect the read-only H5AD validation result."""

    matrix_state: str
    filter_data: bool
    minimum_genes_per_cell: int
    minimum_cells_per_gene: int
    matrix_chunk_rows: int
    raw_count_integer_tolerance: float
    group_analysis: tuple[str, ...]
    correlation_analysis: tuple[str, ...]
    adjust_proportion_column: str | None
    allow_missing_annotation_values: bool


@dataclass(frozen=True)
class CellCovariateEvidence:
    """Immutable parsed identifiers and numeric-check summaries, not a matrix."""

    cells: frozenset[str]
    columns: tuple[str, ...]
    rows: int
    constant_column: str | None
    maximum_constant_deviation: float


@dataclass(frozen=True)
class GeneSetEvidence:
    """Native gene-set membership needed for later atlas compatibility checks."""

    trait: str
    genes: tuple[str, ...]
    weighted_genes: int


def _load_anndata(path: Path, *, error_type):
    try:
        import anndata
    except ImportError as exc:
        raise error_type(
            "scDRS H5AD validation requires the optional single-cell "
            "dependencies; install PostGWAS with the 'single-cell' extra"
        ) from exc
    try:
        return anndata.read_h5ad(path, backed="r")
    except Exception as exc:
        raise error_type(
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
    error_type,
) -> dict[str, dict]:
    summary = {}
    for name in names:
        try:
            values = np.asarray(obs[name], dtype=float)
        except (TypeError, ValueError) as exc:
            raise error_type(
                "scDRS correlation annotation %r must be numeric" % name
            ) from exc
        if np.isinf(values).any():
            raise error_type(
                "scDRS correlation annotation %r contains infinite values" % name
            )
        observed = values[np.isfinite(values)]
        if observed.size < 2 or np.ptp(observed) <= 0:
            raise error_type(
                "scDRS correlation annotation %r requires at least two "
                "finite, nonconstant values %s" % (name, context)
            )
        summary[name] = {
            "observed_cells": int(observed.size),
            "minimum": float(np.min(observed)),
            "maximum": float(np.max(observed)),
        }
    return summary


def _requested_annotation_names(policy) -> tuple[list[str], list[str]]:
    categorical = list(policy.group_analysis)
    if policy.adjust_proportion_column is not None:
        categorical.append(policy.adjust_proportion_column)
    continuous = list(policy.correlation_analysis)
    return categorical, continuous


def _validate_annotation_values(adata, policy, *, error_type) -> dict[str, Any]:
    categorical, continuous = _requested_annotation_names(policy)
    requested = categorical + continuous
    missing_columns = [name for name in requested if name not in adata.obs.columns]
    if missing_columns:
        raise error_type(
            "scDRS annotations are absent from adata.obs: %s"
            % ", ".join(missing_columns)
        )

    missing_counts: dict[str, int] = {}
    for name in requested:
        missing = int(adata.obs[name].isna().sum())
        missing_counts[name] = missing
        if missing and name == policy.adjust_proportion_column:
            raise error_type(
                "scDRS proportion-adjustment annotation %r contains %d missing "
                "values; the upstream cell-weight lookup requires every cell "
                "to have a group" % (name, missing)
            )
        if missing and not policy.allow_missing_annotation_values:
            raise error_type(
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
            error_type=error_type,
        ),
    }


def _inspect_h5ad(
    path: Path, policy: H5adValidationPolicy, *, error_type,
) -> dict[str, Any]:
    """Validate AnnData structure and the exact expression matrix scDRS reads."""
    adata = _load_anndata(path, error_type=error_type)
    try:
        if int(adata.n_obs) <= 0 or int(adata.n_vars) <= 0:
            raise error_type("scDRS H5AD must contain cells and genes")
        if not adata.obs_names.is_unique:
            raise error_type("scDRS H5AD cell identifiers must be unique")
        if not adata.var_names.is_unique:
            raise error_type("scDRS H5AD gene identifiers must be unique")
        if adata.X is None:
            raise error_type(
                "scDRS reads adata.X, but this H5AD file has no X matrix"
            )

        annotation_summary = _validate_annotation_values(
            adata, policy, error_type=error_type,
        )
        gene_cell_counts = np.zeros(int(adata.n_vars), dtype=np.int64)
        passing_cell_mask = np.zeros(int(adata.n_obs), dtype=bool)
        passing_cells = 0
        maximum_fractional_part = 0.0
        chunk_rows = policy.matrix_chunk_rows
        for start in range(0, int(adata.n_obs), chunk_rows):
            stop = min(start + chunk_rows, int(adata.n_obs))
            chunk = adata.X[start:stop]
            values = _matrix_values(chunk)
            if values.size and not np.isfinite(values).all():
                raise error_type(
                    "scDRS H5AD adata.X contains NaN or infinite expression values"
                )
            if values.size and np.any(values < 0):
                raise error_type(
                    "scDRS H5AD adata.X contains negative expression values"
                )
            if policy.matrix_state == "raw_counts" and values.size:
                deviation = float(np.max(np.abs(values - np.rint(values))))
                maximum_fractional_part = max(
                    maximum_fractional_part, deviation,
                )
                if deviation > policy.raw_count_integer_tolerance:
                    raise error_type(
                        "scDRS matrix_state is raw_counts, but adata.X contains "
                        "non-integer values beyond the configured tolerance"
                    )

            cell_gene_counts, _ = _expressed_counts(chunk)
            if policy.filter_data:
                passing = cell_gene_counts >= policy.minimum_genes_per_cell
            else:
                passing = np.ones(stop - start, dtype=bool)
            passing_cell_mask[start:stop] = passing
            passing_cells += int(np.sum(passing))
            if np.any(passing):
                _, filtered_gene_counts = _expressed_counts(chunk[passing])
                gene_cell_counts += filtered_gene_counts.astype(np.int64)

        if passing_cells <= 0:
            raise error_type(
                "No H5AD cells pass the configured scDRS minimum_genes_per_cell"
            )
        if policy.filter_data:
            passing_genes = gene_cell_counts >= policy.minimum_cells_per_gene
        else:
            passing_genes = np.ones(int(adata.n_vars), dtype=bool)
        filtered_genes = int(np.sum(passing_genes))
        if filtered_genes <= 0:
            raise error_type(
                "No H5AD genes pass the configured scDRS minimum_cells_per_gene"
            )
        categorical, continuous = _requested_annotation_names(policy)
        analysis_obs = adata.obs.iloc[np.flatnonzero(passing_cell_mask)]
        if policy.adjust_proportion_column is not None:
            # scDRS 1.0.3 enforces fewer than one group per ten analyzed cells
            # before it computes inverse group-size weights.
            group_count = int(
                analysis_obs[policy.adjust_proportion_column].nunique()
            )
            if group_count >= 0.1 * passing_cells:
                raise error_type(
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
                error_type=error_type,
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
            "matrix_state": policy.matrix_state,
            "maximum_raw_count_fractional_part": maximum_fractional_part,
            "annotation_validation": annotation_summary,
        }
    finally:
        file_manager = getattr(adata, "file", None)
        if file_manager is not None:
            file_manager.close()


def validate_h5ad_expression(
    path: str | Path, *, policy: H5adValidationPolicy, error_type=ValueError,
) -> dict[str, Any]:
    """Inspect the consumer's declared X-matrix contract once per exact policy."""
    source = require_nonempty_file(path, "scDRS H5AD file", error_type=error_type)

    def inspect():
        result = _inspect_h5ad(source, policy, error_type=error_type)
        record_file_validation(
            source, "H5AD expression input",
            checks=(
                "unique cell and gene identifiers", "complete X-matrix value scan",
                "requested annotations", "configured expression-filter feasibility",
            ),
            metrics={
                key: value for key, value in result.items()
                if key not in {"gene_universe", "cell_names"}
            },
        )
        return result

    return deepcopy(validate_once(
        (source,), ("h5ad_expression", asdict(policy)), inspect,
        error_type=error_type,
    ))


def _inspect_cell_covariates(
    path: Path, *, delimiter: str,
    constant_column: str | None, constant_tolerance: float,
    error_type,
) -> CellCovariateEvidence:
    """Validate the table's own identifiers, row widths, numbers and intercept."""
    covariate_file = path
    cells: set[str] = set()
    row_count = 0
    maximum_constant_deviation = 0.0
    try:
        handle = covariate_file.open("r", encoding="utf-8", newline="")
    except OSError as exc:
        raise error_type(
            "Cannot read scDRS covariate file %s: %s" % (covariate_file, exc)
        ) from exc
    with handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader, None)
        if header is None or len(header) < 2 or any(not value.strip() for value in header):
            raise error_type(
                "scDRS covariate file requires a cell-ID column and at least "
                "one named numeric covariate"
            )
        if len(header) != len(set(header)):
            raise error_type("scDRS covariate column names must be unique")
        constant_index = None
        if constant_column is not None:
            if constant_column not in header[1:]:
                raise error_type(
                    "scDRS covariates require the configured constant column %r"
                    % constant_column
                )
            constant_index = header.index(constant_column)
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise error_type(
                    "scDRS covariate line %d has %d fields; expected %d"
                    % (line_number, len(row), len(header))
                )
            cell = row[0].strip()
            if not cell:
                raise error_type(
                    "scDRS covariate line %d has an empty cell identifier"
                    % line_number
                )
            if cell in cells:
                raise error_type(
                    "Duplicate scDRS covariate cell identifier: %s" % cell
                )
            cells.add(cell)
            for column, raw_value in zip(header[1:], row[1:]):
                try:
                    value = float(raw_value)
                except ValueError as exc:
                    raise error_type(
                        "scDRS covariate %r is non-numeric on line %d"
                        % (column, line_number)
                    ) from exc
                if not math.isfinite(value):
                    raise error_type(
                        "scDRS covariate %r is not finite on line %d"
                        % (column, line_number)
                    )
            if constant_index is not None:
                constant_value = float(row[constant_index])
                deviation = abs(constant_value - 1.0)
                maximum_constant_deviation = max(
                    maximum_constant_deviation, deviation,
                )
                if deviation > constant_tolerance:
                    raise error_type(
                        "scDRS covariate constant column %r differs from one "
                        "beyond the configured tolerance on line %d"
                        % (constant_column, line_number)
                    )
            row_count += 1
    return CellCovariateEvidence(
        frozenset(cells), tuple(header[1:]), row_count, constant_column,
        maximum_constant_deviation,
    )


def read_cell_covariates(
    path: str | Path, *, delimiter: str,
    constant_column: str | None, constant_tolerance: float,
    error_type=ValueError,
) -> CellCovariateEvidence:
    """Parse every row; matching these IDs to an atlas remains the consumer's job."""
    source = require_nonempty_file(path, "scDRS covariate file", error_type=error_type)
    contract = {
        "validator": "cell_covariates", "delimiter": delimiter,
        "constant_column": constant_column, "constant_tolerance": constant_tolerance,
    }

    def inspect():
        result = _inspect_cell_covariates(
            source, delimiter=delimiter,
            constant_column=constant_column,
            constant_tolerance=constant_tolerance,
            error_type=error_type,
        )
        record_file_validation(
            source, "Cell covariate table",
            checks=(
                "complete row widths", "unique cell identifiers",
                "finite numeric covariates", "configured intercept constraint",
            ),
            metrics={
                "rows": result.rows, "columns": result.columns,
                "constant_column": result.constant_column,
                "maximum_constant_deviation": result.maximum_constant_deviation,
            },
        )
        return result

    return validate_once((source,), contract, inspect, error_type=error_type)


def normalize_gene_identifier(value: object) -> str | None:
    """Preserve the existing text/whole-number identifier normalization."""
    if value is None or pd.isna(value):
        return None
    if isinstance(value, (float, np.floating)) and float(value).is_integer():
        return str(int(value))
    cleaned = str(value).strip()
    return cleaned or None


def _inspect_identifier_crosswalk(
    path: Path, *, delimiter: str, source_column: str,
    target_column: str, error_type,
) -> tuple[dict[str, str], dict]:
    try:
        frame = pd.read_csv(
            path,
            sep=delimiter,
            dtype=str,
            keep_default_na=False,
        )
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise error_type(
            "Cannot read scDRS gene-identifier crosswalk %s: %s" % (path, exc)
        ) from exc
    required = {
        source_column,
        target_column,
    }
    missing = sorted(required - set(frame.columns))
    if missing or frame.empty:
        raise error_type(
            "scDRS gene-identifier crosswalk is empty or missing columns: %s"
            % ", ".join(missing or sorted(required))
        )
    source = frame[source_column].map(normalize_gene_identifier)
    target = frame[target_column].map(normalize_gene_identifier)
    if source.isna().any() or target.isna().any():
        raise error_type(
            "scDRS gene-identifier crosswalk contains empty source or target IDs"
        )
    pairs = pd.DataFrame({"source": source, "target": target})
    duplicate_sources = pairs["source"].duplicated(keep=False)
    duplicate_targets = pairs["target"].duplicated(keep=False)
    if duplicate_sources.any() or duplicate_targets.any():
        raise error_type(
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


def read_gene_identifier_crosswalk(
    path: str | Path, *, delimiter: str,
    source_column: str, target_column: str, error_type=ValueError,
) -> tuple[dict[str, str], dict]:
    """Read a one-to-one identifier mapping without inferring either namespace."""
    source = require_nonempty_file(
        path, "scDRS gene-identifier crosswalk", error_type=error_type,
    )
    contract = {
        "validator": "one_to_one_gene_crosswalk", "delimiter": delimiter,
        "source_column": source_column, "target_column": target_column,
    }

    def inspect():
        mapping, metrics = _inspect_identifier_crosswalk(
            source, delimiter=delimiter,
            source_column=source_column, target_column=target_column,
            error_type=error_type,
        )
        record_file_validation(
            source, "Gene-identifier crosswalk",
            checks=("required identifier columns", "nonempty identifiers", "one-to-one mapping"),
            metrics=metrics,
        )
        return mapping, metrics

    return deepcopy(validate_once((source,), contract, inspect, error_type=error_type))


def _parse_gene_token(
    token: str, weight_separator: str, *, error_type,
) -> tuple[str, float | None]:
    cleaned = token.strip()
    if not cleaned:
        raise error_type("scDRS gene sets must not contain empty gene entries")
    separator = weight_separator
    if separator not in cleaned:
        return cleaned, None
    gene, raw_weight = cleaned.rsplit(separator, 1)
    gene = gene.strip()
    if not gene or not raw_weight.strip():
        raise error_type("Invalid scDRS gene-weight entry: %r" % cleaned)
    try:
        weight = float(raw_weight)
    except ValueError as exc:
        raise error_type(
            "Invalid scDRS gene weight in entry %r" % cleaned
        ) from exc
    if not math.isfinite(weight):
        raise error_type(
            "scDRS gene weights must be finite: %r" % cleaned
        )
    return gene, weight


def _inspect_gene_sets(
    path: Path, *, trait_column: str, gene_set_column: str,
    delimiter: str, gene_separator: str, weight_separator: str, error_type,
) -> tuple[GeneSetEvidence, ...]:
    """Validate native .gs syntax independently of an atlas or study."""
    gene_set_file = path
    records: list[GeneSetEvidence] = []
    observed_traits: set[str] = set()
    try:
        handle = gene_set_file.open("r", encoding="utf-8", newline="")
    except OSError as exc:
        raise error_type(
            "Cannot read scDRS gene-set file %s: %s" % (gene_set_file, exc)
        ) from exc
    with handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        expected_header = [
            trait_column, gene_set_column,
        ]
        if reader.fieldnames != expected_header:
            raise error_type(
                "scDRS .gs header must be exactly %s"
                % delimiter.join(expected_header)
            )
        for line_number, row in enumerate(reader, start=2):
            trait = validate_filename_component(
                row.get(trait_column, ""),
                "scDRS trait on line %d" % line_number,
                error_type=error_type,
            )
            if "@" in trait:
                raise error_type(
                    "scDRS trait %r contains '@', which is reserved by the "
                    "upstream multi-score filename pattern" % trait
                )
            if trait in observed_traits:
                raise error_type("Duplicate scDRS trait: %s" % trait)
            observed_traits.add(trait)
            raw_gene_set = (row.get(gene_set_column) or "").strip()
            tokens = raw_gene_set.split(gene_separator)
            parsed = [
                _parse_gene_token(token, weight_separator, error_type=error_type)
                for token in tokens
            ]
            genes = [gene for gene, _ in parsed]
            weighted_genes = sum(weight is not None for _, weight in parsed)
            if weighted_genes not in {0, len(parsed)}:
                raise error_type(
                    "scDRS trait %s mixes weighted and unweighted genes; all "
                    "entries must use the same native .gs representation" % trait
                )
            if len(genes) != len(set(genes)):
                raise error_type(
                    "scDRS trait %s contains duplicate genes" % trait
                )
            records.append(GeneSetEvidence(trait, tuple(genes), weighted_genes))
    if not records:
        raise error_type("scDRS .gs file must contain at least one trait")
    return tuple(records)


def read_scdrs_gene_sets(
    path: str | Path, *, trait_column: str, gene_set_column: str,
    delimiter: str, gene_separator: str, weight_separator: str,
    error_type=ValueError,
) -> tuple[GeneSetEvidence, ...]:
    """Parse native .gs syntax; effective atlas overlap is validated downstream."""
    source = require_nonempty_file(path, "scDRS .gs gene-set file", error_type=error_type)
    options = dict(
        trait_column=trait_column, gene_set_column=gene_set_column,
        delimiter=delimiter, gene_separator=gene_separator,
        weight_separator=weight_separator,
    )

    def inspect():
        result = _inspect_gene_sets(source, **options, error_type=error_type)
        record_file_validation(
            source, "scDRS gene-set table",
            checks=(
                "exact native header", "unique trait and gene identifiers",
                "finite weights", "consistent weighted or unweighted representation",
            ),
            metrics={
                "trait_count": len(result),
                "gene_entries": sum(len(row.genes) for row in result),
            },
        )
        return result

    return validate_once(
        (source,), {"validator": "scdrs_gene_sets", **options}, inspect,
        error_type=error_type,
    )
