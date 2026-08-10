"""Validation, normalization, and scientific summaries for GCTA-COJO output."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from postgwas.core.io.tables import read_delimited_table
from postgwas.modules.gcta_cojo.errors import GctaCojoError


def _atomic_write(frame: pl.DataFrame, destination: Path, config) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(
        "%s%s" % (destination.name, config.atomic_output_suffix)
    )
    try:
        frame.write_csv(
            temporary,
            separator=config.normalized_delimiter,
            null_value=config.normalized_null_value,
        )
        temporary.replace(destination)
    except (OSError, pl.exceptions.PolarsError) as exc:
        raise GctaCojoError(
            "Cannot write normalized GCTA-COJO result %s: %s"
            % (destination, exc)
        ) from exc
    finally:
        temporary.unlink(missing_ok=True)


def normalize_cojo_results(
    source: str | Path,
    destination: str | Path,
    module,
) -> tuple[dict, pl.DataFrame]:
    """Validate GCTA's result schema and preserve documented NA estimates."""
    result_config = module.results
    schema = result_config.schemas[module.mode]
    frame, delimiter = read_delimited_table(
        source,
        result_config.delimiter,
        candidates=result_config.delimiter_candidates,
        minimum_columns=len(schema.required_columns),
        maximum_columns=result_config.maximum_columns,
        sample_lines=result_config.sample_lines,
        null_values=result_config.null_values,
        infer_schema_length=result_config.infer_schema_length,
        error_type=GctaCojoError,
        description="GCTA-COJO %s result" % module.mode,
    )
    missing = [column for column in schema.required_columns if column not in frame.columns]
    if missing:
        raise GctaCojoError(
            "GCTA-COJO %s result is missing required columns: %s"
            % (module.mode, ", ".join(missing))
        )
    if frame.is_empty():
        raise GctaCojoError(
            "GCTA-COJO %s produced an empty result; inspect the GCTA log and "
            "reference overlap diagnostics." % module.mode
        )
    if frame.select(
        pl.col(schema.identifier_column).cast(pl.String).is_duplicated().any()
    ).item():
        raise GctaCojoError("GCTA-COJO result contains duplicate SNP identifiers.")

    converted = {
        column: pl.col(column).cast(pl.Float64, strict=False)
        for column in schema.numeric_columns
    }
    for column, expression in converted.items():
        invalid_text = frame.filter(
            pl.col(column).is_not_null() & expression.is_null()
        ).height
        if invalid_text:
            raise GctaCojoError(
                "GCTA-COJO result column %s contains %d non-numeric values."
                % (column, invalid_text)
            )
    frame = frame.with_columns(*converted.values())
    invalid_base = frame.filter(pl.any_horizontal([
        pl.col(column).is_null() | ~pl.col(column).is_finite()
        for column in schema.required_numeric_columns
    ])).height
    if invalid_base:
        raise GctaCojoError(
            "GCTA-COJO result contains %d rows with missing, non-numeric, or "
            "non-finite marginal statistics." % invalid_base
        )
    for column in schema.positive_columns:
        invalid = frame.filter(pl.col(column) <= 0).height
        if invalid:
            raise GctaCojoError(
                "GCTA-COJO result column %s contains %d non-positive values."
                % (column, invalid)
            )
    for column in schema.integer_columns:
        invalid = frame.filter(pl.col(column) != pl.col(column).floor()).height
        if invalid:
            raise GctaCojoError(
                "GCTA-COJO result column %s contains %d non-integer values."
                % (column, invalid)
            )
    for column in schema.closed_unit_interval_columns:
        invalid = frame.filter(
            (pl.col(column) < 0) | (pl.col(column) > 1)
        ).height
        if invalid:
            raise GctaCojoError(
                "GCTA-COJO result column %s contains %d values outside [0, 1]."
                % (column, invalid)
            )
    p_columns = [schema.marginal_p_value_column, schema.p_value_column]
    for column in p_columns:
        invalid = frame.filter(
            pl.col(column).is_not_null()
            & (
                ~pl.col(column).is_finite()
                | (pl.col(column) < 0)
                | (pl.col(column) > 1)
            )
        ).height
        if invalid:
            raise GctaCojoError(
                "GCTA-COJO result column %s contains %d invalid p-values."
                % (column, invalid)
            )

    adjusted = [
        schema.adjusted_effect_column,
        schema.adjusted_standard_error_column,
        schema.p_value_column,
    ]
    null_counts = [
        frame.select(pl.col(column).is_null().sum()).item() for column in adjusted
    ]
    if len(set(null_counts)) != 1:
        raise GctaCojoError(
            "GCTA-COJO adjusted effect, standard error, and p-value columns have "
            "inconsistent missingness."
        )
    not_estimable = int(null_counts[0])
    if not_estimable and module.mode != "cond":
        raise GctaCojoError(
            "GCTA-COJO %s result contains %d missing joint estimates."
            % (module.mode, not_estimable)
        )
    invalid_se = frame.filter(
        pl.col(schema.adjusted_standard_error_column).is_not_null()
        & (
            ~pl.col(schema.adjusted_standard_error_column).is_finite()
            | (pl.col(schema.adjusted_standard_error_column) <= 0)
        )
    ).height
    if invalid_se:
        raise GctaCojoError(
            "GCTA-COJO result contains %d non-positive or non-finite adjusted "
            "standard errors." % invalid_se
        )

    status = result_config.estimation_status_column
    if status in frame.columns:
        raise GctaCojoError(
            "Configured estimation status column conflicts with GCTA output: %s"
            % status
        )
    frame = frame.with_columns(
        pl.when(pl.col(schema.p_value_column).is_null())
        .then(pl.lit(result_config.not_estimable_status))
        .otherwise(pl.lit(result_config.estimated_status))
        .alias(status)
    )
    destination = Path(destination)
    _atomic_write(frame, destination, result_config)

    threshold = module.reporting.finding_threshold
    estimated = frame.filter(pl.col(schema.p_value_column).is_not_null())
    significant = estimated.filter(pl.col(schema.p_value_column) <= threshold)
    newly_significant = 0
    if module.mode == "cond":
        newly_significant = significant.filter(
            pl.col(schema.marginal_p_value_column) > threshold
        ).height
    top_rows = (
        estimated.sort(schema.p_value_column)
        .head(module.reporting.top_result_count)
        .select(
            schema.identifier_column,
            schema.chromosome_column,
            schema.position_column,
            schema.marginal_p_value_column,
            schema.p_value_column,
        )
        .to_dicts()
    )
    top = [{
        "variant_id": row[schema.identifier_column],
        "chromosome": row[schema.chromosome_column],
        "position": row[schema.position_column],
        "marginal_p_value": row[schema.marginal_p_value_column],
        "adjusted_p_value": row[schema.p_value_column],
    } for row in top_rows]
    metrics = {
        "result_rows": frame.height,
        "estimated_rows": estimated.height,
        "not_estimable_rows": not_estimable,
        "significant_rows": significant.height,
        "newly_significant_rows": newly_significant,
        "finding_threshold": threshold,
        "adjusted_p_value_column": schema.p_value_column,
        "delimiter": delimiter.method,
        "normalized_result": str(destination),
        "top_findings": top,
    }
    return metrics, frame


__all__ = ["normalize_cojo_results"]
