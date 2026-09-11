"""Validation and normalization of official GCTA gene-test results."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from postgwas.core.io.reports import write_delimited_report
from postgwas.core.io.tables import read_delimited_table
from postgwas.core.statistics import adjust_p_values
from postgwas.modules.gcta_gene.errors import GctaGeneError


@dataclass(frozen=True)
class ValidatedGctaResults:
    """Parsed raw GCTA results whose schema and scientific values passed."""

    frame: pl.DataFrame
    p_column: str
    primary_p_values: tuple[float, ...]
    input_delimiter: str
    unit_label: str


def multiple_testing_configuration(module_config) -> dict:
    """Return only settings that determine corrected result-table values."""
    reporting = module_config.reporting
    return {
        "nominal_alpha": reporting.nominal_alpha,
        "familywise_alpha": reporting.familywise_alpha,
        "fdr_alpha": reporting.fdr_alpha,
        "correction_columns": reporting.correction_columns.model_dump(
            mode="json"
        ),
        "normalized_schema_version": (
            module_config.results.normalized_schema_version
        ),
    }


def _add_multiple_testing_columns(
    frame: pl.DataFrame,
    p_values: list[float],
    module_config,
) -> tuple[pl.DataFrame, dict]:
    reporting = module_config.reporting
    columns = reporting.correction_columns
    output_names = list(columns.model_dump().values())
    conflicts = [name for name in output_names if name in frame.columns]
    if conflicts:
        raise GctaGeneError(
            "GCTA result already contains configured multiple-testing columns: %s"
            % ", ".join(conflicts)
        )
    family_size = len(p_values)
    bonferroni = adjust_p_values(p_values, "bonferroni").tolist()
    fdr_bh = adjust_p_values(p_values, "fdr_bh").tolist()
    nominal_significant = [
        value <= reporting.nominal_alpha for value in p_values
    ]
    bonferroni_significant = [
        value <= reporting.familywise_alpha for value in bonferroni
    ]
    fdr_bh_significant = [
        value <= reporting.fdr_alpha for value in fdr_bh
    ]
    corrected = frame.with_columns(
        pl.Series(columns.nominal_significant, nominal_significant),
        pl.Series(columns.bonferroni_adjusted_p, bonferroni),
        pl.Series(columns.bonferroni_significant, bonferroni_significant),
        pl.Series(columns.fdr_bh_adjusted_p, fdr_bh),
        pl.Series(columns.fdr_bh_significant, fdr_bh_significant),
    )
    return corrected, {
        "family_size": family_size,
        "nominal_significant": sum(nominal_significant),
        "bonferroni_significant": sum(bonferroni_significant),
        "fdr_bh_significant": sum(fdr_bh_significant),
        "configuration": multiple_testing_configuration(module_config),
    }


def validate_raw_gcta_results(
    raw_result: str | Path,
    module_config,
) -> ValidatedGctaResults:
    """Validate required upstream columns and values without changing results."""
    settings = module_config.results
    schema = settings.schemas[module_config.method]
    required = list(schema.required_columns)
    if module_config.method == "mbat_combo" and module_config.print_component_p_values:
        required.extend(schema.component_p_value_columns)
    p_column = schema.p_value_column
    frame, detected = read_delimited_table(
        raw_result,
        settings.delimiter,
        candidates=settings.delimiter_candidates,
        minimum_columns=len(required),
        maximum_columns=settings.maximum_columns,
        sample_lines=settings.sample_lines,
        null_values=settings.null_values,
        infer_schema_length=settings.infer_schema_length,
        error_type=GctaGeneError,
        description="GCTA %s result" % module_config.method,
    )
    if frame.is_empty():
        raise GctaGeneError("GCTA result contains no tested units: %s" % raw_result)
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise GctaGeneError(
            "GCTA %s result is missing required columns: %s"
            % (module_config.method, ", ".join(missing))
        )
    incomplete_rows = frame.filter(
        pl.any_horizontal(pl.col(required).is_null())
    ).height
    if incomplete_rows:
        raise GctaGeneError(
            "GCTA %s result contains %d rows missing required values."
            % (module_config.method, incomplete_rows)
        )
    p_columns = [p_column]
    if module_config.method == "mbat_combo" and module_config.print_component_p_values:
        p_columns.extend(schema.component_p_value_columns)
    parsed_p_values = {}
    for column in p_columns:
        values = frame[column].cast(pl.Float64, strict=False)
        invalid_p = sum(
            value is None or not math.isfinite(value) or value < 0 or value > 1
            for value in values.to_list()
        )
        if invalid_p:
            raise GctaGeneError(
                "GCTA result contains %d invalid %s values." % (invalid_p, column)
            )
        parsed_p_values[column] = values
    empty_identifiers = frame.filter(pl.any_horizontal(*[
        pl.col(column).cast(pl.String).str.strip_chars() == ""
        for column in schema.identifier_columns
    ])).height
    if empty_identifiers:
        raise GctaGeneError(
            "GCTA result contains %d empty unit identifiers." % empty_identifiers
        )
    duplicate_units = frame.select(schema.identifier_columns).is_duplicated().sum()
    if duplicate_units:
        raise GctaGeneError(
            "GCTA result contains %d duplicate unit identifiers." % duplicate_units
        )
    primary_p_values = tuple(
        float(value) for value in parsed_p_values[p_column].to_list()
    )
    return ValidatedGctaResults(
        frame=frame,
        p_column=p_column,
        primary_p_values=primary_p_values,
        input_delimiter=detected.value,
        unit_label=schema.unit_label,
    )


def add_multiple_testing_results(
    validated: ValidatedGctaResults,
    normalized_result: str | Path,
    module_config,
) -> dict:
    """Add configured correction columns and write the validated normalized TSV."""
    settings = module_config.results
    frame, multiple_testing = _add_multiple_testing_columns(
        validated.frame,
        list(validated.primary_p_values),
        module_config,
    )
    destination = write_delimited_report(
        frame.to_dicts(),
        normalized_result,
        fieldnames=frame.columns,
        delimiter=settings.normalized_delimiter,
        null_value=settings.normalized_null_value,
    )
    minimum_p = min(validated.primary_p_values)
    return {
        "tested_units": frame.height,
        "unit_label": validated.unit_label,
        "minimum_p_value": minimum_p,
        "p_value_column": validated.p_column,
        "input_delimiter": validated.input_delimiter,
        "normalized_result": str(destination),
        "columns": frame.columns,
        "multiple_testing": multiple_testing,
    }


def normalize_gcta_results(
    raw_result: str | Path,
    normalized_result: str | Path,
    module_config,
) -> dict:
    """Validate raw results, add corrections, and write a normalized TSV."""
    validated = validate_raw_gcta_results(raw_result, module_config)
    return add_multiple_testing_results(
        validated, normalized_result, module_config,
    )


__all__ = [
    "ValidatedGctaResults",
    "add_multiple_testing_results",
    "multiple_testing_configuration",
    "normalize_gcta_results",
    "validate_raw_gcta_results",
]
