"""Validation and normalization of native LDSC ``--h2-cts`` results."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any

from postgwas.core.io.reports import write_delimited_report
from postgwas.core.statistics import adjust_p_values
from postgwas.modules.single_cell.errors import SingleCellError


def normalize_ldsc_celltype_results(
    native_results: str | Path,
    output_file: str | Path,
    *,
    expected_cell_types: tuple[str, ...],
    dataset_id: str,
    single_cell_config,
) -> dict[str, Any]:
    """Validate LDSC output and write the shared single-cell result schema."""
    method = single_cell_config.ldsc_celltype
    native = method.result_format
    expected_header = [
        native.name_column,
        native.coefficient_column,
        native.standard_error_column,
        native.p_value_column,
    ]
    try:
        handle = Path(native_results).open("r", encoding="utf-8", newline="")
    except OSError as exc:
        raise SingleCellError(
            "Cannot read LDSC cell-type results %s: %s" % (native_results, exc)
        ) from exc
    rows: list[dict[str, Any]] = []
    observed: set[str] = set()
    with handle:
        reader = csv.DictReader(handle, delimiter=native.delimiter)
        if reader.fieldnames != expected_header:
            raise SingleCellError(
                "LDSC cell-type result header must be exactly %s"
                % native.delimiter.join(expected_header)
            )
        for line_number, row in enumerate(reader, start=2):
            name = (row.get(native.name_column) or "").strip()
            if not name:
                raise SingleCellError(
                    "LDSC cell-type result line %d has an empty name" % line_number
                )
            if name in observed:
                raise SingleCellError(
                    "Duplicate LDSC cell-type result name: %s" % name
                )
            observed.add(name)
            values: dict[str, float] = {}
            for key, label in (
                (native.coefficient_column, "coefficient"),
                (native.standard_error_column, "standard error"),
                (native.p_value_column, "P value"),
            ):
                try:
                    value = float(row.get(key, ""))
                except (TypeError, ValueError) as exc:
                    raise SingleCellError(
                        "LDSC %s is non-numeric on result line %d"
                        % (label, line_number)
                    ) from exc
                if not math.isfinite(value):
                    raise SingleCellError(
                        "LDSC %s is not finite on result line %d"
                        % (label, line_number)
                    )
                values[key] = value
            if values[native.standard_error_column] <= 0:
                raise SingleCellError(
                    "LDSC standard error must be positive on result line %d"
                    % line_number
                )
            p_value = values[native.p_value_column]
            if not 0 <= p_value <= 1:
                raise SingleCellError(
                    "LDSC P value must be within [0, 1] on result line %d"
                    % line_number
                )
            rows.append({"name": name, **values})
    if not rows:
        raise SingleCellError("LDSC cell-type result file contains no tests")
    expected = set(expected_cell_types)
    missing = sorted(expected - observed)
    unexpected = sorted(observed - expected)
    if missing or unexpected:
        details = []
        if missing:
            details.append("missing: %s" % ", ".join(missing))
        if unexpected:
            details.append("unexpected: %s" % ", ".join(unexpected))
        raise SingleCellError(
            "LDSC results do not match the .ldcts hypothesis family (%s)"
            % "; ".join(details)
        )

    multiple_testing = single_cell_config.multiple_testing
    p_values = [row[native.p_value_column] for row in rows]
    adjusted = {
        correction: adjust_p_values(p_values, correction).tolist()
        for correction in multiple_testing.methods
    }
    schema = single_cell_config.result_schema
    records = []
    for index, row in enumerate(rows):
        record = {
            schema.dataset_column: dataset_id,
            schema.method_column: "ldsc_celltype",
            schema.workflow_column: method.workflow,
            schema.cell_type_column: row["name"],
            schema.gene_count_column: None,
            schema.beta_column: row[native.coefficient_column],
            schema.standardized_beta_column: None,
            schema.standard_error_column: row[native.standard_error_column],
            schema.p_value_column: row[native.p_value_column],
        }
        for correction in multiple_testing.methods:
            value = float(adjusted[correction][index])
            record[
                schema.adjusted_p_value_column_pattern.format(method=correction)
            ] = value
            record[
                schema.significance_column_pattern.format(method=correction)
            ] = bool(value <= multiple_testing.significance_threshold)
        records.append(record)
    fieldnames = [
        schema.dataset_column,
        schema.method_column,
        schema.workflow_column,
        schema.cell_type_column,
        schema.gene_count_column,
        schema.beta_column,
        schema.standardized_beta_column,
        schema.standard_error_column,
        schema.p_value_column,
    ]
    for correction in multiple_testing.methods:
        fieldnames.extend((
            schema.adjusted_p_value_column_pattern.format(method=correction),
            schema.significance_column_pattern.format(method=correction),
        ))
    destination = write_delimited_report(
        records,
        output_file,
        fieldnames=fieldnames,
        delimiter=schema.delimiter,
        null_value=schema.null_value,
    )
    threshold = multiple_testing.significance_threshold
    return {
        "output": str(destination),
        "workflow": method.workflow,
        "tested_cell_types": len(records),
        "nominal_significant": sum(value <= threshold for value in p_values),
        "adjusted_significant": {
            correction: int(sum(value <= threshold for value in values))
            for correction, values in adjusted.items()
        },
        "primary_correction": multiple_testing.primary_method,
        "significance_threshold": threshold,
        "normalized_schema_version": schema.normalized_schema_version,
        "coefficient_interpretation": "additional_per_snp_heritability",
        "native_p_value_alternative": "coefficient_greater_than_zero",
    }


__all__ = ["normalize_ldsc_celltype_results"]
