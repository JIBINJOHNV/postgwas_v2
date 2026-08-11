"""Scientific validation and normalization for MAGMA cell typing."""

from __future__ import annotations

from pathlib import Path

from postgwas.core.io.reports import write_delimited_report
from postgwas.core.statistics import adjust_p_values
from postgwas.modules.magmacovar.main import (
    read_magma_covariate_results,
    validate_magma_covariate_table,
)
from postgwas.modules.single_cell.errors import SingleCellError


def validate_magma_celltype_use_case(single_cell_config, magmacovar_config):
    """Resolve and enforce the published FUMA base-model MAGMA contract."""
    method = single_cell_config.magma_celltype
    use_case = magmacovar_config.model_use_cases.get(method.magmacovar_use_case)
    if use_case is None:
        raise SingleCellError(
            "modules.single_cell.magma_celltype.magmacovar_use_case %r is not "
            "defined in modules.magmacovar.model_use_cases"
            % method.magmacovar_use_case
        )
    expected_model = ["condition-hide=%s" % method.average_property]
    if use_case.model != expected_model or use_case.direction != "greater":
        raise SingleCellError(
            "MAGMA cell typing base analysis requires the selected MAGMAcovar "
            "use case to specify model %r and direction 'greater'; observed model "
            "%r and direction %r"
            % (expected_model, use_case.model, use_case.direction)
        )
    return use_case


def validate_magma_celltype_covariates(
    covariates_file: str | Path,
    *,
    average_property: str,
    magmacovar_config,
) -> dict:
    """Validate the expression-property matrix and its hypothesis family."""
    summary = validate_magma_covariate_table(
        covariates_file,
        minimum_genes=magmacovar_config.minimum_genes,
        maximum_missing_fraction=(
            magmacovar_config.input.maximum_missing_fraction
        ),
        missing_genes=magmacovar_config.input.missing_genes,
    )
    properties = list(summary["property_names"])
    if average_property not in properties:
        raise SingleCellError(
            "MAGMA cell-type covariates must contain the configured average "
            "expression property %r" % average_property
        )
    cell_types = [name for name in properties if name != average_property]
    if not cell_types:
        raise SingleCellError(
            "MAGMA cell-type covariates must contain at least one cell-type "
            "property in addition to %r" % average_property
        )
    return {**summary, "cell_type_properties": cell_types}


def normalize_magma_celltype_results(
    gene_property_results: str | Path,
    covariates_file: str | Path,
    output_file: str | Path,
    *,
    dataset_id: str,
    single_cell_config,
    magmacovar_config,
) -> dict:
    """Validate one FUMA-compatible base analysis and write corrected results."""
    method = single_cell_config.magma_celltype
    covariates = validate_magma_celltype_covariates(
        covariates_file,
        average_property=method.average_property,
        magmacovar_config=magmacovar_config,
    )
    rows = read_magma_covariate_results(
        gene_property_results,
        minimum_genes=magmacovar_config.minimum_genes,
    )
    expected = set(covariates["cell_type_properties"])
    observed = {str(row["variable"]) for row in rows}
    missing = sorted(expected - observed)
    unexpected = sorted(observed - expected)
    if missing or unexpected:
        detail = []
        if missing:
            detail.append("missing cell types: %s" % ", ".join(missing))
        if unexpected:
            detail.append("unexpected properties: %s" % ", ".join(unexpected))
        raise SingleCellError(
            "MAGMA gene-property results do not match the configured cell-type "
            "hypothesis family (%s)" % "; ".join(detail)
        )

    multiple_testing = single_cell_config.multiple_testing
    p_values = [float(row["p_value"]) for row in rows]
    adjusted = {
        correction: adjust_p_values(p_values, correction).tolist()
        for correction in multiple_testing.methods
    }
    schema = single_cell_config.result_schema
    records = []
    for index, row in enumerate(rows):
        record = {
            schema.dataset_column: dataset_id,
            schema.method_column: "magma_celltype",
            schema.workflow_column: method.workflow,
            schema.cell_type_column: row["variable"],
            schema.gene_count_column: int(row["genes"]),
            schema.beta_column: float(row["beta"]),
            schema.standardized_beta_column: float(row["beta_standardized"]),
            schema.standard_error_column: float(row["standard_error"]),
            schema.p_value_column: float(row["p_value"]),
        }
        for correction in multiple_testing.methods:
            adjusted_value = float(adjusted[correction][index])
            record[
                schema.adjusted_p_value_column_pattern.format(method=correction)
            ] = adjusted_value
            record[
                schema.significance_column_pattern.format(method=correction)
            ] = bool(adjusted_value <= multiple_testing.significance_threshold)
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
        fieldnames.extend(
            (
                schema.adjusted_p_value_column_pattern.format(method=correction),
                schema.significance_column_pattern.format(method=correction),
            )
        )
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
            correction: int(
                sum(value <= threshold for value in adjusted[correction])
            )
            for correction in multiple_testing.methods
        },
        "primary_correction": multiple_testing.primary_method,
        "significance_threshold": threshold,
        "covariate_genes": covariates["covariate_genes"],
        "average_property": method.average_property,
        "normalized_schema_version": schema.normalized_schema_version,
    }


__all__ = [
    "normalize_magma_celltype_results",
    "validate_magma_celltype_covariates",
    "validate_magma_celltype_use_case",
]
