"""Validated MAGMA gene-property execution."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Mapping, Sequence

from postgwas.core.gene_property_validation import (
    validate_magma_gene_results,
    validate_magma_covariate_table,
)
from postgwas.core.io.reports import write_delimited_report
from postgwas.core.paths import require_nonempty_file
from postgwas.core.processes import run_checked_command
from postgwas.core.statistics import adjust_p_values
from postgwas.modules.magmacovar.errors import MagmaCovarError
from postgwas.modules.magmacovar.stages import MAGMACOVAR_STAGES


def validate_magma_covariate_inputs(
    gene_results_file: str | Path,
    covariates_file: str | Path,
    *,
    minimum_genes: int,
    maximum_missing_fraction: float,
    missing_genes: str,
) -> dict:
    """Validate MAGMA inputs and their gene-ID overlap without a dataframe."""
    gene_results, genes = validate_magma_gene_results(
        gene_results_file, minimum_genes=minimum_genes, error_type=MagmaCovarError,
    )
    summary = validate_magma_covariate_table(
        covariates_file,
        eligible_gene_ids=genes,
        minimum_genes=minimum_genes,
        maximum_missing_fraction=maximum_missing_fraction,
        missing_genes=missing_genes,
        error_type=MagmaCovarError,
    )
    return {
        "gene_results_file": str(gene_results),
        "gene_results_genes": len(genes),
        **summary,
    }


def build_magma_covariate_command(
    magma_bin: str,
    gene_results_file: str | Path,
    covariates_file: str | Path,
    output_prefix: str | Path,
    *,
    model: Sequence[str],
    direction: str,
    missing_values: str,
    maximum_missing_fraction: float,
    missing_genes: str,
) -> list[str]:
    """Build the documented MAGMA command as an argument vector."""
    model_options = []
    for raw_option in model:
        option = str(raw_option).strip()
        if not option or option.startswith("-"):
            raise MagmaCovarError(
                "Each modules.magmacovar.model entry must be a MAGMA --model "
                "modifier without a leading dash"
            )
        model_options.append(option)
    gene_covar = [
        "--gene-covar",
        str(covariates_file),
        "missing-values=%s" % missing_values,
        "max-miss=%.12g" % maximum_missing_fraction,
    ]
    if missing_genes == "fill":
        gene_covar.append("missing-genes=fill")
    return [
        str(magma_bin),
        "--gene-results", str(gene_results_file),
        *gene_covar,
        "--model", *model_options, "direction-covar=%s" % direction,
        "--out", str(output_prefix),
    ]


def read_magma_covariate_results(
    results_file: str | Path, *, minimum_genes: int,
) -> list[dict]:
    """Return validated documented COVAR rows from MAGMA ``.gsa.out`` output."""
    path = require_nonempty_file(
        results_file, "MAGMA gene-property output", error_type=MagmaCovarError,
    )
    header: list[str] | None = None
    rows = []
    required = {"VARIABLE", "TYPE", "NGENES", "BETA", "BETA_STD", "SE", "P"}
    try:
        with path.open("r", encoding="utf-8", errors="strict") as handle:
            for raw_line in handle:
                fields = raw_line.split()
                if not fields:
                    continue
                if header is None and required.issubset(fields):
                    header = fields
                    continue
                if header is None or len(fields) < len(header):
                    continue
                values = dict(zip(header, fields))
                if values.get("TYPE") != "COVAR":
                    continue
                try:
                    genes = int(values["NGENES"])
                    beta = float(values["BETA"])
                    beta_standardized = float(values["BETA_STD"])
                    standard_error = float(values["SE"])
                    p_value = float(values["P"])
                except (KeyError, ValueError) as exc:
                    raise MagmaCovarError(
                        "MAGMA gene-property output contains an invalid COVAR row"
                    ) from exc
                if genes < minimum_genes:
                    raise MagmaCovarError(
                        "MAGMA result %r used %d genes; at least %d are required"
                        % (values.get("VARIABLE", "unknown"), genes, minimum_genes)
                    )
                if (
                    not all(math.isfinite(value) for value in (
                        beta, beta_standardized, standard_error, p_value,
                    ))
                    or standard_error < 0
                    or not 0 <= p_value <= 1
                ):
                    raise MagmaCovarError(
                        "MAGMA result %r contains a non-finite or out-of-range statistic"
                        % values.get("VARIABLE", "unknown")
                    )
                variable = values.get("FULL_NAME") or values["VARIABLE"]
                if not variable.strip():
                    raise MagmaCovarError(
                        "MAGMA gene-property output contains an empty variable name"
                    )
                rows.append(
                    {
                        "variable": variable,
                        "genes": genes,
                        "beta": beta,
                        "beta_standardized": beta_standardized,
                        "standard_error": standard_error,
                        "p_value": p_value,
                    }
                )
    except UnicodeDecodeError as exc:
        raise MagmaCovarError(
            "MAGMA gene-property output is not valid UTF-8 text: %s" % path
        ) from exc
    except OSError as exc:
        raise MagmaCovarError(
            "Cannot read MAGMA gene-property output %s: %s" % (path, exc)
        ) from exc
    if header is None:
        raise MagmaCovarError(
            "MAGMA gene-property output is missing the documented result header"
        )
    if not rows:
        raise MagmaCovarError(
            "MAGMA gene-property output contains no COVAR results"
        )
    variables = [row["variable"] for row in rows]
    if len(variables) != len(set(variables)):
        raise MagmaCovarError(
            "MAGMA gene-property output contains duplicate COVAR variables"
        )
    return rows


def validate_magma_covariate_output(
    results_file: str | Path, *, minimum_genes: int,
) -> dict:
    """Validate and summarize COVAR rows in MAGMA's ``.gsa.out`` output."""
    rows = read_magma_covariate_results(
        results_file, minimum_genes=minimum_genes,
    )
    return {
        "tested_properties": len(rows),
        "minimum_result_genes": min(row["genes"] for row in rows),
        "maximum_result_genes": max(row["genes"] for row in rows),
    }


def _corrected_result_records(
    rows: list[dict], module,
) -> tuple[list[str], list[dict], dict]:
    """Build the configured global correction report from validated COVAR rows."""
    schema = module.result_schema
    correction = module.multiple_testing
    p_values = [row["p_value"] for row in rows]
    try:
        adjusted = {
            method: adjust_p_values(p_values, method)
            for method in correction.methods
        }
    except ValueError as exc:
        raise MagmaCovarError(
            "Cannot correct MAGMA gene-property p-values: %s" % exc
        ) from exc

    adjusted_columns = {
        method: schema.adjusted_p_value_column_pattern.format(method=method)
        for method in correction.methods
    }
    fieldnames = [
        schema.variable_column,
        schema.type_column,
        schema.gene_count_column,
        schema.beta_column,
        schema.standardized_beta_column,
        schema.standard_error_column,
        schema.p_value_column,
        *(adjusted_columns[method] for method in correction.methods),
        schema.primary_correction_method_column,
        schema.primary_adjusted_p_value_column,
        schema.primary_significant_column,
    ]
    primary = correction.primary_method
    records = []
    for index, row in enumerate(rows):
        primary_adjusted = float(adjusted[primary][index])
        record = {
            schema.variable_column: row["variable"],
            schema.type_column: "COVAR",
            schema.gene_count_column: row["genes"],
            schema.beta_column: row["beta"],
            schema.standardized_beta_column: row["beta_standardized"],
            schema.standard_error_column: row["standard_error"],
            schema.p_value_column: row["p_value"],
            schema.primary_correction_method_column: primary,
            schema.primary_adjusted_p_value_column: primary_adjusted,
            schema.primary_significant_column: (
                primary_adjusted <= correction.significance_threshold
            ),
        }
        record.update({
            adjusted_columns[method]: float(adjusted[method][index])
            for method in correction.methods
        })
        records.append(record)

    significant_by_method = {
        method: sum(
            float(adjusted[method][index])
            <= correction.significance_threshold
            for index in range(len(records))
        )
        for method in correction.methods
    }
    ranked_rows = sorted(
        enumerate(rows),
        key=lambda item: (item[1]["p_value"], item[1]["variable"]),
    )[:module.reporting.top_property_count]
    top_properties = [
        {
            "property": row["variable"],
            "genes": row["genes"],
            "beta": row["beta"],
            "standardized_beta": row["beta_standardized"],
            "standard_error": row["standard_error"],
            "p_value": row["p_value"],
            "adjusted_p_values": {
                method: float(adjusted[method][index])
                for method in correction.methods
            },
        }
        for index, row in ranked_rows
    ]
    summary = {
        "tested_properties": len(rows),
        "minimum_result_genes": min(row["genes"] for row in rows),
        "maximum_result_genes": max(row["genes"] for row in rows),
        "primary_correction_method": primary,
        "significance_threshold": correction.significance_threshold,
        "primary_significant_properties": significant_by_method[primary],
        "significant_properties_by_method": significant_by_method,
        "correction_methods": list(correction.methods),
        "top_properties": top_properties,
    }
    return fieldnames, records, summary


def write_corrected_magma_covariate_results(
    results_file: str | Path,
    corrected_results_file: str | Path,
    *,
    module,
    logger,
) -> dict:
    """Write a separate corrected report without modifying MAGMA's native output."""
    rows = read_magma_covariate_results(
        results_file, minimum_genes=module.minimum_genes,
    )
    fieldnames, records, summary = _corrected_result_records(rows, module)
    destination = write_delimited_report(
        records,
        corrected_results_file,
        fieldnames=fieldnames,
        delimiter=module.result_schema.delimiter,
        null_value=module.result_schema.null_value,
    )
    logger.record(
        "RESULT",
        "magmacovar_multiple_testing",
        raw_results_file=str(Path(results_file)),
        corrected_results_file=str(destination),
        **summary,
    )
    return summary


def validate_corrected_magma_covariate_output(
    results_file: str | Path,
    corrected_results_file: str | Path,
    *,
    module,
) -> dict:
    """Prove that a saved corrected report exactly matches its native source."""
    rows = read_magma_covariate_results(
        results_file, minimum_genes=module.minimum_genes,
    )
    fieldnames, expected_records, summary = _corrected_result_records(rows, module)
    corrected = require_nonempty_file(
        corrected_results_file,
        "corrected MAGMA gene-property report",
        error_type=MagmaCovarError,
    )
    try:
        with corrected.open("r", encoding="utf-8", errors="strict", newline="") as handle:
            reader = csv.DictReader(
                handle, delimiter=module.result_schema.delimiter,
            )
            if reader.fieldnames != fieldnames:
                raise MagmaCovarError(
                    "Corrected MAGMA gene-property report columns do not match "
                    "the configured result schema"
                )
            observed_records = list(reader)
    except UnicodeDecodeError as exc:
        raise MagmaCovarError(
            "Corrected MAGMA gene-property report is not valid UTF-8 text: %s"
            % corrected
        ) from exc
    except (OSError, csv.Error) as exc:
        raise MagmaCovarError(
            "Cannot read corrected MAGMA gene-property report %s: %s"
            % (corrected, exc)
        ) from exc

    if len(observed_records) != len(expected_records):
        raise MagmaCovarError(
            "Corrected MAGMA gene-property report has %d rows; %d are required "
            "by the native COVAR results"
            % (len(observed_records), len(expected_records))
        )
    for row_number, (observed, expected) in enumerate(
        zip(observed_records, expected_records), 2,
    ):
        if None in observed:
            raise MagmaCovarError(
                "Corrected MAGMA gene-property report row %d contains more "
                "fields than its configured header" % row_number
            )
        mismatched = [
            field
            for field in fieldnames
            if observed.get(field) != str(expected[field])
        ]
        if mismatched:
            raise MagmaCovarError(
                "Corrected MAGMA gene-property report row %d does not match the "
                "native results and configured corrections in column(s): %s"
                % (row_number, ", ".join(mismatched))
            )
    return summary


def run_magma_covariates(
    *,
    magma_bin: str,
    gene_results_file: str | Path,
    covariates_file: str | Path,
    output_prefix: str | Path,
    results_file: str | Path,
    corrected_results_file: str | Path,
    native_log_file: str | Path,
    module,
    logger,
    timeout_seconds: float | None,
    pipeline_progress=None,
    pipeline_stage_numbers: Mapping[str, int] | None = None,
    finalize_results=None,
) -> dict:
    """Validate inputs, run MAGMA once, and finalize its scientific output."""
    if (pipeline_progress is None) != (pipeline_stage_numbers is None):
        raise MagmaCovarError(
            "MAGMAcovar pipeline progress requires its configured stage mapping"
        )

    def start_pipeline_stage(key: str) -> int | None:
        if pipeline_progress is None:
            return None
        number = int(pipeline_stage_numbers[key])
        pipeline_progress.start(number)
        return number

    def complete_pipeline_stage(number: int | None, fields) -> None:
        if number is not None:
            pipeline_progress.complete(number, outcome_fields=fields)

    total_steps = len(MAGMACOVAR_STAGES)
    pipeline_stage = start_pipeline_stage("gene_results")
    with logger.step(
        1,
        total_steps,
        MAGMACOVAR_STAGES[0],
        "validate_magma_gene_results",
    ) as step:
        gene_results, genes = validate_magma_gene_results(
            gene_results_file, minimum_genes=module.minimum_genes, error_type=MagmaCovarError,
        )
        gene_fields = [
            ("analysis", "MAGMA gene-association results"),
            ("info", "Input file", gene_results.name),
            ("info", "Native format", ".genes.raw"),
            ("genetic", "Eligible gene identifiers", len(genes)),
            ("success", "Required # VERSION metadata", "present"),
            ("success", "File-structure validation", "passed"),
        ]
        step.set_rows(len(genes))
        step.outcome(
            "%s eligible MAGMA gene records passed structural validation."
            % f"{len(genes):,}",
            fields=gene_fields,
            gene_results_file=str(gene_results),
            eligible_genes=len(genes),
            structural_validation="passed",
        )
    complete_pipeline_stage(pipeline_stage, gene_fields)

    pipeline_stage = start_pipeline_stage("covariates")
    with logger.step(
        2,
        total_steps,
        MAGMACOVAR_STAGES[1],
        "validate_magma_covariate_table",
    ) as step:
        covariate_summary = validate_magma_covariate_table(
            covariates_file,
            eligible_gene_ids=genes,
            minimum_genes=module.minimum_genes,
            maximum_missing_fraction=module.input.maximum_missing_fraction,
            missing_genes=module.input.missing_genes,
            error_type=MagmaCovarError,
        )
        overlap = covariate_summary["overlapping_genes"]
        overlap_fraction = overlap / len(genes)
        maximum_observed_missingness = max(
            item["missing_fraction"]
            for item in covariate_summary["property_missingness"]
        )
        covariate_fields = [
            ("analysis", "Gene-property covariate reference"),
            (
                "info", "Input file",
                Path(covariate_summary["covariates_file"]).name,
            ),
            ("count", "Covariate-table genes", covariate_summary["covariate_genes"]),
            ("count", "Gene properties", covariate_summary["properties"]),
            (
                "success", "Gene IDs shared with MAGMA results",
                "%s/%s (%.2f%%)"
                % (f"{overlap:,}", f"{len(genes):,}", overlap_fraction * 100),
            ),
            (
                (
                    "warning"
                    if covariate_summary["gene_results_without_covariates"]
                    else "success"
                ),
                "MAGMA genes without covariates",
                covariate_summary["gene_results_without_covariates"],
            ),
            (
                (
                    "warning"
                    if covariate_summary["covariate_genes_without_gene_results"]
                    else "success"
                ),
                "Covariate genes absent from MAGMA results",
                covariate_summary["covariate_genes_without_gene_results"],
            ),
            ("info", "Missing-value policy", module.input.missing_values),
            ("info", "Missing-gene policy", module.input.missing_genes),
            (
                "info", "Maximum permitted property missingness",
                "%.2f%%" % (module.input.maximum_missing_fraction * 100),
            ),
            (
                "success", "Maximum observed property missingness",
                "%.2f%%" % (maximum_observed_missingness * 100),
            ),
            ("success", "Gene-ID and missingness validation", "passed"),
            ("success", "File-structure validation", "passed"),
        ]
        step.set_rows(overlap)
        step.outcome(
            "%s genes overlap across %s validated gene properties."
            % (f"{overlap:,}", f"{covariate_summary['properties']:,}"),
            fields=covariate_fields,
            **{
                key: value
                for key, value in covariate_summary.items()
                if key not in {"property_names", "property_missingness"}
            },
        )
        logger.record("OBSERVED", "magmacovar_inputs", **{
            "gene_results_file": str(gene_results),
            "gene_results_genes": len(genes),
            **{
                key: value
                for key, value in covariate_summary.items()
                if key != "property_missingness"
            },
        })
        for missingness in covariate_summary["property_missingness"]:
            logger.record(
                "OBSERVED", "magmacovar_property_missingness", **missingness,
            )
    complete_pipeline_stage(pipeline_stage, covariate_fields)

    inputs = {
        "gene_results_file": str(gene_results),
        "gene_results_genes": len(genes),
        **covariate_summary,
    }
    command = build_magma_covariate_command(
        magma_bin,
        inputs["gene_results_file"],
        inputs["covariates_file"],
        output_prefix,
        model=module.model,
        direction=module.direction,
        missing_values=module.input.missing_values,
        maximum_missing_fraction=module.input.maximum_missing_fraction,
        missing_genes=module.input.missing_genes,
    )

    pipeline_stage = start_pipeline_stage("analysis")
    with logger.step(
        3,
        total_steps,
        MAGMACOVAR_STAGES[2],
        "run_magma_gene_property_analysis",
    ) as step:
        run_checked_command(
            command,
            "MAGMA gene-property analysis",
            logger=logger,
            error_type=MagmaCovarError,
            timeout_seconds=timeout_seconds,
            expected_outputs=[results_file, native_log_file],
        )
        analysis_fields = [
            ("analysis", "MAGMA gene-property model"),
            (
                "decision", "Model modifiers",
                ", ".join(module.model) if module.model else "marginal",
            ),
            ("decision", "Property-test direction", module.direction),
            ("info", "Missing-value policy", module.input.missing_values),
            ("info", "Missing-gene policy", module.input.missing_genes),
            ("success", "MAGMA execution", "completed with exit status 0"),
            ("info", "Native result", Path(results_file).name),
        ]
        step.outcome(
            "MAGMA completed the configured gene-property model.",
            fields=analysis_fields,
            model=list(module.model),
            direction=module.direction,
            results_file=str(results_file),
        )
    complete_pipeline_stage(pipeline_stage, analysis_fields)

    pipeline_stage = start_pipeline_stage("results")
    with logger.step(
        4,
        total_steps,
        MAGMACOVAR_STAGES[3],
        "validate_and_correct_magma_covariate_results",
    ) as step:
        output = write_corrected_magma_covariate_results(
            results_file,
            corrected_results_file,
            module=module,
            logger=logger,
        )
        correction_label = module.multiple_testing.reporting_method_labels[
            output["primary_correction_method"]
        ]
        result_fields = [
            ("analysis", "Gene-property results"),
            ("count", "Tested properties", output["tested_properties"]),
            (
                "genetic", "Genes per tested property",
                "%s–%s"
                % (
                    f"{output['minimum_result_genes']:,}",
                    f"{output['maximum_result_genes']:,}",
                ),
            ),
            ("decision", "Primary correction", correction_label),
            (
                "info", "Primary significance threshold",
                output["significance_threshold"],
            ),
            (
                (
                    "success"
                    if output["primary_significant_properties"] else "info"
                ),
                "Significant after primary correction",
                "%s/%s"
                % (
                    f"{output['primary_significant_properties']:,}",
                    f"{output['tested_properties']:,}",
                ),
            ),
        ]
        highlight_method = module.reporting.highlight_method
        highlight_significant = output["significant_properties_by_method"].get(
            highlight_method
        )
        if (
            highlight_significant is not None
            and highlight_method != output["primary_correction_method"]
        ):
            highlight_label = module.multiple_testing.reporting_method_labels[
                highlight_method
            ]
            result_fields.append((
                "success" if highlight_significant else "info",
                "%s significant properties" % highlight_label,
                "%s/%s"
                % (
                    f"{highlight_significant:,}",
                    f"{output['tested_properties']:,}",
                ),
            ))
        result_fields.extend((
            ("success", "Native COVAR-result validation", "passed"),
            (
                "success", "Corrected result report",
                Path(corrected_results_file).name,
            ),
        ))
        finalized = (
            finalize_results(inputs, output)
            if finalize_results is not None
            else None
        )
        step.set_rows(output["tested_properties"])
        step.outcome(
            "%s COVAR results were validated and corrected as one family."
            % f"{output['tested_properties']:,}",
            fields=result_fields,
            **output,
        )
    return {
        "inputs": inputs,
        "output": output,
        "command": command,
        "finalized": finalized,
        "result_fields": result_fields,
    }


__all__ = [
    "build_magma_covariate_command",
    "read_magma_covariate_results",
    "run_magma_covariates",
    "validate_corrected_magma_covariate_output",
    "validate_magma_covariate_inputs",
    "validate_magma_gene_results",
    "validate_magma_covariate_output",
    "validate_magma_covariate_table",
    "write_corrected_magma_covariate_results",
]
