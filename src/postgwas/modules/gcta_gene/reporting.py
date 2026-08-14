"""Scientific summaries for validated GCTA fastBAT and mBAT-combo results."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

import polars as pl

from postgwas.core.contracts import Artifact
from postgwas.core.dataframes import chromosome_expression
from postgwas.core.pipeline_logging import PipelineLogger
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.modules.gcta_gene.errors import GctaGeneError


_METHOD_TITLES = {
    "fastbat_gene": "GCTA fastBAT gene-association summary",
    "fastbat_segment": "GCTA fastBAT segment-association summary",
    "fastbat_set": "GCTA fastBAT set-association summary",
    "mbat_combo": "GCTA mBAT-combo gene-association summary",
}


def _association_identifier(row: dict, columns: list[str]) -> str:
    if len(columns) == 1:
        return str(row[columns[0]])
    return " · ".join("%s %s" % (column, row[column]) for column in columns)


def _format_p_value(value: float, significant_digits: int) -> str:
    return format(float(value), ".%dg" % significant_digits)


def build_gcta_scientific_summary(
    normalized_result: str | Path,
    module_config,
    metrics: Mapping,
) -> dict:
    """Rank validated associations and derive explicitly reporting-only QC."""
    path = Path(normalized_result)
    try:
        frame = pl.read_csv(
            path,
            separator=module_config.results.normalized_delimiter,
            null_values=module_config.results.null_values,
            infer_schema_length=module_config.results.infer_schema_length,
        )
    except (OSError, pl.exceptions.PolarsError) as exc:
        raise GctaGeneError(
            "Cannot read normalized GCTA results for reporting %s: %s"
            % (path, exc)
        ) from exc
    schema = module_config.results.schemas[module_config.method]
    correction_columns = module_config.reporting.correction_columns
    required = [
        schema.p_value_column,
        *schema.identifier_columns,
        *correction_columns.model_dump().values(),
    ]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise GctaGeneError(
            "Normalized GCTA result is missing reporting columns: %s"
            % ", ".join(missing)
        )
    if frame.is_empty():
        raise GctaGeneError("Normalized GCTA result contains no associations: %s" % path)

    p_column = schema.p_value_column
    p_values = frame[p_column].cast(pl.Float64, strict=False)
    if p_values.null_count():
        raise GctaGeneError(
            "Normalized GCTA result contains invalid %s values." % p_column
        )
    frame = frame.with_columns(p_values.alias(p_column))
    numeric_corrections = (
        correction_columns.bonferroni_adjusted_p,
        correction_columns.fdr_bh_adjusted_p,
    )
    boolean_corrections = (
        correction_columns.nominal_significant,
        correction_columns.bonferroni_significant,
        correction_columns.fdr_bh_significant,
    )
    for column in numeric_corrections:
        values = frame[column].cast(pl.Float64, strict=False)
        invalid = frame.filter(
            values.is_null() | ~values.is_finite() | (values < 0) | (values > 1)
        ).height
        if invalid:
            raise GctaGeneError(
                "Normalized GCTA result contains %d invalid %s values."
                % (invalid, column)
            )
        frame = frame.with_columns(values.alias(column))
    for column in boolean_corrections:
        values = frame[column].cast(pl.Boolean, strict=False)
        if values.null_count():
            raise GctaGeneError(
                "Normalized GCTA result contains invalid %s flags." % column
            )
        frame = frame.with_columns(values.alias(column))
    ranked = frame.sort(p_column, maintain_order=True).head(
        module_config.reporting.top_result_count
    )
    component_columns = [
        column for column in schema.component_p_value_columns
        if column in ranked.columns
    ]
    top_associations = []
    for rank, row in enumerate(ranked.to_dicts(), 1):
        top_associations.append({
            "rank": rank,
            "identifier": _association_identifier(row, schema.identifier_columns),
            "p_value": float(row[p_column]),
            "bonferroni_adjusted_p": float(
                row[correction_columns.bonferroni_adjusted_p]
            ),
            "fdr_bh_adjusted_p": float(
                row[correction_columns.fdr_bh_adjusted_p]
            ),
            "component_p_values": {
                column: float(row[column]) for column in component_columns
            },
        })

    tested_units = int(metrics["tested_units"])
    familywise_alpha = module_config.reporting.familywise_alpha
    bonferroni_threshold = familywise_alpha / tested_units
    nominal_significant = frame.filter(
        pl.col(correction_columns.nominal_significant)
    ).height
    bonferroni_significant = frame.filter(
        pl.col(correction_columns.bonferroni_significant)
    ).height
    fdr_bh_significant = frame.filter(
        pl.col(correction_columns.fdr_bh_significant)
    ).height
    input_variants = int(metrics.get("variants", 0))
    resolved_variants = int(
        metrics.get("overlapping_variants", metrics.get("harmonised_input_variants", 0))
    )
    source_units = None
    source_units_label = None
    if module_config.method in {"fastbat_gene", "mbat_combo"}:
        source_units = metrics.get("genes")
        source_units_label = "Gene annotations supplied"
    elif module_config.method == "fastbat_set":
        source_units = metrics.get("input_sets")
        source_units_label = "Custom sets supplied"
    if source_units is not None:
        source_units = int(source_units)

    expected_chromosomes = []
    result_chromosomes = []
    chromosome_column = module_config.reporting.chromosome_columns[
        module_config.method
    ]
    if chromosome_column is not None:
        expected_source = (
            metrics.get("shared_chromosomes")
            if module_config.method in {"fastbat_gene", "mbat_combo"}
            else metrics.get("reference_chromosomes")
        )
        expected_chromosomes = [str(value) for value in (expected_source or [])]
        normalized = frame.select(chromosome_expression(
            chromosome_column,
            frame.schema[chromosome_column],
            strip_chr_prefix=(
                module_config.variant_harmonisation.chromosome_label_policy
                == "strip_chr_prefix"
            ),
        ))
        result_chromosomes = sorted(
            set(normalized[chromosome_column].drop_nulls().to_list())
        )
    result_chromosome_set = set(result_chromosomes)
    missing_result_chromosomes = [
        chromosome for chromosome in expected_chromosomes
        if chromosome not in result_chromosome_set
    ]
    warnings = []
    unresolved = int(metrics.get("unresolved_variants_removed", 0))
    if unresolved:
        warnings.append(
            "%s GWAS variants could not be resolved to exact PLINK BIM IDs and "
            "were excluded before GCTA."
            % format(unresolved, ",")
        )
    direct_input_unmodified = bool(metrics.get("direct_input_unmodified", False))
    absent_from_reference = int(
        metrics.get("input_variants_absent_from_reference", 0)
    )
    if direct_input_unmodified and absent_from_reference:
        variant_label = (
            "variant is" if absent_from_reference == 1 else "variants are"
        )
        warnings.append(
            "%s direct-input GWAS %s absent from PLINK BIM column 2 "
            "and will not be used by GCTA. PostGWAS passed the original "
            "summary-statistics input without rewriting or filtering it."
            % (format(absent_from_reference, ","), variant_label)
        )
    omitted_empty = int(metrics.get("omitted_empty_sets", 0))
    if omitted_empty:
        warnings.append(
            "%s custom sets had no variants shared by the GWAS and LD reference "
            "and were omitted under the configured policy."
            % format(omitted_empty, ",")
        )
    omitted_oversized = int(metrics.get("omitted_oversized_sets", 0))
    if omitted_oversized:
        warnings.append(
            "%s custom sets exceeded GCTA's configured maximum set size and "
            "were omitted under the configured policy."
            % format(omitted_oversized, ",")
        )
    if missing_result_chromosomes:
        warnings.append(
            "No tested %s were present in the GCTA result for shared chromosome(s): "
            "%s."
            % (schema.unit_label, ", ".join(missing_result_chromosomes))
        )
    overlap_fraction = (
        resolved_variants / input_variants if input_variants else None
    )
    return {
        "method": module_config.method,
        "unit_label": schema.unit_label,
        "tested_units": tested_units,
        "source_units": source_units,
        "source_units_label": source_units_label,
        "input_variants": input_variants,
        "resolved_variants": resolved_variants,
        "overlap_fraction": overlap_fraction,
        "direct_input_unmodified": direct_input_unmodified,
        "reference_variants": int(metrics.get("reference_variants", 0)),
        "input_variants_absent_from_reference": absent_from_reference,
        "reference_variants_absent_from_input": int(
            metrics.get("reference_variants_absent_from_input", 0)
        ),
        "expected_chromosomes": expected_chromosomes,
        "result_chromosomes": result_chromosomes,
        "missing_result_chromosomes": missing_result_chromosomes,
        "familywise_alpha": familywise_alpha,
        "nominal_alpha": module_config.reporting.nominal_alpha,
        "fdr_alpha": module_config.reporting.fdr_alpha,
        "bonferroni_threshold": bonferroni_threshold,
        "nominal_significant": nominal_significant,
        "bonferroni_significant": bonferroni_significant,
        "fdr_bh_significant": fdr_bh_significant,
        "correction_columns": correction_columns.model_dump(mode="json"),
        "minimum_p_value": float(frame[p_column].min()),
        "p_value_column": p_column,
        "top_associations": top_associations,
        "warnings": warnings,
    }


def _interpretation_fields(
    summary: Mapping, significant_digits: int,
) -> list[tuple[str, str]]:
    threshold = _format_p_value(
        summary["bonferroni_threshold"],
        significant_digits,
    )
    method = summary["method"]
    if method == "mbat_combo":
        evidence = (
            "mBAT-combo combines signed mBAT and unsigned fastBAT evidence; "
            "component p-values are shown when requested and present."
        )
    else:
        evidence = (
            "fastBAT aggregates single-variant association evidence while "
            "accounting for LD estimated from the configured reference."
        )
    if method in {"fastbat_gene", "mbat_combo"}:
        causal = (
            "A gene association localizes statistical evidence to the configured "
            "gene window; it does not establish a causal gene or variant."
        )
    elif method == "fastbat_segment":
        causal = (
            "A segment association localizes statistical evidence to that interval; "
            "it does not establish a causal gene or variant."
        )
    else:
        causal = (
            "A set association indicates aggregate evidence among included variants; "
            "it does not establish that the set or any member is causal."
        )
    return [
        ("Association evidence", evidence),
        (
            "Nominal significance",
            "Raw p ≤ %s is uncorrected and should not be interpreted as "
            "family-wide significance."
            % _format_p_value(summary["nominal_alpha"], significant_digits),
        ),
        (
            "Benjamini-Hochberg FDR",
            "BH-adjusted p ≤ %s controls the false discovery rate across the "
            "reported family under the procedure's assumptions."
            % _format_p_value(summary["fdr_alpha"], significant_digits),
        ),
        (
            "Bonferroni significance",
            "Raw p ≤ %s controls the configured family-wise alpha across "
            "%s tested %s; results were not filtered by this threshold."
            % (threshold, format(summary["tested_units"], ","), summary["unit_label"]),
        ),
        ("Causal interpretation", causal),
        (
            "Reference dependence",
            "Results depend on the population-matched LD reference, available "
            "variants, reference MAF, and configured window, segment, or set definitions.",
        ),
    ]


def render_gcta_scientific_summary(
    summary: dict,
    dataset: str,
    module_config,
    artifacts: Mapping[str, Artifact],
    log_path: Path,
    label_width: int,
) -> str:
    """Render the validated scientific summary for terminal monitoring."""
    precision = module_config.reporting.p_value_significant_digits
    # Nested report fields are indented four cells farther than the report
    # metadata and output paths. Widen the outer label column by the same
    # amount so every separator and wrapped value begins in one terminal cell.
    outer_label_width = label_width + 4
    status = (
        "COMPLETED WITH SCIENTIFIC WARNINGS"
        if summary["warnings"] else "COMPLETED"
    )
    unit_label = summary["unit_label"]
    overlap = "%s of %s" % (
        format(summary["resolved_variants"], ","),
        format(summary["input_variants"], ","),
    )
    if summary["overlap_fraction"] is not None:
        overlap += " (%.1f%%)" % (100 * summary["overlap_fraction"])
    lines = [
        "",
        screen_line("analysis", _METHOD_TITLES[summary["method"]], indent=2),
        screen_field(
            "info", "Dataset", dataset, indent=6,
            label_width=outer_label_width,
        ),
        screen_field(
            "analysis", "Method", summary["method"], indent=6,
            label_width=outer_label_width,
        ),
        screen_field(
            "warning" if summary["warnings"] else "success",
            "Analysis status", status, indent=6,
            label_width=outer_label_width,
        ),
        "",
        screen_line("genetic", "Scientific findings", indent=6),
        screen_field(
            "count", "%s tested" % unit_label.capitalize(),
            format(summary["tested_units"], ","), indent=10, label_width=label_width,
        ),
    ]
    if summary["source_units"] is not None:
        lines.append(screen_field(
            "count", summary["source_units_label"],
            format(summary["source_units"], ","),
            indent=10, label_width=label_width,
        ))
    if summary["expected_chromosomes"]:
        lines.append(screen_field(
            "count", "Chromosomes represented",
            "%s of %s shared chromosomes" % (
                format(len(summary["result_chromosomes"]), ","),
                format(len(summary["expected_chromosomes"]), ","),
            ),
            indent=10, label_width=label_width,
        ))
    lines.append(screen_field(
        "count", "GWAS variants supplied", format(summary["input_variants"], ","),
        indent=10, label_width=label_width,
    ))
    if summary["direct_input_unmodified"]:
        lines.extend([
            screen_field(
                "count", "PLINK BIM variants supplied",
                format(summary["reference_variants"], ","),
                indent=10, label_width=label_width,
            ),
            screen_field(
                (
                    "warning"
                    if summary["input_variants_absent_from_reference"]
                    else "success"
                ),
                "Summary IDs absent from BIM",
                format(summary["input_variants_absent_from_reference"], ","),
                indent=10, label_width=label_width,
            ),
            screen_field(
                "count", "BIM IDs absent from summary",
                format(summary["reference_variants_absent_from_input"], ","),
                indent=10, label_width=label_width,
            ),
        ])
    lines.extend([
        screen_field(
            "count", "LD-reference-resolved variants", overlap,
            indent=10, label_width=label_width,
        ),
        screen_field(
            "info", "Nominal p-value threshold",
            _format_p_value(summary["nominal_alpha"], precision),
            indent=10, label_width=label_width,
        ),
        screen_field(
            "count", "Nominally significant associations",
            "%s of %s" % (
                format(summary["nominal_significant"], ","),
                format(summary["tested_units"], ","),
            ),
            indent=10, label_width=label_width,
        ),
        screen_field(
            "info", "Benjamini-Hochberg FDR threshold",
            _format_p_value(summary["fdr_alpha"], precision),
            indent=10, label_width=label_width,
        ),
        screen_field(
            "count", "FDR-significant associations",
            "%s of %s" % (
                format(summary["fdr_bh_significant"], ","),
                format(summary["tested_units"], ","),
            ),
            indent=10, label_width=label_width,
        ),
        screen_field(
            "info", "Family-wise alpha", summary["familywise_alpha"],
            indent=10, label_width=label_width,
        ),
        screen_field(
            "info", "Bonferroni threshold",
            _format_p_value(summary["bonferroni_threshold"], precision),
            indent=10, label_width=label_width,
        ),
        screen_field(
            "count", "Bonferroni-significant associations",
            "%s of %s" % (
                format(summary["bonferroni_significant"], ","),
                format(summary["tested_units"], ","),
            ),
            indent=10, label_width=label_width,
        ),
        screen_field(
            "analysis", "Smallest association p-value",
            _format_p_value(summary["minimum_p_value"], precision),
            indent=10, label_width=label_width,
        ),
        "",
        screen_line("decision", "Top associated %s" % unit_label, indent=6),
    ])
    for association in summary["top_associations"]:
        values = [
            "%s %s" % (
                summary["p_value_column"],
                _format_p_value(association["p_value"], precision),
            )
        ]
        values.extend([
            "%s %s" % (
                summary["correction_columns"]["fdr_bh_adjusted_p"],
                _format_p_value(association["fdr_bh_adjusted_p"], precision),
            ),
            "%s %s" % (
                summary["correction_columns"]["bonferroni_adjusted_p"],
                _format_p_value(
                    association["bonferroni_adjusted_p"], precision,
                ),
            ),
        ])
        values.extend(
            "%s %s" % (column, _format_p_value(value, precision))
            for column, value in association["component_p_values"].items()
        )
        lines.append(screen_field(
            "genetic", "%d. %s" % (association["rank"], association["identifier"]),
            " · ".join(values), indent=10, label_width=label_width,
        ))
    lines.extend(["", screen_line("decision", "How to interpret", indent=6)])
    lines.extend(
        screen_field(
            "info", label, value, indent=10, label_width=label_width,
        )
        for label, value in _interpretation_fields(summary, precision)
    )
    if summary["warnings"]:
        lines.extend(["", screen_line("warning", "Scientific warnings", indent=6)])
        lines.extend(
            screen_field(
                "warning", "Warning %d" % number, warning,
                indent=10, label_width=label_width,
            )
            for number, warning in enumerate(summary["warnings"], 1)
        )
    lines.extend([
        "",
        screen_field(
            "success", "Complete normalized results",
            artifacts["normalized_results"].path, indent=6,
            label_width=outer_label_width,
        ),
        screen_field(
            "info", "Original GCTA results", artifacts["raw_results"].path,
            indent=6, label_width=outer_label_width,
        ),
    ])
    if "frequency_qc" in artifacts:
        lines.append(screen_field(
            "info", "GCTA frequency-QC report", artifacts["frequency_qc"].path,
            indent=6, label_width=outer_label_width,
        ))
    lines.extend([
        screen_field(
            "info", "Full log", log_path, indent=6,
            label_width=outer_label_width,
        ),
        "",
    ])
    return "\n".join(lines)


def record_gcta_scientific_summary(
    logger: PipelineLogger,
    summary: Mapping,
) -> None:
    """Record summary findings, ranked results, and warnings in the canonical log."""
    logger.record("OBSERVED", "gcta_scientific_summary", **{
        key: value for key, value in summary.items()
        if key not in {"top_associations", "warnings"}
    })
    for association in summary["top_associations"]:
        logger.record("RESULT", "gcta_top_association", **association)
    for warning in summary["warnings"]:
        logger.warning(warning)


__all__ = [
    "build_gcta_scientific_summary",
    "record_gcta_scientific_summary",
    "render_gcta_scientific_summary",
]
