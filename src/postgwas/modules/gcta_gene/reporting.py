"""Scientific summaries for validated GCTA fastBAT and mBAT-combo results."""

from __future__ import annotations

from html import escape
import json
import os
from pathlib import Path
from typing import Any, Mapping, Sequence

import polars as pl

from postgwas.core.contracts import Artifact
from postgwas.core.dataframes import chromosome_expression
from postgwas.core.io.reports import write_html_report
from postgwas.core.pipeline_logging import PipelineLogger
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.core.values import format_count, format_percentage
from postgwas.core.vcf import vcf_query_field_label
from postgwas.modules.gcta_gene.errors import GctaGeneError


_METHOD_TITLES = {
    "fastbat_gene": "GCTA fastBAT gene-association summary",
    "fastbat_segment": "GCTA fastBAT segment-association summary",
    "fastbat_set": "GCTA fastBAT set-association summary",
    "mbat_combo": "GCTA mBAT-combo gene-association summary",
}


def _formatted_fraction(part, whole) -> str:
    return "%s/%s (%s)" % (
        format_count(part),
        format_count(whole),
        format_percentage(part, whole),
    )


def gcta_ma_input_outcome_fields(
    configuration,
    input_file: str | Path,
    input_metrics: Mapping,
    *,
    formatter_result: Mapping | None = None,
    variant_id_type: str | None = None,
    reference_metrics: Mapping | None = None,
) -> list[tuple]:
    """Describe one validated GCTA ``.ma`` input from configuration evidence."""
    formatting = configuration.modules.formatting
    canonical = formatting.canonical_columns
    schema = formatting.exports["gcta_gene"].outputs["summary_statistics"]
    required_columns = [
        schema.columns[source]
        for source in schema.validation.required_columns
    ]
    identifier_column = schema.columns[canonical.resolved_variant_id]
    allele_columns = (
        schema.columns[canonical.alternate_allele],
        schema.columns[canonical.reference_allele],
    )
    p_sources = [
        source
        for source, transformation in schema.transformations.items()
        if transformation == "negative_log10_to_raw_p"
    ]
    p_source = p_sources[0]
    p_column = schema.columns[p_source]
    path = Path(input_file)
    formatter = dict(formatter_result or {})
    fields: list[tuple] = []
    if formatter:
        vcf_fields = formatting.vcf_fields.root
        identifier_type = variant_id_type or formatter.get("variant_id_type")
        if identifier_type == "unique":
            identifier = formatting.variant_identifiers.unique_id_template.format(
                chromosome=canonical.chromosome,
                position=canonical.position,
                reference_allele=canonical.reference_allele,
                alternate_allele=canonical.alternate_allele,
            )
            identifier_provenance = "%s → %s" % (
                " + ".join((
                    canonical.chromosome,
                    canonical.position,
                    canonical.reference_allele,
                    canonical.alternate_allele,
                )),
                identifier,
            )
        elif identifier_type == "rsid":
            identifier_provenance = "%s → %s" % (
                vcf_query_field_label(vcf_fields[canonical.variant_id]),
                identifier_column,
            )
        else:
            identifier_provenance = (
                "configured BIM-compatible identifier → %s"
                % identifier_column
            )

        column_fields = []
        for source, destination in schema.columns.items():
            if source == canonical.resolved_variant_id:
                kind = "genetic"
                provenance = identifier_provenance
            else:
                vcf_source = vcf_query_field_label(vcf_fields[source])
                transformation = schema.transformations.get(source)
                if transformation == "negative_log10_to_raw_p":
                    kind = "analysis"
                    processing = "raw P = 10⁻ᴸᴾ"
                elif transformation is not None:
                    kind = "analysis"
                    processing = "configured transformation: %s" % transformation
                elif source == canonical.alternate_allele:
                    kind = "genetic"
                    processing = "effect allele"
                elif source == canonical.reference_allele:
                    kind = "genetic"
                    processing = "other allele"
                elif source == canonical.effect_allele_frequency:
                    kind = "genetic"
                    processing = "%s effect-allele frequency" % allele_columns[0]
                elif source == canonical.total_sample_size:
                    kind = "info"
                    processing = "per-variant total N"
                else:
                    kind = "info"
                    processing = "copied without numerical transformation"
                provenance = "%s → %s" % (vcf_source, processing)
            column_fields.append((kind, destination, provenance))

        fields.extend((
            (
                "analysis",
                "Prepare GCTA summary statistics from the harmonised GWAS-VCF",
            ),
            ("info", "GCTA .ma file", path.name),
            ("info", "Input format", ".ma with a header row"),
            ("analysis", "Extract and transform GWAS-VCF fields"),
            ("count", "Variants extracted", formatter.get("rows_in")),
            *column_fields,
            ("analysis", "Write the GCTA .ma input"),
            ("success", "Variants written", formatter.get("rows_out")),
            (
                "loss" if formatter.get("rows_excluded", 0) else "success",
                "Variants excluded during formatting",
                formatter.get("rows_excluded", 0),
            ),
            (
                "genetic",
                "Selected variant-ID format",
                identifier_type or "configured BIM-compatible identifiers",
            ),
            (
                "warning" if formatter.get("p_values_bounded", 0) else "success",
                "Converted P values bounded at configured minimum",
                formatter.get("p_values_bounded", 0),
            ),
        ))
    else:
        fields.extend((
            ("analysis", "Supplied GCTA summary-statistics input"),
            ("info", "GCTA .ma file", path.name),
            ("info", "Input format", ".ma with a header row"),
            ("count", "Variants read", input_metrics["variants"]),
        ))
    fields.extend((
        ("analysis", "Validate the final GCTA .ma input"),
        ("info", "Required columns", ", ".join(required_columns)),
        ("genetic", "Variant-ID column", identifier_column),
        ("info", "Raw P-value column", p_column),
        ("genetic", "Effect/reference allele columns", "/".join(allele_columns)),
        ("success", "Required columns and non-missing values", "passed"),
        ("success", "Scientific numeric ranges", "passed"),
        ("success", "Unique non-empty variant identifiers", "passed"),
        ("success", "Non-empty distinct allele pairs", "passed"),
        ("success", "File-structure validation", "passed"),
    ))
    if reference_metrics is not None:
        reference = dict(reference_metrics)
        total = int(input_metrics["variants"])
        shared = int(reference["overlapping_variants"])
        absent = int(reference["input_variants_absent_from_reference"])
        fields.extend((
            ("analysis", "Compare the final .ma input with PLINK BIM"),
            (
                "success",
                "Exact IDs shared with BIM",
                _formatted_fraction(shared, total),
            ),
            (
                "warning" if absent else "success",
                "Summary IDs absent from BIM",
                "%s; will not be used by GCTA" % format_count(absent)
                if absent else "0",
            ),
            (
                "count",
                "BIM IDs absent from summary",
                reference["reference_variants_absent_from_input"],
            ),
            (
                "success",
                "Compatible shared allele pairs",
                reference["compatible_allele_pairs"],
            ),
            (
                "success",
                "Incompatible shared allele pairs",
                reference["incompatible_allele_pairs"],
            ),
            ("success", "Second reconciled .ma file", "not created"),
        ))
    return fields


def gcta_ld_reference_outcome_fields(
    module,
    reference_prefix: str | Path,
    reference_paths,
    *,
    reference_metrics: Mapping | None = None,
    variant_id_type: str | None = None,
    summary_variants: int | None = None,
    enforce_minimum_overlap: bool = False,
) -> list[tuple]:
    """Describe validated PLINK files and optional GWAS/BIM compatibility."""
    metrics = dict(reference_metrics or {})
    identifier_label = {
        "rsid": "rsID",
        "unique": "chromosome-position-allele ID",
    }.get(variant_id_type, "exact PLINK BIM column 2 values")
    fields: list[tuple] = [
        ("analysis", "PLINK LD reference"),
        ("info", "Reference prefix", Path(reference_prefix).name),
        (
            "success",
            "Required companion files",
            "%s present and non-empty"
            % ", ".join(
                Path(path).suffix.removeprefix(".").upper()
                for path in reference_paths
            ),
        ),
    ]
    if metrics.get("reference_variants") is not None:
        fields.append((
            "count", "PLINK BIM unique IDs", metrics["reference_variants"],
        ))
    fields.extend((
        ("genetic", "Variant identifiers", identifier_label),
        ("info", "BIM columns", ", ".join(module.reference.bim_columns)),
    ))
    if metrics.get("reference_chromosomes"):
        fields.append((
            "genetic",
            "Chromosomes represented",
            ", ".join(metrics["reference_chromosomes"]),
        ))
    analysis_scope = metrics.get("analysis_scope")
    if analysis_scope is not None:
        mhc = analysis_scope.get("mhc_region")
        fields.extend((
            ("analysis", "Configured analysis scope"),
            ("info", "MHC policy", analysis_scope["mhc_policy"]),
            (
                "info",
                "MHC region",
                (
                    "%s:%s-%s"
                    % (mhc["chromosome"], mhc["start"], mhc["end"])
                    if mhc is not None else "not configured"
                ),
            ),
            (
                "info",
                "Excluded chromosomes",
                ", ".join(analysis_scope["exclude_chromosomes"]) or "none",
            ),
            (
                "loss" if metrics["excluded_reference_variants"] else "success",
                "LD-reference variants excluded by scope",
                metrics["excluded_reference_variants"],
            ),
        ))
    fields.extend((
        ("genetic", "Declared genome build", module.genome_build),
        ("genetic", "Declared population", module.reference.population),
        ("success", "BIM structural validation", "passed"),
    ))
    if (
        summary_variants is not None
        and metrics.get("overlapping_variants") is not None
    ):
        shared = int(metrics["overlapping_variants"])
        absent = int(metrics["input_variants_absent_from_reference"])
        fields.extend((
            ("analysis", "GWAS/BIM identifier and allele compatibility"),
            ("count", "Summary-statistic unique IDs", summary_variants),
            (
                "success", "Exact IDs shared",
                _formatted_fraction(shared, summary_variants),
            ),
            (
                "warning" if absent else "success",
                "Summary IDs absent from BIM",
                "%s; will not be used by GCTA" % format_count(absent)
                if absent else "0",
            ),
            (
                "count",
                "BIM IDs absent from summary",
                metrics["reference_variants_absent_from_input"],
            ),
            (
                "success",
                "Compatible shared allele pairs",
                metrics["compatible_allele_pairs"],
            ),
            (
                "success",
                "Incompatible shared allele pairs",
                metrics["incompatible_allele_pairs"],
            ),
            (
                "info",
                "Minimum ID-overlap enforcement",
                (
                    format_percentage(
                        module.variant_harmonisation.minimum_overlap_fraction,
                        1,
                    )
                    if enforce_minimum_overlap
                    else "not applied to supplied direct input; overlap is reported"
                ),
            ),
            ("success", "Identifier and allele compatibility", "passed"),
            ("success", "Summary-statistics input rewritten", "no"),
            (
                "loss" if metrics["excluded_input_variants"] else "success",
                "Shared variants excluded by analysis scope",
                metrics["excluded_input_variants"],
            ),
            (
                "success",
                "Shared variants eligible for GCTA",
                metrics["analyzable_variants"],
            ),
        ))
    fields.append((
        "warning",
        "Build and population provenance",
        "taken from configuration; neither can be inferred from PLINK file contents",
    ))
    return fields


def gcta_gene_reference_outcome_fields(
    module,
    gene_file: str | Path,
    metrics: Mapping,
) -> list[tuple]:
    """Describe one validated GCTA gene-coordinate reference."""
    fields: list[tuple] = [
        ("analysis", "GCTA gene-coordinate reference"),
        ("info", "Reference file", Path(gene_file).name),
        ("info", "Coordinate columns", ", ".join(module.gene_annotation.columns)),
        ("count", "Unique gene identifiers", metrics["genes"]),
        (
            "genetic",
            "Chromosomes represented",
            ", ".join(metrics["gene_chromosomes"]),
        ),
        ("genetic", "Declared genome build", module.genome_build),
        ("info", "Configured gene window", "±%s kb" % module.gene_window_kb),
    ]
    if metrics.get("shared_chromosomes"):
        fields.append((
            "success",
            "Chromosomes shared with BIM",
            ", ".join(metrics["shared_chromosomes"]),
        ))
    fields.extend((
        ("success", "Unique non-empty gene identifiers", "passed"),
        ("success", "Positive ordered coordinates", "passed"),
        ("success", "File-structure validation", "passed"),
    ))
    if metrics.get("input_genes") is not None:
        fields.extend((
            ("analysis", "Coordinate-defined gene scope"),
            ("count", "Genes before analysis-scope exclusions", metrics["input_genes"]),
            (
                "loss" if metrics["excluded_genes"] else "success",
                "Genes excluded by chromosome/MHC policy",
                metrics["excluded_genes"],
            ),
            ("success", "Genes eligible for GCTA", metrics["retained_genes"]),
        ))
    return fields


def gcta_gmt_outcome_fields(
    module,
    gmt_file: str | Path,
    pathway_count: int,
    unique_gene_count: int,
) -> list[tuple]:
    """Describe one validated original pathway GMT resource."""
    conversion = module.set_annotation.conversion
    return [
        ("analysis", "Original pathway GMT"),
        ("info", "Pathway file", Path(gmt_file).name),
        ("info", "Detected format", "GMT"),
        ("count", "Pathway records", pathway_count),
        ("genetic", "Unique GMT gene identifiers", unique_gene_count),
        ("info", "Duplicate-gene policy", conversion.duplicate_gene_policy),
        ("success", "Unique non-empty pathway identifiers", "passed"),
        ("success", "Pathway gene-membership validation", "passed"),
        ("success", "File-structure validation", "passed"),
    ]


def _gcta_gene_overlap_fields(
    module,
    requested: int,
    matched: int,
    missing: int,
    coordinate_reference_genes: int,
) -> list[tuple]:
    conversion = module.set_annotation.conversion
    coordinate_genes_absent = coordinate_reference_genes - matched
    return [
        ("analysis", "Pathway/gene-coordinate compatibility"),
        ("count", "Unique GMT gene identifiers", requested),
        ("success", "Genes represented in coordinate file", matched),
        (
            "warning" if missing else "success",
            "GMT genes absent from coordinate file",
            missing,
        ),
        (
            "genetic",
            "GMT gene mappability",
            format_percentage(matched, requested),
        ),
        (
            "decision",
            "GMT gene mappability criterion",
            "reported only; not used as the stopping criterion",
        ),
        (
            "count",
            "Coordinate-reference genes",
            coordinate_reference_genes,
        ),
        (
            "success",
            "Coordinate-reference genes represented in GMT",
            matched,
        ),
        (
            "warning" if coordinate_genes_absent else "success",
            "Coordinate-reference genes absent from GMT",
            coordinate_genes_absent,
        ),
        (
            "genetic",
            "Coordinate-reference coverage by GMT",
            format_percentage(matched, coordinate_reference_genes),
        ),
        (
            "info",
            "Required minimum coordinate-reference coverage",
            format_percentage(
                conversion.minimum_gene_id_overlap_fraction, 1,
            ),
        ),
        ("info", "Unmapped-gene policy", conversion.unmapped_gene_policy),
        ("decision", "Pathway genes retained", "matched genes only"),
        ("success", "Gene-ID compatibility validation", "passed"),
    ]


def gcta_gene_overlap_outcome_fields(module, preflight) -> list[tuple]:
    """Describe validated compatibility between GMT genes and coordinates."""
    return _gcta_gene_overlap_fields(
        module,
        len(preflight.requested_genes),
        len(preflight.gene_coordinates),
        len(preflight.missing_genes),
        len(preflight.all_gene_coordinates),
    )


def gcta_set_source_outcome_fields(
    set_file: str | Path,
    metrics: Mapping,
) -> list[tuple]:
    """Describe one structurally validated native fastBAT set-list input."""
    return [
        ("analysis", "Prepared fastBAT set-list input"),
        ("info", "Set-list file", Path(set_file).name),
        ("info", "Input format", "GCTA set blocks terminated by END"),
        ("count", "Input sets", metrics["input_sets"]),
        ("count", "Variant memberships", metrics["requested_set_variants"]),
        (
            "count",
            "Unique variant identifiers",
            metrics["unique_requested_set_variants"],
        ),
        ("success", "Unique set identifiers", "passed"),
        ("success", "Within-set variant uniqueness", "passed"),
        ("success", "END terminators and non-empty sets", "passed"),
        ("success", "Set-block validation", "passed"),
    ]


def gcta_gmt_preparation_outcome_fields(
    module,
    gmt_file: str | Path,
    gene_file: str | Path,
    set_file: str | Path,
    manifest: Mapping,
) -> list[tuple]:
    """Describe validated GMT, gene reference, compatibility, and final sets."""
    validation = manifest["validation"]
    requested = validation["input_unique_genes"]
    matched = validation["genes_with_coordinates"]
    missing = validation["unmapped_unique_genes"]
    fields = gcta_gmt_outcome_fields(
        module,
        gmt_file,
        validation["input_pathways"],
        requested,
    )
    fields.extend(gcta_gene_reference_outcome_fields(
        module,
        gene_file,
        {
            "genes": validation["gene_coordinate_reference_genes"],
            "gene_chromosomes": validation.get(
                "gene_coordinate_chromosomes",
                validation["shared_chromosomes"],
            ),
            "shared_chromosomes": validation["shared_chromosomes"],
        },
    ))
    fields.extend(_gcta_gene_overlap_fields(
        module,
        requested,
        matched,
        missing,
        validation["gene_coordinate_reference_genes"],
    ))
    if validation.get("genes_excluded_by_analysis_scope") is not None:
        fields.extend((
            ("analysis", "Pathway genes after analysis-scope exclusions"),
            (
                "loss"
                if validation["genes_excluded_by_analysis_scope"]
                else "success",
                "Coordinate-matched genes excluded by scope",
                validation["genes_excluded_by_analysis_scope"],
            ),
            (
                "success",
                "Coordinate-matched genes eligible for analysis",
                validation["genes_eligible_for_analysis"],
            ),
        ))
    fields.extend([
        ("analysis", "Final fastBAT set resource"),
        ("info", "Set-list file", Path(set_file).name),
        ("success", "Final sets written", validation["pathways_written"]),
        (
            "warning" if validation["pathways_omitted_empty"] else "success",
            "Empty pathways omitted",
            validation["pathways_omitted_empty"],
        ),
        (
            "warning"
            if validation["pathways_omitted_oversized"]
            else "success",
            "Oversized pathways omitted",
            validation["pathways_omitted_oversized"],
        ),
        (
            "count",
            "Retained variant memberships",
            validation["total_pathway_variant_memberships"],
        ),
        (
            "success",
            "Pathways with complete gene-ID mapping",
            validation["pathways_with_complete_gene_id_mapping"],
        ),
        (
            "warning"
            if validation["pathways_with_partial_gene_id_mapping"]
            else "success",
            "Pathways with partial gene-ID mapping",
            validation["pathways_with_partial_gene_id_mapping"],
        ),
        (
            "warning"
            if validation["pathways_without_gene_id_mapping"]
            else "success",
            "Pathways without gene-ID mapping",
            validation["pathways_without_gene_id_mapping"],
        ),
        ("success", "Published resource validation", "passed"),
    ])
    return fields


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
        metrics.get(
            "analyzable_variants",
            metrics.get("overlapping_variants", metrics.get("analysis_input_variants", 0)),
        )
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
            metrics.get(
                "retained_gene_chromosomes",
                metrics.get("shared_chromosomes"),
            )
            if module_config.method in {"fastbat_gene", "mbat_combo"}
            else metrics.get("reference_chromosomes")
        )
        excluded_chromosomes = set(
            metrics.get("analysis_scope", {}).get("exclude_chromosomes", ())
        )
        expected_chromosomes = [
            str(value) for value in (expected_source or [])
            if str(value) not in excluded_chromosomes
        ]
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
    formatter_input_unmodified = bool(
        metrics.get("formatter_input_unmodified", False)
    )
    absent_from_reference = int(
        metrics.get("input_variants_absent_from_reference", 0)
    )
    if (direct_input_unmodified or formatter_input_unmodified) and absent_from_reference:
        variant_label = (
            "variant is" if absent_from_reference == 1 else "variants are"
        )
        source_label = (
            "direct-input" if direct_input_unmodified else "formatter-created"
        )
        warnings.append(
            "%s %s GWAS %s absent from PLINK BIM column 2 and will not be "
            "used by GCTA. PostGWAS passed the validated summary-statistics "
            "input without rewriting it."
            % (
                format(absent_from_reference, ","), source_label,
                variant_label,
            )
        )
    excluded_scope_variants = int(metrics.get("excluded_input_variants", 0))
    if excluded_scope_variants:
        warnings.append(
            "%s GWAS/BIM-shared variants were excluded under the configured "
            "chromosome and MHC analysis scope."
            % format(excluded_scope_variants, ",")
        )
    excluded_scope_genes = int(
        metrics.get(
            "genes_excluded_by_analysis_scope",
            metrics.get("excluded_genes", 0),
        )
    )
    if excluded_scope_genes:
        warnings.append(
            "%s coordinate-defined genes were excluded under the configured "
            "chromosome and MHC analysis scope."
            % format(excluded_scope_genes, ",")
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
    unmapped_unique_genes = int(metrics.get("unmapped_unique_genes", 0))
    if unmapped_unique_genes:
        warnings.append(
            "%s unique GMT genes were absent from the gene-coordinate "
            "reference and could not contribute variants; review the "
            "pathway-level gene-ID mappability audit."
            % format(unmapped_unique_genes, ",")
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
    input_unique_genes = metrics.get("input_unique_genes")
    gene_coordinate_reference_genes = metrics.get(
        "gene_coordinate_reference_genes"
    )
    genes_with_coordinates = metrics.get("genes_with_coordinates")
    return {
        "method": module_config.method,
        "unit_label": schema.unit_label,
        "tested_units": tested_units,
        "source_units": source_units,
        "source_units_label": source_units_label,
        "input_variants": input_variants,
        "resolved_variants": resolved_variants,
        "overlap_fraction": overlap_fraction,
        "input_unique_genes": input_unique_genes,
        "gene_coordinate_reference_genes": gene_coordinate_reference_genes,
        "genes_with_coordinates": genes_with_coordinates,
        "gmt_gene_mappability_fraction": metrics.get(
            "gmt_gene_mappability_fraction"
        ),
        "coordinate_reference_coverage_fraction": metrics.get(
            "coordinate_reference_coverage_fraction"
        ),
        "coordinate_reference_genes_absent_from_gmt": metrics.get(
            "coordinate_reference_genes_absent_from_gmt"
        ),
        "unmapped_unique_genes": unmapped_unique_genes,
        "pathways_with_complete_gene_id_mapping": metrics.get(
            "pathways_with_complete_gene_id_mapping"
        ),
        "pathways_with_partial_gene_id_mapping": metrics.get(
            "pathways_with_partial_gene_id_mapping"
        ),
        "pathways_without_gene_id_mapping": metrics.get(
            "pathways_without_gene_id_mapping"
        ),
        "direct_input_unmodified": direct_input_unmodified,
        "formatter_input_unmodified": formatter_input_unmodified,
        "reference_variants": int(metrics.get("reference_variants", 0)),
        "input_variants_absent_from_reference": absent_from_reference,
        "reference_variants_absent_from_input": int(
            metrics.get("reference_variants_absent_from_input", 0)
        ),
        "analysis_scope": metrics.get("analysis_scope", {}),
        "excluded_scope_variants": excluded_scope_variants,
        "excluded_scope_genes": excluded_scope_genes,
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
    if summary["input_unique_genes"] is not None:
        lines.extend([
            screen_field(
                "genetic", "GMT gene mappability",
                _formatted_fraction(
                    summary["genes_with_coordinates"],
                    summary["input_unique_genes"],
                ),
                indent=10, label_width=label_width,
            ),
            screen_field(
                "decision", "GMT gene mappability criterion",
                "reported only; not used as the stopping criterion",
                indent=10, label_width=label_width,
            ),
            screen_field(
                "genetic", "Coordinate-reference coverage by GMT",
                _formatted_fraction(
                    summary["genes_with_coordinates"],
                    summary["gene_coordinate_reference_genes"],
                ),
                indent=10, label_width=label_width,
            ),
            screen_field(
                "info", "Required minimum coordinate-reference coverage",
                format_percentage(
                    module_config.set_annotation.conversion
                    .minimum_gene_id_overlap_fraction,
                    1,
                ),
                indent=10, label_width=label_width,
            ),
            screen_field(
                "count", "Pathways with complete gene-ID mapping",
                format(
                    summary["pathways_with_complete_gene_id_mapping"], ","
                ),
                indent=10, label_width=label_width,
            ),
            screen_field(
                "warning"
                if summary["pathways_with_partial_gene_id_mapping"]
                else "success",
                "Pathways with partial gene-ID mapping",
                format(
                    summary["pathways_with_partial_gene_id_mapping"], ","
                ),
                indent=10, label_width=label_width,
            ),
            screen_field(
                "warning"
                if summary["pathways_without_gene_id_mapping"]
                else "success",
                "Pathways without gene-ID mapping",
                format(summary["pathways_without_gene_id_mapping"], ","),
                indent=10, label_width=label_width,
            ),
        ])
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
    analysis_scope = summary["analysis_scope"]
    if analysis_scope:
        mhc = analysis_scope.get("mhc_region")
        mhc_label = (
            "%s:%s-%s" % (mhc["chromosome"], mhc["start"], mhc["end"])
            if mhc is not None else "not configured"
        )
        lines.extend([
            screen_field(
                "info", "MHC policy", analysis_scope["mhc_policy"],
                indent=10, label_width=label_width,
            ),
            screen_field(
                "info", "MHC region", mhc_label,
                indent=10, label_width=label_width,
            ),
            screen_field(
                "info", "Excluded chromosomes",
                ", ".join(analysis_scope["exclude_chromosomes"]) or "none",
                indent=10, label_width=label_width,
            ),
            screen_field(
                "count", "Shared variants excluded by scope",
                format(summary["excluded_scope_variants"], ","),
                indent=10, label_width=label_width,
            ),
        ])
        if summary["excluded_scope_genes"]:
            lines.append(screen_field(
                "count", "Coordinate-defined genes excluded by scope",
                format(summary["excluded_scope_genes"], ","),
                indent=10, label_width=label_width,
            ))
    if (
        summary["direct_input_unmodified"]
        or summary["formatter_input_unmodified"]
    ):
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
            "success", "HTML scientific report",
            artifacts["html_report"].path, indent=6,
            label_width=outer_label_width,
        ),
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


def _html_value(value: Any, significant_digits: int) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return format(value, ",")
    if isinstance(value, float):
        return (
            "NA" if value != value
            else _format_p_value(value, significant_digits)
        )
    return str(value)


def _html_metric(
    label: str,
    value: Any,
    significant_digits: int,
    tone: str = "",
) -> str:
    return (
        '<div class="metric %s"><span>%s</span><strong>%s</strong></div>'
        % (
            escape(tone, quote=True),
            escape(label),
            escape(_html_value(value, significant_digits)),
        )
    )


def _html_table(
    rows: Sequence[Mapping[str, Any]],
    columns: Sequence[str],
    significant_digits: int,
) -> str:
    if not rows:
        return '<p class="muted">No records are available.</p>'
    header = "".join("<th>%s</th>" % escape(column) for column in columns)
    body = "".join(
        "<tr>%s</tr>" % "".join(
            "<td>%s</td>"
            % escape(_html_value(row.get(column), significant_digits))
            for column in columns
        )
        for row in rows
    )
    return (
        '<div class="table-wrap"><table><thead><tr>%s</tr></thead>'
        '<tbody>%s</tbody></table></div>' % (header, body)
    )


def _html_file_link(path: str | Path, report_path: Path) -> str:
    source = Path(path)
    href = os.path.relpath(source, report_path.parent)
    return '<a href="%s">%s</a>' % (
        escape(href, quote=True), escape(source.name),
    )


def _gcta_html_result_rows(
    normalized_result: str | Path,
    module_config,
) -> tuple[list[list[str]], list[str]]:
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
            "Cannot read normalized GCTA results for the HTML report %s: %s"
            % (path, exc)
        ) from exc
    schema = module_config.results.schemas[module_config.method]
    configured_columns = module_config.html_report.columns[module_config.method]
    correction_columns = set(
        module_config.reporting.correction_columns.model_dump().values()
    )
    mandatory_columns = set(schema.required_columns) | correction_columns
    missing = [
        column for column in configured_columns
        if column not in frame.columns and column in mandatory_columns
    ]
    if missing:
        raise GctaGeneError(
            "Normalized GCTA result is missing configured HTML-report columns: %s"
            % ", ".join(missing)
        )
    columns = [column for column in configured_columns if column in frame.columns]
    if frame.is_empty():
        raise GctaGeneError(
            "Normalized GCTA result contains no associations for the HTML report: %s"
            % path
        )
    precision = module_config.reporting.p_value_significant_digits
    ranked = frame.select(columns).sort(
        schema.p_value_column, maintain_order=True,
    )
    return [
        [_html_value(value, precision) for value in row]
        for row in ranked.iter_rows()
    ], columns


def _paginated_gcta_table(
    rows: Sequence[Sequence[str]],
    columns: Sequence[str],
    *,
    page_size: int,
) -> str:
    header = "".join(
        '<th><button type="button" data-column="%d">%s</button></th>'
        % (index, escape(column))
        for index, column in enumerate(columns)
    )
    payload = json.dumps(
        rows, ensure_ascii=False, separators=(",", ":"),
    ).replace("</", "<\\/")
    return """
<div class="paginated-results" data-page-size="%d" data-source="gcta-results-data">
  <div class="table-controls">
    <label>Search all rows <input type="search" class="result-search" autocomplete="off"></label>
    <div class="page-controls">
      <button type="button" class="previous-page">Previous</button>
      <span class="page-status" aria-live="polite"></span>
      <button type="button" class="next-page">Next</button>
    </div>
  </div>
  <div class="table-wrap"><table><thead><tr>%s</tr></thead><tbody></tbody></table></div>
  <noscript><p class="notice">Enable JavaScript to browse this table. The complete normalized TSV is linked below.</p></noscript>
</div>
<script type="application/json" id="gcta-results-data">%s</script>""" % (
        page_size, header, payload,
    )


def write_gcta_html_report(
    normalized_result: str | Path,
    destination: str | Path,
    *,
    dataset_id: str,
    module_config,
    summary: Mapping,
    gcta_version: str,
    output_files: Mapping[str, str | Path],
) -> Path:
    """Write a complete standalone report from validated GCTA results."""
    destination = Path(destination)
    rows, columns = _gcta_html_result_rows(normalized_result, module_config)
    if len(rows) != int(summary["tested_units"]):
        raise GctaGeneError(
            "HTML-report row count does not match the validated GCTA result: "
            "%d != %d" % (len(rows), summary["tested_units"])
        )
    precision = module_config.reporting.p_value_significant_digits
    status = (
        "COMPLETED WITH SCIENTIFIC WARNINGS"
        if summary["warnings"] else "COMPLETED"
    )
    title = _METHOD_TITLES[summary["method"]].replace("summary", "report")
    overlap = "%s/%s (%s)" % (
        format_count(summary["resolved_variants"]),
        format_count(summary["input_variants"]),
        format_percentage(
            summary["resolved_variants"], summary["input_variants"],
        ),
    )
    metrics = "".join((
        _html_metric(
            "%s tested" % summary["unit_label"].capitalize(),
            summary["tested_units"], precision,
        ),
        _html_metric(
            "Nominally significant",
            summary["nominal_significant"], precision,
        ),
        _html_metric(
            "BH-FDR significant",
            summary["fdr_bh_significant"], precision, "good",
        ),
        _html_metric(
            "Bonferroni significant",
            summary["bonferroni_significant"], precision, "good",
        ),
        _html_metric("Smallest p-value", summary["minimum_p_value"], precision),
        _html_metric("LD-reference-resolved variants", overlap, precision),
    ))
    scope = summary.get("analysis_scope") or {}
    mhc = scope.get("mhc_region")
    coverage_rows = [
        {"Assessment": "GWAS variants supplied", "Result": summary["input_variants"]},
        {"Assessment": "PLINK BIM variants", "Result": summary["reference_variants"]},
        {"Assessment": "LD-reference-resolved variants", "Result": overlap},
        {
            "Assessment": "Summary IDs absent from BIM",
            "Result": summary["input_variants_absent_from_reference"],
        },
        {
            "Assessment": "BIM IDs absent from summary",
            "Result": summary["reference_variants_absent_from_input"],
        },
        {
            "Assessment": "Shared variants excluded by analysis scope",
            "Result": summary["excluded_scope_variants"],
        },
        {
            "Assessment": "Coordinate-defined genes excluded by analysis scope",
            "Result": summary["excluded_scope_genes"],
        },
    ]
    if summary["source_units"] is not None:
        coverage_rows.insert(0, {
            "Assessment": summary["source_units_label"],
            "Result": summary["source_units"],
        })
    if summary["input_unique_genes"] is not None:
        coverage_rows[0:0] = [
            {
                "Assessment": "GMT gene mappability",
                "Result": _formatted_fraction(
                    summary["genes_with_coordinates"],
                    summary["input_unique_genes"],
                ),
            },
            {
                "Assessment": "GMT gene mappability criterion",
                "Result": "reported only; not used as the stopping criterion",
            },
            {
                "Assessment": "Coordinate-reference coverage by GMT",
                "Result": _formatted_fraction(
                    summary["genes_with_coordinates"],
                    summary["gene_coordinate_reference_genes"],
                ),
            },
            {
                "Assessment": "Required minimum coordinate-reference coverage",
                "Result": format_percentage(
                    module_config.set_annotation.conversion
                    .minimum_gene_id_overlap_fraction,
                    1,
                ),
            },
            {
                "Assessment": "Coordinate-reference genes absent from GMT",
                "Result": summary[
                    "coordinate_reference_genes_absent_from_gmt"
                ],
            },
            {
                "Assessment": "Pathways with complete gene-ID mapping",
                "Result": summary[
                    "pathways_with_complete_gene_id_mapping"
                ],
            },
            {
                "Assessment": "Pathways with partial gene-ID mapping",
                "Result": summary[
                    "pathways_with_partial_gene_id_mapping"
                ],
            },
            {
                "Assessment": "Pathways without gene-ID mapping",
                "Result": summary["pathways_without_gene_id_mapping"],
            },
        ]
    significance_rows = [
        {
            "Procedure": "Nominal",
            "Threshold": summary["nominal_alpha"],
            "Significant associations": summary["nominal_significant"],
        },
        {
            "Procedure": "Benjamini-Hochberg FDR",
            "Threshold": summary["fdr_alpha"],
            "Significant associations": summary["fdr_bh_significant"],
        },
        {
            "Procedure": "Bonferroni family-wise error",
            "Threshold": summary["bonferroni_threshold"],
            "Significant associations": summary["bonferroni_significant"],
        },
    ]
    mhc_label = (
        "%s:%s-%s" % (mhc["chromosome"], mhc["start"], mhc["end"])
        if mhc is not None else "not configured"
    )
    parameters = [
        {"Parameter": "Method", "Value": module_config.method},
        {"Parameter": "Genome build", "Value": module_config.genome_build},
        {
            "Parameter": "LD-reference population",
            "Value": module_config.reference.population,
        },
        {"Parameter": "GCTA version", "Value": gcta_version},
        {"Parameter": "Reference MAF minimum", "Value": module_config.reference_maf_min},
        {
            "Parameter": "Maximum frequency difference",
            "Value": module_config.frequency_difference_max,
        },
        {"Parameter": "fastBAT LD cutoff", "Value": module_config.fastbat_ld_cutoff},
        {"Parameter": "MHC policy", "Value": scope.get("mhc_policy")},
        {"Parameter": "MHC region", "Value": mhc_label},
        {
            "Parameter": "Excluded chromosomes",
            "Value": ", ".join(scope.get("exclude_chromosomes", ())) or "none",
        },
    ]
    if module_config.method in {"fastbat_gene", "mbat_combo"}:
        parameters.append({
            "Parameter": "Gene window (kb)",
            "Value": module_config.gene_window_kb,
        })
    if module_config.method == "fastbat_segment":
        parameters.append({
            "Parameter": "Segment size (kb)",
            "Value": module_config.segment_size_kb,
        })
    if module_config.method == "fastbat_set":
        parameters.extend((
            {
                "Parameter": "Maximum variants per set",
                "Value": module_config.set_annotation.maximum_set_variants,
            },
            {
                "Parameter": "Oversized-set policy",
                "Value": module_config.set_annotation.oversized_set_policy,
            },
        ))
    if module_config.method == "mbat_combo":
        parameters.append({
            "Parameter": "mBAT SVD gamma", "Value": module_config.mbat_svd_gamma,
        })
    warning_section = (
        "<section><h2>Scientific warnings</h2><ol>%s</ol></section>"
        % "".join("<li>%s</li>" % escape(item) for item in summary["warnings"])
        if summary["warnings"] else ""
    )
    interpretation_rows = [
        {"Topic": label, "Interpretation": value}
        for label, value in _interpretation_fields(summary, precision)
    ]
    downloads = "".join(
        "<li><strong>%s:</strong> %s</li>"
        % (escape(label), _html_file_link(path, destination))
        for label, path in output_files.items()
    )
    document = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>%s · %s</title>
<style>
:root{--ink:#172033;--muted:#607084;--line:#dbe2ea;--blue:#075985;--pale:#eef6fb;--green:#08783e;--amber:#9a6700}
*{box-sizing:border-box}body{font-family:Inter,ui-sans-serif,system-ui,sans-serif;margin:0;color:var(--ink);background:#f3f6fa;line-height:1.45}
main{max-width:1500px;margin:auto;padding:2rem}header,section{background:white;padding:1.5rem 1.7rem;margin-bottom:1.2rem;border:1px solid var(--line);border-radius:14px;box-shadow:0 2px 10px #1f293708}
h1{margin:.1rem 0}.subtitle,.muted{color:var(--muted)}.lead{max-width:95ch}.status{font-weight:700;color:%s}.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:.8rem;margin:1rem 0}
.metric{border-left:4px solid var(--blue);background:var(--pale);padding:.8rem 1rem;border-radius:7px}.metric span{display:block;color:var(--muted);font-size:.82rem}.metric strong{font-size:1.2rem}.metric.good{border-color:var(--green)}
.table-wrap{overflow:auto;border:1px solid var(--line);border-radius:8px;margin:.8rem 0 1.2rem}table{border-collapse:collapse;width:100%%;font-size:.88rem;white-space:nowrap}th,td{border-bottom:1px solid var(--line);padding:.55rem .7rem;text-align:left}th{background:var(--pale);position:sticky;top:0}
.table-controls{display:flex;align-items:center;justify-content:space-between;gap:1rem;flex-wrap:wrap;margin:.8rem 0}.table-controls label{font-weight:700}.table-controls input{margin-left:.45rem;min-width:18rem;padding:.48rem .6rem;border:1px solid var(--line);border-radius:6px}.page-controls{display:flex;align-items:center;gap:.65rem}.page-controls button,th button{font:inherit}.page-controls button{padding:.42rem .7rem;border:1px solid var(--line);border-radius:6px;background:white;color:var(--blue);cursor:pointer}.page-controls button:disabled{color:var(--muted);cursor:not-allowed}th button{border:0;background:transparent;color:inherit;font-weight:700;padding:0;cursor:pointer}.page-status{min-width:15rem;text-align:center;color:var(--muted)}
.notice{border-left:4px solid var(--amber);background:#fff8e6;padding:.8rem 1rem}.downloads{background:#f8fafc;padding:.8rem 1rem;border-radius:8px}a{color:var(--blue)}
@media(max-width:700px){main{padding:.7rem}header,section{padding:1rem}table{font-size:.8rem}.table-controls input{min-width:12rem}.page-status{min-width:auto}}
</style></head><body><main>
<header><p class="subtitle">PostGWAS scientific results report</p><h1>%s</h1><p><strong>Dataset:</strong> %s</p><p><strong>Status:</strong> <span class="status">%s</span></p>
<p class="lead">This report summarizes the validated GCTA association results and preserves every tested result row. Association does not establish causality.</p></header>
<section><h2>Results at a glance</h2><div class="metrics">%s</div></section>
<section><h2>Multiple-testing results</h2>%s<p class="muted">The nominal threshold is uncorrected. Bonferroni controls family-wise error and Benjamini-Hochberg controls false discovery rate across this result family.</p></section>
<section><h2>Input compatibility and analysis coverage</h2>%s</section>
%s
<section><h2>Complete association results</h2><p class="muted">All %s validated %s are included and initially ordered by %s. Search, sort, or move between pages without changing the result set.</p>%s</section>
<section><h2>How to interpret</h2>%s</section>
<section><h2>Analysis parameters and provenance</h2>%s<p class="muted">The genome build and reference population are explicit declarations; BED/BIM/FAM contents do not independently encode them.</p><div class="downloads"><strong>Output files</strong><ul>%s</ul></div></section>
</main>
<script>
document.querySelectorAll('.paginated-results').forEach(function(container) {
  const rows = JSON.parse(document.getElementById(container.dataset.source).textContent);
  const pageSize = Number(container.dataset.pageSize);
  const body = container.querySelector('tbody');
  const search = container.querySelector('.result-search');
  const previous = container.querySelector('.previous-page');
  const next = container.querySelector('.next-page');
  const status = container.querySelector('.page-status');
  let filtered = rows.slice(); let page = 0; let sortColumn = null; let ascending = true;
  function compare(left, right) {
    const leftNumber = Number(String(left).replaceAll(',', ''));
    const rightNumber = Number(String(right).replaceAll(',', ''));
    if (Number.isFinite(leftNumber) && Number.isFinite(rightNumber)) return leftNumber - rightNumber;
    return String(left).localeCompare(String(right), undefined, {numeric:true});
  }
  function render() {
    const pageCount = Math.max(1, Math.ceil(filtered.length / pageSize)); page = Math.min(page, pageCount - 1);
    const start = page * pageSize; const end = Math.min(start + pageSize, filtered.length); body.replaceChildren();
    filtered.slice(start, end).forEach(function(row) { const tr = document.createElement('tr'); row.forEach(function(value) { const td = document.createElement('td'); td.textContent = value; tr.appendChild(td); }); body.appendChild(tr); });
    status.textContent = filtered.length ? `${start + 1}–${end} of ${filtered.length.toLocaleString()} · page ${page + 1} of ${pageCount}` : '0 matching rows';
    previous.disabled = page === 0; next.disabled = page + 1 >= pageCount;
  }
  search.addEventListener('input', function() { const term = search.value.trim().toLocaleLowerCase(); filtered = term ? rows.filter(function(row) { return row.some(function(value) { return String(value).toLocaleLowerCase().includes(term); }); }) : rows.slice(); if (sortColumn !== null) filtered.sort(function(left, right) { const order = compare(left[sortColumn], right[sortColumn]); return ascending ? order : -order; }); page = 0; render(); });
  previous.addEventListener('click', function() { page -= 1; render(); }); next.addEventListener('click', function() { page += 1; render(); });
  container.querySelectorAll('th button').forEach(function(button) { button.addEventListener('click', function() { const column = Number(button.dataset.column); ascending = sortColumn === column ? !ascending : true; sortColumn = column; filtered.sort(function(left, right) { const order = compare(left[column], right[column]); return ascending ? order : -order; }); page = 0; render(); }); }); render();
});
</script></body></html>""" % (
        escape(title), escape(dataset_id),
        "var(--amber)" if summary["warnings"] else "var(--green)",
        escape(title), escape(dataset_id), escape(status), metrics,
        _html_table(
            significance_rows,
            ("Procedure", "Threshold", "Significant associations"),
            precision,
        ),
        _html_table(coverage_rows, ("Assessment", "Result"), precision),
        warning_section,
        format_count(len(rows)), escape(summary["unit_label"]),
        escape(summary["p_value_column"]),
        _paginated_gcta_table(
            rows, columns, page_size=module_config.html_report.page_size,
        ),
        _html_table(interpretation_rows, ("Topic", "Interpretation"), precision),
        _html_table(parameters, ("Parameter", "Value"), precision), downloads,
    )
    return write_html_report(document, destination)


__all__ = [
    "build_gcta_scientific_summary",
    "gcta_gene_overlap_outcome_fields",
    "gcta_gene_reference_outcome_fields",
    "gcta_gmt_outcome_fields",
    "gcta_gmt_preparation_outcome_fields",
    "gcta_ld_reference_outcome_fields",
    "gcta_ma_input_outcome_fields",
    "gcta_set_source_outcome_fields",
    "record_gcta_scientific_summary",
    "render_gcta_scientific_summary",
    "write_gcta_html_report",
]
