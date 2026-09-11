"""Machine-readable and human-readable reports for validated LD clumping.

The builders are presentation-only. They consume the counts returned by the
region, standard, and COJO algorithms and never reopen a scientific result file.
Consequently, the terminal, CSV, and HTML surfaces cannot diverge through a
second scientific calculation.
"""

from __future__ import annotations

from html import escape
import json
import os
from pathlib import Path
from typing import Any, Mapping, Sequence

from postgwas.core.io.reports import write_delimited_report, write_html_report
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.modules.ld_clumping.common import (
    format_reference_exclusion_reasons,
    reference_exclusion_analysis_limitation,
    reference_exclusion_reason_items,
    reference_exclusion_reason_label,
)


SUMMARY_COLUMNS = (
    "record_type",
    "dataset_id",
    "genome_build",
    "population",
    "method",
    "status",
    "order",
    "chromosome",
    "detail_name",
    "detail_value",
    "annotated_variants",
    "ld_blocks",
    "significant_blocks",
    "genome_wide_significant_outside_ld_regions",
    "input_significant_variants",
    "analysed_significant_variants",
    "skipped_significant_variants",
    "independent_significant_snps",
    "lead_snps",
    "genomic_risk_loci",
    "cojo_summary_variants",
    "cojo_reference_overlap_variants",
    "cojo_reference_overlap_fraction",
    "cojo_reference_missing_variants",
    "cojo_reference_allele_mismatch_variants",
    "cojo_reference_samples",
    "cojo_selected_signals",
    "cojo_genomic_loci",
    "cojo_warnings",
    "cojo_configured_parallel_output_contract",
    "cojo_parallel_output_contract",
    "cojo_conditional_reconstruction_status",
    "reference_coverage_status",
    "reference_exclusions",
    "reference_exclusions_within_reported_locus_boundaries",
    "reference_exclusions_outside_reported_locus_boundaries",
    "warnings",
    "reference_ld_rows_before_maf",
    "reference_ld_rows_retained",
    "low_maf_partner_rows_excluded",
    "low_maf_partner_variants_excluded",
    "lead_pvalue",
    "candidate_pvalue",
    "clump_r2",
    "lead_r2",
    "window_kb",
    "merge_distance_bp",
    "cojo_merge_distance_bp",
    "cojo_index_pvalue",
    "cojo_model_window_kb",
    "cojo_chromosome",
    "minimum_reference_maf",
    "remove_mhc",
    "missing_index_action",
    "missing_chromosome_action",
    "output_name",
    "output_path",
)


def build_ld_clumping_summary(
    *,
    dataset_id: str,
    input_vcf: str | Path,
    reference_directory: str | Path | None,
    module,
    methods: Sequence[str],
    status: str,
    region_result: Mapping[str, Any] | None,
    standard_result: Mapping[str, Any] | None,
    output_paths: Mapping[str, str | Path | None],
    threads: int,
    memory_gb: float,
    cojo_result: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build the one evidence object shared by every report surface."""
    return {
        "schema_version": 1,
        "dataset_id": str(dataset_id),
        "genome_build": module.genome_build.value,
        "population": module.population.value,
        "methods": [str(method) for method in methods],
        "status": str(status),
        "inputs": {
            "vcf": str(Path(input_vcf).expanduser().resolve()),
            "ld_reference_directory": (
                None
                if reference_directory is None
                else str(Path(reference_directory).expanduser().resolve())
            ),
            "cojo_reference_prefix": (
                (cojo_result or {}).get("reference_prefix")
            ),
        },
        "scientific_policy": {
            "lead_pvalue": module.lead_pvalue,
            "candidate_pvalue": module.candidate_pvalue,
            "clump_r2": module.clump_r2,
            "lead_r2": module.lead_r2,
            "window_kb": module.window_kb,
            "merge_distance_bp": module.merge_distance_bp,
            "minimum_reference_maf": module.minimum_reference_maf,
            "remove_mhc": module.remove_mhc,
            "missing_index_action": module.missing_index_action,
            "missing_chromosome_action": module.missing_chromosome_action,
            "cojo_merge_distance_bp": module.cojo.merge_distance_bp,
            "cojo_index_pvalue": module.cojo.index_pvalue,
            "cojo_parallel_output_contract": (
                module.cojo.parallel_output_contract
            ),
            "reference_orientation": module.reference.orientation,
            "reference_format_version": module.reference.format_version,
        },
        "execution": {
            "threads": int(threads),
            "memory_gb": float(memory_gb),
        },
        "region": dict(region_result or {}),
        "standard": dict(standard_result or {}),
        "cojo": dict(cojo_result or {}),
        "outputs": {
            str(name): None if value is None else str(value)
            for name, value in output_paths.items()
        },
    }


def ld_clumping_summary_records(
    summary: Mapping[str, Any],
) -> list[dict[str, Any]]:
    """Return run, method, chromosome, reference, and output CSV records."""
    policy = summary["scientific_policy"]
    common = {
        "dataset_id": summary["dataset_id"],
        "genome_build": summary["genome_build"],
        "population": summary["population"],
        "lead_pvalue": policy["lead_pvalue"],
        "candidate_pvalue": policy["candidate_pvalue"],
        "clump_r2": policy["clump_r2"],
        "lead_r2": policy["lead_r2"],
        "window_kb": policy["window_kb"],
        "merge_distance_bp": policy["merge_distance_bp"],
        "minimum_reference_maf": policy["minimum_reference_maf"],
        "remove_mhc": policy["remove_mhc"],
        "missing_index_action": policy["missing_index_action"],
        "missing_chromosome_action": policy["missing_chromosome_action"],
        "cojo_merge_distance_bp": policy["cojo_merge_distance_bp"],
        "cojo_index_pvalue": policy["cojo_index_pvalue"],
        "cojo_configured_parallel_output_contract": policy[
            "cojo_parallel_output_contract"
        ],
    }
    region = summary.get("region") or {}
    standard = summary.get("standard") or {}
    cojo = summary.get("cojo") or {}
    records = [{
        **common,
        "record_type": "overall",
        "method": ",".join(summary["methods"]),
        "status": summary["status"],
        "annotated_variants": region.get("annotated_variants"),
        "ld_blocks": region.get("ld_blocks"),
        "significant_blocks": region.get("significant_blocks"),
        "genome_wide_significant_outside_ld_regions": region.get(
            "genome_wide_significant_outside_ld_regions"
        ),
        "input_significant_variants": standard.get(
            "input_significant_variants"
        ),
        "analysed_significant_variants": standard.get(
            "significant_variants"
        ),
        "skipped_significant_variants": standard.get(
            "skipped_significant_variants"
        ),
        "independent_significant_snps": standard.get(
            "independent_significant_snps"
        ),
        "lead_snps": standard.get("lead_snps"),
        "genomic_risk_loci": standard.get("genomic_risk_loci"),
        "reference_coverage_status": standard.get(
            "reference_coverage_status"
        ),
        "reference_exclusions": standard.get("reference_exclusions"),
        "reference_exclusions_within_reported_locus_boundaries": standard.get(
            "reference_exclusions_within_reported_locus_boundaries"
        ),
        "reference_exclusions_outside_reported_locus_boundaries": standard.get(
            "reference_exclusions_outside_reported_locus_boundaries"
        ),
        "warnings": standard.get("warnings"),
        "cojo_summary_variants": cojo.get("summary_variants"),
        "cojo_reference_overlap_variants": cojo.get(
            "reference_overlap_variants"
        ),
        "cojo_reference_overlap_fraction": cojo.get(
            "reference_overlap_fraction"
        ),
        "cojo_reference_missing_variants": cojo.get(
            "reference_missing_variants"
        ),
        "cojo_reference_allele_mismatch_variants": cojo.get(
            "reference_allele_mismatch_variants"
        ),
        "cojo_reference_samples": cojo.get("reference_samples"),
        "cojo_selected_signals": cojo.get("selected_signals"),
        "cojo_genomic_loci": cojo.get("genomic_loci"),
        "cojo_warnings": len(cojo.get("warnings") or ()),
        "cojo_parallel_output_contract": cojo.get(
            "parallel_output_contract"
        ),
        "cojo_conditional_reconstruction_status": cojo.get(
            "conditional_reconstruction_status"
        ),
        "cojo_model_window_kb": cojo.get("gcta_model_window_kb"),
        "cojo_chromosome": cojo.get("gcta_chromosome"),
    }]
    if region:
        records.append({
            **common,
            "record_type": "method",
            "method": "region",
            "status": region.get("status"),
            "annotated_variants": region.get("annotated_variants"),
            "ld_blocks": region.get("ld_blocks"),
            "significant_blocks": region.get("significant_blocks"),
            "genome_wide_significant_outside_ld_regions": region.get(
                "genome_wide_significant_outside_ld_regions"
            ),
        })
    if standard:
        records.append({
            **common,
            "record_type": "method",
            "method": "standard",
            "status": standard.get("status"),
            "input_significant_variants": standard.get(
                "input_significant_variants"
            ),
            "analysed_significant_variants": standard.get(
                "significant_variants"
            ),
            "skipped_significant_variants": standard.get(
                "skipped_significant_variants"
            ),
            "independent_significant_snps": standard.get(
                "independent_significant_snps"
            ),
            "lead_snps": standard.get("lead_snps"),
            "genomic_risk_loci": standard.get("genomic_risk_loci"),
            "reference_coverage_status": standard.get(
                "reference_coverage_status"
            ),
            "reference_exclusions": standard.get("reference_exclusions"),
            "reference_exclusions_within_reported_locus_boundaries": (
                standard.get(
                    "reference_exclusions_within_reported_locus_boundaries"
                )
            ),
            "reference_exclusions_outside_reported_locus_boundaries": (
                standard.get(
                    "reference_exclusions_outside_reported_locus_boundaries"
                )
            ),
            "warnings": standard.get("warnings"),
            "reference_ld_rows_before_maf": standard.get(
                "reference_ld_rows_before_maf"
            ),
            "reference_ld_rows_retained": standard.get(
                "reference_ld_rows_retained"
            ),
            "low_maf_partner_rows_excluded": standard.get(
                "low_maf_partner_rows_excluded"
            ),
            "low_maf_partner_variants_excluded": standard.get(
                "low_maf_partner_variants_excluded"
            ),
        })
    if cojo:
        records.append({
            **common,
            "record_type": "method",
            "method": "cojo-slct",
            "status": cojo.get("status"),
            "cojo_summary_variants": cojo.get("summary_variants"),
            "cojo_reference_overlap_variants": cojo.get(
                "reference_overlap_variants"
            ),
            "cojo_reference_overlap_fraction": cojo.get(
                "reference_overlap_fraction"
            ),
            "cojo_reference_missing_variants": cojo.get(
                "reference_missing_variants"
            ),
            "cojo_reference_allele_mismatch_variants": cojo.get(
                "reference_allele_mismatch_variants"
            ),
            "cojo_reference_samples": cojo.get("reference_samples"),
            "cojo_selected_signals": cojo.get("selected_signals"),
            "cojo_genomic_loci": cojo.get("genomic_loci"),
            "cojo_warnings": len(cojo.get("warnings") or ()),
            "cojo_parallel_output_contract": cojo.get(
                "parallel_output_contract"
            ),
            "cojo_conditional_reconstruction_status": cojo.get(
                "conditional_reconstruction_status"
            ),
            "cojo_model_window_kb": cojo.get("gcta_model_window_kb"),
            "cojo_chromosome": cojo.get("gcta_chromosome"),
        })
        for order, warning in enumerate(cojo.get("warnings") or (), 1):
            records.append({
                **common,
                "record_type": "cojo_warning",
                "method": "cojo-slct",
                "status": "warning",
                "order": order,
                "detail_name": "warning_%d" % order,
                "detail_value": str(warning),
                "cojo_summary_variants": cojo.get("summary_variants"),
                "cojo_reference_overlap_variants": cojo.get(
                    "reference_overlap_variants"
                ),
                "cojo_reference_overlap_fraction": cojo.get(
                    "reference_overlap_fraction"
                ),
                "cojo_reference_missing_variants": cojo.get(
                    "reference_missing_variants"
                ),
                "cojo_reference_allele_mismatch_variants": cojo.get(
                    "reference_allele_mismatch_variants"
                ),
                "cojo_reference_samples": cojo.get("reference_samples"),
                "cojo_warnings": len(cojo.get("warnings") or ()),
            })
    for order, chromosome in enumerate(
        standard.get("chromosome_results") or (), 1,
    ):
        records.append({
            **common,
            "record_type": "chromosome",
            "method": "standard",
            "status": chromosome.get("status"),
            "order": order,
            "chromosome": chromosome.get("chrom"),
            "input_significant_variants": chromosome.get(
                "input_significant", chromosome.get("significant")
            ),
            "analysed_significant_variants": chromosome.get("significant"),
            "independent_significant_snps": chromosome.get("independent"),
            "lead_snps": chromosome.get("lead"),
            "genomic_risk_loci": chromosome.get("loci"),
            "warnings": chromosome.get("warnings"),
            "reference_ld_rows_before_maf": chromosome.get(
                "reference_ld_rows_before_maf"
            ),
            "reference_ld_rows_retained": chromosome.get(
                "reference_ld_rows_retained"
            ),
            "low_maf_partner_rows_excluded": chromosome.get(
                "low_maf_partner_rows_excluded"
            ),
            "low_maf_partner_variants_excluded": chromosome.get(
                "low_maf_partner_variants_excluded"
            ),
        })
    for order, detail in enumerate(
        standard.get("missing_reference_details") or (), 1,
    ):
        records.append({
            **common,
            "record_type": "missing_reference",
            "method": "standard",
            "status": "warning_skip",
            "order": order,
            "chromosome": detail.get("chromosome"),
            "detail_name": "missing_resources",
            "detail_value": json.dumps(
                detail.get("missing_resources") or (), ensure_ascii=False
            ),
            "skipped_significant_variants": detail.get(
                "significant_variants"
            ),
        })
    for order, (reason, count) in enumerate(
        (standard.get("reference_exclusion_reasons") or {}).items(), 1,
    ):
        records.append({
            **common,
            "record_type": "reference_exclusion",
            "method": "standard",
            "order": order,
            "detail_name": reason,
            "detail_value": count,
            "reference_exclusions": count,
        })
    for order, (name, path) in enumerate(summary["outputs"].items(), 1):
        if path:
            records.append({
                **common,
                "record_type": "output",
                "status": "generated",
                "order": order,
                "output_name": name,
                "output_path": path,
            })
    return records


def write_ld_clumping_summary_csv(
    summary: Mapping[str, Any], destination: str | Path,
) -> Path:
    """Write the complete LD-clumping summary atomically as CSV."""
    return write_delimited_report(
        ld_clumping_summary_records(summary),
        destination,
        fieldnames=SUMMARY_COLUMNS,
        delimiter=",",
        null_value="",
    )


def _count(value: Any) -> str:
    if value is None:
        return "not applicable"
    try:
        return f"{int(value):,}"
    except (TypeError, ValueError):
        return str(value)


def _chromosomes(values: Sequence[Any]) -> str:
    labels = ["chr%s" % value for value in values]
    return ", ".join(labels) if labels else "none"


def format_cojo_reference_match(cojo: Mapping[str, Any]) -> str:
    """Explain the validated GWAS/LD-reference overlap without COJO jargon."""
    matched = cojo.get("reference_overlap_variants")
    total = cojo.get("summary_variants")
    fraction = cojo.get("reference_overlap_fraction")
    if matched is None and total is None:
        return "not available"
    if matched is None:
        count_text = "matched count unavailable out of %s" % _count(total)
    elif total is None:
        count_text = "%s matched · total unavailable" % _count(matched)
    else:
        count_text = "%s of %s" % (_count(matched), _count(total))
    if fraction is None:
        return "%s · percentage unavailable" % count_text
    return "%s · %.2f%% matched by SNP ID and allele pair" % (
        count_text,
        100 * float(fraction),
    )


def format_cojo_output_scope(cojo: Mapping[str, Any]) -> str:
    """Describe which scientifically distinct COJO result tables were made."""
    status = cojo.get("conditional_reconstruction_status")
    if status == "not_requested":
        return (
            "selected signals and joint statistics (.jma); "
            "full-genome conditional table (.cma) not generated"
        )
    if status in {"completed", "native_single_command"}:
        return (
            "selected signals and joint statistics (.jma); "
            "full-genome conditional table (.cma) generated"
        )
    if status == "not_applicable_no_signals":
        return "no signals selected; conditional table not applicable"
    return "not reported"


def render_ld_clumping_summary(
    summary: Mapping[str, Any], *, label_width: int,
) -> str:
    """Render one concise terminal summary with one configured value column."""
    indent = 6
    region = summary.get("region") or {}
    standard = summary.get("standard") or {}
    cojo = summary.get("cojo") or {}
    status = summary["status"]
    lines = [
        "",
        screen_line("analysis", "LD clumping summary", indent=2),
        screen_field(
            "info", "Dataset", summary["dataset_id"], indent=indent,
            label_width=label_width,
        ),
        screen_field(
            "genetic", "Genome build", summary["genome_build"], indent=indent,
            label_width=label_width,
        ),
        screen_field(
            "genetic", "Population", summary["population"], indent=indent,
            label_width=label_width,
        ),
        screen_field(
            "analysis", "Methods", ", ".join(summary["methods"]),
            indent=indent, label_width=label_width,
        ),
        screen_field(
            "warning" if status == "partial_reference" else "success",
            "Run status",
            (
                "partial reference · results are not genome-wide"
                if status == "partial_reference" else status
            ),
            indent=indent,
            label_width=label_width,
        ),
    ]
    if region:
        lines.extend((
            "",
            screen_line("genetic", "Annotated-region method", indent=4),
            screen_field(
                "count", "Annotated variants",
                _count(region.get("annotated_variants")), indent=indent,
                label_width=label_width,
            ),
            screen_field(
                "genetic", "LD blocks", _count(region.get("ld_blocks")),
                indent=indent, label_width=label_width,
            ),
            screen_field(
                "success", "Significant blocks",
                _count(region.get("significant_blocks")), indent=indent,
                label_width=label_width,
            ),
            screen_field(
                (
                    "warning"
                    if region.get(
                        "genome_wide_significant_outside_ld_regions"
                    )
                    else "info"
                ),
                "Significant outside LD regions",
                "%s · missing %s_LDblock annotation"
                % (
                    _count(
                        region.get(
                            "genome_wide_significant_outside_ld_regions"
                        )
                    ),
                    summary["population"],
                ),
                indent=indent,
                label_width=label_width,
            ),
        ))
    if standard:
        coverage = standard.get("reference_coverage_status")
        exclusion_reasons = standard.get("reference_exclusion_reasons") or {}
        skipped_chromosomes = standard.get("skipped_chromosomes") or ()
        coverage_value = coverage or "not reported"
        if coverage == "partial" and skipped_chromosomes:
            reference_word = (
                "reference" if len(skipped_chromosomes) == 1 else "references"
            )
            coverage_value += (
                f" · {_chromosomes(skipped_chromosomes)} LD "
                f"{reference_word} unavailable"
            )
        skipped_significant = standard.get("skipped_significant_variants")
        lines.extend((
            "",
            screen_line("genetic", "FUMA-style standard method", indent=4),
            "",
            screen_field(
                "warning" if coverage == "partial" else "success",
                "Reference coverage", coverage_value,
                indent=indent, label_width=label_width,
            ),
            "",
            screen_field(
                "count", "GWS variants supplied",
                _count(standard.get("input_significant_variants")),
                indent=indent, label_width=label_width,
            ),
            screen_field(
                "success", "GWS variants passing LD-reference checks",
                _count(standard.get("significant_variants")), indent=indent,
                label_width=label_width,
            ),
            screen_field(
                "warning" if skipped_significant else "info",
                "GWS variants excluded",
                _count(skipped_significant),
                indent=indent, label_width=label_width,
            ),
        ))
        for reason, count in reference_exclusion_reason_items(exclusion_reasons):
            value = _count(count)
            if reason == "missing_chromosome_reference":
                if skipped_chromosomes:
                    value += f" · {_chromosomes(skipped_chromosomes)}"
            lines.append(
                screen_field(
                    "warning",
                    reference_exclusion_reason_label(reason),
                    value,
                    indent=indent + 2,
                    label_width=label_width - 2,
                )
            )
        lines.extend((
            "",
            screen_field(
                "genetic", "Independent significant SNPs",
                _count(standard.get("independent_significant_snps")),
                indent=indent, label_width=label_width,
            ),
            screen_field(
                "genetic", "Lead SNPs", _count(standard.get("lead_snps")),
                indent=indent, label_width=label_width,
            ),
            screen_field(
                "genetic", "Genomic risk loci",
                _count(standard.get("genomic_risk_loci")), indent=indent,
                label_width=label_width,
            ),
        ))
        if exclusion_reasons:
            lines.extend((
                "",
                screen_line(
                    "analysis",
                    "Excluded variants relative to reported loci",
                    indent=indent,
                ),
                screen_field(
                    "warning",
                    "Outside all locus boundaries",
                    _count(
                        standard.get(
                            "reference_exclusions_outside_reported_locus_boundaries"
                        )
                    ),
                    indent=indent + 2,
                    label_width=label_width - 2,
                ),
                screen_field(
                    "info",
                    "Within reported locus boundaries",
                    "%s · coordinate overlap only"
                    % _count(
                        standard.get(
                            "reference_exclusions_within_reported_locus_boundaries"
                        )
                    ),
                    indent=indent + 2,
                    label_width=label_width - 2,
                ),
            ))
            lines.append("")
            lines.append(
                screen_field(
                    "warning",
                    "Analysis limitation",
                    reference_exclusion_analysis_limitation(
                        skipped_significant
                    ),
                    indent=indent,
                    label_width=label_width,
                )
            )
        if standard.get("other_warnings"):
            lines.append(
                screen_field(
                    "warning",
                    "Other diagnostic warnings",
                    f"{_count(standard['other_warnings'])} · see detailed log",
                    indent=indent,
                    label_width=label_width,
                )
            )
    if cojo:
        warnings = cojo.get("warnings") or ()
        exclusions = cojo.get("exclusions") or {}
        overlap_value = format_cojo_reference_match(cojo)
        overlap_kind = (
            "warning"
            if cojo.get("reference_missing_variants")
            or cojo.get("reference_allele_mismatch_variants")
            else "success"
        )
        lines.extend((
            "",
            screen_line("genetic", "GCTA-COJO stepwise method", indent=4),
            "",
            screen_field(
                "count", "GWAS variants formatted",
                _count(cojo.get("summary_variants")), indent=indent,
                label_width=label_width,
            ),
            screen_field(
                overlap_kind, "GWAS variants matched to LD reference",
                overlap_value,
                indent=indent, label_width=label_width,
            ),
            screen_field(
                "count", "LD-reference samples",
                _count(cojo.get("reference_samples")), indent=indent,
                label_width=label_width,
            ),
            "",
            screen_field(
                "genetic", "COJO-selected signals",
                _count(cojo.get("selected_signals")), indent=indent,
                label_width=label_width,
            ),
            screen_field(
                "genetic", "Physical loci",
                _count(cojo.get("genomic_loci")), indent=indent,
                label_width=label_width,
            ),
            screen_field(
                "info", "Physical grouping distance",
                "%s bp · applied after GCTA model fitting"
                % _count(cojo.get("merge_distance_bp")),
                indent=indent, label_width=label_width,
            ),
            screen_field(
                "info", "GCTA model window",
                "%s kb · not a locus boundary"
                % _count(cojo.get("gcta_model_window_kb")),
                indent=indent, label_width=label_width,
            ),
            screen_field(
                "info", "Locus index statistic",
                "%s p value" % cojo.get("index_pvalue"),
                indent=indent, label_width=label_width,
            ),
            screen_field(
                "info", "COJO output scope",
                format_cojo_output_scope(cojo),
                indent=indent, label_width=label_width,
            ),
        ))
        if cojo.get("gcta_chromosome") is not None:
            lines.append(screen_field(
                "warning", "GCTA chromosome restriction",
                "chr%s · COJO result is chromosome-specific"
                % cojo["gcta_chromosome"],
                indent=indent, label_width=label_width,
            ))
        if exclusions.get("mhc_exclusions"):
            lines.append(screen_field(
                "info", "Reference variants excluded in MHC",
                _count(exclusions.get("mhc_exclusions")), indent=indent,
                label_width=label_width,
            ))
        lines.append(screen_field(
            "warning" if warnings else "info",
            "COJO warnings",
            _count(len(warnings)) if warnings else "none",
            indent=indent, label_width=label_width,
        ))
        for index, warning in enumerate(warnings, 1):
            lines.append(screen_field(
                "warning",
                "Warning %d" % index,
                warning,
                indent=indent + 2,
                label_width=label_width - 2,
            ))
    outputs = summary["outputs"]
    lines.extend((
        "",
        screen_line("analysis", "Saved reports", indent=4),
        screen_field(
            "success", "Summary CSV",
            os.path.basename(str(outputs["summary_csv"])), indent=indent,
            label_width=label_width,
        ),
        screen_field(
            "success", "Detailed HTML report",
            os.path.basename(str(outputs["html_report"])), indent=indent,
            label_width=label_width,
        ),
    ))
    return "\n".join(lines)


def _html_value(value: Any) -> str:
    if value is None or value == "":
        return "Not available"
    if value is True:
        return "Yes"
    if value is False:
        return "No"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        return f"{value:,.8g}"
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return str(value)


def _cell(value: Any) -> str:
    return escape(_html_value(value))


def _path_cell(value: Any) -> str:
    if value is None or value == "":
        return '<span class="muted">Not generated</span>'
    path = Path(str(value)).expanduser()
    try:
        href = path.resolve().as_uri()
    except (OSError, ValueError):
        href = ""
    label = escape(str(value))
    if not href:
        return '<span class="path">%s</span>' % label
    return '<a class="path" href="%s">%s</a>' % (
        escape(href, quote=True), label,
    )


def _table(
    headers: Sequence[str],
    rows: Sequence[Sequence[Any]],
    *,
    table_id: str | None = None,
) -> str:
    head = "".join(
        '<th scope="col">%s</th>' % escape(header) for header in headers
    )
    body = "".join(
        "<tr>%s</tr>" % "".join("<td>%s</td>" % _cell(value) for value in row)
        for row in rows
    )
    if not rows:
        body = '<tr><td colspan="%d" class="muted">No rows</td></tr>' % len(
            headers
        )
    id_attribute = (
        "" if table_id is None
        else ' id="%s"' % escape(table_id, quote=True)
    )
    return (
        '<div class="table-wrap"><table%s><thead><tr>%s</tr></thead>'
        '<tbody>%s</tbody></table></div>' % (id_attribute, head, body)
    )


def _key_values(
    rows: Sequence[tuple[str, Any]], *, path_labels: Sequence[str] = (),
) -> str:
    paths = set(path_labels)
    body = "".join(
        '<tr><th scope="row">%s</th><td>%s</td></tr>'
        % (
            escape(label),
            _path_cell(value) if label in paths else _cell(value),
        )
        for label, value in rows
    )
    return (
        '<div class="table-wrap"><table class="key-value"><tbody>%s</tbody>'
        '</table></div>' % body
    )


def _searchable_table_section(
    *,
    title: str,
    note: str,
    columns: Sequence[str],
    rows: Sequence[Sequence[Any]],
    table_id: str,
    output_label: str,
    output_file: Any,
    total_rows: Any = None,
) -> str:
    """Render one self-contained searchable table with an audit-file link."""
    count_id = "%s-count" % table_id
    return (
        '<section class="result-table"><h2>%s</h2>'
        '<p class="section-note">%s</p>'
        '<div class="result-tools"><label for="%s-filter">Search this '
        'table</label><input id="%s-filter" type="search" '
        'data-table-filter="%s" data-count-target="%s" '
        'placeholder="Chromosome, locus, variant, rsID…">'
        '<span id="%s">%s rows</span></div>%s'
        '<p class="full-result-link"><strong>%s:</strong> %s</p>'
        '</section>'
        % (
            escape(title),
            escape(note),
            table_id,
            table_id,
            table_id,
            count_id,
            count_id,
            _count(len(rows) if total_rows is None else total_rows),
            _table(columns, rows, table_id=table_id),
            escape(output_label),
            _path_cell(output_file),
        )
    )


def _result_table_sections(
    result: Mapping[str, Any], *, table_prefix: str, method_label: str,
) -> list[str]:
    """Render searchable compact views of validated method result tables."""
    sections = []
    for index, (name, table) in enumerate(
        (result.get("_report_tables") or {}).items(), 1,
    ):
        columns = table.get("columns") or ()
        rows = table.get("rows") or ()
        if not columns:
            continue
        title = name.replace("_", " ").title().replace("Snps", "SNPs")
        sections.append(
            _searchable_table_section(
                title=title,
                note=(
                    "Browse all %s %s result rows directly in this report. The "
                    "compact columns are configured in the canonical "
                    "LD-clumping YAML; the linked TSV retains every scientific "
                    "column."
                    % (
                        _count(table.get("total_rows", len(rows))),
                        method_label,
                    )
                ),
                columns=columns,
                rows=rows,
                table_id="%s-result-%d" % (table_prefix, index),
                output_label="Complete TSV",
                output_file=table.get("output_file"),
                total_rows=table.get("total_rows", len(rows)),
            )
        )
    return sections


def _reference_exclusion_locus_section(
    standard: Mapping[str, Any],
) -> str:
    """Render every excluded index and its coordinate-only locus relation."""
    details = standard.get("_reference_exclusion_details") or ()
    if not details:
        return ""
    rows = []
    for detail in details:
        status = detail.get("reported_locus_boundary_status")
        status_label = {
            "inside_reported_locus_boundary": (
                "Inside reported locus boundary"
            ),
            "outside_all_reported_locus_boundaries": (
                "Outside all reported locus boundaries"
            ),
        }.get(status, str(status or "Not available"))
        boundary_chromosome = detail.get("overlapping_locus_chromosome")
        boundary_start = detail.get("overlapping_locus_start")
        boundary_end = detail.get("overlapping_locus_end")
        matching_boundary = "Not applicable"
        if (
            boundary_chromosome
            and boundary_start is not None
            and boundary_end is not None
        ):
            matching_boundary = "chr%s:%s–%s" % (
                boundary_chromosome,
                f"{int(boundary_start):,}",
                f"{int(boundary_end):,}",
            )
        rows.append((
            detail.get("chromosome"),
            detail.get("position"),
            detail.get("canonical_id"),
            detail.get("input_variant_id"),
            reference_exclusion_reason_label(detail.get("reason")),
            status_label,
            detail.get("overlapping_genomic_loci"),
            matching_boundary,
            detail.get("locus_membership_interpretation"),
        ))
    return _searchable_table_section(
        title="Excluded indexes and reported locus boundaries",
        note=(
            "%s excluded indexes are outside every reported locus boundary; "
            "%s fall inside a boundary by coordinate only. Coordinate overlap "
            "does not restore reference validation, LD membership, or the "
            "ability to seed a locus."
            % (
                _count(
                    standard.get(
                        "reference_exclusions_outside_reported_locus_boundaries"
                    )
                ),
                _count(
                    standard.get(
                        "reference_exclusions_within_reported_locus_boundaries"
                    )
                ),
            )
        ),
        columns=(
            "Chromosome",
            "Position",
            "Canonical variant",
            "Input ID",
            "Exclusion reason",
            "Boundary classification",
            "Matching locus ID",
            "Matching locus boundary",
            "Scientific interpretation",
        ),
        rows=rows,
        table_id="reference-exclusion-locus-audit",
        output_label="Complete exclusion audit",
        output_file=standard.get("reference_exclusions_file"),
        total_rows=len(details),
    )


def _region_outside_ld_section(region: Mapping[str, Any]) -> str:
    """Render the complete significant-variant LD-region coverage audit."""
    table = region.get("_significant_outside_ld_region_details") or {}
    columns = table.get("columns") or ()
    if not columns:
        return ""
    display_names = {
        "CHR": "Chromosome",
        "BP": "Position",
        "SNP": "Variant",
        "uniq_id": "Canonical variant",
        "ID": "Input ID",
        "REF": "REF",
        "ALT": "ALT",
        "BETA": "Beta",
        "SE": "SE",
        "AF": "Allele frequency",
        "LP": "-log10(P)",
        "P_value": "P value",
        "LDblock": "LD-block annotation",
        "exclusion_reason": "Reason",
    }
    return _searchable_table_section(
        title="Genome-wide-significant variants outside annotated LD regions",
        note=(
            "%s variants pass the configured genome-wide threshold but lack "
            "the selected population's LD-block annotation after MHC "
            "exclusion. Missing annotation means they were outside the region "
            "method's analysed blocks; it does not prove biological LD "
            "independence."
            % _count(table.get("total_rows", len(table.get("rows") or ())))
        ),
        columns=tuple(display_names.get(column, column) for column in columns),
        rows=table.get("rows") or (),
        table_id="region-significant-outside-ld-regions",
        output_label="Complete outside-region audit",
        output_file=table.get("output_file"),
        total_rows=table.get("total_rows"),
    )


_STYLES = """
:root{--ink:#172033;--muted:#64748b;--line:#dbe3ee;--canvas:#f4f7fb;--panel:#fff;--brand:#4338ca;--brand2:#7c3aed;--ok:#166534;--warn:#92400e}*{box-sizing:border-box}body{margin:0;background:var(--canvas);color:var(--ink);font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}.page-header{padding:36px max(24px,calc((100vw - 1320px)/2));color:#fff;background:linear-gradient(135deg,#312e81,#6d28d9)}.eyebrow{margin:0 0 5px;font-size:12px;font-weight:800;letter-spacing:.14em;text-transform:uppercase;opacity:.86}.page-header h1{margin:0;font-size:clamp(28px,4vw,44px);line-height:1.1}.page-header p{max-width:940px;margin:12px 0 0;color:#ede9fe}.container{max-width:1320px;margin:0 auto;padding:26px 24px 58px}.notice,section{background:var(--panel);border:1px solid var(--line);border-radius:13px;box-shadow:0 7px 24px rgba(15,23,42,.05)}.notice{margin-bottom:20px;padding:14px 16px;border-left:4px solid var(--brand2)}.notice.warn{border-left-color:#d97706;background:#fffbeb}.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(175px,1fr));gap:12px;margin:0 0 22px}.metric{padding:15px;border:1px solid var(--line);border-radius:11px;background:#fff}.metric-label{font-size:12px;font-weight:800;color:var(--muted);letter-spacing:.05em;text-transform:uppercase}.metric-value{margin-top:4px;font-size:22px;font-weight:780;overflow-wrap:anywhere}.badge{display:inline-flex;padding:5px 9px;border-radius:999px;font-size:12px;font-weight:850}.badge.ok{color:var(--ok);background:#dcfce7}.badge.warn{color:var(--warn);background:#fef3c7}section{margin:18px 0;padding:22px}section h2{margin:0 0 5px;font-size:21px}.section-note{margin:0 0 14px;color:var(--muted)}.table-wrap{overflow:auto;border:1px solid var(--line);border-radius:9px}table{width:100%;border-collapse:collapse;background:#fff}th,td{padding:10px 12px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}thead th{background:#f5f3ff;color:#334155;font-size:12px;text-transform:uppercase;letter-spacing:.04em}.key-value th{width:310px;background:#f8fafc;color:#475569}.result-tools{display:flex;align-items:end;gap:10px;flex-wrap:wrap;margin:0 0 12px}.result-tools label{width:100%;font-size:12px;font-weight:800;color:var(--muted);text-transform:uppercase;letter-spacing:.04em}.result-tools input{min-width:min(420px,100%);padding:9px 11px;border:1px solid #b8c4d6;border-radius:8px;background:#fff;color:var(--ink);font:inherit}.result-tools span{color:var(--muted)}.result-table .table-wrap{max-height:68vh}.result-table th,.result-table td{white-space:nowrap}.result-table thead th{position:sticky;top:0;z-index:1}.full-result-link{margin:12px 0 0}.path{overflow-wrap:anywhere;word-break:break-word;color:#5b21b6}.muted{color:var(--muted)}.footer{max-width:1320px;margin:0 auto;padding:0 24px 34px;color:var(--muted);font-size:13px}@media(max-width:700px){.container{padding:16px 10px 40px}section{padding:14px}.key-value th{width:42%}th,td{padding:8px}.page-header{padding:28px 18px}}@media print{body{background:#fff}.page-header{background:#fff;color:var(--ink);padding:0}.page-header p{color:var(--muted)}.container{max-width:none;padding:0}.notice,section{box-shadow:none;break-inside:avoid}.table-wrap{overflow:visible}.result-table .table-wrap{max-height:none}a{color:inherit;text-decoration:none}}
"""


_RESULT_TABLE_SCRIPT = """
<script>
document.querySelectorAll('[data-table-filter]').forEach(function(input) {
  const table = document.getElementById(input.dataset.tableFilter);
  const count = document.getElementById(input.dataset.countTarget);
  if (!table || !count) return;
  const rows = Array.from(table.tBodies[0].rows);
  input.addEventListener('input', function() {
    const query = input.value.trim().toLocaleLowerCase();
    let visible = 0;
    rows.forEach(function(row) {
      const matches = !query || row.textContent.toLocaleLowerCase().includes(query);
      row.hidden = !matches;
      if (matches) visible += 1;
    });
    count.textContent = visible.toLocaleString() + ' of ' +
      rows.length.toLocaleString() + ' rows';
  });
});
</script>
"""


def render_ld_clumping_html(summary: Mapping[str, Any]) -> str:
    """Render a self-contained detailed report from validated result metadata."""
    region = summary.get("region") or {}
    standard = summary.get("standard") or {}
    cojo = summary.get("cojo") or {}
    cards = []
    if region:
        cards.extend((
            ("Annotated variants", _count(region.get("annotated_variants"))),
            ("LD blocks", _count(region.get("ld_blocks"))),
            ("Significant blocks", _count(region.get("significant_blocks"))),
            (
                "Significant outside LD regions",
                _count(
                    region.get(
                        "genome_wide_significant_outside_ld_regions"
                    )
                ),
            ),
        ))
    if standard:
        cards.extend((
            (
                "GWS passing LD-reference checks",
                _count(standard.get("significant_variants")),
            ),
            (
                "GWS excluded",
                _count(standard.get("skipped_significant_variants")),
            ),
            ("Independent significant", _count(standard.get("independent_significant_snps"))),
            ("Lead SNPs", _count(standard.get("lead_snps"))),
            ("Genomic risk loci", _count(standard.get("genomic_risk_loci"))),
            (
                "Excluded outside locus boundaries",
                _count(
                    standard.get(
                        "reference_exclusions_outside_reported_locus_boundaries"
                    )
                ),
            ),
        ))
    if cojo:
        cards.extend((
            (
                "GWAS variants matched to LD reference",
                format_cojo_reference_match(cojo),
            ),
            ("COJO-selected signals", _count(cojo.get("selected_signals"))),
            ("COJO physical loci", _count(cojo.get("genomic_loci"))),
            ("COJO reference samples", _count(cojo.get("reference_samples"))),
        ))
    metrics = '<div class="metrics">%s</div>' % "".join(
        '<div class="metric"><div class="metric-label">%s</div>'
        '<div class="metric-value">%s</div></div>'
        % (escape(label), escape(value)) for label, value in cards
    )
    partial = summary["status"] == "partial_reference"
    if partial:
        reason_text = format_reference_exclusion_reasons(
            standard.get("reference_exclusion_reasons") or {}
        )
        limitation = reference_exclusion_analysis_limitation(
            standard.get("skipped_significant_variants")
        )
        status_notice = (
            '<div class="notice warn"><strong>Partial LD-reference coverage:'
            '</strong> %s. Exclusion reasons: %s. '
            'Reported standard loci are valid for LD-reference-validated content '
            'but are not a genome-wide result.</div>'
            % (
                escape(limitation),
                escape(reason_text),
            )
        )
    elif standard:
        status_notice = (
            '<div class="notice"><strong>Reference coverage:</strong> all '
            'significant chromosomes required by this standard run were '
            'available.</div>'
        )
    elif cojo:
        status_notice = (
            '<div class="notice"><strong>COJO reference:</strong> the GCTA '
            'model used the separately reported PLINK reference. Physical '
            'locus grouping was applied only after signal selection.</div>'
        )
    else:
        status_notice = (
            '<div class="notice"><strong>Annotated-region analysis:</strong> '
            'standard LD-reference clumping was not selected for this run.</div>'
        )
    if cojo.get("conditional_reconstruction_status") == "not_requested":
        status_notice += (
            '<div class="notice"><strong>COJO output scope:</strong> '
            'chromosome-level <code>--cojo-slct</code> produced the selected '
            'signals and their joint statistics. The separate full-genome '
            '<code>--cojo-cond</code> table was not requested; it is not used '
            'to select signals or define the reported physical loci.</div>'
        )
    cojo_warnings = tuple(cojo.get("warnings") or ())
    if cojo_warnings:
        status_notice += (
            '<div class="notice warn"><strong>COJO completed with %s '
            'warning%s:</strong> each reason and its consequence are shown '
            'in the COJO warning table below.</div>'
            % (
                _count(len(cojo_warnings)),
                "" if len(cojo_warnings) == 1 else "s",
            )
        )
    method_rows = []
    if region:
        method_rows.append((
            "region", region.get("status"), region.get("annotated_variants"),
            region.get("significant_blocks"),
            region.get("genome_wide_significant_outside_ld_regions"),
            None, None, None,
        ))
    if standard:
        method_rows.append((
            "standard", standard.get("status"),
            standard.get("significant_variants"), None,
            None,
            standard.get("independent_significant_snps"),
            standard.get("lead_snps"), standard.get("genomic_risk_loci"),
        ))
    if cojo:
        method_rows.append((
            "cojo-slct", cojo.get("status"),
            format_cojo_reference_match(cojo), None, None,
            cojo.get("selected_signals"), None, cojo.get("genomic_loci"),
        ))
    chromosome_rows = [
        (
            row.get("chrom"), row.get("status"),
            row.get("input_significant", row.get("significant")),
            row.get("significant"), row.get("independent"), row.get("lead"),
            row.get("loci"),
            row.get("reference_exclusions", row.get("excluded_indexes", 0)),
            format_reference_exclusion_reasons(
                row.get("reference_exclusion_reasons") or {}
            ),
            row.get("other_warnings", 0),
            row.get("reference_ld_rows_retained"),
        )
        for row in standard.get("chromosome_results") or ()
    ]
    missing_rows = [
        (
            row.get("chromosome"), row.get("significant_variants"),
            row.get("missing_resources"),
        )
        for row in standard.get("missing_reference_details") or ()
    ]
    exclusion_rows = [
        (reference_exclusion_reason_label(reason), count)
        for reason, count in reference_exclusion_reason_items(
            standard.get("reference_exclusion_reasons") or {}
        )
    ]
    policy_rows = tuple(
        (name.replace("_", " ").title(), value)
        for name, value in summary["scientific_policy"].items()
    )
    output_rows = tuple(
        (name.replace("_", " ").title(), path)
        for name, path in summary["outputs"].items()
        if path
    )
    sections = [
        '<section><h2>Method summary</h2><p class="section-note">Region, standard, and COJO are separate scientific definitions; their counts are not pooled. COJO-selected signals are displayed in the independent-signal column but are not FUMA independent significant SNPs.</p>%s</section>'
        % _table(
            (
                "Method", "Status", "Variants analysed / reference-matched",
                "Significant blocks", "Significant outside LD regions",
                "Independent significant", "Lead SNPs", "Risk loci",
            ),
            method_rows,
        ),
    ]
    if region:
        region_outside_section = _region_outside_ld_section(region)
        if region_outside_section:
            sections.append(region_outside_section)
    if standard:
        sections.extend(_result_table_sections(
            standard, table_prefix="standard", method_label="standard",
        ))
        exclusion_locus_section = _reference_exclusion_locus_section(standard)
        if exclusion_locus_section:
            sections.append(exclusion_locus_section)
        sections.extend((
            '<section><h2>Chromosome-level standard clumping</h2><p class="section-note">Counts come directly from each completed chromosome worker. Exclusion reasons count significant index variants, not log messages.</p>%s</section>'
            % _table(
                ("Chromosome", "Status", "GWS input", "GWS passing LD-reference checks", "Independent", "Lead", "Risk loci", "GWS excluded", "Exclusion reasons", "Other diagnostic warnings", "LD rows retained"),
                chromosome_rows,
            ),
            '<section><h2>Reference coverage and exclusions</h2><p class="section-note">Each count is a significant index variant excluded before clumping. A missing chromosome is kept separate from an exact-position, allele, or MAF exclusion.</p>%s<h3>Why significant index variants were skipped</h3>%s<h3>Relationship to final reported loci</h3>%s</section>'
            % (
                _table(
                    (
                        "Chromosome", "Significant variants skipped",
                        "Missing resources",
                    ),
                    missing_rows,
                ),
                _table(("Reason", "Excluded indexes"), exclusion_rows),
                _key_values((
                    (
                        "Outside all reported locus boundaries",
                        standard.get(
                            "reference_exclusions_outside_reported_locus_boundaries"
                        ),
                    ),
                    (
                        "Inside a reported locus boundary",
                        "%s · coordinate overlap only, not LD membership"
                        % _count(
                            standard.get(
                                "reference_exclusions_within_reported_locus_boundaries"
                            )
                        ),
                    ),
                )),
            ),
            '<section><h2>Reference LD filtering</h2>%s</section>'
            % _key_values((
                (
                    "LD rows before runtime MAF",
                    standard.get("reference_ld_rows_before_maf"),
                ),
                (
                    "LD rows retained",
                    standard.get("reference_ld_rows_retained"),
                ),
                (
                    "Low-MAF partner rows excluded",
                    standard.get("low_maf_partner_rows_excluded"),
                ),
                (
                    "Distinct low-MAF partners excluded",
                    standard.get("low_maf_partner_variants_excluded"),
                ),
                (
                    "GWS variants excluded",
                    standard.get("reference_exclusions"),
                ),
                (
                    "Other diagnostic warnings",
                    standard.get("other_warnings"),
                ),
            )),
        ))
    if cojo:
        if cojo_warnings:
            sections.append(
                '<section><h2>COJO warnings</h2><p class="section-note">'
                'The run completed, but these validated warnings affect '
                'coverage or interpretation and should not be hidden behind '
                'a count.</p>%s</section>'
                % _table(
                    ("Warning", "Reason and consequence"),
                    tuple(
                        (index, warning)
                        for index, warning in enumerate(cojo_warnings, 1)
                    ),
                )
            )
        sections.extend(_result_table_sections(
            cojo, table_prefix="cojo", method_label="COJO",
        ))
        sections.append(
            '<section><h2>GCTA-COJO model and physical locus policy</h2>'
            '<p class="section-note">A reference match means that the GWAS '
            'SNP ID exists in the BIM and the complete unordered allele pair '
            'agrees. GCTA selected the signals with joint LD modelling. The '
            'physical grouping below did not refit, add, or remove selected '
            'signals.</p>%s</section>'
            % _key_values((
                ("GWAS variants formatted", cojo.get("summary_variants")),
                (
                    "GWAS variants matched to LD reference",
                    format_cojo_reference_match(cojo),
                ),
                ("LD-reference samples", cojo.get("reference_samples")),
                ("GCTA version", cojo.get("gcta_version")),
                ("COJO output scope", format_cojo_output_scope(cojo)),
                (
                    "Parallel output contract",
                    cojo.get("parallel_output_contract"),
                ),
                ("Executed GCTA command", cojo.get("gcta_command")),
                (
                    "GCTA significance threshold",
                    cojo.get("gcta_significance_threshold"),
                ),
                ("GCTA model window (kb)", cojo.get("gcta_model_window_kb")),
                ("GCTA chromosome restriction", cojo.get("gcta_chromosome")),
                (
                    "GCTA collinearity cutoff",
                    cojo.get("gcta_collinearity_cutoff"),
                ),
                (
                    "GCTA maximum frequency difference",
                    cojo.get("gcta_frequency_difference_max"),
                ),
                (
                    "GCTA reference MAF minimum",
                    cojo.get("gcta_reference_maf_min"),
                ),
                (
                    "Physical grouping distance (bp)",
                    cojo.get("merge_distance_bp"),
                ),
                ("Locus index p value", cojo.get("index_pvalue")),
                ("Warnings", len(cojo.get("warnings") or ())),
                (
                    "Configured SNP exclusions",
                    (cojo.get("exclusions") or {}).get(
                        "configured_exclusions"
                    ),
                ),
                (
                    "Reference variants excluded in MHC",
                    (cojo.get("exclusions") or {}).get("mhc_exclusions"),
                ),
            ))
        )
    sections.extend((
        '<section><h2>Resolved scientific policy</h2><p class="section-note">These schema-validated values determined the analysis. The report does not alter them.</p>%s</section>'
        % _key_values(policy_rows),
        '<section><h2>Inputs and execution provenance</h2>%s</section>'
        % _key_values((
            ("Input GWAS-VCF", summary["inputs"].get("vcf")),
            ("LD reference directory", summary["inputs"].get("ld_reference_directory")),
            ("COJO PLINK reference prefix", summary["inputs"].get("cojo_reference_prefix")),
            ("Threads", summary["execution"].get("threads")),
            ("Memory (GB)", summary["execution"].get("memory_gb")),
        ), path_labels=(
            "Input GWAS-VCF", "LD reference directory",
            "COJO PLINK reference prefix",
        )),
        '<section><h2>Outputs</h2>%s</section>'
        % _key_values(
            output_rows, path_labels=tuple(label for label, _ in output_rows),
        ),
    ))
    section_document = "".join(sections)
    status_class = "warn" if partial else "ok"
    status_label = "PARTIAL REFERENCE" if partial else "COMPLETED"
    return """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>%s — PostGWAS LD clumping report</title>
<style>%s</style>
</head>
<body>
<header class="page-header"><p class="eyebrow">PostGWAS</p><h1>LD clumping report</h1><p>%s · %s · %s · selected methods: %s.</p></header>
<main class="container">
<div class="notice"><strong>Scientific status:</strong> <span class="badge %s">%s</span> Reports reuse validated algorithm outputs and do not recalculate variants, LD, or loci.</div>
%s
%s
%s
</main>
<footer class="footer">PostGWAS LD clumping report · consult the linked CSV, reference-exclusion table, locus tables, detailed method log, and canonical log for auditable evidence.</footer>
%s
</body>
</html>
""" % (
        escape(str(summary["dataset_id"])),
        _STYLES,
        escape(str(summary["dataset_id"])),
        escape(str(summary["genome_build"])),
        escape(str(summary["population"])),
        escape(", ".join(summary["methods"])),
        status_class,
        status_label,
        status_notice,
        metrics,
        section_document,
        _RESULT_TABLE_SCRIPT,
    )


def write_ld_clumping_html_report(
    summary: Mapping[str, Any], destination: str | Path,
) -> Path:
    """Write the self-contained detailed HTML report atomically."""
    return write_html_report(render_ld_clumping_html(summary), destination)


__all__ = [
    "SUMMARY_COLUMNS",
    "build_ld_clumping_summary",
    "format_cojo_output_scope",
    "format_cojo_reference_match",
    "ld_clumping_summary_records",
    "render_ld_clumping_html",
    "render_ld_clumping_summary",
    "write_ld_clumping_html_report",
    "write_ld_clumping_summary_csv",
]
