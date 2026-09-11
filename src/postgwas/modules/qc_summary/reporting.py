"""Shared screen report for genotype-free GWAS-VCF quality control."""

from __future__ import annotations

from contextvars import ContextVar
from html import escape
import json
import os
from pathlib import Path
from typing import Any, Mapping, Sequence

from rich.cells import cell_len

from postgwas.core.validation_reporting import (
    consolidate_validation_fields, flush_file_validation_display,
)
from postgwas.core.ui.screen import (
    SUMMARY_CARD_LABEL_WIDTH,
    SYMBOLS,
    render_action_plan,
    screen_field,
    screen_line,
    screen_section_lines,
    screen_summary_card,
)


_QC_DEFAULT_COLON_INDEX = 66
_QC_BASE_INDENT = 8
_QC_TERMINAL_LABEL_WIDTH: ContextVar[int | None] = ContextVar(
    "qc_terminal_label_width", default=None,
)
QC_SUMMARY_COLUMNS = (
    "record_type",
    "dataset_id",
    "genome_build",
    "order",
    "section",
    "key",
    "detail_key",
    "label",
    "category",
    "purpose",
    "criterion",
    "decision",
    "raw_value",
    "qc_passed_value",
    "variants_matching_raw",
    "fraction_matching_raw",
    "variants_unique_to_rule",
    "variants_overlapping_other_rules",
    "provenance_value",
    "raw_variants",
    "qc_passed_variants",
    "excluded_variants",
    "variants_unique_to_one_rule",
    "variants_matching_multiple_rules",
    "retained_fraction",
    "accounting_balanced",
    "requested_threads",
    "polars_thread_pool_size",
    "thread_budget_enforced",
    "raw_vcf",
    "metric_report",
    "rule_report",
    "assessment_json",
    "summary_csv",
    "html_report",
)

def _qc_label_width(indent: int) -> int:
    """Pad labels so every QC field colon has one terminal column."""
    configured = _QC_TERMINAL_LABEL_WIDTH.get()
    if configured is None:
        return _QC_DEFAULT_COLON_INDEX - int(indent) - 4
    return int(configured) + _QC_BASE_INDENT - int(indent)


def _provenance_label(key: str) -> str:
    """Return a readable label for one reported scientific provenance field."""
    return key.replace("_", " ").title()


def metric_available(value: Any) -> bool:
    if value is None:
        return False
    try:
        return bool(value == value)
    except (TypeError, ValueError):
        return False


def format_metric_count(value: Any) -> str:
    if not metric_available(value):
        return "unavailable"
    try:
        return "{:,}".format(int(value))
    except (TypeError, ValueError):
        return str(value)


def _ratio(value: Any) -> str:
    if not metric_available(value):
        return "unavailable"
    return "%.2f" % float(value)


def _number(value: Any) -> str:
    if not metric_available(value):
        return "unavailable"
    return "{:,.2f}".format(float(value))


def format_metric_percent(
    numerator: Any,
    denominator: Any,
    decimals: int = 2,
) -> str:
    """Format a percentage without hiding an unavailable or zero denominator."""
    if not metric_available(numerator) or not metric_available(denominator):
        return "unavailable"
    try:
        denominator_value = float(denominator)
        if denominator_value <= 0:
            return "unavailable"
        value = 100.0 * float(numerator) / denominator_value
    except (TypeError, ValueError, ZeroDivisionError):
        return "unavailable"
    return ("%%.%df%%%%" % int(decimals)) % value


def metric_warning_kind(value: Any) -> str:
    try:
        return "warning" if int(value) else "success"
    except (TypeError, ValueError):
        return "warning"


def _rule_by_key(
    assessment: Mapping[str, Any],
    key: str,
) -> Mapping[str, Any]:
    return next(
        (rule for rule in assessment.get("rules") or () if rule.get("key") == key),
        {},
    )


def _rule_detail(rule: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    return next(
        (detail for detail in rule.get("details") or () if detail.get("key") == key),
        {},
    )


def _field(kind: str, label: str, value: Any, indent: int) -> list[str]:
    return screen_field(
        kind, label, value, indent=indent, label_width=_qc_label_width(indent),
    ).split("\n")


def _render_qc_summary_lines(
    assessment: dict[str, Any],
    *,
    pre_vcf: dict[str, Any] | None = None,
    title: str = "GWAS-VCF quality-control summary",
) -> list[str]:
    """Render the shared raw-rule and combined final QC assessment."""
    include_pre_vcf = pre_vcf is not None
    pre_vcf = pre_vcf or {}
    raw = assessment["raw"]
    passed = assessment["qc_passed"]
    cutoff = assessment["af_difference_cutoff"]
    fields = assessment["field_labels"]
    format_af = fields["study_format_af"]
    format_info = fields["imputation_format"]
    format_neff = fields["effective_sample_size_format"]
    study_info_af = fields["study_info_af"]
    external_info_af = fields["external_info_af"]
    reports = assessment.get("reports") or {}
    header_validation = assessment.get("vcf_header_validation") or {}
    provenance = assessment.get("vcf_provenance") or {}
    provenance_values = provenance.get("values") or {}
    aggregation = assessment.get("aggregation") or {}
    lines = [
        screen_line("analysis", title, indent=4),
        screen_field(
            "info", "Assessment model",
            "the input GWAS-VCF is extracted once; every active configured QC "
            "condition is evaluated independently against the same records",
            indent=6, label_width=_qc_label_width(6),
        ),
        screen_field(
            "info", "VCF output",
            "the input GWAS-VCF is unchanged; no QC-filtered VCF is created",
            indent=6, label_width=_qc_label_width(6),
        ),
    ]
    if include_pre_vcf:
        lines.extend([
            "",
            screen_line("genetic", "1. Before VCF creation", indent=6),
        ])
    for label, key, kind in (
        ("Input rows", "total_variant_infile", "count"),
        ("Rows read", "total_variant_read", "count"),
        (
            "Missing read-stage mandatory values",
            "total_variant_removed_missing_values",
            "loss",
        ),
        ("Duplicate rows", "total_variant_removed_duplicates", "loss"),
        ("Invalid coordinates", "total_variant_removed_null_coords", "loss"),
        ("Non-standard alleles", "total_variant_removed_non_standard_alleles", "loss"),
        ("Passed initial input validation", "total_variant_remaining_for_harmonisation", "success"),
        ("Ready SNPs", "total_variant_ready_snps", "genetic"),
        ("Ready indels / other variants", "total_variant_ready_indels_or_other", "genetic"),
        (
            "Palindromic orientation unavailable",
            "total_variant_removed_palindromic_orientation_unavailable",
            "loss",
        ),
        (
            "Palindromic frequency discordant",
            "total_variant_removed_palindromic_frequency_discordant",
            "loss",
        ),
        ("Palindromic ambiguous removed", "total_variant_removed_palindromic_ambiguous", "loss"),
        ("Reference unmatched removed", "total_variant_removed_reference_unmatched", "loss"),
        ("Reference ambiguous removed", "total_variant_removed_reference_ambiguous", "loss"),
        ("Missing effect frequency", "total_variant_with_missing_eaf", "warning"),
        ("Invalid effect statistics", "total_variant_with_invalid_beta_se", "loss"),
        ("All chromosome-stage removals", "total_variant_removed_chromosome_harmonisation", "loss"),
        (
            "Passed chromosome harmonisation",
            "total_variant_passed_chromosome_harmonisation",
            "success",
        ),
        ("Sent to GWAS-to-VCF", "total_variant_in_vcf_input", "count"),
    ):
        if key in pre_vcf:
            lines.extend(_field(
                kind, label, format_metric_count(pre_vcf.get(key)), 8,
            ))

    initial = pre_vcf.get("total_variant_remaining_for_harmonisation")
    removed = pre_vcf.get("total_variant_removed_chromosome_harmonisation")
    passed_chromosomes = pre_vcf.get(
        "total_variant_passed_chromosome_harmonisation"
    )
    sent = pre_vcf.get("total_variant_in_vcf_input")
    if all(
        metric_available(value)
        for value in (initial, removed, passed_chromosomes, sent)
    ):
        balanced = (
            int(initial) - int(removed) == int(passed_chromosomes)
            and int(passed_chromosomes) == int(sent)
        )
        lines.extend(_field(
            "success" if balanced else "loss",
            "Pre-VCF count reconciliation",
            "%s initial = %s sent to GWAS-to-VCF + %s removed%s"
            % (
                format_metric_count(initial), format_metric_count(sent),
                format_metric_count(removed),
                "" if balanced else " — DOES NOT BALANCE",
            ),
            8,
        ))

    validation_section_number = 2 if include_pre_vcf else 1
    profile_section_number = validation_section_number + 1
    lines.extend([
        "",
        screen_line(
            "genetic",
            "%d. %s" % (
                validation_section_number,
                (
                    "Validated raw merged GWAS-VCF"
                    if include_pre_vcf
                    else "Validated input summary-statistics GWAS-VCF"
                ),
            ),
            indent=6,
        ),
    ])
    lines.extend(_field(
        "info", "Input VCF", os.path.basename(assessment["raw_vcf"]), 8,
    ))
    if header_validation:
        lines.extend(_field(
            "genetic",
            "Genome build",
            header_validation.get("genome_build") or "unavailable",
            8,
        ))
        lines.extend(_field(
            "count",
            "Declared contigs",
            format_metric_count(header_validation.get("declared_contig_count")),
            8,
        ))
    lines.extend(_field(
        "count", "Total variants", format_metric_count(raw["num_records"]), 8,
    ))
    if header_validation:
        required_fields = list(header_validation.get("required_fields") or ())
        required_count = header_validation.get("required_field_count")
        lines.extend(_field(
            (
                "success"
                if header_validation.get("status") == "passed"
                else "loss"
            ),
            "VCF header field contract",
            "%s · %s required INFO/FORMAT fields declared"
            % (
                str(header_validation.get("status", "unavailable")).upper(),
                format_metric_count(required_count),
            ),
            8,
        ))
        lines.extend(_field(
            "info",
            "Required declared fields",
            ", ".join(required_fields) if required_fields else "unavailable",
            8,
        ))
    if provenance:
        provenance_status = provenance.get("status", "not_checked")
        lines.extend(_field(
            "success" if provenance_status == "available" else "info",
            "PostGWAS scientific provenance",
            provenance_status.replace("_", " "),
            8,
        ))
        for label, key in (
            ("Allele-frequency source", "af_source"),
            ("Allele-frequency meaning", "af_input_type"),
            ("Imputation quality score source", "info_source"),
            ("Imputation quality score interpretation", "info_interpretation"),
            ("Imputation quality score output field", "info_output_field"),
        ):
            if provenance_values.get(key):
                lines.extend(_field("info", label, provenance_values[key], 8))
    lines.extend([
        "",
        screen_line(
            "analysis",
            "%d. %s" % (
                profile_section_number,
                (
                    "Profile raw merged GWAS-VCF"
                    if include_pre_vcf
                    else "Profile input GWAS-VCF"
                ),
            ),
            indent=6,
        ),
    ])
    if aggregation:
        lines.extend(_field(
            "info",
            "Aggregation",
            "%s temporary-table scans using %s; Polars %s"
            % (
                format_metric_count(aggregation.get("temporary_table_scans")),
                aggregation.get("streaming_collection_api", "unavailable"),
                aggregation.get("polars_version", "unavailable"),
            ),
            8,
        ))
    lines.extend(_field("count", "SNPs", format_metric_count(raw["num_snps"]), 8))
    lines.extend(_field(
        "count", "Indels / other variants",
        format_metric_count(raw["num_non_snps"]), 8,
    ))
    lines.extend(_field("analysis", "Transition / transversion", _ratio(raw["ts_tv_ratio"]), 8))
    lines.append("")
    lines.append(screen_line("analysis", "Effective sample-size distribution", indent=10))
    lines.extend(_field(
        "count", "Usable %s" % format_neff,
        format_metric_count(raw["effective_sample_size_available"]), 12,
    ))
    lines.extend(_field(
        metric_warning_kind(raw["effective_sample_size_missing_or_invalid"]),
        "Missing or invalid %s" % format_neff,
        format_metric_count(raw["effective_sample_size_missing_or_invalid"]), 12,
    ))
    for label, key in (
        ("Minimum", "effective_sample_size_minimum"),
        ("Maximum", "effective_sample_size_maximum"),
        ("Mean", "effective_sample_size_mean"),
        ("Sample standard deviation", "effective_sample_size_standard_deviation"),
        (
            "Upper outlier threshold (mean + %g SD)"
            % assessment["sample_size_outlier_standard_deviations"],
            "effective_sample_size_outlier_threshold",
        ),
    ):
        lines.extend(_field("analysis", label, _number(raw[key]), 12))
    lines.extend(_field(
        "analysis",
        "Upper-tail diagnostic count",
        format_metric_count(raw["effective_sample_size_above_outlier_threshold"]), 12,
    ))
    lines.extend(_field(
        "analysis",
        "Low-Neff reference (raw q=%.3f)"
        % assessment["sample_size_reference_quantile"],
        _number(raw["effective_sample_size_reference_quantile_value"]),
        12,
    ))
    lines.extend(_field(
        "analysis",
        "Low-Neff threshold (%.2f%% of reference)"
        % (100.0 * assessment["sample_size_minimum_fraction_of_reference"]),
        _number(raw["effective_sample_size_minimum_threshold"]),
        12,
    ))
    lines.extend(_field(
        metric_warning_kind(
            raw["effective_sample_size_below_minimum_threshold"]
        ),
        "Below low-Neff threshold",
        "%s (%s of usable values)"
        % (
            format_metric_count(
                raw["effective_sample_size_below_minimum_threshold"]
            ),
            format_metric_percent(
                raw["effective_sample_size_below_minimum_threshold"],
                raw["effective_sample_size_available"],
            ),
        ),
        12,
    ))
    lines.append("")
    lines.append(screen_line("genetic", "Frequency and imputation fields", indent=10))
    for label, key in (
        ("Missing %s" % format_af, "format_af_missing"),
        ("Missing %s" % format_info, "format_si_missing"),
        ("Missing study AF (%s)" % study_info_af, "study_af_missing"),
        ("Missing %s" % external_info_af, "external_af_missing"),
        ("Study/reference AF pairs", "af_comparable"),
        ("|%s − %s| > %s" % (study_info_af, external_info_af, cutoff),
         "af_difference_above_cutoff"),
    ):
        value = raw[key]
        kind = (
            metric_warning_kind(value)
            if "Missing" in label or key == "af_difference_above_cutoff"
            else "count"
        )
        lines.extend(_field(kind, label, format_metric_count(value), 12))

    lines.extend([
        "",
        screen_line(
            "analysis", "%d. QC conditions assessed on the input GWAS-VCF"
            % (validation_section_number + 2),
            indent=6,
        ),
        screen_field(
            "info", "Assessment basis",
            "all %s input variants are tested against every active condition. Rule "
            "counts may overlap and must not be added together"
            % format_metric_count(raw["num_records"]),
            indent=8, label_width=_qc_label_width(8),
        ),
    ])
    for rule in assessment["rules"]:
        lines.extend([
            "",
            screen_line(
                "analysis",
                "Rule %d · %s" % (rule["number"], rule["label"]),
                indent=10,
            ),
        ])
        lines.extend(_field("info", "Condition", rule["criterion"], 12))
        lines.extend(_field("info", "Category", rule["category"], 12))
        lines.extend(_field("info", "Purpose", rule["purpose"], 12))
        lines.extend(_field(
            "loss" if rule["failed_raw"] else "success",
            "Would be excluded",
            "%s (%s of input)"
            % (
                format_metric_count(rule["failed_raw"]),
                format_metric_percent(rule["failed_raw"], raw["num_records"]),
            ),
            12,
        ))
        lines.extend(_field(
            "analysis", "Unique to this rule",
            format_metric_count(rule["unique_only_raw"]), 12,
        ))
        lines.extend(_field(
            "analysis", "Also matched another rule",
            format_metric_count(rule["overlap_raw"]), 12,
        ))
        for detail in rule["details"]:
            decision = "excluded from the virtual subset"
            if detail["decision"] == "retain_by_configured_policy":
                decision = "retained by configured policy"
            lines.extend(_field(
                (
                    "loss"
                    if detail["matched_raw"]
                    and detail["decision"] == "exclude_from_virtual_subset"
                    else "info"
                ),
                detail["label"],
                "%s (%s of input); %s"
                % (
                    format_metric_count(detail["matched_raw"]),
                    format_metric_percent(
                        detail["matched_raw"], raw["num_records"],
                    ),
                    decision,
                ),
                12,
            ))

    if assessment.get("inactive_rules"):
        lines.extend([
            "",
            screen_line("info", "Inactive QC conditions", indent=10),
        ])
        for inactive in assessment["inactive_rules"]:
            lines.extend(_field(
                "info", inactive["label"], inactive["reason"], 12,
            ))

    retained_percent = 100.0 * assessment["retained_fraction"]
    lines.extend([
        "",
        screen_line(
            "success", "%d. Final virtual QC assessment"
            % (validation_section_number + 3), indent=6,
        ),
        screen_field(
            "info", "Meaning of QC-passed", assessment["definition"],
            indent=8, label_width=_qc_label_width(8),
        ),
        screen_field(
            "decision", "Combined decision",
            "all %s active conditions are applied together to every raw record; "
            "a variant is retained only if it passes every condition"
            % format_metric_count(assessment["active_rule_count"]),
            indent=8, label_width=_qc_label_width(8),
        ),
    ])
    lines.extend(_field(
        "count", "Input GWAS-VCF variants", format_metric_count(raw["num_records"]), 8,
    ))
    lines.extend(_field(
        "success", "QC-passed variants",
        "%s (%.2f%% of input)"
        % (format_metric_count(passed["num_records"]), retained_percent),
        8,
    ))
    lines.extend(_field(
        "loss", "Excluded by combined virtual policy",
        format_metric_count(assessment["excluded_total"]),
        8,
    ))
    lines.extend(_field(
        "analysis", "Unique to one condition",
        format_metric_count(assessment["unique_rule_only_total"]),
        8,
    ))
    lines.extend(_field(
        "analysis", "Matched multiple conditions",
        format_metric_count(assessment["overlap_variants"]),
        8,
    ))
    lines.extend(_field(
        "success", "Count reconciliation",
        "%s input = %s virtual QC-passed + %s excluded" % (
            format_metric_count(raw["num_records"]),
            format_metric_count(passed["num_records"]),
            format_metric_count(assessment["excluded_total"]),
        ),
        8,
    ))
    lines.extend([
        "",
        screen_line("genetic", "Final virtual subset metrics", indent=10),
    ])
    for label, key, kind, formatter in (
        ("SNPs", "num_snps", "count", format_metric_count),
        ("Indels / other variants", "num_non_snps", "count", format_metric_count),
        ("Transition / transversion", "ts_tv_ratio", "analysis", _ratio),
        (
            "Missing %s" % format_af, "format_af_missing", "warning",
            format_metric_count,
        ),
        (
            "Missing %s" % format_info, "format_si_missing", "warning",
            format_metric_count,
        ),
        (
            "Missing %s" % external_info_af, "external_af_missing", "warning",
            format_metric_count,
        ),
        (
            "AF differences above cutoff", "af_difference_above_cutoff",
            "warning", format_metric_count,
        ),
        (
            "Usable %s" % format_neff, "effective_sample_size_available", "count",
            format_metric_count,
        ),
        (
            "Missing or invalid %s" % format_neff,
            "effective_sample_size_missing_or_invalid", "warning",
            format_metric_count,
        ),
        ("Minimum %s" % format_neff, "effective_sample_size_minimum", "analysis", _number),
        ("Maximum %s" % format_neff, "effective_sample_size_maximum", "analysis", _number),
        ("Mean %s" % format_neff, "effective_sample_size_mean", "analysis", _number),
        (
            "Sample SD %s" % format_neff,
            "effective_sample_size_standard_deviation", "analysis", _number,
        ),
        (
            "%s upper outlier threshold" % format_neff,
            "effective_sample_size_outlier_threshold", "analysis", _number,
        ),
        (
            "%s upper-tail diagnostic count" % format_neff,
            "effective_sample_size_above_outlier_threshold", "analysis",
            format_metric_count,
        ),
        (
            "%s raw q=%.3f reference"
            % (format_neff, assessment["sample_size_reference_quantile"]),
            "effective_sample_size_reference_quantile_value", "analysis",
            _number,
        ),
        (
            "%s low-Neff threshold" % format_neff,
            "effective_sample_size_minimum_threshold", "analysis", _number,
        ),
        (
            "%s values below low-Neff threshold" % format_neff,
            "effective_sample_size_below_minimum_threshold", "warning",
            lambda value: "%s (%s of usable values)" % (
                format_metric_count(value),
                format_metric_percent(
                    value,
                    passed["effective_sample_size_available"],
                ),
            ),
        ),
    ):
        resolved_kind = (
            metric_warning_kind(passed[key]) if kind == "warning" else kind
        )
        lines.extend(_field(resolved_kind, label, formatter(passed[key]), 12))
    lines.extend(_field(
        "info", "QC-passed VCF", "not created; use the unchanged input GWAS-VCF",
        8,
    ))
    if reports.get("summary"):
        lines.extend(_field(
            "info", "Metric report", os.path.basename(reports["summary"]), 8,
        ))
    if reports.get("rules"):
        lines.extend(_field(
            "info", "Rule report", os.path.basename(reports["rules"]), 8,
        ))
    if reports.get("json"):
        lines.extend(_field(
            "info", "Assessment JSON", os.path.basename(reports["json"]), 8,
        ))
    if reports.get("csv"):
        lines.extend(_field(
            "success", "Summary CSV", os.path.basename(reports["csv"]), 8,
        ))
    if reports.get("html"):
        lines.extend(_field(
            "success", "Detailed HTML report",
            os.path.basename(reports["html"]), 8,
        ))
    return lines


def qc_summary_lines(
    assessment: dict[str, Any],
    *,
    pre_vcf: dict[str, Any] | None = None,
    title: str = "GWAS-VCF quality-control summary",
    label_width: int | None = None,
) -> list[str]:
    """Render QC fields in one configured value column at every hierarchy."""
    if label_width is not None and int(label_width) < 1:
        raise ValueError("label_width must be a positive integer")
    token = _QC_TERMINAL_LABEL_WIDTH.set(
        None if label_width is None else int(label_width)
    )
    try:
        return _render_qc_summary_lines(
            assessment, pre_vcf=pre_vcf, title=title,
        )
    finally:
        _QC_TERMINAL_LABEL_WIDTH.reset(token)


def _summary_card_options(label_width: int | None) -> dict[str, int]:
    if label_width is not None and int(label_width) < 1:
        raise ValueError("label_width must be a positive integer")
    return {} if label_width is None else {"label_width": int(label_width)}


def qc_input_validation_screen_lines(
    evidence: Mapping[str, Any],
    *,
    label_width: int | None = None,
    include_total: bool = True,
) -> list[str]:
    """Render validation using the filtering module's screen hierarchy."""
    card_options = _summary_card_options(label_width)
    raw = evidence.get("raw") or {}
    fields = evidence.get("field_labels") or {}
    header = evidence.get("vcf_header_validation") or {}
    provenance = evidence.get("vcf_provenance") or {}
    provenance_values = provenance.get("values") or {}
    header_passed = header.get("status") == "passed"
    provenance_available = provenance.get("status") == "available"
    score_source = provenance_values.get("info_source") or (
        "source not declared in VCF provenance"
    )
    card_fields = [
        (
            "info",
            "Input VCF",
            os.path.basename(str(evidence.get("raw_vcf") or "unavailable")),
        ),
        (
            "genetic",
            "Genome build",
            header.get("genome_build")
            or evidence.get("genome_build")
            or "unavailable",
        ),
        (
            "count",
            "Declared contigs",
            format_metric_count(header.get("declared_contig_count")),
        ),
    ]
    if include_total:
        card_fields.append((
            "count", "Total variants", format_metric_count(raw.get("num_records")),
        ))
    card_fields.extend([
        (
            "success" if header_passed else "warning",
            "Header contract",
            "%s · %s/%s required INFO/FORMAT fields declared"
            % (
                str(header.get("status") or "unavailable").upper(),
                format_metric_count(header.get("required_field_count")),
                format_metric_count(header.get("required_field_count")),
            ),
        ),
        (
            "success" if provenance_available else "info",
            "Field provenance",
            str(provenance.get("status") or "not checked").replace(
                "_", " ",
            ).upper(),
        ),
        (
            "genetic",
            "Allele-frequency fields",
            "%s for MAF; %s versus %s"
            % (
                fields.get("study_format_af", "unavailable"),
                fields.get("study_info_af", "unavailable"),
                fields.get("external_info_af", "unavailable"),
            ),
        ),
        (
            "analysis",
            "Imputation quality score",
            "%s · %s"
            % (fields.get("imputation_format", "unavailable"), score_source),
        ),
    ])
    card_fields = consolidate_validation_fields(card_fields)
    flush_file_validation_display()
    if not card_fields:
        return []
    lines = screen_summary_card(
        "success" if header_passed else "warning",
        "Input VCF validation",
        card_fields,
        title_indent=4,
        field_indent=8,
        **card_options,
    )
    lines.insert(2, "")
    lines.append("")
    return lines


def qc_rules_applied_screen_lines(
    rule_contract: Mapping[str, Any],
    *,
    label_width: int | None = None,
) -> list[str]:
    """Render the virtual-subset plan with the shared action-plan formatter."""
    _summary_card_options(label_width)
    decision_plan = list(rule_contract.get("decision_plan") or ())
    if not decision_plan:
        decision_plan = [
            {
                "action": (
                    "EXCLUDE"
                    if detail.get("decision") == "exclude_from_virtual_subset"
                    else "KEEP"
                ),
                "label": detail.get("label") or "Unnamed QC condition",
            }
            for rule in rule_contract.get("rules") or ()
            for detail in rule.get("details") or ()
        ]

    inactive_labels = [
        str(rule.get("label") or rule.get("key") or "unknown rule")
        for rule in rule_contract.get("inactive_rules") or ()
    ]
    return render_action_plan(
        "QC assessment plan · conditions evaluated independently",
        decision_plan,
        empty_message="No QC exclusion condition is active",
        footer=(
            (
                "info",
                "EXCLUDE affects only the reported virtual QC-passed subset; "
                "the input VCF remains unchanged",
            ),
            (
                "info",
                "Inactive optional rules · %s"
                % (", ".join(inactive_labels) if inactive_labels else "none"),
            ),
        ),
    ).split("\n")


def qc_rule_results_screen_lines(
    assessment: Mapping[str, Any],
    *,
    label_width: int | None = None,
) -> list[str]:
    """Render independent QC outcomes using filtering-style rule groups."""
    card_options = _summary_card_options(label_width)
    resolved_label_width = card_options.get(
        "label_width", SUMMARY_CARD_LABEL_WIDTH,
    )
    separator_column = (
        8 + cell_len(SYMBOLS["count"]) + 2 + resolved_label_width
    )
    raw_total = (assessment.get("raw") or {}).get("num_records")
    rules = list(assessment.get("rules") or ())
    lines = screen_section_lines("analysis", "QC results by rule")
    current_group = None
    for index, rule in enumerate(rules):
        group = str(
            rule.get("display_group")
            or rule.get("category")
            or "QC rule"
        )
        group_kind = str(rule.get("display_group_kind") or "analysis")
        if group_kind not in SYMBOLS:
            group_kind = "analysis"
        if group != current_group:
            lines.extend((screen_line(group_kind, group, indent=8), ""))
            current_group = group

        failed = int(rule.get("failed_raw") or 0)
        variant_word = "variant" if failed == 1 else "variants"
        lines.extend((
            screen_line(
                "loss" if failed else "success",
                str(rule.get("label") or rule.get("key") or "QC rule"),
                indent=12,
            ),
            screen_line(
                "info",
                "%s %s failed this rule (%s of input):"
                % (
                    format_metric_count(failed),
                    variant_word,
                    format_metric_percent(failed, raw_total),
                ),
                indent=16,
            ),
            "%s• %s failed only this QC rule"
            % (" " * 20, format_metric_count(rule.get("unique_only_raw"))),
            "%s• %s also failed one or more other QC rules"
            % (" " * 20, format_metric_count(rule.get("overlap_raw"))),
        ))
        details = list(rule.get("details") or ())
        if details:
            lines.append(
                "%s• Trigger counts: %s"
                % (
                    " " * 20,
                    " · ".join(
                        "%s (%s)"
                        % (
                            detail.get("label")
                            or detail.get("key")
                            or "Trigger",
                            format_metric_count(detail.get("matched_raw")),
                        )
                        for detail in details
                    ),
                )
            )
        if index < len(rules) - 1:
            lines.append("")

    inactive = list(assessment.get("inactive_rules") or ())
    if inactive:
        lines.extend(("", screen_line("info", "Inactive QC rules", indent=8)))
        lines.extend(
            screen_field(
                "info",
                str(rule.get("label") or rule.get("key") or "QC rule"),
                rule.get("reason") or "inactive by configuration",
                indent=12,
                separator_column=separator_column,
            )
            for rule in inactive
        )
    return lines


def qc_final_screen_lines(
    assessment: Mapping[str, Any],
    *,
    label_width: int | None = None,
) -> list[str]:
    """Render the filtering-style final QC summary after report publication."""
    card_options = _summary_card_options(label_width)
    resolved_label_width = card_options.get(
        "label_width", SUMMARY_CARD_LABEL_WIDTH,
    )
    raw = assessment["raw"]
    passed = assessment["qc_passed"]
    reports = assessment.get("reports") or {}
    raw_total = raw.get("num_records")
    passed_total = passed.get("num_records")
    excluded_total = assessment.get("excluded_total")

    low_neff = raw.get("effective_sample_size_below_minimum_threshold")
    missing_neff = raw.get("effective_sample_size_missing_or_invalid")
    lines = [
        "",
        screen_line("analysis", "GWAS-VCF quality-control summary", indent=4),
        "",
        screen_field(
            "info", "Dataset", assessment.get("dataset_id") or "unavailable",
            indent=8, label_width=resolved_label_width,
        ),
        screen_field(
            "genetic", "Genome build",
            assessment.get("genome_build") or "unavailable",
            indent=8, label_width=resolved_label_width,
        ),
        screen_field(
            "info", "Assessment method",
            "Every active QC rule is evaluated independently; a variant enters "
            "the virtual QC-passed subset only when no active rule excludes it",
            indent=8, label_width=resolved_label_width,
        ),
    ]
    lines.extend(qc_rule_results_screen_lines(
        assessment, label_width=label_width,
    ))

    lines.extend(screen_section_lines("count", "Exact virtual-subset accounting"))
    lines.extend((
        screen_field(
            "analysis", "Variants failing two or more QC rules",
            format_metric_count(assessment.get("overlap_variants")),
            indent=8, label_width=resolved_label_width,
        ),
        screen_field(
            "loss" if excluded_total else "success",
            "Excluded by combined QC policy", format_metric_count(excluded_total),
            indent=8, label_width=resolved_label_width,
        ),
        screen_field(
            "success" if assessment.get("accounting_balanced") else "error",
            "Virtual-subset count check",
            "%s passed + %s excluded %s %s input"
            % (
                format_metric_count(passed_total),
                format_metric_count(excluded_total),
                "=" if assessment.get("accounting_balanced") else "!=",
                format_metric_count(raw_total),
            ),
            indent=8, label_width=resolved_label_width,
        ),
    ))

    lines.extend(screen_section_lines("analysis", "Saved reports"))
    for kind, label, key in (
        ("info", "Metric report", "summary"),
        ("info", "Rule report", "rules"),
        ("info", "Assessment JSON", "json"),
        ("success", "Summary CSV", "csv"),
        ("success", "Detailed HTML report", "html"),
    ):
        value = reports.get(key)
        if value:
            lines.append(screen_field(
                kind, label, os.path.basename(str(value)),
                indent=8, label_width=resolved_label_width,
            ))

    retained_fraction = assessment.get("retained_fraction")
    retained_percent = (
        "unavailable"
        if retained_fraction is None
        else "%.2f%%" % (100.0 * float(retained_fraction))
    )
    lines.extend(screen_section_lines("count", "Final virtual QC outcome"))
    lines.extend((
        screen_field(
            "count", "Input GWAS-VCF", format_metric_count(raw_total),
            indent=8, label_width=resolved_label_width,
        ),
        screen_field(
            "count", "Virtual QC-passed subset", format_metric_count(passed_total),
            indent=8, label_width=resolved_label_width,
        ),
        screen_field(
            "success", "Variants retained", retained_percent,
            indent=8, label_width=resolved_label_width,
        ),
        screen_field(
            (
                "warning"
                if int(low_neff or 0) or int(missing_neff or 0)
                else "success"
            ),
            "Effective sample size",
            "%s usable · %s missing/invalid · %s below %s"
            % (
                format_metric_count(raw.get("effective_sample_size_available")),
                format_metric_count(missing_neff),
                format_metric_count(low_neff),
                _number(raw.get("effective_sample_size_minimum_threshold")),
            ),
            indent=8, label_width=resolved_label_width,
        ),
        screen_field(
            "info", "VCF output", "unchanged · no filtered VCF created",
            indent=8, label_width=resolved_label_width,
        ),
        "",
    ))
    return lines


def qc_summary_screen_lines(
    assessment: Mapping[str, Any],
    *,
    label_width: int | None = None,
) -> list[str]:
    """Render validation, plan, and final summary when stages are not visible."""
    _summary_card_options(label_width)
    lines = qc_input_validation_screen_lines(
        assessment, label_width=label_width,
    )
    lines.extend(qc_rules_applied_screen_lines(
        assessment, label_width=label_width,
    ))
    lines.extend(qc_final_screen_lines(
        assessment, label_width=label_width,
    ))
    return lines


def qc_summary_csv_records(
    assessment: Mapping[str, Any],
) -> list[dict[str, Any]]:
    """Return one overall row plus validation, metric, rule, and detail rows."""
    raw = assessment["raw"]
    passed = assessment["qc_passed"]
    aggregation = assessment.get("aggregation") or {}
    reports = assessment.get("reports") or {}
    common = {
        "dataset_id": assessment.get("dataset_id"),
        "genome_build": assessment.get("genome_build"),
    }
    records = [{
        **common,
        "record_type": "overall",
        "section": "assessment",
        "raw_variants": raw.get("num_records"),
        "qc_passed_variants": passed.get("num_records"),
        "excluded_variants": assessment.get("excluded_total"),
        "variants_unique_to_one_rule": assessment.get("unique_rule_only_total"),
        "variants_matching_multiple_rules": assessment.get("overlap_variants"),
        "retained_fraction": assessment.get("retained_fraction"),
        "accounting_balanced": assessment.get("accounting_balanced"),
        "requested_threads": aggregation.get("requested_threads"),
        "polars_thread_pool_size": aggregation.get("polars_thread_pool_size"),
        "thread_budget_enforced": aggregation.get("thread_budget_enforced"),
        "raw_vcf": assessment.get("raw_vcf"),
        "metric_report": reports.get("summary"),
        "rule_report": reports.get("rules"),
        "assessment_json": reports.get("json"),
        "summary_csv": reports.get("csv"),
        "html_report": reports.get("html"),
    }]
    for order, (key, value) in enumerate(
        (assessment.get("vcf_header_validation") or {}).items(), 1,
    ):
        records.append({
            **common,
            "record_type": "validation",
            "order": order,
            "section": "vcf_header_validation",
            "key": key,
            "label": key.replace("_", " ").title(),
            "raw_value": (
                ",".join(str(item) for item in value)
                if isinstance(value, (list, tuple)) else value
            ),
        })
    metric_keys = list(raw)
    metric_keys.extend(key for key in passed if key not in raw)
    for order, key in enumerate(metric_keys, 1):
        records.append({
            **common,
            "record_type": "metric",
            "order": order,
            "section": "raw_vs_qc_passed",
            "key": key,
            "label": key.replace("_", " ").title(),
            "raw_value": raw.get(key),
            "qc_passed_value": passed.get(key),
        })
    for order, rule in enumerate(assessment.get("rules") or (), 1):
        records.append({
            **common,
            "record_type": "rule",
            "order": order,
            "section": "active_qc_rules",
            "key": rule.get("key"),
            "label": rule.get("label"),
            "category": rule.get("category"),
            "purpose": rule.get("purpose"),
            "criterion": rule.get("criterion"),
            "decision": rule.get("decision"),
            "variants_matching_raw": rule.get("failed_raw"),
            "fraction_matching_raw": rule.get("failed_fraction_raw"),
            "variants_unique_to_rule": rule.get("unique_only_raw"),
            "variants_overlapping_other_rules": rule.get("overlap_raw"),
        })
        for detail_order, detail in enumerate(rule.get("details") or (), 1):
            records.append({
                **common,
                "record_type": "rule_detail",
                "order": "%s.%s" % (order, detail_order),
                "section": "active_qc_rules",
                "key": rule.get("key"),
                "detail_key": detail.get("key"),
                "label": detail.get("label"),
                "category": rule.get("category"),
                "purpose": rule.get("purpose"),
                "criterion": rule.get("criterion"),
                "decision": detail.get("decision"),
                "variants_matching_raw": detail.get("matched_raw"),
                "fraction_matching_raw": detail.get("matched_fraction_raw"),
            })
    for order, inactive in enumerate(assessment.get("inactive_rules") or (), 1):
        records.append({
            **common,
            "record_type": "inactive_rule",
            "order": order,
            "section": "inactive_qc_rules",
            "key": inactive.get("key"),
            "label": inactive.get("label"),
            "category": inactive.get("category"),
            "purpose": inactive.get("reason"),
            "decision": "inactive",
        })
    provenance = assessment.get("vcf_provenance") or {}
    records.append({
        **common,
        "record_type": "vcf_provenance",
        "order": 1,
        "section": "vcf_scientific_provenance",
        "key": "status",
        "label": "Scientific provenance status",
        "provenance_value": provenance.get("status"),
    })
    for order, (key, value) in enumerate(
        (provenance.get("values") or {}).items(), 2,
    ):
        records.append({
            **common,
            "record_type": "vcf_provenance",
            "order": order,
            "section": "vcf_scientific_provenance",
            "key": key,
            "label": _provenance_label(key),
            "provenance_value": value,
        })
    for order, (key, value) in enumerate(aggregation.items(), 1):
        records.append({
            **common,
            "record_type": "provenance",
            "order": order,
            "section": "aggregation",
            "key": key,
            "label": key.replace("_", " ").title(),
            "raw_value": value,
        })
    return records


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


def _html_cell(value: Any) -> str:
    return escape(_html_value(value))


def _html_path(value: Any) -> str:
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


def _html_table(
    headers: Sequence[str], rows: Sequence[Sequence[Any]],
) -> str:
    heading = "".join(
        '<th scope="col">%s</th>' % escape(label) for label in headers
    )
    body = "".join(
        "<tr>%s</tr>" % "".join(
            "<td>%s</td>" % _html_cell(value) for value in row
        )
        for row in rows
    )
    if not rows:
        body = '<tr><td colspan="%d" class="muted">No rows</td></tr>' % len(
            headers
        )
    return (
        '<div class="table-wrap"><table><thead><tr>%s</tr></thead>'
        '<tbody>%s</tbody></table></div>' % (heading, body)
    )


def _html_key_values(
    rows: Sequence[tuple[str, Any]], *, path_labels: Sequence[str] = (),
) -> str:
    paths = set(path_labels)
    body = "".join(
        '<tr><th scope="row">%s</th><td>%s</td></tr>'
        % (
            escape(label),
            _html_path(value) if label in paths else _html_cell(value),
        )
        for label, value in rows
    )
    return (
        '<div class="table-wrap"><table class="key-value"><tbody>%s</tbody>'
        '</table></div>' % body
    )


_QC_HTML_STYLES = """
:root{--ink:#172033;--muted:#64748b;--line:#dbe3ee;--canvas:#f4f7fb;--panel:#fff;--brand:#155e75;--brand2:#0891b2;--ok:#166534;--warn:#9a3412;--bad:#991b1b}*{box-sizing:border-box}body{margin:0;background:var(--canvas);color:var(--ink);font:15px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}.page-header{padding:38px max(24px,calc((100vw - 1160px)/2));color:#fff;background:linear-gradient(135deg,#164e63,#0e7490)}.eyebrow{margin:0 0 5px;font-size:12px;font-weight:800;letter-spacing:.14em;text-transform:uppercase;opacity:.86}.page-header h1{margin:0;font-size:clamp(28px,4vw,44px);line-height:1.1}.page-header p{max-width:900px;margin:12px 0 0;color:#cffafe}.container{max-width:1160px;margin:0 auto;padding:26px 24px 58px}.notice,section{background:var(--panel);border:1px solid var(--line);border-radius:13px;box-shadow:0 7px 24px rgba(15,23,42,.05)}.notice{margin-bottom:20px;padding:14px 16px;border-left:4px solid var(--brand2)}.notice.warning,section.warning{border-left:4px solid var(--warn);background:#fff7ed}.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px;margin:0 0 22px}.metric{padding:15px;border:1px solid var(--line);border-radius:11px;background:#fff}.metric-label{font-size:12px;font-weight:800;color:var(--muted);letter-spacing:.05em;text-transform:uppercase}.metric-value{margin-top:4px;font-size:22px;font-weight:780;overflow-wrap:anywhere}.badge,.category{display:inline-flex;padding:4px 9px;border-radius:999px;font-size:12px;font-weight:800}.badge.ok{color:var(--ok);background:#dcfce7}.badge.fail{color:var(--bad);background:#fee2e2}.category{color:#075985;background:#e0f2fe}section{margin:18px 0;padding:22px}section h2{margin:0 0 5px;font-size:22px}.section-note{margin:0 0 16px;color:var(--muted)}.findings{margin:10px 0 0;padding-left:22px}.findings li{margin:8px 0}.rule-list{display:grid;gap:14px}.rule-card{border:1px solid var(--line);border-radius:11px;padding:18px;background:#fff}.rule-head{display:flex;align-items:flex-start;justify-content:space-between;gap:18px}.rule-head h3{margin:5px 0 0;font-size:19px}.rule-impact{text-align:right;white-space:nowrap}.rule-impact strong{display:block;font-size:21px}.rule-impact span{color:var(--muted);font-size:13px}.rule-purpose{margin:12px 0;color:#334155}.rule-meta{display:grid;grid-template-columns:minmax(150px,1fr) minmax(180px,2fr);gap:5px 16px;margin:12px 0}.rule-meta dt{font-weight:750;color:#475569}.rule-meta dd{margin:0}.impact-bar{height:8px;margin:12px 0;background:#e2e8f0;border-radius:999px;overflow:hidden}.impact-bar span{display:block;height:100%;width:var(--impact);background:linear-gradient(90deg,#0891b2,#155e75)}details{margin-top:12px}summary{cursor:pointer;font-weight:750;color:#075985}.table-wrap{overflow:auto;border:1px solid var(--line);border-radius:9px}table{width:100%;border-collapse:collapse;background:#fff}th,td{padding:10px 12px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}thead th{background:#ecfeff;color:#334155;font-size:12px;text-transform:uppercase;letter-spacing:.04em}.key-value th{width:300px;background:#f8fafc;color:#475569}.path{overflow-wrap:anywhere;word-break:break-word;color:#075985}.muted{color:var(--muted)}.footer{max-width:1160px;margin:0 auto;padding:0 24px 34px;color:var(--muted);font-size:13px}@media(max-width:700px){.container{padding:16px 10px 40px}section{padding:14px}.rule-head{display:block}.rule-impact{text-align:left;margin-top:8px}.rule-meta{grid-template-columns:1fr}.rule-meta dd{margin-bottom:7px}.key-value th{width:42%}th,td{padding:8px}.page-header{padding:28px 18px}}@media print{body{background:#fff}.page-header{background:#fff;color:var(--ink);padding:0}.page-header p{color:var(--muted)}.container{max-width:none;padding:0}.notice,section,.rule-card{box-shadow:none;break-inside:avoid}.table-wrap{overflow:visible}details{display:block}details>*{display:block}a{color:inherit;text-decoration:none}}
"""


def _decision_label(value: Any) -> str:
    return {
        "exclude_from_virtual_subset": "Exclude from virtual subset",
        "retain_by_configured_policy": "Retain by configured policy",
        "inactive": "Inactive",
    }.get(str(value), str(value).replace("_", " ").capitalize())


def _missing_policy_label(value: Any) -> str:
    return {
        "remove": "Exclude when missing",
        "keep": "Retain when missing",
    }.get(str(value), _decision_label(value))


def _html_percent(numerator: Any, denominator: Any) -> str:
    return format_metric_percent(numerator, denominator)


def _html_rule_cards(assessment: Mapping[str, Any]) -> str:
    raw_total = (assessment.get("raw") or {}).get("num_records")
    cards = []
    for rule in assessment.get("rules") or ():
        fraction = 100.0 * float(rule.get("failed_fraction_raw") or 0.0)
        width = min(100.0, max(0.0, fraction))
        detail_rows = [
            (
                detail.get("label"),
                _decision_label(detail.get("decision")),
                "%s (%s)" % (
                    format_metric_count(detail.get("matched_raw")),
                    _html_percent(detail.get("matched_raw"), raw_total),
                ),
            )
            for detail in rule.get("details") or ()
        ]
        cards.append(
            '<article class="rule-card"><div class="rule-head"><div>'
            '<span class="category">%s</span><h3>Rule %s · %s</h3></div>'
            '<div class="rule-impact"><strong>%s</strong><span>%s of input</span>'
            '</div></div><p class="rule-purpose">%s</p>'
            '<div class="impact-bar" aria-hidden="true"><span style="--impact:%.6g%%">'
            '</span></div><dl class="rule-meta"><dt>Condition</dt><dd>%s</dd>'
            '<dt>Decision</dt><dd>%s</dd><dt>Total matched</dt><dd>%s</dd>'
            '<dt>Unique to this rule</dt><dd>%s</dd>'
            '<dt>Also matched another rule</dt><dd>%s</dd></dl>'
            '<details open><summary>Detailed trigger counts</summary>%s</details></article>'
            % (
                escape(str(rule.get("category") or "Unclassified")),
                _html_cell(rule.get("number")),
                escape(str(rule.get("label") or "Unnamed rule")),
                escape(format_metric_count(rule.get("failed_raw"))),
                escape(_html_percent(rule.get("failed_raw"), raw_total)),
                escape(str(rule.get("purpose") or "")),
                width,
                escape(str(rule.get("criterion") or "")),
                escape(_decision_label(rule.get("decision"))),
                escape(format_metric_count(rule.get("failed_raw"))),
                escape(format_metric_count(rule.get("unique_only_raw"))),
                escape(format_metric_count(rule.get("overlap_raw"))),
                _html_table(("Trigger", "Decision", "Input matches"), detail_rows),
            )
        )
    return '<div class="rule-list">%s</div>' % "".join(cards)


def _html_key_findings(assessment: Mapping[str, Any]) -> str:
    raw = assessment["raw"]
    rules = list(assessment.get("rules") or ())
    findings = [
        "%s of %s input variants (%s) would remain in the virtual QC-passed "
        "subset; %s (%s) would be excluded by at least one active condition."
        % (
            format_metric_count(assessment["qc_passed"].get("num_records")),
            format_metric_count(raw.get("num_records")),
            _html_percent(
                assessment["qc_passed"].get("num_records"), raw.get("num_records"),
            ),
            format_metric_count(assessment.get("excluded_total")),
            _html_percent(assessment.get("excluded_total"), raw.get("num_records")),
        ),
    ]
    if rules:
        largest = max(rules, key=lambda rule: int(rule.get("failed_raw") or 0))
        findings.append(
            "The largest individual rule-match count is %s: %s variants (%s of "
            "input). Rule counts overlap, so this is not an independent removal count."
            % (
                largest.get("label"),
                format_metric_count(largest.get("failed_raw")),
                _html_percent(largest.get("failed_raw"), raw.get("num_records")),
            )
        )
    score_rule = next(
        (rule for rule in rules if rule.get("key") == "imputation_quality"), {},
    )
    if score_rule:
        missing = _rule_detail(score_rule, "missing").get("matched_raw")
        low = _rule_detail(score_rule, "below_minimum").get("matched_raw")
        high = _rule_detail(score_rule, "above_maximum").get("matched_raw")
        findings.append(
            "For the imputation quality score, %s variants are missing, %s are "
            "below the minimum, and %s are above the maximum. Interpret these "
            "counts using the VCF provenance shown below."
            % tuple(format_metric_count(value) for value in (missing, low, high))
        )
    findings.append(
        "%s variants match more than one active condition; per-rule totals must "
        "not be added together."
        % format_metric_count(assessment.get("overlap_variants"))
    )
    return '<ul class="findings">%s</ul>' % "".join(
        "<li>%s</li>" % escape(item) for item in findings
    )


def _html_provenance(assessment: Mapping[str, Any]) -> str:
    provenance = assessment.get("vcf_provenance") or {}
    values = provenance.get("values") or {}
    rows = (
        ("Provenance status", provenance.get("status")),
        ("PostGWAS version", values.get("postgwas_version")),
        ("PostGWAS dataset", values.get("dataset_id")),
        ("VCF status", values.get("vcf_status")),
        ("Allele-frequency source", values.get("af_source")),
        ("Allele-frequency input column", values.get("af_input_column")),
        ("Allele-frequency meaning", values.get("af_input_type")),
        ("Allele-frequency decision", values.get("af_type_decision")),
        ("Allele-frequency harmonisation", values.get("af_harmonisation")),
        ("Allele-frequency output field", values.get("af_output_field")),
        ("Imputation quality score source", values.get("info_source")),
        ("Imputation quality score input column", values.get("info_input_column")),
        ("Imputation quality score interpretation", values.get("info_interpretation")),
        ("Imputation quality score output field", values.get("info_output_field")),
    )
    note = (
        "These values explain what the assessed AF and imputation-quality "
        "fields mean, including whether the score is a study measurement or "
        "an external proxy."
        if provenance.get("status") == "available"
        else "Complete scientific provenance was not available in the configured "
        "VCF metadata. Interpret the AF and score rules cautiously."
    )
    return '<p class="section-note">%s</p>%s' % (
        escape(note), _html_key_values(rows),
    )


def render_qc_summary_html(assessment: Mapping[str, Any]) -> str:
    """Render a self-contained report from the reconciled assessment evidence."""
    raw = assessment["raw"]
    passed = assessment["qc_passed"]
    retained = 100.0 * float(assessment.get("retained_fraction") or 0.0)
    cards = (
        ("Input variants", format_metric_count(raw.get("num_records"))),
        ("Virtual QC-passed", format_metric_count(passed.get("num_records"))),
        ("Excluded by ≥1 rule", format_metric_count(assessment.get("excluded_total"))),
        ("Retained", "%.2f%%" % retained),
        ("Active rules", format_metric_count(assessment.get("active_rule_count"))),
        ("Matched multiple rules", format_metric_count(assessment.get("overlap_variants"))),
    )
    metric_cards = '<div class="metrics">%s</div>' % "".join(
        '<div class="metric"><div class="metric-label">%s</div>'
        '<div class="metric-value">%s</div></div>'
        % (escape(label), escape(value)) for label, value in cards
    )
    header = assessment.get("vcf_header_validation") or {}
    header_rows = (
        ("Validation status", header.get("status")),
        ("Required field count", header.get("required_field_count")),
        ("Required declared fields", header.get("required_fields")),
        ("Declared INFO fields", header.get("declared_info_field_count")),
        ("Declared FORMAT fields", header.get("declared_format_field_count")),
    )
    input_validation_rows = (
        ("Dataset", assessment.get("dataset_id")),
        (
            "Input VCF",
            os.path.basename(str(assessment.get("raw_vcf") or "")),
        ),
        ("Genome build", assessment.get("genome_build")),
        ("Declared contigs", header.get("declared_contig_count")),
        ("Total variants", raw.get("num_records")),
        ("Reference-frequency tag", assessment.get("external_af_name")),
        *header_rows,
    )
    accounting_rows = (
        ("Input variants", raw.get("num_records")),
        ("Virtual QC-passed variants", passed.get("num_records")),
        ("Variants excluded by ≥1 rule", assessment.get("excluded_total")),
        ("Variants unique to one rule", assessment.get("unique_rule_only_total")),
        ("Variants matching multiple rules", assessment.get("overlap_variants")),
        ("Retained percentage", "%.2f%%" % retained),
        ("Accounting reconciled", assessment.get("accounting_balanced")),
    )
    metric_rows = [
        (label, raw.get(key), passed.get(key))
        for label, key in (
            ("Variant records", "num_records"),
            ("SNPs", "num_snps"),
            ("Indels and other non-SNPs", "num_non_snps"),
            ("Transition/transversion ratio", "ts_tv_ratio"),
            ("Missing FORMAT/AF", "format_af_missing"),
            ("Missing imputation quality score", "format_si_missing"),
            ("Missing study INFO/AF", "study_af_missing"),
            ("Missing reference AF", "external_af_missing"),
            ("Comparable study/reference AF pairs", "af_comparable"),
            ("AF differences above cutoff", "af_difference_above_cutoff"),
            ("Usable effective sample sizes", "effective_sample_size_available"),
            ("Missing/invalid effective sample sizes", "effective_sample_size_missing_or_invalid"),
            ("Minimum effective sample size", "effective_sample_size_minimum"),
            ("Maximum effective sample size", "effective_sample_size_maximum"),
            ("Mean effective sample size", "effective_sample_size_mean"),
            ("Effective sample-size sample SD", "effective_sample_size_standard_deviation"),
            ("Low-Neff threshold", "effective_sample_size_minimum_threshold"),
            ("Below low-Neff threshold", "effective_sample_size_below_minimum_threshold"),
            ("Below low-Neff fraction", "effective_sample_size_below_minimum_threshold_fraction"),
            ("Upper-tail diagnostic threshold", "effective_sample_size_outlier_threshold"),
            ("Upper-tail diagnostic count", "effective_sample_size_above_outlier_threshold"),
        )
    ]
    policy = assessment.get("variant_qc_policy") or {}
    policy_rows = (
        ("Minimum −log10(P)", "Disabled" if policy.get("minimum_neglog10_p") is None else policy.get("minimum_neglog10_p")),
        ("Missing P-value decision", _missing_policy_label(policy.get("missing_pvalue_action"))),
        ("Minimum MAF", policy.get("maf_min")),
        ("Missing AF decision", _missing_policy_label(policy.get("missing_af_action"))),
        ("Minimum imputation quality score", policy.get("info_min")),
        ("Maximum imputation quality score", policy.get("info_max")),
        ("Missing imputation quality score decision", _missing_policy_label(policy.get("missing_info_action"))),
        ("Maximum absolute AF difference", policy.get("maximum_af_difference")),
        ("Include indels/non-SNPs", policy.get("include_indels")),
        ("Exclude frequency-ambiguous palindromic SNPs", policy.get("remove_palindromic")),
        ("Palindromic AF interval", "%s–%s" % (policy.get("palindromic_lower"), policy.get("palindromic_upper"))),
        ("Exclude configured MHC region", policy.get("remove_mhc")),
        ("MHC interval", "%s:%s–%s" % (policy.get("mhc_chromosome"), policy.get("mhc_start"), policy.get("mhc_end")) if policy.get("remove_mhc") else "Inactive"),
    )
    fields = assessment.get("field_labels") or {}
    field_rows = (
        ("Study AF used for concordance", fields.get("study_info_af")),
        ("Reference AF used for concordance", fields.get("external_info_af")),
        ("Study AF used for MAF/palindromic rules", fields.get("study_format_af")),
        ("Imputation quality score", fields.get("imputation_format")),
        ("Association −log10(P)", fields.get("log_pvalue_format")),
        ("Effective sample size", fields.get("effective_sample_size_format")),
    )
    inactive_rows = [
        (item.get("label"), item.get("category"), item.get("reason"))
        for item in assessment.get("inactive_rules") or ()
    ]
    aggregation = assessment.get("aggregation") or {}
    execution_rows = (
        ("Aggregation strategy", aggregation.get("aggregation_strategy")),
        ("Temporary-table scans", aggregation.get("temporary_table_scans")),
        ("Streaming collection API", aggregation.get("streaming_collection_api")),
        ("Polars version", aggregation.get("polars_version")),
        ("Requested threads", aggregation.get("requested_threads")),
        ("Effective Polars threads", aggregation.get("polars_thread_pool_size")),
        ("Thread budget enforced", aggregation.get("thread_budget_enforced")),
    )
    reports = assessment.get("reports") or {}
    output_rows = (
        ("Input GWAS-VCF", assessment.get("raw_vcf")),
        ("Metric TSV", reports.get("summary")),
        ("Rule TSV", reports.get("rules")),
        ("Assessment JSON", reports.get("json")),
        ("Summary CSV", reports.get("csv")),
        ("Detailed HTML report", reports.get("html")),
    )
    status_class = "ok" if assessment.get("accounting_balanced") else "fail"
    status = "RECONCILED" if assessment.get("accounting_balanced") else "FAILED"
    provenance_status = (assessment.get("vcf_provenance") or {}).get("status")
    provenance_class = "" if provenance_status == "available" else " warning"
    sections = "".join((
        (
            '<section class="%s"><h2>Validated input summary-statistics '
            'GWAS-VCF</h2><p class="section-note">The build and configured '
            'INFO/FORMAT declarations were validated before extraction, '
            'declared contigs were counted from the header, and the total '
            'variant count was measured during the single extraction used by '
            'this assessment.</p>%s<h3>Scientific field provenance</h3>%s'
            '</section>'
        )
        % (
            provenance_class.strip(),
            _html_key_values(input_validation_rows),
            _html_provenance(assessment),
        ),
        '<section><h2>What these results mean</h2><p class="section-note">A deterministic interpretation of the reconciled counts; no additional filtering was performed.</p>%s</section>'
        % _html_key_findings(assessment),
        '<section><h2>Active QC rules</h2><p class="section-note">Each card separates the scientific purpose, configured decision, total matches, variants unique to that rule, and overlap with other rules. “Excluded” means excluded only from the report’s virtual subset.</p>%s</section>'
        % _html_rule_cards(assessment),
        '<section><h2>Inactive QC rules</h2><p class="section-note">These configurable conditions do not contribute to the combined virtual-subset decision.</p>%s</section>'
        % _html_table(("Rule", "Category", "Why inactive"), inactive_rows),
        '<section><h2>Exact QC accounting</h2><p class="section-note">Rule matches may overlap. Only the combined mask defines the unique excluded and virtual QC-passed totals.</p>%s</section>'
        % _html_key_values(accounting_rows),
        '<section><h2>Input versus virtual QC-passed metrics</h2><p class="section-note">Both columns come from the same two-pass streaming assessment; the VCF is not reread to build this report.</p>%s</section>'
        % _html_table(("Metric", "Input GWAS-VCF", "Virtual QC-passed subset"), metric_rows),
        '<section><h2>Resolved scientific policy</h2>%s</section>' % _html_key_values(policy_rows),
        '<section><h2>Resolved VCF fields</h2>%s</section>' % _html_key_values(field_rows),
        '<section><h2>Execution provenance</h2>%s</section>' % _html_key_values(execution_rows),
        '<section><h2>Outputs</h2>%s</section>' % _html_key_values(output_rows, path_labels=tuple(label for label, _ in output_rows)),
    ))
    return """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>%s — PostGWAS QC report</title>
<style>%s</style>
</head>
<body>
<header class="page-header"><p class="eyebrow">PostGWAS</p><h1>GWAS-VCF quality-control report</h1><p>%s · %s · input evidence, virtual-subset decisions, rule overlap, and scientific provenance.</p></header>
<main class="container">
<div class="notice"><strong>Report-only assessment:</strong> the input GWAS-VCF is unchanged. “QC-passed” describes a virtual subset; no filtered VCF was created.</div>
%s
<div class="notice"><strong>Count reconciliation:</strong> <span class="badge %s">%s</span> %s input = %s virtual QC-passed + %s excluded.</div>
%s
</main>
<footer class="footer">PostGWAS QC report · consult the linked CSV, TSV, JSON, and canonical log for complete machine-readable evidence.</footer>
</body>
</html>
""" % (
        escape(str(assessment.get("dataset_id") or "dataset")),
        _QC_HTML_STYLES,
        escape(str(assessment.get("dataset_id") or "dataset")),
        escape(str(assessment.get("genome_build") or "unknown build")),
        metric_cards,
        status_class,
        status,
        _html_cell(raw.get("num_records")),
        _html_cell(passed.get("num_records")),
        _html_cell(assessment.get("excluded_total")),
        sections,
    )


__all__ = [
    "QC_SUMMARY_COLUMNS",
    "format_metric_count",
    "format_metric_percent",
    "metric_available",
    "metric_warning_kind",
    "qc_final_screen_lines",
    "qc_input_validation_screen_lines",
    "qc_rule_results_screen_lines",
    "qc_rules_applied_screen_lines",
    "qc_summary_csv_records",
    "qc_summary_lines",
    "qc_summary_screen_lines",
    "render_qc_summary_html",
]
