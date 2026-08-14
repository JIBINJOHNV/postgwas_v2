"""Shared screen report for genotype-free GWAS-VCF quality control."""

from __future__ import annotations

import os
from typing import Any

from postgwas.core.ui.screen import screen_field, screen_line


_QC_COLON_INDEX = 66


def _qc_label_width(indent: int) -> int:
    """Pad labels so every QC field colon has one terminal column."""
    return _QC_COLON_INDEX - int(indent) - 4


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


def _field(kind: str, label: str, value: Any, indent: int) -> list[str]:
    return screen_field(
        kind, label, value, indent=indent, label_width=_qc_label_width(indent),
    ).splitlines()


def qc_summary_lines(
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
    lines = [
        screen_line("analysis", title, indent=4),
        screen_field(
            "info", "Assessment model",
            "the merged VCF is extracted once; every active configured QC condition "
            "is evaluated independently against the same raw records",
            indent=6, label_width=_qc_label_width(6),
        ),
        screen_field(
            "info", "VCF output",
            "the raw merged VCF is unchanged; no QC-filtered VCF is created",
            indent=6, label_width=_qc_label_width(6),
        ),
        "",
        screen_line("genetic", "1. Before VCF creation", indent=6),
    ]
    for label, key, kind in (
        ("Input rows", "total_variant_infile", "count"),
        ("Rows read", "total_variant_read", "count"),
        ("Missing required fields", "total_variant_removed_missing_values", "loss"),
        ("Duplicate rows", "total_variant_removed_duplicates", "loss"),
        ("Invalid coordinates", "total_variant_removed_null_coords", "loss"),
        ("Non-standard alleles", "total_variant_removed_non_standard_alleles", "loss"),
        ("Passed initial input validation", "total_variant_remaining_for_harmonisation", "success"),
        ("Ready SNPs", "total_variant_ready_snps", "genetic"),
        ("Ready indels / other variants", "total_variant_ready_indels_or_other", "genetic"),
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

    lines.extend([
        "",
        screen_line("genetic", "2. Raw merged VCF", indent=6),
    ])
    lines.extend(_field("info", "File", os.path.basename(assessment["raw_vcf"]), 8))
    lines.extend(_field(
        "count", "Variant records", format_metric_count(raw["num_records"]), 8,
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
        metric_warning_kind(raw["effective_sample_size_above_outlier_threshold"]),
        "Above upper outlier threshold",
        format_metric_count(raw["effective_sample_size_above_outlier_threshold"]), 12,
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
        screen_line("analysis", "3. QC conditions assessed on the raw VCF", indent=6),
        screen_field(
            "info", "Assessment basis",
            "all %s raw variants are tested against every active condition. Rule "
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
        lines.extend(_field(
            "loss" if rule["failed_raw"] else "success",
            "Raw variants failing", format_metric_count(rule["failed_raw"]), 12,
        ))
        for detail in rule["details"]:
            action = "fails QC"
            if detail["action"] == "pass_qc":
                action = "passes by configured policy"
            lines.extend(_field(
                (
                    "loss"
                    if detail["matched_raw"] and detail["action"] == "fail_qc"
                    else "info"
                ),
                detail["label"],
                "%s in raw VCF; %s"
                % (format_metric_count(detail["matched_raw"]), action),
                12,
            ))

    retained_percent = 100.0 * assessment["retained_fraction"]
    lines.extend([
        "",
        screen_line("success", "4. Final QC-passed assessment", indent=6),
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
        "count", "Raw VCF variants", format_metric_count(raw["num_records"]), 8,
    ))
    lines.extend(_field(
        "success", "QC-passed variants",
        "%s (%.2f%% of raw VCF)"
        % (format_metric_count(passed["num_records"]), retained_percent),
        8,
    ))
    lines.extend(_field(
        "loss", "Failed combined QC",
        format_metric_count(assessment["excluded_total"]),
        8,
    ))
    lines.extend(_field(
        "analysis", "Failed multiple conditions",
        format_metric_count(assessment["overlap_variants"]),
        8,
    ))
    lines.extend(_field(
        "analysis", "Rule failure matches",
        "%s total = %s unique failed variants + %s additional overlapping matches"
        % (
            format_metric_count(assessment["rule_match_total"]),
            format_metric_count(assessment["excluded_total"]),
            format_metric_count(assessment["extra_rule_matches"]),
        ),
        8,
    ))
    lines.extend(_field(
        "success", "Count reconciliation",
        "%s raw = %s QC-passed + %s failed combined QC" % (
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
            "%s values above threshold" % format_neff,
            "effective_sample_size_above_outlier_threshold", "warning",
            format_metric_count,
        ),
    ):
        resolved_kind = (
            metric_warning_kind(passed[key]) if kind == "warning" else kind
        )
        lines.extend(_field(resolved_kind, label, formatter(passed[key]), 12))
    lines.extend(_field(
        "info", "QC-passed VCF", "not created; use the unchanged raw merged VCF",
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
    if not include_pre_vcf:
        pre_start = next(
            index for index, line in enumerate(lines)
            if "1. Before VCF creation" in line
        )
        raw_start = next(
            index for index, line in enumerate(lines)
            if "2. Raw merged VCF" in line
        )
        del lines[pre_start:raw_start]
        lines = [
            line.replace("2. Raw merged VCF", "1. Raw merged VCF")
            .replace(
                "3. QC conditions assessed on the raw VCF",
                "2. QC conditions assessed on the raw VCF",
            )
            .replace(
                "4. Final QC-passed assessment",
                "3. Final QC-passed assessment",
            )
            for line in lines
        ]
    return lines


__all__ = [
    "format_metric_count",
    "format_metric_percent",
    "metric_available",
    "metric_warning_kind",
    "qc_summary_lines",
]
