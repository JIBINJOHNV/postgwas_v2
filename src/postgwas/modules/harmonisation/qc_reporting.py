"""Hierarchical screen report for harmonisation and virtual VCF QC assessment."""

from __future__ import annotations

import os
from typing import Any

from postgwas.core.ui.screen import screen_field, screen_line


_QC_COLON_INDEX = 66


def _qc_label_width(indent: int) -> int:
    """Pad labels so every QC field colon has one terminal column."""
    return _QC_COLON_INDEX - int(indent) - 4


def _available(value: Any) -> bool:
    if value is None:
        return False
    try:
        return bool(value == value)
    except (TypeError, ValueError):
        return False


def _count(value: Any) -> str:
    if not _available(value):
        return "unavailable"
    try:
        return "{:,}".format(int(value))
    except (TypeError, ValueError):
        return str(value)


def _ratio(value: Any) -> str:
    if not _available(value):
        return "unavailable"
    return "%.2f" % float(value)


def _number(value: Any) -> str:
    if not _available(value):
        return "unavailable"
    return "{:,.2f}".format(float(value))


def _percent(numerator: Any, denominator: Any, decimals: int = 2) -> str:
    """Format a percentage without hiding an unavailable or zero denominator."""
    if not _available(numerator) or not _available(denominator):
        return "unavailable"
    try:
        denominator_value = float(denominator)
        if denominator_value <= 0:
            return "unavailable"
        value = 100.0 * float(numerator) / denominator_value
    except (TypeError, ValueError, ZeroDivisionError):
        return "unavailable"
    return ("%%.%df%%%%" % int(decimals)) % value


def _warning_if_nonzero(value: Any) -> str:
    try:
        return "warning" if int(value) else "success"
    except (TypeError, ValueError):
        return "warning"


def _field(kind: str, label: str, value: Any, indent: int) -> list[str]:
    return screen_field(
        kind, label, value, indent=indent, label_width=_qc_label_width(indent),
    ).splitlines()


def harmonisation_qc_summary_lines(
    pre_vcf: dict[str, Any],
    assessment: dict[str, Any],
) -> list[str]:
    """Render pre-VCF, raw-rule, and combined final QC outcomes."""
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
        screen_line("analysis", "Harmonisation quality-control summary", indent=4),
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
        ("Palindromic ambiguous removed", "total_variant_removed_palindromic_ambiguous", "loss"),
        ("Reference unmatched removed", "total_variant_removed_reference_unmatched", "loss"),
        ("Reference ambiguous removed", "total_variant_removed_reference_ambiguous", "loss"),
        ("Missing effect frequency", "total_variant_with_missing_eaf", "warning"),
        ("Invalid effect statistics", "total_variant_with_invalid_beta_se", "loss"),
        ("All chromosome-stage removals", "total_variant_removed_chromosome_harmonisation", "loss"),
        ("Passed chromosome harmonisation", "total_variant_passed_chromosome_harmonisation", "success"),
        ("Sent to GWAS-to-VCF", "total_variant_in_vcf_input", "count"),
    ):
        if key in pre_vcf:
            lines.extend(_field(kind, label, _count(pre_vcf.get(key)), 8))

    initial = pre_vcf.get("total_variant_remaining_for_harmonisation")
    removed = pre_vcf.get("total_variant_removed_chromosome_harmonisation")
    passed_chromosomes = pre_vcf.get(
        "total_variant_passed_chromosome_harmonisation"
    )
    sent = pre_vcf.get("total_variant_in_vcf_input")
    if all(
        _available(value)
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
                _count(initial), _count(sent), _count(removed),
                "" if balanced else " — DOES NOT BALANCE",
            ),
            8,
        ))

    lines.extend([
        "",
        screen_line("genetic", "2. Raw merged VCF", indent=6),
    ])
    lines.extend(_field("info", "File", os.path.basename(assessment["raw_vcf"]), 8))
    lines.extend(_field("count", "Variant records", _count(raw["num_records"]), 8))
    lines.extend(_field("count", "SNPs", _count(raw["num_snps"]), 8))
    lines.extend(_field("count", "Indels / other variants", _count(raw["num_non_snps"]), 8))
    lines.extend(_field("analysis", "Transition / transversion", _ratio(raw["ts_tv_ratio"]), 8))
    lines.append("")
    lines.append(screen_line("analysis", "Effective sample-size distribution", indent=10))
    lines.extend(_field(
        "count", "Usable %s" % format_neff,
        _count(raw["effective_sample_size_available"]), 12,
    ))
    lines.extend(_field(
        _warning_if_nonzero(raw["effective_sample_size_missing_or_invalid"]),
        "Missing or invalid %s" % format_neff,
        _count(raw["effective_sample_size_missing_or_invalid"]), 12,
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
        _warning_if_nonzero(raw["effective_sample_size_above_outlier_threshold"]),
        "Above upper outlier threshold",
        _count(raw["effective_sample_size_above_outlier_threshold"]), 12,
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
            _warning_if_nonzero(value)
            if "Missing" in label or key == "af_difference_above_cutoff"
            else "count"
        )
        lines.extend(_field(kind, label, _count(value), 12))

    lines.extend([
        "",
        screen_line("analysis", "3. QC conditions assessed on the raw VCF", indent=6),
        screen_field(
            "info", "Assessment basis",
            "all %s raw variants are tested against every active condition. Rule "
            "counts may overlap and must not be added together"
            % _count(raw["num_records"]),
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
            "Raw variants failing", _count(rule["failed_raw"]), 12,
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
                "%s in raw VCF; %s" % (_count(detail["matched_raw"]), action),
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
            % _count(assessment["active_rule_count"]),
            indent=8, label_width=_qc_label_width(8),
        ),
    ])
    lines.extend(_field("count", "Raw VCF variants", _count(raw["num_records"]), 8))
    lines.extend(_field(
        "success", "QC-passed variants",
        "%s (%.2f%% of raw VCF)" % (_count(passed["num_records"]), retained_percent),
        8,
    ))
    lines.extend(_field(
        "loss", "Failed combined QC", _count(assessment["excluded_total"]),
        8,
    ))
    lines.extend(_field(
        "analysis", "Failed multiple conditions", _count(assessment["overlap_variants"]),
        8,
    ))
    lines.extend(_field(
        "analysis", "Rule failure matches",
        "%s total = %s unique failed variants + %s additional overlapping matches"
        % (
            _count(assessment["rule_match_total"]),
            _count(assessment["excluded_total"]),
            _count(assessment["extra_rule_matches"]),
        ),
        8,
    ))
    lines.extend(_field(
        "success", "Count reconciliation",
        "%s raw = %s QC-passed + %s failed combined QC" % (
            _count(raw["num_records"]),
            _count(passed["num_records"]),
            _count(assessment["excluded_total"]),
        ),
        8,
    ))
    lines.extend([
        "",
        screen_line("genetic", "Final virtual subset metrics", indent=10),
    ])
    for label, key, kind, formatter in (
        ("SNPs", "num_snps", "count", _count),
        ("Indels / other variants", "num_non_snps", "count", _count),
        ("Transition / transversion", "ts_tv_ratio", "analysis", _ratio),
        ("Missing %s" % format_af, "format_af_missing", "warning", _count),
        ("Missing %s" % format_info, "format_si_missing", "warning", _count),
        ("Missing %s" % external_info_af, "external_af_missing", "warning", _count),
        ("AF differences above cutoff", "af_difference_above_cutoff", "warning", _count),
        ("Usable %s" % format_neff, "effective_sample_size_available", "count", _count),
        (
            "Missing or invalid %s" % format_neff,
            "effective_sample_size_missing_or_invalid", "warning", _count,
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
            "effective_sample_size_above_outlier_threshold", "warning", _count,
        ),
    ):
        resolved_kind = _warning_if_nonzero(passed[key]) if kind == "warning" else kind
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
    return lines


def harmonisation_qc_takeaway_lines(
    pre_vcf: dict[str, Any],
    assessment: dict[str, Any],
    *,
    genome_build_info: dict[str, Any],
    study_decisions: dict[str, Any],
    population_frequency: dict[str, Any],
    dataset_summary: dict[str, Any],
    merge_summary: dict[str, Any],
    primary_outputs: dict[str, str],
    final_status: str,
) -> list[str]:
    """Render six compact QC cards using only previously calculated results."""
    raw = assessment["raw"]
    passed = assessment["qc_passed"]
    lines = [screen_line("decision", "Final QC takeaways", indent=4)]

    def add_card(
        kind: str,
        title: str,
        details: list[tuple[str, str, str]],
    ) -> None:
        """Append one visually separated takeaway and its supporting evidence."""
        lines.append("")
        lines.append(screen_line(kind, title, indent=8))
        for detail_kind, label, value in details:
            lines.extend(screen_field(
                detail_kind,
                label,
                value,
                indent=12,
                label_width=20,
            ).splitlines())

    input_rows = pre_vcf.get("total_variant_infile")
    rows_read = pre_vcf.get("total_variant_read")
    sent = pre_vcf.get("total_variant_in_vcf_input")
    final_records = raw.get("num_records")
    total_snps = raw.get("num_snps")
    total_non_snps = raw.get("num_non_snps")
    passed_snps = passed.get("num_snps")
    failed_snps = (
        max(int(total_snps) - int(passed_snps), 0)
        if _available(total_snps) and _available(passed_snps)
        else None
    )
    variant_balanced = (
        all(
            _available(value)
            for value in (sent, final_records, total_snps, passed_snps)
        )
        and int(sent) == int(final_records)
        and 0 <= int(passed_snps) <= int(total_snps)
        and bool(assessment.get("accounting_balanced", False))
    )
    add_card(
        "success" if variant_balanced else "warning",
        "Variant flow",
        [
            (
                "info",
                "Input / read",
                "%s / %s variants" % (_count(input_rows), _count(rows_read)),
            ),
            (
                "info",
                "VCF creation",
                "%s variants used → %s variants in final VCF"
                % (_count(sent), _count(final_records)),
            ),
            (
                "genetic",
                "VCF content",
                "%s SNPs · %s indels / other variants"
                % (_count(total_snps), _count(total_non_snps)),
            ),
            (
                "analysis",
                "QC passed",
                "%s / %s SNPs passed all %s active rules (%s)"
                % (
                    _count(passed_snps), _count(total_snps),
                    _count(assessment.get("active_rule_count")),
                    _percent(passed_snps, total_snps),
                ),
            ),
            (
                (
                    "warning"
                    if _available(failed_snps) and int(failed_snps)
                    else "success"
                ),
                "QC failed",
                "%s SNPs failed ≥1 active rule (%s)"
                % (
                    _count(failed_snps),
                    _percent(failed_snps, total_snps),
                ),
            ),
        ],
    )

    genome_build = str(
        genome_build_info.get("inferred_build") or "unavailable"
    )
    build_matches = (genome_build_info.get("matches") or {}).get(genome_build)
    testable = genome_build_info.get("testable_variants")
    strand = str(study_decisions.get("strand") or "unresolved")
    strand_consensus = study_decisions.get("strand_consensus") or {}
    dominant_fraction = strand_consensus.get("dominant_fraction")
    strand_text = strand.capitalize()
    if _available(dominant_fraction):
        strand_text += " · %s dominant" % _percent(dominant_fraction, 1)
    build_resolved = (
        genome_build != "unavailable"
        and _available(build_matches)
        and _available(testable)
        and strand in {"forward", "reverse"}
    )
    add_card(
        "success" if build_resolved else "warning",
        "Reference alignment",
        [
            (
                "genetic",
                "Genome build",
                "%s · %s / %s testable variants matched (%s)"
                % (
                    genome_build, _count(build_matches), _count(testable),
                    _percent(build_matches, testable),
                ),
            ),
            ("genetic", "Strand", strand_text),
        ],
    )

    frequency_type = str(
        study_decisions.get("frequency_type") or "unresolved"
    ).replace("_", " ")
    closest_population = population_frequency.get("closest_population")
    selected_check = population_frequency.get("selected_population_check") or {}
    selected_population = (
        selected_check.get("normalized_population")
        or assessment.get("external_af_name")
        or "unavailable"
    )
    selected_status = selected_check.get("status")
    if closest_population:
        if selected_status == "match":
            population_text = (
                "%s identified as closest population and selected as reference"
                % closest_population
            )
        elif selected_status == "mismatch":
            population_text = (
                "%s identified as closest population, but %s selected as reference"
                % (closest_population, selected_population)
            )
        elif selected_status:
            population_text = (
                "%s identified as closest population; selected %s not comparable"
                % (closest_population, selected_population)
            )
        else:
            population_text = "%s identified as closest population" % closest_population
    else:
        population_text = "Closest population inconclusive"

    filename_checks = population_frequency.get("external_file_checks") or []
    if filename_checks:
        filename_matches = sum(
            check.get("status") == "match" for check in filename_checks
        )
        filename_mismatches = sum(
            check.get("status") in {"mismatch", "ambiguous"}
            for check in filename_checks
        )
        filename_unresolved = len(filename_checks) - filename_matches - filename_mismatches
        population_text += " · filenames: %s match, %s mismatch, %s unresolved" % (
            _count(filename_matches), _count(filename_mismatches),
            _count(filename_unresolved),
        )

    comparable = raw.get("af_comparable")
    mismatched = raw.get("af_difference_above_cutoff")
    concordant = (
        max(int(comparable) - int(mismatched), 0)
        if _available(comparable) and _available(mismatched)
        else None
    )
    study_missing = raw.get("study_af_missing")
    reference_missing = raw.get("external_af_missing")
    cutoff = assessment.get("af_difference_cutoff")
    frequency_warning = (
        frequency_type == "unresolved"
        or closest_population is None
        or selected_status == "mismatch"
        or any(
            check.get("status") in {"mismatch", "ambiguous"}
            for check in filename_checks
        )
        or any(
            not _available(value) or int(value) != 0
            for value in (mismatched, study_missing, reference_missing)
        )
    )
    add_card(
        "warning" if frequency_warning else "success",
        "Allele frequency",
        [
            ("genetic", "Study frequency", frequency_type.capitalize()),
            ("decision", "Population", population_text),
            (
                "analysis",
                "Reference AF",
                "%s / %s concordant (%s) · %s mismatched (%s)"
                % (
                    _count(concordant), _count(comparable),
                    _percent(concordant, comparable, 4), _count(mismatched),
                    _percent(mismatched, comparable, 4),
                ),
            ),
            (
                "info",
                "AF check",
                "%s missing study AF · %s missing reference AF · "
                "|AF difference| ≤ %s"
                % (
                    _count(study_missing), _count(reference_missing),
                    "unavailable" if not _available(cutoff) else "%g" % float(cutoff),
                ),
            ),
        ],
    )

    effect_type = str(study_decisions.get("effect_type") or "unresolved")
    effect_text = (
        "odds ratio harmonised to log-odds BETA"
        if effect_type == "odds_ratio"
        else "BETA-scale effect"
        if effect_type == "beta"
        else "effect type unresolved"
    )
    invalid_effects = pre_vcf.get("total_variant_with_invalid_beta_se")
    missing_neff = passed.get("effective_sample_size_missing_or_invalid")
    missing_info = passed.get("format_si_missing")
    neff_outliers = passed.get("effective_sample_size_above_outlier_threshold")
    statistical_warning = any(
        not _available(value) or int(value) != 0
        for value in (invalid_effects, missing_neff, missing_info, neff_outliers)
    ) or effect_type not in {"beta", "odds_ratio"}
    add_card(
        "warning" if statistical_warning else "success",
        "Statistical quality",
        [
            ("analysis", "Effect scale", effect_text),
            (
                "analysis",
                "Pre-VCF",
                "%s invalid effect-statistic removals" % _count(invalid_effects),
            ),
            (
                "analysis",
                "QC-passed variants",
                "%s missing/invalid Neff · %s missing imputation score · "
                "%s Neff upper outliers"
                % (
                    _count(missing_neff), _count(missing_info),
                    _count(neff_outliers),
                ),
            ),
        ],
    )

    reconciliation = dataset_summary.get("reconciliation") or {}
    rejected_rows = dataset_summary.get("rejected_variants_rows")
    rejection_file = dataset_summary.get("rejected_variants_file")
    reconciled = reconciliation.get("balanced") is True
    accounting_balanced = bool(assessment.get("accounting_balanced", False))
    provenance_complete = (
        rejection_file is not None
        and _available(rejected_rows)
        and reconciled
        and accounting_balanced
    )
    add_card(
        "success" if provenance_complete else "warning",
        "Audit trail",
        [
            (
                "loss",
                "Harmonisation",
                "%s variants removed with reasons recorded" % _count(rejected_rows),
            ),
            (
                "analysis",
                "QC assessment",
                "%s merged-VCF variants failed ≥1 of %s active rules"
                % (
                    _count(assessment.get("excluded_total")),
                    _count(assessment.get("active_rule_count")),
                ),
            ),
            (
                "success" if reconciled and accounting_balanced else "warning",
                "Accounting",
                (
                    "Variant and chromosome counts reconciled"
                    if reconciled and accounting_balanced
                    else "Variant or chromosome counts did not reconcile"
                ),
            ),
        ],
    )

    completed = list(dataset_summary.get("completed") or [])
    failed = list(dataset_summary.get("failed") or [])
    chromosome_total = len(dataset_summary.get("chromosomes") or {})
    if chromosome_total == 0:
        chromosome_total = len(completed) + len(failed)
    merge_status = str(merge_summary.get("merge_status") or "unavailable").upper()
    required_failures = list(merge_summary.get("required_merge_failures") or [])
    integrity_ok = (
        str(final_status).upper() == "OK"
        and not failed
        and len(completed) == chromosome_total
        and merge_status == "OK"
        and not required_failures
        and bool(primary_outputs)
    )
    merged_vcf_ok = merge_status == "OK" and not required_failures
    merge_text = "Indexed and validated" if merged_vcf_ok else (
        "Required merge failures: %s"
        % (", ".join(required_failures) if required_failures else merge_status)
    )
    add_card(
        "success" if integrity_ok else "warning",
        "VCF integrity",
        [
            (
                "genetic",
                "Chromosomes",
                "%s / %s completed"
                % (_count(len(completed)), _count(chromosome_total)),
            ),
            (
                "success" if merged_vcf_ok else "warning",
                "Final merged VCF",
                merge_text,
            ),
            (
                "success" if integrity_ok else "warning",
                "Outputs",
                "%s primary outputs completed · dataset status %s"
                % (
                    _count(len(primary_outputs)), str(final_status).upper(),
                ),
            ),
        ],
    )
    return lines
