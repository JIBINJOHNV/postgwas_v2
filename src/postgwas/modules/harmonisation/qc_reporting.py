"""Hierarchical screen report for harmonisation and virtual VCF QC assessment."""

from __future__ import annotations

import math
from typing import Any

from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.modules.qc_summary.reporting import (
    format_metric_count,
    format_metric_percent,
    metric_available,
    qc_summary_lines,
)


def effect_scale_description(effect_type: Any) -> str:
    """Describe the final harmonised effect scale without reinterpreting data."""
    normalized = str(effect_type or "").strip().lower()
    if normalized == "odds_ratio":
        return "odds ratio harmonised to log-odds BETA"
    if normalized == "beta":
        return "BETA-scale effect"
    return "effect type unresolved"


_ORIENTED_STRAND_ACTIONS = (
    "forward",
    "forward_swapped",
    "reverse_complement",
    "reverse_complement_swapped",
)


def _nonnegative_integer(value: Any) -> int | None:
    """Return one internal count without accepting booleans or negative values."""
    if value is None or isinstance(value, bool):
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    if not math.isfinite(numeric) or numeric < 0 or not numeric.is_integer():
        return None
    return int(numeric)


def summarise_strand_orientation(
    dataset_summary: dict[str, Any],
) -> dict[str, Any]:
    """Combine completed-chromosome strand counters without rescanning variants."""
    completed = [str(value) for value in dataset_summary.get("completed") or []]
    chromosome_state = dataset_summary.get("chromosomes") or {}
    expected = len(chromosome_state) or (
        len(completed) + len(dataset_summary.get("failed") or [])
    )
    summaries = dataset_summary.get("chromosome_summaries") or {}
    actions = {name: 0 for name in (*_ORIENTED_STRAND_ACTIONS, "disabled")}
    removed = {
        "reference_unmatched": 0,
        "palindromic_ambiguous": 0,
        "reference_ambiguous": 0,
    }
    evaluated = 0
    retained = 0
    summarized = 0
    counts_complete = True
    statuses: set[str] = set()
    populations: set[str] = set()
    consensuses: set[str] = set()
    reference_files: set[str] = set()
    population_observations = 0
    consensus_observations = 0
    reference_file_observations = 0

    for chromosome in completed:
        summary = summaries.get(chromosome) or {}
        orientation = (
            ((summary.get("stage_qc") or {}).get("eaf_qc") or {}).get(
                "strand_orientation"
            )
            or {}
        )
        if not orientation:
            continue
        summarized += 1
        orientation_status = str(
            orientation.get("status") or "unknown"
        ).strip().lower()
        statuses.add(orientation_status)
        population = orientation.get("reference_population_column")
        consensus = orientation.get("study_strand_consensus")
        reference_file = orientation.get("reference_file")
        if population:
            populations.add(str(population))
            population_observations += 1
        if consensus:
            consensuses.add(str(consensus))
            consensus_observations += 1
        if reference_file:
            reference_files.add(str(reference_file))
            reference_file_observations += 1

        for name in actions:
            value = _nonnegative_integer(
                (orientation.get("actions") or {}).get(name, 0)
            )
            if value is None:
                counts_complete = False
            else:
                actions[name] += value
        for name in removed:
            value = _nonnegative_integer(
                orientation.get(name, 0 if orientation_status == "disabled" else None)
            )
            if value is None:
                counts_complete = False
            else:
                removed[name] += value
        initial = _nonnegative_integer(orientation.get("initial_variants"))
        final = _nonnegative_integer(orientation.get("final_variants"))
        if initial is None or final is None:
            counts_complete = False
        else:
            evaluated += initial
            retained += final

    oriented_matched = sum(actions[name] for name in _ORIENTED_STRAND_ACTIONS)
    removed_total = sum(removed.values())
    disabled = statuses == {"disabled"}
    coverage_complete = summarized == len(completed) == expected
    metadata_mixed = len(populations) > 1 or len(consensuses) > 1
    metadata_complete = disabled or (
        population_observations == summarized
        and consensus_observations == summarized
        and reference_file_observations == summarized
    )
    accounting_balanced = (
        counts_complete
        and evaluated == retained + removed_total
        and retained == sum(actions.values())
        and (disabled or retained == oriented_matched)
    )
    if summarized == 0:
        status = "unavailable"
    elif disabled:
        status = "disabled" if coverage_complete else "partial"
    elif "disabled" in statuses or len(statuses) > 1 or metadata_mixed:
        status = "mixed"
    elif not coverage_complete:
        status = "partial"
    elif statuses != {"success"} or not metadata_complete:
        status = "inconsistent"
    elif not accounting_balanced:
        status = "inconsistent"
    else:
        status = "success"

    def one_or_mixed(values: set[str]) -> str | None:
        if not values:
            return None
        return next(iter(values)) if len(values) == 1 else "mixed"

    return {
        "status": status,
        "chromosomes_expected": expected,
        "chromosomes_completed": len(completed),
        "chromosomes_summarized": summarized,
        "complete_chromosome_coverage": coverage_complete,
        "metadata_complete": metadata_complete,
        "reference_population": one_or_mixed(populations),
        "reference_file_count": len(reference_files),
        "reference_files": sorted(reference_files),
        "study_consensus": one_or_mixed(consensuses),
        "variants_evaluated": evaluated if counts_complete and summarized else None,
        "variants_retained": retained if counts_complete and summarized else None,
        "variants_matched": (
            None if disabled or not counts_complete or not summarized
            else oriented_matched
        ),
        "forward": actions["forward"],
        "forward_swapped": actions["forward_swapped"],
        "reverse_complement": actions["reverse_complement"],
        "reverse_complement_swapped": actions["reverse_complement_swapped"],
        "disabled": actions["disabled"],
        "reference_unmatched": removed["reference_unmatched"],
        "palindromic_ambiguous": removed["palindromic_ambiguous"],
        "reference_ambiguous": removed["reference_ambiguous"],
        "removed_total": removed_total if counts_complete and summarized else None,
        "accounting_balanced": accounting_balanced,
    }


def harmonisation_qc_summary_lines(
    pre_vcf: dict[str, Any],
    assessment: dict[str, Any],
) -> list[str]:
    """Add harmonisation context to the shared GWAS-VCF QC report."""
    return qc_summary_lines(
        assessment,
        pre_vcf=pre_vcf,
        title="Harmonisation quality-control summary",
    )


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
    reference_source: str | None = None,
    reference_population: str | None = None,
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
        if metric_available(total_snps) and metric_available(passed_snps)
        else None
    )
    variant_balanced = (
        all(
            metric_available(value)
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
                "%s / %s variants"
                % (format_metric_count(input_rows), format_metric_count(rows_read)),
            ),
            (
                "info",
                "VCF creation",
                "%s variants used → %s variants in final VCF"
                % (format_metric_count(sent), format_metric_count(final_records)),
            ),
            (
                "genetic",
                "VCF content",
                "%s SNPs · %s indels / other variants"
                % (
                    format_metric_count(total_snps),
                    format_metric_count(total_non_snps),
                ),
            ),
            (
                "analysis",
                "QC passed",
                "%s / %s SNPs passed all %s active rules (%s)"
                % (
                    format_metric_count(passed_snps),
                    format_metric_count(total_snps),
                    format_metric_count(assessment.get("active_rule_count")),
                    format_metric_percent(passed_snps, total_snps),
                ),
            ),
            (
                (
                    "warning"
                    if metric_available(failed_snps) and int(failed_snps)
                    else "success"
                ),
                "QC failed",
                "%s SNPs failed ≥1 active rule (%s)"
                % (
                    format_metric_count(failed_snps),
                    format_metric_percent(failed_snps, total_snps),
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
    if metric_available(dominant_fraction):
        strand_text += " · %s dominant" % format_metric_percent(
            dominant_fraction, 1,
        )
    build_resolved = (
        genome_build != "unavailable"
        and metric_available(build_matches)
        and metric_available(testable)
        and strand in {"forward", "reverse"}
    )
    build_source = (
        "configured explicitly"
        if genome_build_info.get("forced") is True
        else "automatically detected"
    )
    orientation = summarise_strand_orientation(dataset_summary)
    orientation_status = str(orientation["status"])
    reference_population_value = (
        orientation.get("reference_population") or reference_population
    )
    reference_text = " · ".join(
        str(value)
        for value in (reference_source, reference_population_value)
        if value
    ) or "unavailable"
    if orientation.get("reference_file_count"):
        reference_text += " · %s chromosome reference file%s" % (
            format_metric_count(orientation["reference_file_count"]),
            "" if orientation["reference_file_count"] == 1 else "s",
        )
    alignment_details = [
        (
            "genetic",
            "Genome build",
            "%s · %s · %s / %s testable variants matched (%s)"
            % (
                genome_build, build_source,
                format_metric_count(build_matches),
                format_metric_count(testable),
                format_metric_percent(build_matches, testable),
            ),
        ),
        ("genetic", "Strand consensus", strand_text),
    ]
    if orientation_status != "unavailable":
        alignment_details.append(("info", "Reference panel", reference_text))
        if orientation_status == "disabled":
            alignment_details.append((
                "warning",
                "Strand orientation",
                "disabled · %s / %s chromosomes summarized"
                % (
                    format_metric_count(orientation["chromosomes_summarized"]),
                    format_metric_count(orientation["chromosomes_expected"]),
                ),
            ))
        else:
            alignment_details.extend([
                (
                    "success" if orientation_status == "success" else "warning",
                    "Strand matched",
                    "%s total · forward %s · forward-swapped %s · "
                    "reverse-complement %s · reverse-complement-swapped %s"
                    % (
                        format_metric_count(orientation["variants_matched"]),
                        format_metric_count(orientation["forward"]),
                        format_metric_count(orientation["forward_swapped"]),
                        format_metric_count(orientation["reverse_complement"]),
                        format_metric_count(
                            orientation["reverse_complement_swapped"]
                        ),
                    ),
                ),
                (
                    "loss" if orientation["removed_total"] else "success",
                    "Strand removed",
                    "%s total · unmatched %s · palindromic ambiguous %s · "
                    "reference ambiguous %s"
                    % (
                        format_metric_count(orientation["removed_total"]),
                        format_metric_count(orientation["reference_unmatched"]),
                        format_metric_count(orientation["palindromic_ambiguous"]),
                        format_metric_count(orientation["reference_ambiguous"]),
                    ),
                ),
                (
                    "success" if orientation["accounting_balanced"] else "warning",
                    "Strand accounting",
                    "%s evaluated = %s retained + %s removed · %s / %s "
                    "chromosomes summarized"
                    % (
                        format_metric_count(orientation["variants_evaluated"]),
                        format_metric_count(orientation["variants_retained"]),
                        format_metric_count(orientation["removed_total"]),
                        format_metric_count(orientation["chromosomes_summarized"]),
                        format_metric_count(orientation["chromosomes_expected"]),
                    ),
                ),
            ])
    add_card(
        (
            "success"
            if build_resolved and orientation_status in {"success", "unavailable"}
            else "warning"
        ),
        "Reference alignment",
        alignment_details,
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
            format_metric_count(filename_matches),
            format_metric_count(filename_mismatches),
            format_metric_count(filename_unresolved),
        )

    comparable = raw.get("af_comparable")
    mismatched = raw.get("af_difference_above_cutoff")
    concordant = (
        max(int(comparable) - int(mismatched), 0)
        if metric_available(comparable) and metric_available(mismatched)
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
            not metric_available(value) or int(value) != 0
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
                    format_metric_count(concordant),
                    format_metric_count(comparable),
                    format_metric_percent(concordant, comparable, 4),
                    format_metric_count(mismatched),
                    format_metric_percent(mismatched, comparable, 4),
                ),
            ),
            (
                "info",
                "AF check",
                "%s missing study AF · %s missing reference AF · "
                "|AF difference| ≤ %s"
                % (
                    format_metric_count(study_missing),
                    format_metric_count(reference_missing),
                    "unavailable"
                    if not metric_available(cutoff) else "%g" % float(cutoff),
                ),
            ),
        ],
    )

    effect_type = str(study_decisions.get("effect_type") or "unresolved")
    effect_source = str(study_decisions.get("effect_type_source") or "unavailable")
    effect_detected = study_decisions.get("effect_type_detected")
    pvalue_type = str(study_decisions.get("pvalue_type") or "unresolved")
    pvalue_source = str(study_decisions.get("pvalue_type_source") or "unavailable")
    pvalue_detected = study_decisions.get("pvalue_type_detected")

    def decision_text(value: str, source: str, detected: Any) -> str:
        if source == "detector":
            return "%s · automatically detected" % value
        if detected is not None:
            agreement = "detector agrees" if str(detected) == value else (
                "detector found %s" % detected
            )
            return "%s · source %s; %s" % (value, source, agreement)
        return "%s · source %s" % (value, source)

    effect_text = effect_scale_description(effect_type)
    invalid_effects = pre_vcf.get("total_variant_with_invalid_beta_se")
    missing_neff = passed.get("effective_sample_size_missing_or_invalid")
    missing_info = passed.get("format_si_missing")
    neff_outliers = passed.get("effective_sample_size_above_outlier_threshold")
    statistical_warning = any(
        not metric_available(value) or int(value) != 0
        for value in (invalid_effects, missing_neff, missing_info, neff_outliers)
    ) or effect_type not in {"beta", "odds_ratio"}
    add_card(
        "warning" if statistical_warning else "success",
        "Statistical quality",
        [
            (
                "analysis", "Effect type",
                decision_text(effect_type, effect_source, effect_detected),
            ),
            (
                "analysis", "P-value type",
                decision_text(pvalue_type, pvalue_source, pvalue_detected),
            ),
            ("analysis", "Effect scale", effect_text),
            (
                "analysis",
                "Pre-VCF",
                "%s invalid effect-statistic removals"
                % format_metric_count(invalid_effects),
            ),
            (
                "analysis",
                "QC-passed variants",
                "%s missing/invalid Neff · %s missing imputation score · "
                "%s Neff upper outliers"
                % (
                    format_metric_count(missing_neff),
                    format_metric_count(missing_info),
                    format_metric_count(neff_outliers),
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
        and metric_available(rejected_rows)
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
                "%s variants removed with reasons recorded"
                % format_metric_count(rejected_rows),
            ),
            (
                "analysis",
                "QC assessment",
                "%s merged-VCF variants failed ≥1 of %s active rules"
                % (
                    format_metric_count(assessment.get("excluded_total")),
                    format_metric_count(assessment.get("active_rule_count")),
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
                % (
                    format_metric_count(len(completed)),
                    format_metric_count(chromosome_total),
                ),
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
                    format_metric_count(len(primary_outputs)),
                    str(final_status).upper(),
                ),
            ),
        ],
    )
    return lines
