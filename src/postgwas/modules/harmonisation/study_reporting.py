"""Present saved study evidence without reading inputs or repeating analyses.

This module consumes the decision dictionaries produced by preparation. It does
not infer a type, evaluate acceptance thresholds or modify scientific values.
Display percentages are scalar formatting of recorded counts/fractions only.
The 0.5 frequency label is the minor-allele probability boundary, not a new
configurable QC threshold.
"""

from __future__ import annotations

import json

from postgwas.core.ui import screen_line, screen_summary_card
from postgwas.core.values import (
    format_count, format_fraction_percentage, format_number, format_percentage,
)

from .shared.runtime import resolve_policies


def _text(value):
    return "not available" if value is None else str(value)


def _words(value):
    return _text(value).replace("_", " ")


def _input_column(decisions, key):
    if key not in decisions:
        return "Not recorded"
    return str(decisions[key]) if decisions[key] is not None else "Not supplied"


def _pair(part, whole):
    if part is None or whole is None or whole == 0:
        return "not available"
    return "%s / %s (%s)" % (
        format_count(part), format_count(whole), format_percentage(part, whole),
    )


def _source(value):
    return {
        "detector": "Automatically inferred",
        "policy": "YAML policy declaration",
        "sample_sheet": "Sample-sheet declaration",
        "dataset_z_pvalue_cross_check": "Dataset Z/P-value comparison",
    }.get(value, _words(value))


def _range(evidence):
    return "%s to %s; median %s" % tuple(
        format_number(evidence.get(key)) for key in ("min", "max", "median")
    )


def _numeric_fields(evidence, detector):
    """Keep count/range evidence on one recorded subset, never mixed casts.

    The canonical detectors record a finite numeric subset. Generic cast
    evidence can also include infinities/NaN, so it cannot supply a numerator
    or range for the detector's finite denominator. Exclusion counts below are
    scalar differences of already-recorded counts, not new input analyses.
    """
    finite = "n_usable" in detector
    usable = detector.get("n_usable") if finite else evidence.get("non_null")
    fields = [
        ("count", "Finite numeric values" if finite else "Numeric values before finite QC", format_count(usable, missing="not available")),
        ("info", "Observed range and median", _range(detector if finite else evidence)),
    ]
    non_numeric = evidence.get("unparseable")
    cast_null = evidence.get("null")
    if cast_null is not None and non_numeric is not None:
        fields.append(("count", "Missing input", format_count(cast_null - non_numeric)))
    if non_numeric is not None:
        fields.append(("count", "Non-numeric input", format_count(non_numeric)))
    if finite:
        cast_non_null = detector.get("n_non_null", evidence.get("non_null"))
        rows = detector.get("n_rows", evidence.get("rows"))
        if usable is not None and cast_non_null is not None:
            fields.append(("count", "Non-finite input", format_count(cast_non_null - usable)))
        if usable is not None and rows is not None:
            fields.append(("count", "Excluded from inference", format_count(rows - usable)))
    return fields, usable, finite


def _build_fields(info, references):
    selected = info.get("inferred_build")
    # "Ambiguous" is the established genome-build decision sentinel, never a
    # resource build name. An accompanying ambiguity reason also bars selection.
    ambiguous = selected == "Ambiguous" or bool(info.get("ambiguous_reason"))
    fields = [("genetic", "Selected build", "No build selected" if ambiguous else _text(selected))]
    if info.get("acceptance_path"):
        fields.append(("info", "Acceptance path", _words(info["acceptance_path"])))
    fallback = info.get("low_coverage_fallback") or {}
    if fallback:
        fields.append(("info", "Low-coverage fallback", _words(fallback.get("status"))))
        if fallback.get("accepted"):
            fields.append(("warning", "Experimental inference", "Dominant-build inference; match fraction is not a posterior probability and small mixed components can pass."))
        if fallback.get("failed_gates") and fallback.get("attempted"):
            fields.append(("warning", "Fallback failed gates", ", ".join(fallback["failed_gates"])))
    for build, count in ((info.get("distinct_loci") or {}).get("exclusive_snp_loci") or {}).items():
        fields.append(("count", "%s exclusive SNP loci" % build, format_count(count)))
    for chromosome, item in (info.get("per_chromosome") or {}).items():
        fields.append(("genetic", "Chromosome %s build evidence" % chromosome,
                       "%s; %s testable; matches %s; exclusive SNP loci %s" % (
                           _words(item.get("assessment")), format_count(item.get("testable_variants")),
                           _text(item.get("matches")), _text(item.get("exclusive_snp_loci")),
                       )))
    fields.extend(("count", label, format_count(info.get(key), missing="not available"))
                  for label, key in (
                      ("Input variants", "input_variants"),
                      ("Coordinate-testable", "testable_variants"),
                      ("Not build-testable", "untestable_variants"),
                  ))
    for build, count in (info.get("matches") or {}).items():
        fields.append(("genetic", "%s matches" % build,
                       "%s (%s%% input; %s%% testable)" % (
                           format_count(count),
                           format_number((info.get("input_percentages") or {}).get(build)),
                           format_number((info.get("percentages") or {}).get(build)),
                       )))
        if build in (info.get("reference_marker_counts") or {}):
            fields.append(("genetic", "%s reference coverage" % build,
                           "%s / %s markers (%s%%)" % (
                               format_count((info.get("matched_reference_marker_counts") or {}).get(build), missing="not available"),
                               format_count(info["reference_marker_counts"][build]),
                               format_number((info.get("reference_percentages") or {}).get(build)),
                           )))
    if selected in references:
        fields.append(("info", "Reference", str(references[selected])))
    reason = info.get("ambiguous_reason") or info.get("decision_reason")
    if reason:
        fields.append(("info", "Decision reason", str(reason)))
    if ambiguous:
        fields.append(("error", "Planned action", "Build-specific resources are not selected; resolve the recorded genome-build ambiguity."))
    elif selected is not None:
        fields.append(("info", "Planned action", "Use the selected %s source-build resources." % selected))
    return fields


def _effect_fields(decisions, policies):
    evidence = decisions.get("effect_evidence") or {}
    detector = evidence.get("detector_stats") or {}
    source = decisions.get("effect_type_source")
    effect = decisions.get("effect_type")
    fields = [
        ("info", "Input column", _input_column(decisions, "effect_column")),
        ("info", "Decision source", _source(source)),
    ]
    if evidence:
        numeric_fields, usable, finite = _numeric_fields(evidence, detector)
        fields.extend(numeric_fields)
        non_positive = detector.get("non_positive_value_count") if finite else evidence.get("non_positive")
        if finite and non_positive is None and all(key in detector for key in ("n_negative", "n_zero")):
            non_positive = detector["n_negative"] + detector["n_zero"]
        fields.append(("info", "Zero or negative", _pair(non_positive, usable)))
    if source == "detector" or decisions.get("effect_type_detected") is not None:
        low, high = policies.get("effect.or_detection_require_median_range")
        fields.extend([
            ("info", "OR non-positive limit", format_fraction_percentage(policies.get("effect.or_detection_max_non_positive_fraction"))),
            ("info", "OR median range", "%s to %s inclusive" % (format_number(low), format_number(high))),
        ])
    fields.append(("analysis", "Selected type", _words(effect)))
    if "effect_type_matches_declaration" in decisions:
        fields.append(("info", "Declaration cross-check", _agreement(decisions["effect_type_matches_declaration"])))
    if source == "detector":
        fields.append(("info", "Decision reason", "The recorded OR detection rule was satisfied." if effect == "odds_ratio" else "The recorded detector selected beta; the OR rule was not satisfied."))
        fields.append(("warning", "Interpretation warning", "Numeric patterns support this selection; they cannot establish the column's meaning. Check the study documentation."))
    if effect is not None and decisions.get("effect_column") is not None:
        fields.append(("info", "Planned action", "Convert using BETA = ln(OR)." if effect == "odds_ratio" else "Keep the effect on its supplied beta scale."))
    elif "effect_column" in decisions and decisions["effect_column"] is None:
        fields.append(("info", "Planned action", "No effect-size column supplied; later effect recovery remains pending."))
    return fields


def _agreement(value):
    return "Agrees" if value is True else "Disagrees; see recorded warning/policy action" if value is False else "Not independently established"


def _se_fields(decisions, policies):
    if ("se_column" in decisions and decisions["se_column"] is None) or decisions.get("se_scale_source") == "no_standard_error_column":
        return [("info", "Check status", "No standard-error column supplied; SE recovery remains pending under the existing policies.")]
    if decisions.get("effect_type") == "beta":
        return [("info", "Check status", "Not applicable to beta; beta-scale SE is unchanged.")]
    evidence = decisions.get("se_scale_evidence") or {}
    fields = [
        ("info", "Input column", decisions.get("se_column") or "Not recorded"),
        ("info", "Decision source", _source(decisions.get("se_scale_source"))),
    ]
    if evidence:
        fields.append(("info", "Comparison evidence", _words(evidence.get("comparison_source"))))
        informative = evidence.get("informative_variants")
        fields.append(("count", "Informative variants", format_count(informative, missing="not available")))
        if informative:
            fields.extend([
                ("info", "Log-odds agreement", _pair(evidence.get("log_odds_matches"), informative)),
                ("info", "Raw-OR agreement", _pair(evidence.get("raw_odds_ratio_matches"), informative)),
                ("info", "Exponentiated-SE agreement", _pair(evidence.get("exponentiated_log_odds_matches"), informative)),
            ])
        fields.extend([
            ("info", "Minimum informative count", format_count(evidence.get("minimum_informative_variants", policies.get("effect.se_scale_min_variants")))),
            ("info", "Required agreement", format_fraction_percentage(evidence.get("minimum_agreement_fraction", policies.get("effect.se_scale_min_agreement_fraction")))),
            ("info", "Required agreement margin", format_fraction_percentage(evidence.get("minimum_agreement_margin", policies.get("effect.se_scale_min_agreement_margin")))),
            ("info", "Comparison status", _words(evidence.get("status"))),
        ])
        if evidence.get("error"):
            fields.append(("warning", "Comparison unavailable", str(evidence["error"])))
    scale = decisions.get("se_scale")
    fields.append(("analysis", "Selected scale", scale or "Not decided"))
    if "se_scale_matches_declaration" in decisions:
        fields.append(("info", "Declaration cross-check", _agreement(decisions["se_scale_matches_declaration"])))
    if scale in ("log_odds", "as_given"):
        fields.append(("info", "Planned action", "Keep the supplied SE on the log-odds scale." if scale == "log_odds" else "Divide raw OR-scale SE by the original OR before allele orientation."))
    elif scale == "exponentiated_log_odds":
        fields.append(("info", "Planned action", "Decode SE_logOR = ln(stored_SE) in a separate working column; preserve the supplied values."))
    return fields


def _pvalue_fields(decisions, policies):
    evidence = decisions.get("pvalue_evidence") or {}
    detector = evidence.get("detector_stats") or {}
    source = decisions.get("pvalue_type_source")
    kind = decisions.get("pvalue_type")
    fields = [
        ("info", "Input column", _input_column(decisions, "pvalue_column")),
        ("info", "Decision source", _source(source)),
    ]
    if evidence:
        numeric_fields, usable, finite = _numeric_fields(evidence, detector)
        fields.extend(numeric_fields)
        negative = detector.get("n_negative") if finite else evidence.get("n_negative", evidence.get("negative"))
        fields.append(("info", "Negative values", _pair(negative, usable) if finite else format_count(negative, missing="not available")))
        fields.append(("info", "Negative fraction limit", format_fraction_percentage(evidence.get("max_negative_fraction", policies.get("pvalue.max_negative_fraction")))))
    if source == "detector" or decisions.get("pvalue_type_detected") is not None:
        finite = "n_usable" in detector
        threshold = detector.get("detect_threshold", evidence.get("detect_threshold", policies.get("pvalue.mlogp_detect_threshold")))
        fields.extend([
            ("info", "Above %s" % format_number(threshold), _pair(detector.get("n_above_threshold") if finite else evidence.get("above_threshold"), detector.get("n_usable") if finite else evidence.get("non_null"))),
            ("info", "Negative-log10 rule", "At least %s values AND at least %s must exceed %s; the configured median rule also applies when both gates pass." % (
                format_count(policies.get("pvalue.mlogp_detect_min_count")),
                format_fraction_percentage(policies.get("pvalue.mlogp_detect_proportion")),
                format_number(threshold),
            )),
            ("info", "Median detection rule", "Within %s of %s for negative-log10 detection; not used to select raw scale." % (
                format_number(detector.get("median_tolerance", policies.get("pvalue.mlogp_median_tolerance"))),
                format_number(detector.get("expected_median_neglog10", policies.get("pvalue.mlogp_expected_median"))),
            )),
        ])
    fields.append(("analysis", "Selected type", {"raw": "Raw p-values", "neglog10": "Negative-log10 p-values"}.get(kind, "Not decided")))
    if "pvalue_type_matches_declaration" in decisions:
        fields.append(("info", "Declaration cross-check", _agreement(decisions["pvalue_type_matches_declaration"])))
    if source == "detector":
        fields.append(("info", "Decision reason", "The negative-log10 count and fraction gates were not both met." if kind == "raw" else "The recorded negative-log10 count, fraction and median requirements were met."))
        fields.append(("warning", "Interpretation warning", "Numeric patterns cannot prove column meaning or one-sided versus two-sided testing. Check the study documentation."))
    if kind is not None and decisions.get("pvalue_column") is not None:
        fields.append(("info", "Planned action", "Keep raw probabilities; chromosome P-value QC still applies." if kind == "raw" else "Convert negative-log10 values to canonical P; chromosome P-value QC still applies."))
    elif "pvalue_column" in decisions and decisions["pvalue_column"] is None:
        fields.append(("info", "Planned action", "No p-value column supplied; later p-value recovery remains pending."))
    return fields


def _beta_z_fields(decisions):
    evidence = decisions.get("beta_z_sign_consistency") or {}
    if evidence.get("applies") is False:
        return [("info", "Check status", "Not applicable — %s." % evidence.get("reason", "required inputs are absent"))]
    comparable = evidence.get("comparable_missing_se_rows")
    if comparable == 0:
        return [("info", "Check status", "Not assessed — no comparable rows needing SE recovery.")]
    if comparable is None:
        return [("info", "Check status", "Not recorded")]
    return [
        ("count", "Comparable missing-SE rows", format_count(comparable)),
        ("info", "Incompatible directions", _pair(evidence.get("discordant_rows"), comparable)),
        ("info", "Maximum allowed fraction", format_fraction_percentage(evidence.get("maximum_fraction"))),
        ("info", "Decision", _words(evidence.get("decision"))),
        ("info", "Affected-row action", _words(evidence.get("row_action"))),
    ]


def replacement_frequency_fields(decisions):
    """Present final orchestration evidence, never re-evaluate numeric gates.

    A numeric ``passed`` alone is not an adoption decision: preparation issues
    its row/source-bound contract only after the final identity checks succeed.
    """
    if decisions.get("eaf_selected_source") != "replacement_reference":
        return []
    evidence = decisions.get("eaf_replacement_validation") or {}
    source = evidence.get("source") or decisions.get("eaf_replacement") or {}
    status = evidence.get("status")
    accepted = status == "passed" and bool(evidence.get("contract"))
    outcome = {
        "inconclusive": "Inconclusive; dataset stopped before chromosome analysis",
        "invalid": "Invalid; dataset stopped before chromosome analysis",
    }.get(status)
    if outcome is None:
        outcome = ("Passed; replacement EAF accepted for assignment" if accepted else
                   "Numeric gates passed; replacement approval not recorded" if status == "passed" else
                   "Not validated; dataset stopped before replacement approval" if decisions.get("report_failure") else
                   "Provisionally selected; replacement validation pending; assignment pending.")
    fields = [
        ("analysis" if accepted else "error" if status in ("inconclusive", "invalid") else "info", "Replacement validation", outcome),
        ("info", "Replacement source", "%s; column %s" % (_text(source.get("file")), _text(source.get("column")))),
        ("info", "Original study AF preserved", "Original study MAF remains unchanged; working EAF is a population-reference proxy."),
    ]
    if not evidence:
        return fields
    fields.extend([
        ("count", "Study rows entering AF preflight (fixed N)", format_count(evidence.get("study_variant_count"), missing="not available")),
        ("count", "Replacement paired input rows", format_count(evidence.get("input_pair_rows"), missing="not available")),
        ("count", "Replacement comparable variants", "%s; required >= %s" % (format_count(evidence.get("comparable_variants"), missing="not available"), format_count(evidence.get("minimum_overlap"), missing="not available"))),
        ("info", "Replacement diagnostic coverage", "%s; required >= %s" % (_pair(evidence.get("comparable_variants"), evidence.get("study_variant_count")), format_fraction_percentage(evidence.get("minimum_overlap_fraction")))),
        ("info", "Replacement correlation method", _text(evidence.get("correlation_method"))),
        ("info", "Replacement correlation", "%s; required >= %s" % (format_number(evidence.get("correlation")), format_number(evidence.get("minimum_correlation")))),
        ("info", "Replacement mean absolute difference", "%s; required <= %s" % (format_number(evidence.get("mean_absolute_difference")), format_number(evidence.get("maximum_mean_absolute_difference")))),
        ("info", "Replacement major frequency (B > 0.5)", "%s; required > %s" % (_pair(evidence.get("major_count"), evidence.get("comparable_variants")), format_fraction_percentage(evidence.get("minimum_major_fraction")))),
        ("info", "Replacement shared evidence gate", "Passed" if evidence.get("shared_evidence_gate_passed") is True else "Failed" if evidence.get("shared_evidence_gate_passed") is False else "Not recorded"),
        ("info", "Replacement failed requirements", "; ".join(_words(item) for item in evidence.get("failed_requirements") or []) or "None recorded"),
        ("info", "Replacement decision reason", _words(evidence.get("failure_reason", evidence.get("decision_reason")))),
    ])
    fields.extend(("count", "Replacement excluded: %s" % _words(key), format_count(value, missing="not available"))
                  for key, value in (evidence.get("excluded_counts") or {}).items())
    for key, label in (("mapping", "Replacement declared allele mapping"),
                       ("resolved_files", "Replacement resolved source files"),
                       ("identities", "Replacement source identities")):
        if key in source:
            fields.append(("info", label, json.dumps(source[key], sort_keys=True, ensure_ascii=False)))
    reused = source.get("same_source_evidence_reused")
    fields.extend([
        ("info", "Reused default-reference evidence", "Yes" if reused is True else "No" if reused is False else "Not recorded"),
        ("info", "Chromosomes reusing default evidence", ", ".join(map(str, source.get("reused_default_chromosomes") or [])) or "None recorded"),
        ("info", "Replacement interpretation", "Folded-frequency and major-profile agreement supports consistency under the declared allele mapping; it does not prove allele direction or measured study EAF. Reused source evidence is not independent evidence."),
    ])
    return fields


def render_replacement_frequency_validation(study_decisions):
    """Render the later preparation outcome without repeating initial decisions."""
    fields = replacement_frequency_fields(study_decisions or {})
    return "\n".join(screen_summary_card("genetic", "Whole-study replacement AF validation", fields)) if fields else ""


def _frequency_fields(decisions, policies, resources):
    evidence = decisions.get("eaf_evidence") or {}
    comparison = decisions.get("eaf_reference_validation") or {}
    status = comparison.get("status")
    source = decisions.get("eaf_selected_source")
    column = decisions.get("eaf_column")
    if column is None and decisions.get("eaf_is_maf_source") == "no_study_frequency_column":
        return [
            ("info", "Input frequency", "No study frequency column supplied"),
            ("info", "Selected source", "External EAF selected; per-chromosome validation and assignment pending."),
        ]
    fields = [("info", "Input column", column or "Not recorded")]
    if evidence:
        fields.extend([
            ("count", "Usable frequencies", format_count(evidence.get("usable"), missing="not available")),
            ("count", "Missing or invalid", format_count(evidence.get("unusable"), missing="not available")),
            ("info", "Observed range and median", _range(evidence)),
            ("info", "At or below 0.5", _pair(evidence.get("at_or_below_0_5"), evidence.get("usable"))),
            ("info", "MAF-like screen cutoff", "strictly above %s" % format_fraction_percentage(decisions.get("eaf_maf_decision_cutoff", policies.get("eaf.maf_decision_cutoff")))),
        ])
        fields.extend(("count", label, format_count(evidence[key], missing="not available"))
                      for label, key in (("Missing input", "input_missing"),
                                         ("Non-numeric input", "unparseable"),
                                         ("Non-finite input", "non_finite"),
                                         ("Out-of-range input", "out_of_range"))
                      if key in evidence)
    initial = decisions.get("eaf_is_maf_initial", decisions.get("eaf_is_maf"))
    fields.append(("info", "Distribution screen", "MAF-like; reference confirmation required" if initial is True else "Not MAF-like" if initial is False else "Not recorded"))
    fields.append(("info", "Reference comparison", {
        "not_required": "Not required — distribution screen did not trigger comparison",
        "passed": "Allele-aligned reference agreement passed",
        "replacement_selected": "Input MAF supported; replacement source selected provisionally",
        "failed": "FAILED",
        "not_applicable": "Not applicable",
    }.get(status, "Pending")))
    if comparison.get("comparable_variants") is not None:
        references = dict.fromkeys(
            (item.get("default_eaf_file"), item.get("default_eaf_reference_column"))
            for item in resources.values() if item.get("default_eaf_file")
        )
        for path, population in references:
            fields.append(("info", "Comparison reference", "%s; population %s" % (path, _text(population))))
        minor_fraction_text = format_fraction_percentage(comparison.get("reference_effect_allele_minor_fraction"))
        if comparison.get("minimum_reference_effect_allele_minor_fraction") is not None:
            minor_fraction_text += "; required > %s" % format_fraction_percentage(
                comparison["minimum_reference_effect_allele_minor_fraction"]
            )
        fields.extend([
            ("count", "Comparable variants", "%s; required >= %s" % (format_count(comparison["comparable_variants"]), format_count(comparison.get("minimum_overlap"), missing="not available"))),
            ("info", "Correlation method", _text(comparison.get("correlation_method", comparison.get("maf_correlation_method")))),
            ("info", "Direct EAF correlation", "%s; required >= %s" % (format_number(comparison.get("direct_eaf_correlation", comparison.get("correlation"))), format_number(comparison.get("minimum_study_correlation", comparison.get("minimum_correlation", comparison.get("minimum_maf_correlation")))))),
            ("info", "Mean absolute EAF difference", "%s; required <= %s" % (format_number(comparison.get("mean_eaf_error", comparison.get("mean_absolute_difference"))), format_number(comparison.get("maximum_mean_absolute_difference", comparison.get("maximum_mean_absolute_maf_difference"))))),
            ("info", "Reference EA at or below 0.5", minor_fraction_text),
        ])
        if comparison.get("study_variant_count") is not None:
            fields.extend([
                ("info", "Usable diagnostic coverage", "%s; required >= %s" % (
                    _pair(comparison["comparable_variants"], comparison["study_variant_count"]),
                    format_fraction_percentage(comparison.get("minimum_overlap_fraction")),
                )),
            ])
        if "rule_1_passed" in comparison:
            fields.extend([
                ("info", "Shared evidence gate", "Passed" if comparison["shared_evidence_gate_passed"] else "Failed"),
                ("info", "Both rules evaluated", "Yes" if comparison["whole_study_rules_evaluated"] else "No — insufficient shared evidence"),
                ("info", "Supplied AF versus reference MAF correlation", "%s; required >= %s" % (
                    format_number(comparison.get("supplied_to_reference_maf_correlation")),
                    format_number(comparison.get("minimum_study_correlation")),
                )),
                ("info", "Mean absolute supplied-AF / reference-MAF difference", "%s; required <= %s" % (
                    format_number(comparison.get("mean_maf_error")),
                    format_number(comparison.get("maximum_mean_absolute_maf_difference")),
                )),
                ("info", "Reference EA above 0.5", "%s; required > %s" % (
                    format_fraction_percentage(comparison.get("reference_effect_allele_major_fraction")),
                    format_fraction_percentage(comparison.get("minimum_reference_effect_allele_major_fraction")),
                )),
            ])
            for rule, label in ((1, "Rule 1 (retain original AF)"), (2, "Rule 2 (select replacement)")):
                fields.append(("info", label, "Passed" if comparison["rule_%s_passed" % rule] else "Failed"))
                failures = comparison.get("rule_%s_failed_requirements" % rule) or []
                fields.append(("info", "Rule %s failed requirements" % rule,
                               "; ".join(_words(item) for item in failures) if failures else "None"))
            fields.append(("info", "Combined AF decision", _words(comparison.get("decision_reason"))))
    if comparison.get("maf_reference_decision") == "maf":
        fields.append(("info", "Input interpretation", "Input MAF supported by reference comparison; the supplied column is not accepted as EAF."))
    if status == "failed":
        fields.extend([
            ("error", "Frequency decision", "FAILED — %s" % _words(comparison.get("failure_reason", comparison.get("decision_reason")))),
            ("error", "Planned action", "No chromosome analysis; no frequency replacement or EAF acceptance was authorised."),
        ])
    elif status == "replacement_selected" and source == "replacement_reference":
        fields.extend([
            ("analysis", "Frequency decision", "Replacement EAF selected provisionally by original AF rules."),
        ])
        fields.extend(replacement_frequency_fields(decisions))
    elif status in ("passed", "not_required"):
        fields.extend([
            ("analysis", "Frequency decision", "Study EAF retained" if status == "passed" else "Declared study EAF retained"),
            ("info", "Planned action", "Keep the supplied EAF; normal chromosome QC still applies."),
        ])
    return fields


def _strand_fields(decisions, policies):
    evidence = {**(decisions.get("strand_consensus") or {}), **decisions}
    mode = evidence.get("strand_mode")
    strand = evidence.get("strand")
    if mode == "reference_aligned" or strand == "reference_aligned":
        return [
            ("info", "Decision source", "Explicit reference-aligned declaration"),
            ("genetic", "Selected strand mode", "reference_aligned"),
            ("info", "Planned action", "Validate direct/swapped REF/ALT order per variant; no strand consensus was inferred."),
        ]
    fields = [("info", "Decision source", "Inferred from selected-build allele matches" if mode == "auto" else _text(mode))]
    for label, key in (("Forward matches", "forward"), ("Reverse matches", "reverse"),
                       ("Ambiguous matches excluded", "ambiguous"), ("Informative matches", "informative")):
        fields.append(("count", label, format_count(evidence.get(key), missing="not available")))
    fields.extend([
        ("info", "Minimum informative count", format_count(evidence.get("minimum_informative_variants", policies.get("strand.min_informative_variants")))),
        ("info", "Dominant fraction", format_fraction_percentage(evidence.get("dominant_fraction"))),
        ("info", "Required dominant fraction", format_fraction_percentage(evidence.get("consensus_threshold", policies.get("strand.consensus_threshold")))),
        ("genetic", "Selected consensus", strand or "Not decided"),
        ("info", "Decision reason", _text(evidence.get("reason"))),
    ])
    if strand is not None:
        fields.append(("info", "Planned action", "Use the recorded consensus during per-variant reference orientation; unresolved palindromes still require decisive evidence."))
    return fields


def render_study_decisions(
    study_decisions, *, genome_build_info=None, build_reference_files=None,
    policies=None, resource_maps=None,
) -> str:
    """Render only recorded evidence, keeping pending and failed states explicit.

    Optional context contains already-resolved metadata. Reference identifiers
    may be existing ``File N`` labels or literal full paths; they are never
    opened, inspected or reduced to potentially ambiguous basenames here.
    """
    decisions = study_decisions or {}
    policy = resolve_policies(policies)
    failure = decisions.get("report_failure")
    sections = []
    if genome_build_info:
        sections.append(("genetic", "Genome build", _build_fields(genome_build_info, build_reference_files or {})))
    sections.extend([
        ("analysis", "Effect type", _effect_fields(decisions, policy)),
        ("analysis", "Standard-error scale", _se_fields(decisions, policy)),
        ("analysis", "P-value type", _pvalue_fields(decisions, policy)),
        ("analysis", "BETA/Z direction for missing-SE recovery", _beta_z_fields(decisions)),
        ("genetic", "Allele frequency", _frequency_fields(decisions, policy, resource_maps or {})),
        ("genetic", "Strand consensus", _strand_fields(decisions, policy)),
    ])
    lines = [screen_line("error" if failure else "decision", "Study-wide decisions and evidence", indent=4)]
    for kind, title, fields in sections:
        lines.extend(screen_summary_card(kind, title, fields))
    if failure:
        lines.extend(screen_summary_card("error", "Dataset stopped", [("error", "Reason", str(failure))]))
    return "\n".join(lines)


__all__ = ["render_study_decisions", "render_replacement_frequency_validation", "replacement_frequency_fields"]
