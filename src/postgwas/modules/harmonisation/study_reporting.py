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
        if effect is not None:
            fields.append(("info", "Recommendation", "Check the study documentation. If confirmed, set effect_type = %s in the sample sheet, or effect.type in YAML." % effect))
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
        if kind is not None:
            fields.append(("info", "Recommendation", "Check the study documentation. If confirmed, set p_value_type = %s in the sample sheet, or pvalue.type in YAML." % kind))
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


def _effective_frequency_screen_fields(screen, label, *, reason=None):
    """Format the saved whole-study screen without reclassifying its values."""
    if not screen:
        return []
    result = ("No usable frequencies; replacement rejected" if reason == "no_usable_replacement_frequencies" else
        "MAF-like; replacement rejected" if reason == "replacement_frequency_maf_like" else
        "No usable frequencies" if screen.get("usable") == 0 else
        "MAF-like" if screen.get("maf_like") is True else
        "Not MAF-like" if screen.get("maf_like") is False else "Not recorded")
    return [
        ("info", label + " screen scope", {
            "merged_gwas_rows": "All merged GWAS rows; one contribution per row, pooled across chromosomes",
            "gwas_rows": "All GWAS rows; one contribution per row, pooled across chromosomes",
        }.get(screen.get("scope"), _words(screen.get("scope")))),
        ("info", label + " allele direction", {
            "original_gwas_effect_allele": "Frequency of the original listed GWAS effect allele",
        }.get(screen.get("allele_direction"), _words(screen.get("allele_direction")))),
        ("count", label + " screen rows", format_count(screen.get("rows"), missing="not available")),
        ("count", label + " usable frequencies", format_count(screen.get("usable"), missing="not available")),
        ("count", label + " missing or invalid frequencies", format_count(screen.get("unusable"), missing="not available")),
        ("info", label + " AF at or below 0.5", _pair(screen.get("at_or_below_0_5"), screen.get("usable"))),
        ("info", label + " MAF-like cutoff", "strictly above %s" % format_fraction_percentage(screen.get("decision_cutoff"))),
        ("info", label + " screen result", result),
    ]


def native_frequency_screen_fields(decisions, *, policies=None):
    """Share saved native-A screen presentation with terminal and HTML reports."""
    if decisions.get("eaf_column") is None:
        return []
    screen = decisions.get("eaf_effective_frequency_screen")
    if not screen and decisions.get("eaf_evidence"):
        screen = dict(decisions["eaf_evidence"])
        screen.setdefault("decision_cutoff", decisions.get("eaf_maf_decision_cutoff",
            resolve_policies(policies).get("eaf.maf_decision_cutoff")))
        screen.setdefault("maf_like", decisions.get("eaf_is_maf_initial", decisions.get("eaf_is_maf")))
        screen.setdefault("scope", "gwas_rows")
        screen.setdefault("allele_direction", "original_gwas_effect_allele")
    return _effective_frequency_screen_fields(screen, "Initial selected AF A")


def _frequency_source_fields(source, label, *, include_column=True):
    """Keep declarations readable; complete file/identity maps remain in evidence."""
    declared = tuple(dict.fromkeys(path for path in (source.get("declared_files") or {}).values() if path))
    resolved = tuple(dict.fromkeys(path for path in (source.get("resolved_files") or {}).values() if path))
    declaration = source.get("file") or (declared[0] if len(declared) == 1 else
        "%s declared files" % format_count(len(declared)) if declared else None)
    fields = [("info", label + " source", "%s; column %s" % (_text(declaration), _text(source.get("column")))
        if include_column else _text(declaration))]
    if declared and source.get("file"):
        fields.append(("info", label + " declared source files", "%s declared files" % format_count(len(declared))))
    if resolved:
        fields.append(("info", label + " resolved source files", "%s resolved files; full identities retained in structured evidence" % format_count(len(resolved))))
    if source.get("mapping"):
        fields.append(("info", label + " declared allele mapping", json.dumps(source["mapping"], sort_keys=True, ensure_ascii=False)))
    return fields


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
        *_frequency_source_fields(source, "Replacement"),
        ("info", "Original study AF preserved" if decisions.get("eaf_column") is not None else "Selected input AF A preserved",
            "Original study MAF remains unchanged; working EAF is a population-reference proxy."
            if decisions.get("eaf_column") is not None else
            "Initially selected external AF A remains unchanged; no native GWAS frequency column mapped. Working EAF B is a population-reference proxy."),
    ]
    if not evidence:
        return fields
    reason = evidence.get("failure_reason", evidence.get("decision_reason"))
    fields.extend(_effective_frequency_screen_fields(evidence.get("effective_frequency_screen"), "Replacement B", reason=reason))
    fields.append(("info", "Replacement decision reason", _words(reason)))
    if reason in ("replacement_frequency_maf_like", "no_usable_replacement_frequencies"):
        fields.append(("error", "Replacement stop action", "No further default comparison, replacement or fallback; dataset stopped before chromosome analysis."))
    fields.append(
        ("count", "Study rows entering AF preflight (fixed N)", format_count(evidence.get("study_variant_count"), missing="not available")))
    if "comparable_variants" in evidence:
        fields.extend([
        ("count", "Replacement paired input rows", format_count(evidence.get("input_pair_rows"), missing="not available")),
        ("count", "Replacement comparable variants", "%s; required >= %s" % (format_count(evidence.get("comparable_variants"), missing="not available"), format_count(evidence.get("minimum_overlap"), missing="not available"))),
        ("info", "Replacement diagnostic coverage", "%s; required >= %s" % (_pair(evidence.get("comparable_variants"), evidence.get("study_variant_count")), format_fraction_percentage(evidence.get("minimum_overlap_fraction")))),
        ("info", "Replacement correlation method", _text(evidence.get("correlation_method"))),
        ("info", "Replacement correlation", "%s; required >= %s" % (format_number(evidence.get("correlation")), format_number(evidence.get("minimum_correlation")))),
        ("info", "Replacement mean absolute difference", "%s; required <= %s" % (format_number(evidence.get("mean_absolute_difference")), format_number(evidence.get("maximum_mean_absolute_difference")))),
        ("info", "Replacement major frequency (B > 0.5)", "%s; required > %s" % (_pair(evidence.get("major_count"), evidence.get("comparable_variants")), format_fraction_percentage(evidence.get("minimum_major_fraction")))),
        ("info", "Replacement shared evidence gate", "Passed" if evidence.get("shared_evidence_gate_passed") is True else "Failed" if evidence.get("shared_evidence_gate_passed") is False else "Not recorded"),
        ("info", "Replacement failed requirements", "; ".join(_words(item) for item in evidence.get("failed_requirements") or []) or "None recorded"),
        ])
    fields.extend(("count", "Replacement excluded: %s" % _words(key), format_count(value, missing="not available"))
                  for key, value in (evidence.get("excluded_counts") or {}).items())
    reused = source.get("same_source_evidence_reused")
    fields.extend([
        ("info", "Reused default-reference evidence", "Yes" if reused is True else "No" if reused is False else "Not recorded"),
        ("count", "Chromosomes reusing default evidence", format_count(len(source.get("reused_default_chromosomes") or []))),
        ("info", "Replacement interpretation", "Folded-frequency and major-profile agreement supports consistency under the declared allele mapping; it does not prove allele direction or measured study EAF. Reused source evidence is not independent evidence."),
    ])
    return fields


def render_replacement_frequency_validation(study_decisions):
    """Render the later preparation outcome without repeating initial decisions."""
    fields = frequency_summary_fields(study_decisions) if study_decisions and study_decisions.get("eaf_selected_source") == "replacement_reference" else []
    return "\n".join(screen_summary_card("genetic", "Whole-study replacement AF validation", fields)) if fields else ""


def external_frequency_fields(decisions):
    """Present initial external-source evidence without calling it study AF."""
    if not decisions.get("eaf_external_validation_required"):
        return []
    evidence = decisions.get("eaf_external_validation") or {}
    status = evidence.get("status")
    accepted = status == "passed" and bool(evidence.get("contract")) and decisions.get("eaf_selected_source") != "replacement_reference"
    outcome = ("Passed; external EAF accepted for assignment" if accepted else
        "Numeric gates passed; external approval not recorded" if status == "passed" else
        "Initial external A has MAF support; explicit replacement selected provisionally" if status == "replacement_selected" else
        "Inconclusive; dataset stopped before chromosome analysis" if status == "inconclusive" else
        "Invalid; dataset stopped before chromosome analysis" if status == "invalid" else
        "Not validated; dataset stopped before external approval" if decisions.get("report_failure") else
        "External EAF selected; whole-study validation and assignment pending.")
    fields = [
        ("analysis" if accepted else "error" if status in ("inconclusive", "invalid") else "info", "External validation", outcome),
        ("info", "Frequency provenance", "Supplied external frequency; no native GWAS frequency column mapped. Population-reference proxy."),
    ]
    if not evidence:
        return fields
    screen = evidence.get("raw_frequency_screen") or {}
    fields.extend([
        ("count", "Study rows entering AF preflight (fixed N)", format_count(evidence.get("study_variant_count"), missing="not available")),
        ("info", "External decision reason", _words(evidence.get("failure_reason", evidence.get("decision_reason")))),
    ])
    if screen:
        fields.extend([
            ("info", "Source-panel frequencies (descriptive only)", "Raw source-panel values do not decide the GWAS-row screen."),
            ("info", "Raw external AF at or below 0.5", _pair(screen.get("at_or_below_0_5"), screen.get("usable"))),
        ])
    fields.extend(_effective_frequency_screen_fields(evidence.get("effective_frequency_screen"), "Initial selected AF A"))
    if evidence.get("decision_reason") == "external_frequency_not_predominantly_minor":
        fields.append(("info", "Default semantic comparison", "Skipped; selected A was not MAF-like. Normal strand and chromosome QC still apply."))
    partial_screens = evidence.get("raw_frequency_screens_by_chromosome") or {}
    if partial_screens and "raw_frequency_screen" not in evidence:
        fields.append(("info", "Partial source-panel description",
            "%s chromosome(s) described; these panel counts do not decide the GWAS-row screen." % len(partial_screens)))
    if "comparable_variants" in evidence:
        fields.extend([
            ("count", "External comparable variants", "%s; required >= %s" % (format_count(evidence.get("comparable_variants")), format_count(evidence.get("minimum_overlap")))),
            ("info", "External diagnostic coverage", "%s; required >= %s" % (_pair(evidence.get("comparable_variants"), evidence.get("study_variant_count")), format_fraction_percentage(evidence.get("minimum_overlap_fraction")))),
            ("info", "External correlation method", _text(evidence.get("maf_correlation_method"))),
            ("info", "External shared evidence gate", "Passed" if evidence.get("shared_evidence_gate_passed") else "Failed"),
            ("info", "External direct correlation", "%s; required >= %s" % (format_number(evidence.get("direct_eaf_correlation")), format_number(evidence.get("minimum_study_correlation")))),
            ("info", "External direct mean absolute difference", "%s; required <= %s" % (format_number(evidence.get("mean_eaf_error")), format_number(evidence.get("maximum_mean_absolute_maf_difference")))),
            ("info", "Reference frequency of original GWAS effect allele <= 0.5", "%s; required > %s" % (format_fraction_percentage(evidence.get("reference_effect_allele_minor_fraction")), format_fraction_percentage(evidence.get("minimum_reference_effect_allele_minor_fraction")))),
            ("info", "External AF versus reference MAF correlation", format_number(evidence.get("supplied_to_reference_maf_correlation"))),
            ("info", "External AF versus reference MAF mean absolute difference", format_number(evidence.get("mean_maf_error"))),
            ("info", "Reference frequency of original GWAS effect allele > 0.5", "%s; required > %s" % (format_fraction_percentage(evidence.get("reference_effect_allele_major_fraction")), format_fraction_percentage(evidence.get("minimum_reference_effect_allele_major_fraction")))),
        ])
        for rule, label in ((1, "Retain declared external EAF"), (2, "External MAF support; select explicit replacement")):
            evaluated = "rule_%s_passed" % rule in evidence
            fields.append(("info", "External rule %s: %s" % (rule, label),
                _text(evidence["rule_%s_passed" % rule]) if evaluated else "Not evaluated"))
            fields.append(("info", "External rule %s failed requirements" % rule,
                ("; ".join(_words(item) for item in evidence.get("rule_%s_failed_requirements" % rule) or []) or "None")
                if evaluated else "Not evaluated"))
    for key, label in (("estimated_rank_memory_bytes", "Spearman estimated rank memory (bytes)"),
                       ("available_rank_memory_bytes", "Spearman available rank memory (bytes)")):
        if key in evidence:
            fields.append(("info", label, format_number(evidence[key])))
    if evidence.get("failure_message"):
        fields.append(("error", "External failure details", evidence["failure_message"]))
    fields.extend(("count", "External excluded: %s" % _words(key), format_count(value))
        for key, value in (evidence.get("excluded_counts") or {}).items())
    fields.extend(_frequency_source_fields(evidence.get("source") or decisions.get("eaf_external_source") or {}, "External"))
    fields.append(("info", "External interpretation", "Configured QC heuristics support this source decision; they do not prove column semantics or measured study EAF. Normal chromosome QC and output concordance still apply."))
    return fields


def render_external_frequency_validation(study_decisions):
    fields = frequency_summary_fields(study_decisions) if study_decisions and study_decisions.get("eaf_external_validation_required") else []
    return "\n".join(screen_summary_card("genetic", "Whole-study external AF validation", fields)) if fields else ""


def _frequency_approval(evidence, label, failure=False):
    """Describe recorded approval; numeric success alone grants no permission."""
    status = evidence.get("status")
    if status == "passed":
        return (label + " approved" if evidence.get("contract") else
                "Numeric gates passed; " + label.lower() + " approval not recorded")
    if status in ("failed", "invalid", "inconclusive"):
        return _words(status).capitalize() + "; dataset stopped before chromosome analysis"
    return ("Not validated; dataset stopped before %s approval" % label.lower() if failure else
            label + " validation pending; assignment pending")


def _frequency_comparison_fields(evidence):
    """Summarise recorded semantic comparisons without evaluating their rules."""
    if evidence.get("comparable_variants") is None:
        return []
    fields = [("count", "Comparable variants", _pair(evidence.get("comparable_variants"),
        evidence.get("study_variant_count")) if evidence.get("study_variant_count") else
        format_count(evidence["comparable_variants"]))]
    for rule, label, correlation, error in (
        (1, "Agreement with EAF", "direct_eaf_correlation", "mean_eaf_error"),
        (2, "Agreement with MAF", "supplied_to_reference_maf_correlation", "mean_maf_error"),
    ):
        passed = evidence.get("rule_%s_passed" % rule)
        if rule == 1 and passed is None:
            passed = evidence.get("direct_eaf_accepted")
        outcome = ("Not evaluated" if evidence.get("whole_study_rules_evaluated") is False else
            "Passed" if passed is True else "Failed" if passed is False else "Outcome not recorded")
        fields.append(("info", label, "%s; correlation %s; mean absolute difference %s" % (
            outcome, format_number(evidence.get(correlation)), format_number(evidence.get(error)))))
    return fields


def _frequency_failure_fields(evidence, label):
    fields = []
    for key in ("failure_reason", "decision_reason"):
        if evidence.get(key) and (key == "failure_reason" or not evidence.get("failure_reason")):
            if evidence.get("status") in ("failed", "invalid", "inconclusive"):
                fields.append(("error", label + " stop reason", _words(evidence[key])))
    failures = list(dict.fromkeys([item for key in (
        "failed_requirements", "shared_evidence_failed_requirements",
        "rule_1_failed_requirements", "rule_2_failed_requirements")
        for item in evidence.get(key) or []]))
    if failures and evidence.get("status") in ("failed", "invalid", "inconclusive"):
        fields.append(("error", label + " failed requirements", "; ".join(_words(item) for item in failures)))
    if evidence.get("status") in ("failed", "invalid", "inconclusive"):
        limits = {**evidence,
            "minimum_correlation": evidence.get("minimum_study_correlation", evidence.get("minimum_correlation")),
            "maximum_mean_absolute_difference": evidence.get("maximum_mean_absolute_difference",
                evidence.get("maximum_mean_absolute_maf_difference"))}
        for key, name in (("minimum_overlap", "Minimum comparable variants"),
            ("minimum_overlap_fraction", "Minimum study coverage"),
            ("minimum_correlation", "Minimum correlation"),
            ("maximum_mean_absolute_difference", "Maximum mean absolute difference"),
            ("minimum_major_fraction", "Minimum major-allele fraction"),
            ("estimated_rank_memory_bytes", "Estimated rank memory (bytes)"),
            ("available_rank_memory_bytes", "Available rank memory (bytes)")):
            if limits.get(key) is not None:
                fields.append(("info", label + " " + name.lower(), format_count(limits[key]) if key == "minimum_overlap" else format_number(limits[key])))
    if evidence.get("failure_message"):
        fields.append(("error", label + " failure details", evidence["failure_message"]))
    return fields


def _input_frequency_summary_fields(decisions, initial, external, replacement):
    """Describe input quality and the saved screen in original effect-allele direction."""
    source = initial.get("source") or decisions.get("eaf_external_source") or {}
    fields = [("info", "Input source", "Supplied external frequency; no native GWAS frequency column mapped; %s" %
        _frequency_source_fields(source, "Input")[0][2] if external else "GWAS column %s" % _text(decisions.get("eaf_column")))]
    screen = (initial.get("effective_frequency_screen") if external else
        decisions.get("eaf_effective_frequency_screen") or decisions.get("eaf_evidence")) or {}
    if screen:
        fields.extend([
            ("info", "Screen scope", "Whole-study merged GWAS rows, original listed GWAS effect allele"),
            ("count", "Variants assessed", format_count(screen.get("rows"), missing="not available")),
            ("info", "Input values at or below 0.5", _pair(screen.get("at_or_below_0_5"), screen.get("usable"))),
            ("info", "Input distribution", "No usable frequencies" if screen.get("usable") == 0 else
                ("MAF-like" if initial.get("status") in ("passed", "replacement_selected") else
                    "MAF-like; reference confirmation required") if screen.get("maf_like",
                    decisions.get("eaf_is_maf_initial", decisions.get("eaf_is_maf"))) is True else
                "Not MAF-like" if screen.get("maf_like", decisions.get("eaf_is_maf_initial",
                    decisions.get("eaf_is_maf"))) is False else "Not recorded"),
        ])
        if screen.get("unusable") or screen.get("usable") == 0:
            fields.extend([
                ("count", "Input usable frequencies", format_count(screen.get("usable"), missing="not available")),
                ("count", "Input missing or invalid frequencies", format_count(screen.get("unusable"), missing="not available")),
            ])
        cutoff = screen.get("decision_cutoff", decisions.get("eaf_maf_decision_cutoff"))
        if cutoff is not None:
            fields = [(kind, label, value + "; cutoff strictly above %s" % format_fraction_percentage(cutoff)
                if label == "Input distribution" else value) for kind, label, value in fields]
    quality = decisions.get("eaf_evidence") or {}
    quality_counts = [(name, quality.get(key)) for name, key in (
        ("Missing input", "input_missing"), ("Non-numeric input", "unparseable"),
        ("Non-finite input", "non_finite"), ("Out-of-range input", "out_of_range"))]
    if any(count for _, count in quality_counts):
        fields.append(("info", "Input quality", "; ".join("%s: %s" % (name, format_count(count))
            for name, count in quality_counts if count is not None)))
    if replacement or decisions.get("eaf_input_interpretation") == "maf":
        fields.append(("info", "Input interpretation", "Input MAF supported by reference comparison"))
        fields.append(("info", "Original column", "Preserved unchanged; working EAF is a population-reference proxy"))
    return fields


def _replacement_summary_fields(decisions, validation):
    """Describe replacement preparation independently of final chromosome losses."""
    fields = []
    selected = validation.get("source") or decisions.get("eaf_replacement") or {}
    fields.extend([
        ("info", "Replacement source", _frequency_source_fields(selected, "Replacement", include_column=False)[0][2]),
        ("info", "Frequency column", _text(selected.get("column"))),
    ])
    screen = validation.get("effective_frequency_screen") or {}
    if screen:
        fields.extend([
            ("info", "Replacement distribution", "%s; %s of usable values at or below 0.5" % (
                "No usable frequencies; replacement rejected" if screen.get("usable") == 0 else
                "MAF-like; replacement rejected" if screen.get("maf_like") is True else
                "Not MAF-like" if screen.get("maf_like") is False else "Not recorded",
                _pair(screen.get("at_or_below_0_5"), screen.get("usable")))),
            ("count", "Unresolved at preparation", "%s rows; not a final removal count" % format_count(screen.get("unusable"), missing="not available")),
        ])
    if validation.get("comparable_variants") is not None:
        fields.extend([
            ("count", "Replacement comparable variants", _pair(validation.get("comparable_variants"), validation.get("study_variant_count"))),
            ("info", "Replacement consistency", "%s; correlation %s; mean absolute difference %s" % (
                "Passed" if validation.get("status") == "passed" else "Not approved",
                format_number(validation.get("correlation")), format_number(validation.get("mean_absolute_difference")))),
        ])
    fields.extend(_frequency_failure_fields(validation, "Replacement"))
    if validation.get("failure_reason", validation.get("decision_reason")) in (
        "replacement_frequency_maf_like", "no_usable_replacement_frequencies"):
        fields.append(("error", "Stop action", "No further default comparison, replacement or fallback"))
    return fields


def frequency_summary_fields(decisions, *, policies=None, resources=None):
    """One concise study summary shared by terminal and HTML; evidence stays intact."""
    resources = resources or {}
    external = bool(decisions.get("eaf_external_validation_required")
        or "eaf_external_validation" in decisions or decisions.get("eaf_selected_source") == "external_eaf"
        or (decisions.get("eaf_initial_frequency") or {}).get("source_kind") == "external_eaf")
    initial = decisions.get("eaf_external_validation" if external else "eaf_reference_validation") or {}
    replacement = decisions.get("eaf_selected_source") == "replacement_reference"
    validation = (decisions.get("eaf_replacement_validation") or {}) if replacement else initial
    if replacement:
        outcome = _frequency_approval(validation, "Replacement", decisions.get("report_failure"))
    elif external:
        outcome = _frequency_approval(initial, "External EAF", decisions.get("report_failure"))
    else:
        outcome = {"passed": "Study EAF retained", "not_required": "Declared study EAF retained",
            "failed": "Failed; dataset stopped before chromosome analysis",
            "invalid": "Invalid; dataset stopped before chromosome analysis",
            "inconclusive": "Inconclusive; dataset stopped before chromosome analysis"}.get(initial.get("status"),
            "Input frequency assessment pending")
    if not decisions.get("eaf_column") and not external and not replacement:
        return [("info", "Frequency decision", "Not recorded; no AF source decision available")]
    fields = [("analysis", "Frequency decision", outcome)]
    fields.extend(_input_frequency_summary_fields(decisions, initial, external, replacement))
    references = {(r.get("default_eaf_file"), r.get("default_eaf_reference_column"))
        for r in resources.values() if r.get("default_eaf_file")}
    if references:
        populations = ", ".join(dict.fromkeys(_text(population) for _, population in sorted(references, key=lambda item: str(item))))
        fields.append(("info", "Reference", "%s chromosome files; population %s; full paths in run manifest" % (
            format_count(len(references)), populations) if len(references) > 1 else
            "%s; population %s" % next(iter(references))))
    fields.extend(_frequency_comparison_fields(initial))
    if initial.get("status") in ("failed", "invalid", "inconclusive"):
        fields.append(("error", "Action", "No chromosome analysis; correct the recorded AF validation failure"))
    if initial.get("status") == "not_required" or initial.get("decision_reason") == "external_frequency_not_predominantly_minor":
        fields.append(("info", "Reference comparison", "Not required; input was not MAF-like"))
    fields.extend(_frequency_failure_fields(initial, "Input"))
    if replacement:
        fields.extend(_replacement_summary_fields(decisions, validation))
    if external:
        fields.append(("info", "Frequency provenance", "Population-reference proxy; supplied external AF"))
    if replacement and (validation.get("source") or {}).get("same_source_evidence_reused"):
        fields.append(("info", "Source reuse", "Prepared reference reused; this is not independent confirmation"))
    if decisions.get("report_failure"):
        fields.append(("error", "Failure", _text(decisions["report_failure"])))
    fields.append(("info", "Detailed evidence", "Thresholds, allele mappings, source identities and diagnostic exclusions are recorded in the run manifest"))
    return fields


def chromosome_frequency_source_fields(qc):
    """Render the worker's recorded source separately from original input MAF."""
    replacement = (qc.get("eaf_selected_source") or qc.get("decision_source")) == "replacement_reference"
    if not replacement:
        return None
    evidence = qc.get("eaf_replacement_validation") or {}
    approved = evidence.get("status") == "passed" and bool(evidence.get("contract"))
    return [
        ("Source", "Approved replacement EAF" if approved else _frequency_approval(evidence, "Replacement")),
        ("Source file", _text(qc.get("eaf_replacement_source_file"))),
        ("Frequency column", _text(qc.get("eaf_replacement_source_column"))),
        ("Working column", _text(qc.get("final_eaf_col"))),
        ("Frequency interpretation", "EAF; whole-study decision reused" if approved else "Approval not established"),
        ("Assignment", "Prepared frequencies selected after allele orientation; original input retained unchanged"
            if approved else "Assignment not authorised"),
    ]


def _frequency_fields(decisions, policies, resources):
    return frequency_summary_fields(decisions, policies=policies, resources=resources)


def _info_fields(decisions):
    """Describe the saved source-scale decision without implying row-level QC."""
    evidence = decisions.get("info_score_type_evidence") or {}
    selected = decisions.get("info_score_type", evidence.get("resolved_type"))
    if not evidence and selected is None:
        return [("info", "Check status", "Not recorded; INFO source-scale resolution has not been established.")]
    fields = [
        ("info", "Selected source", _text(evidence.get("source"))),
        ("info", "Decision source", _source(decisions.get("info_score_type_source", evidence.get("decision_source")))),
        ("info", "Configured score type", _words(evidence.get("configured_type"))),
    ]
    if evidence:
        finite = evidence.get("finite_values")
        standard_max = format_number(evidence.get("standard_info_maximum"))
        mach_max = format_number(evidence.get("mach_rsq_maximum"))
        fields.extend([
            ("count", "Selected source rows", format_count(evidence.get("source_rows"), missing="not available")),
            ("count", "Finite numeric values", format_count(finite, missing="not available")),
            ("count", "Missing or non-finite values", format_count(evidence.get("unusable_values"), missing="not available")),
            ("info", "Observed source range", "%s to %s" % (format_number(evidence.get("minimum")), format_number(evidence.get("maximum")))),
            ("count", "Below configured minimum", format_count(evidence.get("below_minimum"), missing="not available")),
            ("count", "Standard range values (through %s)" % standard_max, format_count(evidence.get("standard_range_values"), missing="not available")),
            ("info", "MaCH range values (above %s through %s)" % (standard_max, mach_max), _pair(evidence.get("mach_range_values"), finite)),
            ("info", "Automatic MaCH fraction threshold", "at least %s" % format_fraction_percentage(evidence.get("auto_mach_rsq_fraction"))),
            ("info", "Above MaCH maximum %s" % mach_max, _pair(evidence.get("above_mach_maximum"), finite)),
            ("info", "Maximum invalid fraction", format_fraction_percentage(evidence.get("maximum_invalid_fraction"))),
        ])
    fields.append(("analysis", "Selected score type", _words(selected) if selected is not None else "Not resolved"))
    if selected is not None:
        fields.append(("info", "Planned action", "Use the recorded source score type; normal chromosome INFO QC still applies."))
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
    policies=None, resource_maps=None, only=None, include_title=True,
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
        ("analysis", "Imputation-quality score", _info_fields(decisions)),
    ])
    lines = ([screen_line("error" if failure else "decision", "Study-wide decisions and evidence", indent=4)]
        if include_title else [])
    for kind, title, fields in sections:
        if only is None or title in only:
            lines.extend(screen_summary_card(kind, title, fields))
    if failure:
        lines.extend(screen_summary_card("error", "Dataset stopped", [("error", "Reason", str(failure))]))
    return "\n".join(lines)


class StudyDecisionReporter:
    """Publish each completed section once, using invocation-local scalar state."""

    sections = ("Genome build", "Effect type", "Standard-error scale", "P-value type",
        "BETA/Z direction for missing-SE recovery", "Allele frequency", "Strand consensus",
        "Imputation-quality score")

    def __init__(self, emit, render=render_study_decisions):
        self.emit = emit
        self.render = render
        self.published = set()
        self.pending_checks = set()
        self.started = False
        self.failed = False

    def publish(self, decisions, names, **context):
        unknown = set(names) - set(self.sections)
        if unknown:
            raise ValueError("Unknown study report sections: %s" % sorted(unknown))
        pending = tuple(name for name in names if name not in self.published)
        if not pending and not decisions.get("report_failure"):
            return
        text = self.render(decisions, only=pending, include_title=not self.started, **context)
        self.emit(text)
        self.started = True
        self.published.update(pending)

    def pending(self, name, message):
        if name not in self.pending_checks:
            self.emit("\n".join(screen_summary_card("decision", "Pending study-wide check",
                [("info", name, message)])))
            self.pending_checks.add(name)

    def fail(self, decisions, error, **context):
        if self.failed:
            return
        self.failed = True
        self.publish({**decisions, "report_failure": str(error)}, self.sections, **context)


__all__ = ["render_study_decisions", "render_replacement_frequency_validation", "replacement_frequency_fields",
    "render_external_frequency_validation", "external_frequency_fields", "native_frequency_screen_fields"]
