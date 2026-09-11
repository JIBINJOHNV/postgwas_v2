"""Presentation-only workflow cards from persisted harmonisation evidence.

The stage order below is the orchestration contract in ``service.py`` (8 study,
16 chromosome and 5 post-merge stages), not a second configuration source.
These functions neither open scientific inputs nor infer missing outcomes.
In particular, an overall successful dataset is not evidence that an optional
check ran, and a logger's VCF-stage row count is not a VCF record count.
"""

from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
from math import isfinite
from typing import Any


# Only workflow identities and their purpose are invariant; thresholds, sources
# and actions come exclusively from the recorded run.
_DATASET = (
    ("Configuration validation", "Check the resolved sample-sheet mappings, required inputs and configuration before analysis.", ("columns.", "input.")),
    ("Header validation", "Check that the configured columns exist in the study header.", ("columns.", "input.")),
    ("Read and prepare the study", "Parse missing and numeric values, normalize coordinates and alleles, apply supported-chromosome rules and handle duplicates.", ("input.", "columns.", "chromosome.", "allele.", "duplicates.")),
    ("Content validation and recovery plan", "Check usable inputs, plan missing sample-size handling and enforce chromosome-X Z-only reconstruction safety before workers start.", ("sample_size.", "effect_from_z.x_")),
    ("Genome-build determination", "Compare testable study variants with configured build references, or validate the declared build.", ("build.",)),
    ("Study-wide strand consensus", "Use informative non-palindromic reference matches to determine the study-wide strand evidence.", ("strand.",)),
    ("Study-property decisions", "Resolve effect type, SE scale, p-value scale and initial frequency interpretation once for the study; compare declarations with evidence.", ("effect.", "pvalue.", "eaf.", "study_properties.")),
    ("Reference preflight and chromosome preparation", "Validate exact build/chromosome resources, stage shared references, resolve INFO metric type and prepare chromosome partitions with bounded scheduling.", ("execution.", "info.", "input.chromosome_")),
)

_CHROMOSOME = (
    ("Load chromosome", "Read the typed chromosome partition prepared at study level.", (), None),
    ("Reference resources", "Recheck the configured resources for this chromosome and genome build.", (), None),
    ("Effect-scale conversion", "Apply the study-level effect and SE-scale decisions before any allele swap.", ("effect.",), "beta_or_oddsratio_qc"),
    ("Allele alignment and effect frequency", "Match reference alleles, resolve strand and palindromic orientation, validate frequency interpretation and apply EAF policies. Strand removals are part of this step, not an additional total.", ("strand.", "eaf.", "duplicates."), "eaf_qc"),
    ("Sample size", "Use the configured counts and trait type, then apply the recorded missingness plan and sample-size checks.", ("sample_size.",), "sample_size_qc"),
    ("Recover BETA or SE from Z", "Use the applicable supplied-BETA/SE/Z combination; reconstruct both only when needed and permitted, with the recorded method and assumptions.", ("effect_from_z.",), "effect_from_z_qc"),
    ("P-value normalization", "Apply the study-level p-value scale decision and validate or transform values under the resolved policies.", ("pvalue.",), "pval_detection_and_conversion_qc"),
    ("Recover remaining SE from BETA and P", "Keep usable SE and recover missing SE only where the recorded policy and p-value evidence permit it.", ("pvalue.", "validation.se_from_clipped_pval"), "se_from_beta_pvalue_qc"),
    ("Z-score processing", "Process Z using the available BETA and SE while preserving the recorded source and validity checks.", ("validation.",), "z_from_beta_se_qc"),
    ("Effect-statistic validation", "Check BETA, SE and Z validity and their configured statistical agreement, including p-values.", ("validation.",), "effect_statistics_validation_qc"),
    ("Imputation quality", "Use the selected internal, external or explicitly fixed INFO source and apply its metric-specific range and missingness policies.", ("info.",), "info_qc"),
    ("Variant identifiers", "Preserve and clean study identifiers; fill missing identifiers from chromosome, position and alleles after applying the configured missing tokens.", ("input.null_values",), None),
    ("Final completeness and row accounting", "Enforce the configured final required fields; reconcile harmonised and rejected rows before export when rejection provenance is enabled.", ("final_check.", "rejects."), "final_completeness_qc"),
    ("Harmonised table export", "Export the checked scientific columns and audit-only strand action for GWAS-to-VCF.", (), None),
    ("GWAS-to-VCF conversion", "Convert the harmonised table against the source-build FASTA; adapter losses are separate from earlier harmonisation rejections.", (), None),
    ("VCF normalization, annotation and liftover", "Normalize variants, annotate exact-allele IDs and reference frequencies, add consequences, lift over and normalize/index the target VCF. Separate failed liftover from configured exclusions of successfully lifted variants.", ("vcf.",), None),
)

_POST_MERGE = (
    ("Merge and validate chromosome VCFs", "Concatenate required output groups and validate chromosome coverage, merged outputs and indexes, using the configured retry/failure policy.", ("vcf.",)),
    ("Population-frequency similarity", "Compare study AF with the configured populations on the unfiltered input-build merged VCF; similarity is not proof of ancestry.", ()),
    ("Combine audit files and clean intermediates", "Merge side files and independently reconcile adapter-input counts before configured cleanup.", ("rejects.", "vcf.")),
    ("Merged-VCF quality assessment", "Assess the raw VCF and the virtual subset passing all active QC rules. This assessment does not itself filter the delivered merged VCF; rule failures may overlap.", ()),
    ("Finalize outputs, reports and provenance", "Validate delivered outputs, save summaries and logs, and retain or remove the raw GWAS-to-VCF intermediate according to the recorded policy.", ("vcf.keep_gwas2vcf_intermediate",)),
)

# These are recorded field identities, not new scientific calculations.
_OBSERVED = {
    "status": "Recorded outcome", "final_status": "Recorded outcome",
    "source": "Input source", "effect_col": "Effect input column",
    "effect_column": "Effect column", "effect_type": "Effect type",
    "effect_decision_source": "Decision source", "conversion": "Transformation",
    "se_input_scale": "Supplied SE scale", "se_rescaled_by_or": "SE converted from OR scale",
    "decision_source": "Decision source", "eaf_provenance": "Frequency source",
    "final_eaf_col": "Final frequency column", "maf_reference_decision": "Reference MAF/EAF decision",
    "ncase_source": "Case-count source", "ncontrol_source": "Control or total-count source",
    "trait_type": "Trait type", "Neff_status": "Sample-size handling",
    "se_column": "SE column", "z_column": "Z column",
    "beta_computed": "BETA derived", "se_computed": "SE derived",
    "effect_estimate_method_label": "Reconstruction method",
    "effect_estimate_citation": "Reconstruction reference",
    "effect_estimate_scale": "Reconstructed effect scale",
    "detected_scale": "P-value input scale", "pvalue_column": "P-value column",
    "pvalue_decision_source": "P-value decision source",
    "pvalue_out_of_range_action": "Out-of-range P action",
    "z_source": "Z source", "info_column": "INFO column",
    "info_score_type": "INFO metric type", "info_on_missing_action": "Missing INFO action",
    "on_missing": "Missing-value action", "resource_source": "Resource validation",
    "gwas2vcf_exit_code": "GWAS-to-VCF exit code", "target_vcf": "Target-build VCF",
    "inferred_build": "Selected genome build", "mode": "Build mode",
    "strand": "Study strand", "reason": "Reason", "pvalue_type": "P-value type",
    "se_scale": "SE scale", "frequency_type": "Frequency interpretation",
    "eaf_is_maf_initial": "Initial MAF suspicion", "effect_type_source": "Effect decision source",
    "pvalue_type_source": "P-value decision source", "action": "Configured action",
    "merge_status": "Merge outcome", "closest_population": "Most similar reference population",
    "decision_reason": "Decision explanation", "definition": "QC definition",
    "keep_requested": "Keep intermediate requested", "retained": "Intermediate retained",
    "retention_reason": "Retention reason",
    "total_adapter_input_rows": "Rows supplied to GWAS-to-VCF",
    "raw_variants": "Raw merged-VCF records", "qc_passed_variants": "Virtual QC-passed records",
    "excluded": "Records failing at least one virtual QC rule",
    "testable_variants": "Build-testable study variants",
    "untestable_variants": "Study variants not build-testable",
    "informative": "Strand-informative variants", "forward": "Forward-strand evidence",
    "reverse": "Reverse-strand evidence", "ambiguous": "Ambiguous strand evidence",
    "missing_variants": "Variants missing sample-size inputs",
    "se_recovered_from_z": "Missing SE recovered from Z",
    "se_unavailable_from_z": "Missing SE not recovered from Z",
    "beta_z_sign_mismatches": "BETA/Z direction disagreements",
    "se_missing_recovered_from_pvalue": "Missing SE recovered from P",
    "se_from_zero_p_approximation": "SE reconstructed with explicit zero-P approximation",
    "identifiers_rewritten": "Identifiers with separators rewritten",
    "missing_identifiers": "Missing identifiers detected",
    "identifiers_filled": "Missing identifiers filled",
    "rescaled_within_tolerance": "INFO values corrected within rounding tolerance",
    "rejected_out_of_range": "Out-of-range INFO variants rejected",
    "rejected_missing": "Missing INFO variants rejected",
    "variants_removed": "Total removed within this combined stage",
    "variants_removed_total": "Variants removed within this stage",
}


def _mapping(value: Any) -> dict:
    return dict(value) if isinstance(value, Mapping) else {}


def _value(value: Any) -> str:
    if value is None:
        return "not recorded"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return f"{value:,}"
    return str(value)


def _observed(payload: Mapping[str, Any]) -> list[str]:
    # Translate only recorded outcome enums, never arbitrary column/path text.
    display = dict(payload)
    for key, translations in {
        "conversion": {"OR_to_Beta_log_transform_applied": "odds ratios converted to BETA using ln(OR)"},
        "Neff_status": {"calculated_from_case_control": "effective sample size calculated as 4 / (1/Ncase + 1/Ncontrol)",
                        "not_computed_no_inputs": "not calculated because no sample-size inputs were available"},
        "z_source": {"calculated_from_beta_se": "calculated as BETA / SE"},
        "status": {"skipped_no_input_effect_column": "no supplied effect column; no OR conversion needed"},
    }.items():
        if isinstance(display.get(key), str):
            display[key] = translations.get(display[key], display[key])
    neff = payload.get("Neff_status")
    if isinstance(neff, str) and neff.startswith("fallback_ncontrol_only:"):
        display["Neff_status"] = "control-count input used as total sample size; column " + neff.split(":", 1)[1]
    return [f"{label}: {_value(display[key])}." for key, label in _OBSERVED.items()
            if key in payload and not isinstance(payload[key], (Mapping, list, tuple))]


def _chromosome_groups(number: int, payload: Mapping, actions: list) -> list[dict]:
    """Expose existing sub-stage evidence; never sum or reinterpret row losses."""
    groups = []

    def add(title, values, fields, note):
        metrics = {label: deepcopy(values[key]) for key, label in fields if key in values}
        if metrics:
            groups.append({"title": title, "metrics": metrics, "note": note})

    if number == 4:
        strand = _mapping(payload.get("strand_orientation"))
        orientation = {**strand, **_mapping(strand.get("actions"))}
        add("Strand and allele orientation", orientation, (
            ("reference_file", "Reference file"), ("reference_population_column", "Reference population"),
            ("initial_variants", "Variants evaluated"), ("forward", "Alleles unchanged"),
            ("forward_swapped", "Alleles swapped"), ("reverse_complement", "Reverse-complement"),
            ("reverse_complement_swapped", "Complement and swap"),
            ("reference_unmatched", "Unmatched removed"), ("reference_unmatched_retained", "Unmatched retained"),
            ("palindromic_frequency_conflict", "Palindromic consensus/AF conflicts removed"),
            ("palindromic_orientation_unavailable", "Palindromic orientation unavailable"),
            ("palindromic_frequency_discordant", "Palindromic frequency failures removed"),
            ("palindromic_ambiguous", "Palindromic ambiguous removed"),
            ("reference_ambiguous", "Reference ambiguous removed"), ("final_variants", "Retained after strand alignment"),
        ), "Palindromic conflicts mean study-wide strand and allele-frequency evidence supported opposite orientations.")
        fields = [("final_eaf_col", "Final EAF column"), ("strand_af_tolerance", "Maximum absolute AF difference")]
        notes = []
        for label, prefix, action_key in (
            ("Non-palindromic", "strand_non_palindromic_af", "strand_af_discordance_action"),
            ("Palindromic", "strand_palindromic_af", "strand_palindromic_af_discordance_action"),
        ):
            fields.extend(((prefix + "_comparable", label + " comparable"), (prefix + "_discordant", label + " discordant")))
            fields.append((action_key, label + " disagreement policy"))
            action_text = {"warn": "retained with warning; no frequencies changed by this check",
                           "reject": "removed by this check", "fail": "configured to stop on disagreement"}.get(payload.get(action_key))
            if action_text and prefix + "_discordant" in payload:
                notes.append(label + " disagreements: " + action_text + ".")
        fields.extend((("variants_removed", "Combined strand/EAF removals"), ("final_variants", "Retained after the complete strand/EAF step")))
        add("Frequency agreement", payload, fields, " ".join(notes) + " Combined strand/EAF removals already include strand losses; do not add the subgroup counts again.")
    if number == 16:
        by_check = {str(action.get("check")): action for action in actions}
        add("Source-build normalization", _mapping(by_check.get("NORM")), (("before", "Input VCF records"), ("after", "Normalized VCF records")), "Counts include the recorded normalization and exact-duplicate handling.")
        annotation = {name: by_check[key].get("after") for key, name in (("ID", "After ID annotation"), ("EAF", "After reference-frequency annotation"), ("CSQ", "After consequence annotation")) if key in by_check}
        add("VCF annotations", annotation, [(name, name) for name in annotation], "These are successive record counts, not separate removal totals.")
        add("Target-build liftover", _mapping(payload.get("liftover_accounting")), (
            ("input", "Variants entering liftover"), ("rejected", "Plugin rejections"),
            ("swap_excluded", "Successfully lifted variants excluded by swap policy"),
            ("swap_policy", "Swap policy"), ("final", "Target-build VCF records"),
        ), "Successfully lifted variants excluded by policy are not failed liftover records.")
    return groups


def _stage_status(record: Mapping, payload: Mapping) -> str:
    statuses = [str(value).lower() for value in (
        record.get("status"), payload.get("status"), payload.get("merge_status"),
    ) if value is not None]
    if any(status in {"failed", "fail", "error", "frequency_inversion_suspected"} for status in statuses):
        return "Failed"
    if "failed_warning" in statuses:
        return "Failed — continued by policy"
    if any(status in {"running", "in_progress"} for status in statuses):
        return "Running"
    if "skipped_incomplete_input_build_vcf" in statuses:
        return "Not run — incomplete input"
    if any(status in {"disabled", "skipped", "not_applicable"} or status.startswith("skipped_") for status in statuses):
        return "Not needed"
    if "partial" in statuses:
        return "Partial"
    if "warning" in statuses:
        return "Completed with warnings"
    if any(status in {"ok", "passed", "success", "complete", "completed"} for status in statuses):
        return "Completed"
    return "Evidence recorded" if record or payload else "Not recorded"


def _card(number: int, definition: tuple, record: Mapping, payload: Mapping,
          policies: Mapping, *, row_counts: bool = True) -> dict:
    title, purpose, prefixes = definition[:3]
    metrics = {label: record[key] for key, label in (
        ("rows_in", "Rows entering step"), ("rows_out", "Rows leaving step"),
        ("removed", "Rows removed at this step"), ("elapsed", "Elapsed seconds"),
    ) if key in record and (row_counts or key == "elapsed")}
    messages = _observed(payload)
    if record.get("error"):
        messages.insert(0, str(record["error"]))
    outcome = _mapping(_mapping(record.get("extra")).get("outcome"))
    if outcome.get("message"):
        messages.insert(0, str(outcome["message"]))
    if not messages:
        messages = ["The step outcome was recorded; detailed evidence is below." if record
                    else "Detailed evidence is available below; a separate completion outcome was not recorded." if payload
                    else "No stage evidence was retained in this manifest; this does not mean zero changes or a passed check."]
    return {
        "number": number, "title": title, "purpose": purpose,
        "status": _stage_status(record, payload), "summary": messages,
        "metrics": deepcopy(metrics),
        "policies": deepcopy({key: value for key, value in policies.items()
                              if any(str(key).startswith(prefix) for prefix in prefixes)}),
        "details": deepcopy({"Recorded stage": dict(record), "Observed evidence": dict(payload)}),
    }


def _records(container: Mapping, total: int) -> dict:
    # The total separates dataset/post-merge steps, which reuse step numbers.
    return {step["number"]: _mapping(step) for step in container.get("steps") or []
            if isinstance(step, Mapping) and step.get("total") == total
            and isinstance(step.get("number"), int)}


def _evidence_coverage(expected: list[str], recorded: list[str]) -> dict:
    """Coverage describes evidence, never successful scientific execution."""
    missing = [chromosome for chromosome in expected if chromosome not in recorded]
    return {"expected_chromosomes": len(expected), "recorded_chromosomes": len(recorded),
            "missing_chromosomes": missing, "complete": bool(expected) and not missing}


def _recorded_number(value: Any, *, count: bool = False) -> int | float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not isfinite(value):
        return None
    if count:
        return int(value) if value >= 0 and value == int(value) else None
    return value


def _count_difference(before: Any, after: Any) -> int | None:
    before, after = (_recorded_number(value, count=True) for value in (before, after))
    return before - after if before is not None and after is not None and before >= after else None


def _evidence_group(expected: list[str], payloads: Mapping[str, Mapping],
                    counts: tuple[str, ...], sources: tuple[str, ...]) -> dict:
    """Sum a metric across chromosomes, not across possibly overlapping checks."""
    result = {"counts": {}, "sources": {}, "source_coverage": {},
              "coverage": _evidence_coverage(expected, [c for c, p in payloads.items() if p])}
    for key in counts:
        values = {c: _recorded_number(p.get(key), count=True) for c, p in payloads.items()}
        values = {c: value for c, value in values.items() if value is not None}
        result["counts"][key] = {
            "value": sum(values.values()) if values else None,
            **_evidence_coverage(expected, list(values)),
        }
    for key in sources:
        values = {c: p[key] for c, p in payloads.items()
                  if key in p and isinstance(p[key], (str, bool, int, float))
                  and (not isinstance(p[key], float) or isfinite(p[key]))}
        result["sources"][key] = sorted(set(values.values()), key=str)
        result["source_coverage"][key] = _evidence_coverage(expected, list(values))
    return result


def _effect_evidence(stages: Mapping[int, Mapping]) -> dict:
    effect = stages[6]
    invalid_reconstructed = _recorded_number(effect.get("invalid_reconstructed_statistics"), count=True)
    # BETA=Z*SE can leave null results; whole-frame row counts alone therefore
    # cannot establish the number reconstructed. Z-only reconstruction records
    # an explicit finite/positive validation before returning its row count.
    beta_count = None
    if effect.get("beta_computed") is False:
        beta_count = 0
    elif effect.get("beta_computed") is True:
        if "beta_undefined_from_z" in effect:
            beta_count = _count_difference(effect.get("variants_remaining"), effect["beta_undefined_from_z"])
        elif invalid_reconstructed == 0:
            beta_count = _recorded_number(effect.get("variants_remaining"), count=True)
    reconstructed_se = None
    if effect.get("has_beta") is True or effect.get("se_computed") is False:
        reconstructed_se = 0
    elif (effect.get("se_computed") is True and effect.get("has_se") is False
          and invalid_reconstructed == 0):
        reconstructed_se = _recorded_number(effect.get("variants_remaining"), count=True)
    result = {**stages[3], **effect, **stages[8], **stages[9]}
    result.update(beta_derived_from_z=beta_count, se_reconstructed_from_z=reconstructed_se)
    if "se_unavailable_from_z" not in result and "se_undefined_from_z" in effect:
        result["se_unavailable_from_z"] = effect["se_undefined_from_z"]
    return result if any(stages[number] for number in (3, 6, 8, 9)) else {}


def summarise_chromosome_evidence(manifest: Mapping[str, Any]) -> dict:
    """Aggregate saved overview evidence without opening input or VCF files.

    Every count retains its chromosome coverage: ``None`` means unavailable,
    zero means explicitly recorded, and partial sums are not whole-study totals.
    Sources are unique recorded values, including an unresolved ``auto`` trait.
    Sample-size ranges use extrema only, never averages of chromosome medians.
    Retained discordances mean retained *by that check*, not necessarily present
    in final outputs. Counts from different checks must not be added together.
    """
    dataset = _mapping(manifest.get("dataset"))
    summaries = {str(c): _mapping(s) for c, s in _mapping(dataset.get("chromosome_summaries")).items()}
    chromosome_ids = set(map(str, _mapping(dataset.get("chromosomes")))) | set(summaries)
    for key in ("completed", "failed"):
        if isinstance(dataset.get(key), (list, tuple)):
            chromosome_ids.update(str(value) for value in dataset[key]
                                  if isinstance(value, (str, int)) and not isinstance(value, bool))
    expected = sorted(chromosome_ids)
    payloads = {name: {} for name in ("effects", "frequency", "pvalues", "validation", "sample_size", "info", "vcf")}
    for chromosome, summary in summaries.items():
        records = _records(summary, len(_CHROMOSOME))
        qc = _mapping(summary.get("stage_qc"))
        stages = {number: {**_mapping(_mapping(records.get(number)).get("extra")),
                           **_mapping(qc.get(definition[3]))}
                  for number, definition in enumerate(_CHROMOSOME, 1)}
        payloads["effects"][chromosome] = _effect_evidence(stages)
        for name, number in (("frequency", 4), ("pvalues", 7), ("validation", 10), ("info", 11)):
            payloads[name][chromosome] = dict(stages[number])
        payloads["validation"][chromosome].update(_mapping(stages[10].get("basic_matched_by_reason")))
        if stages[11]:
            payloads["info"][chromosome]["removed_total"] = _count_difference(
                stages[11].get("initial_variants"), stages[11].get("final_variants"))
        sample = dict(stages[5])
        missing = _mapping(sample.get("missing_sample_size"))
        if missing:
            sample.update({key: missing[key] for key in ("variants_affected", "variants_removed", "variants_imputed") if key in missing})
            sample["missing_action"] = missing.get("action")
        payloads["sample_size"][chromosome] = sample
        actions = {str(action.get("check")): action for action in summary.get("qc_actions") or []
                   if isinstance(action, Mapping) and str(action.get("step", "")).split(" ", 1)[0] == "16"}
        normal = _mapping(actions.get("NORM"))
        lifted = _mapping(stages[16].get("liftover_accounting"))
        adapter_input = _mapping(records.get(15)).get("rows_in")
        payloads["vcf"][chromosome] = {
            "adapter_input": adapter_input, "adapter_output": normal.get("before"),
            "adapter_removed": _count_difference(adapter_input, normal.get("before")),
            "normalization_input": normal.get("before"), "normalization_output": normal.get("after"),
            "normalization_removed": normal.get("removed"),
            **{label: lifted.get(key) for key, label in (("input", "liftover_input"), ("rejected", "liftover_rejected"),
                ("swap_excluded", "liftover_swap_excluded"), ("final", "target_output"), ("swap_policy", "swap_policy"))},
        } if normal or lifted or adapter_input is not None else {}

    # All field names below are persisted metric identities, not configurable
    # scientific policies or column defaults. Their producers own calculations.
    fields = {
        "effects": (("beta_derived_from_z", "se_reconstructed_from_z", "se_recovered_from_z", "se_unavailable_from_z",
                     "se_from_pvalue", "se_missing_recovered_from_pvalue", "se_from_clipped_pval", "se_from_zero_p_approximation",
                     "beta_z_sign_mismatches", "removed_beta_z_sign_mismatches"),
                    ("effect_type", "effect_col", "conversion", "se_input_scale", "has_beta", "has_se", "has_zscore",
                     "beta_computed", "se_computed", "effect_estimate_method_label", "effect_estimate_beta_formula",
                     "effect_estimate_se_formula", "effect_estimate_scale", "effect_estimate_citation", "z_source")),
        "frequency": (("strand_non_palindromic_af_comparable", "strand_non_palindromic_af_discordant",
                       "strand_palindromic_af_comparable", "strand_palindromic_af_discordant", "missing_eaf_count", "eaf_degenerate_removed"),
                      ("eaf_provenance", "final_eaf_col", "strand_af_discordance_action", "strand_palindromic_af_discordance_action")),
        "pvalues": (("variants_with_pvalues_clipped_low", "variants_with_pvalues_clipped_high", "variants_with_zero_pvalues_replaced",
                     "pvalues_preserved_below_float64", "variants_with_null_pvalues_removed", "variants_with_lt0_pvalues_removed",
                     "variants_with_non_finite_pvalues_removed", "variants_with_gt1_05_pvalues_removed"),
                    ("detected_scale", "pvalue_column", "pvalue_decision_source", "pvalue_out_of_range_action")),
        "validation": (("removed_total", "variants_removed_invalid_beta_se", "se_null", "se_non_finite", "se_non_positive", "beta_invalid", "beta_zero", "z_invalid"),
                       ("beta_zero_action",)),
        "sample_size": (("initial_variants", "final_variants", "variants_affected", "variants_removed", "variants_imputed"),
                        ("trait_type", "Neff_status", "ncase_source", "ncontrol_source", "missing_action")),
        "info": (("initial_variants", "missing_info", "missing_info_after", "rescaled_within_tolerance", "rejected_out_of_range",
                  "rejected_missing", "final_variants", "removed_total"),
                 ("source", "info_column", "info_score_type", "info_on_missing_action", "info_out_of_range_action")),
        "vcf": (("adapter_input", "adapter_output", "adapter_removed", "normalization_input", "normalization_output", "normalization_removed",
                 "liftover_input", "liftover_rejected", "liftover_swap_excluded", "target_output"), ("swap_policy",)),
    }
    result = {name: _evidence_group(expected, payloads[name], *spec) for name, spec in fields.items()}
    result["coverage"] = _evidence_coverage(expected, list(summaries))
    for name, key in (("beta_se_z", "beta_se_z_concordance"), ("z_p", "z_pval_concordance")):
        checks = {c: _mapping(p.get(key)) for c, p in payloads["validation"].items()}
        for check in checks.values():
            if check.get("status") == "ran" and check.get("action") in {"warn", "reject"}:
                check["retained_discordant"] = _count_difference(check.get("discordant"), check.get("removed"))
        result["validation"][name] = _evidence_group(
            expected, checks, ("checked", "discordant", "removed", "retained_discordant"), ("status", "action"))
    sample_stats = {c: _mapping(p.get("sample_size_stats")) for c, p in payloads["sample_size"].items()}
    ranges = {}
    for column in sorted({key for stats in sample_stats.values() for key in stats}):
        values = {c: (_recorded_number(_mapping(stats.get(column)).get("min")),
                      _recorded_number(_mapping(stats.get(column)).get("max"))) for c, stats in sample_stats.items()}
        values = {c: (low, high) for c, (low, high) in values.items() if low is not None and high is not None and low <= high}
        if values:
            ranges[column] = {"min": min(pair[0] for pair in values.values()), "max": max(pair[1] for pair in values.values()),
                              "coverage": _evidence_coverage(expected, list(values))}
    result["sample_size"].update(ranges=ranges, missing_plan=deepcopy(_mapping(dataset.get("missing_sample_size"))))
    return result


def dataset_steps(manifest: Mapping[str, Any]) -> list[dict]:
    """Return all eight study stages, retaining absent completion evidence."""
    dataset = _mapping(manifest.get("dataset"))
    study = _mapping(dataset.get("study_decisions"))
    strand = _mapping(study.get("strand_consensus"))
    if "strand" in study:
        strand["strand"] = study["strand"]
    pre_vcf = _mapping(manifest.get("pre_vcf_summary"))
    input_counts = {key: pre_vcf[key] for key in (
        "total_variant_infile", "total_variant_read", "total_variant_removed_missing_values",
        "total_variant_removed_null_coords", "total_variant_removed_unsupported_chromosomes",
        "total_variant_removed_non_standard_alleles", "total_variant_removed_duplicates",
        "total_variant_remaining_for_harmonisation", "total_variant_ready_snps",
        "total_variant_ready_indels_or_other",
    ) if key in pre_vcf}
    payloads = [
        {"Resolved inputs": _mapping(manifest.get("config"))} if manifest.get("config") else {},
        {},
        {key: value for key, value in {
            "Input accounting": input_counts or None,
            "Field completeness": dataset.get("field_completeness"),
        }.items() if value is not None},
        _mapping(dataset.get("missing_sample_size")),
        _mapping(dataset.get("genome_build")), strand,
        {key: value for key, value in study.items() if not key.startswith("info_")},
        {key: value for key, value in {
            "Resource preflight": dataset.get("resource_preflight"),
            "Parallelism": dataset.get("parallelism"),
            "INFO metric decision": study.get("info_score_type_evidence"),
        }.items() if value is not None},
    ]
    if study.get("x_chromosome_z_only_policy") is not None:
        payloads[3]["Chromosome-X reconstruction safety"] = study["x_chromosome_z_only_policy"]
    records = _records(manifest, len(_DATASET)) | _records(dataset, len(_DATASET))
    cards = [_card(number, definition, records.get(number, {}), payload,
                   _mapping(manifest.get("policies")))
             for number, (definition, payload) in enumerate(zip(_DATASET, payloads), 1)]
    if manifest.get("sumstat_file"):
        cards[2]["summary"].insert(0, f"Study input: {manifest['sumstat_file']}.")
    cards[2]["metrics"].update({label: input_counts[key] for key, label in (
        ("total_variant_infile", "Input variants"), ("total_variant_read", "Variants read"),
        ("total_variant_remaining_for_harmonisation", "Ready for chromosome analysis"),
    ) if key in input_counts})
    return cards


def chromosome_steps(summary: Mapping[str, Any], policies: Mapping | None = None) -> list[dict]:
    """Return sixteen stages without adding overlapping exclusion counts."""
    records = _records(summary, len(_CHROMOSOME))
    stages = _mapping(summary.get("stage_qc"))
    cards = []
    for number, definition in enumerate(_CHROMOSOME, 1):
        record = records.get(number, {})
        payload = _mapping(record.get("extra"))
        payload.update(_mapping(stages.get(definition[3])))
        card = _card(number, definition, record, payload, _mapping(policies), row_counts=number < 15)
        actions = [deepcopy(action) for action in summary.get("qc_actions") or []
                   if isinstance(action, Mapping)
                   and str(action.get("step", "")).split(" ", 1)[0] == f"{number:02d}"]
        if actions:
            card["details"]["Recorded checks (not additive)"] = actions
        card["groups"] = _chromosome_groups(number, payload, actions)
        if number == 13:
            card["details"]["Pre-export reconciliation"] = {
                key: deepcopy(summary[key]) for key in ("rows_in", "rows_out", "rejected", "balanced", "reject_counts")
                if key in summary
            }
        if number == 15 and record:
            card["summary"].append("The step row count describes the harmonised input, not the number of VCF records created.")
        if number == 16:
            accounting = _mapping(payload.get("liftover_accounting"))
            card["metrics"].update({label: accounting[key] for key, label in (
                ("input", "Variants entering liftover"), ("rejected", "Liftover plugin rejections"),
                ("swap_excluded", "Successfully lifted variants excluded by swap policy"),
                ("final", "Target-build VCF records"),
            ) if key in accounting})
        cards.append(card)
    return cards


def post_merge_steps(manifest: Mapping[str, Any]) -> list[dict]:
    """Return five post-merge stages; optional concordance follows separately."""
    payloads = [
        _mapping(manifest.get("merged_vcfs")),
        _mapping(manifest.get("population_frequency_qc")),
        _mapping(manifest.get("gwas2vcf_summary_audit")),
        _mapping(manifest.get("qc_assessment")),
        _mapping(manifest.get("gwas2vcf_intermediate")),
    ]
    records = _records(manifest, len(_POST_MERGE))
    cards = [_card(number, definition, records.get(number, {}), payload,
                   _mapping(manifest.get("policies")))
             for number, (definition, payload) in enumerate(zip(_POST_MERGE, payloads), 1)]
    for key in ("combined_log", "completed_at", "vcf_provenance"):
        if key in manifest:
            cards[-1]["details"][key] = deepcopy(manifest[key])
    return cards
