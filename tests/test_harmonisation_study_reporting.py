"""Detailed reports preserve the distinction between evidence and decisions."""

from copy import deepcopy
import importlib
import os

import pytest
import polars as pl
from rich.cells import cell_len

from postgwas.config import load_configuration
from postgwas.core.ui.screen import terminal_presentation
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.study_properties import resolve_study_properties


def render(decisions, **kwargs):
    module = importlib.import_module("postgwas.modules.harmonisation.study_reporting")
    return module.render_study_decisions(decisions, **kwargs)


def normal(text):
    return " ".join(text.split())


@pytest.fixture
def decisions():
    return {
        "effect_type": "odds_ratio", "effect_type_source": "detector",
        "effect_column": "OR", "se_column": "SE", "pvalue_column": "P",
        "effect_evidence": {"non_null": 1200, "non_positive": 0,
                            "median": 1.01, "min": 0.5, "max": 2.0},
        "pvalue_type": "raw", "pvalue_type_source": "sample_sheet",
        "pvalue_type_detected": "raw", "pvalue_type_matches_declaration": True,
        "pvalue_evidence": {"non_null": 1200, "null": 3, "unparseable": 2,
                            "min": 1e-100, "max": 1.0, "median": 0.4,
                            "above_threshold": 0, "negative": 0},
        "se_scale": "log_odds", "se_scale_source": "dataset_z_pvalue_cross_check",
        "se_scale_evidence": {"comparison_source": "supplied_z", "status": "decided",
                              "informative_variants": 1100, "log_odds_matches": 1089,
                              "raw_odds_ratio_matches": 220, "minimum_informative_variants": 1000,
                              "minimum_agreement_fraction": 0.95, "minimum_agreement_margin": 0.1},
        "beta_z_sign_consistency": {"applies": False,
                                    "reason": "effect or signed Z column is not supplied"},
        "eaf_column": "AF", "eaf_is_maf": False, "eaf_is_maf_initial": True,
        "eaf_selected_source": "study_supplied", "eaf_input_interpretation": "eaf",
        "eaf_evidence": {"usable": 1200, "unusable": 5, "input_missing": 3,
                         "unparseable": 1, "non_finite": 1, "out_of_range": 0,
                         "at_or_below_0_5": 1190, "low_fraction_of_usable": 1190/1200,
                         "min": 0.001, "max": 0.8, "median": 0.2},
        "eaf_maf_decision_cutoff": 0.95,
        "eaf_reference_validation": {
            "status": "passed", "decision_reason": "allele_aligned_frequency_agreement",
            "direct_eaf_accepted": True, "comparable_variants": 1100,
            "minimum_overlap": 1000, "correlation_method": "pearson",
            "direct_eaf_correlation": 0.98, "minimum_correlation": 0.9,
            "mean_eaf_error": 0.02, "maximum_mean_absolute_difference": 0.03,
            "reference_effect_allele_minor_fraction": 0.99,
        },
        "strand": "forward", "strand_mode": "auto", "forward": 1090,
        "reverse": 10, "ambiguous": 100, "informative": 1100,
        "dominant_fraction": 1090/1100, "minimum_informative_variants": 1000,
        "consensus_threshold": 0.99, "reason": "dominant strand fraction 0.9909",
    }


def test_report_uses_existing_evidence_and_actual_policy_limits_without_mutation(decisions, monkeypatch):
    policies = default_policies()
    render({}, policies=policies)  # Import/resolve presentation dependencies first.
    original = deepcopy(decisions)
    def forbidden(*args, **kwargs):
        pytest.fail("report rendering must not read study or reference files")
    monkeypatch.setattr("builtins.open", forbidden)
    text = normal(render(decisions, policies=policies))
    assert "Study-wide decisions and evidence" in text
    assert "Automatically inferred" in text
    assert "Sample-sheet declaration" in text
    assert "1,089 / 1,100 (99.00%)" in text
    assert "Direct EAF correlation" in text and "0.98; required >= 0.9" in text
    assert "0.02; required <= 0.03" in text
    assert "Study EAF retained" in text
    assert "99.09%" in text and "99.00%" in text
    assert decisions == original


def test_reference_comparison_not_required_is_not_presented_as_a_pass(decisions):
    decisions.update(eaf_is_maf_initial=False)
    decisions["eaf_reference_validation"] = {
        "status": "not_required", "decision_reason": "study_frequency_not_predominantly_minor",
    }
    text = normal(render(decisions))
    assert "Not required" in text
    assert "Declared study EAF retained" in text
    assert "Direct EAF correlation" not in text
    assert "Allele-aligned reference agreement passed" not in text


def test_maf_distribution_without_completed_reference_check_stays_provisional(decisions):
    decisions.update(eaf_is_maf=True, eaf_is_maf_source="study_level_statistic")
    decisions.pop("eaf_reference_validation")
    decisions.pop("eaf_selected_source")
    text = normal(render(decisions))
    assert "MAF-like; reference confirmation required" in text
    assert "Pending" in text
    assert "Study EAF retained" not in text


def test_selected_replacement_is_a_planned_source_not_a_completed_output(decisions):
    decisions.update(eaf_is_maf=True, eaf_selected_source="replacement_reference",
                     eaf_input_interpretation="maf",
                     eaf_replacement={"file": "/references/replacement [EUR].tsv.gz", "column": "EUR"})
    decisions["eaf_reference_validation"].update(
        status="replacement_selected", direct_eaf_accepted=False,
        maf_reference_decision="maf", maf_correlation=0.97,
        mean_absolute_maf_difference=0.01, mean_maf_error=0.01,
        mean_eaf_error=0.3, minimum_error_margin=0.05,
        decision_reason="maf_interpretation_supported",
    )
    text = normal(render(decisions))
    assert "Input MAF supported" in text
    assert "Replacement EAF selected" in text
    assert "/references/replacement [EUR].tsv.gz" in text
    assert "EUR" in text
    assert "assignment pending" in text
    assert "Study EAF retained" not in text
    assert "completed" not in text.lower()


def test_failed_or_inconclusive_reference_evidence_never_reports_accepted_eaf(decisions):
    decisions["eaf_reference_validation"].update(
        status="failed", direct_eaf_accepted=False,
        maf_reference_decision="inconclusive", decision_reason="insufficient_overlap",
        comparable_variants=20,
    )
    text = normal(render(decisions))
    assert "FAILED" in text
    assert "insufficient overlap" in text.lower()
    assert "20; required >= 1,000" in text
    assert "No chromosome analysis" in text
    assert "Study EAF retained" not in text


def test_absent_frequency_has_no_fabricated_distribution_or_comparison():
    text = normal(render({"eaf_column": None, "eaf_is_maf": None,
                          "eaf_is_maf_source": "no_study_frequency_column",
                          "eaf_selected_source": "external_eaf",
                          "eaf_reference_validation": {"status": "not_applicable"}}))
    assert "Supplied external frequency; no native GWAS frequency column mapped" in text
    assert "External EAF selected" in text
    assert "Usable frequencies" not in text
    assert "Direct EAF correlation" not in text


@pytest.mark.parametrize("strand", ["mixed", "unresolved", "reverse"])
def test_strand_outcome_is_reported_without_turning_it_into_forward(decisions, strand):
    decisions["strand"] = strand
    text = normal(render(decisions))
    strand_text = text.split("Strand consensus", 1)[1]
    assert strand in strand_text.lower()
    assert "Forward matches" in strand_text and "Reverse matches" in strand_text
    assert "1,090" in strand_text and "1,100" in strand_text


def test_reference_aligned_mode_is_a_declaration_not_an_inferred_consensus():
    text = normal(render({"strand": "reference_aligned", "strand_mode": "reference_aligned",
                          "forward": 0, "reverse": 0, "informative": 0,
                          "dominant_fraction": None,
                          "reason": "explicit reference-aligned declaration"}))
    assert "Explicit reference-aligned declaration" in text
    assert "Forward matches" not in text
    assert "Dominant fraction" not in text
    assert "0 / 0" not in text


def test_no_comparable_missing_se_rows_is_not_a_pass(decisions):
    decisions["beta_z_sign_consistency"] = {"applies": True,
                                           "comparable_missing_se_rows": 0,
                                           "discordant_rows": 0, "decision": "continue"}
    text = normal(render(decisions))
    section = text.split("BETA/Z direction for missing-SE recovery", 1)[1].split("Allele frequency", 1)[0]
    assert "Not assessed" in section
    assert "no comparable rows needing SE recovery" in section
    assert "0 / 0" not in section


def test_missing_se_and_beta_scales_are_not_fabricated(decisions):
    decisions.update(effect_type="beta", se_scale=None, se_scale_source="not_applicable_to_beta")
    text = normal(render(decisions))
    section = text.split("Standard-error scale", 1)[1].split("P-value type", 1)[0]
    assert "not applicable to beta" in section.lower()
    assert "Log-odds agreement" not in section


def test_unknown_values_are_unavailable_rather_than_zero_or_success():
    text = normal(render({"eaf_column": "AF", "eaf_evidence": {"usable": None},
                          "eaf_reference_validation": {"status": "failed", "direct_eaf_correlation": float('nan')}}))
    assert "not available" in text.lower()
    assert ": nan" not in text.lower()
    assert "Study EAF retained" not in text


def test_partial_failure_keeps_unvisited_properties_unknown_instead_of_absent():
    text = normal(render({"report_failure": "Effect type could not be established"}))
    assert "❌ Study-wide decisions and evidence" in text
    assert "Dataset stopped" in text
    assert "Effect type could not be established" in text
    assert "Input column : Not recorded" in text
    assert "Input column : Not supplied" not in text


def test_evidence_explains_type_choices_and_frequency_exclusions(decisions):
    decisions.update(pvalue_type_source="detector")
    text = normal(render(decisions))
    assert "The recorded OR detection rule was satisfied" in text
    assert "The negative-log10 count and fraction gates were not both met" in text
    assert "Median detection rule" in text
    assert "Missing input : 3" in text
    assert "Non-numeric input : 1" in text
    assert "Non-finite input : 1" in text


def test_type_evidence_uses_one_finite_subset_with_distinct_exclusions():
    frame = pl.DataFrame({
        "B": ["-0.1", "0.2", "inf", "-inf", None, "bad", "NaN", "0.3"],
        "P": ["0.01", "0.2", "inf", "inf", None, "bad", "NaN", "0.03"],
    })
    resolved = resolve_study_properties(frame, {
        "beta_or_col": "B", "declared_effect_type": "beta",
        "pval_col": "P", "declared_pvalue_type": "raw",
    })
    original = deepcopy(resolved)
    text = normal(render(resolved))
    effect = text.split("Effect type", 1)[1].split("Standard-error scale", 1)[0]
    pvalue = text.split("P-value type", 1)[1].split("BETA/Z direction", 1)[0]
    assert "Above 1.05 : 0 / 3 (0.00%)" in pvalue
    assert "-0.1 to 0.3; median 0.2" in effect
    assert "1 / 3 (33.33%)" in effect
    assert "Finite numeric values : 3" in effect
    assert "Finite numeric values : 3" in pvalue
    assert "0.01 to 0.2; median 0.03" in pvalue
    for section in (effect, pvalue):
        assert "Missing input : 1" in section
        assert "Non-numeric input : 1" in section
        assert "Non-finite input : 3" in section
        assert "Excluded from inference : 5" in section
    assert resolved == original


def test_explicit_types_do_not_authorise_actions_on_absent_columns():
    policies = default_policies().with_overrides({"effect.type": "odds_ratio", "pvalue.type": "raw"})
    resolved = resolve_study_properties(pl.DataFrame({"Z": [1.0, 2.0]}), {"imp_z_col": "Z"}, policies=policies)
    text = normal(render(resolved, policies=policies))
    effect = text.split("Effect type", 1)[1].split("Standard-error scale", 1)[0]
    pvalue = text.split("P-value type", 1)[1].split("BETA/Z direction", 1)[0]
    assert "odds ratio" in effect
    assert "Raw p-values" in pvalue
    assert "BETA = ln(OR)" not in effect
    assert "Keep raw probabilities" not in pvalue
    assert "recovery" in effect.lower()
    assert "recovery" in pvalue.lower()


def test_absent_beta_se_requires_recovery_rather_than_unchanged_se():
    resolved = resolve_study_properties(pl.DataFrame({"B": [0.1], "Z": [1.0]}), {
        "beta_or_col": "B", "imp_z_col": "Z", "declared_effect_type": "beta",
    })
    text = normal(render(resolved))
    section = text.split("Standard-error scale", 1)[1].split("P-value type", 1)[0]
    assert "recovery" in section.lower()
    assert "unchanged" not in section.lower()


def test_ambiguous_build_never_selects_resources():
    text = normal(render({"report_failure": "The leading builds tied"}, genome_build_info={
        "inferred_build": "Ambiguous", "ambiguous_reason": "The leading builds tied",
        "input_variants": 20, "testable_variants": 10, "untestable_variants": 10,
    }))
    assert "No build selected" in text
    assert "Use the selected Ambiguous" not in text
    assert "resources are not selected" in text


def test_build_names_and_reference_identity_come_from_actual_context(decisions):
    text = normal(render(decisions, genome_build_info={
        "inferred_build": "AssemblyA", "input_variants": 1200, "testable_variants": 1000,
        "untestable_variants": 200, "matches": {"AssemblyA": 990, "AssemblyB": 10},
        "percentages": {"AssemblyA": 99.0, "AssemblyB": 1.0},
        "input_percentages": {"AssemblyA": 82.5, "AssemblyB": 0.833},
        "matched_reference_marker_counts": {"AssemblyA": 990, "AssemblyB": 10},
        "reference_marker_counts": {"AssemblyA": 1100, "AssemblyB": 1100},
        "reference_percentages": {"AssemblyA": 90.0, "AssemblyB": 0.909},
    }, build_reference_files={"AssemblyA": "/references/A/check.tsv", "AssemblyB": "/references/B/check.tsv"}))
    assert "AssemblyA" in text and "AssemblyB matches" in text
    assert "/references/A/check.tsv" in text
    assert "GRCh37" not in text


@pytest.mark.parametrize("path,selected", [
    ("normal", "AssemblyA"), ("declared", "AssemblyA"),
    ("low_coverage_fallback", "AssemblyA"), ("ambiguous", "Ambiguous"),
])
def test_build_report_is_pooled_and_preserves_complete_chromosome_evidence(path, selected):
    info = {
        "inferred_build": selected, "acceptance_path": path,
        "input_variants": 2200, "testable_variants": 2000, "untestable_variants": 200,
        "matches": {"AssemblyA": 1900, "AssemblyB": 100},
        "percentages": {"AssemblyA": 95.0, "AssemblyB": 5.0},
        "input_percentages": {"AssemblyA": 1900 / 22, "AssemblyB": 100 / 22},
        "distinct_loci": {"exclusive_snp_loci": {"AssemblyA": 1800, "AssemblyB": 100}},
        "per_chromosome": {str(chromosome): {
            "assessment": "competitor_veto", "testable_variants": 100,
            "matches": {"AssemblyA": 95, "AssemblyB": 5},
            "exclusive_snp_loci": {"AssemblyA": 90, "AssemblyB": 5},
            "competitor_veto": True,
        } for chromosome in range(1, 23)},
        "low_coverage_fallback": {"status": "accepted" if path == "low_coverage_fallback" else "not_needed",
            "accepted": path == "low_coverage_fallback"},
    }
    if selected == "Ambiguous":
        info["ambiguous_reason"] = "The leading builds tied"
    before = deepcopy(info)
    text = normal(render({}, genome_build_info=info))
    assert text.count("Selected build") == 1
    assert "2,200" in text and "1,900" in text and "1,800" in text
    assert "build evidence" not in text
    assert "Chromosome 1" not in text
    assert info == before
    if path == "normal":
        assert "normal" in text
        assert "No build selected" not in text


@pytest.mark.parametrize("effect,pvalue", [("beta", "raw"), ("odds_ratio", "neglog10")])
def test_automatic_inference_recommendations_share_the_structured_report(decisions, effect, pvalue):
    decisions.update(effect_type=effect, pvalue_type=pvalue, pvalue_type_source="detector",
        effect_type_inference_warning="AUTOMATIC EFFECT-TYPE INFERENCE: repeated full evidence",
        pvalue_type_inference_warning="AUTOMATIC P-VALUE-TYPE INFERENCE: repeated full evidence")
    text = normal(render(decisions))
    assert "effect_type = %s in the sample sheet, or effect.type in YAML" % effect in text
    assert "p_value_type = %s in the sample sheet, or pvalue.type in YAML" % pvalue in text
    assert text.count("Interpretation warning") == 2
    assert "repeated full evidence" not in text
    assert "type automatically inferred" not in text


def test_info_subsection_uses_saved_source_evidence_and_keeps_row_qc_pending(decisions):
    decisions.update(info_score_type="mach_rsq", info_score_type_source="automatic_full_dataset",
        info_score_type_evidence={"source": "external INFO column 'Rsq'", "source_rows": 1205,
            "finite_values": 1200, "unusable_values": 5, "below_minimum": 1,
            "standard_range_values": 1000, "mach_range_values": 198, "above_mach_maximum": 1,
            "minimum": -.1, "maximum": 2.5, "configured_type": "auto",
            "resolved_type": "mach_rsq", "decision_source": "automatic_full_dataset",
            "mach_range_fraction": .165, "invalid_fraction": 1 / 1200,
            "auto_mach_rsq_fraction": .01, "maximum_invalid_fraction": .1,
            "standard_info_maximum": 1.0, "mach_rsq_maximum": 2.0})
    before = deepcopy(decisions)
    text = normal(render(decisions))
    section = text.split("Imputation-quality score", 1)[1]
    assert "Selected score type : mach rsq" in section
    assert "external INFO column 'Rsq'" in section
    assert "Selected source rows : 1,205" in section
    assert "Finite numeric values : 1,200" in section
    assert "198 / 1,200 (16.50%)" in section
    assert "1.00%" in section and "10.00%" in section
    assert "normal chromosome INFO QC still applies" in section
    assert decisions == before


def test_unvisited_info_check_is_not_presented_as_resolved():
    text = normal(render({"report_failure": "External AF preparation failed"}))
    section = text.split("Imputation-quality score", 1)[1].split("Dataset stopped", 1)[0]
    assert "Not recorded" in section
    assert "Selected score type" not in section
    assert "normal chromosome INFO QC" not in section


@pytest.mark.parametrize("width", [72, 120, 160])
def test_report_uses_shared_wrapping_and_one_value_column(decisions, monkeypatch, width):
    monkeypatch.setattr("shutil.get_terminal_size", lambda fallback: os.terminal_size((width, 30)))
    settings = load_configuration().logging
    with terminal_presentation(settings):
        text = render(decisions)
    assert max(map(cell_len, text.splitlines())) <= width
    colon_columns = {cell_len(line.split(" : ", 1)[0]) for line in text.splitlines() if " : " in line}
    assert len(colon_columns) == 1
