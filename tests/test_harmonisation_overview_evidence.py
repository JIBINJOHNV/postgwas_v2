"""Overview totals retain the coverage and meaning of persisted metrics."""

from copy import deepcopy

import pytest

from postgwas.modules.harmonisation.report_steps import summarise_chromosome_evidence


def _manifest(*summaries):
    return {"dataset": {"chromosomes": {str(i): {} for i in range(1, len(summaries) + 1)},
                        "chromosome_summaries": {str(i): value for i, value in enumerate(summaries, 1)}}}


def _count(result, group, key):
    return result[group]["counts"][key]


def test_recorded_zero_is_distinct_from_absent_and_partial_evidence():
    manifest = _manifest({"stage_qc": {"info_qc": {"missing_info": 0}}}, {})
    manifest["dataset"]["chromosomes"]["X"] = {"status": "failed"}
    result = summarise_chromosome_evidence(manifest)
    assert result["coverage"] == {"expected_chromosomes": 3, "recorded_chromosomes": 2,
                                  "missing_chromosomes": ["X"], "complete": False}
    assert _count(result, "info", "missing_info") == {
        "value": 0, "expected_chromosomes": 3, "recorded_chromosomes": 1,
        "missing_chromosomes": ["2", "X"], "complete": False,
    }
    assert _count(result, "info", "rejected_missing")["value"] is None
    assert result["info"]["sources"]["source"] == []
    assert result["info"]["source_coverage"]["source"]["complete"] is False


def test_empty_manifest_never_implies_completed_checks_or_zero_counts():
    result = summarise_chromosome_evidence({"status": "OK"})
    assert result["coverage"]["complete"] is False
    assert result["coverage"]["expected_chromosomes"] == 0
    for name in ("effects", "frequency", "pvalues", "validation", "sample_size", "info", "vcf"):
        assert all(item["value"] is None for item in result[name]["counts"].values())
        assert result[name]["coverage"]["complete"] is False
    assert result["sample_size"]["ranges"] == {}


def test_sparse_manifest_keeps_explicit_completed_and_failed_chromosomes():
    manifest = {"dataset": {"completed": ["1", "2"], "failed": ["X"],
                            "chromosome_summaries": {"1": {"stage_qc": {"info_qc": {"missing_info": 0}}}}}}
    result = summarise_chromosome_evidence(manifest)
    assert result["coverage"]["expected_chromosomes"] == 3
    assert result["coverage"]["missing_chromosomes"] == ["2", "X"]
    assert _count(result, "info", "missing_info")["complete"] is False


def test_multiple_chromosomes_sum_each_metric_once_and_preserve_sources():
    first = {"stage_qc": {
        "info_qc": {"initial_variants": 100, "final_variants": 96, "missing_info": 2,
                    "rescaled_within_tolerance": 3, "rejected_out_of_range": 4,
                    "rejected_missing": 0, "source": "study_column:INFO"},
        "eaf_qc": {"strand_non_palindromic_af_comparable": 90, "strand_non_palindromic_af_discordant": 2,
                   "strand_af_discordance_action": "warn", "eaf_provenance": "internal_eaf"},
        "pval_detection_and_conversion_qc": {"variants_with_pvalues_clipped_low": 2,
                                             "variants_with_zero_pvalues_replaced": 1},
    }}
    second = deepcopy(first)
    second["stage_qc"]["info_qc"]["source"] = "fixed_cli:1"
    result = summarise_chromosome_evidence(_manifest(first, second))
    assert _count(result, "info", "removed_total")["value"] == 8
    assert _count(result, "info", "rescaled_within_tolerance")["value"] == 6
    assert _count(result, "frequency", "strand_non_palindromic_af_discordant")["value"] == 4
    assert result["info"]["sources"]["source"] == ["fixed_cli:1", "study_column:INFO"]
    assert result["info"]["source_coverage"]["source"]["complete"] is True
    # Zero-P replacements can also be low-clipped: no combined total is made.
    assert _count(result, "pvalues", "variants_with_pvalues_clipped_low")["value"] == 4
    assert _count(result, "pvalues", "variants_with_zero_pvalues_replaced")["value"] == 2
    assert "removed_total" not in result["frequency"]["counts"]


def test_stage_qc_overrides_duplicate_step_extra_and_supports_extra_only():
    summary = {"steps": [{"number": 11, "total": 16, "extra": {"missing_info": 2}},
                         {"number": 5, "total": 16, "extra": {"Neff_status": "not_computed_no_inputs"}}],
               "stage_qc": {"info_qc": {"missing_info": 3}}}
    result = summarise_chromosome_evidence(_manifest(summary))
    assert _count(result, "info", "missing_info")["value"] == 3
    assert result["sample_size"]["sources"]["Neff_status"] == ["not_computed_no_inputs"]


@pytest.mark.parametrize("invalid", [None, True, False, "3", -1, 1.5, float("nan"), float("inf")])
def test_invalid_counts_do_not_become_zero(invalid):
    result = summarise_chromosome_evidence(_manifest({"stage_qc": {"info_qc": {"missing_info": invalid}}}))
    assert _count(result, "info", "missing_info")["value"] is None
    assert _count(result, "info", "missing_info")["recorded_chromosomes"] == 0


def test_partial_se_recovery_is_not_whole_frame_reconstruction():
    result = summarise_chromosome_evidence(_manifest({"stage_qc": {
        "effect_from_z_qc": {"has_beta": True, "has_se": True, "has_zscore": True,
                             "beta_computed": False, "se_computed": True,
                             "se_recovered_from_z": 2, "se_unavailable_from_z": 1,
                             "variants_remaining": 100},
        "se_from_beta_pvalue_qc": {"se_from_pvalue": 1, "se_missing_recovered_from_pvalue": 1,
                                   "se_from_zero_p_approximation": 1},
        "z_from_beta_se_qc": {"z_source": "study_column:Z", "total_variants": 100},
    }}))
    assert _count(result, "effects", "se_recovered_from_z")["value"] == 2
    assert _count(result, "effects", "se_reconstructed_from_z")["value"] == 0
    assert _count(result, "effects", "beta_derived_from_z")["value"] == 0
    assert _count(result, "effects", "se_from_pvalue")["value"] == 1
    assert result["effects"]["sources"]["has_zscore"] == [True]
    assert "z_derived" not in result["effects"]["counts"]


def test_beta_derivation_excludes_undefined_and_requires_supporting_counts():
    known = {"stage_qc": {"effect_from_z_qc": {"beta_computed": True,
              "beta_undefined_from_z": 2, "variants_remaining": 10}}}
    unknown = {"stage_qc": {"effect_from_z_qc": {"beta_computed": True, "variants_remaining": 90}}}
    value = _count(summarise_chromosome_evidence(_manifest(known, unknown)), "effects", "beta_derived_from_z")
    assert value["value"] == 8
    assert value["complete"] is False
    assert value["missing_chromosomes"] == ["2"]


def test_z_only_reconstruction_reuses_recorded_formulas_not_trait_inference():
    result = summarise_chromosome_evidence(_manifest({"stage_qc": {
        "effect_from_z_qc": {"has_beta": False, "has_se": False,
                             "beta_computed": True, "se_computed": True,
                             "variants_remaining": 80, "invalid_reconstructed_statistics": 0,
                             "effect_estimate_beta_formula": "saved beta formula",
                             "effect_estimate_se_formula": "saved SE formula"},
        "sample_size_qc": {"trait_type": "auto", "Neff_status": "calculated_from_case_control"},
    }}))
    assert _count(result, "effects", "beta_derived_from_z")["value"] == 80
    assert _count(result, "effects", "se_reconstructed_from_z")["value"] == 80
    assert result["effects"]["sources"]["effect_estimate_beta_formula"] == ["saved beta formula"]
    assert result["sample_size"]["sources"]["trait_type"] == ["auto"]


def test_reconstruction_boolean_is_not_a_recorded_zero_invalid_count():
    result = summarise_chromosome_evidence(_manifest({"stage_qc": {"effect_from_z_qc": {
        "has_beta": False, "has_se": False, "beta_computed": True, "se_computed": True,
        "variants_remaining": 10, "invalid_reconstructed_statistics": False,
    }}}))
    assert _count(result, "effects", "beta_derived_from_z")["value"] is None
    assert _count(result, "effects", "se_reconstructed_from_z")["value"] is None


def test_concordance_retention_is_check_local_and_not_inferred_from_fail_or_off():
    summaries = [{"stage_qc": {"effect_statistics_validation_qc": {
        "removed_total": 3, "variants_removed_invalid_beta_se": 1,
        "basic_matched_by_reason": {"se_null": 1, "beta_zero": 2},
        "beta_se_z_concordance": {"checked": 90, "discordant": 3, "removed": 0, "status": "ran", "action": "warn"},
        "z_pval_concordance": {"checked": 90, "discordant": 2, "removed": 2, "status": "ran", "action": "reject"},
    }}}, {"stage_qc": {"effect_statistics_validation_qc": {
        "beta_se_z_concordance": {"checked": 10, "discordant": 1, "removed": 0, "status": "ran", "action": "fail"},
        "z_pval_concordance": {"checked": 0, "discordant": 0, "removed": 0, "status": "off", "action": "off"},
    }}}]
    result = summarise_chromosome_evidence(_manifest(*summaries))
    beta = result["validation"]["beta_se_z"]["counts"]
    assert beta["discordant"]["value"] == 4
    assert beta["retained_discordant"]["value"] == 3
    assert beta["retained_discordant"]["complete"] is False
    assert result["validation"]["z_p"]["counts"]["retained_discordant"]["value"] == 0
    assert _count(result, "validation", "removed_total")["value"] == 3


def test_sample_size_ranges_use_recorded_extrema_and_preserve_missing_plan():
    manifest = _manifest({"stage_qc": {"sample_size_qc": {"sample_size_stats": {
        "Neff": {"min": 100, "max": 1000, "median": 900, "mean": 500},
        "Ncase": {"median": 99}, "ncontrol_col": "missing",
    }, "missing_sample_size": {"action": "median", "variants_affected": 1, "variants_imputed": 1}}}},
        {"stage_qc": {"sample_size_qc": {"sample_size_stats": {
            "Neff": {"min": 200, "max": 2000, "median": 300, "mean": 400},
            "Ncase": {"min": 20, "max": 100},
        }}}})
    manifest["dataset"]["missing_sample_size"] = {"action": "median", "fill_values": {"Ncase": 100}}
    original = deepcopy(manifest)
    result = summarise_chromosome_evidence(manifest)
    assert result["sample_size"]["ranges"]["Neff"] == {
        "min": 100, "max": 2000, "coverage": {"expected_chromosomes": 2,
        "recorded_chromosomes": 2, "missing_chromosomes": [], "complete": True},
    }
    assert result["sample_size"]["ranges"]["Ncase"]["coverage"]["complete"] is False
    assert "ncontrol_col" not in result["sample_size"]["ranges"]
    result["sample_size"]["missing_plan"]["fill_values"]["Ncase"] = 999
    assert manifest == original


def test_vcf_counts_do_not_reuse_step_rows_as_vcf_outputs_or_sum_liftover_losses():
    result = summarise_chromosome_evidence(_manifest({"steps": [
        {"number": 15, "total": 16, "rows_in": 100, "rows_out": 100},
        {"number": 16, "total": 16, "rows_out": 100, "extra": {"liftover_accounting": {
            "input": 95, "rejected": 1, "swap_excluded": 2, "final": 92, "swap_policy": "exclude"}}},
    ], "qc_actions": [{"step": "16 bcftools_annotate_liftover", "check": "NORM", "before": 98, "after": 95, "removed": 3}]}))
    counts = result["vcf"]["counts"]
    assert counts["adapter_input"]["value"] == 100
    assert counts["adapter_output"]["value"] == 98
    assert counts["adapter_removed"]["value"] == 2
    assert counts["normalization_output"]["value"] == 95
    assert counts["liftover_rejected"]["value"] == 1
    assert counts["liftover_swap_excluded"]["value"] == 2
    assert counts["target_output"]["value"] == 92
    missing = summarise_chromosome_evidence(_manifest({"steps": [
        {"number": 15, "total": 16, "rows_in": 100, "rows_out": 100},
        {"number": 16, "total": 16, "rows_out": 100},
    ]}))
    assert _count(missing, "vcf", "adapter_output")["value"] is None
    assert _count(missing, "vcf", "target_output")["value"] is None
