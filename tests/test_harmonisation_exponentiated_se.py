"""Decode exponentiated uncertainty without changing supplied study values."""

import math

import polars as pl
import pytest

from postgwas.modules.harmonisation.effect_type import (
    EffectTypeError,
    detect_or_standard_error_scale,
    harmonise_effect_estimates,
)
from postgwas.modules.harmonisation.policies import default_policies


def _frame(scale, rows=120):
    odds_ratios = [0.5, 2.0] * (rows // 2) + ([0.5] if rows % 2 else [])
    log_se = [0.15] * rows
    stored = {
        "log_odds": log_se,
        "as_given": [value * ratio for value, ratio in zip(log_se, odds_ratios)],
        "exponentiated_log_odds": [math.exp(value) for value in log_se],
    }[scale]
    z = [math.log(ratio) / 0.15 for ratio in odds_ratios]
    return pl.DataFrame({
        "OR": odds_ratios, "SE": stored, "Z": z,
        "P": [math.erfc(abs(value) / math.sqrt(2.0)) for value in z],
    })


@pytest.mark.parametrize("scale", ["log_odds", "as_given", "exponentiated_log_odds"])
@pytest.mark.parametrize("comparison", ["z", "raw_p", "neglog10_p"])
def test_auto_scale_distinguishes_all_supported_uncertainty_representations(scale, comparison):
    frame = _frame(scale)
    if comparison == "neglog10_p":
        frame = frame.with_columns((-pl.col("P").log10()).alias("P"))
    original = frame.clone()
    detected, evidence = detect_or_standard_error_scale(
        frame, "OR", "SE", z_col="Z" if comparison == "z" else None,
        pvalue_col="P", pvalue_type="neglog10" if comparison == "neglog10_p" else "raw",
        policies=default_policies(),
    )
    assert detected == scale
    assert evidence["informative_variants"] == 120
    prefix = "raw_odds_ratio" if scale == "as_given" else scale
    assert evidence[prefix + "_matches"] == 120
    assert evidence[prefix + "_agreement_fraction"] == pytest.approx(1.0)
    assert frame.equals(original)


def test_reported_nmosd_example_decodes_into_a_collision_safe_working_column():
    frame = pl.DataFrame({
        "OR": [0.893986], "standard_error": [1.130373],
        "standard_error_": [99.0], "beta": [7.0],
    })
    original = frame.clone()
    result, qc, mapping = harmonise_effect_estimates(
        "1", frame,
        {"beta_or_col": "OR", "se_col": "standard_error", "imp_z_col": "standard_error__"},
        policies=default_policies(),
        decision={"effect_type": "odds_ratio", "se_scale": "exponentiated_log_odds"},
    )
    destination = mapping["se_col"]
    assert destination not in {*frame.columns, "standard_error__"}
    assert result.select(frame.columns).equals(original)
    assert frame.equals(original)
    assert result[destination].to_list() == pytest.approx([0.1225476667526027])
    beta = result[mapping["beta_col"]].item()
    decoded_p = math.erfc(abs(beta / result[destination].item()) / math.sqrt(2.0))
    assert beta == pytest.approx(-0.11206516388649862)
    assert decoded_p == pytest.approx(0.360477, abs=3e-6)
    assert qc["se_exponentiated_decoded"] is True
    assert qc["se_values_decoded"] == 1
    assert qc["se_input_column"] == "standard_error"
    assert qc["se_output_column"] == destination


@pytest.mark.parametrize("invalid", ["0", "-1", "1", "inf", "nan", "1e309", "not-a-number"])
def test_selected_exponentiated_scale_refuses_noninvertible_provided_values(invalid):
    frame = pl.DataFrame({"OR": [2.0], "SE": [invalid]})
    original = frame.clone()
    with pytest.raises(EffectTypeError, match=r"SE.*finite.*greater than 1"):
        harmonise_effect_estimates(
            "1", frame, {"beta_or_col": "OR", "se_col": "SE"},
            policies=default_policies(),
            decision={"effect_type": "odds_ratio", "se_scale": "exponentiated_log_odds"},
        )
    assert frame.equals(original)


def test_missing_encoded_values_remain_missing_and_large_finite_values_decode():
    frame = pl.DataFrame({
        "OR": [2.0, 0.5, 1.0], "SE": [math.exp(0.15), None, 1.7976931348623157e308],
    })
    result, qc, mapping = harmonise_effect_estimates(
        "1", frame, {"beta_or_col": "OR", "se_col": "SE"},
        policies=default_policies(),
        decision={"effect_type": "odds_ratio", "se_scale": "exponentiated_log_odds"},
    )
    assert result[mapping["se_col"]].to_list() == pytest.approx([0.15, None, 709.782712893384])
    assert result["SE"].equals(frame["SE"])
    assert qc["se_values_decoded"] == 2


@pytest.mark.parametrize("rows,expected", [(99, None), (100, "exponentiated_log_odds")])
def test_exponentiated_auto_scale_keeps_the_existing_minimum_evidence_floor(rows, expected):
    detected, evidence = detect_or_standard_error_scale(
        _frame("exponentiated_log_odds", rows), "OR", "SE", z_col="Z",
        policies=default_policies(),
    )
    assert detected == expected
    assert evidence["informative_variants"] == rows


def test_null_effect_and_censored_p_do_not_confirm_exponentiated_se():
    frame = _frame("exponentiated_log_odds", 99).drop("Z")
    extra = pl.DataFrame({"OR": [1.0, 2.0, None], "SE": [math.exp(0.15)] * 3, "P": [1.0, 0.0, 0.1]})
    detected, evidence = detect_or_standard_error_scale(
        pl.concat([frame, extra]), "OR", "SE", pvalue_col="P", pvalue_type="raw",
        policies=default_policies(),
    )
    assert detected is None
    assert evidence["informative_variants"] == 99


def test_invalid_decoded_candidate_never_confirms_exponentiated_scale():
    frame = _frame("log_odds")
    detected, evidence = detect_or_standard_error_scale(
        frame, "OR", "SE", z_col="Z", policies=default_policies(),
    )
    assert detected == "log_odds"
    assert evidence["exponentiated_log_odds_matches"] == 0
    assert evidence["exponentiated_log_odds_agreement_fraction"] == 0.0


@pytest.mark.parametrize("encoded_rows,expected", [(95, "exponentiated_log_odds"), (94, None), (50, None)])
def test_mixed_scale_evidence_uses_a_common_denominator_and_existing_agreement_gate(encoded_rows, expected):
    frame = _frame("exponentiated_log_odds", 100)
    ordinary_z = frame["OR"].log() / frame["SE"]
    frame = frame.with_columns(pl.Series("Z", [
        frame["Z"][i] if i < encoded_rows else ordinary_z[i] for i in range(100)
    ]))
    detected, evidence = detect_or_standard_error_scale(
        frame, "OR", "SE", z_col="Z", policies=default_policies(),
    )
    assert detected == expected
    assert evidence["informative_variants"] == 100
    assert evidence["exponentiated_log_odds_matches"] == encoded_rows


def test_equal_best_candidates_remain_ambiguous_even_with_a_zero_margin_override():
    frame = _frame("exponentiated_log_odds", 100)
    ordinary_z = frame["OR"].log() / frame["SE"]
    frame = frame.with_columns(pl.Series("Z", [
        frame["Z"][i] if i < 50 else ordinary_z[i] for i in range(100)
    ]))
    policies = default_policies().with_overrides({
        "effect.se_scale_min_agreement_fraction": 0.5,
        "effect.se_scale_min_agreement_margin": 0.0,
    })
    detected, evidence = detect_or_standard_error_scale(
        frame, "OR", "SE", z_col="Z", policies=policies,
    )
    assert detected is None
    assert evidence["status"] == "ambiguous_or_inconsistent"


def test_ambiguous_supplied_z_is_not_overruled_by_pvalue_evidence():
    frame = _frame("exponentiated_log_odds").with_columns(pl.lit(0.0).alias("Z"))
    detected, evidence = detect_or_standard_error_scale(
        frame, "OR", "SE", z_col="Z", pvalue_col="P", pvalue_type="raw",
        policies=default_policies(),
    )
    assert detected is None
    assert evidence["comparison_source"] == "supplied_z"
    assert evidence["pvalue_attempt"] is None


@pytest.mark.parametrize("scale,ratio,stored,true_se", [
    ("log_odds", 2.0, 1.5, 1.5),
    ("as_given", 8.0, 1.2, 0.15),
])
def test_a_stored_se_above_one_is_not_itself_evidence_of_exponentiation(scale, ratio, stored, true_se):
    frame = pl.DataFrame({"OR": [ratio] * 100, "SE": [stored] * 100,
                          "Z": [math.log(ratio) / true_se] * 100})
    detected, evidence = detect_or_standard_error_scale(
        frame, "OR", "SE", z_col="Z", policies=default_policies(),
    )
    assert detected == scale
    assert evidence["exponentiated_log_odds_matches"] == 0


@pytest.mark.parametrize("shared_rows,expected", [(85, "log_odds"), (86, None)])
def test_agreement_margin_is_inclusive_and_uses_the_best_competing_scale(shared_rows, expected):
    ratios = [1.005] * shared_rows + [2.0] * (100 - shared_rows)
    z = [math.log(ratio) / 1.1 for ratio in ratios]
    z[-5:] = [100.0] * 5
    frame = pl.DataFrame({"OR": ratios, "SE": [1.1] * 100, "Z": z})
    detected, evidence = detect_or_standard_error_scale(
        frame, "OR", "SE", z_col="Z", policies=default_policies(),
    )
    assert evidence["informative_variants"] == 100
    assert evidence["log_odds_matches"] == 95
    assert evidence["raw_odds_ratio_matches"] == shared_rows
    assert detected == expected


def test_all_null_effects_provide_no_scale_confirmation():
    frame = pl.DataFrame({"OR": [1.0] * 100, "SE": [math.exp(0.15)] * 100,
                          "Z": [0.0] * 100, "P": [1.0] * 100})
    detected, evidence = detect_or_standard_error_scale(
        frame, "OR", "SE", z_col="Z", pvalue_col="P", pvalue_type="raw",
        policies=default_policies(),
    )
    assert detected is None
    assert evidence["informative_variants"] == 0
