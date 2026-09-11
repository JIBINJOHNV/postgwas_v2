"""Scientific contracts for effect type and OR standard-error scale."""

import math
from pathlib import Path

import polars as pl
import pytest

from postgwas.core.pipeline_logging import PipelineLogger
from postgwas.modules.harmonisation.effect_type import (
    EffectTypeError,
    detect_effect_type,
    detect_or_standard_error_scale,
    harmonise_effect_estimates,
)
from postgwas.modules.harmonisation.policies import PolicyError, default_policies
from postgwas.modules.harmonisation.rejects import (
    RejectCollector,
    SOURCE_INPUT_ROW_COLUMN,
)
from postgwas.modules.harmonisation.study_properties import (
    StudyPropertyError,
    resolve_study_properties,
)


def _detect(values):
    return detect_effect_type(
        pl.DataFrame({"EFFECT": values}), "EFFECT", default_policies()
    )


def test_or_detection_defaults_are_conservative_and_validated():
    policies = default_policies()
    assert policies.get("effect.or_detection_max_non_positive_fraction") == 0.005
    assert policies.get("effect.or_detection_require_median_range") == [0.8, 1.25]
    assert policies.get("effect.or_non_positive") == "reject"
    assert policies.get("effect.se_scale") == "auto"
    assert policies.get("effect.se_scale_min_variants") == 100
    assert policies.get("effect.se_scale_min_agreement_fraction") == 0.95
    assert policies.get("effect.se_scale_min_agreement_margin") == 0.1
    assert policies.get("effect.cast_strict") is False
    assert "effect.or_detection_negative_fraction" not in policies
    with pytest.raises(PolicyError, match="or_detection_max_non_positive_fraction"):
        policies.with_overrides({
            "effect.or_detection_max_non_positive_fraction": 1.01,
        })
    with pytest.raises(PolicyError, match="se_scale_min_variants"):
        policies.with_overrides({"effect.se_scale_min_variants": 0})
    with pytest.raises(PolicyError, match="se_scale_min_agreement_fraction"):
        policies.with_overrides({"effect.se_scale_min_agreement_fraction": 1.01})


def test_direct_effect_harmonisation_nulls_and_counts_unparseable_values():
    result, qc, mapping = harmonise_effect_estimates(
        "1",
        pl.DataFrame({"BETA": ["0.1", "not-a-number", "-0.2"]}),
        {"beta_or_col": "BETA", "se_col": None},
        policies=default_policies(),
        decision={"effect_type": "beta"},
    )

    assert result["BETA"].to_list() == [0.1, None, -0.2]
    assert qc["effect_cast_strict"] is False
    assert qc["effect_unparseable_values"] == 1
    assert mapping["beta_col"] == "BETA"


def test_explicit_strict_effect_cast_still_fails_on_unparseable_values():
    policies = default_policies().with_overrides({"effect.cast_strict": True})

    with pytest.raises(
        EffectTypeError,
        match=(
            "Chromosome 1 effect column 'BETA' contains 1 non-missing value.*"
            "effect.cast_strict=true"
        ),
    ):
        harmonise_effect_estimates(
            "1",
            pl.DataFrame({"BETA": ["0.1", "not-a-number"]}),
            {"beta_or_col": "BETA", "se_col": None},
            policies=policies,
            decision={"effect_type": "beta"},
        )


def _or_se_frame(scale, *, include_z=True, include_p=False):
    odds_ratios = [2.0, 0.5] * 60
    log_se = [0.15] * len(odds_ratios)
    supplied_se = (
        log_se
        if scale == "log_odds"
        else [se * odds_ratio for se, odds_ratio in zip(log_se, odds_ratios)]
    )
    z_scores = [
        math.log(odds_ratio) / se
        for odds_ratio, se in zip(odds_ratios, log_se)
    ]
    columns = {"OR": odds_ratios, "SE": supplied_se}
    if include_z:
        columns["Z"] = z_scores
    if include_p:
        columns["P"] = [
            math.erfc(abs(z_score) / math.sqrt(2.0)) for z_score in z_scores
        ]
    return pl.DataFrame(columns)


@pytest.mark.parametrize("scale", ["log_odds", "as_given"])
def test_or_se_scale_is_decided_once_from_supplied_z(scale):
    detected, evidence = detect_or_standard_error_scale(
        _or_se_frame(scale),
        "OR",
        "SE",
        z_col="Z",
        policies=default_policies(),
    )

    assert detected == scale
    assert evidence["comparison_source"] == "supplied_z"
    assert evidence["informative_variants"] == 120
    assert evidence[
        "log_odds_agreement_fraction"
        if scale == "log_odds"
        else "raw_odds_ratio_agreement_fraction"
    ] == pytest.approx(1.0)


@pytest.mark.parametrize(
    "decision",
    [None, {"effect_type": "odds_ratio", "se_scale": "as_given"}],
    ids=["automatically_detected_or", "explicit_or"],
)
def test_raw_or_scale_se_uses_original_or_when_effect_column_is_named_beta(
    decision,
):
    result, qc, mapping = harmonise_effect_estimates(
        "1",
        pl.DataFrame({
            "beta": [0.5, 1.0, 2.0],
            "SE": [0.2, 0.1, 0.3],
        }),
        {"beta_or_col": "beta", "se_col": "SE"},
        policies=default_policies().with_overrides({
            "effect.se_scale": "as_given",
        }),
        decision=decision,
    )

    assert result["beta"].to_list() == pytest.approx([
        math.log(0.5), 0.0, math.log(2.0),
    ])
    assert result["SE"].to_list() == pytest.approx([0.4, 0.1, 0.15])
    assert mapping["beta_col"] == "beta"
    assert mapping["beta_or_col"] == "beta"
    assert qc["se_rescaled_by_or"] is True
    assert qc["beta_output_name_collision_avoided"] is False


def test_log_odds_se_is_unchanged_when_or_column_is_named_beta():
    result, qc, _mapping = harmonise_effect_estimates(
        "1",
        pl.DataFrame({
            "beta": [0.5, 1.0, 2.0],
            "SE": [0.4, 0.1, 0.15],
        }),
        {"beta_or_col": "beta", "se_col": "SE"},
        policies=default_policies(),
        decision={"effect_type": "odds_ratio", "se_scale": "log_odds"},
    )

    assert result["beta"].to_list() == pytest.approx([
        math.log(0.5), 0.0, math.log(2.0),
    ])
    assert result["SE"].to_list() == pytest.approx([0.4, 0.1, 0.15])
    assert qc["se_rescaled_by_or"] is False


def test_unrelated_beta_column_is_preserved_with_collision_safe_output():
    result, qc, mapping = harmonise_effect_estimates(
        "1",
        pl.DataFrame({"OR": [2.0], "beta": [99.0], "SE": [0.3]}),
        {"beta_or_col": "OR", "se_col": "SE"},
        policies=default_policies(),
        decision={"effect_type": "odds_ratio", "se_scale": "as_given"},
    )

    assert result.item(0, "beta") == pytest.approx(99.0)
    assert result.item(0, "beta_") == pytest.approx(math.log(2.0))
    assert result.item(0, "SE") == pytest.approx(0.15)
    assert mapping["beta_col"] == "beta_"
    assert mapping["beta_or_col"] == "beta_"
    assert qc["beta_output_column"] == "beta_"
    assert qc["beta_output_name_collision_avoided"] is True


def test_effect_and_se_cannot_map_to_the_same_input_column():
    with pytest.raises(
        EffectTypeError,
        match=(
            "maps the effect size and standard error to the same input column "
            "'STAT'.*Correct effect_column and standard_error_column"
        ),
    ):
        harmonise_effect_estimates(
            "1",
            pl.DataFrame({"STAT": [1.0]}),
            {"beta_or_col": "STAT", "se_col": "STAT"},
            policies=default_policies(),
            decision={"effect_type": "beta"},
        )


@pytest.mark.parametrize("pvalue_type", ["raw", "neglog10"])
def test_or_se_scale_falls_back_to_two_sided_pvalue_when_z_is_absent(pvalue_type):
    frame = _or_se_frame("as_given", include_z=False, include_p=True)
    if pvalue_type == "neglog10":
        frame = frame.with_columns((-pl.col("P").log10()).alias("P"))
    detected, evidence = detect_or_standard_error_scale(
        frame,
        "OR",
        "SE",
        pvalue_col="P",
        pvalue_type=pvalue_type,
        policies=default_policies(),
    )

    assert detected == "as_given"
    assert evidence["comparison_source"] == "reported_two_sided_pvalue"
    assert evidence["informative_variants"] == 120


def test_reported_raw_p_zero_is_not_used_as_exact_se_scale_evidence():
    frame = _or_se_frame("log_odds", include_z=False).with_columns(
        pl.lit(0.0).alias("P")
    )
    detected, evidence = detect_or_standard_error_scale(
        frame,
        "OR",
        "SE",
        pvalue_col="P",
        pvalue_type="raw",
        policies=default_policies(),
    )

    assert detected is None
    assert evidence["comparable_variants"] == 0
    assert evidence["status"] == "insufficient_informative_variants"


def test_auto_or_se_scale_is_recorded_as_one_study_wide_decision():
    decisions = resolve_study_properties(
        _or_se_frame("log_odds"),
        {"beta_or_col": "OR", "se_col": "SE", "imp_z_col": "Z"},
        policies=default_policies().with_overrides({"effect.type": "odds_ratio"}),
    )

    assert decisions["se_scale"] == "log_odds"
    assert decisions["se_scale_source"] == "dataset_z_pvalue_cross_check"
    assert decisions["se_scale_detected"] == "log_odds"
    assert decisions["se_scale_matches_declaration"] is None
    assert decisions["se_scale_evidence"]["informative_variants"] == 120


def test_dataset_beta_z_guard_allows_exact_one_percent_and_reports(capsys):
    beta = [0.2] * 100
    z_score = [2.0] * 99 + [-2.0]

    decisions = resolve_study_properties(
        pl.DataFrame({"BETA": beta, "Z": z_score}),
        {"beta_or_col": "BETA", "se_col": None, "imp_z_col": "Z"},
        policies=default_policies().with_overrides({"effect.type": "beta"}),
    )

    assessment = decisions["beta_z_sign_consistency"]
    assert assessment["comparable_missing_se_rows"] == 100
    assert assessment["discordant_rows"] == 1
    assert assessment["discordant_fraction"] == pytest.approx(0.01)
    assert assessment["maximum_fraction"] == pytest.approx(0.01)
    assert assessment["decision"] == "continue_and_reject"
    screen = capsys.readouterr().out
    assert "BETA/Z SIGN DISCORDANCE" in screen
    assert "1 of 100 comparable missing-SE rows (1.0000%)" in screen
    assert "Every affected variant that reaches chromosome step 06" in screen
    assert "beta_z_sign_discordant" in screen


def test_dataset_beta_z_guard_stops_strictly_above_one_percent(capsys):
    beta = [0.2] * 100
    z_score = [2.0] * 98 + [-2.0, -2.0]

    with pytest.raises(
        StudyPropertyError,
        match=(
            r"2 discordant row\(s\) among 100 comparable row\(s\).*"
            r"2\.0000%.*max_beta_z_sign_mismatch_fraction=0\.01.*"
            r"No chromosome processing was started"
        ),
    ):
        resolve_study_properties(
            pl.DataFrame({"BETA": beta, "Z": z_score}),
            {"beta_or_col": "BETA", "se_col": None, "imp_z_col": "Z"},
            policies=default_policies().with_overrides({
                "effect.type": "beta",
            }),
        )

    screen = capsys.readouterr().out
    assert "BETA/Z SIGN DISCORDANCE LIMIT EXCEEDED" in screen
    assert "No chromosome processing was started" in screen


def test_explicit_beta_z_fail_policy_stops_on_any_dataset_mismatch():
    with pytest.raises(
        StudyPropertyError,
        match=r"beta_z_sign_mismatch='fail' stops on any discordant row",
    ):
        resolve_study_properties(
            pl.DataFrame({"BETA": [0.2, 0.2], "Z": [2.0, -2.0]}),
            {"beta_or_col": "BETA", "se_col": None, "imp_z_col": "Z"},
            policies=default_policies().with_overrides({
                "effect.type": "beta",
                "effect_from_z.beta_z_sign_mismatch": "fail",
                "effect_from_z.max_beta_z_sign_mismatch_fraction": 1.0,
            }),
        )


def test_beta_z_guard_denominator_excludes_unusable_and_populated_se_rows():
    decisions = resolve_study_properties(
        pl.DataFrame({
            "BETA": [0.2, 0.2, 0.2, None, 0.2],
            "Z": [2.0, -2.0, None, -2.0, -2.0],
            "SE": [None, None, None, None, 0.1],
        }),
        {"beta_or_col": "BETA", "se_col": "SE", "imp_z_col": "Z"},
        policies=default_policies().with_overrides({
            "effect.type": "beta",
            "effect_from_z.max_beta_z_sign_mismatch_fraction": 0.5,
        }),
    )

    assessment = decisions["beta_z_sign_consistency"]
    assert assessment["comparable_missing_se_rows"] == 2
    assert assessment["discordant_rows"] == 1
    assert assessment["discordant_fraction"] == pytest.approx(0.5)


def test_dataset_sign_guard_uses_log_odds_direction_for_odds_ratios():
    decisions = resolve_study_properties(
        pl.DataFrame({
            "OR": [2.0, 0.5, 2.0, 0.5, 1.0],
            "Z": [2.0, -2.0, -2.0, 2.0, 3.0],
        }),
        {"beta_or_col": "OR", "se_col": None, "imp_z_col": "Z"},
        policies=default_policies().with_overrides({
            "effect.type": "odds_ratio",
            "effect_from_z.max_beta_z_sign_mismatch_fraction": 0.6,
        }),
    )

    assessment = decisions["beta_z_sign_consistency"]
    assert assessment["effect_type"] == "odds_ratio"
    assert assessment["comparable_missing_se_rows"] == 5
    assert assessment["discordant_rows"] == 3
    assert assessment["discordant_fraction"] == pytest.approx(0.6)


def test_auto_or_se_scale_fails_before_split_when_evidence_is_insufficient():
    frame = _or_se_frame("log_odds").head(50)

    with pytest.raises(StudyPropertyError, match="No chromosome processing was started"):
        resolve_study_properties(
            frame,
            {"beta_or_col": "OR", "se_col": "SE", "imp_z_col": "Z"},
            policies=default_policies().with_overrides({"effect.type": "odds_ratio"}),
        )


def test_explicit_or_se_scale_is_cross_checked_but_remains_authoritative(capsys):
    decisions = resolve_study_properties(
        _or_se_frame("as_given"),
        {"beta_or_col": "OR", "se_col": "SE", "imp_z_col": "Z"},
        policies=default_policies().with_overrides({
            "effect": {"type": "odds_ratio", "se_scale": "log_odds"},
        }),
    )

    assert decisions["se_scale"] == "log_odds"
    assert decisions["se_scale_detected"] == "as_given"
    assert decisions["se_scale_matches_declaration"] is False
    assert "DECLARATION MISMATCH" in capsys.readouterr().out


def test_beta_study_does_not_run_or_se_scale_detection():
    decisions = resolve_study_properties(
        pl.DataFrame({"BETA": [-0.1, 0.0, 0.1], "SE": [0.1, 0.1, 0.1]}),
        {"beta_or_col": "BETA", "se_col": "SE"},
        policies=default_policies().with_overrides({"effect.type": "beta"}),
    )

    assert decisions["se_scale"] is None
    assert decisions["se_scale_source"] == "not_applicable_to_beta"


def test_or_accepts_at_most_half_percent_combined_non_positive_contamination():
    decision, evidence = _detect([0.0] + [1.0] * 199)

    assert decision == "odds_ratio"
    assert evidence["non_positive_value_count"] == 1
    assert evidence["non_positive_fraction"] == pytest.approx(0.005)


def test_beta_requires_high_non_positive_fraction_and_non_or_median():
    decision, evidence = _detect([-0.1, 0.0] + [0.05] * 198)

    assert decision == "beta"
    assert evidence["non_positive_fraction"] == pytest.approx(0.01)
    assert evidence["median_range_confirmed"] is False


def test_all_positive_beta_is_ambiguous_instead_of_log_transformed():
    with pytest.raises(EffectTypeError, match="all-positive beta"):
        _detect([0.05] * 200)


def test_or_like_median_with_excess_non_positive_values_is_conflicting():
    with pytest.raises(EffectTypeError, match="conflicting evidence"):
        _detect([-0.1, 0.0] + [1.0] * 198)


def test_missing_values_do_not_dilute_the_non_positive_fraction():
    decision, evidence = _detect([-0.1] * 2 + [0.05] * 198 + [None] * 800)

    assert decision == "beta"
    assert evidence["n_usable"] == 200
    assert evidence["non_positive_fraction"] == pytest.approx(0.01)


def test_applied_auto_effect_type_warns_on_screen_and_in_canonical_log(
    tmp_path, capsys,
):
    policies = default_policies()
    logger = PipelineLogger(
        "study", "dataset", str(tmp_path), policies=policies,
        screen_level="ERROR",
    )
    try:
        decisions = resolve_study_properties(
            pl.DataFrame({"EFFECT": [0.9, 1.0, 1.1]}),
            {"beta_or_col": "EFFECT"},
            policies=policies,
            logger=logger,
        )
    finally:
        logger.close()

    warning = decisions["effect_type_inference_warning"]
    assert decisions["effect_type"] == "odds_ratio"
    assert decisions["effect_type_source"] == "detector"
    assert "AUTOMATIC EFFECT-TYPE INFERENCE" in warning
    assert "BETA = ln(effect)" in warning
    assert "sample-sheet effect_type" in warning
    screen = " ".join(capsys.readouterr().out.split())
    assert "Effect type automatically inferred" in screen
    assert "Evidence used" in screen
    assert "OR rule: at most" in screen
    assert "Continuing with the inferred type" in screen
    log_text = Path(logger.log_path).read_text(encoding="utf-8")
    assert "AUTOMATIC EFFECT-TYPE INFERENCE" in log_text
    assert "BETA = ln(effect)" in log_text
    assert "sample-sheet" in log_text
    assert "effect_type field" in log_text


def test_inference_screen_uses_observed_counts_and_overridden_thresholds(capsys):
    policies = default_policies().with_overrides({
        "effect.or_detection_max_non_positive_fraction": 0.02,
        "effect.or_detection_require_median_range": [0.7, 1.4],
    })
    decisions = resolve_study_properties(
        pl.DataFrame({"EFFECT": [0.0] + [0.9999] * 99 + [None, float("inf")]}),
        {"beta_or_col": "EFFECT"}, policies=policies,
    )
    assert decisions["effect_type"] == "odds_ratio"
    screen = " ".join(capsys.readouterr().out.split())
    assert "1 / 100 (1.000%); OR rule: at most 2%" in screen
    assert "0.9999; OR rule: 0.7 to 1.4 inclusive" in screen
    assert "effect_type = odds_ratio" in screen
    assert "BETA = ln(effect)" in screen


def test_explicit_sample_sheet_effect_type_does_not_emit_inference_warning(capsys):
    decisions = resolve_study_properties(
        pl.DataFrame({"EFFECT": [0.9, 1.0, 1.1]}),
        {"beta_or_col": "EFFECT", "declared_effect_type": "odds_ratio"},
        policies=default_policies().with_overrides({"effect.type": "odds_ratio"}),
    )

    assert decisions["effect_type_source"] == "sample_sheet"
    assert "effect_type_inference_warning" not in decisions
    assert "AUTOMATIC EFFECT-TYPE INFERENCE" not in capsys.readouterr().out


def test_inferred_beta_warning_states_that_values_are_not_log_transformed(capsys):
    decisions = resolve_study_properties(
        pl.DataFrame({"EFFECT": [-0.1, 0.0] + [0.05] * 198}),
        {"beta_or_col": "EFFECT"},
        policies=default_policies(),
    )

    warning = decisions["effect_type_inference_warning"]
    assert decisions["effect_type"] == "beta"
    assert "treated directly as beta coefficients" in warning
    screen = " ".join(capsys.readouterr().out.split())
    assert "Beta coefficient" in screen
    assert "treated directly as beta coefficients" in screen


def test_confirmed_or_rejects_tolerated_non_positive_rows_before_log(tmp_path):
    frame = pl.DataFrame({"OR": [0.0] + [1.0] * 199}).with_row_index(
        SOURCE_INPUT_ROW_COLUMN, offset=1,
    )
    collector = RejectCollector(
        source_snapshot=frame,
        logger=None,
        out_path=str(tmp_path / "rejected.tsv"),
        delimiter="\t",
        compress=False,
    )

    result, qc, mapping = harmonise_effect_estimates(
        "1",
        frame,
        {"beta_or_col": "OR", "se_col": None},
        policies=default_policies(),
        rejects=collector,
    )

    assert result.height == 199
    assert result["beta"].to_list() == [0.0] * 199
    assert collector.counts() == {"effect_non_positive_or": 1}
    assert qc["or_non_positive_action"] == "reject"
    assert "effect_type_inference_warning" in qc
    assert mapping["beta_col"] == "beta"
