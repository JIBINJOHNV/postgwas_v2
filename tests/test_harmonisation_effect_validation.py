"""Scientific regression tests for chromosome effect-statistic concordance."""

import inspect
import math

import polars as pl
import pytest

from postgwas.core.statistics import normal_z_magnitude_from_ln_p
from postgwas.modules.harmonisation import effect_validation
from postgwas.modules.harmonisation.effect_validation import (
    EffectValidationError,
    FIELD_COLUMN_KEYS,
    FIELD_REJECT_REASONS,
    final_completeness_check,
    validate_effect_statistics,
)
from postgwas.modules.harmonisation.policies import PolicyError, default_policies
from postgwas.modules.harmonisation.p_values import harmonise_p_values
from postgwas.modules.harmonisation.rejects import (
    REASONS,
    RejectCollector,
    SOURCE_INPUT_ROW_COLUMN,
)
from postgwas.modules.harmonisation.z_score import (
    derive_z_score_from_effect_and_standard_error,
)


COLUMNS = {
    "beta_col": "BETA",
    "se_col": "SE",
    "imp_z_col": "Z",
    "pval_col": "P",
}


def _frame():
    return pl.DataFrame({
        "BETA": [0.2, 0.2, 0.2],
        "SE": [0.1, 0.1, 0.1],
        "Z": [2.0, 2.02, 5.0],
        "P": [0.0455002639, 0.043383, 0.0455002639],
    })


def test_concordance_checks_are_enabled_as_non_destructive_warnings_by_default():
    policies = default_policies()

    result, qc, _ = validate_effect_statistics(
        "1", _frame(), dict(COLUMNS), policies=policies,
    )

    assert policies.get("validation.beta_se_z_concordance") == "warn"
    assert policies.get("validation.z_pval_concordance") == "warn"
    assert result.height == 3
    assert qc["beta_se_z_concordance"]["checked"] == 3
    assert qc["beta_se_z_concordance"]["discordant"] == 1
    assert qc["z_pval_concordance"]["checked"] == 3
    assert qc["z_pval_concordance"]["discordant"] == 1


def test_z_pvalue_concordance_uses_exact_log_for_float64_underflow():
    policies = default_policies()
    z = normal_z_magnitude_from_ln_p(
        [-400.0 * math.log(10.0)], tail=2,
    )[0]
    frame, _p_qc, columns = harmonise_p_values(
        "1",
        pl.DataFrame({
            "BETA": [0.2], "SE": [0.2 / z], "Z": [z], "P": ["1e-400"],
        }),
        dict(COLUMNS),
        decision="raw",
        policies=policies,
    )

    result, qc, _ = validate_effect_statistics(
        "1", frame, columns, policies=policies,
    )

    assert result.height == 1
    assert qc["z_pval_concordance"]["checked"] == 1
    assert qc["z_pval_concordance"]["discordant"] == 0


def test_beta_se_z_combined_tolerance_allows_rounding_but_rejects_discordance():
    policies = default_policies().with_overrides({
        "validation.beta_se_z_concordance": "reject",
        "validation.z_pval_concordance": "off",
    })

    result, qc, _ = validate_effect_statistics(
        "1", _frame(), dict(COLUMNS), policies=policies,
    )

    assert result["Z"].to_list() == [2.0, 2.02]
    assert qc["beta_se_z_concordance"]["removed"] == 1
    assert qc["removed_by_reason"] == {"beta_se_z_discordant": 1}


def test_beta_se_z_fail_policy_stops_the_chromosome():
    policies = default_policies().with_overrides({
        "validation.beta_se_z_concordance": "fail",
        "validation.z_pval_concordance": "off",
    })

    with pytest.raises(EffectValidationError, match="BETA/SE versus Z"):
        validate_effect_statistics("1", _frame(), dict(COLUMNS), policies=policies)


def test_basic_effect_checks_run_once_in_step_10_without_a_certificate(monkeypatch):
    policies = default_policies().with_overrides({
        "validation.beta_se_z_concordance": "off",
        "validation.z_pval_concordance": "off",
    })
    frame = pl.DataFrame({"BETA": [0.0, 0.2], "SE": [0.1, 0.1]})
    columns = {"beta_col": "BETA", "se_col": "SE", "imp_z_col": None}
    frame, _z_qc, columns = derive_z_score_from_effect_and_standard_error(
        "1", frame, columns, policies=policies,
    )
    applied_reasons = []
    original_apply_rule = effect_validation._apply_rule

    def recording_apply_rule(*args, **kwargs):
        applied_reasons.append(args[2])
        return original_apply_rule(*args, **kwargs)

    monkeypatch.setattr(effect_validation, "_apply_rule", recording_apply_rule)
    result, qc, _ = validate_effect_statistics(
        "1", frame, columns, policies=policies,
    )

    assert result.height == 2
    assert applied_reasons == [
        "se_null", "se_non_finite", "se_non_positive",
        "beta_invalid", "beta_zero", "z_invalid",
    ]
    assert qc["basic_checks_executed"] == applied_reasons
    assert "checks_reused" not in qc
    assert "validation_state" not in inspect.signature(
        derive_z_score_from_effect_and_standard_error
    ).parameters
    assert "upstream_validation" not in inspect.signature(
        validate_effect_statistics
    ).parameters


def test_step_10_applies_the_resolved_basic_policy_after_z_calculation():
    initial_policies = default_policies().with_overrides({
        "validation.beta_se_z_concordance": "off",
        "validation.z_pval_concordance": "off",
    })
    frame = pl.DataFrame({"BETA": [0.0, 0.2], "SE": [0.1, 0.1]})
    columns = {"beta_col": "BETA", "se_col": "SE", "imp_z_col": None}
    frame, _z_qc, columns = derive_z_score_from_effect_and_standard_error(
        "1", frame, columns, policies=initial_policies,
    )
    changed_policies = initial_policies.with_overrides({
        "validation.beta_zero": "reject",
    })

    result, qc, _ = validate_effect_statistics(
        "1", frame, columns, policies=changed_policies,
    )

    assert result["BETA"].to_list() == [0.2]
    assert qc["removed_by_reason"] == {"beta_zero": 1}


def test_basic_keep_policy_counts_six_disjoint_invalidity_classes_once():
    policies = default_policies().with_overrides({
        "validation.se_invalid": "keep",
        "validation.beta_invalid": "keep",
        "validation.beta_zero": "keep",
        "validation.z_invalid": "keep",
        "validation.beta_se_z_concordance": "off",
        "validation.z_pval_concordance": "off",
    })
    frame = pl.DataFrame({
        "BETA": [0.2, 0.2, 0.2, None, 0.0, 0.2],
        "SE": [None, float("inf"), 1.0e-13, 0.1, 0.1, 0.1],
        "Z": [2.0, 2.0, 2.0, 2.0, 0.0, float("inf")],
        "P": [0.05] * 6,
    })

    result, qc, _ = validate_effect_statistics(
        "1", frame, dict(COLUMNS), policies=policies,
    )

    assert result.height == 6
    assert qc["basic_matched_by_reason"] == {
        "se_null": 1,
        "se_non_finite": 1,
        "se_non_positive": 1,
        "beta_invalid": 1,
        "beta_zero": 1,
        "z_invalid": 1,
    }
    assert qc["variants_removed_invalid_beta_se"] == 0


def test_concordance_tolerances_are_schema_validated():
    with pytest.raises(PolicyError, match="beta_se_z_absolute_tolerance"):
        default_policies().with_overrides({
            "validation.beta_se_z_absolute_tolerance": -0.01,
        })


def _frame_with_missing_standard_errors(missing, total=100):
    return pl.DataFrame({
        "BETA": [0.2] * total,
        "SE": [None] * missing + [0.1] * (total - missing),
        "Z": [2.0] * total,
        "P": [0.0455002639] * total,
    })


def test_effect_rejection_guard_defaults_to_ninety_five_percent():
    policies = default_policies()

    assert policies.get("validation.warn_reject_fraction") == 0.2
    assert policies.get("validation.max_reject_fraction") == 0.95


def test_effect_rejection_guard_allows_exactly_ninety_five_percent():
    policies = default_policies().with_overrides({
        "validation.beta_se_z_concordance": "off",
        "validation.z_pval_concordance": "off",
    })

    result, qc, _ = validate_effect_statistics(
        "1", _frame_with_missing_standard_errors(95), dict(COLUMNS),
        policies=policies,
    )

    assert result.height == 5
    assert qc["removed_fraction"] == pytest.approx(0.95)


def test_effect_rejection_guard_fails_above_ninety_five_percent():
    policies = default_policies().with_overrides({
        "validation.beta_se_z_concordance": "off",
        "validation.z_pval_concordance": "off",
    })

    with pytest.raises(
        EffectValidationError,
        match=r"above validation\.max_reject_fraction \(0\.95\)",
    ):
        validate_effect_statistics(
            "1", _frame_with_missing_standard_errors(96), dict(COLUMNS),
            policies=policies,
        )


def test_effect_rejection_guard_can_be_explicitly_disabled_with_one():
    policies = default_policies().with_overrides({
        "validation.max_reject_fraction": 1.0,
        "validation.beta_se_z_concordance": "off",
        "validation.z_pval_concordance": "off",
    })

    result, qc, _ = validate_effect_statistics(
        "1", _frame_with_missing_standard_errors(100), dict(COLUMNS),
        policies=policies,
    )

    assert result.height == 0
    assert qc["removed_fraction"] == pytest.approx(1.0)


def test_warning_rejection_fraction_cannot_exceed_failure_fraction():
    with pytest.raises(PolicyError, match="warn_reject_fraction"):
        default_policies().with_overrides({
            "validation.warn_reject_fraction": 0.96,
            "validation.max_reject_fraction": 0.95,
        })


def test_every_configurable_final_required_field_has_a_rejection_reason():
    assert set(FIELD_COLUMN_KEYS) == set(FIELD_REJECT_REASONS)
    assert set(FIELD_REJECT_REASONS.values()).issubset(REASONS)


@pytest.mark.parametrize("action", ["reject", "keep", "fail"])
def test_absent_required_column_always_fails(action):
    policies = default_policies().with_overrides({
        "final_check.require": ["se"],
        "final_check.on_missing": action,
    })

    with pytest.raises(EffectValidationError, match="absent required column"):
        final_completeness_check(
            "1",
            pl.DataFrame({"BETA": [0.1, 0.2]}),
            {"beta_col": "BETA", "se_col": None},
            policies=policies,
        )


def test_optional_required_field_is_rejected_with_registered_reason(tmp_path):
    frame = pl.DataFrame({"SNP": ["rs1", None]}).with_row_index(
        SOURCE_INPUT_ROW_COLUMN, offset=1,
    )
    collector = RejectCollector(
        source_snapshot=frame,
        logger=None,
        out_path=str(tmp_path / "rejected.tsv"),
        delimiter="\t",
        compress=False,
    )
    policies = default_policies().with_overrides({
        "final_check.require": ["snp"],
        "final_check.on_missing": "reject",
    })

    result, qc, _mapping = final_completeness_check(
        "1",
        frame,
        {"snp_id_col": "SNP"},
        policies=policies,
        rejects=collector,
    )

    assert result["SNP"].to_list() == ["rs1"]
    assert collector.counts() == {"final_missing_snp": 1}
    assert qc["missing_by_field"] == {"snp": 1}
