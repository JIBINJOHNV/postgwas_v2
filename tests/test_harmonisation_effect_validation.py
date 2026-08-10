"""Scientific regression tests for chromosome effect-statistic concordance."""

import polars as pl
import pytest

from postgwas.modules.harmonisation import effect_validation
from postgwas.modules.harmonisation.effect_validation import (
    EffectValidationError,
    FIELD_COLUMN_KEYS,
    FIELD_REJECT_REASONS,
    final_completeness_check,
    validate_effect_statistics,
)
from postgwas.modules.harmonisation.policies import PolicyError, default_policies
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


def test_basic_effect_checks_are_reused_only_for_the_exact_z_step_frame(monkeypatch):
    policies = default_policies().with_overrides({
        "validation.beta_se_z_concordance": "off",
        "validation.z_pval_concordance": "off",
    })
    frame = pl.DataFrame({"BETA": [0.0, 0.2], "SE": [0.1, 0.1]})
    columns = {"beta_col": "BETA", "se_col": "SE", "imp_z_col": None}
    state = {}

    frame, _z_qc, columns = derive_z_score_from_effect_and_standard_error(
        "1", frame, columns, policies=policies, validation_state=state,
    )
    numeric_calls = []
    original_numeric = effect_validation._numeric

    def counting_numeric(dataframe, column):
        numeric_calls.append(column)
        return original_numeric(dataframe, column)

    monkeypatch.setattr(effect_validation, "_numeric", counting_numeric)
    result, qc, _ = validate_effect_statistics(
        "1", frame, columns, policies=policies, upstream_validation=state,
    )

    assert result.height == 2
    assert qc["checks_reused"] == [
        "beta_invalid", "beta_zero", "se_invalid", "z_invalid",
    ]
    assert numeric_calls == []

    _result, cloned_qc, _ = validate_effect_statistics(
        "1", frame.clone(), columns, policies=policies, upstream_validation=state,
    )
    assert cloned_qc["checks_reused"] == []
    assert numeric_calls == ["SE", "BETA", "imp_z_col"]


def test_changed_policy_reruns_only_its_basic_effect_check():
    initial_policies = default_policies().with_overrides({
        "validation.beta_se_z_concordance": "off",
        "validation.z_pval_concordance": "off",
    })
    frame = pl.DataFrame({"BETA": [0.0, 0.2], "SE": [0.1, 0.1]})
    columns = {"beta_col": "BETA", "se_col": "SE", "imp_z_col": None}
    state = {}
    frame, _z_qc, columns = derive_z_score_from_effect_and_standard_error(
        "1", frame, columns, policies=initial_policies, validation_state=state,
    )
    changed_policies = initial_policies.with_overrides({
        "validation.beta_zero": "reject",
    })

    result, qc, _ = validate_effect_statistics(
        "1", frame, columns, policies=changed_policies,
        upstream_validation=state,
    )

    assert result["BETA"].to_list() == [0.2]
    assert qc["checks_reused"] == ["beta_invalid", "se_invalid", "z_invalid"]
    assert qc["removed_by_reason"] == {"beta_zero": 1}


def test_concordance_tolerances_are_schema_validated():
    with pytest.raises(PolicyError, match="beta_se_z_absolute_tolerance"):
        default_policies().with_overrides({
            "validation.beta_se_z_absolute_tolerance": -0.01,
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
