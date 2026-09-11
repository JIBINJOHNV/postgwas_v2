"""Regression tests for p=0 when SE cannot be reconstructed exactly."""

import math

import polars as pl
import pytest
from scipy.special import ndtri_exp
from scipy.stats import norm

from postgwas.modules.harmonisation.effect_from_z import (
    SE_UNAVAILABLE_FROM_Z_COLUMN,
    derive_effect_and_standard_error_from_z,
)
from postgwas.modules.harmonisation.effect_validation import (
    validate_effect_statistics,
)
from postgwas.modules.harmonisation.p_values import (
    CLIPPED_LOW_COLUMN,
    FLOAT_UNDERFLOW_COLUMN,
    REPORTED_ZERO_COLUMN,
    harmonise_p_values,
)
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.rejects import (
    RejectCollector,
    SOURCE_INPUT_ROW_COLUMN,
)
from postgwas.modules.harmonisation.standard_error import (
    SE_FROM_CLIPPED_PVAL_COL,
    SE_FROM_ZERO_P_APPROXIMATION_COL,
    derive_standard_error_from_effect_and_p_value,
)
from postgwas.modules.harmonisation.z_score import (
    derive_z_score_from_effect_and_standard_error,
)


def _raw_pvalues(frame, mapping, policies=None):
    return harmonise_p_values(
        "1",
        frame,
        mapping,
        decision="raw",
        policies=policies or default_policies(),
    )


def test_partial_supplied_se_fills_only_null_cells_from_beta_and_pvalue():
    policies = default_policies()
    assert policies.get("pvalue.derive_partial_missing_se") is True
    frame, _qc, mapping = _raw_pvalues(
        pl.DataFrame({
            "BETA": [0.1, 0.2, 0.3],
            "SE": [0.02, None, 0.05],
            "P": [0.05, 0.01, 0.001],
        }),
        {"beta_col": "BETA", "se_col": "SE", "pval_col": "P"},
        policies,
    )

    result, qc, result_mapping = derive_standard_error_from_effect_and_p_value(
        "1", frame, mapping, policies=policies
    )

    expected = abs(0.2 / norm.isf(0.01 / policies.get("pvalue.se_tail")))
    assert result["SE"].to_list() == pytest.approx([0.02, expected, 0.05])
    assert result_mapping["se_col"] == "SE"
    assert qc["se_missing_on_entry"] == 1
    assert qc["se_missing_recovered_from_pvalue"] == 1
    assert qc["se_from_pvalue"] == 1

    with_z, _z_qc, z_mapping = derive_z_score_from_effect_and_standard_error(
        "1", result, result_mapping, policies=policies
    )
    validated, effect_qc, _ = validate_effect_statistics(
        "1", with_z, z_mapping, policies=policies
    )
    assert validated.height == 3
    assert effect_qc["variants_removed_due_to_missing_se"] == 0


@pytest.mark.parametrize(
    ("decision", "pvalues"),
    [
        ("raw", ["1e-300", "1e-320", "1e-400"]),
        ("neglog10", ["300", "320", "400"]),
    ],
)
def test_extreme_pvalues_produce_distinct_log_domain_standard_errors(
    decision, pvalues,
):
    policies = default_policies()
    frame, _qc, mapping = harmonise_p_values(
        "1",
        pl.DataFrame({"BETA": ["0.2", "0.2", "0.2"], "P": pvalues}),
        {"beta_col": "BETA", "pval_col": "P", "se_col": None},
        decision=decision,
        policies=policies,
    )

    result, qc, _ = derive_standard_error_from_effect_and_p_value(
        "1", frame, mapping, policies=policies
    )

    log_p = -pl.Series([300.0, 320.0, 400.0]) * math.log(10.0)
    z = -ndtri_exp(log_p.to_numpy() - math.log(2.0))
    expected = abs(0.2 / z)
    assert result["SE"].to_list() == pytest.approx(expected, rel=1e-12)
    assert len(set(result["SE"].to_list())) == 3
    assert frame[CLIPPED_LOW_COLUMN].to_list() == [False, False, False]
    assert frame[FLOAT_UNDERFLOW_COLUMN].to_list() == [False, False, True]
    assert qc["se_from_clipped_pval"] == 0
    assert qc["se_missing_recovered_from_pvalue"] == 3


def test_partial_supplied_se_recovery_can_be_disabled_in_yaml_policy():
    policies = default_policies().with_overrides({
        "pvalue.derive_partial_missing_se": False,
    })
    frame, _qc, mapping = _raw_pvalues(
        pl.DataFrame({"BETA": [0.1, 0.2], "SE": [0.02, None], "P": [0.05, 0.01]}),
        {"beta_col": "BETA", "se_col": "SE", "pval_col": "P"},
        policies,
    )

    result, qc, _ = derive_standard_error_from_effect_and_p_value(
        "1", frame, mapping, policies=policies
    )

    assert result["SE"].to_list() == [0.02, None]
    assert qc["status"] == "partial missing SE recovery disabled"
    assert qc["se_missing_recovered_from_pvalue"] == 0


def test_partial_supplied_se_does_not_bypass_zero_pvalue_safety_policy():
    frame, _qc, mapping = _raw_pvalues(
        pl.DataFrame({"BETA": [0.1, 0.2], "SE": [0.02, None], "P": [0.05, 0.0]}),
        {"beta_col": "BETA", "se_col": "SE", "pval_col": "P"},
    )

    with pytest.raises(ValueError, match="1 highly significant variant"):
        derive_standard_error_from_effect_and_p_value(
            "1", frame, mapping, policies=default_policies()
        )


def test_failed_beta_z_recovery_falls_back_to_beta_and_pvalue():
    policies = default_policies()
    frame, _qc, mapping = derive_effect_and_standard_error_from_z(
        "1",
        pl.DataFrame({"BETA": [0.1, 0.2], "Z": [2.0, 0.0], "P": [0.01, 0.05]}),
        {
            "beta_or_col": "BETA",
            "beta_col": "BETA",
            "se_col": None,
            "imp_z_col": "Z",
            "pval_col": "P",
        },
        policies=policies,
    )
    frame, _qc, mapping = _raw_pvalues(frame, mapping, policies)

    result, qc, _ = derive_standard_error_from_effect_and_p_value(
        "1", frame, mapping, policies=policies
    )

    expected = abs(0.2 / norm.isf(0.05 / policies.get("pvalue.se_tail")))
    assert result["SE"].to_list() == pytest.approx([0.05, expected])
    assert result[SE_UNAVAILABLE_FROM_Z_COLUMN].to_list() == [False, True]
    assert qc["se_missing_recovered_from_pvalue"] == 1


def test_default_fails_clearly_without_discarding_zero_p_variant():
    frame, _qc, mapping = _raw_pvalues(
        pl.DataFrame({"BETA": [0.1, 0.2], "P": [0.0, 0.05]}),
        {"beta_col": "BETA", "pval_col": "P", "se_col": None},
    )

    with pytest.raises(ValueError) as error:
        derive_standard_error_from_effect_and_p_value(
            "1", frame, mapping, policies=default_policies()
        )

    message = str(error.value)
    assert "Chromosome 1: 1 highly significant variant(s)" in message
    assert "No variants were silently discarded" in message
    assert "--zero-p-se-action approximate" in message
    assert frame.height == 2
    assert frame[REPORTED_ZERO_COLUMN].to_list() == [True, False]


def test_supplied_se_or_usable_z_preserves_zero_p_variant():
    supplied, _qc, supplied_mapping = _raw_pvalues(
        pl.DataFrame({"BETA": [0.1], "SE": [0.02], "P": [0.0]}),
        {"beta_col": "BETA", "se_col": "SE", "pval_col": "P"},
    )
    supplied, supplied_qc, _ = derive_standard_error_from_effect_and_p_value(
        "1", supplied, supplied_mapping, policies=default_policies()
    )
    assert supplied["SE"].to_list() == [0.02]
    assert supplied_qc["zero_p_missing_se"] == 0

    z_frame, _z_qc, z_mapping = derive_effect_and_standard_error_from_z(
        "1",
        pl.DataFrame({"BETA": [0.1], "Z": [2.0], "P": [0.0]}),
        {
            "beta_or_col": "BETA",
            "beta_col": "BETA",
            "se_col": None,
            "imp_z_col": "Z",
            "pval_col": "P",
        },
        policies=default_policies(),
    )
    z_frame, _p_qc, z_mapping = _raw_pvalues(z_frame, z_mapping)
    z_frame, z_qc, _ = derive_standard_error_from_effect_and_p_value(
        "1", z_frame, z_mapping, policies=default_policies()
    )
    assert z_frame["SE"].to_list() == [0.05]
    assert z_frame[SE_UNAVAILABLE_FROM_Z_COLUMN].to_list() == [False]
    assert z_qc["zero_p_missing_se"] == 0


def test_approximate_override_marks_only_unresolved_zero_p_rows():
    policies = default_policies().with_overrides(
        {"pvalue.zero_missing_se": "approximate"}
    )
    frame, _z_qc, mapping = derive_effect_and_standard_error_from_z(
        "1",
        pl.DataFrame({"BETA": [0.1, 0.2], "Z": [2.0, 0.0], "P": [0.0, 0.0]}),
        {
            "beta_or_col": "BETA",
            "beta_col": "BETA",
            "se_col": None,
            "imp_z_col": "Z",
            "pval_col": "P",
        },
        policies=policies,
    )
    frame, _p_qc, mapping = _raw_pvalues(frame, mapping, policies)
    result, qc, _ = derive_standard_error_from_effect_and_p_value(
        "1", frame, mapping, policies=policies
    )

    assert result["SE"][0] == pytest.approx(0.05)
    assert result["SE"][1] is not None
    assert result[SE_FROM_ZERO_P_APPROXIMATION_COL].to_list() == [False, True]
    assert result[SE_FROM_CLIPPED_PVAL_COL].to_list() == [False, False]
    assert qc["se_from_zero_p_approximation"] == 1
    assert qc["zero_p_missing_se_action"] == "approximate"


def test_approximate_override_uses_provenance_after_range_policy_blanks_zero():
    policies = default_policies().with_overrides(
        {
            "pvalue.zero_missing_se": "approximate",
            "pvalue.out_of_range": "null",
        }
    )
    frame, _qc, mapping = _raw_pvalues(
        pl.DataFrame({"BETA": [0.1], "P": [0.0]}),
        {"beta_col": "BETA", "pval_col": "P", "se_col": None},
        policies,
    )
    assert frame["PVAL"].to_list() == [None]
    assert frame[REPORTED_ZERO_COLUMN].to_list() == [True]

    result, qc, _ = derive_standard_error_from_effect_and_p_value(
        "1", frame, mapping, policies=policies
    )

    assert result["SE"][0] is not None
    assert result[SE_FROM_ZERO_P_APPROXIMATION_COL].to_list() == [True]
    assert qc["se_from_zero_p_approximation"] == 1


def test_reject_override_writes_specific_provenance(tmp_path):
    policies = default_policies().with_overrides(
        {"pvalue.zero_missing_se": "reject"}
    )
    frame, _qc, mapping = _raw_pvalues(
        pl.DataFrame({"ID": ["zero", "ordinary"], "BETA": [0.1, 0.2], "P": [0.0, 0.05]}),
        {"beta_col": "BETA", "pval_col": "P", "se_col": None},
        policies,
    )
    frame = frame.with_row_index(SOURCE_INPUT_ROW_COLUMN, offset=1)
    collector = RejectCollector(
        source_snapshot=frame,
        logger=None,
        out_path=str(tmp_path / "rejected.tsv"),
        delimiter="\t",
        compress=False,
    )

    result, qc, _ = derive_standard_error_from_effect_and_p_value(
        "1", frame, mapping, policies=policies, rejects=collector
    )

    assert result["ID"].to_list() == ["ordinary"]
    assert collector.counts() == {"zero_p_missing_se": 1}
    assert collector.frame()["reject_reason"].to_list() == ["zero_p_missing_se"]
    assert qc["zero_p_missing_se_rejected"] == 1
