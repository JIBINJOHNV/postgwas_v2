"""Regression tests for p=0 when SE cannot be reconstructed exactly."""

import polars as pl
import pytest

from postgwas.modules.harmonisation.effect_from_z import (
    SE_UNAVAILABLE_FROM_Z_COLUMN,
    derive_effect_and_standard_error_from_z,
)
from postgwas.modules.harmonisation.p_values import (
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


def _raw_pvalues(frame, mapping, policies=None):
    return harmonise_p_values(
        "1",
        frame,
        mapping,
        decision="raw",
        policies=policies or default_policies(),
    )


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
    assert frame["LP"].to_list() == [None]
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
