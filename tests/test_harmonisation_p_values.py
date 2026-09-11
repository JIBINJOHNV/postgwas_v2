"""Scientific contracts for study-wide p-value representation detection."""

from decimal import Decimal
import math
from pathlib import Path

import polars as pl
import pytest
import yaml

from postgwas.modules.harmonisation.p_values import (
    AmbiguousPValueTypeError,
    CLIPPED_LOW_COLUMN,
    EXACT_RAW_P_COLUMN,
    ExcessiveNegativePValueError,
    FLOAT_UNDERFLOW_COLUMN,
    LOG_P_COLUMN,
    PValueTypeError,
    REPORTED_ZERO_COLUMN,
    SOURCE_TEXT_COLUMN,
    convert_negative_log10_to_p_value,
    detect_p_value_type,
    harmonise_p_values,
)
from postgwas.modules.harmonisation.policies import (
    PolicyError,
    default_policies,
    registry_path,
)
from postgwas.modules.harmonisation.rejects import (
    RejectCollector,
    SOURCE_INPUT_ROW_COLUMN,
)
from postgwas.modules.harmonisation.study_properties import (
    StudyPropertyError,
    resolve_study_properties,
)
from postgwas.modules.harmonisation.summary_statistics_io import (
    load_summary_statistics_table,
)


def _uniform_quantiles(size=1000):
    return [(index + 0.5) / size for index in range(size)]


def test_only_auto_raw_and_neglog10_are_schema_supported():
    registry = yaml.safe_load(Path(registry_path()).read_text(encoding="utf-8"))
    pvalue = registry["policies"]["pvalue"]

    assert pvalue["type"]["validator"]["members"] == [
        "auto",
        "raw",
        "neglog10",
    ]
    assert pvalue["mlogp_expected_median"]["default"] == pytest.approx(
        0.3010299956639812
    )
    assert default_policies().get("pvalue.mlogp_expected_median") == pytest.approx(
        pvalue["mlogp_expected_median"]["default"]
    )
    assert pvalue["output_column"]["default"] == "PVAL"
    assert default_policies().get("pvalue.output_column") == "PVAL"
    assert pvalue["max_negative_fraction"]["default"] == pytest.approx(0.001)
    assert default_policies().get("pvalue.max_negative_fraction") == pytest.approx(
        0.001
    )
    assert "mlogp_max" not in pvalue
    with pytest.raises(PolicyError, match="pvalue.type"):
        default_policies().with_overrides({"pvalue.type": "negln"})


def test_raw_pvalue_source_text_preserves_float64_underflow_and_literal_zero():
    source = ["0.05", "1e-300", "1e-320", "1e-400", "0"]
    frame, qc, _mapping = harmonise_p_values(
        chromosome="1",
        df=pl.DataFrame({"P": source}),
        sample_column_dict={"pval_col": "P"},
        decision="raw",
        policies=default_policies(),
    )

    assert frame[SOURCE_TEXT_COLUMN].to_list() == source
    assert frame[FLOAT_UNDERFLOW_COLUMN].to_list() == [
        False, False, False, True, False,
    ]
    assert frame[REPORTED_ZERO_COLUMN].to_list() == [
        False, False, False, False, True,
    ]
    assert frame[CLIPPED_LOW_COLUMN].to_list() == [
        False, False, False, False, True,
    ]
    assert frame[EXACT_RAW_P_COLUMN].to_list() == [None, None, None, "1e-400", None]
    raw = frame["PVAL"].to_list()
    assert raw[0] == pytest.approx(0.05)
    assert raw[1] == pytest.approx(1e-300, rel=1e-12, abs=0.0)
    assert raw[2] > 0.0
    assert raw[3] is None
    assert raw[4] == pytest.approx(1e-300, rel=1e-12, abs=0.0)
    expected_log = [math.log(0.05), math.log(1e-300), math.log(1e-320),
                    -400.0 * math.log(10.0), math.log(1e-300)]
    assert frame[LOG_P_COLUMN].to_list() == pytest.approx(expected_log)
    assert qc["pvalue_float_underflow"] == 1


def test_neglog10_extremes_keep_distinct_log_magnitudes_without_clipping():
    frame, qc, _mapping = harmonise_p_values(
        chromosome="1",
        df=pl.DataFrame({"P": ["1.30103", "300", "320", "400"]}),
        sample_column_dict={"pval_col": "P"},
        decision="neglog10",
        policies=default_policies(),
    )

    assert frame[LOG_P_COLUMN].to_list() == pytest.approx([
        -1.30103 * math.log(10.0),
        -300.0 * math.log(10.0),
        -320.0 * math.log(10.0),
        -400.0 * math.log(10.0),
    ])
    assert frame[FLOAT_UNDERFLOW_COLUMN].to_list() == [False, False, False, True]
    assert frame[CLIPPED_LOW_COLUMN].to_list() == [False, False, False, False]
    assert frame[EXACT_RAW_P_COLUMN][3] is not None
    assert Decimal(frame[EXACT_RAW_P_COLUMN][3]) == Decimal("1e-400")
    assert frame["PVAL"][1] == pytest.approx(1e-300, rel=1e-12, abs=0.0)
    assert frame["PVAL"][2] > 0.0
    assert frame["PVAL"][3] is None
    assert "mlogp_values_clipped_at_max" not in qc
    assert qc["pvalue_float_underflow"] == 1


def test_negative_raw_underflow_is_not_misclassified_as_zero_or_positive():
    frame, qc, _mapping = harmonise_p_values(
        chromosome="1",
        df=pl.DataFrame({"ID": ["negative", "positive"],
                         "P": ["-1e-400", "1e-400"]}),
        sample_column_dict={"pval_col": "P"},
        decision="raw",
        policies=default_policies(),
    )

    assert frame["ID"].to_list() == ["positive"]
    assert frame[REPORTED_ZERO_COLUMN].to_list() == [False]
    assert qc["variants_with_lt0_pvalues_removed"] == 1


def test_initial_table_read_keeps_extreme_pvalue_token_as_text(tmp_path):
    path = tmp_path / "study.tsv"
    path.write_text("SNP\tP\na\t1e-400\nb\t0\n", encoding="utf-8")

    frame, _mapping = load_summary_statistics_table(
        str(path),
        {"snp_id_col": "SNP", "pval_col": "P"},
        policies=default_policies(),
    )

    assert frame.schema["P"] == pl.String
    assert frame["P"].to_list() == ["1e-400", "0"]


@pytest.mark.parametrize(
    ("decision", "supplied"),
    [("raw", 0.05), ("neglog10", -math.log10(0.05))],
)
def test_every_supported_input_scale_leaves_the_same_canonical_raw_pvalue(
    decision, supplied,
):
    policies = default_policies()
    frame, qc, mapping = harmonise_p_values(
        chromosome="1",
        df=pl.DataFrame({"P": [supplied]}),
        sample_column_dict={"pval_col": "P"},
        decision=decision,
        policies=policies,
    )

    output_column = policies.get("pvalue.output_column")
    assert mapping["pval_col"] == output_column
    assert frame[output_column].to_list() == pytest.approx([0.05])
    assert "LP" not in frame.columns
    assert qc["pvalue_input_column"] == "P"
    assert qc["pvalue_output_column"] == output_column


def test_canonical_raw_pvalue_column_comes_from_validated_yaml_policy():
    policies = default_policies().with_overrides({
        "pvalue.output_column": "RAW_ASSOCIATION_P",
    })
    frame, _qc, mapping = harmonise_p_values(
        chromosome="1",
        df=pl.DataFrame({"P": [0.05]}),
        sample_column_dict={"pval_col": "P"},
        decision="raw",
        policies=policies,
    )

    assert mapping["pval_col"] == "RAW_ASSOCIATION_P"
    assert frame["RAW_ASSOCIATION_P"].to_list() == [0.05]
    assert "PVAL" not in frame.columns
    with pytest.raises(PolicyError, match="pvalue.output_column"):
        policies.with_overrides({"pvalue.output_column": "invalid column"})
    with pytest.raises(PolicyError, match="pvalue.output_column"):
        policies.with_overrides({"pvalue.output_column": "LP"})


def test_canonical_raw_pvalue_column_never_overwrites_unrelated_input_data():
    with pytest.raises(PValueTypeError, match="Refusing to overwrite unrelated data"):
        harmonise_p_values(
            chromosome="1",
            df=pl.DataFrame({"P": [0.05], "PVAL": [0.75]}),
            sample_column_dict={"pval_col": "P"},
            decision="raw",
            policies=default_policies(),
        )


def test_detection_policy_relationships_are_schema_validated():
    policies = default_policies()

    with pytest.raises(PolicyError, match="pvalue.mlogp_detect_min_count"):
        policies.with_overrides({"pvalue.mlogp_detect_min_count": 1})
    with pytest.raises(PolicyError, match="pvalue.mlogp_detect_proportion"):
        policies.with_overrides({"pvalue.mlogp_detect_proportion": 0.0})
    with pytest.raises(PolicyError, match="pvalue.mlogp_median_tolerance"):
        policies.with_overrides({
            "pvalue.mlogp_expected_median": 0.2,
            "pvalue.mlogp_median_tolerance": 0.2,
        })
    with pytest.raises(PolicyError, match="pvalue.mlogp_detect_threshold"):
        policies.with_overrides({"pvalue.mlogp_detect_threshold": 1.01})
    with pytest.raises(PolicyError, match="pvalue.max_negative_fraction"):
        policies.with_overrides({"pvalue.max_negative_fraction": 1.01})


def test_exactly_point_one_percent_negative_warns_but_does_not_stop(capsys):
    values = [-0.2] + [0.5] * 999

    decisions = resolve_study_properties(
        pl.DataFrame({"P": values}),
        {"pval_col": "P"},
        policies=default_policies(),
    )

    assert decisions["pvalue_type"] == "raw"
    assert decisions["pvalue_evidence"]["negative_fraction"] == pytest.approx(0.001)
    assert "NEGATIVE P-VALUES" in decisions["pvalue_negative_warning"]
    assert "will be rejected" in capsys.readouterr().out


@pytest.mark.parametrize("configured_type", ["auto", "raw", "neglog10"])
def test_more_than_point_one_percent_negative_stops_before_chromosomes(
    configured_type,
):
    values = [-0.2, -0.1] + [0.5] * 998
    policies = default_policies().with_overrides({"pvalue.type": configured_type})

    with pytest.raises(
        StudyPropertyError,
        match=(
            r"2 negative usable values out of 1,000 \(0\.2000%\).*"
            r"pvalue\.max_negative_fraction=0\.001.*No chromosome processing"
        ),
    ):
        resolve_study_properties(
            pl.DataFrame({"P": values}),
            {"pval_col": "P"},
            policies=policies,
        )


def test_sample_sheet_declaration_cannot_bypass_negative_fraction_gate():
    values = [-0.2, -0.1] + [0.5] * 998

    with pytest.raises(StudyPropertyError, match="incorrect p-value-column mapping"):
        resolve_study_properties(
            pl.DataFrame({"P": values}),
            {"pval_col": "P", "declared_pvalue_type": "raw"},
            policies=default_policies(),
        )


def test_missing_values_cannot_dilute_the_negative_fraction():
    values = [-0.2, 0.5] + [None] * 98

    with pytest.raises(
        ExcessiveNegativePValueError,
        match=r"1 negative usable value out of 2 \(50\.0000%\)",
    ):
        detect_p_value_type(
            pl.DataFrame({"P": values}), "P", default_policies()
        )


@pytest.mark.parametrize("decision", ["raw", "neglog10"])
def test_allowed_negative_cells_are_rejected_by_shared_canonical_gate(decision):
    valid_value = 0.5 if decision == "raw" else -math.log10(0.5)
    frame, qc, _mapping = harmonise_p_values(
        chromosome="1",
        df=pl.DataFrame({"ID": ["invalid", "valid"], "P": [-0.2, valid_value]}),
        sample_column_dict={"pval_col": "P"},
        decision=decision,
        policies=default_policies(),
    )

    assert frame["ID"].to_list() == ["valid"]
    assert frame["PVAL"].to_list() == pytest.approx([0.5])
    assert qc["pvalue_source_negative"] == 1
    assert qc["variants_with_lt0_pvalues_removed"] == 1


@pytest.mark.parametrize("decision", ["raw", "neglog10"])
def test_missing_and_nonfinite_sources_have_step_seven_reject_reasons(
    decision, tmp_path,
):
    valid = "0.5" if decision == "raw" else str(-math.log10(0.5))
    source = pl.DataFrame({
        "ID": ["valid", "missing", "unparseable", "posinf", "neginf", "nan"],
        "P": [valid, None, "not-a-number", "inf", "-inf", "NaN"],
    }).with_row_index(SOURCE_INPUT_ROW_COLUMN, offset=1)
    collector = RejectCollector(
        source_snapshot=source,
        logger=None,
        out_path=str(tmp_path / (decision + "_rejected.tsv")),
        delimiter="\t",
        compress=False,
    )

    frame, qc, _mapping = harmonise_p_values(
        chromosome="1",
        df=source,
        sample_column_dict={"pval_col": "P"},
        decision=decision,
        policies=default_policies(),
        rejects=collector,
    )

    assert frame["ID"].to_list() == ["valid"]
    assert collector.counts() == {"pval_null": 2, "pval_out_of_range": 3}
    assert qc["pvalue_source_missing_or_unparseable"] == 2
    assert qc["pvalue_source_non_finite"] == 3
    if decision == "neglog10":
        assert qc["mlogp_missing_or_unparseable"] == 2
        assert qc["mlogp_non_finite"] == 3
        assert "mlogp_invalid_or_missing" not in qc


@pytest.mark.parametrize("decision", ["raw", "neglog10"])
@pytest.mark.parametrize(
    ("action", "expected_ids", "expected_invalid_p"),
    [
        ("clip", ["valid"], None),
        ("reject", ["valid"], None),
        ("null", ["invalid", "valid"], None),
    ],
)
def test_negative_source_values_follow_out_of_range_policy(
    decision, action, expected_ids, expected_invalid_p,
):
    valid_value = 0.5 if decision == "raw" else -math.log10(0.5)
    policies = default_policies().with_overrides({
        "pvalue.out_of_range": action,
    })

    frame, qc, _mapping = harmonise_p_values(
        chromosome="1",
        df=pl.DataFrame({"ID": ["invalid", "valid"], "P": [-0.2, valid_value]}),
        sample_column_dict={"pval_col": "P"},
        decision=decision,
        policies=policies,
    )

    assert frame["ID"].to_list() == expected_ids
    if action == "null":
        assert frame["PVAL"].to_list()[0] is expected_invalid_p
        assert frame["PVAL"].to_list()[1] == pytest.approx(0.5)
    else:
        assert frame["PVAL"].to_list() == pytest.approx([0.5])
    assert qc["pvalue_source_negative"] == 1


@pytest.mark.parametrize("decision", ["raw", "neglog10"])
def test_out_of_range_fail_stops_for_negative_source_value(decision):
    valid_value = 0.5 if decision == "raw" else -math.log10(0.5)
    policies = default_policies().with_overrides({
        "pvalue.out_of_range": "fail",
    })

    with pytest.raises(PValueTypeError, match="out_of_range.*fail"):
        harmonise_p_values(
            chromosome="1",
            df=pl.DataFrame({"P": [-0.2, valid_value]}),
            sample_column_dict={"pval_col": "P"},
            decision=decision,
            policies=policies,
        )


def test_negative_neglog10_below_float64_is_never_converted_to_p_one():
    frame, qc, _mapping = harmonise_p_values(
        chromosome="1",
        df=pl.DataFrame({
            "ID": ["negative_underflow", "valid"],
            "P": ["-1e-400", str(-math.log10(0.5))],
        }),
        sample_column_dict={"pval_col": "P"},
        decision="neglog10",
        policies=default_policies(),
    )

    assert frame["ID"].to_list() == ["valid"]
    assert frame["PVAL"].to_list() == pytest.approx([0.5])
    assert qc["pvalue_source_negative"] == 1


def test_low_level_neglog10_converter_never_floors_negative_value_to_p_one():
    with pytest.raises(PValueTypeError, match="Cannot convert 1 negative value"):
        convert_negative_log10_to_p_value(
            pl.DataFrame({"P": [-0.2]}),
            {"pval_col": "P"},
            output_col="PVAL",
            policies=default_policies(),
        )


def test_one_extreme_sentinel_cannot_convert_a_raw_p_column():
    values = _uniform_quantiles(10_000) + [999.0]

    detected, evidence = detect_p_value_type(
        pl.DataFrame({"P": values}), "P", default_policies()
    )

    assert detected == "raw"
    assert evidence["n_above_threshold"] == 1
    assert evidence["count_supports_neglog10"] is False
    assert evidence["on_log_scale"] is False
    assert "within_range_fraction" not in evidence


def test_small_file_with_one_outlier_is_refused_instead_of_converted():
    values = _uniform_quantiles(100) + [999.0]

    with pytest.raises(AmbiguousPValueTypeError, match="single sentinel"):
        detect_p_value_type(
            pl.DataFrame({"P": values}), "P", default_policies()
        )


def test_standard_neglog10_distribution_requires_count_fraction_and_yaml_median():
    values = [-math.log10(value) for value in _uniform_quantiles()]

    detected, evidence = detect_p_value_type(
        pl.DataFrame({"P": values}), "P", default_policies()
    )

    assert detected == "neglog10"
    assert evidence["count_supports_neglog10"] is True
    assert evidence["fraction_supports_neglog10"] is True
    assert evidence["median"] == pytest.approx(
        default_policies().get("pvalue.mlogp_expected_median"), rel=1e-5
    )


def test_detector_uses_configured_median_instead_of_a_code_constant():
    values = [-math.log10(value) for value in _uniform_quantiles()]
    policies = default_policies().with_overrides({
        "pvalue.mlogp_expected_median": 0.4,
        "pvalue.mlogp_median_tolerance": 0.01,
    })

    with pytest.raises(AmbiguousPValueTypeError, match="configured -log10 median"):
        detect_p_value_type(pl.DataFrame({"P": values}), "P", policies)


@pytest.mark.parametrize(
    "transform",
    [math.log10, math.log],
    ids=["signed_log10", "signed_ln"],
)
def test_signed_log_columns_fail_before_row_filtering(transform):
    values = [transform(value) for value in _uniform_quantiles()]

    with pytest.raises(PValueTypeError, match="negative usable values"):
        detect_p_value_type(
            pl.DataFrame({"P": values}), "P", default_policies()
        )


def test_positive_negative_natural_log_is_not_misread_as_neglog10():
    values = [-math.log(value) for value in _uniform_quantiles()]

    with pytest.raises(AmbiguousPValueTypeError, match="does not agree"):
        detect_p_value_type(
            pl.DataFrame({"P": values}), "P", default_policies()
        )


def test_unsupported_signed_scale_fails_during_study_resolution_before_fanout():
    values = [math.log10(value) for value in _uniform_quantiles()]

    with pytest.raises(StudyPropertyError, match="No chromosome processing was started"):
        resolve_study_properties(
            pl.DataFrame({"P": values}),
            {"pval_col": "P"},
            policies=default_policies(),
        )
