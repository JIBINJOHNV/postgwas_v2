"""Scientific and configuration contracts for missing per-variant sample size."""

import polars as pl
import pytest

from postgwas.modules.harmonisation.policies import PolicyError, default_policies
from postgwas.modules.harmonisation.rejects import (
    RejectCollector,
    SOURCE_INPUT_ROW_COLUMN,
)
from postgwas.modules.harmonisation.sample_size import (
    harmonise_sample_sizes,
    prepare_missing_sample_sizes,
)


def _policies(action="remove", maximum=0.01):
    return default_policies().with_overrides({
        "sample_size.missing_action": action,
        "sample_size.max_missing_fraction": maximum,
    })


def _mapping(**updates):
    values = {
        "ncase_col": None,
        "ncontrol_col": "N",
        "ncase": None,
        "ncontrol": None,
    }
    values.update(updates)
    return values


def test_missing_sample_size_policy_defaults_and_validation():
    policies = default_policies()
    assert policies.get("sample_size.missing_action") == "remove"
    assert policies.get("sample_size.max_missing_fraction") == pytest.approx(0.01)
    assert policies.get("sample_size.min_value") == 1
    assert policies.get("sample_size.cases_only") == "fail"

    with pytest.raises(PolicyError, match="sample_size.missing_action"):
        _policies(action="guess")
    with pytest.raises(PolicyError, match="sample_size.max_missing_fraction"):
        _policies(maximum=1.01)
    with pytest.raises(PolicyError, match="sample_size.min_value"):
        policies.with_overrides({"sample_size.min_value": 0})
    with pytest.raises(PolicyError, match="sample_size.cases_only"):
        policies.with_overrides({"sample_size.cases_only": "use_ncase"})


def test_case_only_policy_stops_with_clear_scientific_reason():
    with pytest.raises(ValueError, match="cannot determine case/control effective sample size"):
        harmonise_sample_sizes(
            "1",
            pl.DataFrame({"NCASE": [10_000]}),
            _mapping(ncase_col="NCASE", ncontrol_col=None),
            policies=default_policies(),
        )


@pytest.mark.parametrize(
    ("frame", "mapping"),
    [
        (
            pl.DataFrame({"NCASE": [10_000], "NCONTROL": [20_000]}),
            _mapping(ncase_col="NCASE", ncontrol_col="NCONTROL"),
        ),
        (
            pl.DataFrame({"NCONTROL": [20_000]}),
            _mapping(ncase=10_000, ncontrol_col="NCONTROL"),
        ),
    ],
)
def test_quantitative_trait_rejects_case_count_inputs(frame, mapping):
    policies = default_policies().with_overrides({
        "sample_size.trait_type": "quantitative",
    })

    with pytest.raises(ValueError, match="quantitative.*case-count inputs are not allowed"):
        harmonise_sample_sizes("1", frame, mapping, policies=policies)


def test_quantitative_trait_uses_control_field_as_total_sample_size():
    policies = default_policies().with_overrides({
        "sample_size.trait_type": "quantitative",
    })

    result, qc, mapping = harmonise_sample_sizes(
        "1",
        pl.DataFrame({"N": [20_000]}),
        _mapping(),
        policies=policies,
    )

    assert result["Neff"].to_list() == [20_000]
    assert mapping["neff_col"] == "Neff"
    assert qc["Neff_status"] == "fallback_ncontrol_only:N"


def test_exact_maximum_is_allowed_but_greater_fraction_stops_clearly():
    frame = pl.DataFrame({"N": [100, 100, 100, None]})

    plan = prepare_missing_sample_sizes(
        frame, _mapping(), policies=_policies(maximum=0.25), scope="dataset 'study'",
    )
    assert plan["missing_variants"] == 1
    assert plan["missing_fraction"] == pytest.approx(0.25)

    with pytest.raises(ValueError) as error:
        prepare_missing_sample_sizes(
            frame, _mapping(), policies=_policies(maximum=0.249), scope="dataset 'study'",
        )
    message = str(error.value)
    assert "1 of 4 variants (25.00%)" in message
    assert "exceeding sample_size.max_missing_fraction = 24.90%" in message
    assert "Provide a more complete sample-size column" in message


def test_binary_missing_fraction_counts_each_affected_variant_once():
    frame = pl.DataFrame({
        "NCASE": [100, None, None, 100],
        "NCONTROL": [200, 200, None, None],
    })
    plan = prepare_missing_sample_sizes(
        frame,
        _mapping(ncase_col="NCASE", ncontrol_col="NCONTROL"),
        policies=_policies(maximum=0.75),
    )

    assert plan["missing_by_column"] == {"NCASE": 2, "NCONTROL": 2}
    assert plan["missing_variants"] == 3
    assert plan["missing_fraction"] == pytest.approx(0.75)


def test_imputation_statistic_excludes_counts_below_configured_floor():
    policies = default_policies().with_overrides({
        "sample_size.missing_action": "median",
        "sample_size.max_missing_fraction": 0.25,
        "sample_size.min_value": 1,
    })
    plan = prepare_missing_sample_sizes(
        pl.DataFrame({"N": [-100, 100, 300, None]}),
        _mapping(),
        policies=policies,
    )

    assert plan["fill_values"] == {"N": 200}


@pytest.mark.parametrize(
    ("action", "expected"),
    [("median", 200), ("mean", 3433)],
)
def test_dataset_statistic_is_reused_when_chromosome_values_are_filled(action, expected):
    full_dataset = pl.DataFrame({"N": [100, 200, 10_000, None]})
    policies = _policies(action=action, maximum=0.25)
    plan = prepare_missing_sample_sizes(full_dataset, _mapping(), policies=policies)

    chromosome, qc, mapping = harmonise_sample_sizes(
        "2",
        pl.DataFrame({"N": [10_000, None]}),
        _mapping(),
        policies=policies,
        missing_plan=plan,
    )

    assert plan["fill_values"] == {"N": expected}
    assert chromosome["N"].to_list() == [10_000, expected]
    assert chromosome["Neff"].to_list() == [10_000, expected]
    assert qc["missing_sample_size"]["variants_imputed"] == 1
    assert qc["missing_sample_size"]["fill_values"] == {"N": expected}
    assert mapping["neff_col"] == "Neff"


def test_remove_action_records_missing_sample_size_in_reject_output(tmp_path):
    frame = pl.DataFrame({
        "variant": ["kept", "missing"], "N": [100, None],
    }).with_row_index(SOURCE_INPUT_ROW_COLUMN, offset=1)
    policies = _policies(action="remove", maximum=0.5)
    plan = prepare_missing_sample_sizes(frame, _mapping(), policies=policies)
    rejects = RejectCollector(
        source_snapshot=frame,
        logger=None,
        out_path=str(tmp_path / "rejected.tsv"),
        delimiter="\t",
        compress=False,
    )

    result, qc, _mapping_result = harmonise_sample_sizes(
        "1", frame, _mapping(), policies=policies, rejects=rejects,
        missing_plan=plan,
    )

    assert result["variant"].to_list() == ["kept"]
    assert qc["removed_by_reason"]["sample_size_invalid"] == 1
    rejected = rejects.frame()
    assert rejected["variant"].to_list() == ["missing"]
    assert rejected["reject_reason"].to_list() == ["sample_size_invalid"]


def test_fail_action_stops_on_any_missing_value_with_clear_message():
    with pytest.raises(ValueError) as error:
        prepare_missing_sample_sizes(
            pl.DataFrame({"N": [100, None]}),
            _mapping(),
            policies=_policies(action="fail", maximum=1.0),
            scope="dataset 'study'",
        )

    assert "1 of 2 variants (50.00%)" in str(error.value)
    assert "sample_size.missing_action is 'fail'" in str(error.value)


def test_fixed_sample_size_has_no_per_variant_missingness():
    plan = prepare_missing_sample_sizes(
        pl.DataFrame({"variant": ["a", "b"]}),
        _mapping(ncontrol_col=None, ncontrol=50_000),
        policies=default_policies(),
    )

    assert plan["columns"] == []
    assert plan["missing_variants"] == 0


def test_default_minimum_rejects_zero_sample_size():
    result, qc, _mapping_result = harmonise_sample_sizes(
        "1",
        pl.DataFrame({"variant": ["zero", "one"], "N": [0, 1]}),
        _mapping(),
        policies=default_policies(),
    )

    assert result["variant"].to_list() == ["one"]
    assert result["Neff"].to_list() == [1]
    assert qc["removed_by_reason"]["sample_size_invalid"] == 1
