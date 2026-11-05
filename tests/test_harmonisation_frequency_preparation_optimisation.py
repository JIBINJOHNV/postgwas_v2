"""AF preparation reuses work while preserving row and duplicate semantics."""

from copy import deepcopy

import polars as pl
from polars.testing import assert_frame_equal
import pytest

from postgwas.modules.harmonisation import allele_frequency as af
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.shared import allele_join
from postgwas.modules.harmonisation.shared.variant_columns import canonicalize_variant_frame
from test_harmonisation_external_af_preflight import fixture, validate


KEYS = {role: role for role in ("chr", "pos", "ea", "oa")}


@pytest.mark.parametrize("exact_action,discarded,positions", [
    ("keep_one", 1, [1, 3, 4, 5, 6, 7]),
    ("discard_all", 4, [5, 6, 7]),
])
def test_duplicate_statistics_use_one_reduction_without_changing_invalid_values(
        monkeypatch, exact_action, discarded, positions):
    reference = pl.DataFrame({"chr": ["1"] * 11,
        "pos": [1, 1, 2, 2, 3, 3, 4, 4, 5, 6, 7],
        "ea": ["A"] * 11, "oa": ["G"] * 11,
        "AF": [.2, .2, .1, None, float("nan"), float("nan"),
               float("inf"), float("inf"), -.1, None, .4]})
    before = reference.clone()
    select = pl.DataFrame.select
    reductions = []

    def measured_select(frame, *args, **kwargs):
        result = select(frame, *args, **kwargs)
        if (any(name.startswith("__postgwas_duplicate_count__") for name in frame.columns)
                and result.height == 1 and not set(KEYS).intersection(result.columns)):
            reductions.append(result.columns)
        return result

    monkeypatch.setattr(pl.DataFrame, "select", measured_select)
    result, statistics = allele_join.deduplicate_reference_values(
        reference, list(KEYS), "AF", exact_action=exact_action,
        non_identical_action="discard_all", reference_label="audit", error_type=RuntimeError)
    assert statistics == {"reference_duplicate_groups": 4, "reference_exact_duplicate_groups": 3,
        "reference_non_identical_duplicate_groups": 1, "reference_duplicate_groups_discarded": discarded,
        "reference_exact_duplicate_rows": 6, "reference_non_identical_duplicate_rows": 2}
    expected = before.filter(pl.col("pos").is_in(positions)).unique(
        subset=list(KEYS), keep="first", maintain_order=True)
    assert_frame_equal(result, expected)
    assert_frame_equal(reference, before)
    assert len(reductions) == 1


def _effect_inputs():
    frame = pl.DataFrame({af.SOURCE_INPUT_ROW_COLUMN: [5, 0, 3, 1, 4, 2],
        "chr": ["chr01"] * 6, "pos": ["106.0", "101", "104", "102", "105", "1.03e2"],
        "ea": [" a "] * 6, "oa": [" g "] * 6,
        "__postgwas_effect_key_ea": ["original source"] * 6})
    matches = pl.DataFrame({af.SOURCE_INPUT_ROW_COLUMN: [0, 1, 2, 2, 4, 5],
        "__strand_candidates": [1, 1, 2, 2, 1, 1],
        "__strand_swapped": [False, True, False, True, False, True],
        "__strand_af": [.3, .3, .1, .2, None, .3],
        "__strand_ref_af": [.6, .6, .7, .7, None, .6],
        "__strand_alt": ["A"] * 6, "__strand_ref": ["G"] * 6,
        af.STRAND_ACTION_COLUMN: ["forward", "forward_swapped", "forward", "forward_swapped", "forward", "forward_swapped"]})
    return frame, matches


def test_same_source_effect_view_avoids_unused_membership_join_and_keeps_all_source_rows(monkeypatch):
    frame, matches = _effect_inputs()
    before = frame.clone()
    join = pl.DataFrame.join
    joins = []

    def measured_join(frame, *args, **kwargs):
        joins.append(kwargs.get("on"))
        return join(frame, *args, **kwargs)

    monkeypatch.setattr(pl.DataFrame, "join", measured_join)
    result = af._effect_frequency_view(frame, matches, KEYS, default_policies())
    assert result[af.SOURCE_INPUT_ROW_COLUMN].to_list() == before[af.SOURCE_INPUT_ROW_COLUMN].to_list()
    assert result["effect_frequency"].to_list() == [.6, .3, None, .6, None, None]
    assert result["default_frequency"].to_list() == result["effect_frequency"].to_list()
    assert_frame_equal(frame, before)
    assert len(joins) == 1


@pytest.mark.parametrize("helper", ["effect", "candidates"])
def test_helpers_can_reuse_an_explicitly_canonical_full_key_frame(monkeypatch, helper):
    frame, matches = _effect_inputs()
    policies = default_policies().with_overrides({"chromosome.rename_map": {"1": "2", "2": "3"}})
    canonical, _ = canonicalize_variant_frame(frame, KEYS, policies)
    if helper == "effect":
        function = af._effect_frequency_view
    else:
        function = af._replacement_candidate_keys
    expected = function(frame, matches, KEYS, policies)
    monkeypatch.setattr(af, "canonicalize_variant_frame",
        lambda *a, **k: pytest.fail("canonical key projection was normalized again"))
    actual = function(canonical, matches, KEYS, policies, study_columns_canonical=True)
    assert_frame_equal(actual, expected)
    assert actual[af.SOURCE_INPUT_ROW_COLUMN].n_unique() == frame.height


@pytest.mark.parametrize("same_source", [False, True])
def test_real_external_preparation_normalizes_original_keys_once_and_keeps_unmatched_rows(
        tmp_path, monkeypatch, same_source):
    values = [.6 + i % 30 / 100 for i in range(1200)]

    def remove_default_rows(source, default, external):
        if not same_source:
            pl.read_csv(default, separator="\t").filter(pl.col("POS") < 1100).write_csv(default, separator="\t")
        return source

    parts = list(fixture(tmp_path, values=values, default_values=values,
        modify=remove_default_rows, policies={"strand.unmatched_action": "retain"}))
    if same_source:
        parts[3] = {chromosome: dict(resources, user_eaf_file=resources["default_eaf_file"])
                    for chromosome, resources in parts[3].items()}
    normalize = af.canonicalize_variant_frame
    labels = []

    def measured_normalization(*args, **kwargs):
        labels.append(kwargs.get("label"))
        return normalize(*args, **kwargs)

    monkeypatch.setattr(af, "canonicalize_variant_frame", measured_normalization)
    original = deepcopy(parts[2])
    selected = validate(parts)
    screen = selected["eaf_external_validation"]["effective_frequency_screen"]
    assert screen["rows"] == screen["usable"] == 1200
    assert screen["at_or_below_0_5"] == 0
    assert selected["eaf_external_validation"]["status"] == "passed"
    assert parts[2] == original
    assert labels.count("replacement frequency diagnostic keys") == 2
    assert "external AF assignment keys" not in labels
    assert "effect-frequency screen keys" not in labels
