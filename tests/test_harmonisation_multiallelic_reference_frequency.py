"""Complete reference allele sets supply REF AF without corrupting partial sets."""

import polars as pl
import pytest

from postgwas.modules.harmonisation import allele_frequency as frequency_module
from postgwas.modules.harmonisation.allele_frequency import (
    AlleleFrequencyError,
    merge_external_allele_frequencies,
)
from postgwas.modules.harmonisation.effect_validation import final_completeness_check
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.shared.allele_join import allele_oriented_left_join
from postgwas.modules.harmonisation.shared import allele_join as join_module


STUDY = {"chr": "CHR", "pos": "POS", "ea": "EA", "oa": "OA"}
REFERENCE = {"chr": "CHROM", "pos": "POS", "ea": "ALT", "oa": "REF"}
FILE_MAPPING = {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF", "delimiter": "tab"}


def _reference(rows):
    return pl.DataFrame(
        rows, schema={"CHROM": pl.String, "POS": pl.Int64,
                      "REF": pl.String, "ALT": pl.String, "AF": pl.Float64},
        orient="row",
    )


def _study(ea="A", oa="G"):
    return pl.DataFrame({"CHR": ["1"], "POS": [10], "EA": [ea], "OA": [oa],
                         "RAW_MAF": [0.05], "BETA": [0.2], "ROW": [7]})


def _join(reference, study=None, **overrides):
    options = dict(
        study_columns=STUDY, reference_columns=REFERENCE,
        policies=default_policies(), value_column="AF", output_column="EAF",
        duplicate_exact_action="keep_one", duplicate_non_identical_action="discard_all",
        swapped_value="one_minus", prefer_non_null_value=True,
        screen_frequency=False, require_complete_frequency_site=True,
        frequency_reason_column="ASSIGNMENT_REASON",
    )
    options.update(overrides)
    return allele_oriented_left_join(_study() if study is None else study, reference, **options)


@pytest.mark.parametrize("alternates, expected", [
    ([("G", .3), ("C", .1)], .6),
    ([("G", .2), ("C", .1), ("T", .05)], .65),
    ([("G", .4), ("C", .6)], 0.0),
    ([("G", .4791), ("C", 0.0)], .5209),
])
def test_complete_snv_alleles_reconstruct_reference_frequency(alternates, expected):
    source = _study()
    reference = _reference([("1", 10, "A", allele, value) for allele, value in alternates])
    source_before, reference_before = source.clone(), reference.clone()

    result, orientation, stats = _join(reference, source)

    assert result["EAF"].to_list() == pytest.approx([expected])
    assert result[orientation].to_list() == [1]
    assert stats["swapped_value_matches"] == 1
    assert stats["missing_values"] == 0
    assert stats["reference_frequency_reconstructed_rows"] == 1
    assert stats["reference_frequency_unresolved_rows"] == 0
    assert stats["reference_multiallelic_sites"] == 1
    assert "ASSIGNMENT_REASON" not in result.columns
    assert source.equals(source_before)
    assert reference.equals(reference_before)
    assert result.select(source.columns).equals(source_before)


def test_mixed_direct_biallelic_and_multiallelic_rows_preserve_order_and_values():
    source = pl.DataFrame({"CHR": ["1"] * 3, "POS": [30, 10, 20],
                           "EA": ["A", "G", "T"], "OA": ["G", "A", "C"],
                           "RAW_MAF": [.1, .2, .25], "ROW": [3, 1, 2]})
    reference = _reference([
        ("1", 10, "A", "G", .2), ("1", 20, "T", "C", .25),
        ("1", 30, "A", "G", .3), ("1", 30, "A", "C", .1),
    ])

    result, orientation, stats = _join(reference, source)

    assert result["EAF"].to_list() == pytest.approx([.6, .2, .75])
    assert result[orientation].to_list() == [1, 0, 1]
    assert result.select(source.columns).equals(source)
    assert stats["direct_value_matches"] == 1
    assert stats["swapped_value_matches"] == 2
    assert stats["reference_frequency_reconstructed_rows"] == 1
    assert stats["reference_frequency_unresolved_rows"] == 0


def test_exact_duplicate_alt_is_counted_once_in_reference_sum():
    reference = _reference([
        ("1", 10, "A", "G", .3), ("1", 10, "A", "G", .3),
        ("1", 10, "A", "C", .1),
    ])

    result, _, stats = _join(reference)

    assert result["EAF"].to_list() == pytest.approx([.6])
    assert stats["reference_exact_duplicate_groups"] == 1
    assert stats["reference_duplicate_rows_removed"] == 1
    assert stats["reference_frequency_reconstructed_rows"] == 1


@pytest.mark.parametrize("extra", [None, float("nan"), float("inf"), -.01, 1.01])
def test_invalid_alt_component_never_becomes_partial_reference_sum(extra):
    reference = _reference([("1", 10, "A", "G", .3), ("1", 10, "A", "C", extra)])

    result, orientation, stats = _join(reference)

    assert result["EAF"].to_list() == [None]
    assert result[orientation].to_list() == [1]
    assert stats["swapped_key_matches"] == 1
    assert stats["swapped_value_matches"] == 0
    assert stats["missing_values"] == 1
    assert result["RAW_MAF"].to_list() == [.05]
    assert stats["reference_frequency_reconstructed_rows"] == 0
    assert stats["reference_frequency_unresolved_rows"] == 1
    assert result["ASSIGNMENT_REASON"].to_list() == [
        "missing_or_invalid_reference_frequency_component",
    ]


@pytest.mark.parametrize("rows, reason", [
    ([("1", 10, "A", "G", .7), ("1", 10, "A", "C", .4)],
     "reference_alt_frequency_sum_out_of_range"),
    ([("1", 10, "A", "G", .3), ("1", 10, "T", "C", .1)],
     "unsupported_reference_allele_set"),
    ([("1", 10, "A", "G", .3), ("1", 10, "A", "AT", .1)],
     "unsupported_reference_allele_set"),
    ([("1", 10, "A", "G", .3), ("1", 10, "A", "*", .1)],
     "unsupported_reference_allele_set"),
    ([("1", 10, "A", "G", .3), ("1", 10, "A", "A", .1)],
     "unsupported_reference_allele_set"),
])
def test_incompatible_reference_sets_leave_swapped_frequency_unresolved(rows, reason):
    result, orientation, stats = _join(_reference(rows))

    assert result["EAF"].to_list() == [None]
    assert result[orientation].to_list() == [1]
    assert stats["missing_values"] == 1
    assert stats["reference_frequency_unresolved_rows"] == 1
    assert stats["reference_frequency_unresolved_reasons"] == {reason: 1}
    assert result["ASSIGNMENT_REASON"].to_list() == [reason]


@pytest.mark.parametrize("exact_action, extra_values", [
    ("keep_one", [.1, .2]),
    ("discard_all", [.1, .1]),
])
def test_discarded_alt_component_invalidates_remaining_site(exact_action, extra_values):
    # The study matches the unchanged G row; discarded C rows must not turn
    # 1 - G into an apparently valid substitute for the full reference sum.
    reference = _reference([
        ("1", 10, "A", "G", .3),
        *[("1", 10, "A", "C", value) for value in extra_values],
    ])

    result, orientation, stats = _join(reference, duplicate_exact_action=exact_action)

    assert result[orientation].to_list() == [1]
    assert result["EAF"].to_list() == [None]
    assert stats["reference_duplicate_groups_discarded"] == 1
    assert stats["missing_values"] == 1
    assert stats["reference_frequency_unresolved_rows"] == 1
    assert result["ASSIGNMENT_REASON"].to_list() == [
        "discarded_or_conflicting_reference_component",
    ]


def test_direct_alt_does_not_require_reconstructing_missing_other_allele():
    reference = _reference([("1", 10, "A", "G", .3), ("1", 10, "A", "C", None)])

    result, orientation, _ = _join(reference, _study("G", "A"))

    assert result["EAF"].to_list() == [.3]
    assert result[orientation].to_list() == [0]


def test_existing_biallelic_indel_inversion_stays_supported():
    result, orientation, _ = _join(
        _reference([("1", 10, "A", "AT", .2)]), _study("A", "AT"),
    )

    assert result["EAF"].to_list() == pytest.approx([.8])
    assert result[orientation].to_list() == [1]


def test_nonfrequency_annotation_keeps_existing_swapped_value():
    reference = _reference([("1", 10, "A", "G", .3), ("1", 10, "A", "C", .1)])

    result, _, _ = _join(reference, swapped_value="same", require_complete_frequency_site=False)

    assert result["EAF"].to_list() == [.3]


def test_direct_only_frequency_matching_skips_unneeded_site_aggregation(monkeypatch):
    reference = _reference([("1", 10, "A", "G", .3), ("1", 10, "A", "C", .1)])

    def unneeded_aggregation(*args, **kwargs):
        raise AssertionError("Direct-only matching must not reconstruct REF frequencies")

    monkeypatch.setattr(join_module, "_reference_frequency_sites", unneeded_aggregation)

    result, orientation, stats = _join(reference, _study("G", "A"), orientations=["direct"])

    assert result["EAF"].to_list() == [.3]
    assert result[orientation].to_list() == [0]
    assert stats["direct_value_matches"] == 1
    assert "ASSIGNMENT_REASON" not in result.columns


def test_ordinary_external_semantics_exclude_multiallelic_frequency_pairs(tmp_path, monkeypatch):
    reference = tmp_path / "ordinary_external.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tAF\n"
        "1\t10\tA\tG\t0.3\n1\t10\tA\tC\t0.1\n"
        "1\t20\tA\tG\t0.2\n1\t30\tA\tC\t0.3\n"
        "1\t40\tG\tA\t0.4\n1\t50\tT\tC\t0.1\n"
    )
    source = pl.DataFrame({
        "CHR": ["1"] * 6, "POS": [10, 10, 20, 30, 40, 50],
        "EA": ["A", "G", "G", "C", "A", "T"],
        "OA": ["G", "A", "A", "A", "G", "C"],
        "DEFAULT_EAF": [.6, .3, .2, .3, .4, .9],
        "RAW_MAF": [.3, .3, .2, .3, .4, .1],
    })
    original_source = source.clone()
    observed = []
    confirm = frequency_module.confirm_eaf_with_reference

    def record_real_confirmation(frame, *args, **kwargs):
        observed.append(frame.clone())
        return confirm(frame, *args, **kwargs)

    monkeypatch.setattr(frequency_module, "confirm_eaf_with_reference", record_real_confirmation)
    result, selected, stats = merge_external_allele_frequencies(
        source, str(reference), "AF", FILE_MAPPING, STUDY, chromosome="1",
        validation_reference_col="DEFAULT_EAF", declared_alt_frequency=False,
        policies=default_policies().with_overrides({"eaf.maf_reference_min_overlap": 3}),
    )

    assert len(observed) == 1
    # A valid multiallelic attachment does not make 1 - REF_AF a raw ALT AF:
    # at position 10 that would fabricate 0.4 instead of the source G AF 0.3.
    # Both its direct and swapped rows must be excluded from biallelic-only
    # MAF/EAF diagnostic evidence, while retaining their valid assigned EAF.
    assert observed[0]["POS"].to_list() == [20, 30, 40, 50]
    assert stats["maf_reference_validation"]["comparable_variants"] == 4
    assert result[selected].to_list() == pytest.approx([.6, .3, .2, .3, .4, .9])
    assert result.select(source.columns).equals(original_source)


def _coverage_fixture(tmp_path):
    reference = tmp_path / "replacement.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tAF\n"
        "1\t10\tA\tG\t0.3\n1\t10\tA\tC\t0.1\n"
        "1\t20\tA\tG\t0.2\n1\t30\tA\tG\t0.2\n1\t40\tA\tG\t0.2\n"
        "1\t50\tA\tG\t0.3\n1\t50\tA\tC\tNA\n"
    )
    source = pl.DataFrame({"CHR": ["1"] * 5, "POS": [10, 20, 30, 40, 50],
                           "EA": ["A"] * 5, "OA": ["G"] * 5,
                           "RAW_MAF": [.1] * 5})
    return reference, source


def test_unresolved_assignment_stays_in_coverage_denominator_then_missing_policy(tmp_path):
    reference, source = _coverage_fixture(tmp_path)
    result, selected, stats = merge_external_allele_frequencies(
        source, str(reference), "AF", FILE_MAPPING, STUDY,
        chromosome="1", declared_alt_frequency=True,
    )

    assert result.height == 5
    assert result[selected].to_list() == pytest.approx([.6, .8, .8, .8, None])
    assert stats["key_match_fraction"] == 1.0
    assert stats["usable_match_rows"] == 4
    assert stats["usable_match_fraction"] == .8
    assert stats["matched_missing_frequency_rows"] == 1
    assert stats["reference_frequency_reconstructed_rows"] == 1
    assert stats["reference_frequency_unresolved_rows"] == 1
    reasons = [name for name in result.columns if name.startswith("__postgwas_external_eaf_assignment_reason")]
    assert len(reasons) == 1
    assert result[reasons[0]].to_list() == [
        None, None, None, None, "missing_or_invalid_reference_frequency_component",
    ]
    policies = default_policies().with_overrides({"final_check.require": ["eaf"]})
    kept, qc, _ = final_completeness_check(
        "1", result, {"eaf_col": selected}, policies=policies,
    )
    assert kept["POS"].to_list() == [10, 20, 30, 40]
    assert kept["RAW_MAF"].to_list() == [.1] * 4
    assert qc["removed_by_reason"]["final_missing_eaf"] == 1


def test_unresolved_frequency_can_still_fail_existing_coverage_gate(tmp_path):
    reference, source = _coverage_fixture(tmp_path)

    with pytest.raises(AlleleFrequencyError, match="usable EAF.*only 4.*80"):
        merge_external_allele_frequencies(
            source, str(reference), "AF", FILE_MAPPING, STUDY,
            chromosome="1", declared_alt_frequency=True,
            policies=default_policies().with_overrides({"eaf.external_min_match_fraction": .81}),
        )
