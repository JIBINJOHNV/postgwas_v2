"""Focused tests for reference confirmation of MAF-like frequency columns."""

from unittest.mock import patch

import polars as pl
import pytest

from postgwas.modules.harmonisation.allele_frequency import (
    AlleleFrequencyError,
    confirm_eaf_with_reference,
    harmonise_allele_frequency,
    merge_external_allele_frequencies,
    validate_allele_frequency,
)
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.shared.statistics import (
    frequency_maf_screen,
)


MAPPING = {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF", "delimiter": "tab"}
STANDARD = {"chr": "CHR", "pos": "POS", "ea": "EA", "oa": "OA"}


def _reference(tmp_path):
    path = tmp_path / "GRCh37_1000G_freq_chr1.tsv"
    path.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.8\n"
        "1\t200\tC\tT\t0.3\n",
        encoding="utf-8",
    )
    return str(path)


def _policies(minimum):
    return default_policies().with_overrides(
        {"eaf": {"maf_reference_min_overlap": minimum}}
    )


def _raw_frequency_screen(*values):
    return frequency_maf_screen(
        pl.DataFrame({"AF": list(values)}),
        "AF",
        float(default_policies().get("eaf.maf_decision_cutoff")),
    )


def _external_main_contract(tmp_path):
    return (
        {
            "gwas_outputname": "study",
            "output_folder": str(tmp_path),
            "chr_col": "CHR",
            "pos_col": "POS",
            "ea_col": "EA",
            "oa_col": "OA",
            "eaf_col": None,
            "eafcolumn": "EUR",
        },
        {
            "frequency_qc_directory": "eaf_qc",
            "missing_eaf": "eaf_qc/{dataset_id}_missing_chr{chromosome}.tsv",
            "out_of_range_eaf": "eaf_qc/{dataset_id}_range_chr{chromosome}.tsv",
            "post_orientation_duplicates": (
                "qc/{dataset_id}_duplicates_chr{chromosome}.tsv"
            ),
        },
    )


def test_out_of_range_eaf_is_rejected_by_default_without_silent_clipping():
    frame = pl.DataFrame({
        "EAF": [-0.1, 0.0, 0.2, 1.0, 1.02, float("nan"), float("inf")],
    })

    valid, _maf_like, count, missing, result, _low, stats = (
        validate_allele_frequency(frame, "EAF", policies=default_policies())
    )

    assert not valid
    assert count == 4
    assert missing == 0
    assert result["EAF"].to_list() == [0.0, 0.2, 1.0]
    assert stats["out_of_range_action"] == "reject"
    assert stats["initial"]["non_finite_rows"] == 2
    assert stats["initial"]["usable_rows"] == 3
    assert stats["initial"]["le_0_5_rows"] == 2
    assert stats["initial"]["le_0_5_fraction_among_usable"] == pytest.approx(2 / 3)
    assert stats["final"]["invalid_total_rows"] == 0


def test_missing_eaf_is_not_misclassified_as_out_of_range():
    frame = pl.DataFrame({"EAF": [None, -0.1, 0.2]})

    _valid, _maf_like, count, missing, result, _low, stats = (
        validate_allele_frequency(frame, "EAF", policies=default_policies())
    )

    assert count == 1
    assert missing == 1
    assert result["EAF"].to_list() == [None, 0.2]
    assert stats["initial"]["null_rows"] == 1
    assert stats["initial"]["invalid_total_rows"] == 1


def test_null_out_of_range_policy_blanks_invalid_values_and_does_not_clip_endpoints():
    policies = default_policies().with_overrides({"eaf.out_of_range": "null"})

    _valid, _maf_like, count, missing, result, _low, _stats = (
        validate_allele_frequency(
            pl.DataFrame({"EAF": [-0.1, 0.0, 0.2, 1.0, 1.02]}),
            "EAF",
            policies=policies,
        )
    )

    assert count == 2
    assert missing == 2
    assert result["EAF"].to_list() == [None, 0.0, 0.2, 1.0, None]


def test_clipping_eaf_requires_explicit_policy():
    policies = default_policies().with_overrides({"eaf.out_of_range": "clip"})

    _valid, _maf_like, count, missing, result, _low, _stats = (
        validate_allele_frequency(
            pl.DataFrame({"EAF": [-0.1, 0.2, 1.02, float("inf")]}),
            "EAF",
            policies=policies,
        )
    )

    assert count == 3
    assert missing == 1
    assert result["EAF"].to_list() == [0.0, 0.2, 1.0, None]


def test_reference_check_confirms_unaligned_maf_after_direct_and_swapped_matches(tmp_path):
    study = pl.DataFrame({
        "CHR": ["1", "1"], "POS": [100, 200],
        "EA": ["A", "C"], "OA": ["G", "T"], "AF": [0.2, 0.3],
    })

    decision, stats = confirm_eaf_with_reference(
        study, "AF", _reference(tmp_path), "EUR", MAPPING, STANDARD,
        policies=_policies(2),
    )

    assert decision == "maf"
    assert stats["direct_matches"] == 1
    assert stats["swapped_matches"] == 1
    assert stats["reference_effect_allele_minor_fraction"] == pytest.approx(0.0)
    assert stats["mean_maf_error"] == pytest.approx(0.0)
    assert stats["mean_eaf_error"] == pytest.approx(0.5)
    assert stats["maf_correlation"] == pytest.approx(1.0)
    assert stats["mean_absolute_maf_difference"] == pytest.approx(0.0)
    assert stats["reference_suitability_status"] == "suitable"


def test_reference_check_is_inconclusive_when_reference_effect_alleles_are_all_minor(tmp_path):
    study = pl.DataFrame({
        "CHR": ["1", "1"], "POS": [100, 200],
        "EA": ["G", "T"], "OA": ["A", "C"], "AF": [0.2, 0.3],
    })

    decision, stats = confirm_eaf_with_reference(
        study, "AF", _reference(tmp_path), "EUR", MAPPING, STANDARD,
        policies=_policies(2),
    )

    assert decision == "inconclusive"
    assert stats["reference_effect_allele_minor_fraction"] == pytest.approx(1.0)
    assert stats["mean_eaf_error"] == pytest.approx(0.0)
    assert stats["mean_maf_error"] == pytest.approx(0.0)
    assert stats["decision_reason"] == "reference_effect_allele_mostly_minor"


def test_reference_check_confirms_eaf_only_with_informative_separation(tmp_path):
    study = pl.DataFrame({
        "CHR": ["1", "1"], "POS": [100, 200],
        "EA": ["A", "C"], "OA": ["G", "T"], "AF": [0.8, 0.7],
    })

    decision, stats = confirm_eaf_with_reference(
        study, "AF", _reference(tmp_path), "EUR", MAPPING, STANDARD,
        policies=_policies(2),
    )

    assert decision == "eaf"
    assert stats["reference_effect_allele_minor_fraction"] == pytest.approx(0.0)
    assert stats["mean_eaf_error"] == pytest.approx(0.0)
    assert stats["mean_maf_error"] == pytest.approx(0.5)
    assert stats["decision_reason"] == "eaf_error_lower_by_required_margin"


def test_reference_check_is_inconclusive_when_errors_are_too_close(tmp_path):
    reference = tmp_path / "reference.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.51\n"
        "1\t200\tC\tT\t0.49\n"
        "1\t300\tT\tC\t0.52\n",
        encoding="utf-8",
    )
    study = pl.DataFrame({
        "CHR": ["1", "1", "1"], "POS": [100, 200, 300],
        "EA": ["A", "T", "C"], "OA": ["G", "C", "T"],
        "AF": [0.50, 0.49, 0.48],
    })

    decision, stats = confirm_eaf_with_reference(
        study, "AF", str(reference), "EUR", MAPPING, STANDARD,
        policies=_policies(3),
    )

    assert decision == "inconclusive"
    assert stats["mean_eaf_error"] == pytest.approx(0.05 / 3)
    assert stats["mean_maf_error"] == pytest.approx(0.01 / 3)
    assert stats["absolute_error_difference"] == pytest.approx(0.04 / 3)
    assert stats["minimum_error_margin"] == pytest.approx(0.02)
    assert stats["decision_reason"] == "error_difference_below_required_margin"

    strict_policies = _policies(3).with_overrides({
        "eaf.maf_reference_min_correlation": 0.90,
    })
    strict_decision, strict_stats = confirm_eaf_with_reference(
        study,
        "AF",
        str(reference),
        "EUR",
        MAPPING,
        STANDARD,
        policies=strict_policies,
    )
    assert strict_decision == "inconclusive"
    assert strict_stats["reference_suitability_reason"] == (
        "maf_correlation_below_minimum"
    )


def test_reference_check_stops_when_folded_maf_correlation_is_too_low(tmp_path):
    reference = tmp_path / "reference.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.30\n"
        "1\t200\tC\tT\t0.10\n"
        "1\t300\tT\tC\t0.20\n",
        encoding="utf-8",
    )
    study = pl.DataFrame({
        "CHR": ["1", "1", "1"],
        "POS": [100, 200, 300],
        "EA": ["A", "T", "C"],
        "OA": ["G", "C", "T"],
        "AF": [0.10, 0.20, 0.30],
    })

    decision, stats = confirm_eaf_with_reference(
        study,
        "AF",
        str(reference),
        "EUR",
        MAPPING,
        STANDARD,
        policies=_policies(3),
    )

    assert decision == "inconclusive"
    assert stats["maf_correlation"] == pytest.approx(-0.5)
    assert stats["reference_suitability_status"] == "inconclusive"
    assert stats["reference_suitability_reason"] == (
        "maf_correlation_below_minimum"
    )


def test_reference_check_stops_when_folded_maf_difference_is_too_large(
    tmp_path,
):
    reference = tmp_path / "reference.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.25\n"
        "1\t200\tC\tT\t0.35\n"
        "1\t300\tT\tC\t0.45\n",
        encoding="utf-8",
    )
    study = pl.DataFrame({
        "CHR": ["1", "1", "1"],
        "POS": [100, 200, 300],
        "EA": ["A", "T", "C"],
        "OA": ["G", "C", "T"],
        "AF": [0.10, 0.20, 0.30],
    })

    decision, stats = confirm_eaf_with_reference(
        study,
        "AF",
        str(reference),
        "EUR",
        MAPPING,
        STANDARD,
        policies=_policies(3),
    )

    assert decision == "inconclusive"
    assert stats["maf_correlation"] == pytest.approx(1.0)
    assert stats["mean_absolute_maf_difference"] == pytest.approx(0.15)
    assert stats["reference_suitability_reason"] == (
        "maf_difference_above_maximum"
    )

    relaxed_policies = _policies(3).with_overrides({
        "eaf.maf_reference_max_mean_absolute_difference": 0.16,
    })
    _decision, relaxed_stats = confirm_eaf_with_reference(
        study,
        "AF",
        str(reference),
        "EUR",
        MAPPING,
        STANDARD,
        policies=relaxed_policies,
    )
    assert relaxed_stats["reference_suitability_status"] == "suitable"
    assert relaxed_stats["reference_suitability_reason"] == (
        "maf_correlation_and_difference_passed"
    )


def test_reference_check_uses_configured_spearman_correlation(tmp_path):
    policies = default_policies().with_overrides({
        "eaf": {
            "maf_reference_min_overlap": 2,
            "maf_reference_correlation_method": "spearman",
        },
    })
    study = pl.DataFrame({
        "CHR": ["1", "1"],
        "POS": [100, 200],
        "EA": ["A", "C"],
        "OA": ["G", "T"],
        "AF": [0.20, 0.30],
    })

    decision, stats = confirm_eaf_with_reference(
        study,
        "AF",
        _reference(tmp_path),
        "EUR",
        MAPPING,
        STANDARD,
        policies=policies,
    )

    assert decision == "maf"
    assert stats["maf_correlation_method"] == "spearman"
    assert stats["maf_correlation"] == pytest.approx(1.0)


def test_reference_check_is_inconclusive_below_configured_overlap(tmp_path):
    study = pl.DataFrame({
        "CHR": ["1"], "POS": [100], "EA": ["A"], "OA": ["G"], "AF": [0.2],
    })

    decision, stats = confirm_eaf_with_reference(
        study, "AF", _reference(tmp_path), "EUR", MAPPING, STANDARD,
        policies=_policies(2),
    )

    assert decision == "inconclusive"
    assert stats["comparable_variants"] == 1


def test_reference_check_excludes_palindromic_variants(tmp_path):
    reference = tmp_path / "reference.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n1\t100\tT\tA\t0.8\n",
        encoding="utf-8",
    )
    study = pl.DataFrame({
        "CHR": ["1"], "POS": [100], "EA": ["A"], "OA": ["T"], "AF": [0.2],
    })

    decision, stats = confirm_eaf_with_reference(
        study, "AF", str(reference), "EUR", MAPPING, STANDARD,
        policies=_policies(2),
    )

    assert decision == "inconclusive"
    assert stats["comparable_variants"] == 0


def test_reference_check_rejects_invalid_reference_frequency(tmp_path):
    reference = tmp_path / "reference.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n1\t100\tG\tA\t1.2\n",
        encoding="utf-8",
    )
    study = pl.DataFrame({
        "CHR": ["1"], "POS": [100], "EA": ["A"], "OA": ["G"], "AF": [0.2],
    })

    with pytest.raises(AlleleFrequencyError, match="outside 0 to 1"):
        confirm_eaf_with_reference(
            study, "AF", str(reference), "EUR", MAPPING, STANDARD,
            policies=_policies(2),
        )


def test_external_eaf_duplicate_missing_and_finite_values_discards_both(tmp_path):
    reference = tmp_path / "external_eaf.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\tNA\n"
        "1\t100\tG\tA\t0.25\n",
        encoding="utf-8",
    )
    study = pl.DataFrame({
        "CHR": ["1"], "POS": [100], "EA": ["A"], "OA": ["G"],
    })

    result, column, stats = merge_external_allele_frequencies(
        study,
        str(reference),
        "EUR",
        MAPPING,
        STANDARD,
        policies=default_policies().with_overrides({
            "eaf.external_min_match_fraction": 0.0,
        }),
        study_columns_canonical=True,
    )

    assert result[column].to_list() == [None]
    assert stats["reference_duplicate_rows_removed"] == 2
    assert stats["reference_non_identical_duplicate_groups"] == 1


def test_external_eaf_matches_policy_mapped_reference_chromosome(tmp_path):
    reference = tmp_path / "external_eaf_x.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n23\t100\tG\tA\t0.25\n",
        encoding="utf-8",
    )
    study = pl.DataFrame({
        "CHR": ["X"], "POS": [100], "EA": ["A"], "OA": ["G"],
    })

    result, column, stats = merge_external_allele_frequencies(
        study,
        str(reference),
        "EUR",
        MAPPING,
        STANDARD,
        policies=default_policies(),
        study_columns_canonical=True,
    )

    assert result[column].to_list() == pytest.approx([0.25])
    assert stats["usable_direct_match_rows"] == 1
    assert stats["unmatched_rows"] == 0


def test_external_eaf_non_identical_duplicates_are_all_discarded(tmp_path):
    reference = tmp_path / "external_eaf.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.25\n"
        "1\t100\tG\tA\t0.45\n",
        encoding="utf-8",
    )
    study = pl.DataFrame({
        "CHR": ["1"], "POS": [100], "EA": ["A"], "OA": ["G"],
    })

    result, column, stats = merge_external_allele_frequencies(
        study,
        str(reference),
        "EUR",
        MAPPING,
        STANDARD,
        policies=default_policies().with_overrides({
            "eaf.external_min_match_fraction": 0.0,
        }),
        study_columns_canonical=True,
    )

    assert result[column].to_list() == [None]
    assert stats["reference_duplicate_rows_removed"] == 2
    assert stats["reference_non_identical_duplicate_groups"] == 1


def test_external_eaf_exact_duplicates_keep_one(tmp_path):
    reference = tmp_path / "external_eaf.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.25\n"
        "1\t100\tG\tA\t0.25\n",
        encoding="utf-8",
    )
    study = pl.DataFrame({
        "CHR": ["1"], "POS": [100], "EA": ["A"], "OA": ["G"],
    })

    result, column, stats = merge_external_allele_frequencies(
        study,
        str(reference),
        "EUR",
        MAPPING,
        STANDARD,
        policies=default_policies(),
        study_columns_canonical=True,
    )

    assert result[column].to_list() == pytest.approx([0.25])
    assert stats["reference_duplicate_rows_removed"] == 1
    assert stats["reference_exact_duplicate_groups"] == 1


def test_external_palindromic_alt_frequency_uses_declared_mapping_without_changing_effect(tmp_path):
    reference = tmp_path / "external_eaf.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n1\t100\tT\tA\t0.90\n",
        encoding="utf-8",
    )
    study = pl.DataFrame({
        "CHR": ["1"], "POS": [100], "EA": ["A"], "OA": ["T"],
        "BETA": [0.25], "Z": [2.5], "strand_action": ["forward"],
        "strand_reference_af": [0.10],
    })

    result, column, stats = merge_external_allele_frequencies(
        study,
        str(reference),
        "EUR",
        MAPPING,
        STANDARD,
        policies=default_policies(),
        study_columns_canonical=True,
    )

    assert result[column].to_list() == pytest.approx([0.90])
    assert result.item(0, "BETA") == pytest.approx(0.25)
    assert result.item(0, "Z") == pytest.approx(2.5)
    assert result.item(0, "strand_action") == "forward"
    assert "eaf_palindromic_flag" not in result.columns
    assert stats["palindromic_evaluated_rows"] == 1
    assert stats["palindromic_orientation_resolved_rows"] == 1
    assert stats["palindromic_orientation_rejected_rows"] == 0
    assert stats["palindromic_orientation_state_counts"] == {
        "resolved_by_study_orientation": 1
    }


def test_external_palindromic_eaf_near_half_is_retained_after_orientation(tmp_path):
    reference = tmp_path / "external_eaf.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n1\t100\tT\tA\t0.49\n",
        encoding="utf-8",
    )
    study = pl.DataFrame({
        "CHR": ["1"], "POS": [100], "EA": ["A"], "OA": ["T"],
        "BETA": [0.25], "strand_action": ["forward"],
        "strand_reference_af": [0.48],
    })
    policies = default_policies().with_overrides({
        "eaf.external_min_match_fraction": 0.0,
    })

    result, _column, stats = merge_external_allele_frequencies(
        study,
        str(reference),
        "EUR",
        MAPPING,
        STANDARD,
        policies=policies,
        study_columns_canonical=True,
    )

    assert result.height == 1
    assert result.item(0, "EUR") == pytest.approx(0.49)
    assert result.item(0, "BETA") == pytest.approx(0.25)
    assert stats["palindromic_orientation_state_counts"] == {
        "resolved_by_study_orientation": 1
    }
    assert stats["palindromic_orientation_rejected_rows"] == 0


def test_external_palindromic_eaf_cannot_replace_missing_study_orientation(tmp_path):
    reference = tmp_path / "external_eaf.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n1\t100\tT\tA\t0.10\n",
        encoding="utf-8",
    )
    study = pl.DataFrame({
        "CHR": ["1"], "POS": [100], "EA": ["A"], "OA": ["T"],
        "BETA": [0.25], "strand_reference_af": [0.10],
    })
    policies = default_policies().with_overrides({
        "eaf.external_min_match_fraction": 0.0,
    })

    result, _column, stats = merge_external_allele_frequencies(
        study,
        str(reference),
        "EUR",
        MAPPING,
        STANDARD,
        policies=policies,
        study_columns_canonical=True,
    )

    assert result.is_empty()
    assert stats["palindromic_orientation_state_counts"] == {
        "orientation_unavailable": 1
    }


def test_external_eaf_match_fraction_counts_only_usable_frequencies(tmp_path):
    reference = tmp_path / "external_eaf.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.2\n"
        "1\t200\tT\tC\tNA\n"
        "1\t300\tG\tA\t0.3\n"
        "1\t400\tC\tT\tNaN\n"
        "1\t500\tA\tC\tinf\n"
        "1\t600\tT\tG\t1.2\n",
        encoding="utf-8",
    )
    study = pl.DataFrame({
        "CHR": ["1"] * 7,
        "POS": [100, 200, 300, 400, 500, 600, 700],
        "EA": ["A", "C", "G", "T", "C", "G", "A"],
        "OA": ["G", "T", "A", "C", "A", "T", "C"],
    })
    policies = default_policies().with_overrides({
        "eaf": {"external_min_match_fraction": 0.0},
    })

    result, column, stats = merge_external_allele_frequencies(
        study,
        str(reference),
        "EUR",
        MAPPING,
        STANDARD,
        policies=policies,
        study_columns_canonical=True,
    )

    assert result.height == study.height
    assert result.get_column("POS").to_list() == study.get_column("POS").to_list()
    assert column == "EUR"
    assert stats["key_match_rows"] == 6
    assert stats["key_match_fraction"] == pytest.approx(6 / 7)
    assert stats["usable_direct_match_rows"] == 1
    assert stats["usable_flip_match_rows"] == 1
    assert stats["usable_match_rows"] == 2
    assert stats["usable_match_fraction"] == pytest.approx(2 / 7)
    assert stats["match_fraction"] == pytest.approx(2 / 7)
    assert stats["matched_missing_frequency_rows"] == 1
    assert stats["matched_non_finite_frequency_rows"] == 2
    assert stats["matched_out_of_range_frequency_rows"] == 1
    assert stats["unmatched_rows"] == 1


def test_external_eaf_default_fails_below_eighty_percent_usable(tmp_path):
    reference = tmp_path / "external_eaf.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.2\n"
        "1\t200\tT\tC\tNA\n"
        "1\t300\tA\tG\t0.3\n"
        "1\t400\tC\tT\t0.4\n"
        "1\t500\tA\tC\tNA\n",
        encoding="utf-8",
    )
    study = pl.DataFrame({
        "CHR": ["1"] * 5,
        "POS": [100, 200, 300, 400, 500],
        "EA": ["A", "C", "G", "T", "C"],
        "OA": ["G", "T", "A", "C", "A"],
    })

    assert default_policies().get("eaf.external_min_match_fraction") == pytest.approx(0.8)
    with pytest.raises(
        AlleleFrequencyError,
        match=r"alleles for 5 of 5.*usable EAF for only 3 variants \(60.00%\)",
    ):
        merge_external_allele_frequencies(
            study,
            str(reference),
            "EUR",
            MAPPING,
            STANDARD,
            policies=default_policies(),
            study_columns_canonical=True,
        )


def test_external_eaf_exactly_eighty_percent_usable_passes(tmp_path):
    reference = tmp_path / "external_eaf.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.2\n"
        "1\t200\tT\tC\t0.3\n"
        "1\t300\tA\tG\t0.4\n"
        "1\t400\tC\tT\t0.5\n",
        encoding="utf-8",
    )
    study = pl.DataFrame({
        "CHR": ["1"] * 5,
        "POS": [100, 200, 300, 400, 500],
        "EA": ["A", "C", "G", "T", "C"],
        "OA": ["G", "T", "A", "C", "A"],
    })

    _result, _column, stats = merge_external_allele_frequencies(
        study,
        str(reference),
        "EUR",
        MAPPING,
        STANDARD,
        policies=default_policies(),
        study_columns_canonical=True,
    )

    assert stats["usable_match_rows"] == 4
    assert stats["usable_match_fraction"] == pytest.approx(0.8)
    assert stats["match_fraction"] == pytest.approx(0.8)


def test_external_maf_screen_uses_raw_values_before_swapped_alignment(tmp_path):
    reference = tmp_path / "external_eaf.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.2\n"
        "1\t200\tT\tC\t0.2\n",
        encoding="utf-8",
    )
    study = pl.DataFrame({
        "CHR": ["1", "1"],
        "POS": [100, 200],
        "EA": ["A", "T"],
        "OA": ["G", "C"],
    })

    result, column, stats = merge_external_allele_frequencies(
        study,
        str(reference),
        "EUR",
        MAPPING,
        STANDARD,
        policies=_policies(2),
        study_columns_canonical=True,
    )

    assert result[column].to_list() == pytest.approx([0.2, 0.8])
    assert stats["raw_frequency_screen"]["at_or_below_0_5"] == 2
    assert stats["raw_frequency_screen"]["low_fraction_of_usable"] == pytest.approx(1.0)
    assert stats["raw_frequency_screen"]["maf_like"] is True


def test_main_reuses_identical_explicit_external_eaf_and_strand_reference(tmp_path):
    reference = tmp_path / "reference.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.25\n"
        "1\t200\tT\tC\t0.75\n",
        encoding="utf-8",
    )
    external_alias = tmp_path / "external_alias.tsv"
    external_alias.symlink_to(reference)
    oriented = pl.DataFrame({
        "CHR": ["1", "1"],
        "POS": [100, 200],
        "EA": ["A", "C"],
        "OA": ["G", "T"],
        "strand_action": ["forward", "forward_swapped"],
        "strand_reference_af": [0.25, 0.75],
    })
    columns, layout = _external_main_contract(tmp_path)
    strand_qc = {
        "status": "mocked",
        "actions": {"forward": 1, "forward_swapped": 1},
        "raw_frequency_screen": _raw_frequency_screen(0.25, 0.75),
    }

    with patch(
        "postgwas.modules.harmonisation.allele_frequency.harmonise_strand_orientation",
        return_value=(oriented, strand_qc, dict(columns)),
    ), patch(
        "postgwas.modules.harmonisation.allele_frequency.merge_external_allele_frequencies"
    ) as merge:
        result, qc, resolved_columns = harmonise_allele_frequency(
            "1",
            oriented.drop(["strand_action", "strand_reference_af"]),
            dict(columns),
            layout,
            "\t",
            eaffile=str(external_alias),
            external_eaf_colmap={**MAPPING, "delimiter": "auto"},
            default_eaf_file=str(reference),
            default_eaf_column="EUR",
            default_eaf_colmap=dict(MAPPING),
            study_decision={"eaf_is_maf": None},
        )

    merge.assert_not_called()
    assert result["strand_reference_af"].to_list() == pytest.approx([0.25, 0.75])
    assert resolved_columns["eaf_col"] == "strand_reference_af"
    assert qc["external_merge_reused_aligned_strand_reference"] is True
    assert qc["external_merge_usable_match_fraction"] == pytest.approx(1.0)
    assert qc["decision_source"] == "external_eaf"


def test_reused_identical_external_eaf_keeps_usable_coverage_gate(tmp_path):
    reference = tmp_path / "reference.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.25\n"
        "1\t999\tG\tA\t0.75\n",
        encoding="utf-8",
    )
    oriented = pl.DataFrame({
        "CHR": ["1", "1"],
        "POS": [100, 200],
        "EA": ["A", "C"],
        "OA": ["G", "T"],
        "strand_action": ["forward", "forward"],
        "strand_reference_af": [0.25, None],
    })
    columns, layout = _external_main_contract(tmp_path)

    with patch(
        "postgwas.modules.harmonisation.allele_frequency.harmonise_strand_orientation",
        return_value=(
            oriented,
            {
                "status": "mocked",
                "actions": {"forward": 2},
                "raw_frequency_screen": _raw_frequency_screen(0.25, 0.75),
            },
            dict(columns),
        ),
    ), pytest.raises(
        AlleleFrequencyError,
        match=r"usable EAF for only 1 of 2 retained variants.*50.00%.*below the 80.00%",
    ):
        harmonise_allele_frequency(
            "1",
            oriented.drop(["strand_action", "strand_reference_af"]),
            dict(columns),
            layout,
            "\t",
            eaffile=str(reference),
            external_eaf_colmap={**MAPPING, "delimiter": "auto"},
            default_eaf_file=str(reference),
            default_eaf_column="EUR",
            default_eaf_colmap=dict(MAPPING),
            study_decision={"eaf_is_maf": None},
        )


def test_same_external_and_default_file_stops_only_when_raw_column_is_maf_like(
    tmp_path,
):
    reference = tmp_path / "reference.tsv"
    reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.20\n"
        "1\t200\tT\tC\t0.30\n",
        encoding="utf-8",
    )
    oriented = pl.DataFrame({
        "CHR": ["1", "1"],
        "POS": [100, 200],
        "EA": ["A", "C"],
        "OA": ["G", "T"],
        "strand_action": ["forward", "forward"],
        "strand_reference_af": [0.20, 0.30],
    })
    columns, layout = _external_main_contract(tmp_path)
    strand_qc = {
        "status": "mocked",
        "actions": {"forward": 2},
        "raw_frequency_screen": _raw_frequency_screen(0.20, 0.30),
    }

    with patch(
        "postgwas.modules.harmonisation.allele_frequency.harmonise_strand_orientation",
        return_value=(oriented, strand_qc, dict(columns)),
    ), pytest.raises(
        AlleleFrequencyError,
        match=(
            r"external AF column 'EUR' is MAF-like.*resolve to the same file.*"
            r"modules\.harmonisation\.default_eaf\.source"
        ),
    ):
        harmonise_allele_frequency(
            "1",
            oriented.drop(["strand_action", "strand_reference_af"]),
            dict(columns),
            layout,
            "\t",
            eaffile=str(reference),
            external_eaf_colmap={**MAPPING, "delimiter": "auto"},
            default_eaf_file=str(reference),
            default_eaf_column="EUR",
            default_eaf_colmap=dict(MAPPING),
            study_decision={"eaf_is_maf": None},
        )


def test_independent_default_reference_confirms_external_maf_and_stops(tmp_path):
    external = tmp_path / "external.tsv"
    external.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.20\n"
        "1\t200\tT\tC\t0.30\n",
        encoding="utf-8",
    )
    default_reference = tmp_path / "default.tsv"
    default_reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.80\n"
        "1\t200\tT\tC\t0.70\n",
        encoding="utf-8",
    )
    oriented = pl.DataFrame({
        "CHR": ["1", "1"],
        "POS": [100, 200],
        "EA": ["A", "C"],
        "OA": ["G", "T"],
        "strand_action": ["forward", "forward"],
        "strand_reference_af": [0.80, 0.70],
    })
    columns, layout = _external_main_contract(tmp_path)

    with patch(
        "postgwas.modules.harmonisation.allele_frequency.harmonise_strand_orientation",
        return_value=(oriented, {"status": "mocked"}, dict(columns)),
    ), pytest.raises(
        AlleleFrequencyError,
        match=r"external AF column 'EUR' was confirmed as minor allele frequency",
    ):
        harmonise_allele_frequency(
            "1",
            oriented.drop(["strand_action", "strand_reference_af"]),
            dict(columns),
            layout,
            "\t",
            eaffile=str(external),
            external_eaf_colmap=dict(MAPPING),
            default_eaf_file=str(default_reference),
            default_eaf_column="EUR",
            default_eaf_colmap=dict(MAPPING),
            study_decision={"eaf_is_maf": None},
            policies=_policies(2),
        )


def test_independent_default_reference_confirms_external_eaf_and_continues(
    tmp_path,
):
    row_count = 100
    positions = list(range(1, row_count + 1))
    minor_frequencies = [
        0.05 + (index % 10) * 0.02 for index in range(row_count)
    ]
    minor_frequencies[95] = 0.20
    minor_frequencies[96:] = [0.05] * 4
    candidate = minor_frequencies[:96] + [
        1.0 - value for value in minor_frequencies[96:]
    ]
    default = minor_frequencies[:95] + [
        1.0 - value for value in minor_frequencies[95:]
    ]
    external = tmp_path / "external.tsv"
    default_reference = tmp_path / "default.tsv"

    def reference_text(values):
        rows = ["CHROM\tPOS\tREF\tALT\tEUR"]
        rows.extend(
            "1\t%s\tG\tA\t%s" % (position, value)
            for position, value in zip(positions, values)
        )
        return "\n".join(rows) + "\n"

    external.write_text(reference_text(candidate), encoding="utf-8")
    default_reference.write_text(reference_text(default), encoding="utf-8")
    oriented = pl.DataFrame({
        "CHR": ["1"] * row_count,
        "POS": positions,
        "EA": ["A"] * row_count,
        "OA": ["G"] * row_count,
        "strand_action": ["forward"] * row_count,
        "strand_reference_af": default,
    })
    columns, layout = _external_main_contract(tmp_path)

    with patch(
        "postgwas.modules.harmonisation.allele_frequency.harmonise_strand_orientation",
        return_value=(oriented, {"status": "mocked"}, dict(columns)),
    ):
        result, qc, resolved_columns = harmonise_allele_frequency(
            "1",
            oriented.drop(["strand_action", "strand_reference_af"]),
            dict(columns),
            layout,
            "\t",
            eaffile=str(external),
            external_eaf_colmap=dict(MAPPING),
            default_eaf_file=str(default_reference),
            default_eaf_column="EUR",
            default_eaf_colmap=dict(MAPPING),
            study_decision={"eaf_is_maf": None},
            policies=_policies(row_count),
        )

    assert result[resolved_columns["eaf_col"]].to_list() == pytest.approx(
        candidate
    )
    assert qc["external_raw_af_maf_like"] is True
    assert qc["external_maf_reference_decision"] == "eaf"
    assert qc["external_maf_reference_reference_effect_allele_minor_fraction"] == pytest.approx(
        0.95
    )
    assert qc["external_maf_reference_maf_correlation"] == pytest.approx(1.0)
    assert qc[
        "external_maf_reference_mean_absolute_maf_difference"
    ] == pytest.approx(0.0)


def test_independent_default_reference_stops_when_external_maf_screen_is_inconclusive(
    tmp_path,
):
    external = tmp_path / "external.tsv"
    external.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.20\n"
        "1\t200\tT\tC\t0.20\n",
        encoding="utf-8",
    )
    default_reference = tmp_path / "default.tsv"
    default_reference.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.20\n"
        "1\t200\tT\tC\t0.20\n",
        encoding="utf-8",
    )
    oriented = pl.DataFrame({
        "CHR": ["1", "1"],
        "POS": [100, 200],
        "EA": ["A", "C"],
        "OA": ["G", "T"],
        "strand_action": ["forward", "forward"],
        "strand_reference_af": [0.20, 0.20],
    })
    columns, layout = _external_main_contract(tmp_path)

    with patch(
        "postgwas.modules.harmonisation.allele_frequency.harmonise_strand_orientation",
        return_value=(oriented, {"status": "mocked"}, dict(columns)),
    ), pytest.raises(
        AlleleFrequencyError,
        match=r"comparison with the independent default EAF reference was inconclusive",
    ):
        harmonise_allele_frequency(
            "1",
            oriented.drop(["strand_action", "strand_reference_af"]),
            dict(columns),
            layout,
            "\t",
            eaffile=str(external),
            external_eaf_colmap=dict(MAPPING),
            default_eaf_file=str(default_reference),
            default_eaf_column="EUR",
            default_eaf_colmap=dict(MAPPING),
            study_decision={"eaf_is_maf": None},
            policies=_policies(2),
        )


def test_main_calls_reference_check_only_for_study_wide_maf_flag(tmp_path):
    frame = pl.DataFrame({
        "CHR": ["1", "1"], "POS": [100, 200],
        "EA": ["A", "G"], "OA": ["G", "A"], "AF": [0.2, 0.3],
    })
    columns = {
        "gwas_outputname": "study", "output_folder": str(tmp_path),
        "chr_col": "CHR", "pos_col": "POS", "ea_col": "EA", "oa_col": "OA",
        "eaf_col": "AF",
    }
    layout = {
        "frequency_qc_directory": "eaf_qc",
        "missing_eaf": "eaf_qc/{dataset_id}_missing_chr{chromosome}.tsv",
        "out_of_range_eaf": "eaf_qc/{dataset_id}_range_chr{chromosome}.tsv",
        "post_orientation_duplicates": "qc/{dataset_id}_duplicates_chr{chromosome}.tsv",
    }

    with patch(
        "postgwas.modules.harmonisation.allele_frequency.confirm_eaf_with_reference",
        return_value=("eaf", {"comparable_variants": 2}),
    ) as confirm, patch(
        "postgwas.modules.harmonisation.allele_frequency.harmonise_strand_orientation",
        side_effect=lambda chromosome, df, sample_column_dict, **kwargs: (
            df, {"status": "mocked"}, sample_column_dict,
        ),
    ):
        _, qc, _ = harmonise_allele_frequency(
            "1", frame, dict(columns), layout, "\t",
            external_eaf_colmap=MAPPING,
            default_eaf_file="reference.tsv",
            default_eaf_column="EUR",
            default_eaf_colmap=MAPPING,
            study_decision={"eaf_is_maf": True},
        )
        assert qc["decision_source"] == "internal_eaf_reference_validation"
        confirm.assert_called_once()

        harmonise_allele_frequency(
            "1", frame, dict(columns), layout, "\t",
            external_eaf_colmap=MAPPING,
            study_decision={"eaf_is_maf": False},
        )
        confirm.assert_called_once()


def test_main_keeps_opt_in_unmatched_variant_with_internal_eaf(tmp_path):
    frame = pl.DataFrame({
        "CHR": ["1"], "POS": [999], "EA": ["A"], "OA": ["G"],
        "AF": [0.2],
    })
    columns = {
        "gwas_outputname": "study", "output_folder": str(tmp_path),
        "chr_col": "CHR", "pos_col": "POS", "ea_col": "EA",
        "oa_col": "OA", "eaf_col": "AF",
    }
    layout = {
        "frequency_qc_directory": "eaf_qc",
        "missing_eaf": "eaf_qc/{dataset_id}_missing_chr{chromosome}.tsv",
        "out_of_range_eaf": "eaf_qc/{dataset_id}_range_chr{chromosome}.tsv",
        "post_orientation_duplicates": (
            "qc/{dataset_id}_duplicates_chr{chromosome}.tsv"
        ),
    }
    policies = default_policies().with_overrides({
        "strand.unmatched_action": "retain",
    })

    result, qc, _ = harmonise_allele_frequency(
        "1", frame, columns, layout, "\t",
        external_eaf_colmap=MAPPING,
        default_eaf_file=_reference(tmp_path),
        default_eaf_column="EUR",
        default_eaf_colmap=MAPPING,
        study_decision={
            "strand": "forward",
            "eaf_is_maf": False,
            "effect_type": "beta",
        },
        policies=policies,
    )

    assert result.height == 1
    assert result.item(0, "AF") == pytest.approx(0.2)
    assert result.item(0, "strand_reference_af") is None
    assert result.item(0, "strand_action") == "reference_unmatched_retained"
    assert qc["strand_orientation"]["reference_unmatched_retained"] == 1
    assert qc["strand_af_comparable"] == 0


def test_main_raises_when_reference_confirms_unaligned_maf(tmp_path):
    frame = pl.DataFrame({
        "CHR": ["1"], "POS": [100], "EA": ["A"], "OA": ["G"], "AF": [0.2],
    })
    columns = {
        "gwas_outputname": "study", "output_folder": str(tmp_path),
        "chr_col": "CHR", "pos_col": "POS", "ea_col": "EA", "oa_col": "OA",
        "eaf_col": "AF",
    }
    layout = {
        "frequency_qc_directory": "eaf_qc",
        "missing_eaf": "eaf_qc/{dataset_id}_missing_chr{chromosome}.tsv",
        "out_of_range_eaf": "eaf_qc/{dataset_id}_range_chr{chromosome}.tsv",
        "post_orientation_duplicates": "qc/{dataset_id}_duplicates_chr{chromosome}.tsv",
    }

    with patch(
        "postgwas.modules.harmonisation.allele_frequency.confirm_eaf_with_reference",
        return_value=("maf", {"comparable_variants": 1}),
    ), patch(
        "postgwas.modules.harmonisation.allele_frequency.harmonise_strand_orientation",
        side_effect=lambda chromosome, df, sample_column_dict, **kwargs: (
            df, {"status": "mocked"}, sample_column_dict,
        ),
    ), pytest.raises(AlleleFrequencyError, match="confirmed as minor allele frequency"):
        harmonise_allele_frequency(
            "1", frame, columns, layout, "\t",
            external_eaf_colmap=MAPPING,
            default_eaf_file="reference.tsv",
            default_eaf_column="EUR",
            default_eaf_colmap=MAPPING,
            study_decision={"eaf_is_maf": True},
        )


def test_main_raises_when_maf_reference_comparison_is_inconclusive(tmp_path):
    frame = pl.DataFrame({
        "CHR": ["1"], "POS": [100], "EA": ["A"], "OA": ["G"], "AF": [0.2],
    })
    columns = {
        "gwas_outputname": "study", "output_folder": str(tmp_path),
        "chr_col": "CHR", "pos_col": "POS", "ea_col": "EA", "oa_col": "OA",
        "eaf_col": "AF",
    }
    layout = {
        "frequency_qc_directory": "eaf_qc",
        "missing_eaf": "eaf_qc/{dataset_id}_missing_chr{chromosome}.tsv",
        "out_of_range_eaf": "eaf_qc/{dataset_id}_range_chr{chromosome}.tsv",
        "post_orientation_duplicates": "qc/{dataset_id}_duplicates_chr{chromosome}.tsv",
    }
    evidence = {
        "comparable_variants": 1,
        "minimum_overlap": 1000,
        "mean_eaf_error": None,
        "mean_maf_error": None,
        "minimum_error_margin": 0.02,
        "reference_effect_allele_minor_fraction": None,
        "reference_minor_fraction_cutoff": 0.95,
        "decision_reason": "insufficient_overlap",
    }

    with patch(
        "postgwas.modules.harmonisation.allele_frequency.confirm_eaf_with_reference",
        return_value=("inconclusive", evidence),
    ), patch(
        "postgwas.modules.harmonisation.allele_frequency.harmonise_strand_orientation",
        side_effect=lambda chromosome, df, sample_column_dict, **kwargs: (
            df, {"status": "mocked"}, sample_column_dict,
        ),
    ), pytest.raises(
        AlleleFrequencyError,
        match=r"reference comparison was inconclusive.*PostGWAS will not guess",
    ):
        harmonise_allele_frequency(
            "1", frame, columns, layout, "\t",
            external_eaf_colmap=MAPPING,
            default_eaf_file="reference.tsv",
            default_eaf_column="EUR",
            default_eaf_colmap=MAPPING,
            study_decision={"eaf_is_maf": True},
        )
