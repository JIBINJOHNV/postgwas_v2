"""Focused tests for reference confirmation of MAF-like frequency columns."""

from unittest.mock import patch

import polars as pl
import pytest

from postgwas.modules.harmonisation.allele_frequency import (
    AlleleFrequencyError,
    confirm_eaf_with_reference,
    harmonise_allele_frequency,
    merge_external_allele_frequencies,
)
from postgwas.modules.harmonisation.policies import default_policies


MAPPING = {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF", "delimiter": "tab"}
STANDARD = {"chr": "CHR", "pos": "POS", "ea": "EA", "oa": "OA"}


def _reference(tmp_path):
    path = tmp_path / "GRCh37_1000G_freq_chr1.tsv"
    path.write_text(
        "CHROM\tPOS\tREF\tALT\tEUR\n"
        "1\t100\tG\tA\t0.8\n"
        "1\t200\tC\tT\t0.2\n",
        encoding="utf-8",
    )
    return str(path)


def _policies(minimum):
    return default_policies().with_overrides(
        {"eaf": {"maf_reference_min_overlap": minimum}}
    )


def test_reference_check_confirms_unaligned_maf_after_direct_and_swapped_matches(tmp_path):
    study = pl.DataFrame({
        "CHR": ["1", "1"], "POS": [100, 200],
        "EA": ["A", "C"], "OA": ["G", "T"], "AF": [0.2, 0.2],
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
    assert stats["mean_eaf_error"] == pytest.approx(0.6)


def test_reference_check_accepts_frequency_aligned_to_minor_effect_allele(tmp_path):
    study = pl.DataFrame({
        "CHR": ["1", "1"], "POS": [100, 200],
        "EA": ["G", "T"], "OA": ["A", "C"], "AF": [0.2, 0.2],
    })

    decision, stats = confirm_eaf_with_reference(
        study, "AF", _reference(tmp_path), "EUR", MAPPING, STANDARD,
        policies=_policies(2),
    )

    assert decision == "eaf"
    assert stats["reference_effect_allele_minor_fraction"] == pytest.approx(1.0)
    assert stats["mean_eaf_error"] == pytest.approx(0.0)
    assert stats["mean_maf_error"] == pytest.approx(0.0)


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
        policies=_policies(1),
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
            policies=_policies(1),
        )


def test_external_eaf_duplicate_prefers_valid_frequency_over_missing(tmp_path):
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
        policies=default_policies(),
        study_columns_canonical=True,
    )

    assert result[column].to_list() == pytest.approx([0.25])
    assert stats["reference_duplicate_rows_removed"] == 1
    assert stats["reference_duplicate_groups_preferred_usable_value"] == 1


def test_external_eaf_duplicate_rejects_conflicting_frequencies(tmp_path):
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

    with pytest.raises(AlleleFrequencyError, match="conflicting finite values"):
        merge_external_allele_frequencies(
            study,
            str(reference),
            "EUR",
            MAPPING,
            STANDARD,
            policies=default_policies(),
            study_columns_canonical=True,
        )


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
