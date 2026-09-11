"""Scientific and performance contracts for shared allele-oriented joins."""

import polars as pl
import pytest

from postgwas.modules.harmonisation.imputation_quality import (
    ImputationQualityError,
    harmonise_imputation_quality,
)
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.shared.allele_join import (
    allele_oriented_left_join,
)


STUDY_COLUMNS = {"chr": "CHR", "pos": "POS", "ea": "EA", "oa": "OA"}
REFERENCE_COLUMNS = {"chr": "CHROM", "pos": "BP", "ea": "A1", "oa": "A2"}
DUPLICATE_ACTIONS = {
    "duplicate_exact_action": "keep_one",
    "duplicate_non_identical_action": "discard_all",
}


def test_shared_join_preserves_order_and_applies_only_requested_swap_transform():
    study = pl.DataFrame({
        "CHR": ["1", "1", "1"],
        "POS": [20, 10, 30],
        "EA": ["C", "A", "T"],
        "OA": ["T", "G", "C"],
        "ROW": ["swapped", "direct", "unmatched"],
    })
    reference = pl.DataFrame({
        "CHROM": ["chr1", "chr1"],
        "BP": ["10", "20"],
        "A1": [" a ", " t "],
        "A2": [" g ", " c "],
        "AF": [0.2, 0.3],
    })

    result, orientation, stats = allele_oriented_left_join(
        study,
        reference,
        study_columns=STUDY_COLUMNS,
        reference_columns=REFERENCE_COLUMNS,
        policies=default_policies(),
        value_column="AF",
        output_column="EAF",
        **DUPLICATE_ACTIONS,
        swapped_value="one_minus",
        study_columns_canonical=True,
    )

    assert result["ROW"].to_list() == ["swapped", "direct", "unmatched"]
    assert result["EAF"].to_list() == pytest.approx([0.7, 0.2, None])
    assert result[orientation].to_list() == [1, 0, None]
    assert stats["direct_key_matches"] == 1
    assert stats["swapped_key_matches"] == 1
    assert stats["unmatched_rows"] == 1


def test_shared_join_applies_same_chromosome_aliases_to_reference_keys():
    study = pl.DataFrame({
        "CHR": ["1", "X"],
        "POS": [10, 20],
        "EA": ["A", "C"],
        "OA": ["G", "T"],
    })
    reference = pl.DataFrame({
        "CHROM": ["01", "23"],
        "BP": [10, 20],
        "A1": ["A", "C"],
        "A2": ["G", "T"],
        "VALUE": [0.8, 0.9],
    })

    result, orientation, stats = allele_oriented_left_join(
        study,
        reference,
        study_columns=STUDY_COLUMNS,
        reference_columns=REFERENCE_COLUMNS,
        policies=default_policies(),
        value_column="VALUE",
        output_column="VALUE",
        **DUPLICATE_ACTIONS,
    )

    assert result["VALUE"].to_list() == pytest.approx([0.8, 0.9])
    assert result[orientation].to_list() == [0, 0]
    assert stats["direct_value_matches"] == 2
    assert stats["unmatched_rows"] == 0


def test_reference_aliases_are_canonicalized_before_duplicate_resolution():
    study = pl.DataFrame({
        "CHR": ["X"], "POS": [10], "EA": ["A"], "OA": ["G"],
    })
    reference = pl.DataFrame({
        "CHROM": ["23", "X"],
        "BP": [10, 10],
        "A1": ["A", "A"],
        "A2": ["G", "G"],
        "VALUE": [0.8, 0.9],
    })

    result, orientation, stats = allele_oriented_left_join(
        study,
        reference,
        study_columns=STUDY_COLUMNS,
        reference_columns=REFERENCE_COLUMNS,
        policies=default_policies(),
        value_column="VALUE",
        output_column="VALUE",
        **DUPLICATE_ACTIONS,
    )

    assert result["VALUE"].to_list() == [None]
    assert result[orientation].to_list() == [None]
    assert stats["reference_non_identical_duplicate_groups"] == 1
    assert stats["reference_duplicate_rows_removed"] == 2


def test_shared_join_prefers_direct_and_can_fall_back_from_an_empty_value():
    study = pl.DataFrame({
        "CHR": ["1"], "POS": [10], "EA": ["A"], "OA": ["G"],
    })
    reference = pl.DataFrame({
        "CHROM": ["1", "1"],
        "BP": [10, 10],
        "A1": ["A", "G"],
        "A2": ["G", "A"],
        "VALUE": [None, 0.8],
    })

    direct, direct_orientation, _ = allele_oriented_left_join(
        study,
        reference,
        study_columns=STUDY_COLUMNS,
        reference_columns=REFERENCE_COLUMNS,
        policies=default_policies(),
        value_column="VALUE",
        output_column="VALUE",
        **DUPLICATE_ACTIONS,
        swapped_value="same",
        prefer_non_null_value=False,
    )
    fallback, fallback_orientation, _ = allele_oriented_left_join(
        study,
        reference,
        study_columns=STUDY_COLUMNS,
        reference_columns=REFERENCE_COLUMNS,
        policies=default_policies(),
        value_column="VALUE",
        output_column="VALUE",
        **DUPLICATE_ACTIONS,
        swapped_value="same",
        prefer_non_null_value=True,
    )

    assert direct["VALUE"].to_list() == [None]
    assert direct[direct_orientation].to_list() == [0]
    assert fallback["VALUE"].to_list() == [0.8]
    assert fallback[fallback_orientation].to_list() == [1]


def test_shared_join_discards_missing_and_finite_non_identical_duplicates():
    study = pl.DataFrame({
        "CHR": ["1"], "POS": [10], "EA": ["A"], "OA": ["G"],
    })
    reference = pl.DataFrame({
        "CHROM": ["1", "1", "1"],
        "BP": [10, 10, 10],
        "A1": ["A", "A", "A"],
        "A2": ["G", "G", "G"],
        "VALUE": [None, float("nan"), 0.9],
    })

    result, _, stats = allele_oriented_left_join(
        study,
        reference,
        study_columns=STUDY_COLUMNS,
        reference_columns=REFERENCE_COLUMNS,
        policies=default_policies(),
        value_column="VALUE",
        output_column="VALUE",
        **DUPLICATE_ACTIONS,
    )

    assert result["VALUE"].to_list() == [None]
    assert stats["reference_duplicate_rows_removed"] == 3
    assert stats["reference_duplicate_groups"] == 1
    assert stats["reference_exact_duplicate_groups"] == 0
    assert stats["reference_non_identical_duplicate_groups"] == 1
    assert stats["reference_duplicate_groups_discarded"] == 1


def test_shared_join_discards_conflicting_finite_reference_duplicates():
    study = pl.DataFrame({
        "CHR": ["1"], "POS": [10], "EA": ["A"], "OA": ["G"],
    })
    reference = pl.DataFrame({
        "CHROM": ["1", "1"], "BP": [10, 10],
        "A1": ["A", "A"], "A2": ["G", "G"],
        "VALUE": [0.7, 0.9],
    })

    result, orientation, stats = allele_oriented_left_join(
        study,
        reference,
        study_columns=STUDY_COLUMNS,
        reference_columns=REFERENCE_COLUMNS,
        policies=default_policies(),
        value_column="VALUE",
        output_column="VALUE",
        **DUPLICATE_ACTIONS,
    )

    assert result["VALUE"].to_list() == [None]
    assert result[orientation].to_list() == [None]
    assert stats["reference_non_identical_duplicate_groups"] == 1
    assert stats["reference_duplicate_rows_removed"] == 2


def test_shared_join_can_fail_on_non_identical_reference_duplicates():
    study = pl.DataFrame({
        "CHR": ["1"], "POS": [10], "EA": ["A"], "OA": ["G"],
    })
    reference = pl.DataFrame({
        "CHROM": ["1", "1"], "BP": [10, 10],
        "A1": ["A", "A"], "A2": ["G", "G"],
        "VALUE": [0.7, 0.9],
    })

    with pytest.raises(
        RuntimeError,
        match=r"non-identical duplicate values.*CHR='1'.*POS=10.*\[0.7, 0.9\]",
    ):
        allele_oriented_left_join(
            study,
            reference,
            study_columns=STUDY_COLUMNS,
            reference_columns=REFERENCE_COLUMNS,
            policies=default_policies(),
            value_column="VALUE",
            output_column="VALUE",
            duplicate_exact_action="keep_one",
            duplicate_non_identical_action="fail",
        )


def test_shared_join_collapses_equivalent_finite_reference_duplicates():
    study = pl.DataFrame({
        "CHR": ["1"], "POS": [10], "EA": ["A"], "OA": ["G"],
    })
    reference = pl.DataFrame({
        "CHROM": ["1", "1"], "BP": [10, 10],
        "A1": ["A", "A"], "A2": ["G", "G"],
        "VALUE": [0.9, 0.9],
    })

    result, _, stats = allele_oriented_left_join(
        study,
        reference,
        study_columns=STUDY_COLUMNS,
        reference_columns=REFERENCE_COLUMNS,
        policies=default_policies(),
        value_column="VALUE",
        output_column="VALUE",
        **DUPLICATE_ACTIONS,
    )

    assert result["VALUE"].to_list() == pytest.approx([0.9])
    assert stats["reference_duplicate_rows_removed"] == 1
    assert stats["reference_exact_duplicate_groups"] == 1
    assert stats["reference_non_identical_duplicate_groups"] == 0


def test_shared_join_can_discard_all_exact_reference_duplicates():
    study = pl.DataFrame({
        "CHR": ["1"], "POS": [10], "EA": ["A"], "OA": ["G"],
    })
    reference = pl.DataFrame({
        "CHROM": ["1", "1"], "BP": [10, 10],
        "A1": ["A", "A"], "A2": ["G", "G"],
        "VALUE": [0.9, 0.9],
    })

    result, orientation, stats = allele_oriented_left_join(
        study,
        reference,
        study_columns=STUDY_COLUMNS,
        reference_columns=REFERENCE_COLUMNS,
        policies=default_policies(),
        value_column="VALUE",
        output_column="VALUE",
        duplicate_exact_action="discard_all",
        duplicate_non_identical_action="discard_all",
    )

    assert result["VALUE"].to_list() == [None]
    assert result[orientation].to_list() == [None]
    assert stats["reference_duplicate_groups_discarded"] == 1
    assert stats["reference_duplicate_rows_removed"] == 2


def _info_inputs(tmp_path, duplicate=False):
    rows = [
        "1\t10\tA\tG\t0.9",
        "1\t20\tT\tC\t0.8",
    ]
    if duplicate:
        rows.insert(1, "1\t10\tA\tG\t0.7")
    reference = tmp_path / "info.tsv"
    reference.write_text(
        "CHROM\tPOS\tA1\tA2\tINFO\n" + "\n".join(rows) + "\n",
        encoding="utf-8",
    )
    study = pl.DataFrame({
        "CHR": ["1", "1", "1"],
        "POS": [10, 20, 30],
        "EA": ["A", "C", "A"],
        "OA": ["G", "T", "C"],
    })
    columns = {
        "chr_col": "CHR", "pos_col": "POS", "ea_col": "EA", "oa_col": "OA",
        "imp_info_col": None,
    }
    mapping = {
        "chr": "CHROM", "pos": "POS", "a1": "A1", "a2": "A2",
        "delimiter": "tab",
    }
    return reference, study, columns, mapping


def _info_policies(**duplicate_overrides):
    return default_policies().with_overrides({
        "info": {"source": "reference", "score_type": "standard_info"},
        "external_reference": duplicate_overrides,
    })


def test_external_info_uses_one_join_for_direct_swapped_and_unmatched(
    tmp_path, monkeypatch,
):
    reference, study, columns, mapping = _info_inputs(tmp_path)
    original_join = pl.DataFrame.join
    joins = 0

    def counting_join(self, *args, **kwargs):
        nonlocal joins
        joins += 1
        return original_join(self, *args, **kwargs)

    monkeypatch.setattr(pl.DataFrame, "join", counting_join)
    result, qc, resolved = harmonise_imputation_quality(
        "1",
        study,
        columns,
        info_file=str(reference),
        info_column="INFO",
        external_info_colmap=mapping,
        policies=_info_policies(),
    )

    assert joins == 1
    assert result[resolved["imp_info_col"]].to_list() == pytest.approx(
        [0.9, 0.8, None]
    )
    assert qc["matched_direct"] == 1
    assert qc["matched_flipped"] == 1
    assert qc["missing_info"] == 1


def test_external_info_matches_policy_mapped_reference_chromosome(tmp_path):
    reference = tmp_path / "info_x.tsv"
    reference.write_text(
        "CHROM\tPOS\tA1\tA2\tINFO\n23\t10\tA\tG\t0.9\n",
        encoding="utf-8",
    )
    study = pl.DataFrame({
        "CHR": ["X"], "POS": [10], "EA": ["A"], "OA": ["G"],
    })
    columns = {
        "chr_col": "CHR", "pos_col": "POS", "ea_col": "EA", "oa_col": "OA",
        "imp_info_col": None,
    }
    mapping = {
        "chr": "CHROM", "pos": "POS", "a1": "A1", "a2": "A2",
        "delimiter": "tab",
    }

    result, qc, resolved = harmonise_imputation_quality(
        "X",
        study,
        columns,
        info_file=str(reference),
        info_column="INFO",
        external_info_colmap=mapping,
        policies=_info_policies(),
    )

    assert result[resolved["imp_info_col"]].to_list() == pytest.approx([0.9])
    assert qc["matched_direct"] == 1
    assert qc["missing_info"] == 0


def test_external_info_non_identical_duplicates_are_all_discarded(tmp_path):
    reference, study, columns, mapping = _info_inputs(tmp_path, duplicate=True)

    result, qc, resolved = harmonise_imputation_quality(
        "1",
        study,
        columns,
        info_file=str(reference),
        info_column="INFO",
        external_info_colmap=mapping,
        policies=_info_policies(),
    )

    assert result[resolved["imp_info_col"]].to_list() == pytest.approx(
        [None, 0.8, None]
    )
    assert qc["matched_direct"] == 0
    assert qc["matched_flipped"] == 1
    assert qc["missing_info"] == 2


def test_external_info_non_identical_duplicate_policy_can_fail(tmp_path):
    reference, study, columns, mapping = _info_inputs(tmp_path, duplicate=True)

    with pytest.raises(
        ImputationQualityError,
        match="external_reference.non_identical_duplicate_action",
    ):
        harmonise_imputation_quality(
            "1",
            study,
            columns,
            info_file=str(reference),
            info_column="INFO",
            external_info_colmap=mapping,
            policies=_info_policies(non_identical_duplicate_action="fail"),
        )


def test_external_info_duplicate_missing_and_finite_values_discards_both(tmp_path):
    reference, study, columns, mapping = _info_inputs(tmp_path)
    reference.write_text(
        "CHROM\tPOS\tA1\tA2\tINFO\n"
        "1\t10\tA\tG\tNA\n"
        "1\t10\tA\tG\t0.9\n"
        "1\t20\tT\tC\t0.8\n",
        encoding="utf-8",
    )

    result, qc, resolved = harmonise_imputation_quality(
        "1",
        study,
        columns,
        info_file=str(reference),
        info_column="INFO",
        external_info_colmap=mapping,
        policies=_info_policies(),
    )

    assert result[resolved["imp_info_col"]].to_list() == pytest.approx(
        [None, 0.8, None]
    )
    assert qc["missing_info"] == 2
    assert qc["matched_direct"] == 0
    assert qc["matched_flipped"] == 1


def test_external_info_exact_duplicates_keep_one(tmp_path):
    reference, study, columns, mapping = _info_inputs(tmp_path)
    reference.write_text(
        "CHROM\tPOS\tA1\tA2\tINFO\n"
        "1\t10\tA\tG\t0.9\n"
        "1\t10\tA\tG\t0.9\n"
        "1\t20\tT\tC\t0.8\n",
        encoding="utf-8",
    )

    result, qc, resolved = harmonise_imputation_quality(
        "1",
        study,
        columns,
        info_file=str(reference),
        info_column="INFO",
        external_info_colmap=mapping,
        policies=_info_policies(),
    )

    assert result[resolved["imp_info_col"]].to_list() == pytest.approx(
        [0.9, 0.8, None]
    )
    assert qc["matched_direct"] == 1
    assert qc["matched_flipped"] == 1


def test_eaf_and_info_share_one_external_reference_duplicate_policy():
    policies = default_policies()

    assert policies.get("external_reference.exact_duplicate_action") == "keep_one"
    assert (
        policies.get("external_reference.non_identical_duplicate_action")
        == "discard_all"
    )
    with pytest.raises(KeyError):
        policies.get("eaf.deduplicate_reference")
    with pytest.raises(KeyError):
        policies.get("info.deduplicate_reference")
