"""Scientific regression tests for deterministic duplicate handling."""

from pathlib import Path

import polars as pl
import pytest

from postgwas.modules.harmonisation.policies import PolicyError, default_policies
from postgwas.modules.harmonisation.rejects import (
    RejectCollector,
    SOURCE_INPUT_ROW_COLUMN,
)
from postgwas.modules.harmonisation.summary_statistics_io import (
    DUPLICATE_REPORT_COLUMNS,
    read_summary_statistics,
    resolve_duplicate_variants,
)


MAPPING = {
    "gwas_outputname": "study",
    "chr_col": "CHR",
    "pos_col": "POS",
    "chr_pos_col": None,
    "snp_id_col": "SNP",
    "ea_col": "EA",
    "oa_col": "OA",
    "eaf_col": "EAF",
    "beta_or_col": "BETA",
    "se_col": "SE",
    "imp_z_col": "Z",
    "pval_col": "P",
    "ncontrol_col": "NCONTROL",
    "ncase_col": "NCASE",
    "imp_info_col": "INFO",
}


def _collector(frame, tmp_path):
    return RejectCollector(
        source_snapshot=frame,
        logger=None,
        out_path=str(tmp_path / "rejected.tsv"),
        delimiter="\t",
        compress=False,
    )


def test_default_policy_replaces_arbitrary_keep_modes():
    policies = default_policies()

    assert "duplicates.keep" not in policies
    assert policies.get("duplicates.consistency_fields") == [
        "beta", "se", "zscore", "pval", "eaf",
    ]
    assert policies.get("duplicates.selection_order") == [
        "completeness", "sample_size", "info", "input_order",
    ]
    assert policies.get("duplicates.conflicting_action") == "remove_all"


@pytest.mark.parametrize(
    "selection_order",
    [
        ["completeness", "sample_size", "info"],
        ["input_order", "completeness"],
    ],
)
def test_input_order_is_required_as_the_final_tie_breaker(selection_order):
    with pytest.raises(PolicyError, match="input_order"):
        default_policies().with_overrides({
            "duplicates": {"selection_order": selection_order},
        })


def test_consistent_duplicates_use_quality_ranking_and_restore_input_order(tmp_path):
    frame = pl.DataFrame({
        "SNP": [
            "a", "b", "c", "d", "e", "f", "g", "h", "unique",
            "imbalanced", "balanced",
        ],
        "CHR": ["1"] * 11,
        "POS": [100, 100, 101, 101, 102, 102, 103, 103, 104, 105, 105],
        "EA": ["a", "A", "C", "C", "G", "G", "T", "T", "A", "A", "A"],
        "OA": ["g", "G", "T", "T", "A", "A", "C", "C", "C", "G", "G"],
        "BETA": [0.1] * 11,
        "SE": [
            0.02, 0.02, 0.02, 0.02, 0.02, 0.02, None, 0.02, 0.02, 0.02, 0.02,
        ],
        "Z": [5.0] * 11,
        "P": [1e-6] * 11,
        "EAF": [0.2] * 11,
        "NCONTROL": [
            100, 200, 200, 200, 200, 200, 5000, 100, 100, 9900, 1000,
        ],
        "NCASE": [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 1000],
        "INFO": [
            0.9, 0.8, 0.8, 0.9, 0.9, 0.9, 0.99, 0.5, 0.9, 0.9, 0.9,
        ],
    })
    frame = frame.with_row_index(SOURCE_INPUT_ROW_COLUMN, offset=1)
    collector = _collector(frame, tmp_path)

    retained, report, statistics = resolve_duplicate_variants(
        frame,
        dict(MAPPING),
        policies=default_policies(),
        rejects=collector,
    )

    # POS 100: larger N; POS 101: higher INFO; POS 102: earlier input row;
    # POS 103: greater completeness wins before its competitor's larger N.
    # POS 105: balanced N=2,000 has greater effective N than 100+9,900.
    assert retained["SNP"].to_list() == [
        "b", "d", "e", "h", "unique", "balanced",
    ]
    assert report.filter(pl.col("duplicate_action") == "kept")["SNP"].to_list() == [
        "b", "d", "e", "h", "balanced",
    ]
    assert report["duplicate_input_row"].to_list() == [
        1, 2, 3, 4, 5, 6, 7, 8, 10, 11,
    ]
    assert report["duplicate_class"].unique().to_list() == ["consistent"]
    assert statistics == {
        "duplicate_groups": 5,
        "duplicate_rows": 10,
        "consistent_groups": 5,
        "conflicting_groups": 0,
        "consistent_rows_removed": 5,
        "conflicting_rows_removed": 0,
        "rows_removed": 5,
    }
    assert collector.counts() == {"duplicate_variant": 5}
    assert report.columns == [
        column for column in frame.columns if column != SOURCE_INPUT_ROW_COLUMN
    ] + list(DUPLICATE_REPORT_COLUMNS)


def test_conflicting_scientific_values_remove_the_entire_group(tmp_path):
    frame = pl.DataFrame({
        "SNP": ["beta_a", "beta_b", "eaf_a", "eaf_b", "unique"],
        "CHR": ["1"] * 5,
        "POS": [100, 100, 101, 101, 102],
        "EA": ["A", "A", "C", "C", "G"],
        "OA": ["G", "G", "T", "T", "A"],
        "BETA": [0.1, 0.2, 0.3, 0.3, 0.4],
        "SE": [0.02] * 5,
        "Z": [5.0] * 5,
        "P": [1e-6] * 5,
        "EAF": [0.2, 0.2, 0.3, 0.4, 0.5],
        "NCONTROL": [1000] * 5,
        "NCASE": [1000] * 5,
        "INFO": [0.9] * 5,
    })
    frame = frame.with_row_index(SOURCE_INPUT_ROW_COLUMN, offset=1)
    collector = _collector(frame, tmp_path)

    retained, report, statistics = resolve_duplicate_variants(
        frame,
        dict(MAPPING),
        policies=default_policies(),
        rejects=collector,
    )

    assert retained["SNP"].to_list() == ["unique"]
    assert report.filter(pl.col("POS") == 100)[
        "duplicate_conflicting_fields"
    ].unique().to_list() == ["beta"]
    assert report.filter(pl.col("POS") == 101)[
        "duplicate_conflicting_fields"
    ].unique().to_list() == ["eaf"]
    assert report["duplicate_action"].unique().to_list() == ["remove_all"]
    assert statistics["conflicting_groups"] == 2
    assert statistics["conflicting_rows_removed"] == 4
    assert collector.counts() == {"conflicting_duplicate": 4}


def test_duplicate_validation_never_silently_skips_an_unresolved_key():
    frame = pl.DataFrame({
        "CHR": ["1", "1"],
        "POS": [100, 100],
        "EA": ["A", "A"],
        "BETA": [0.1, 0.1],
        "P": [0.1, 0.1],
    })
    mapping = dict(MAPPING)
    mapping["oa_col"] = None

    with pytest.raises(ValueError, match="unresolved fields: oa"):
        resolve_duplicate_variants(frame, mapping, policies=default_policies())


def test_fail_dataset_writes_duplicate_report_before_stopping(tmp_path):
    source = tmp_path / "study.tsv"
    source.write_text(
        "SNP\tCHR\tPOS\tEA\tOA\tBETA\tSE\tZ\tP\tEAF\tNCONTROL\tNCASE\tINFO\n"
        "a\t1\t100\tA\tG\t0.1\t0.02\t5\t1e-6\t0.2\t1000\t1000\t0.9\n"
        "b\t1\t100\tA\tG\t0.2\t0.02\t5\t1e-6\t0.2\t1000\t1000\t0.9\n",
        encoding="utf-8",
    )
    policies = default_policies().with_overrides({
        "duplicates": {"conflicting_action": "fail_dataset"},
    })
    output_layout = {
        "rejected_directory": "rejected",
        "input_reject": "rejected/{dataset_id}_input.tsv",
        "duplicates": "{dataset_id}_duplicates.tsv",
    }

    with pytest.raises(ValueError, match="conflicting duplicate group"):
        read_summary_statistics(
            sumstat_file=str(source),
            output_dir=str(tmp_path),
            sample_column_dict=dict(MAPPING),
            output_layout=output_layout,
            table_delimiter="\t",
            policies=policies,
            input_line_count=2,
        )

    report_path = Path(tmp_path) / "study_duplicates.tsv"
    assert report_path.is_file()
    report = pl.read_csv(report_path, separator="\t")
    assert report.height == 2
    assert report["duplicate_class"].unique().to_list() == ["conflicting"]
    assert report["duplicate_action"].unique().to_list() == ["fail_dataset"]


def test_duplicate_input_row_survives_earlier_input_rejection(tmp_path):
    source = tmp_path / "study.tsv"
    source.write_text(
        "SNP\tCHR\tPOS\tEA\tOA\tBETA\tSE\tZ\tP\tEAF\tNCONTROL\tNCASE\tINFO\n"
        "invalid\t1\t99\tA\tG\t0.1\t0.02\t5\tNA\t0.2\t100\t100\t0.9\n"
        "lower_n\t1\t100\tA\tG\t0.1\t0.02\t5\t1e-6\t0.2\t100\t100\t0.9\n"
        "higher_n\t1\t100\tA\tG\t0.1\t0.02\t5\t1e-6\t0.2\t200\t200\t0.9\n",
        encoding="utf-8",
    )
    output_layout = {
        "rejected_directory": "rejected",
        "input_reject": "rejected/{dataset_id}_input.tsv",
        "duplicates": "{dataset_id}_duplicates.tsv",
    }

    retained, *_ = read_summary_statistics(
        sumstat_file=str(source), output_dir=str(tmp_path),
        sample_column_dict=dict(MAPPING), output_layout=output_layout,
        table_delimiter="\t", policies=default_policies(), input_line_count=3,
    )

    report = pl.read_csv(tmp_path / "study_duplicates.tsv", separator="\t")
    assert retained["SNP"].to_list() == ["higher_n"]
    assert report["duplicate_input_row"].to_list() == [2, 3]
    assert report.filter(pl.col("duplicate_action") == "kept")[
        "duplicate_input_row"
    ].to_list() == [3]
    assert retained[SOURCE_INPUT_ROW_COLUMN].to_list() == [3]
    assert SOURCE_INPUT_ROW_COLUMN not in report.columns


def test_input_reject_restores_values_from_before_coordinate_normalization(tmp_path):
    source = tmp_path / "study.tsv"
    source.write_text(
        "SNP\tCHR\tPOS\tEA\tOA\tBETA\tSE\tZ\tP\tEAF\tNCONTROL\tNCASE\tINFO\n"
        "raw_lower\tchr01\t100.0\ta\tg\t0.1\t0.02\t5\t1e-6\t0.2\t100\t100\t0.9\n"
        "retained\t1\t100\tA\tG\t0.1\t0.02\t5\t1e-6\t0.2\t200\t200\t0.9\n",
        encoding="utf-8",
    )
    output_layout = {
        "rejected_directory": "rejected",
        "input_reject": "rejected/{dataset_id}_input.tsv",
        "duplicates": "{dataset_id}_duplicates.tsv",
    }

    retained, *_ = read_summary_statistics(
        sumstat_file=str(source),
        output_dir=str(tmp_path),
        sample_column_dict=dict(MAPPING),
        output_layout=output_layout,
        table_delimiter="\t",
        policies=default_policies(),
        input_line_count=2,
    )

    rejected = pl.read_csv(
        tmp_path / "rejected" / "study_input.tsv.gz",
        separator="\t",
        infer_schema_length=0,
    )
    assert retained["SNP"].to_list() == ["retained"]
    assert rejected.select(["SNP", "CHR", "POS", "EA", "OA"]).row(0) == (
        "raw_lower", "chr01", "100.0", "a", "g",
    )
    assert rejected["reject_reason"].to_list() == ["duplicate_variant"]
