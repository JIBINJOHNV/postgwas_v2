"""Regression tests for behavior-preserving Harmonisation optimizations."""

import csv
import json
from pathlib import Path

import polars as pl
import pytest

from postgwas.config import load_configuration
from postgwas.modules.harmonisation.chromosome_partition import (
    write_chromosome_partitions,
)
from postgwas.modules.harmonisation.cleanup import finalise_harmonisation_outputs
from postgwas.modules.harmonisation.coordinates import (
    harmonise_coordinates_and_alleles,
)
from postgwas.modules.harmonisation.effect_type import harmonise_effect_estimates
from postgwas.modules.harmonisation.effect_type import detect_effect_type
from postgwas.modules.harmonisation.effect_from_z import (
    derive_effect_and_standard_error_from_z,
)
from postgwas.modules.harmonisation.gwas2vcf_export import (
    summarise_gwas2vcf_columns,
)
from postgwas.modules.harmonisation.imputation_quality import (
    harmonise_imputation_quality,
)
from postgwas.modules.harmonisation.input_validation import validate_content
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.p_values import detect_p_value_type
from postgwas.modules.harmonisation.qc_results import qc_results_to_dataframe
from postgwas.modules.harmonisation.rejects import (
    RejectCollector,
    SOURCE_INPUT_ROW_COLUMN,
)
from postgwas.modules.harmonisation.sample_size import harmonise_sample_sizes
from postgwas.modules.harmonisation.service import _save_qc_results
from postgwas.modules.harmonisation.study_properties import resolve_study_properties
from postgwas.modules.harmonisation.summary_statistics_io import (
    check_file_truncation,
    inspect_summary_statistics_file,
    load_summary_statistics_table,
    normalise_summary_statistics_values,
    read_summary_statistics,
)
from postgwas.modules.harmonisation.shared.variant_columns import (
    CANONICAL_VARIANT_COLUMNS_KEY,
    canonical_variant_schema,
    has_canonical_variant_columns,
    mark_canonical_variant_columns,
    study_string_schema,
)


def _cleanup_layout():
    return {
        "qc_directory": "qc_summary",
        "adapter_merged_mapping": "qc_summary/{dataset_id}_column_mapping.json",
        "gwas2vcf_summary": "qc_summary/{dataset_id}_summary.tsv",
        "adapter_mapping": "{dataset_id}_chr{chromosome}.dict",
        "adapter_summary": "{dataset_id}_chr{chromosome}_summary.tsv",
        "adapter_input_archive": "qc_summary/{dataset_id}_inputs",
        "adapter_input": "{dataset_id}_chr{chromosome}_vcf_input.tsv",
        "chromosome_table": "{dataset_id}_chr{chromosome}.tsv",
        "chromosome_source_snapshot": (
            "rejected/{dataset_id}_chr{chromosome}_source.parquet"
        ),
    }


def test_input_inspection_counts_data_during_the_truncation_pass(tmp_path):
    source = tmp_path / "study.tsv"
    source.write_text(
        "##metadata\nCHR\tPOS\n1\t10\n1\t20\n",
        encoding="utf-8",
    )

    warnings, data_rows = inspect_summary_statistics_file(
        source, policies=default_policies(),
    )

    assert warnings == []
    assert data_rows == 2
    assert check_file_truncation(source, policies=default_policies()) == []


def test_study_property_detection_uses_one_aggregation(monkeypatch):
    frame = pl.DataFrame({
        "EFFECT": [-0.1] * 2 + [0.05] * 198,
        "P": [0.5] * 200,
        "EAF": [0.2] * 200,
    })
    policies = default_policies()
    expected_effect, effect_evidence = detect_effect_type(
        frame, "EFFECT", policies,
    )
    expected_pvalue, pvalue_evidence = detect_p_value_type(frame, "P", policies)

    original_select = pl.DataFrame.select
    calls = 0

    def counting_select(self, *expressions, **named_expressions):
        nonlocal calls
        calls += 1
        return original_select(self, *expressions, **named_expressions)

    monkeypatch.setattr(pl.DataFrame, "select", counting_select)
    decisions = resolve_study_properties(
        frame,
        {"beta_or_col": "EFFECT", "pval_col": "P", "eaf_col": "EAF"},
        policies=policies,
    )

    assert calls == 1
    assert decisions["effect_type"] == expected_effect == "beta"
    assert decisions["pvalue_type"] == expected_pvalue == "raw"
    assert decisions["eaf_is_maf"] is True
    assert (
        decisions["effect_evidence"]["detector_stats"]["median"]
        == effect_evidence["median"]
    )
    assert (
        decisions["pvalue_evidence"]["detector_stats"]["median"]
        == pvalue_evidence["median"]
    )


def test_chromosome_partition_preserves_rows_and_columns(tmp_path):
    frame = pl.DataFrame({
        "CHR": ["2", "1", "2", "X"],
        "POS": [20, 10, 21, 30],
        "EA": ["A", "C", "G", "T"],
    })
    paths = write_chromosome_partitions(
        frame,
        {"output_folder": str(tmp_path), "chr_col": "CHR", "gwas_outputname": "study"},
        {"chromosome_table": "{dataset_id}_chr{chromosome}.tsv"},
        "\t",
    )

    assert set(paths) == {"1", "2", "X"}
    chromosome_two = pl.read_csv(paths["2"], separator="\t")
    assert chromosome_two.columns == frame.columns
    assert chromosome_two["POS"].to_list() == [20, 21]
    assert sum(pl.read_csv(path, separator="\t").height for path in paths.values()) == frame.height


def test_chromosome_partition_writes_original_source_rows_by_stable_id(tmp_path):
    source = pl.DataFrame({
        "CHR": ["chr01", "2"],
        "POS": ["100.0", "200"],
        "EA": ["a", "C"],
        "OA": ["g", "T"],
        "EAF": [0.2, 0.3],
    }).with_row_index(SOURCE_INPUT_ROW_COLUMN, offset=1)
    working = source.with_columns([
        pl.Series("CHR", ["1", "2"]),
        pl.Series("POS", [100, 200]),
        pl.col("EA").str.to_uppercase(),
        pl.col("OA").str.to_uppercase(),
    ])
    layout = {
        "chromosome_table": "{dataset_id}_chr{chromosome}.tsv",
        "chromosome_source_snapshot": (
            "rejected/{dataset_id}_chr{chromosome}_source.parquet"
        ),
    }

    write_chromosome_partitions(
        working,
        {
            "output_folder": str(tmp_path),
            "chr_col": "CHR",
            "gwas_outputname": "study",
        },
        layout,
        "\t",
        source_snapshot=source,
    )

    chromosome_one = pl.read_parquet(
        tmp_path / "rejected" / "study_chr1_source.parquet"
    )
    assert chromosome_one.select(["CHR", "POS", "EA", "OA", "EAF"]).row(0) == (
        "chr01", "100.0", "a", "g", 0.2,
    )


def test_successful_cleanup_removes_internal_source_snapshots(tmp_path):
    chromosome_table = tmp_path / "study_chr1.tsv"
    source_snapshot = tmp_path / "rejected" / "study_chr1_source.parquet"
    adapter_mapping = tmp_path / "study_chr1.dict"
    adapter_summary = tmp_path / "study_chr1_summary.tsv"
    source_snapshot.parent.mkdir()
    chromosome_table.write_text("CHR\tPOS\n1\t100\n", encoding="utf-8")
    pl.DataFrame({"CHR": ["1"], "POS": [100]}).write_parquet(source_snapshot)
    adapter_mapping.write_text('{"chr_col": 0}\n', encoding="utf-8")
    adapter_summary.write_text(
        "chromosome\tkey\tnum_rows\n1\tsnp_id_col\t1\n",
        encoding="utf-8",
    )
    layout = _cleanup_layout()

    result = finalise_harmonisation_outputs(
        output_dir=str(tmp_path),
        gwas_outputname="study",
        output_layout=layout,
        threads=1,
        compression_executable=None,
    )

    assert not chromosome_table.exists()
    assert not source_snapshot.exists()
    assert set(result["removed"]) == {
        str(chromosome_table),
        str(source_snapshot),
        str(adapter_mapping),
        str(adapter_summary),
    }
    assert json.loads(Path(result["dict"]).read_text(encoding="utf-8")) == {
        "chr_col": 0,
    }
    assert "snp_id_col" in Path(result["summary"]).read_text(encoding="utf-8")


def test_cleanup_rejects_inconsistent_chromosome_column_mappings(tmp_path):
    (tmp_path / "study_chr1.dict").write_text(
        '{"chr_col": 0}\n', encoding="utf-8",
    )
    (tmp_path / "study_chr2.dict").write_text(
        '{"chr_col": 1}\n', encoding="utf-8",
    )
    for chromosome in ("1", "2"):
        (tmp_path / ("study_chr%s_summary.tsv" % chromosome)).write_text(
            "chromosome\tkey\tnum_rows\n%s\tsnp_id_col\t1\n" % chromosome,
            encoding="utf-8",
        )

    with pytest.raises(RuntimeError, match="column mappings are inconsistent"):
        finalise_harmonisation_outputs(
            output_dir=str(tmp_path),
            gwas_outputname="study",
            output_layout=_cleanup_layout(),
            threads=1,
            compression_executable=None,
        )

    assert (tmp_path / "study_chr1.dict").is_file()
    assert (tmp_path / "study_chr2.dict").is_file()


def test_cleanup_rejects_stale_merged_side_outputs(tmp_path):
    chromosome_table = tmp_path / "study_chr1.tsv"
    chromosome_table.write_text("CHR\tPOS\n1\t100\n", encoding="utf-8")
    qc_dir = tmp_path / "qc_summary"
    qc_dir.mkdir()
    stale_mapping = qc_dir / "study_column_mapping.json"
    stale_summary = qc_dir / "study_summary.tsv"
    stale_mapping.write_text("previous mapping\n", encoding="utf-8")
    stale_summary.write_text("previous summary\n", encoding="utf-8")
    layout = _cleanup_layout()

    with pytest.raises(RuntimeError, match="may belong to an earlier run"):
        finalise_harmonisation_outputs(
            output_dir=str(tmp_path),
            gwas_outputname="study",
            output_layout=layout,
            threads=1,
            compression_executable=None,
        )

    assert stale_mapping.read_text(encoding="utf-8") == "previous mapping\n"
    assert stale_summary.read_text(encoding="utf-8") == "previous summary\n"
    assert chromosome_table.exists()


def test_first_read_keeps_variant_identity_columns_as_strings(tmp_path):
    source = tmp_path / "study.tsv"
    source.write_text(
        "CHR\tPOS\tSNP\tEA\tOA\tBETA\n"
        "1\t100\t000123\tA\tG\t0.1\n"
        "Y\t200\t9007199254740993\tC\tT\t0.2\n",
        encoding="utf-8",
    )
    mapping = {
        "chr_col": "CHR", "pos_col": "POS", "snp_id_col": "SNP",
        "ea_col": "EA", "oa_col": "OA", "beta_or_col": "BETA",
    }
    policies = default_policies().with_overrides({"input.schema_inference_rows": 1})

    frame, _ = load_summary_statistics_table(
        str(source), mapping, policies=policies,
    )

    assert frame.schema["CHR"] == pl.String
    assert frame.schema["SNP"] == pl.String
    assert frame.schema["EA"] == pl.String
    assert frame.schema["OA"] == pl.String
    assert frame["CHR"].to_list() == ["1", "Y"]
    assert frame["SNP"].to_list() == ["000123", "9007199254740993"]


def test_first_read_keeps_combined_coordinate_as_string(tmp_path):
    source = tmp_path / "combined.tsv"
    source.write_text(
        "VARIANT\tSNP\tEA\tOA\n1:100000\t001\tA\tG\n",
        encoding="utf-8",
    )
    mapping = {
        "chr_pos_col": "VARIANT", "snp_id_col": "SNP",
        "ea_col": "EA", "oa_col": "OA",
    }

    frame, _ = load_summary_statistics_table(
        str(source), mapping, policies=default_policies(),
    )

    assert frame.schema["VARIANT"] == pl.String
    assert frame.schema["SNP"] == pl.String
    assert frame.row(0)[:2] == ("1:100000", "001")


def test_study_variant_columns_are_canonicalised_once_and_partition_types_restore(tmp_path):
    source = tmp_path / "study.tsv"
    source.write_text(
        "CHR\tPOS\tEA\tOA\tSNP\tBETA\n"
        "chr01\t100.0\ta\tg\t000123\t0.25\n",
        encoding="utf-8",
    )
    (
        frame,
        _source_snapshot,
        _input_rows,
        _read_rows,
        mapping,
        _removed_coordinates,
        _removed_alleles,
        _removed_duplicates,
        _removed_missing,
    ) = read_summary_statistics(
        sumstat_file=str(source),
        output_dir=str(tmp_path),
        sample_column_dict={
            "chr_col": "CHR",
            "pos_col": "POS",
            "ea_col": "EA",
            "oa_col": "OA",
            "snp_id_col": "SNP",
            "beta_or_col": "BETA",
            "gwas_outputname": "study",
        },
        output_layout={
            "rejected_directory": "rejected",
            "input_reject": "rejected/{dataset_id}_input.tsv",
            "duplicates": "{dataset_id}_duplicates.tsv",
        },
        table_delimiter="\t",
        policies=default_policies(),
        input_line_count=1,
    )

    assert mapping[CANONICAL_VARIANT_COLUMNS_KEY] is True
    assert has_canonical_variant_columns(frame, mapping)
    assert frame.select(["CHR", "POS", "EA", "OA"]).row(0) == (
        "1", 100, "A", "G",
    )
    assert frame.schema["BETA"] == pl.Float64

    mapping.update({"output_folder": str(tmp_path), "gwas_outputname": "study"})
    paths = write_chromosome_partitions(
        frame, mapping, {"chromosome_table": "{dataset_id}_chr{chromosome}.tsv"}, "\t",
    )
    restored = pl.read_csv(
        paths["1"], separator="\t", schema_overrides=canonical_variant_schema(mapping),
    )
    assert has_canonical_variant_columns(restored, mapping)
    assert restored.schema["CHR"] == pl.String
    assert restored.schema["POS"] == pl.Int64
    assert restored.schema["SNP"] == pl.String
    assert restored.item(0, "SNP") == "000123"


def test_content_validation_skips_repeated_scan_for_canonical_variants(monkeypatch):
    frame, mapping = harmonise_coordinates_and_alleles(
        chromosome="All_Chrs",
        df=pl.DataFrame({
            "CHR": ["1", "2"],
            "POS": [100, 200],
            "EA": ["A", "C"],
            "OA": ["G", "T"],
        }),
        sample_column_dict={
            "chr_col": "CHR", "pos_col": "POS",
            "ea_col": "EA", "oa_col": "OA",
        },
        policies=default_policies(),
    )
    mark_canonical_variant_columns(mapping)
    original_select = pl.DataFrame.select
    calls = 0

    def counting_select(self, *expressions, **named_expressions):
        nonlocal calls
        calls += 1
        return original_select(self, *expressions, **named_expressions)

    monkeypatch.setattr(pl.DataFrame, "select", counting_select)
    valid, problems = validate_content(
        mapping, frame, policies=default_policies(),
    )

    assert valid
    assert problems == []
    assert calls == 1  # The configured-column null-count aggregation only.


def test_content_validation_keeps_coordinate_scan_for_direct_callers():
    frame, mapping = harmonise_coordinates_and_alleles(
        chromosome="All_Chrs",
        df=pl.DataFrame({
            "CHR": ["1"],
            "POS": [100],
            "EA": ["not-an-allele"],
            "OA": ["G"],
        }),
        sample_column_dict={
            "chr_col": "CHR", "pos_col": "POS",
            "ea_col": "EA", "oa_col": "OA",
        },
        policies=default_policies(),
    )

    valid, problems = validate_content(
        mapping, frame, policies=default_policies(),
    )

    assert CANONICAL_VARIANT_COLUMNS_KEY not in mapping
    assert not valid
    assert any(problem.category == "no_variants" for problem in problems)


def test_partition_schema_keeps_snp_and_combined_coordinate_as_strings():
    mapping = {
        "chr_col": "CHR", "pos_col": "POS", "chr_pos_col": "CHR_POS",
        "snp_id_col": "SNP", "ea_col": "EA", "oa_col": "OA",
    }

    assert study_string_schema(mapping) == {
        "CHR": pl.String,
        "CHR_POS": pl.String,
        "SNP": pl.String,
        "EA": pl.String,
        "OA": pl.String,
    }
    assert canonical_variant_schema(mapping) == {
        "CHR": pl.String,
        "CHR_POS": pl.String,
        "SNP": pl.String,
        "EA": pl.String,
        "OA": pl.String,
        "POS": pl.Int64,
    }


def test_only_configured_numeric_columns_are_retyped():
    frame = normalise_summary_statistics_values(
        pl.DataFrame({"BETA": ["0.25"], "EXTRA": [123], "SNP": [9007199254740993]}),
        "CHR",
        policies=default_policies(),
        numeric_columns=["BETA"],
    )

    assert frame.schema == {
        "BETA": pl.Float64,
        "EXTRA": pl.Int64,
        "SNP": pl.Int64,
    }
    assert frame.item(0, "SNP") == 9007199254740993


def test_numeric_conversion_keeps_valid_cells_and_nulls_only_failures():
    frame = normalise_summary_statistics_values(
        pl.DataFrame({
            "P": ["0.1", ".", "bad"],
            "BETA": ["0.2", "bad", "0.4"],
            "N": ["1,000", "bad", "100.9"],
            "INFO": ["0.9", "bad", None],
            "SNP": ["001", "002", "003"],
        }),
        "CHR",
        policies=default_policies(),
        numeric_columns=["P", "BETA", "N", "INFO"],
        sample_count_columns=["N"],
        compound_numeric_columns=["INFO"],
    )

    assert frame.schema == {
        "P": pl.Float64,
        "BETA": pl.Float64,
        "N": pl.Int64,
        "INFO": pl.Float64,
        "SNP": pl.String,
    }
    assert frame["P"].to_list() == [0.1, None, None]
    assert frame["BETA"].to_list() == [0.2, None, 0.4]
    assert frame["N"].to_list() == [1000, None, 100]
    assert frame["INFO"].to_list() == [0.9, None, None]
    assert frame["SNP"].to_list() == ["001", "002", "003"]


def test_compound_info_remains_text_for_its_field_specific_parser():
    frame = normalise_summary_statistics_values(
        pl.DataFrame({"INFO": ["0.9,0.8", "bad"]}),
        "CHR",
        policies=default_policies(),
        numeric_columns=["INFO"],
        compound_numeric_columns=["INFO"],
    )

    assert frame.schema["INFO"] == pl.String
    assert frame["INFO"].to_list() == ["0.9,0.8", "bad"]


def test_dot_is_a_default_missing_token_and_required_p_is_rejected(tmp_path):
    source = tmp_path / "study.tsv"
    source.write_text(
        "SNP\tCHR\tPOS\tEA\tOA\tBETA\tP\n"
        "missing_p\t1\t100\tA\tG\t0.1\t.\n"
        "retained\t1\t101\tC\tT\t0.2\t1e-6\n",
        encoding="utf-8",
    )
    policies = default_policies()
    assert "." in policies.get("input.null_values")

    result = read_summary_statistics(
        sumstat_file=str(source),
        output_dir=str(tmp_path),
        sample_column_dict={
            "gwas_outputname": "study",
            "chr_col": "CHR",
            "pos_col": "POS",
            "snp_id_col": "SNP",
            "ea_col": "EA",
            "oa_col": "OA",
            "beta_or_col": "BETA",
            "pval_col": "P",
        },
        output_layout={
            "rejected_directory": "rejected",
            "input_reject": "rejected/{dataset_id}_input.tsv",
            "duplicates": "{dataset_id}_duplicates.tsv",
        },
        table_delimiter="\t",
        policies=policies,
        input_line_count=2,
    )
    retained = result[0]
    removed_missing = result[-1]
    rejected = pl.read_csv(
        tmp_path / "rejected" / "study_input.tsv.gz",
        separator="\t",
        infer_schema_length=0,
    )

    assert retained["SNP"].to_list() == ["retained"]
    assert removed_missing == 1
    assert rejected["SNP"].to_list() == ["missing_p"]
    assert rejected["reject_reason"].to_list() == ["missing_required_columns"]


@pytest.mark.parametrize(
    "combined",
    [
        "1:100000",
        "1_100000",
        "1:100000:C:T",
        "1:100000.0",
        "1:100000-100001",
    ],
)
def test_combined_coordinate_default_preserves_first_position(combined):
    frame, mapping = harmonise_coordinates_and_alleles(
        chromosome="All_Chrs",
        df=pl.DataFrame({"VARIANT": [combined]}),
        sample_column_dict={
            "chr_col": None,
            "pos_col": None,
            "chr_pos_col": "VARIANT",
        },
        policies=default_policies(),
    )

    assert default_policies().get("position.extraction") == "leading_digits"
    assert mapping["chr_col"] == "Chr"
    assert mapping["pos_col"] == "Pos"
    assert frame["Chr"].to_list() == ["1"]
    assert frame["Pos"].to_list() == [100000]


def test_separate_fractional_position_keeps_existing_truncation_behavior():
    frame, _mapping = harmonise_coordinates_and_alleles(
        chromosome="All_Chrs",
        df=pl.DataFrame({
            "CHR": ["1"], "POS": ["100000.9"], "EA": ["A"], "OA": ["G"],
        }),
        sample_column_dict={
            "chr_col": "CHR", "pos_col": "POS", "ea_col": "EA", "oa_col": "OA",
        },
        policies=default_policies(),
    )

    assert frame.schema["POS"] == pl.Int64
    assert frame.item(0, "POS") == 100000


def test_export_summary_retains_numeric_and_text_statistics():
    frame = pl.DataFrame({
        "CHR": ["1", "1", "1"],
        "BETA": [0.1, None, -0.2],
        "SNP": ["a", "a", None],
    })
    summary = summarise_gwas2vcf_columns(
        frame,
        {"chr_col": "CHR", "beta_col": "BETA", "snp_id_col": "SNP"},
        "1",
    )
    rows = {row["key"]: row for row in summary.to_dicts()}

    assert "n_unique" not in summary.columns
    assert rows["beta_col"]["n_missing"] == 1
    assert rows["beta_col"]["min"] == pytest.approx(-0.2)
    assert rows["beta_col"]["max"] == pytest.approx(0.1)
    assert rows["beta_col"]["mean"] == pytest.approx(-0.05)
    assert rows["beta_col"]["median"] == pytest.approx(-0.05)
    assert rows["beta_col"]["std"] == pytest.approx(0.21213203435596428)
    assert rows["snp_id_col"]["n_missing"] == 1
    assert rows["snp_id_col"]["mean"] is None


def test_reject_collector_keeps_exact_complement_and_original_order(tmp_path):
    frame = pl.DataFrame({
        "variant": ["a", "b", "c", "d"], "value": [1, None, -1, 2],
    }).with_row_index(SOURCE_INPUT_ROW_COLUMN, offset=1)
    collector = RejectCollector(
        source_snapshot=frame,
        logger=None,
        out_path=str(tmp_path / "rejects.tsv"),
        delimiter="\t",
        compress=False,
    )

    kept = collector.reject(
        frame,
        pl.col("value") > 0,
        "test",
        "invalid_position",
    )

    assert kept["variant"].to_list() == ["c"]
    assert collector.frame()["variant"].to_list() == ["a", "b", "d"]
    assert kept.height + collector.total() == frame.height


def test_odds_ratio_conversion_qc_still_describes_log_odds():
    frame, qc, mapping = harmonise_effect_estimates(
        "1",
        pl.DataFrame({"OR": [0.5, 1.0, 2.0]}),
        {"beta_or_col": "OR", "se_col": None},
        decision="odds_ratio",
        policies=default_policies(),
    )

    assert frame["beta"].to_list() == pytest.approx([
        -0.6931471805599453, 0.0, 0.6931471805599453,
    ])
    assert mapping["beta_col"] == "beta"
    assert mapping["beta_or_col"] == "beta"
    assert qc["post_conversion_min"] == pytest.approx(qc["min_beta"])
    assert qc["post_conversion_max"] == pytest.approx(qc["max_beta"])
    assert qc["post_conversion_mean"] == pytest.approx(qc["mean_beta"])
    assert qc["post_conversion_std"] == pytest.approx(qc["std_beta"])


def test_odds_ratio_is_converted_before_standard_error_is_derived_from_z():
    expected_se = 0.1
    log_or = 0.6931471805599453
    frame = pl.DataFrame({"OR": [2.0], "Z": [log_or / expected_se]})
    mapping = {"beta_or_col": "OR", "se_col": None, "imp_z_col": "Z"}

    frame, effect_qc, mapping = harmonise_effect_estimates(
        "1",
        frame,
        mapping,
        decision="odds_ratio",
        policies=default_policies(),
    )
    frame, z_qc, mapping = derive_effect_and_standard_error_from_z(
        "1",
        frame,
        mapping,
        policies=default_policies(),
    )

    assert effect_qc["conversion"] == "OR_to_Beta_log_transform_applied"
    assert frame.item(0, "beta") == pytest.approx(log_or)
    assert frame.item(0, "SE") == pytest.approx(expected_se)
    assert mapping["beta_or_col"] == "beta"
    assert mapping["beta_col"] == "beta"
    assert z_qc["status"] == "SE derived from BETA and Z"


def test_info_aggregate_counts_match_configured_actions():
    policies = default_policies().with_overrides({
        "info": {
            "source": "study",
            "out_of_range": "clip",
            "on_missing": "keep",
        }
    })
    frame = pl.DataFrame({
        "CHR": ["1"] * 5,
        "POS": [1, 2, 3, 4, 5],
        "EA": ["A"] * 5,
        "OA": ["G"] * 5,
        "INFO": [None, -0.1, 0.5, 1.02, 1.2],
    })
    result, qc, _mapping = harmonise_imputation_quality(
        "1",
        frame,
        {"chr_col": "CHR", "pos_col": "POS", "ea_col": "EA", "oa_col": "OA", "imp_info_col": "INFO"},
        policies=policies,
    )

    assert qc["missing_info"] == 1
    assert qc["lt_zero"] == 1
    assert qc["gt_one"] == 2
    assert qc["gt_1_05"] == 1
    assert qc["rescaled_within_tolerance"] == 1
    assert qc["out_of_range"] == 2
    assert qc["missing_info_after"] == 1
    assert result["INFO"].to_list() == [None, 0.0, 0.5, 1.0, 1.0]


def test_explicit_fixed_info_populates_every_variant_with_provenance():
    frame = pl.DataFrame({
        "CHR": ["1", "1"],
        "POS": [10, 20],
        "EA": ["A", "C"],
        "OA": ["G", "T"],
    })
    result, qc, mapping = harmonise_imputation_quality(
        "1",
        frame,
        {
            "chr_col": "CHR",
            "pos_col": "POS",
            "ea_col": "EA",
            "oa_col": "OA",
            "imp_info_col": None,
            "fixed_info": 0.99,
            "fixed_info_column": "__postgwas_fixed_info",
        },
        policies=default_policies(),
    )

    assert result["__postgwas_fixed_info"].to_list() == [0.99, 0.99]
    assert mapping["imp_info_col"] == "__postgwas_fixed_info"
    assert qc["source"] == "fixed_cli:0.99"
    assert qc["fixed_info"] == 0.99
    assert qc["fixed_info_user_assigned"] is True
    assert qc["missing_info"] == 0


def test_resolved_internal_info_precedence_is_final_for_cli_datasets():
    frame = pl.DataFrame({
        "CHR": ["1"],
        "POS": [10],
        "EA": ["A"],
        "OA": ["G"],
        "INFO": [0.93],
    })
    policies = default_policies().with_overrides({
        "info": {"source": "reference"}
    })
    result, qc, mapping = harmonise_imputation_quality(
        "1",
        frame,
        {
            "chr_col": "CHR",
            "pos_col": "POS",
            "ea_col": "EA",
            "oa_col": "OA",
            "imp_info_col": "INFO",
            "info_source": "internal",
        },
        policies=policies,
    )

    assert result["INFO"].to_list() == [0.93]
    assert mapping["imp_info_col"] == "INFO"
    assert qc["source"] == "study_column:INFO"


def test_sample_size_aggregate_retains_neff_and_per_column_qc():
    frame, qc, mapping = harmonise_sample_sizes(
        "1",
        pl.DataFrame({"NCASE": [100, 200], "NCONTROL": [300, 400]}),
        {
            "ncase_col": "NCASE",
            "ncontrol_col": "NCONTROL",
            "ncase": None,
            "ncontrol": None,
        },
        policies=default_policies(),
    )

    assert frame["Neff"].to_list() == [300, 533]
    assert mapping["neff_col"] == "Neff"
    assert qc["sample_size_stats"]["NCASE"]["min"] == 100
    assert qc["sample_size_stats"]["NCONTROL"]["max"] == 400
    assert qc["sample_size_stats"]["Neff"]["median"] == pytest.approx(416.5)


def test_qc_table_from_memory_matches_saved_json(tmp_path):
    data = {
        "1": {"effect_qc": {"final_variants": 10}},
        "total_variant_read": 10,
    }
    path = tmp_path / "qc.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    allowed = default_policies().get("chromosome.allowed")

    from_memory = qc_results_to_dataframe(data=data, allowed_chromosomes=allowed)
    from_disk = qc_results_to_dataframe(path, allowed_chromosomes=allowed)

    assert from_memory.equals(from_disk)


def test_chromosomewise_harmonisation_metrics_are_written_as_tsv(tmp_path):
    path = tmp_path / "study_chromosomewise_harmonisation_metrics.tsv"
    table = _save_qc_results(
        {
            "1": {"effect_qc": {"final_variants": 10}},
            "total_variant_read": 10,
        },
        path,
        allowed_chromosomes=["1", "2"],
        delimiter="\t",
        null_output="",
    )

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert list(table.columns) == ["section", "metric", "chr1", "GLOBAL_VALUE"]
    assert rows == [
        {
            "section": "effect_qc",
            "metric": "final_variants",
            "chr1": "10",
            "GLOBAL_VALUE": "",
        },
        {
            "section": "GLOBAL",
            "metric": "total_variant_read",
            "chr1": "",
            "GLOBAL_VALUE": "10",
        },
    ]


def test_resolved_executable_preflight_is_recorded_in_engine_defaults():
    from postgwas.modules.harmonisation.cli import _engine_defaults

    resolved = {
        "bash": "/tools/bash",
        "bcftools": "/tools/bcftools",
        "python": "/tools/python",
        "tabix": "/tools/tabix",
    }
    defaults = _engine_defaults(
        load_configuration(), resolved_executables=resolved,
    )

    assert defaults["executables"] == resolved
    assert defaults["executables_prevalidated"] is True
