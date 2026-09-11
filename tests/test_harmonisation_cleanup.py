"""Regression tests for retry-safe harmonisation cleanup."""

import csv
from pathlib import Path

import polars as pl
import pytest

from postgwas.config import load_configuration
from postgwas.core.paths import configured_output_path
from postgwas.modules.harmonisation import cleanup as cleanup_module
from postgwas.modules.harmonisation.cleanup import (
    _merge_adapter_summaries,
    remove_merged_gwas2vcf_intermediate,
    remove_partial_chromosome_outputs,
)
from postgwas.modules.harmonisation.gwas2vcf_export import (
    GWAS2VCF_SUMMARY_COLUMNS,
)


def _write_adapter_summary(
    path, chromosome, num_rows, delimiter="\t", status="success",
):
    row = {
        "chromosome": str(chromosome),
        "key": "snp_id_col",
        "column_name": "SNP",
        "dtype": "String",
        "n_missing": 0,
        "min": "",
        "max": "",
        "mean": "",
        "median": "",
        "std": "",
        "num_rows": num_rows,
        "num_cols": 12,
        "status": status,
    }
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(GWAS2VCF_SUMMARY_COLUMNS),
            delimiter=delimiter,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(row)


def test_raw_gwas2vcf_cleanup_removes_only_the_exact_vcf_and_indexes(tmp_path):
    raw_vcf = tmp_path / "study_gwas2vcf_GRCh37_merged.vcf.gz"
    raw_tbi = Path(str(raw_vcf) + ".tbi")
    final_vcf = tmp_path / "study_GRCh37_merged.vcf.gz"
    other_study = tmp_path / "other_gwas2vcf_GRCh37_merged.vcf.gz"
    for path in (raw_vcf, raw_tbi, final_vcf, other_study):
        path.write_bytes(b"validated\n")

    removed = remove_merged_gwas2vcf_intermediate(
        raw_vcf, output_directory=tmp_path,
    )

    assert removed == [str(raw_vcf), str(raw_tbi)]
    assert not raw_vcf.exists()
    assert not raw_tbi.exists()
    assert final_vcf.is_file()
    assert other_study.is_file()


def test_raw_gwas2vcf_cleanup_refuses_paths_outside_output_directory(tmp_path):
    output_directory = tmp_path / "dataset"
    output_directory.mkdir()
    outside = tmp_path / "study_gwas2vcf_GRCh37_merged.vcf.gz"
    outside.write_bytes(b"validated\n")

    with pytest.raises(RuntimeError, match="outside the dataset output directory"):
        remove_merged_gwas2vcf_intermediate(
            outside, output_directory=output_directory,
        )

    assert outside.is_file()


def test_raw_gwas2vcf_cleanup_validates_indexes_before_deleting_vcf(tmp_path):
    raw_vcf = tmp_path / "study_gwas2vcf_GRCh37_merged.vcf.gz"
    raw_vcf.write_bytes(b"validated\n")
    Path(str(raw_vcf) + ".tbi").mkdir()

    with pytest.raises(RuntimeError, match="invalid index path"):
        remove_merged_gwas2vcf_intermediate(
            raw_vcf, output_directory=tmp_path,
        )

    assert raw_vcf.is_file()


def test_retry_cleanup_preserves_inputs_and_removes_partial_outputs(tmp_path):
    layout = dict(load_configuration().modules.harmonisation.output_layout.root)
    values = {
        "dataset_id": "study",
        "chromosome": "1",
        "build": "GRCh37",
        "target_build": "GRCh38",
    }
    chromosome_table = configured_output_path(
        tmp_path, layout["chromosome_table"], **values,
    )
    source_snapshot = configured_output_path(
        tmp_path, layout["chromosome_source_snapshot"], **values,
    )
    rejected = configured_output_path(
        tmp_path, layout["chromosome_reject"], **values,
    )
    post_orientation_duplicates = configured_output_path(
        tmp_path, layout["post_orientation_duplicates"], **values,
    )
    adapter_input = configured_output_path(
        tmp_path, layout["adapter_input"], **values,
    )
    partial_vcf = configured_output_path(
        tmp_path, layout["chromosome_raw_vcf"], **values,
    )
    partial_index = Path(str(partial_vcf) + ".tbi")

    for path in (
        source_snapshot,
        rejected,
        post_orientation_duplicates,
        adapter_input,
        partial_vcf,
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"partial\n")
    chromosome_table.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame({"CHR": ["1"], "POS": [100]}).write_parquet(chromosome_table)
    partial_index.write_bytes(b"index\n")

    removed = remove_partial_chromosome_outputs(
        tmp_path, "study", "1", layout,
    )

    assert removed == 5
    assert chromosome_table.is_file()
    retained = pl.read_parquet(chromosome_table)
    assert retained.height == 1
    assert source_snapshot.is_file()
    assert not rejected.exists()
    assert not post_orientation_duplicates.exists()
    assert not adapter_input.exists()
    assert not partial_vcf.exists()
    assert not partial_index.exists()


def test_adapter_summary_merge_writes_one_header_and_uses_configured_delimiter(
    tmp_path,
):
    chromosome_one = tmp_path / "study_chr1_summary.txt"
    chromosome_two = tmp_path / "study_chr2_summary.txt"
    destination = tmp_path / "study_summary.txt"
    _write_adapter_summary(chromosome_one, "1", 3, delimiter="|")
    _write_adapter_summary(chromosome_two, "2", 5, delimiter="|")

    result = _merge_adapter_summaries(
        [("1", chromosome_one), ("2", chromosome_two)],
        destination,
        "|",
        {"1": 3, "2": 5},
    )

    lines = destination.read_text(encoding="utf-8").splitlines()
    expected_header = "|".join(GWAS2VCF_SUMMARY_COLUMNS)
    assert lines.count(expected_header) == 1
    with destination.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="|"))
    assert [row["chromosome"] for row in rows] == ["1", "2"]
    assert result["rows_by_chromosome"] == {"1": 3, "2": 5}
    assert result["total_rows"] == 8


def test_adapter_summary_merge_preserves_audit_when_counts_disagree(tmp_path):
    source = tmp_path / "study_chr1_summary.tsv"
    destination = tmp_path / "study_summary.tsv"
    _write_adapter_summary(source, "1", 7)
    destination.write_text("previous validated audit\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="row accounting disagree"):
        _merge_adapter_summaries(
            [("1", source)], destination, "\t", {"1": 6},
        )

    assert destination.read_text(encoding="utf-8") == (
        "previous validated audit\n"
    )
    assert source.is_file()
    assert list(tmp_path.glob("study_summary.tsv.part*")) == []


def test_adapter_summary_merge_preserves_audit_when_atomic_replace_fails(
    monkeypatch, tmp_path,
):
    source = tmp_path / "study_chr1_summary.tsv"
    destination = tmp_path / "study_summary.tsv"
    _write_adapter_summary(source, "1", 4)
    destination.write_text("previous validated audit\n", encoding="utf-8")

    def fail_replace(_source, _destination):
        raise OSError("simulated replace failure")

    monkeypatch.setattr(cleanup_module.os, "replace", fail_replace)
    with pytest.raises(RuntimeError, match="previous merged audit"):
        _merge_adapter_summaries(
            [("1", source)], destination, "\t", {"1": 4},
        )

    assert destination.read_text(encoding="utf-8") == (
        "previous validated audit\n"
    )
    assert source.is_file()
    assert list(tmp_path.glob("study_summary.tsv.part*")) == []
