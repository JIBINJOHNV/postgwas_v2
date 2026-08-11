"""Failure-path tests for required rejected-variant provenance."""

import polars as pl
import pytest

from postgwas.modules.harmonisation import rejects as reject_module
from postgwas.modules.harmonisation.rejects import (
    RejectCollector,
    RejectOutputError,
    SOURCE_INPUT_ROW_COLUMN,
    concat_reject_files,
)
from postgwas.modules.harmonisation.service import _finalize_chromosome_rejects


class _Logger:
    def __init__(self):
        self.errors = []

    def error(self, message):
        self.errors.append(message)

    def qc(self, *_args, **_kwargs):
        pass


def test_reject_flush_failure_is_fatal_even_with_a_logger(tmp_path, monkeypatch):
    logger = _Logger()
    source = pl.DataFrame({"variant": ["rs1"]}).with_row_index(
        SOURCE_INPUT_ROW_COLUMN, offset=1,
    )
    collector = RejectCollector(
        source_snapshot=source,
        logger=logger,
        out_path=str(tmp_path / "rejected.tsv"),
        delimiter="\t",
        compress=False,
    )
    collector.reject(
        source,
        pl.lit(True),
        "test",
        "invalid_position",
    )

    def fail_write(*_args, **_kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(reject_module, "_write_table", fail_write)

    with pytest.raises(RejectOutputError, match="disk full"):
        collector.flush()

    assert collector.summary()["flushed"] is False
    assert logger.errors


def test_rejected_rows_restore_immutable_values_after_allele_transformations(tmp_path):
    source = pl.DataFrame({
        "EA": ["A"],
        "OA": ["G"],
        "EAF": [0.2],
        "OR": [2.0],
        "SE": [0.1],
    }).with_row_index(SOURCE_INPUT_ROW_COLUMN, offset=1)
    snapshot_path = tmp_path / "source.parquet"
    source.write_parquet(snapshot_path)
    collector = RejectCollector(
        source_snapshot=snapshot_path,
        logger=None,
        out_path=str(tmp_path / "rejected.tsv"),
        delimiter="\t",
        compress=False,
    )
    transformed = source.with_columns([
        pl.lit("G").alias("EA"),
        pl.lit("A").alias("OA"),
        pl.lit(0.8).alias("EAF"),
        pl.lit(0.5).alias("OR"),
        pl.lit(0.05).alias("SE"),
    ])

    collector.reject(
        transformed, pl.lit(True), "10 effect validation", "se_non_positive",
    )
    rejected = collector.frame()

    assert rejected.select(["EA", "OA", "EAF", "OR", "SE"]).row(0) == (
        "A", "G", "0.2", "2.0", "0.1",
    )
    assert SOURCE_INPUT_ROW_COLUMN not in rejected.columns


def test_rejection_fails_when_working_row_loses_stable_source_id(tmp_path):
    source = pl.DataFrame({"variant": ["rs1"]}).with_row_index(
        SOURCE_INPUT_ROW_COLUMN, offset=1,
    )
    collector = RejectCollector(
        source_snapshot=source,
        logger=None,
        out_path=str(tmp_path / "rejected.tsv"),
        delimiter="\t",
        compress=False,
    )

    with pytest.raises(RejectOutputError, match="no stable source-row ID"):
        collector.reject(
            source.drop(SOURCE_INPUT_ROW_COLUMN),
            pl.lit(True),
            "test",
            "invalid_position",
        )


def test_successful_chromosome_becomes_failed_when_reject_flush_fails():
    class FailingCollector:
        def flush(self):
            raise RejectOutputError("permission denied")

        def summary(self):
            return {"path": "/unwritten/rejected.tsv", "flushed": False}

    qc = {}
    status = _finalize_chromosome_rejects(FailingCollector(), "ok", qc)

    assert status == "failed"
    assert "permission denied" in qc["error"]
    assert qc["reject_output_error"] == qc["error"]
    assert qc["rejects"]["flushed"] is False


def test_reject_flush_does_not_replace_an_existing_analysis_error():
    class FailingCollector:
        def flush(self):
            raise RejectOutputError("permission denied")

        def summary(self):
            return {"path": "/unwritten/rejected.tsv", "flushed": False}

    qc = {"error": "SchemaError: invalid numeric column"}
    status = _finalize_chromosome_rejects(FailingCollector(), "failed", qc)

    assert status == "failed"
    assert qc["error"] == "SchemaError: invalid numeric column"
    assert "permission denied" in qc["reject_output_error"]


def test_missing_reject_input_fails_dataset_concatenation(tmp_path):
    with pytest.raises(RejectOutputError, match="missing or unreadable"):
        concat_reject_files(
            [tmp_path / "missing.tsv.gz"],
            tmp_path / "combined.tsv.gz",
            delimiter="\t",
        )

    assert not (tmp_path / "combined.tsv.gz").exists()


def test_combined_reject_write_failure_is_fatal(tmp_path, monkeypatch):
    source = tmp_path / "chromosome.tsv"
    pl.DataFrame({
        "variant": ["rs1"],
        "reject_step": ["03"],
        "reject_detail": [None],
        "reject_reason": ["reference_unmatched"],
    }).write_csv(source, separator="\t")

    def fail_write(*_args, **_kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(reject_module, "_write_table", fail_write)

    with pytest.raises(RejectOutputError, match="disk full"):
        concat_reject_files(
            [source],
            tmp_path / "combined.tsv.gz",
            delimiter="\t",
        )


def test_validated_combined_reject_removes_all_source_shards(tmp_path):
    input_reject = tmp_path / "input.tsv"
    chromosome_reject = tmp_path / "chromosome.tsv"
    output = tmp_path / "rejected" / "study_rejected_variants.tsv.gz"
    pl.DataFrame({
        "variant": ["rs1"],
        "reject_step": ["01"],
        "reject_detail": ["missing P"],
        "reject_reason": ["missing_required_columns"],
    }).write_csv(input_reject, separator="\t")
    pl.DataFrame({
        "variant": ["rs2"],
        "reject_step": ["04"],
        "reject_detail": ["ambiguous alleles"],
        "reject_reason": ["palindromic_ambiguous"],
    }).write_csv(chromosome_reject, separator="\t")

    result = concat_reject_files(
        [input_reject, chromosome_reject],
        output,
        delimiter="\t",
        remove_sources=True,
    )

    combined = pl.read_csv(output, separator="\t", infer_schema_length=0)
    assert combined["variant"].to_list() == ["rs1", "rs2"]
    assert result["rows"] == 2
    assert set(result["source_files_removed"]) == {
        str(input_reject), str(chromosome_reject),
    }
    assert result["source_files_retained"] == []
    assert not input_reject.exists()
    assert not chromosome_reject.exists()


def test_failed_combined_reject_validation_preserves_source_shards(
    tmp_path, monkeypatch,
):
    source = tmp_path / "chromosome.tsv"
    output = tmp_path / "rejected" / "study_rejected_variants.tsv.gz"
    pl.DataFrame({
        "variant": ["rs1"],
        "reject_step": ["04"],
        "reject_detail": [None],
        "reject_reason": ["reference_unmatched"],
    }).write_csv(source, separator="\t")

    def fail_validation(*_args, **_kwargs):
        raise RejectOutputError("row-count validation failed")

    monkeypatch.setattr(
        reject_module, "_validate_reject_output", fail_validation,
    )

    with pytest.raises(RejectOutputError, match="row-count validation failed"):
        concat_reject_files(
            [source], output, delimiter="\t", remove_sources=True,
        )

    assert source.is_file()
    assert not output.exists()
    assert not list(output.parent.glob("*.part*"))
