"""Failure-path tests for required rejected-variant provenance."""

import gzip

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
from postgwas.modules.harmonisation.service import (
    PipelineError,
    _reconcile_dataset_rows,
)


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
            batch_rows=2,
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

    monkeypatch.setattr(reject_module, "_write_reject_batches", fail_write)

    with pytest.raises(RejectOutputError, match="disk full"):
        concat_reject_files(
            [source],
            tmp_path / "combined.tsv.gz",
            delimiter="\t",
            batch_rows=2,
        )


def test_validated_combined_reject_removes_all_source_shards(tmp_path):
    input_reject = tmp_path / "input.tsv"
    chromosome_reject = tmp_path / "chromosome.tsv"
    output = tmp_path / "rejected" / "study_rejected_variants.tsv.gz"
    pl.DataFrame({
        "variant": ["rs1"],
        "reject_step": ["01"],
        "reject_detail": ["missing P"],
        "reject_reason": ["missing_read_mandatory_value"],
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
        batch_rows=2,
        remove_sources=True,
    )

    combined = pl.read_csv(output, separator="\t", infer_schema_length=0)
    assert combined["variant"].to_list() == ["rs1", "rs2"]
    assert result["rows"] == 2
    assert result["rows_by_source"] == {
        str(input_reject): 1,
        str(chromosome_reject): 1,
    }
    assert result["reject_counts_by_source"] == {
        str(input_reject): {"missing_read_mandatory_value": 1},
        str(chromosome_reject): {"palindromic_ambiguous": 1},
    }
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
            [source], output, delimiter="\t", batch_rows=2,
            remove_sources=True,
        )

    assert source.is_file()
    assert not output.exists()
    assert not list(output.parent.glob("*.part*"))


def test_combined_reject_streams_diagonal_schema_in_bounded_batches(
    tmp_path, monkeypatch,
):
    first = tmp_path / "input.tsv.gz"
    second = tmp_path / "chromosome.tsv"
    output = tmp_path / "combined.tsv.gz"
    input_frame = pl.DataFrame({
        "input_only": ["a", "b", "c"],
        "reject_step": ["01", "01", "01"],
        "reject_detail": ["first", "line one\nline two", None],
        "reject_reason": ["invalid_position"] * 3,
    })
    with gzip.open(first, "wb") as handle:
        input_frame.write_csv(handle, separator="\t")
    pl.DataFrame({
        "chromosome_only": ["x", "y"],
        "reject_step": ["03", "03"],
        "reject_detail": [None, "last"],
        "reject_reason": ["reference_unmatched"] * 2,
    }).write_csv(second, separator="\t")

    batch_sizes = []
    normalise = reject_module._normalise_reject_batch

    def record_batch_size(frame, columns):
        batch_sizes.append(frame.height)
        return normalise(frame, columns)

    def forbid_full_materialisation(*_args, **_kwargs):
        raise AssertionError("full reject materialisation must not be used")

    monkeypatch.setattr(
        reject_module, "_normalise_reject_batch", record_batch_size,
    )
    monkeypatch.setattr(
        reject_module.pl, "concat", forbid_full_materialisation,
    )

    result = concat_reject_files(
        [first, second], output, delimiter="\t", batch_rows=2,
    )

    combined = pl.read_csv(output, separator="\t", infer_schema_length=0)
    assert combined.columns == [
        "input_only",
        "chromosome_only",
        "reject_step",
        "reject_detail",
        "reject_reason",
    ]
    assert combined["input_only"].to_list() == ["a", "b", "c", None, None]
    assert combined["chromosome_only"].to_list() == [None, None, None, "x", "y"]
    assert combined["reject_detail"].to_list()[1] == "line one\nline two"
    assert combined["reject_reason"].to_list() == [
        "invalid_position",
        "invalid_position",
        "invalid_position",
        "reference_unmatched",
        "reference_unmatched",
    ]
    assert result["rows"] == 5
    assert result["rows_by_source"] == {str(first): 3, str(second): 2}
    assert result["batch_rows"] == 2
    assert batch_sizes
    assert max(batch_sizes) <= 2


def test_combined_reject_preserves_union_header_when_every_shard_is_empty(
    tmp_path,
):
    first = tmp_path / "input.tsv"
    second = tmp_path / "chromosome.tsv"
    output = tmp_path / "combined.tsv.gz"
    pl.DataFrame(schema={
        "input_only": pl.String,
        "reject_step": pl.String,
        "reject_detail": pl.String,
        "reject_reason": pl.String,
    }).write_csv(first, separator="\t")
    pl.DataFrame(schema={
        "chromosome_only": pl.String,
        "reject_step": pl.String,
        "reject_detail": pl.String,
        "reject_reason": pl.String,
    }).write_csv(second, separator="\t")

    result = concat_reject_files(
        [first, second], output, delimiter="\t", batch_rows=2,
    )

    combined = pl.read_csv(output, separator="\t", infer_schema_length=0)
    assert combined.height == 0
    assert combined.columns == [
        "input_only",
        "chromosome_only",
        "reject_step",
        "reject_detail",
        "reject_reason",
    ]
    assert result["rows"] == 0
    assert result["rows_by_source"] == {str(first): 0, str(second): 0}
    assert result["reject_counts_by_source"] == {
        str(first): {}, str(second): {},
    }


@pytest.mark.parametrize("reason", [None, "not_a_registered_reason"])
def test_combined_reject_refuses_unattributable_rows(tmp_path, reason):
    malformed = tmp_path / "malformed.tsv"
    pl.DataFrame({
        "variant": ["rs1"],
        "reject_step": ["03"],
        "reject_detail": [None],
        "reject_reason": pl.Series([reason], dtype=pl.String),
    }).write_csv(malformed, separator="\t")

    with pytest.raises(RejectOutputError, match="reject_reason"):
        concat_reject_files(
            [malformed], tmp_path / "combined.tsv.gz",
            delimiter="\t", batch_rows=2,
        )


def test_combined_reject_rejects_malformed_source_header(tmp_path):
    malformed = tmp_path / "malformed.tsv"
    pl.DataFrame({
        "variant": ["rs1"],
        "reject_reason": ["invalid_position"],
    }).write_csv(malformed, separator="\t")

    with pytest.raises(RejectOutputError, match="missing required column"):
        concat_reject_files(
            [malformed], tmp_path / "combined.tsv.gz",
            delimiter="\t", batch_rows=2,
        )


@pytest.mark.parametrize(
    "batch_rows", [None, False, 0, -1, 1.5, "not-an-integer"],
)
def test_combined_reject_requires_positive_batch_rows(tmp_path, batch_rows):
    with pytest.raises(RejectOutputError, match="positive integer"):
        concat_reject_files(
            [], tmp_path / "combined.tsv.gz",
            delimiter="\t", batch_rows=batch_rows,
        )


def test_complete_dataset_row_reconciliation_covers_every_parsed_row():
    result = _reconcile_dataset_rows(
        parsed_rows=1000,
        ready_rows=900,
        partition_rows={"1": 500, "2": 400},
        completed=["1", "2"],
        failed=[],
        chromosome_summaries={
            "1": {"rows_in": 500, "rejected": 20, "rows_out": 480},
            "2": {"rows_in": 400, "rejected": 30, "rows_out": 370},
        },
        reject_rows_by_source={"input": 100, "1": 20, "2": 30},
        combined_reject_rows=150,
    )

    assert result["dataset_stage_rejected"] == 100
    assert result["chromosome_stage_rejected"] == 50
    assert result["rejected"] == 150
    assert result["rows_exported"] == 850
    assert result["unprocessed_failed_chromosome_rows"] == 0
    assert result["accounted_rows"] == 1000
    assert result["balanced"] is True
    assert result["complete"] is True


def test_partial_dataset_labels_failed_chromosome_survivors_as_unprocessed():
    result = _reconcile_dataset_rows(
        parsed_rows=1000,
        ready_rows=900,
        partition_rows={"1": 500, "2": 400},
        completed=["1"],
        failed=["2"],
        chromosome_summaries={
            "1": {"rows_in": 500, "rejected": 20, "rows_out": 480},
            "2": {"rows_in": 400, "rejected": 10, "status": "failed"},
        },
        reject_rows_by_source={"input": 100, "1": 20, "2": 10},
        combined_reject_rows=130,
    )

    assert result["rows_exported"] == 480
    assert result["rejected"] == 130
    assert result["unprocessed_rows_by_chromosome"] == {"2": 390}
    assert result["unprocessed_failed_chromosome_rows"] == 390
    assert result["accounted_rows"] == 1000
    assert result["balanced"] is True
    assert result["complete"] is False


def test_dataset_row_reconciliation_refuses_missing_chromosome_evidence():
    with pytest.raises(
        PipelineError,
        match=r"chromosome_summaries\[1\]\.rows_out",
    ):
        _reconcile_dataset_rows(
            parsed_rows=10,
            ready_rows=9,
            partition_rows={"1": 9},
            completed=["1"],
            failed=[],
            chromosome_summaries={
                "1": {"rows_in": 9, "rejected": 1},
            },
            reject_rows_by_source={"input": 1, "1": 1},
            combined_reject_rows=2,
        )


def test_dataset_row_reconciliation_refuses_incomplete_status_coverage():
    with pytest.raises(PipelineError, match="cover exactly"):
        _reconcile_dataset_rows(
            parsed_rows=10,
            ready_rows=9,
            partition_rows={"1": 5, "2": 4},
            completed=["1"],
            failed=[],
            chromosome_summaries={
                "1": {"rows_in": 5, "rejected": 0, "rows_out": 5},
            },
            reject_rows_by_source={"input": 1, "1": 0, "2": 0},
            combined_reject_rows=1,
        )
