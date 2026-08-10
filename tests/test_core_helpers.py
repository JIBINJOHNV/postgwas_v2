from io import StringIO
import sys

import polars as pl
import pytest
from rich.console import Console

from postgwas.core.io.delimiters import resolve_delimiter
from postgwas.core.io.tables import read_delimited_table
from postgwas.core.paths import configured_output_matches
from postgwas.core.processes import build_container_command, run_checked_command
from postgwas.core.resource_preparation import (
    JsonLinesRunLogger,
    ResourcePreparationError,
    resumable_download,
    sha256,
)
from postgwas.core.dataframes import (
    chromosome_expression,
    count_matching_rows,
    numeric_column,
    position_expression,
)
from postgwas.core.values import (
    format_fraction_percentage,
    format_percentage,
    optional_text,
    parse_integer,
)
from postgwas.core.ui import StageProgress
from postgwas.modules.harmonisation.shared.runtime import reject_rows


def test_optional_values_have_one_normalization_rule():
    assert optional_text("  beta  ") == "beta"
    assert optional_text("NaN") is None
    assert optional_text(None) is None
    assert parse_integer("38,691") == 38691
    assert parse_integer("not-a-count") is None


def test_stage_progress_shows_current_completed_and_failed_stages():
    stream = StringIO()
    progress = StageProgress(
        "Test analysis progress",
        enabled=True,
        outcome_label_width=42,
        console=Console(
            file=stream,
            force_terminal=True,
            color_system=None,
            width=120,
        ),
    )

    progress.start_step(1, 3, "Validate inputs")
    progress.complete_step(
        1,
        3,
        "Validate inputs",
        "10 records validated; none excluded.",
        [
            ("count", "Records read", 10),
            ("success", "Records retained", 10),
        ],
    )
    progress.start_step(2, 3, "Run analysis")
    progress.complete_step(
        2,
        3,
        "Run analysis",
        outcome_fields=[("count", "A", 8)],
    )
    progress.start_step(3, 3, "Validate results")
    progress.fail_step(3, 3, "Validate results")

    text = stream.getvalue()
    assert "Test analysis progress" in text
    assert "Current 1/3 · Validate inputs" in text
    assert "Current 2/3 · Run analysis" in text
    assert "Completed 1/3 · Validate inputs" in text
    assert "Outcome" in text
    assert all(
        line == line.rstrip()
        for line in text.splitlines()
        if "Outcome" in line
    )
    assert "Records read" in text
    assert "Records retained" in text
    outcome_lines = [
        line
        for line in text.splitlines()
        if "Records read" in line or "Records retained" in line or "🧮  A " in line
    ]
    assert len({line.index(" : ") for line in outcome_lines}) == 1
    assert "Failed 3/3 · Validate results" in text
    assert "2/3" in text


def test_percentage_formatters_handle_missing_and_zero_denominators():
    assert format_percentage(1, 4) == "25.00%"
    assert format_percentage(1, 0) == "0.00%"
    assert format_fraction_percentage(0.125) == "12.50%"
    assert format_fraction_percentage(None) == "not available"


def test_dataframe_helpers_make_null_polarity_explicit():
    frame = pl.DataFrame({"value": [1, None, -1]})
    condition = pl.col("value") > 0
    assert count_matching_rows(frame, condition) == 1
    assert count_matching_rows(frame, condition, null_is_match=True) == 2
    converted = frame.select(numeric_column(frame, "value").alias("numeric"))
    assert converted.schema["numeric"] == pl.Float64


def test_coordinate_normalization_preserves_integer_coordinates():
    frame = pl.DataFrame({"chrom": ["chr01", "X"], "pos": [" 123 ", "bad"]})
    normalized = frame.with_columns(
        chromosome_expression(
            "chrom", frame.schema["chrom"], strip_leading_zero=True
        ),
        position_expression("pos", frame.schema["pos"]),
    )
    assert normalized["chrom"].to_list() == ["1", "X"]
    assert normalized["pos"].to_list() == [123, None]


def test_shared_rejection_uses_null_as_rejected():
    frame = pl.DataFrame({"value": [1, None, -1]})
    kept, removed = reject_rows(
        frame,
        pl.col("value") > 0,
        step_label="test",
        reason="invalid_position",
    )
    assert removed == 2
    assert kept["value"].to_list() == [-1]


def test_container_command_uses_same_path_mounts_without_a_shell(tmp_path):
    inputs = tmp_path / "input"
    output = tmp_path / "output"
    inputs.mkdir()
    output.mkdir()
    command = build_container_command(
        "docker",
        "example/tool:1.0",
        ["python", "/tools/tool.py", "--input", str(inputs / "study.tsv")],
        mount_directories=[inputs, output, inputs],
        platform="linux/amd64",
        work_directory=output,
        use_host_user=False,
        validate=False,
    )

    assert command[:5] == ["docker", "run", "--rm", "--platform", "linux/amd64"]
    assert command.count("--volume") == 2
    assert "%s:%s" % (inputs, inputs) in command
    assert "%s:%s" % (output, output) in command
    assert command[-5:] == [
        "example/tool:1.0", "python", "/tools/tool.py", "--input",
        str(inputs / "study.tsv"),
    ]


def test_streamed_data_output_has_no_process_metadata(tmp_path):
    destination = tmp_path / "data.tsv"
    run_checked_command(
        ["printf", "value\\n"],
        "Writing data fixture",
        stdout_path=destination,
    )
    assert destination.read_text(encoding="utf-8") == "value\n"


def test_resumable_download_reuses_only_checksum_validated_archives(tmp_path):
    destination = tmp_path / "reference.zip"
    destination.write_bytes(b"verified archive")
    expected = sha256(destination)
    logger = JsonLinesRunLogger(tmp_path / "download.jsonl")

    observed = resumable_download(
        curl="not-used-for-a-valid-cache",
        url="https://example.invalid/reference.zip",
        destination=destination,
        retries=1,
        partial_suffix=".partial",
        logger=logger,
        expected_sha256=expected,
    )

    assert observed == destination
    assert '"event": "download_reused"' in logger.path.read_text(encoding="utf-8")
    destination.write_bytes(b"corrupt archive")
    with pytest.raises(ResourcePreparationError, match="failed SHA-256"):
        resumable_download(
            curl="not-used-for-an-invalid-cache",
            url="https://example.invalid/reference.zip",
            destination=destination,
            retries=1,
            partial_suffix=".partial",
            logger=logger,
            expected_sha256=expected,
        )


def test_process_transcript_records_command_outcome(tmp_path):
    destination = tmp_path / "tool.log"
    run_checked_command(
        ["printf", "tool output\\n"],
        "Writing process transcript",
        stdout_path=destination,
        stdout_header="tool=fixture\n",
    )
    assert destination.read_text(encoding="utf-8") == (
        "tool=fixture\ntool output\n\nexit_code=0\n"
    )


def test_failed_process_reports_stdout_and_stderr():
    with pytest.raises(RuntimeError) as error:
        run_checked_command(
            [
                sys.executable,
                "-c",
                "import sys; print('stdout diagnostic'); "
                "print('stderr diagnostic', file=sys.stderr); raise SystemExit(7)",
            ],
            "Failing fixture",
        )

    message = str(error.value)
    assert "exit status 7" in message
    assert "stdout diagnostic" in message
    assert "stderr diagnostic" in message


def test_configured_output_matches_uses_the_canonical_pattern(tmp_path):
    first = tmp_path / "study_chr1.vcf.gz"
    second = tmp_path / "study_chr2.vcf.gz"
    first.touch()
    second.touch()
    assert configured_output_matches(
        tmp_path,
        "{dataset_id}_chr{chromosome}.vcf.gz",
        dataset_id="study",
        chromosome="*",
    ) == [first, second]


def test_delimiter_resolution_uses_configured_candidates_and_bounds(tmp_path):
    table = tmp_path / "study.txt"
    table.write_text(
        "##study=example\nCHR|BP|P\n1|100|0.5\n1|200|0.1\n",
        encoding="utf-8",
    )
    result = resolve_delimiter(
        table,
        "auto",
        candidates=["tab", "pipe", "comma"],
        minimum_columns=2,
        maximum_columns=10,
        sample_lines=3,
    )
    assert result.value == "|"
    assert result.column_count == 3


def test_delimiter_resolution_never_silently_falls_back(tmp_path):
    table = tmp_path / "inconsistent.txt"
    table.write_text("A,B,C\n1,2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="No configured delimiter"):
        resolve_delimiter(
            table,
            "auto",
            candidates=["comma", "tab"],
            minimum_columns=2,
            maximum_columns=10,
            sample_lines=2,
        )


def test_shared_table_reader_uses_configured_delimiter_and_null_values(tmp_path):
    table = tmp_path / "reference.csv"
    table.write_text("marker,score\nrs1,NA\nrs2,0.9\n", encoding="utf-8")

    frame, detected = read_delimited_table(
        table,
        "comma",
        candidates=["comma"],
        minimum_columns=2,
        maximum_columns=2,
        sample_lines=2,
        null_values=["NA"],
        infer_schema_length=10,
    )

    assert detected.method == "configured"
    assert frame.columns == ["marker", "score"]
    assert frame["score"].to_list() == [None, 0.9]
