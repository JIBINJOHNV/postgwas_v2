from io import StringIO
import sys

import polars as pl
import pytest
from rich.console import Console
from rich.cells import cell_len

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
from postgwas.core.ui import MeasuredProgress, StageProgress, run_with_progress
from postgwas.core.ui.screen import screen_field
from postgwas.modules.harmonisation.shared.runtime import reject_rows


def test_optional_values_have_one_normalization_rule():
    assert optional_text("  beta  ") == "beta"
    assert optional_text("NaN") is None
    assert optional_text(None) is None
    assert parse_integer("38,691") == 38691
    assert parse_integer("not-a-count") is None


def test_screen_fields_align_unicode_symbols_and_wrap_long_values():
    count = screen_field(
        "count", "Count label", "10", width=80, indent=10, label_width=42,
    )
    warning = screen_field(
        "warning", "Warning label", "1", width=80, indent=10, label_width=42,
    )
    assert cell_len(count.split(":", 1)[0]) == cell_len(
        warning.split(":", 1)[0]
    )

    rendered = screen_field(
        "info", "Long path",
        "/Users/example/" + "a" * 100 + "/report.yaml",
        width=80,
        indent=10,
        label_width=42,
        break_long_values=True,
    )
    lines = rendered.splitlines()
    assert len(lines) > 1
    assert all(cell_len(line) <= 80 for line in lines)
    assert all(line.startswith(" " * 59) for line in lines[1:])

    long_label = screen_field(
        "loss",
        "Palindromic SNPs with ambiguous frequency (0.4 ≤ AF ≤ 0.6)",
        "matched 153,472; removed for this reason 150,313",
        width=120,
        indent=8,
        label_width=42,
    )
    long_lines = long_label.splitlines()
    assert len(long_lines) == 2
    assert " : " not in long_lines[0]
    aligned_reference = screen_field(
        "count", "Count label", "10", width=120, indent=8, label_width=42,
    )
    assert cell_len(long_lines[1].split(" : ", 1)[0]) == cell_len(
        aligned_reference.split(" : ", 1)[0]
    )


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
    assert "Started 1/3 · Validate inputs" not in text
    assert "Current 1/3 · Validate inputs" in text
    assert "\n      🔬  Current 1/3 · Validate inputs\n" not in text
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


def test_measured_progress_does_not_leave_a_static_current_line():
    stream = StringIO()
    progress = MeasuredProgress(
        "Measured analysis progress",
        enabled=True,
        console=Console(
            file=stream,
            force_terminal=True,
            color_system=None,
            width=120,
        ),
    )

    progress.start("Run external analysis", total=10)
    progress.complete(10, total=10)

    text = stream.getvalue()
    assert "Current 0/10 · Run external analysis" in text
    assert "\n      🔬  Current 0/10 · Run external analysis\n" not in text


def test_measured_progress_records_an_observed_phase_change():
    stream = StringIO()
    progress = MeasuredProgress(
        "External phases",
        enabled=True,
        console=Console(
            file=stream,
            force_terminal=False,
            color_system=None,
            width=120,
        ),
    )

    progress.start("Load inputs")
    progress.set_phase("Match variants")
    progress.update(1, total=2)
    progress.complete(2, total=2)

    text = stream.getvalue()
    assert "Started 0/? · Load inputs" in text
    assert "Current 0/? · Match variants" in text
    assert "Progress 1/2 · Match variants · 50%" in text


def test_stage_progress_context_marks_success_and_failure_truthfully():
    stream = StringIO()
    progress = StageProgress(
        "Context-managed stages",
        enabled=True,
        console=Console(
            file=stream,
            force_terminal=False,
            color_system=None,
            width=120,
        ),
    )

    with progress.step(1, 2, "Prepare inputs"):
        pass
    with pytest.raises(ValueError, match="invalid result"):
        with progress.step(2, 2, "Validate results"):
            raise ValueError("invalid result")

    text = stream.getvalue()
    assert "Completed 1/2 · Prepare inputs" in text
    assert "Failed 2/2 · Validate results" in text
    assert "All 2 stages completed" not in text


def test_top_level_progress_reports_success_failure_and_disabled_runs():
    stream = StringIO()
    console = Console(
        file=stream,
        force_terminal=True,
        color_system=None,
        width=120,
    )

    assert run_with_progress(
        lambda: 0,
        label="Module progress",
        title="Run module",
        enabled=True,
        console=console,
    ) == 0
    assert "Completed 1/1 · Run module" in stream.getvalue()

    with pytest.raises(ValueError, match="failed operation"):
        run_with_progress(
            lambda: (_ for _ in ()).throw(ValueError("failed operation")),
            label="Module progress",
            title="Run module",
            enabled=True,
            console=console,
        )
    assert "Failed 1/1 · Run module" in stream.getvalue()

    with pytest.raises(SystemExit) as successful_exit:
        run_with_progress(
            lambda: (_ for _ in ()).throw(SystemExit(0)),
            label="Module progress",
            title="Run module",
            enabled=True,
            console=console,
        )
    assert successful_exit.value.code == 0
    assert stream.getvalue().count("Completed 1/1 · Run module") >= 2

    hidden = StringIO()
    assert run_with_progress(
        lambda: "result",
        label="Module progress",
        title="Run module",
        enabled=False,
        console=Console(file=hidden, force_terminal=True, color_system=None),
    ) == "result"
    assert hidden.getvalue() == ""


def test_detailed_progress_can_run_inside_mandatory_module_progress():
    stream = StringIO()
    console = Console(
        file=stream,
        force_terminal=True,
        color_system=None,
        width=120,
    )
    outer = StageProgress("Module progress", enabled=True, console=console)
    inner = StageProgress("Scientific stages", enabled=True, console=console)

    outer.start_step(1, 1, "Run module")
    inner.start_step(1, 2, "Validate inputs")
    inner.complete_step(1, 2, "Validate inputs")
    inner.start_step(2, 2, "Calculate associations")
    inner.complete_step(2, 2, "Calculate associations")
    outer.complete_step(1, 1, "Run module")

    text = stream.getvalue()
    assert "Scientific stages" in text
    assert "Completed 2/2 · Calculate associations" in text
    assert "Completed 1/1 · Run module" in text
    assert not any(
        "Current 1/1 · Run module" in line
        and "Current 1/2 · Validate inputs" in line
        for line in text.splitlines()
    )


def test_nested_progress_uses_deterministic_milestones_without_a_terminal():
    stream = StringIO()
    console = Console(
        file=stream,
        force_terminal=False,
        color_system=None,
        width=120,
    )
    outer = StageProgress("Module progress", enabled=True, console=console)
    inner = StageProgress("Scientific stages", enabled=True, console=console)

    outer.start_step(1, 1, "Run module")
    inner.start_step(1, 2, "Validate inputs")
    inner.complete_step(1, 2, "Validate inputs")
    inner.start_step(2, 2, "Calculate associations")
    inner.fail_step(2, 2, "Calculate associations")
    outer.fail_step(1, 1, "Run module")

    text = stream.getvalue()
    assert "None" not in text
    assert "━" not in text
    assert text.count("Started 1/1 · Run module") == 1
    assert text.count("Started 1/2 · Validate inputs") == 1
    assert text.count("Completed 1/2 · Validate inputs") == 1
    assert text.count("Started 2/2 · Calculate associations") == 1
    assert text.count("Failed 2/2 · Calculate associations") == 1
    assert text.count("Failed 1/1 · Run module") == 1


def test_failed_nested_progress_keeps_enclosing_live_region_paused():
    stream = StringIO()
    console = Console(
        file=stream,
        force_terminal=True,
        color_system=None,
        width=120,
    )
    outer = StageProgress("Module progress", enabled=True, console=console)
    inner = StageProgress("Scientific stages", enabled=True, console=console)
    outer.start_step(1, 1, "Run module")
    inner.start_step(1, 1, "Validate results")
    resume_calls = []
    outer._resume = lambda: resume_calls.append(True)

    inner.fail_step(1, 1, "Validate results")

    assert resume_calls == []
    outer.fail_step(1, 1, "Run module")


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


def test_checked_process_calls_progress_while_command_is_active():
    updates = []

    output = run_checked_command(
        [
            sys.executable,
            "-c",
            "import time; time.sleep(0.08); print('finished')",
        ],
        "Progress callback fixture",
        progress_callback=lambda: updates.append(len(updates)),
        progress_refresh_seconds=0.01,
    )

    assert output.strip() == "finished"
    assert len(updates) >= 3


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
