from io import StringIO
from pathlib import Path
import shutil
import sys
import time

import polars as pl
import pytest
from rich.console import Console
from rich.cells import cell_len
from rich.text import Text
from postgwas.config import load_configuration
from postgwas.config.models.logging import TerminalStyleConfig
from pydantic import ValidationError

from postgwas.core.io import tables as table_module
from postgwas.core.io.artifacts import publish_artifact_set
from postgwas.core.io.delimiters import pandas_csv_engine, resolve_delimiter
from postgwas.core.io.reports import write_html_report
from postgwas.core.io.tables import read_delimited_table, read_parquet_table
from postgwas.core.paths import configured_output_matches
from postgwas.core.preflight import (
    PreflightLogBuffer,
    capture_preflight_file_identities,
    replay_preflight_log,
    require_unchanged_preflight_files,
)
from postgwas.core.input_validation import InputValidationSession
from postgwas.core.processes import (
    build_container_command,
    run_checked_command,
    run_checked_line_command,
    run_checked_pipeline,
)
from postgwas.core.resource_preparation import (
    JsonLinesRunLogger,
    ResourcePreparationError,
    resumable_download,
    sha256,
)
from postgwas.core.errors import MissingRequiredArgumentsError
from postgwas.core.required_arguments import (
    RequiredArgument,
    require_resolved_arguments,
)
from postgwas.core.dataframes import (
    chromosome_expression,
    collect_streaming,
    count_matching_rows,
    numeric_column,
    position_expression,
)
from postgwas.core.execution.runtime import safe_thread_count
from postgwas.core.values import (
    format_fraction_percentage,
    format_percentage,
    optional_text,
    parse_integer,
)
from postgwas.core.ui import (
    MeasuredProgress,
    StageProgress,
    print_screen_block,
    run_with_progress,
)
from postgwas.core.ui.progress import _OperationColumn
from postgwas.core.ui.screen import (
    SYMBOLS,
    render_action_plan,
    screen_field,
    screen_line,
    screen_section_lines,
    screen_summary_card,
    style_screen_block,
    terminal_presentation,
)
from postgwas.core.vcf import (
    count_indexed_vcf_records,
    extract_vcf_table,
    select_vcf_sample,
    validate_indexed_vcf,
)
from postgwas.modules.harmonisation.shared.runtime import reject_rows


def test_html_report_writer_adds_shared_column_sorting_and_filtering(tmp_path):
    destination = tmp_path / "report.html"
    write_html_report(
        """<!doctype html><html><head><title>Report</title></head><body>
        <table><thead><tr><th>Gene</th><th>P</th></tr></thead>
        <tbody><tr><td>GENE2</td><td>0.2</td></tr>
        <tr><td>GENE1</td><td>1e-8</td></tr></tbody></table>
        <table><thead><tr><th>Empty column</th></tr></thead><tbody></tbody></table>
        <table class="key-value"><tbody><tr><th>Variants</th><td>2</td></tr></tbody></table>
        </body></html>""",
        destination,
    )

    report = destination.read_text(encoding="utf-8")
    assert report.count('id="postgwas-table-interactivity"') == 1
    assert 'className = \'postgwas-column-filter\'' in report
    assert "stableSort(filtered, column, ascending)" in report
    assert "header.cloneNode(true)" in report
    assert "const columnCount = header && header.cells.length" in report
    assert "['Field', 'Value']" in report
    assert report.index("postgwas-table-interactivity-style") < report.index(
        "</head>"
    )
    assert report.index('id="postgwas-table-interactivity"') < report.index(
        "</body>"
    )


def test_html_report_writer_does_not_duplicate_shared_interactivity(tmp_path):
    first = tmp_path / "first.html"
    second = tmp_path / "second.html"
    write_html_report(
        "<html><body><table><tr><td>value</td></tr></table></body></html>",
        first,
    )
    write_html_report(first.read_text(encoding="utf-8"), second)

    assert second.read_text(encoding="utf-8").count(
        'id="postgwas-table-interactivity"'
    ) == 1


def test_shared_preflight_file_identity_rejects_a_changed_resource(tmp_path):
    resource = tmp_path / "reference.txt"
    resource.write_text("first\n", encoding="utf-8")
    identities = capture_preflight_file_identities([resource])
    resource.write_text("changed resource\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="changed after pipeline preflight"):
        require_unchanged_preflight_files(identities)


def test_shared_preflight_log_buffer_replays_in_order():
    class Logger:
        def __init__(self):
            self.events = []

        def record(self, *args, **kwargs):
            self.events.append(("record", args, kwargs))

        def info(self, *args, **kwargs):
            self.events.append(("info", args, kwargs))

        def warn(self, *args, **kwargs):
            self.events.append(("warn", args, kwargs))

    buffer = PreflightLogBuffer()
    buffer.record("VALIDATE", "resource", status="passed")
    buffer.info("version output")
    buffer.warning("scientific warning", indent=1)
    logger = Logger()

    replay_preflight_log(buffer.events, logger)

    assert logger.events == [
        ("record", ("VALIDATE", "resource"), {"status": "passed"}),
        ("info", ("version output",), {"indent": 0}),
        ("warn", ("scientific warning",), {"indent": 1}),
    ]


def test_indexed_vcf_record_count_is_validated(monkeypatch):
    monkeypatch.setattr(
        "postgwas.core.vcf.run_checked_command",
        lambda *args, **kwargs: "42\n",
    )
    assert count_indexed_vcf_records("study.vcf.gz", "bcftools") == 42


def test_indexed_vcf_record_count_rejects_non_numeric_output(monkeypatch):
    monkeypatch.setattr(
        "postgwas.core.vcf.run_checked_command",
        lambda *args, **kwargs: "unknown\n",
    )
    with pytest.raises(RuntimeError, match="invalid VCF record count"):
        count_indexed_vcf_records(
            "study.vcf.gz", "bcftools", error_type=RuntimeError,
        )


def test_vcf_sample_contract_accepts_one_matching_sample(monkeypatch):
    monkeypatch.setattr(
        "postgwas.core.vcf._read_vcf_validation_command",
        lambda *_args, **_kwargs: "STUDY\n",
    )

    with InputValidationSession(reuse_checks=False) as session:
        sample = select_vcf_sample(
            "study.vcf.gz",
            "STUDY",
            "bcftools",
            error_type=RuntimeError,
        )

    assert sample == "STUDY"
    record = next(item for item in session.records if item.role == "VCF sample")
    assert record.status == "passed"
    assert record.checks == (
        "exactly one sample column",
        "dataset identifier match",
    )
    assert record.metrics == {
        "sample_count": 1,
        "sample": "STUDY",
        "dataset_id": "STUDY",
    }


@pytest.mark.parametrize(
    ("listed_samples", "expected_count"),
    (("", 0), ("STUDY\nOTHER\n", 2)),
)
def test_vcf_sample_contract_rejects_zero_or_multiple_samples(
    monkeypatch,
    listed_samples,
    expected_count,
):
    monkeypatch.setattr(
        "postgwas.core.vcf._read_vcf_validation_command",
        lambda *_args, **_kwargs: listed_samples,
    )

    with InputValidationSession(reuse_checks=False) as session:
        with pytest.raises(
            RuntimeError,
            match=(
                r"must contain exactly one sample column; bcftools reported %d"
                % expected_count
            ),
        ):
            select_vcf_sample(
                "study.vcf.gz",
                "STUDY",
                "bcftools",
                error_type=RuntimeError,
            )

    record = next(item for item in session.records if item.role == "VCF sample")
    assert record.status == "failed"
    assert record.metrics == {"sample_count": expected_count}


def test_vcf_sample_contract_warns_for_one_mismatching_sample(monkeypatch):
    monkeypatch.setattr(
        "postgwas.core.vcf._read_vcf_validation_command",
        lambda *_args, **_kwargs: "VCF_STUDY\n",
    )

    class Logger:
        def __init__(self):
            self.messages = []

        def warning(self, message):
            self.messages.append(message)

    logger = Logger()
    with InputValidationSession(reuse_checks=False) as session:
        sample = select_vcf_sample(
            "study.vcf.gz",
            "RUN_STUDY",
            "bcftools",
            logger=logger,
            error_type=RuntimeError,
        )

    assert sample == "VCF_STUDY"
    assert len(logger.messages) == 1
    assert "VCF sample 'VCF_STUDY' differs from dataset ID 'RUN_STUDY'" in (
        logger.messages[0]
    )
    record = next(item for item in session.records if item.role == "VCF sample")
    assert record.status == "warning"
    assert record.metrics == {
        "sample_count": 1,
        "sample": "VCF_STUDY",
        "dataset_id": "RUN_STUDY",
    }


def test_real_bcftools_sample_contract_rejects_requested_sample_in_multisample_vcf(
    tmp_path,
):
    bcftools = shutil.which("bcftools")
    if bcftools is None:
        pytest.skip("bcftools is required for integration test")
    vcf = tmp_path / "multisample.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1,length=1000>\n"
        "##FORMAT=<ID=ES,Number=1,Type=Float,Description=\"Effect size\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTUDY\tOTHER\n"
        "1\t101\trs1\tA\tG\t.\tPASS\t.\tES\t0.1\t0.2\n",
        encoding="utf-8",
    )

    with pytest.raises(
        RuntimeError,
        match=r"must contain exactly one sample column; bcftools reported 2",
    ):
        select_vcf_sample(
            vcf,
            "STUDY",
            bcftools,
            error_type=RuntimeError,
        )


def test_indexed_vcf_validation_reuses_only_unchanged_in_run_evidence(
    monkeypatch, tmp_path,
):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"first")
    Path(str(vcf) + ".tbi").write_bytes(b"index")
    header = (
        "##fileformat=VCFv4.2\n"
        "##genome_build=GRCh37\n"
        "##contig=<ID=1>\n"
        "##FORMAT=<ID=ES,Number=1,Type=Float,Description=\"effect\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTUDY\n"
    )
    calls = {"header": 0, "count": 0, "sample": 0}

    def read_header(*_args, **_kwargs):
        calls["header"] += 1
        return header

    def count_records(*_args, **_kwargs):
        calls["count"] += 1
        return 8

    def select_sample(*_args, **_kwargs):
        calls["sample"] += 1
        return "STUDY"

    monkeypatch.setattr("postgwas.core.vcf.read_vcf_header", read_header)
    monkeypatch.setattr(
        "postgwas.core.vcf.count_indexed_vcf_records", count_records,
    )
    monkeypatch.setattr("postgwas.core.vcf.select_vcf_sample", select_sample)

    first = validate_indexed_vcf(
        vcf,
        "STUDY",
        "bcftools",
        genome_build_header="##genome_build={build}",
        supported_genome_builds=["GRCh37", "GRCh38"],
        required_fields=["FORMAT/ES"],
        expected_genome_build="GRCh37",
        error_type=RuntimeError,
    )
    reused = validate_indexed_vcf(
        vcf,
        "STUDY",
        "bcftools",
        genome_build_header="##genome_build={build}",
        supported_genome_builds=["GRCh37", "GRCh38"],
        required_fields=["FORMAT/ES"],
        expected_genome_build="GRCh37",
        cached=first,
        error_type=RuntimeError,
    )

    assert reused.reused is True
    assert calls == {"header": 1, "count": 1, "sample": 1}

    vcf.write_bytes(b"changed-size")
    refreshed = validate_indexed_vcf(
        vcf,
        "STUDY",
        "bcftools",
        genome_build_header="##genome_build={build}",
        supported_genome_builds=["GRCh37", "GRCh38"],
        required_fields=["FORMAT/ES"],
        cached=reused,
        error_type=RuntimeError,
    )

    assert refreshed.reused is False
    assert calls == {"header": 2, "count": 2, "sample": 2}


def test_vcf_extraction_reuses_a_validated_sample_without_listing_samples(
    monkeypatch, tmp_path,
):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"vcf")
    destination = tmp_path / "projection.tsv"
    commands = []

    monkeypatch.setattr(
        "postgwas.core.vcf.select_vcf_sample",
        lambda *_args, **_kwargs: pytest.fail(
            "a preflight-validated sample must not be queried again"
        ),
    )

    def run(command, *_args, **kwargs):
        commands.append(command)
        kwargs["stdout_path"].write_text("1\t100\n", encoding="utf-8")
        return ""

    monkeypatch.setattr("postgwas.core.vcf.run_checked_command", run)
    extract_vcf_table(
        vcf,
        destination,
        "STUDY",
        {"CHROM": "%CHROM", "POS": "%POS"},
        "bcftools",
        delimiter="\t",
        io_buffer_bytes=128,
        validated_sample="VCF_STUDY",
    )

    assert destination.read_text(encoding="utf-8") == "CHROM\tPOS\n1\t100\n"
    sample_index = commands[0].index("--samples")
    assert commands[0][sample_index + 1] == "VCF_STUDY"


def test_required_arguments_support_cli_only_values_and_deduplicate_messages():
    requirements = (
        RequiredArgument("--vcf", None, None),
        RequiredArgument("--vcf", None, ""),
        RequiredArgument("--reference", "modules.example.reference", None),
    )

    with pytest.raises(MissingRequiredArgumentsError) as captured:
        require_resolved_arguments(requirements)

    assert str(captured.value).splitlines() == [
        "Required argument not provided: --vcf. Provide --vcf VALUE.",
        "Required argument not provided: --reference. Provide --reference VALUE "
        "or set modules.example.reference in the run configuration.",
    ]


def test_publish_artifact_set_restores_existing_outputs_on_failure(
    tmp_path, monkeypatch,
):
    first_temporary = tmp_path / "first.tmp"
    second_temporary = tmp_path / "second.tmp"
    first_final = tmp_path / "first.txt"
    second_final = tmp_path / "second.txt"
    first_temporary.write_text("new first", encoding="utf-8")
    second_temporary.write_text("new second", encoding="utf-8")
    first_final.write_text("old first", encoding="utf-8")
    second_final.write_text("old second", encoding="utf-8")

    original_replace = type(first_temporary).replace

    def fail_second_publish(path, target):
        if path == second_temporary:
            raise OSError("simulated publication failure")
        return original_replace(path, target)

    monkeypatch.setattr(type(first_temporary), "replace", fail_second_publish)

    with pytest.raises(RuntimeError, match="Existing outputs were restored"):
        publish_artifact_set((
            (first_temporary, first_final),
            (second_temporary, second_final),
        ))

    assert first_final.read_text(encoding="utf-8") == "old first"
    assert second_final.read_text(encoding="utf-8") == "old second"
    assert second_temporary.read_text(encoding="utf-8") == "new second"


def test_optional_values_have_one_normalization_rule():
    assert optional_text("  beta  ") == "beta"
    assert optional_text("NaN") is None
    assert optional_text(None) is None
    assert parse_integer("38,691") == 38691
    assert parse_integer("not-a-count") is None


@pytest.mark.parametrize(
    ("separator", "expected"),
    [
        (r"\s+", "c"),
        ("\t", "c"),
        (",", "c"),
        (r"\t+", "python"),
        (r"[|;]", "python"),
    ],
)
def test_pandas_csv_engine_preserves_fast_supported_separators_and_regexes(
    separator, expected,
):
    assert pandas_csv_engine(separator) == expected


def test_pandas_csv_engine_rejects_an_empty_separator():
    with pytest.raises(ValueError, match="separator must not be empty"):
        pandas_csv_engine("")


def test_safe_thread_count_respects_an_explicit_memory_budget():
    assert safe_thread_count(
        16,
        gb_per_thread=20,
        available_ram_gb=64,
        reporter=None,
    ) == 3
    assert safe_thread_count(
        2,
        gb_per_thread=20,
        available_ram_gb=64,
        reporter=None,
    ) == 2


def test_safe_thread_count_uses_scheduler_aware_detection_without_a_budget(
    monkeypatch,
):
    monkeypatch.setattr(
        "postgwas.core.execution.runtime.auto_detect_ram_gb",
        lambda: 45,
    )

    assert safe_thread_count(8, gb_per_thread=20, reporter=None) == 2


def test_streaming_collection_uses_only_a_declared_polars_keyword():
    class LegacyFrame:
        def __init__(self):
            self.streaming = False

        def collect(self, *, streaming=False, **_kwargs):
            self.streaming = streaming
            return "legacy"

    class CurrentFrame:
        def __init__(self):
            self.engine = "auto"

        def collect(self, *, engine="auto"):
            self.engine = engine
            return "current"

    class SilentFallbackFrame:
        def collect(self, **_kwargs):
            raise AssertionError("an unsupported keyword must not be attempted")

    legacy = LegacyFrame()
    current = CurrentFrame()
    assert collect_streaming(legacy) == "legacy"
    assert legacy.streaming is True
    assert collect_streaming(current) == "current"
    assert current.engine == "streaming"
    with pytest.raises(RuntimeError, match="does not expose"):
        collect_streaming(SilentFallbackFrame())


def test_screen_fields_align_unicode_symbols_and_wrap_long_values():
    assert len(SYMBOLS["warning"]) == 1
    assert cell_len(SYMBOLS["warning"]) == cell_len(SYMBOLS["count"]) == 2
    count = screen_field(
        "count", "Count label", "10", width=80, indent=10, label_width=42,
    )
    warning = screen_field(
        "warning", "Warning label", "1", width=80, indent=10, label_width=42,
    )
    assert cell_len(count.split(":", 1)[0]) == cell_len(
        warning.split(":", 1)[0]
    )
    assert len(count.split(":", 1)[0]) == len(warning.split(":", 1)[0])

    nested = screen_field(
        "success",
        "Nested label",
        "valid",
        width=80,
        indent=14,
        separator_column=cell_len(count.split(" : ", 1)[0]),
    )
    assert cell_len(nested.split(" : ", 1)[0]) == cell_len(
        count.split(" : ", 1)[0]
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

    path = (
        "/Users/example/analysis/results/"
        "STUDY_formatter_report_with_a_long_name.html"
    )
    path_rendered = screen_field(
        "success", "Detailed HTML report", path,
        width=80, indent=10, label_width=42, path_value=True,
    )
    path_lines = path_rendered.splitlines()
    assert len(path_lines) > 1
    assert all(cell_len(line) <= 80 for line in path_lines)
    assert "".join([
        path_lines[0].split(" : ", 1)[1],
        *(line[59:] for line in path_lines[1:]),
    ]) == path
    assert path_lines[0].split(" : ", 1)[1]
    assert all(line.startswith(" " * 59) for line in path_lines[1:])

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


def test_summary_cards_keep_shared_long_labels_on_one_aligned_row():
    lines = screen_summary_card(
        "analysis",
        "Statistical quality",
        [
            ("analysis", "Effective sample size", "12,519,186 usable"),
            ("analysis", "Upper-tail diagnostic", "0 informational matches"),
            (
                "warning",
                "Reference unmatched retained",
                "1 variant retained by policy",
            ),
        ],
    )

    fields = [line for line in lines if " : " in line]
    assert len(fields) == 3
    assert all(label in "\n".join(fields) for label in (
        "Effective sample size",
        "Upper-tail diagnostic",
        "Reference unmatched retained",
    ))
    assert len({cell_len(line.split(" : ", 1)[0]) for line in fields}) == 1


def test_shared_result_section_matches_filtering_spacing():
    assert screen_section_lines("analysis", "QC results by rule") == [
        "",
        "",
        "      🔬  QC results by rule",
        "",
    ]


def test_shared_action_plan_preserves_filtering_order_and_footer():
    rendered = render_action_plan(
        "QC assessment plan · conditions evaluated independently",
        [
            {"action": "exclude", "label": "Missing FORMAT/AF"},
            {"action": "keep", "label": "Missing FORMAT/SI"},
        ],
        empty_message="No active conditions",
        footer=(("info", "The input VCF remains unchanged"),),
    )

    assert rendered.splitlines() == [
        "",
        "    🔬  QC assessment plan · conditions evaluated independently",
        "",
        "        1. EXCLUDE · Missing FORMAT/AF",
        "        2. KEEP · Missing FORMAT/SI",
        "",
        "        🔹  The input VCF remains unchanged",
    ]


def test_semantic_screen_styles_preserve_plain_text_and_keep_values_neutral():
    block = "\n".join((
        screen_line("analysis", "Filtering summary", indent=4),
        screen_field(
            "info", "Dataset", "SCZ_2026", width=120,
            indent=8, label_width=24,
        ),
        screen_field(
            "decision", "Action", "Rerun without --samp-prev", width=120,
            indent=8, label_width=24,
        ),
        screen_field(
            "attention", "Important", "LDSC used the provided value", width=120,
            indent=8, label_width=24,
        ),
        "",
        screen_line("loss", "Variants removed by MAF", indent=8),
        screen_line("success", "Removal accounting passed", indent=8),
    ))

    styled = style_screen_block(block)

    assert styled.plain == block
    styled_fragments = [
        (styled.plain[span.start:span.end], str(span.style))
        for span in styled.spans
    ]
    assert any(
        "Filtering summary" in fragment and style == "bold cyan"
        for fragment, style in styled_fragments
    )
    assert any(
        "Dataset" in fragment and style == "default"
        for fragment, style in styled_fragments
    )
    assert not any(
        "SCZ_2026" in fragment for fragment, _style in styled_fragments
    )
    assert any(
        "Action" in fragment and style == "bold magenta"
        for fragment, style in styled_fragments
    )
    assert not any(
        "Rerun without --samp-prev" in fragment
        for fragment, _style in styled_fragments
    )
    assert any(
        "Important" in fragment and style == "bold red"
        for fragment, style in styled_fragments
    )
    assert any(
        "LDSC used the provided value" in fragment
        and style == "bold red"
        for fragment, style in styled_fragments
    )
    assert any(
        "Variants removed by MAF" in fragment and style == "bold yellow"
        for fragment, style in styled_fragments
    )
    assert any(
        "Removal accounting passed" in fragment and style == "bold green"
        for fragment, style in styled_fragments
    )


def test_semantic_screen_styles_take_heading_style_from_run_yaml():
    block = "\n".join((
        screen_line("analysis", "QC results by rule", indent=6),
        screen_field(
            "info", "Input variants", "12,519,186",
            width=120, indent=8, label_width=24,
        ),
    ))

    logging = load_configuration(cli_overrides={
        "logging.terminal_style.styles.analysis": "bold #d8b4fe",
    }).logging
    with terminal_presentation(logging):
        styled = style_screen_block(block)
    styled_fragments = [
        (styled.plain[span.start:span.end], str(span.style))
        for span in styled.spans
    ]

    assert styled.plain == block
    assert any(
        "QC results by rule" in fragment and style == "bold #d8b4fe"
        for fragment, style in styled_fragments
    )
    assert not any(
        "12,519,186" in fragment
        for fragment, _style in styled_fragments
    )


def test_semantic_screen_styles_allow_scoped_neutral_yaml_and_restore_default():
    block = "\n".join((
        screen_line("info", "Neutral information", indent=8),
        screen_line("success", "Validated result", indent=8),
    ))

    logging = load_configuration(cli_overrides={
        "logging.terminal_style.styles.info": "default",
    }).logging
    with terminal_presentation(logging):
        styled = style_screen_block(block)
    styled_fragments = [
        (styled.plain[span.start:span.end], str(span.style))
        for span in styled.spans
    ]

    assert styled.plain == block
    assert not any(
        "Neutral information" in fragment and style != "default"
        for fragment, style in styled_fragments
    )
    assert any(
        "Validated result" in fragment and style == "bold green"
        for fragment, style in styled_fragments
    )


def test_semantic_screen_styles_reject_unknown_yaml_kinds():
    values = load_configuration().logging.terminal_style.model_dump()
    values["styles"]["unknown_kind"] = "red"
    with pytest.raises(ValidationError, match="unknown_kind"):
        TerminalStyleConfig.model_validate(values)


def test_print_screen_block_emits_attention_color_to_a_terminal():
    stream = StringIO()
    console = Console(
        file=stream,
        force_terminal=True,
        color_system="standard",
        width=120,
    )
    block = screen_field(
        "attention",
        "Important — value selected",
        "LDSC USED the provided prevalence",
        width=120,
        indent=10,
        label_width=42,
    )

    print_screen_block(block, console=console)

    terminal_output = stream.getvalue()
    assert "\x1b[" in terminal_output
    assert Text.from_ansi(terminal_output).plain.rstrip("\n") == block


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
        outcome_fields=[
            ("analysis", "Compatibility assessment"),
            ("count", "Reference identifiers", 8),
        ],
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
    top_level_outcome_lines = [
        line
        for line in text.splitlines()
        if "Records read" in line
        or "Records retained" in line
    ]
    assert len({line.index(" : ") for line in top_level_outcome_lines}) == 1
    nested_line = next(
        line for line in text.splitlines() if "Reference identifiers" in line
    )
    assert cell_len(nested_line.split(" : ", 1)[0]) == cell_len(
        top_level_outcome_lines[0].split(" : ", 1)[0]
    )
    assert "🔬  Compatibility assessment" in text
    assert "Failed 3/3 · Validate results" in text
    assert "2/3" in text


def test_stage_progress_live_row_distinguishes_active_from_completed_stage():
    from rich.progress import Task

    task = Task(
        0,
        "Current 5/6 · Fit models and calculate scores",
        total=6,
        completed=4,
        _get_time=time.monotonic,
        fields={"active_stage": 5, "active_stage_total": 6},
    )

    rendered = _OperationColumn().render(task).plain
    assert "5/6 · Fit models and calculate scores" in rendered


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

    with progress.step(1, 2, "Prepare inputs") as step:
        step.outcome(
            "Inputs are ready.",
            fields=(
                ("count", "Variants", 10),
                ("success", "Reference", "validated"),
            ),
        )
    with pytest.raises(ValueError, match="invalid result"):
        with progress.step(2, 2, "Validate results"):
            raise ValueError("invalid result")

    text = stream.getvalue()
    assert "Completed 1/2 · Prepare inputs" in text
    assert "Outcome" in text
    assert "Variants" in text
    assert "Reference" in text
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


def test_shared_rejection_respects_an_explicit_domain_null_policy():
    frame = pl.DataFrame({"value": [1, None, -1]})
    kept, removed = reject_rows(
        frame,
        (pl.col("value") > 0).fill_null(False),
        step_label="test",
        reason="invalid_position",
    )
    assert removed == 1
    assert kept["value"].to_list() == [None, -1]


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


def test_checked_process_reports_peak_process_tree_memory():
    resources = {}

    run_checked_command(
        [
            sys.executable,
            "-c",
            "import time; payload = bytearray(8 * 1024 * 1024); "
            "time.sleep(0.3); assert payload",
        ],
        "Resource measurement fixture",
        resource_metrics=resources,
        resource_poll_seconds=0.005,
    )

    assert resources["memory_samples"] >= 1
    assert resources["peak_rss_bytes"] >= 8 * 1024 * 1024
    assert (
        resources["process_tree_samples"] + resources["root_only_samples"]
        == resources["memory_samples"]
    )


def test_checked_process_retains_root_memory_when_tree_enumeration_is_denied(
    monkeypatch,
):
    import psutil

    def deny_descendant_enumeration(_process, recursive=False):
        del recursive
        raise psutil.AccessDenied(pid=_process.pid)

    monkeypatch.setattr(psutil.Process, "children", deny_descendant_enumeration)
    resources = {}

    run_checked_command(
        [
            sys.executable,
            "-c",
            "import time; payload = bytearray(8 * 1024 * 1024); "
            "time.sleep(0.3); assert payload",
        ],
        "Restricted resource measurement fixture",
        resource_metrics=resources,
        resource_poll_seconds=0.005,
    )

    assert resources["memory_samples"] >= 1
    assert resources["root_only_samples"] == resources["memory_samples"]
    assert resources["peak_rss_bytes"] >= 8 * 1024 * 1024


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


def test_checked_line_command_consumes_large_outputs_incrementally():
    lines = []
    line_count = run_checked_line_command(
        [
            sys.executable,
            "-c",
            "import sys; "
            "[print('row-%d' % value) for value in range(1000)]; "
            "print('checked stderr', file=sys.stderr)",
        ],
        "Streaming line fixture",
        line_consumer=lambda line: lines.append(line.rstrip("\n")),
    )

    assert line_count == 1000
    assert lines[0] == "row-0"
    assert lines[-1] == "row-999"


def test_checked_line_command_reports_failure_and_terminates_on_consumer_error():
    with pytest.raises(RuntimeError, match="exit status 7.*stream failure"):
        run_checked_line_command(
            [
                sys.executable,
                "-c",
                "import sys; print('stream failure', file=sys.stderr); "
                "raise SystemExit(7)",
            ],
            "Failing line fixture",
            line_consumer=lambda _line: None,
        )

    started = time.monotonic()
    with pytest.raises(ValueError, match="consumer rejected row"):
        run_checked_line_command(
            [
                sys.executable,
                "-c",
                "import time; print('row', flush=True); time.sleep(10)",
            ],
            "Rejected line fixture",
            line_consumer=lambda _line: (_ for _ in ()).throw(
                ValueError("consumer rejected row")
            ),
        )
    assert time.monotonic() - started < 5


def test_checked_pipeline_preserves_binary_streams_between_stages():
    results = run_checked_pipeline(
        [
            [
                sys.executable,
                "-c",
                "import sys; sys.stdout.buffer.write(b'a\\x00b')",
            ],
            [
                sys.executable,
                "-c",
                "import sys; data=sys.stdin.buffer.read(); "
                "sys.stdout.buffer.write(data.upper())",
            ],
        ],
        "Binary pipeline fixture",
    )

    assert len(results) == 2
    assert all(result.returncode == 0 for result in results)
    assert results[-1].stdout == "A\x00B"


def test_checked_pipeline_drains_stderr_and_terminates_other_stages():
    started = time.monotonic()
    with pytest.raises(RuntimeError) as captured:
        run_checked_pipeline(
            [
                [sys.executable, "-c", "import time; time.sleep(10)"],
                [
                    sys.executable,
                    "-c",
                    "import sys; sys.stderr.write('root-cause-' * 25000); "
                    "raise SystemExit(7)",
                ],
                [sys.executable, "-c", "import time; time.sleep(10)"],
            ],
            "Failing pipeline fixture",
        )

    assert time.monotonic() - started < 5
    message = str(captured.value)
    assert "stage 2 failed with exit status 7" in message
    assert "root-cause" in message


def test_checked_pipeline_requires_declared_outputs(tmp_path):
    expected = tmp_path / "missing.vcf.gz"
    with pytest.raises(RuntimeError, match="expected output is missing or empty"):
        run_checked_pipeline(
            [[sys.executable, "-c", "pass"]],
            "Missing pipeline output fixture",
            expected_outputs=(expected,),
        )


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


def test_shared_delimited_reader_pushes_projection_into_polars(
    monkeypatch, tmp_path,
):
    table = tmp_path / "wide_reference.csv"
    table.write_text(
        "marker,score,unused\nrs1,0.8,large\nrs2,0.9,payload\n",
        encoding="utf-8",
    )
    observed = {}
    original = table_module.pl.read_csv

    def tracked_read_csv(*args, **kwargs):
        observed["columns"] = kwargs.get("columns")
        return original(*args, **kwargs)

    monkeypatch.setattr(table_module.pl, "read_csv", tracked_read_csv)
    frame, _detected = read_delimited_table(
        table,
        "comma",
        candidates=["comma"],
        minimum_columns=2,
        maximum_columns=10,
        sample_lines=2,
        null_values=["NA"],
        infer_schema_length=10,
        columns=["score", "marker"],
    )

    assert observed["columns"] == ["score", "marker"]
    assert frame.columns == ["score", "marker"]
    assert "unused" not in frame.columns


def test_shared_delimited_reader_rejects_absent_projection_before_data_read(
    monkeypatch, tmp_path,
):
    table = tmp_path / "wide_reference.csv"
    table.write_text(
        "marker,score,unused\nrs1,0.8,large\n",
        encoding="utf-8",
    )

    def unexpected_read(*_args, **_kwargs):
        pytest.fail("data reader must not run after header validation fails")

    monkeypatch.setattr(table_module.pl, "read_csv", unexpected_read)
    with pytest.raises(
        RuntimeError,
        match=(
            r"missing configured column\(s\): population\. "
            r"Available columns: marker, score, unused"
        ),
    ):
        read_delimited_table(
            table,
            "comma",
            candidates=["comma"],
            minimum_columns=2,
            maximum_columns=10,
            sample_lines=2,
            null_values=["NA"],
            infer_schema_length=10,
            columns=["marker", "population"],
        )


def test_shared_whitespace_reader_pushes_projection_into_pandas(
    monkeypatch, tmp_path,
):
    import pandas as pd

    table = tmp_path / "wide_reference.txt"
    table.write_text(
        "marker score unused\nrs1 0.8 large\nrs2 0.9 payload\n",
        encoding="utf-8",
    )
    observed = {}
    original = pd.read_csv

    def tracked_read_csv(*args, **kwargs):
        observed["usecols"] = kwargs.get("usecols")
        return original(*args, **kwargs)

    monkeypatch.setattr(pd, "read_csv", tracked_read_csv)
    frame, _detected = read_delimited_table(
        table,
        "whitespace",
        candidates=["whitespace"],
        minimum_columns=2,
        maximum_columns=10,
        sample_lines=2,
        null_values=["NA"],
        infer_schema_length=10,
        columns=["score", "marker"],
    )

    assert observed["usecols"] == ["score", "marker"]
    assert frame.columns == ["score", "marker"]


def test_shared_parquet_reader_projects_after_metadata_validation(
    monkeypatch, tmp_path,
):
    table = tmp_path / "wide_reference.parquet"
    pl.DataFrame({
        "marker": ["rs1", "rs2"],
        "score": [0.8, 0.9],
        "unused": ["large", "payload"],
    }).write_parquet(table)
    observed = {}
    original = table_module.pl.read_parquet

    def tracked_read_parquet(*args, **kwargs):
        observed["columns"] = kwargs.get("columns")
        return original(*args, **kwargs)

    monkeypatch.setattr(table_module.pl, "read_parquet", tracked_read_parquet)
    frame = read_parquet_table(table, columns=["score", "marker"])

    assert observed["columns"] == ["score", "marker"]
    assert frame.columns == ["score", "marker"]


def test_shared_parquet_reader_rejects_absent_projection_before_data_read(
    monkeypatch, tmp_path,
):
    table = tmp_path / "wide_reference.parquet"
    pl.DataFrame({"marker": ["rs1"], "score": [0.8]}).write_parquet(table)

    def unexpected_read(*_args, **_kwargs):
        pytest.fail("data reader must not run after metadata validation fails")

    monkeypatch.setattr(table_module.pl, "read_parquet", unexpected_read)
    with pytest.raises(
        RuntimeError,
        match=r"missing configured column\(s\): population",
    ):
        read_parquet_table(table, columns=["marker", "population"])
