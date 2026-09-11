"""Presentation contracts; no GWAS analysis or network services are executed."""

from argparse import Namespace
import importlib
from io import StringIO
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

from pydantic import ValidationError
import pytest
from rich.cells import cell_len
from rich.console import Console
from rich.text import Text

from postgwas.config import load_configuration
from postgwas.config.models.logging import TerminalStyleConfig
from postgwas.core.ui.terminal_text import TerminalTextDecoder
from postgwas.core.ui import print_screen_message, StageProgress
from postgwas.core.ui.screen import (
    SYMBOLS, screen_field, screen_line, style_screen_block,
    terminal_presentation, terminal_settings, terminal_style,
)
from postgwas.pipeline.registry import REGISTRY


def test_theme_has_exactly_the_shared_semantic_roles():
    theme = load_configuration().logging.terminal_style
    assert set(theme.styles) == set(SYMBOLS) | {"text"}
    assert theme.color == "auto"


@pytest.mark.parametrize("change", ("missing", "unknown", "bad_style"))
def test_invalid_theme_fails_before_analysis(change):
    theme = load_configuration().logging.terminal_style.model_dump()
    if change == "missing":
        theme["styles"].pop("warning")
    elif change == "unknown":
        theme["styles"]["unexpected"] = "red"
    else:
        theme["styles"]["warning"] = "not-a-rich-colour"
    with pytest.raises(ValidationError):
        TerminalStyleConfig.model_validate(theme)


def test_yaml_theme_is_scoped_and_restored_even_after_failure(tmp_path):
    path = tmp_path / "run.yaml"
    path.write_text("logging:\n  terminal_style:\n    styles:\n      error: bold blue\n")
    original = terminal_settings()
    configured = load_configuration(path).logging
    with pytest.raises(RuntimeError):
        with terminal_presentation(configured):
            assert terminal_style("error") == "bold blue"
            assert terminal_style("success") == "bold green"
            raise RuntimeError("test")
    assert terminal_settings() is original


@pytest.mark.parametrize("width", (60, 80, 120, 180))
def test_fields_use_one_value_column_and_wrap_without_losing_paths(width):
    logging = load_configuration().logging
    path = "/reference/" + "long-segment/" * 9 + "study[1].vcf.gz"
    with terminal_presentation(logging):
        fields = [screen_field("info", label, path, width=width, indent=indent,
                               label_width=local_width, path_value=True)
                  for label, indent, local_width in (
                      ("Input", 6, 20), ("Output directory", 8, 28),
                      ("Very long nested annotation reference file", 18, 42),
                  )]
    assert len({cell_len(field.split(" : ")[0].splitlines()[-1]) for field in fields}) == 1
    for field in fields:
        assert max(map(cell_len, field.splitlines())) <= width
        value = field.split(" : ", 1)[1]
        assert "".join(line.strip() for line in value.splitlines()) == path


def test_semantic_styles_never_interpret_markup_or_colour_unrelated_text():
    value = "  ❌  Failure [blue]literal[/blue]\n      More detail\nnative continuation\n"
    styled = style_screen_block(Text(value, style="white on magenta"))
    assert styled.plain == value
    assert not any("native continuation" in value[s.start:s.end] for s in styled.spans)
    assert not any(str(s.style) == "white on magenta" for s in styled.spans)
    assert any("More detail" in value[s.start:s.end] and str(s.style) == "bold red"
               for s in styled.spans)


@pytest.mark.parametrize("role", ("warning", "error", "success", "analysis"))
def test_messages_wrap_with_one_prefix_and_literal_square_brackets(role):
    stream = StringIO()
    console = Console(file=stream, width=60, force_terminal=True, color_system="standard")
    message = "File [red]literal[/red] " + "description " * 8 + "\nSecond source line"
    print_screen_message(role, message, console=console)
    rendered = Text.from_ansi(stream.getvalue()).plain
    assert rendered.count(SYMBOLS[role]) == 1
    assert "[red]literal[/red]" in rendered
    assert " ".join(message.split()) == " ".join(rendered.replace(SYMBOLS[role], "", 1).split())
    assert all(cell_len(line) <= 60 for line in rendered.splitlines())


def test_never_colour_applies_even_to_an_explicit_forced_console():
    logging = load_configuration(cli_overrides={"logging.terminal_style.color": "never"}).logging
    stream = StringIO()
    console = Console(file=stream, force_terminal=True, color_system="standard", no_color=False)
    with terminal_presentation(logging):
        print_screen_message("warning", "Uncoloured warning", console=console)
    assert "\x1b" not in stream.getvalue()


@pytest.mark.parametrize("block_size", (1, 2, 3, 7, 64, 65536))
def test_decoder_preserves_utf8_and_diagnostics_across_every_boundary(block_size):
    source = ("\x1b[31m❌ failed\x1b[0m\r\n"
              "native [A/G] β=−0.25\rnext\n"
              "\x1b]8;;https://example.invalid\x1b\\literal-link\x1b]8;;\x1b\\\n"
              "\x1b]malformed control\nvisible diagnostic\nend").encode()
    decoder = TerminalTextDecoder("utf-8")
    actual = "".join(decoder.decode(source[i:i + block_size])
                     for i in range(0, len(source), block_size))
    actual += decoder.decode(b"", final=True)
    assert actual == "❌ failed\nnative [A/G] β=−0.25\nnext\nliteral-link\n\nvisible diagnostic\nend"


_CAPTURE_SCRIPT = r'''
from pathlib import Path
import os, sys
from unittest.mock import patch
from postgwas.config import load_configuration
from postgwas.core.screen_logging import ScreenSettings, record_screen
from postgwas.core.ui import StageProgress, print_screen_message
config = load_configuration(cli_overrides={"logging.terminal_style.color": sys.argv[3]})
with patch("postgwas.core.screen_logging.os.isatty", return_value=sys.argv[2] == "tty"):
    with record_screen(ScreenSettings(sys.argv[2] != "hidden", Path(sys.argv[1]), config)):
        progress = StageProgress("Shared run", enabled=True)
        progress.start_step(1, 2, "Read input")
        print_screen_message("warning", "Warning [literal]")
        os.write(1, b"\x1b[35mnative [A/G]\x1b[0m\r\n")
        progress.complete_step(1, 2, "Read input")
        progress.start_step(2, 2, "Validate result")
        progress.fail_step(2, 2, "Validate result")
        print_screen_message("error", "Result invalid [literal]", stderr=True)
        os.write(2, b"final without newline")
'''


@pytest.mark.parametrize("mode", ("tty", "redirected", "hidden"))
@pytest.mark.parametrize("colour", ("auto", "never", "no_color"))
def test_colours_do_not_leak_to_logs_or_hide_native_diagnostics(tmp_path, mode, colour):
    log = tmp_path / "screen.log"
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    if colour == "no_color":
        env["NO_COLOR"] = "1"
    else:
        env.pop("NO_COLOR", None)
    result = subprocess.run(
        [sys.executable, "-B", "-c", _CAPTURE_SCRIPT, str(log), mode,
         "auto" if colour == "no_color" else colour],
        capture_output=True, text=True, env=env, timeout=30,
    )
    assert result.returncode == 0, result.stderr
    saved = log.read_text()
    for message in ("Warning [literal]", "native [A/G]", "Failed 2/2",
                    "Result invalid [literal]", "final without newline"):
        assert saved.count(message) == 1
    assert "\x1b" not in saved and "\r" not in saved and "━" not in saved
    assert "All 2 stages completed" not in saved
    terminal = result.stdout + result.stderr
    if mode == "hidden":
        assert not terminal
    elif mode == "redirected":
        assert "\x1b" not in terminal
        assert "native [A/G]" in terminal
    else:
        assert "native [A/G]" in Text.from_ansi(terminal).plain
        # Cursor-control ANSI remains necessary for a live bar even without colour.
        assert ("\x1b[1;33m" in terminal) is (colour == "auto")
        assert "\x1b[35m" not in terminal  # native colour must not override the theme


_DIRECT_COMMANDS = tuple(name for name, spec in REGISTRY.commands().items()
                         if spec.direct_checkpoint != "not_applicable")


@pytest.mark.parametrize("command", _DIRECT_COMMANDS)
@pytest.mark.parametrize("status", (0, 1))
def test_every_public_direct_command_uses_the_shared_presentation_boundary(
    tmp_path, command, status,
):
    # Exercise the real dispatcher/recorder, replacing only analysis and audit IO.
    log = tmp_path / "screen.log"
    source = '''
from pathlib import Path
import sys
from unittest.mock import patch
from postgwas import __main__ as main
from postgwas.config import load_configuration
from postgwas.core.screen_logging import ScreenSettings
from postgwas.core.ui import print_screen_message
from postgwas.core.ui.screen import terminal_style
log, command, status = Path(sys.argv[1]), sys.argv[2], int(sys.argv[3])
config = load_configuration(cli_overrides={"logging.terminal_style.styles.warning": "bold blue"})
def operation():
    assert terminal_style("warning") == "bold blue"
    print_screen_message("warning", "Shared literal warning [file]")
    return status
sys.argv = ["postgwas", command, "--dataset-id", "STUDY"]
with patch.object(main, "resolve_screen_settings", return_value=ScreenSettings(True, log, config)):
    with patch.object(main, "resolve_reference", return_value=operation):
        sys.exit(main.main())
'''
    result = subprocess.run([sys.executable, "-B", "-c", source, str(log), command, str(status)],
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == status, result.stderr
    assert "Shared literal warning [file]" in result.stdout
    text = log.read_text()
    assert "Shared literal warning [file]" in text
    assert ("Failed 1/1" in text) is bool(status)
    assert ("All 1 stages completed" in text) is (status == 0)


@pytest.mark.parametrize("package,runner", (
    ("formatting", "run_formatter_direct"), ("single_cell", "run_single_cell_direct"),
    ("imputation", "run_sumstat_imputation_direct"), ("ldsc", "run_ldsc_direct"),
    ("magmacovar", "run_magma_covar_direct"), ("flames", "run_flames_direct"),
    ("magma", "run_magma_direct"), ("pops", "run_pops_direct"),
    ("mixer", "run_mixer_direct"), ("kpops", "run_kpops_direct"),
    ("caldera", "run_caldera_direct"), ("gcta_gene", "run_gcta_gene_direct"),
))
def test_module_error_adapters_preserve_literal_reason_and_nonzero_status(
    monkeypatch, capsys, package, runner,
):
    from postgwas.core.errors import ConfigurationError

    cli = importlib.import_module("postgwas.modules.%s.cli" % package)
    service = importlib.import_module("postgwas.modules.%s.service" % package)
    parser = SimpleNamespace(parse_args=lambda *_: Namespace(), print_help=lambda **_: None)
    monkeypatch.setattr(cli, "build_parser", lambda *_: parser)

    def fail(*_):
        raise ConfigurationError("Invalid [red]literal[/red] file")

    monkeypatch.setattr(service, runner, fail)
    if hasattr(cli, runner):
        monkeypatch.setattr(cli, runner, fail)
    assert cli.main(["--dataset-id", "STUDY"]) == 1
    captured = capsys.readouterr()
    assert "❌" in captured.err
    assert "Invalid [red]literal[/red] file" in captured.err
    assert not captured.out


def test_pipeline_error_adapter_uses_the_same_literal_theme():
    from postgwas.pipeline import cli
    from unittest.mock import patch

    stream = StringIO()
    console = Console(file=stream, force_terminal=True, color_system="standard", no_color=False)
    with patch.object(cli, "console", console):
        cli._print_error("Pipeline failed:", ValueError("Invalid [blue]file[/blue]"))
    terminal = stream.getvalue()
    assert "\x1b[1;31m" in terminal
    assert "Invalid [blue]file[/blue]" in Text.from_ansi(terminal).plain


@pytest.mark.parametrize("fail", (False, True))
def test_pipeline_registered_stage_rendering_preserves_order_and_failure(
    tmp_path, monkeypatch, fail,
):
    """Exercise orchestration with deterministic stub analyses/checkpoints."""
    from postgwas.core.errors import ModuleExecutionError
    from postgwas.pipeline import executor
    from postgwas.pipeline.planner import PipelinePlan

    stages = tuple(name for name in REGISTRY.names()
                   if REGISTRY.get(name).pipeline_enabled)
    plan = PipelinePlan(stages, stages, stages)
    config = load_configuration()
    args = Namespace(output_directory=str(tmp_path), resume=True, overwrite=False)
    stream = StringIO()
    console = Console(file=stream, force_terminal=False, width=120)
    monkeypatch.setattr(executor, "console", console)
    visited = []

    def resolve(reference):
        def run(*arguments):
            if len(arguments) == 1:  # module-specific title factory
                return "Selected analysis"
            visited.append(reference)
            print_screen_message("warning", "Shared stage warning", console=console)
            if fail and len(visited) == 2:
                raise ValueError("Invalid fixture [path]")
        return run

    checkpoint = SimpleNamespace(
        prepare=lambda: SimpleNamespace(action="run"),
        write=lambda **_: None, should_record_partial=True,
    )
    monkeypatch.setattr(executor, "resolve_reference", resolve)
    monkeypatch.setattr(executor, "_pipeline_checkpoint",
                        lambda **_: (checkpoint, tmp_path / "stub-checkpoint.yaml"))
    with terminal_presentation(config.logging):
        if fail:
            with pytest.raises(ModuleExecutionError, match="Invalid fixture"):
                executor.execute_pipeline(args, plan, config)
        else:
            executor.execute_pipeline(args, plan, config)
    text = stream.getvalue()
    assert visited == [REGISTRY.get(name).runner for name in stages[:len(visited)]]
    if fail:
        assert len(visited) == 2
        assert "Failed 2/%d" % len(stages) in text
        assert "All tasks completed successfully" not in text
    else:
        assert len(visited) == len(stages)
        assert "Completed %d/%d" % (len(stages), len(stages)) in text
        assert "All tasks completed successfully" in text
    assert text.count("Shared stage warning") == len(visited)


def test_literal_controls_cannot_bypass_the_shared_renderer():
    value = "❌  Input \x1b[35mBAD\x1b[0m \x1b]8;;url\x07file\x1b]8;;\x07"
    rendered = style_screen_block(value)
    assert rendered.plain == "❌  Input BAD file"
    assert all(str(span.style) == "bold red" for span in rendered.spans)


def test_configured_layout_and_column_styles_survive_background_refresh():
    from rich.progress import Task

    logging = load_configuration(cli_overrides={
        "logging.terminal_style.heading_indent": 0,
        "logging.terminal_style.field_indent": 4,
        "logging.terminal_style.indent_step": 2,
        "logging.terminal_style.styles.text": "blue",
    }).logging
    with terminal_presentation(logging):
        assert screen_line("analysis", "Heading", indent=2).startswith("🔬")
        assert screen_line("info", "Field", indent=6).startswith("    🔹")
        assert screen_line("info", "Nested", indent=10).startswith("      🔹")
        columns = StageProgress._columns()
    task = Task(0, "Fixture", total=10, completed=3, _get_time=lambda: 0.0)
    for column in columns[2:]:
        assert str(column.render(task).style) == "blue"


def test_multiline_field_labels_keep_the_same_label_column():
    field = screen_field("info", "one two three four five six seven", "1", label_width=6)
    lines = field.splitlines()
    label_column = cell_len(lines[0].split("one", 1)[0])
    assert len(lines) > 3
    for line in lines[1:]:
        assert cell_len(line) - cell_len(line.lstrip()) == label_column


@pytest.mark.parametrize("width", (60, 80, 120, 180))
def test_long_live_descriptions_never_hide_counts_percentage_or_elapsed(width):
    from rich.progress import Progress

    output = StringIO()
    console = Console(file=output, width=width, color_system=None)
    progress = Progress(*StageProgress._columns(), console=console, expand=True)
    progress.add_task("Current 6,652/18,569 · Test mapped genes and validate "
                      "every required output before reporting completion", total=18569, completed=6652)
    console.print(progress)
    text = output.getvalue()
    assert "6652/18569" in text
    assert "36%" in text
    assert "0:00:00" in text
    if width >= 80:
        assert "Test mapped genes" in text
    assert all(cell_len(line) <= width for line in text.splitlines())
