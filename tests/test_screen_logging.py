"""Integration tests for the shared terminal display and transcript controls."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def _test_environment() -> dict[str, str]:
    environment = os.environ.copy()
    source_directory = str(REPOSITORY_ROOT / "src")
    environment["PYTHONPATH"] = os.pathsep.join(
        value
        for value in (source_directory, environment.get("PYTHONPATH"))
        if value
    )
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    return environment


def _run_postgwas(*arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "postgwas", *arguments],
        cwd=REPOSITORY_ROOT,
        env=_test_environment(),
        check=False,
        capture_output=True,
        text=True,
    )


def _run_python(source: str, *arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-c", source, *arguments],
        cwd=REPOSITORY_ROOT,
        env=_test_environment(),
        check=False,
        capture_output=True,
        text=True,
    )


@pytest.mark.parametrize(
    ("visibility_argument", "visible"),
    ((None, True), ("--show-screen", True), ("--hide-screen", False)),
)
def test_pipeline_screen_stream_is_always_saved(
    tmp_path: Path,
    visibility_argument: str | None,
    visible: bool,
):
    output = tmp_path / (visibility_argument or "default")
    arguments = ["pipeline", "--output-directory", str(output)]
    if visibility_argument is not None:
        arguments.append(visibility_argument)

    completed = _run_postgwas(*arguments)

    assert completed.returncode == 0
    terminal = completed.stdout + completed.stderr
    assert ("Choose the final analysis" in terminal) is visible
    transcript = output / "run_metadata" / "screen.log"
    assert transcript.is_file()
    assert "Choose the final analysis" in transcript.read_text(encoding="utf-8")


def test_show_screen_overrides_a_hidden_yaml_default(tmp_path: Path):
    output = tmp_path / "configured-output"
    run_config = tmp_path / "run.yaml"
    run_config.write_text(
        "run:\n"
        "  output_directory: %s\n"
        "logging:\n"
        "  show_screen: false\n"
        "  screen_log_file: transcripts/terminal.log\n" % output,
        encoding="utf-8",
    )

    hidden = _run_postgwas("pipeline", "--run-config", str(run_config))
    shown = _run_postgwas(
        "pipeline",
        "--run-config",
        str(run_config),
        "--show-screen",
    )

    assert hidden.returncode == 0
    assert hidden.stdout == ""
    assert hidden.stderr == ""
    assert shown.returncode == 0
    assert "Choose the final analysis" in shown.stdout
    transcript = output / "transcripts" / "terminal.log"
    text = transcript.read_text(encoding="utf-8")
    assert text.count("Choose the final analysis") == 2


def test_mandatory_progress_ignores_the_optional_detail_setting(tmp_path: Path):
    output = tmp_path / "mandatory-progress"
    run_config = tmp_path / "run.yaml"
    run_config.write_text(
        "run:\n"
        "  output_directory: %s\n"
        "logging:\n"
        "  show_screen: false\n"
        "  show_progress: false\n" % output,
        encoding="utf-8",
    )

    completed = _run_postgwas(
        "magma",
        "--run-config",
        str(run_config),
        "--output-directory",
        str(output),
        "--not-a-postgwas-option",
    )

    assert completed.returncode == 2
    transcript = output / "run_metadata" / "screen.log"
    text = transcript.read_text(encoding="utf-8")
    assert "Magma progress" in text
    assert "Failed 1/1" in text
    assert "━" not in text
    assert "\x1b" not in text


def test_screen_recorder_uses_a_terminal_preserving_progress_channel(
    tmp_path: Path,
):
    transcript = tmp_path / "screen.log"
    completed = _run_python(
        "\n".join((
            "from pathlib import Path",
            "import sys",
            "from unittest.mock import patch",
            "from rich.console import Console",
            "from postgwas.core.screen_logging import ScreenSettings, record_screen",
            "from postgwas.core.ui import StageProgress",
            "transcript = Path(sys.argv[1])",
            "with patch('postgwas.core.screen_logging.os.isatty', return_value=True):",
            "    with record_screen(ScreenSettings(True, transcript)):",
            "        progress = StageProgress('Live progress', enabled=True, console=Console())",
            "        progress.start_step(1, 1, 'Run analysis')",
            "        progress.complete_step(1, 1, 'Run analysis')",
        )),
        str(transcript),
    )

    assert completed.returncode == 0, completed.stderr
    terminal = completed.stdout
    assert "Current 1/1 · Run analysis" in terminal
    assert "Started 1/1 · Run analysis" not in terminal
    assert "All 1 stages completed" in terminal
    log = transcript.read_text(encoding="utf-8")
    assert "Live progress" in log
    assert "Started 1/1 · Run analysis" in log
    assert "Completed 1/1 · Run analysis" in log
    assert "━" not in log
    assert "\x1b" not in log


def test_nested_live_progress_has_one_plain_record_per_stage_event(tmp_path: Path):
    transcript = tmp_path / "nested-screen.log"
    completed = _run_python(
        "\n".join((
            "from pathlib import Path",
            "import sys",
            "from unittest.mock import patch",
            "from rich.console import Console",
            "from postgwas.core.screen_logging import ScreenSettings, record_screen",
            "from postgwas.core.ui import StageProgress",
            "transcript = Path(sys.argv[1])",
            "with patch('postgwas.core.screen_logging.os.isatty', return_value=True):",
            "    with record_screen(ScreenSettings(True, transcript)):",
            "        console = Console()",
            "        outer = StageProgress('Module progress', enabled=True, console=console)",
            "        inner = StageProgress('Scientific stages', enabled=True, console=console)",
            "        outer.start_step(1, 1, 'Run module')",
            "        inner.start_step(1, 2, 'Prepare inputs')",
            "        print('warning while nested progress is active')",
            "        inner.complete_step(1, 2, 'Prepare inputs')",
            "        inner.start_step(2, 2, 'Validate results')",
            "        inner.fail_step(2, 2, 'Validate results')",
            "        outer.fail_step(1, 1, 'Run module')",
        )),
        str(transcript),
    )

    assert completed.returncode == 0, completed.stderr
    assert "None" not in completed.stdout
    assert "warning while nested progress is active" in completed.stdout
    log = transcript.read_text(encoding="utf-8")
    assert log.count("warning while nested progress is active") == 1
    assert log.count("Started 1/1 · Run module") == 1
    assert log.count("Started 1/2 · Prepare inputs") == 1
    assert log.count("Completed 1/2 · Prepare inputs") == 1
    assert log.count("Started 2/2 · Validate results") == 1
    assert log.count("Failed 2/2 · Validate results") == 1
    assert log.count("Failed 1/1 · Run module") == 1
    assert "━" not in log
    assert "\x1b" not in log


def test_hidden_measured_progress_preserves_unit_counts(tmp_path: Path):
    transcript = tmp_path / "measured-screen.log"
    completed = _run_python(
        "\n".join((
            "from pathlib import Path",
            "import sys",
            "from rich.console import Console",
            "from postgwas.core.screen_logging import ScreenSettings, record_screen",
            "from postgwas.core.ui import MeasuredProgress",
            "transcript = Path(sys.argv[1])",
            "with record_screen(ScreenSettings(False, transcript)):",
            "    progress = MeasuredProgress('Measured genes', enabled=True, console=Console())",
            "    progress.start('Load and match variants')",
            "    progress.set_phase('Test mapped genes')",
            "    progress.update(1, total=3)",
            "    progress.update(2, total=3)",
            "    progress.update(3, total=3)",
            "    progress.complete(3, total=3, title='Validate GCTA result rows')",
        )),
        str(transcript),
    )

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout == ""
    assert completed.stderr == ""
    log = transcript.read_text(encoding="utf-8")
    assert "Started 0/? · Load and match variants" in log
    assert "Current 0/? · Test mapped genes" in log
    assert "Progress 1/3 · Test mapped genes · 33%" in log
    assert "Progress 2/3 · Test mapped genes · 66%" in log
    assert "Progress 3/3" not in log
    assert "Completed 3/3 · Validate GCTA result rows" in log
    assert "━" not in log
    assert "\x1b" not in log


@pytest.mark.parametrize(
    ("command", "progress_label"),
    (("manhattan", "Manhattan progress"), ("magma", "Magma progress")),
)
def test_hidden_module_parse_errors_are_saved_to_the_screen_log(
    tmp_path: Path,
    command: str,
    progress_label: str,
):
    output = tmp_path / command

    completed = _run_postgwas(
        command,
        "--output-directory",
        str(output),
        "--hide-screen",
        "--not-a-postgwas-option",
    )

    assert completed.returncode == 2
    assert completed.stdout == ""
    assert completed.stderr == ""
    transcript = output / "run_metadata" / "screen.log"
    transcript_text = transcript.read_text(encoding="utf-8")
    assert "unrecognized arguments" in transcript_text
    assert progress_label in transcript_text
    assert "Failed 1/1" in transcript_text


def test_hidden_harmonisation_router_still_reaches_shared_transcript(
    tmp_path: Path,
):
    transcript = tmp_path / "run_metadata" / "screen.log"
    dataset_report = tmp_path / "STUDY_screen_report.txt"

    completed = _run_python(
        "\n".join((
            "from pathlib import Path",
            "import sys",
            "from postgwas.core.screen_logging import ScreenSettings, record_screen",
            "from postgwas.modules.harmonisation.cli import _DatasetScreenRouter",
            "transcript, report = map(Path, sys.argv[1:])",
            "with record_screen(ScreenSettings(False, transcript)):",
            "    with _DatasetScreenRouter({'STUDY': report}, display=False):",
            "        print('hidden harmonisation progress')",
        )),
        str(transcript),
        str(dataset_report),
    )

    assert completed.returncode == 0
    assert completed.stdout == ""
    assert completed.stderr == ""
    assert "hidden harmonisation progress" in transcript.read_text(encoding="utf-8")
    assert "hidden harmonisation progress" in dataset_report.read_text(
        encoding="utf-8"
    )
