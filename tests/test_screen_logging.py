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
    ((None, True), ("--hide-screen", False)),
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


def test_hidden_yaml_setting_remains_authoritative_without_cli_override(
    tmp_path: Path,
):
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

    assert hidden.returncode == 0
    assert hidden.stdout == ""
    assert hidden.stderr == ""
    transcript = output / "transcripts" / "terminal.log"
    text = transcript.read_text(encoding="utf-8")
    assert text.count("Choose the final analysis") == 1


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


def test_parent_stage_does_not_resume_while_measured_child_is_active(
    tmp_path: Path,
):
    transcript = tmp_path / "measured-child-screen.log"
    completed = _run_python(
        "\n".join((
            "from pathlib import Path",
            "import sys",
            "from unittest.mock import patch",
            "from rich.console import Console",
            "from postgwas.core.screen_logging import ScreenSettings, record_screen",
            "from postgwas.core.ui import MeasuredProgress, PipelineStageController, StageProgress",
            "transcript = Path(sys.argv[1])",
            "with patch('postgwas.core.screen_logging.os.isatty', return_value=True):",
            "    with record_screen(ScreenSettings(True, transcript)):",
            "        console = Console()",
            "        outer = StageProgress('Module progress', enabled=True, console=console)",
            "        stages = PipelineStageController(",
            "            'Scientific stages', ('Run tool', 'Validate results'),",
            "            console=console,",
            "        )",
            "        outer.start_step(1, 1, 'Run module')",
            "        stages.start(1)",
            "        measured = MeasuredProgress(",
            "            'Native measured progress', enabled=True, console=console,",
            "        )",
            "        measured.start('Test genes', total=10)",
            "        measured.update(9, total=10)",
            "        stages.complete(1)",
            "        stages.start(2)",
            "        measured.complete(10, total=10, title='Validate result rows')",
            "        stages.complete(2)",
            "        outer.complete_step(1, 1, 'Run module')",
        )),
        str(transcript),
    )

    assert completed.returncode == 0, completed.stderr
    assert "RecursionError" not in completed.stdout + completed.stderr
    assert "Completed 2/2 · Validate results" in completed.stdout
    assert "Completed 1/1 · Run module" in completed.stdout
    log = transcript.read_text(encoding="utf-8")
    assert log.count("Completed 1/2 · Run tool") == 1
    assert log.count("Completed 10/10 · Validate result rows") == 1
    assert log.count("Completed 2/2 · Validate results") == 1
    assert "━" not in log
    assert "\x1b" not in log


def test_multiline_result_block_is_preserved_while_parent_progress_is_live(
    tmp_path: Path,
):
    transcript = tmp_path / "summary-screen.log"
    completed = _run_python(
        "\n".join((
            "from pathlib import Path",
            "import sys",
            "from unittest.mock import patch",
            "from rich.console import Console",
            "from postgwas.core.screen_logging import ScreenSettings, record_screen",
            "from postgwas.core.ui import StageProgress, style_screen_block",
            "transcript = Path(sys.argv[1])",
            "with patch('postgwas.core.screen_logging.os.isatty', return_value=True):",
            "    with record_screen(ScreenSettings(True, transcript)):",
            "        console = Console()",
            "        outer = StageProgress('Module progress', enabled=True, console=console)",
            "        inner = StageProgress('Scientific stages', enabled=True, console=console)",
            "        outer.start_step(1, 1, 'Run module')",
            "        inner.start_step(1, 1, 'Save reports')",
            "        inner.complete_step(1, 1, 'Save reports')",
            "        inner.print_block(style_screen_block('  🔬  Summary\\n\\n    🔻  Variants in the MHC region\\n        🔹  failed this rule: 54,381; primary removals assigned here: 43,066; overlapping earlier removal rules: 11,315\\n'))",
            "        outer.complete_step(1, 1, 'Run module')",
        )),
        str(transcript),
    )

    assert completed.returncode == 0, completed.stderr
    for expected in (
        "Summary",
        "Variants in the MHC region",
        "failed this rule: 54,381; primary removals assigned here: 43,066",
        "Completed 1/1 · Run module",
    ):
        assert expected in completed.stdout
    log = transcript.read_text(encoding="utf-8")
    assert log.count("Variants in the MHC region") == 1
    assert log.count(
        "failed this rule: 54,381; primary removals assigned here: 43,066"
    ) == 1
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


def test_hidden_harmonisation_progress_reaches_both_required_transcripts(
    tmp_path: Path,
):
    transcript = tmp_path / "run_metadata" / "screen.log"
    dataset_report = tmp_path / "STUDY_screen_report.txt"
    dataset_log = tmp_path / "STUDY_dataset.log"

    completed = _run_python(
        "\n".join((
            "from pathlib import Path",
            "import sys",
            "from postgwas.core.pipeline_logging import PipelineLogger",
            "from postgwas.core.screen_logging import ScreenSettings, record_screen",
            "from postgwas.modules.harmonisation.cli import _DatasetScreenRouter",
            "from postgwas.modules.harmonisation.service import _HarmonisationProgress",
            "transcript, report, dataset_log = map(Path, sys.argv[1:])",
            "with record_screen(ScreenSettings(False, transcript)):",
            "    with _DatasetScreenRouter({'STUDY': report}, display=False) as router:",
            "        router.select('STUDY')",
            "        progress = _HarmonisationProgress(enabled=True, outcome_label_width=24)",
            "        logger = PipelineLogger(",
            "            'STUDY', 'dataset', str(dataset_log.parent),",
            "            log_path=str(dataset_log),",
            "            stage_progress=progress.dataset_stages,",
            "        )",
            "        with logger.step(1, 1, 'Validate dataset', 'test.validate'):",
            "            pass",
            "        progress.start_chromosomes(('1', '2'))",
            "        progress.record_chromosome_result('1', 'ok', 1)",
            "        progress.record_chromosome_result('2', 'ok', 1)",
            "        progress.finish_chromosomes()",
            "        logger.close()",
            "        progress.close()",
        )),
        str(transcript),
        str(dataset_report),
        str(dataset_log),
    )

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout == ""
    assert completed.stderr == ""
    for path in (transcript, dataset_report):
        text = path.read_text(encoding="utf-8")
        assert "Harmonisation dataset stages" in text
        assert "Completed 1/1 · Validate dataset" in text
        assert "Harmonisation chromosome progress" in text
        assert (
            "Completed 2/2 · All chromosome outputs and required counts validated"
            in text
        )


def test_screen_recorder_does_not_wait_for_inherited_pipe_writers(tmp_path: Path):
    transcript = tmp_path / "multiprocessing-screen.log"
    completed = _run_python(
        "\n".join((
            "from pathlib import Path",
            "import subprocess",
            "import sys",
            "import time",
            "from postgwas.core.screen_logging import ScreenSettings, record_screen",
            "transcript = Path(sys.argv[1])",
            "child = None",
            "started = time.monotonic()",
            "try:",
            "    with record_screen(ScreenSettings(False, transcript)):",
            "        print('recorded before inherited writer shutdown', flush=True)",
            "        child = subprocess.Popen([",
            "            sys.executable, '-c', 'import time; time.sleep(3)'",
            "        ])",
            "    print(time.monotonic() - started)",
            "finally:",
            "    if child is not None and child.poll() is None:",
            "        child.terminate()",
            "        child.wait(timeout=5)",
        )),
        str(transcript),
    )

    assert completed.returncode == 0, completed.stderr
    assert float(completed.stdout.strip()) < 1.0
    assert completed.stderr == ""
    assert transcript.read_text(encoding="utf-8") == (
        "recorded before inherited writer shutdown\n"
    )
