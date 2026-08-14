"""Integration tests for the global legacy direct-command checkpoint boundary."""

from io import StringIO
from pathlib import Path

import pytest
from rich.console import Console

from postgwas.config import load_configuration
from postgwas.core.direct_execution import run_direct_with_checkpoint


def _configuration(output: Path, overrides=None):
    return load_configuration(cli_overrides={
        "run.output_directory": str(output),
        **dict(overrides or {}),
    })


def _run(
    operation,
    output: Path,
    *,
    arguments=(),
    stream=None,
    configuration_overrides=None,
):
    return run_direct_with_checkpoint(
        operation,
        module_name="manhattan",
        public_command="manhattan",
        arguments=("--output-directory", str(output), *arguments),
        configuration=_configuration(output, configuration_overrides),
        output_directory=output,
        console=Console(file=stream or StringIO(), force_terminal=False),
    )


def test_completed_direct_command_resumes_without_calling_operation(tmp_path):
    output = tmp_path / "results"
    calls = []

    def operation():
        calls.append(1)
        output.mkdir(parents=True, exist_ok=True)
        (output / "plot.png").write_bytes(b"validated plot")
        return 0

    assert _run(operation, output) == 0
    assert _run(operation, output) == 0

    assert calls == [1]
    assert (
        output / "run_metadata/checkpoints/manhattan_direct.yaml"
    ).is_file()


def test_changed_direct_parameter_warns_and_rewrites_owned_outputs(tmp_path):
    output = tmp_path / "results"
    calls = []
    stream = StringIO()

    def operation():
        calls.append(len(calls) + 1)
        output.mkdir(parents=True, exist_ok=True)
        (output / "plot.png").write_text(
            "attempt=%d\n" % calls[-1], encoding="utf-8",
        )
        return 0

    _run(operation, output, arguments=("--minimum-p", "0.05"))
    _run(
        operation,
        output,
        arguments=("--minimum-p", "0.01"),
        stream=stream,
    )

    assert calls == [1, 2]
    assert (output / "plot.png").read_text(encoding="utf-8") == "attempt=2\n"
    assert "Resolved parameters changed" in stream.getvalue()
    assert "changed_parameters" in (
        output / "run_metadata/checkpoints/checkpoint_events.log"
    ).read_text(encoding="utf-8")


def test_failed_direct_command_records_partial_and_restarts_safely(tmp_path):
    output = tmp_path / "results"
    calls = []
    stream = StringIO()

    def operation():
        calls.append(len(calls) + 1)
        output.mkdir(parents=True, exist_ok=True)
        (output / "plot.png").write_text(
            "attempt=%d\n" % calls[-1], encoding="utf-8",
        )
        if len(calls) == 1:
            raise RuntimeError("simulated interruption")
        return 0

    with pytest.raises(RuntimeError, match="simulated interruption"):
        _run(operation, output)
    _run(operation, output, stream=stream)

    assert calls == [1, 2]
    assert (output / "plot.png").read_text(encoding="utf-8") == "attempt=2\n"
    assert "real partial result exists" in stream.getvalue()
    assert "partial_results" in (
        output / "run_metadata/checkpoints/checkpoint_events.log"
    ).read_text(encoding="utf-8")


def test_modified_direct_output_is_preserved_and_operation_is_not_called(tmp_path):
    output = tmp_path / "results"
    calls = []

    def operation():
        calls.append(1)
        output.mkdir(parents=True, exist_ok=True)
        (output / "plot.png").write_bytes(b"validated plot")
        return 0

    _run(operation, output)
    changed = output / "plot.png"
    changed.write_bytes(b"user replacement")

    with pytest.raises(RuntimeError, match="modified after validation"):
        _run(operation, output)

    assert calls == [1]
    assert changed.read_bytes() == b"user replacement"


def test_no_resume_failure_does_not_replace_prior_completed_manifest(tmp_path):
    output = tmp_path / "results"

    def completed():
        output.mkdir(parents=True, exist_ok=True)
        (output / "plot.png").write_bytes(b"validated plot")
        return 0

    _run(completed, output)
    manifest = output / "run_metadata/checkpoints/manhattan_direct.yaml"
    before = manifest.read_bytes()

    def failed():
        raise RuntimeError("rerun failed before output mutation")

    with pytest.raises(RuntimeError, match="rerun failed"):
        _run(failed, output, arguments=("--no-resume",))

    assert manifest.read_bytes() == before
    assert (output / "plot.png").read_bytes() == b"validated plot"


def test_run_controls_do_not_invalidate_completed_scientific_content(tmp_path):
    output = tmp_path / "results"
    calls = []

    def operation():
        calls.append(len(calls) + 1)
        output.mkdir(parents=True, exist_ok=True)
        (output / "plot.png").write_text(
            "attempt=%d\n" % calls[-1], encoding="utf-8",
        )
        return 0

    _run(
        operation,
        output,
        arguments=("--overwrite",),
        configuration_overrides={"run.overwrite": True},
    )
    _run(operation, output, arguments=("--resume",))

    assert calls == [1]
    assert (output / "plot.png").read_text(encoding="utf-8") == "attempt=1\n"
