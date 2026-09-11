"""Explicit child environments must not mutate or disclose parent state."""

import os
import sys
from types import MappingProxyType
from unittest.mock import Mock

import pytest

from postgwas.core.processes import run_checked_command


@pytest.fixture(params=["blocking", "progress", "resources"])
def execution_options(request):
    if request.param == "progress":
        return {"progress_callback": lambda: None, "progress_refresh_seconds": 0.01}
    if request.param == "resources":
        return {"resource_metrics": {}, "resource_poll_seconds": 0.01}
    return {}


def test_native_child_gets_override_without_mutating_parent_or_mapping(monkeypatch, execution_options):
    monkeypatch.setenv("POSTGWAS_TEST_CHILD_VALUE", "parent-value")
    parent = dict(os.environ)
    child = MappingProxyType({**parent, "POSTGWAS_TEST_CHILD_VALUE": "child-value"})
    child_before = dict(child)
    command = [sys.executable, "-I", "-c", "import os; print(os.environ['POSTGWAS_TEST_CHILD_VALUE']); os.environ['POSTGWAS_TEST_CHILD_VALUE']='child-mutation'"]
    assert run_checked_command(command, "Child environment probe", env=child, **execution_options).strip() == "child-value"
    assert dict(child) == child_before
    assert dict(os.environ) == parent
    assert run_checked_command(command, "Default environment probe", **execution_options).strip() == "parent-value"
    assert run_checked_command(command, "Explicit None environment probe", env=None, **execution_options).strip() == "parent-value"
    assert dict(os.environ) == parent


def test_empty_child_environment_is_not_implicitly_merged(monkeypatch, execution_options):
    monkeypatch.setenv("POSTGWAS_TEST_PARENT_ONLY", "must-not-be-inherited")
    output = run_checked_command(
        [sys.executable, "-I", "-c", "import os; print('POSTGWAS_TEST_PARENT_ONLY' in os.environ)"],
        "Empty child environment probe", env={}, **execution_options,
    )
    assert output.strip() == "False"
    assert os.environ["POSTGWAS_TEST_PARENT_ONLY"] == "must-not-be-inherited"


def test_environment_is_not_logged_and_success_output_logging_is_preserved(execution_options):
    logger = Mock()
    private_marker = "environment-private-marker-not-for-audit"
    command = [sys.executable, "-I", "-c", "import sys; print('normal stdout'); print('normal stderr', file=sys.stderr)"]
    output = run_checked_command(
        command, "Logged environment probe", logger=logger,
        env={**os.environ, "POSTGWAS_TEST_PRIVATE": private_marker}, **execution_options,
    )
    assert output.strip() == "normal stdout"
    logger.record.assert_called_once_with("INPUT", "external_command", purpose="Logged environment probe", command=command)
    logger.info.assert_any_call("Logged environment probe stdout: normal stdout")
    logger.info.assert_any_call("Logged environment probe stderr: normal stderr")
    assert private_marker not in str(logger.mock_calls)


def test_nonzero_child_preserves_checked_error_and_process_transcript(tmp_path, execution_options):
    logger = Mock()
    destination = tmp_path / "child.log"
    private_marker = "environment-private-marker-not-for-audit"
    command = [sys.executable, "-I", "-c", "import sys; print('stdout diagnostic'); print('stderr diagnostic', file=sys.stderr); raise SystemExit(7)"]
    with pytest.raises(ValueError, match="Environment failure probe failed with exit status 7") as caught:
        run_checked_command(
            command, "Environment failure probe", logger=logger, error_type=ValueError,
            env={**os.environ, "POSTGWAS_TEST_PRIVATE": private_marker},
            stdout_path=destination, stderr_to_stdout=True, stdout_header="tool=fixture\n",
            **execution_options,
        )
    transcript = destination.read_text()
    assert transcript.startswith("tool=fixture\n") and transcript.endswith("\nexit_code=7\n")
    for message in ("stdout diagnostic", "stderr diagnostic"):
        assert message in transcript and message in str(caught.value)
    logger.record.assert_called_once_with("INPUT", "external_command", purpose="Environment failure probe", command=command)
    assert private_marker not in transcript + str(logger.mock_calls) + str(caught.value)


def test_missing_required_output_still_fails_with_explicit_environment(tmp_path, execution_options):
    with pytest.raises(RuntimeError, match="expected output is missing or empty"):
        run_checked_command(
            [sys.executable, "-I", "-c", "pass"], "Missing output probe", env={},
            expected_outputs=[tmp_path / "missing.txt"], **execution_options,
        )
