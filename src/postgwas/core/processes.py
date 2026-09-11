"""Checked external-process execution shared by PostGWAS modules."""

from __future__ import annotations

import os
import queue
import signal
import shutil
import subprocess
import tempfile
import threading
import time
from pathlib import Path
from typing import Callable, Mapping, Sequence, Type


class SupervisedProcessTimeout(TimeoutError):
    """A supervised external process exceeded its configured wall-clock limit."""

    def __init__(
        self,
        command,
        timeout_seconds,
        elapsed_seconds,
        stdout,
        stderr,
        forced_termination,
    ):
        self.command = tuple(str(value) for value in command)
        self.timeout_seconds = float(timeout_seconds)
        self.elapsed_seconds = float(elapsed_seconds)
        self.stdout = stdout or ""
        self.stderr = stderr or ""
        self.forced_termination = bool(forced_termination)
        super().__init__(
            "Command exceeded %.12g seconds after %.3f seconds"
            % (self.timeout_seconds, self.elapsed_seconds)
        )


def _signal_process_tree(process, sig):
    """Signal the isolated process tree, using its process group on POSIX."""
    if os.name == "posix":
        try:
            os.killpg(process.pid, sig)
        except ProcessLookupError:
            pass
        except PermissionError:
            # A watcher may reap the group leader between poll() and killpg().
            # Fall back to the direct child when it is still alive so cleanup
            # never hides the original external-tool failure.
            if process.poll() is None:
                try:
                    process.send_signal(sig)
                except (ProcessLookupError, PermissionError):
                    pass
        return
    if sig == signal.SIGTERM:
        process.terminate()
    else:
        process.kill()


def _process_tree_exists(process):
    """Return whether the isolated POSIX process group still has members."""
    if os.name != "posix":
        return process.poll() is None
    try:
        os.killpg(process.pid, 0)
    except ProcessLookupError:
        return False
    return True


def run_supervised_process(
    arguments: Sequence[str],
    *,
    cwd: str | Path | None = None,
    timeout_seconds: float | None = None,
    termination_grace_seconds: float,
) -> subprocess.CompletedProcess:
    """Run a process in an isolated group and terminate its tree on timeout."""
    command = [str(value) for value in arguments]
    if timeout_seconds is not None and float(timeout_seconds) <= 0:
        raise ValueError("timeout_seconds must be greater than zero")
    if float(termination_grace_seconds) <= 0:
        raise ValueError("termination_grace_seconds must be greater than zero")

    popen_options = {
        "cwd": None if cwd is None else str(cwd),
        "stdout": subprocess.PIPE,
        "stderr": subprocess.PIPE,
        "text": True,
    }
    if os.name == "posix":
        popen_options["start_new_session"] = True
    elif os.name == "nt":
        popen_options["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP

    started = time.monotonic()
    process = subprocess.Popen(command, **popen_options)
    try:
        stdout, stderr = process.communicate(timeout=timeout_seconds)
    except subprocess.TimeoutExpired:
        _signal_process_tree(process, signal.SIGTERM)
        forced = False
        try:
            stdout, stderr = process.communicate(
                timeout=float(termination_grace_seconds)
            )
            if _process_tree_exists(process):
                forced = True
                _signal_process_tree(process, signal.SIGKILL)
        except subprocess.TimeoutExpired:
            forced = True
            _signal_process_tree(process, signal.SIGKILL)
            stdout, stderr = process.communicate()
        raise SupervisedProcessTimeout(
            command=command,
            timeout_seconds=timeout_seconds,
            elapsed_seconds=time.monotonic() - started,
            stdout=stdout,
            stderr=stderr,
            forced_termination=forced,
        ) from None
    return subprocess.CompletedProcess(
        command, process.returncode, stdout=stdout, stderr=stderr
    )


def run_checked_command(
    arguments: Sequence[str],
    purpose: str,
    *,
    logger=None,
    error_type: Type[Exception] = RuntimeError,
    stdout_path: str | Path | None = None,
    append_stdout: bool = False,
    stderr_to_stdout: bool = False,
    stdout_header: str | None = None,
    timeout_seconds: float | None = None,
    expected_outputs: Sequence[str | Path] = (),
    dry_run: bool = False,
    progress_callback: Callable[[], None] | None = None,
    progress_refresh_seconds: float | None = None,
    resource_metrics: dict[str, int] | None = None,
    resource_poll_seconds: float | None = None,
    env: Mapping[str, str] | None = None,
) -> str:
    """Run one command safely, optionally streaming stdout to a file.

    Large VCF queries must not be captured in memory.  ``stdout_path`` keeps
    the same checked execution and audit behaviour while streaming output.
    ``env`` follows the native subprocess contract: a supplied mapping replaces
    the child's environment; ``None`` preserves normal parent inheritance.
    Neither the parent environment nor the supplied mapping is modified, and
    environment values are not added to command audit records. Callers needing
    overrides plus inherited variables must supply their own merged mapping.
    See https://docs.python.org/3/library/subprocess.html#subprocess.run.
    """
    command = [str(value) for value in arguments]
    if progress_callback is not None:
        if progress_refresh_seconds is None or float(progress_refresh_seconds) <= 0:
            raise ValueError(
                "progress_refresh_seconds must be greater than zero when a "
                "progress callback is supplied"
            )
    if resource_metrics is not None:
        if resource_poll_seconds is None or float(resource_poll_seconds) <= 0:
            raise ValueError(
                "resource_poll_seconds must be greater than zero when resource "
                "metrics are requested"
            )
        resource_metrics.clear()
        resource_metrics.update({
            "peak_rss_bytes": 0,
            "memory_samples": 0,
            "process_tree_samples": 0,
            "root_only_samples": 0,
            "sampling_failures": 0,
        })
    if logger is not None and hasattr(logger, "record"):
        logger.record("INPUT", "external_command", purpose=purpose, command=command)
    if dry_run:
        if logger is not None:
            logger.record("SKIP", "external_command", purpose=purpose, reason="dry_run")
        return ""
    output_handle = None
    try:
        if stdout_path is not None:
            destination = Path(stdout_path)
            destination.parent.mkdir(parents=True, exist_ok=True)
            output_handle = destination.open(
                "a" if append_stdout else "w", encoding="utf-8", newline="",
            )
            if stdout_header:
                output_handle.write(stdout_header)
                output_handle.flush()
        process_stdout = (
            output_handle if output_handle is not None else subprocess.PIPE
        )
        process_stderr = (
            subprocess.STDOUT if stderr_to_stdout else subprocess.PIPE
        )
        environment_options = {} if env is None else {"env": env}
        if progress_callback is None and resource_metrics is None:
            result = subprocess.run(
                command,
                stdout=process_stdout,
                stderr=process_stderr,
                text=True,
                timeout=timeout_seconds,
                **environment_options,
            )
        else:
            started = time.monotonic()
            process = subprocess.Popen(
                command,
                stdout=process_stdout,
                stderr=process_stderr,
                text=True,
                **environment_options,
            )
            callback = progress_callback

            def sample_process_memory() -> None:
                if resource_metrics is None:
                    return
                try:
                    import psutil

                    root = psutil.Process(process.pid)
                    processes = [root]
                    process_tree_complete = True
                    try:
                        processes.extend(root.children(recursive=True))
                    except (
                        OSError,
                        psutil.AccessDenied,
                        psutil.NoSuchProcess,
                        psutil.ZombieProcess,
                    ):
                        # Some managed macOS environments allow the direct
                        # child RSS query but deny process-tree enumeration.
                        # Retain the useful root measurement and expose that
                        # the descendant component was unavailable.
                        process_tree_complete = False
                    rss_bytes = 0
                    for observed in processes:
                        try:
                            rss_bytes += int(observed.memory_info().rss)
                        except (
                            psutil.AccessDenied,
                            psutil.NoSuchProcess,
                            psutil.ZombieProcess,
                        ):
                            continue
                    resource_metrics["peak_rss_bytes"] = max(
                        int(resource_metrics["peak_rss_bytes"]), rss_bytes,
                    )
                    resource_metrics["memory_samples"] = (
                        int(resource_metrics["memory_samples"]) + 1
                    )
                    sample_kind = (
                        "process_tree_samples"
                        if process_tree_complete
                        else "root_only_samples"
                    )
                    resource_metrics[sample_kind] = (
                        int(resource_metrics[sample_kind]) + 1
                    )
                except Exception:
                    # Resource sampling is advisory. Checked command execution
                    # and its scientific outputs remain authoritative. A process
                    # that exits between communicate() and this sample is a
                    # normal race, not an execution failure.
                    resource_metrics["sampling_failures"] = (
                        int(resource_metrics["sampling_failures"]) + 1
                    )
                    return

            def refresh_progress() -> None:
                nonlocal callback
                if callback is None:
                    return
                try:
                    callback()
                except Exception as exc:
                    callback = None
                    if logger is not None and hasattr(logger, "record"):
                        logger.record(
                            "WARNING",
                            "external_progress_disabled",
                            purpose=purpose,
                            error="%s: %s" % (type(exc).__name__, exc),
                        )

            refresh_progress()
            sample_process_memory()
            try:
                while True:
                    elapsed = time.monotonic() - started
                    remaining = (
                        None
                        if timeout_seconds is None
                        else float(timeout_seconds) - elapsed
                    )
                    if remaining is not None and remaining <= 0:
                        process.kill()
                        stdout, stderr = process.communicate()
                        raise subprocess.TimeoutExpired(
                            command,
                            timeout_seconds,
                            output=stdout,
                            stderr=stderr,
                        )
                    wait_candidates = []
                    if progress_callback is not None:
                        wait_candidates.append(float(progress_refresh_seconds))
                    if resource_metrics is not None:
                        wait_candidates.append(float(resource_poll_seconds))
                    wait_seconds = min(wait_candidates)
                    if remaining is not None:
                        wait_seconds = min(wait_seconds, remaining)
                    try:
                        stdout, stderr = process.communicate(timeout=wait_seconds)
                    except subprocess.TimeoutExpired:
                        refresh_progress()
                        sample_process_memory()
                        continue
                    refresh_progress()
                    sample_process_memory()
                    result = subprocess.CompletedProcess(
                        command,
                        process.returncode,
                        stdout=stdout,
                        stderr=stderr,
                    )
                    break
            except BaseException:
                if process.poll() is None:
                    process.kill()
                    process.communicate()
                raise
        if output_handle is not None and stdout_header is not None:
            output_handle.write("\nexit_code=%s\n" % result.returncode)
            output_handle.flush()
    except subprocess.TimeoutExpired as exc:
        if output_handle is not None and stdout_header is not None:
            output_handle.write("\nexecution_error=timeout\n")
            output_handle.flush()
        raise error_type(
            "%s exceeded the configured timeout of %s seconds"
            % (purpose, timeout_seconds)
        ) from exc
    except OSError as exc:
        if output_handle is not None and stdout_header is not None:
            output_handle.write("\nexecution_error=%s\n" % exc)
            output_handle.flush()
        raise error_type("%s could not start: %s" % (purpose, exc)) from exc
    finally:
        if output_handle is not None:
            output_handle.close()
    if result.returncode != 0:
        diagnostic = "\n".join(
            value for value in (result.stdout, result.stderr) if value
        )
        if stderr_to_stdout and stdout_path is not None:
            try:
                diagnostic = Path(stdout_path).read_text(
                    encoding="utf-8", errors="replace",
                )
            except OSError:
                diagnostic = ""
        stderr = diagnostic.strip().splitlines()
        raise error_type(
            "%s failed with exit status %d: %s"
            % (purpose, result.returncode, "\n".join(stderr[-20:]))
        )
    if result.stdout and logger is not None:
        logger.info("%s stdout: %s" % (purpose, result.stdout.strip()))
    if result.stderr and result.stderr.strip() and logger is not None:
        logger.info("%s stderr: %s" % (purpose, result.stderr.strip()))
    missing = [
        str(Path(value))
        for value in expected_outputs
        if not Path(value).is_file() or Path(value).stat().st_size <= 0
    ]
    if missing:
        raise error_type(
            "%s returned success but expected output is missing or empty: %s"
            % (purpose, ", ".join(missing))
        )
    return "" if stdout_path is not None else (result.stdout or "")


def run_checked_line_command(
    arguments: Sequence[str],
    purpose: str,
    *,
    line_consumer: Callable[[str], None],
    logger=None,
    error_type: Type[Exception] = RuntimeError,
) -> int:
    """Consume a checked command's stdout incrementally without buffering it.

    This is intended for genome-scale textual streams where capturing stdout or
    materialising an intermediate table would be unsafe. Stderr is spooled to a
    temporary file so verbose diagnostics cannot deadlock the child process.
    """
    command = [str(value) for value in arguments]
    if not command:
        raise ValueError("arguments must contain a nonempty command")
    if not callable(line_consumer):
        raise TypeError("line_consumer must be callable")
    if logger is not None and hasattr(logger, "record"):
        logger.record("INPUT", "external_command", purpose=purpose, command=command)

    popen_options = {}
    if os.name == "posix":
        popen_options["start_new_session"] = True
    elif os.name == "nt":
        popen_options["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP

    stderr_handle = tempfile.TemporaryFile(mode="w+b")
    try:
        try:
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=stderr_handle,
                text=True,
                encoding="utf-8",
                **popen_options,
            )
        except OSError as exc:
            raise error_type("%s could not start: %s" % (purpose, exc)) from exc

        line_count = 0
        try:
            if process.stdout is None:
                raise error_type("%s did not expose a readable stdout stream" % purpose)
            try:
                for line in process.stdout:
                    line_consumer(line)
                    line_count += 1
            except UnicodeError as exc:
                if process.poll() is None:
                    _signal_process_tree(process, signal.SIGKILL)
                process.wait()
                raise error_type(
                    "%s produced output that is not valid UTF-8" % purpose
                ) from exc
            except BaseException:
                if process.poll() is None:
                    _signal_process_tree(process, signal.SIGKILL)
                process.wait()
                raise
            return_code = process.wait()
        finally:
            if process.stdout is not None:
                process.stdout.close()

        stderr_handle.seek(0)
        stderr = stderr_handle.read().decode("utf-8", errors="replace")
        if return_code != 0:
            detail = "\n".join(stderr.strip().splitlines()[-20:])
            raise error_type(
                "%s failed with exit status %d: %s"
                % (purpose, return_code, detail or "no stderr was produced")
            )
        if stderr.strip() and logger is not None:
            logger.info("%s stderr: %s" % (purpose, stderr.strip()))
        return line_count
    finally:
        stderr_handle.close()


def run_checked_pipeline(
    commands: Sequence[Sequence[str]],
    purpose: str,
    *,
    logger=None,
    error_type: Type[Exception] = RuntimeError,
    expected_outputs: Sequence[str | Path] = (),
) -> tuple[subprocess.CompletedProcess, ...]:
    """Run a binary-safe process pipeline and check every stage.

    Intermediate stdout is connected directly to the next stage. Each stderr
    stream and the final stdout stream are spooled to temporary files so an
    external tool cannot deadlock the pipeline by filling an unread pipe. If
    any stage fails, every still-running process tree is terminated before the
    collected stage diagnostics are reported.
    """
    pipeline = [[str(value) for value in command] for command in commands]
    if not pipeline or any(not command for command in pipeline):
        raise ValueError("commands must contain one or more nonempty commands")
    if logger is not None and hasattr(logger, "record"):
        logger.record(
            "INPUT",
            "external_pipeline",
            purpose=purpose,
            commands=pipeline,
        )

    popen_group_options = {}
    if os.name == "posix":
        popen_group_options["start_new_session"] = True
    elif os.name == "nt":
        popen_group_options["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP

    processes: list[subprocess.Popen] = []
    stderr_handles = []
    previous_stdout = None
    final_stdout_handle = tempfile.TemporaryFile(mode="w+b")

    def terminate_running_processes() -> None:
        for process in processes:
            if process.poll() is None:
                _signal_process_tree(process, signal.SIGKILL)

    try:
        for index, command in enumerate(pipeline):
            stderr_handle = tempfile.TemporaryFile(mode="w+b")
            stderr_handles.append(stderr_handle)
            process = subprocess.Popen(
                command,
                stdin=previous_stdout,
                stdout=(
                    final_stdout_handle
                    if index == len(pipeline) - 1
                    else subprocess.PIPE
                ),
                stderr=stderr_handle,
                **popen_group_options,
            )
            processes.append(process)
            if previous_stdout is not None:
                previous_stdout.close()
            previous_stdout = process.stdout
    except (OSError, ValueError) as exc:
        if previous_stdout is not None:
            previous_stdout.close()
        terminate_running_processes()
        for process in processes:
            process.wait()
        final_stdout_handle.close()
        for handle in stderr_handles:
            handle.close()
        raise error_type("%s could not start: %s" % (purpose, exc)) from exc

    completed_queue: queue.SimpleQueue[tuple[int, int]] = queue.SimpleQueue()
    watchers = []

    def wait_for_stage(index: int, process: subprocess.Popen) -> None:
        completed_queue.put((index, process.wait()))

    try:
        for index, process in enumerate(processes):
            watcher = threading.Thread(
                target=wait_for_stage,
                args=(index, process),
                daemon=True,
            )
            watcher.start()
            watchers.append(watcher)

        return_codes: dict[int, int] = {}
        while len(return_codes) < len(processes):
            index, return_code = completed_queue.get()
            return_codes[index] = return_code
            if return_code != 0:
                terminate_running_processes()
        for watcher in watchers:
            watcher.join()
    except BaseException:
        terminate_running_processes()
        for process in processes:
            process.wait()
        final_stdout_handle.close()
        for handle in stderr_handles:
            handle.close()
        raise
    finally:
        if previous_stdout is not None:
            previous_stdout.close()

    final_stdout_handle.seek(0)
    final_stdout = final_stdout_handle.read().decode("utf-8", errors="replace")
    stderr_values = []
    for handle in stderr_handles:
        handle.seek(0)
        stderr_values.append(handle.read().decode("utf-8", errors="replace"))
        handle.close()
    final_stdout_handle.close()

    results = tuple(
        subprocess.CompletedProcess(
            command,
            return_codes[index],
            stdout=final_stdout if index == len(pipeline) - 1 else "",
            stderr=stderr_values[index],
        )
        for index, command in enumerate(pipeline)
    )
    failed = [
        (index, result)
        for index, result in enumerate(results, start=1)
        if result.returncode != 0
    ]
    if logger is not None and hasattr(logger, "record"):
        for index, result in enumerate(results, start=1):
            logger.record(
                "ERROR" if result.returncode else "OUTPUT",
                "external_pipeline_stage",
                purpose=purpose,
                stage=index,
                exit_code=result.returncode,
                stderr=result.stderr.strip(),
            )
    if failed:
        diagnostics = []
        for index, result in failed:
            lines = result.stderr.strip().splitlines()
            detail = "\n".join(lines[-20:]) or "no stderr was produced"
            diagnostics.append(
                "stage %d failed with exit status %d (%s): %s"
                % (index, result.returncode, " ".join(result.args), detail)
            )
        raise error_type("%s failed: %s" % (purpose, "\n".join(diagnostics)))

    missing = [
        str(Path(value))
        for value in expected_outputs
        if not Path(value).is_file() or Path(value).stat().st_size <= 0
    ]
    if missing:
        raise error_type(
            "%s returned success but expected output is missing or empty: %s"
            % (purpose, ", ".join(missing))
        )
    return results


def build_container_command(
    runtime: str,
    image: str,
    arguments: Sequence[str],
    *,
    mount_directories: Sequence[str | Path],
    platform: str | None = None,
    work_directory: str | Path | None = None,
    use_host_user: bool = True,
    validate: bool = True,
) -> list[str]:
    """Build a shell-free container command with same-path read/write mounts."""
    executable = shutil.which(str(runtime)) if validate else str(runtime)
    if validate and not executable:
        raise RuntimeError("Container runtime was not found: %s" % runtime)

    mounts = []
    for value in mount_directories:
        path = Path(value).expanduser().resolve()
        if validate and not path.is_dir():
            raise RuntimeError("Container mount directory does not exist: %s" % path)
        if path == Path(path.anchor):
            raise RuntimeError("Refusing to mount a filesystem root into a tool container")
        if path not in mounts:
            mounts.append(path)

    command = [str(executable), "run", "--rm"]
    if platform:
        command.extend(["--platform", str(platform)])
    if use_host_user and hasattr(os, "getuid") and hasattr(os, "getgid"):
        command.extend(["--user", "%d:%d" % (os.getuid(), os.getgid())])
    for path in mounts:
        command.extend(["--volume", "%s:%s" % (path, path)])
    if work_directory is not None:
        command.extend(["--workdir", str(Path(work_directory).expanduser().resolve())])
    command.append(str(image))
    command.extend(str(value) for value in arguments)
    return command


__all__ = [
    "SupervisedProcessTimeout",
    "build_container_command",
    "run_checked_command",
    "run_checked_line_command",
    "run_checked_pipeline",
    "run_supervised_process",
]
