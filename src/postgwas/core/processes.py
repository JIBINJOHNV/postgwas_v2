"""Checked external-process execution shared by PostGWAS modules."""

from __future__ import annotations

import os
import signal
import shutil
import subprocess
import time
from pathlib import Path
from typing import Callable, Sequence, Type


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
) -> str:
    """Run one command safely, optionally streaming stdout to a file.

    Large VCF queries must not be captured in memory.  ``stdout_path`` keeps
    the same checked execution and audit behaviour while streaming output.
    """
    command = [str(value) for value in arguments]
    if progress_callback is not None:
        if progress_refresh_seconds is None or float(progress_refresh_seconds) <= 0:
            raise ValueError(
                "progress_refresh_seconds must be greater than zero when a "
                "progress callback is supplied"
            )
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
        if progress_callback is None:
            result = subprocess.run(
                command,
                stdout=process_stdout,
                stderr=process_stderr,
                text=True,
                timeout=timeout_seconds,
            )
        else:
            started = time.monotonic()
            process = subprocess.Popen(
                command,
                stdout=process_stdout,
                stderr=process_stderr,
                text=True,
            )
            callback = progress_callback

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
                    wait_seconds = float(progress_refresh_seconds)
                    if remaining is not None:
                        wait_seconds = min(wait_seconds, remaining)
                    try:
                        stdout, stderr = process.communicate(timeout=wait_seconds)
                    except subprocess.TimeoutExpired:
                        refresh_progress()
                        continue
                    refresh_progress()
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
    "run_supervised_process",
]
