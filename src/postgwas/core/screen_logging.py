"""Record the complete PostGWAS screen stream while optionally displaying it."""

from __future__ import annotations

from contextlib import contextmanager, nullcontext
from contextvars import copy_context
from dataclasses import dataclass
import os
from pathlib import Path
import select
import sys
import threading
from typing import Any, Iterator, Sequence, TextIO

from rich.console import Console

from postgwas.core.ui.terminal_text import TerminalTextDecoder

from postgwas.config import load_configuration, load_run_configuration_for_module
from postgwas.config.loader import canonical_module_name
from postgwas.core.ui.screen import (
    default_terminal_settings,
    style_screen_block,
    terminal_presentation,
    terminal_settings,
)


_NON_SCIENTIFIC_COMMANDS = {"config", "resources"}
_ACTIVE_RECORDERS = 0
_PROGRESS_DISPLAY_CONSOLE: Console | None = None
_PROGRESS_DISPLAY_STREAM: TextIO | None = None
_PROGRESS_SUMMARY_CONSOLE: Console | None = None

# Bounded transport buffer, not an analysis or user-facing compute setting.
_PIPE_READ_BYTES = 65536


@dataclass(frozen=True)
class ScreenSettings:
    """Resolved display decision and destination for one public command."""

    show_screen: bool
    log_file: Path
    configuration: Any | None = None
    output_directory: Path | None = None


def screen_recording_active() -> bool:
    """Return whether the process-level screen transcript currently owns stdout."""
    return _ACTIVE_RECORDERS > 0


def progress_display_console() -> Console | None:
    """Return the terminal-preserving console owned by the active recorder."""
    return _PROGRESS_DISPLAY_CONSOLE


def progress_summary_console() -> Console | None:
    """Return the plain screen-log console paired with live terminal progress."""
    return _PROGRESS_SUMMARY_CONSOLE


def _write_all(file_descriptor: int, data: bytes) -> None:
    remaining = memoryview(data)
    while remaining:
        written = os.write(file_descriptor, remaining)
        if written <= 0:
            raise OSError("screen transcript write returned no progress")
        remaining = remaining[written:]


class _LockedTranscriptStream:
    """Write progress summaries atomically beside the recorder worker."""

    encoding = "utf-8"

    def __init__(self, file_descriptor: int, lock: threading.Lock):
        self.file_descriptor = file_descriptor
        self.lock = lock

    def write(self, value: str) -> int:
        data = str(value).encode(self.encoding, errors="replace")
        with self.lock:
            _write_all(self.file_descriptor, data)
        return len(value)

    def flush(self) -> None:
        return None

    def isatty(self) -> bool:
        return False


def _option_value(arguments: Sequence[str], option: str) -> str | None:
    """Return the last explicit option value without pre-empting real parsing."""
    value = None
    prefix = option + "="
    for index, argument in enumerate(arguments):
        if argument.startswith(prefix):
            value = argument[len(prefix):]
        elif argument == option and index + 1 < len(arguments):
            value = arguments[index + 1]
    return value


def _screen_override(arguments: Sequence[str]) -> bool | None:
    return False if "--hide-screen" in arguments else None


def resolve_screen_settings(
    command: str,
    arguments: Sequence[str],
) -> ScreenSettings | None:
    """Resolve screen settings before importing or running a scientific module."""
    if (
        command in _NON_SCIENTIFIC_COMMANDS
        or not arguments
        or "--help" in arguments
        or "-h" in arguments
    ):
        return None

    config_file = _option_value(arguments, "--run-config")
    show_screen = _screen_override(arguments)
    global_overrides = {}
    if show_screen is not None:
        global_overrides["logging.show_screen"] = show_screen

    resolved_configuration = None
    try:
        if command == "pipeline":
            configuration = load_configuration(
                config_file,
                cli_overrides=global_overrides,
            )
            module_output = None
        else:
            canonical_module = canonical_module_name(command)
            configuration = load_run_configuration_for_module(
                canonical_module,
                config_file,
                global_overrides=global_overrides,
            )
            module = getattr(configuration.modules, canonical_module)
            module_output = getattr(module, "output_directory", None)
        resolved_configuration = configuration
    except Exception:
        # The real command owns configuration error reporting. Start recording
        # with packaged presentation defaults so that failure is still saved.
        configuration = load_configuration(cli_overrides=global_overrides)
        module_output = None

    output_override = _option_value(arguments, "--output-directory")
    output_directory = Path(
        output_override
        or module_output
        or configuration.run.output_directory
    ).expanduser()
    return ScreenSettings(
        show_screen=configuration.logging.show_screen,
        log_file=output_directory / configuration.logging.screen_log_file,
        configuration=resolved_configuration,
        output_directory=output_directory,
    )


def _copy_pipe(
    read_fd: int,
    stop_fd: int,
    display_fd: int,
    log_fd: int,
    show_screen: bool,
    log_lock: threading.Lock,
    display_console: Console | None,
    display_encoding: str,
) -> None:
    """Copy one redirected stream until EOF or an explicit recorder stop.

    Multiprocessing helper processes can inherit file descriptors 1 and 2 and
    remain alive until their parent exits.  Waiting only for pipe EOF therefore
    deadlocks recorder shutdown: the parent waits for this worker while the
    helper waits for the parent.  A private, non-inherited control pipe lets the
    parent request shutdown after it has flushed and restored the real streams.
    The worker drains every byte already queued before returning.
    """
    decoder = TerminalTextDecoder(display_encoding)
    pending = ""

    def copy(data: bytes, *, final: bool = False) -> None:
        nonlocal pending
        plain = decoder.decode(data, final=final)
        if plain:
            with log_lock:
                _write_all(log_fd, plain.encode("utf-8"))
        if not show_screen:
            return
        if display_console is None:
            if plain:
                _write_all(display_fd, plain.encode(display_encoding, errors="replace"))
            return
        pending += plain
        # Complete lines keep semantic labels and UTF-8 intact across pipe
        # reads. An unusually long native line is streamed, never accumulated
        # without a bound or truncated. No word-based severity guessing.
        boundary = pending.rfind("\n") + 1
        if final or len(pending) >= _PIPE_READ_BYTES:
            boundary = len(pending)
        if boundary:
            display_console.print(style_screen_block(pending[:boundary]), end="", soft_wrap=True)
            pending = pending[boundary:]

    try:
        while True:
            readable, _, _ = select.select((read_fd, stop_fd), (), ())
            if read_fd in readable:
                block = os.read(read_fd, _PIPE_READ_BYTES)
                if not block:
                    break
                copy(block)
            if stop_fd in readable:
                os.read(stop_fd, 1)
                while select.select((read_fd,), (), (), 0)[0]:
                    block = os.read(read_fd, _PIPE_READ_BYTES)
                    if not block:
                        break
                    copy(block)
                break
        copy(b"", final=True)
    finally:
        os.close(read_fd)
        os.close(stop_fd)


@contextmanager
def record_screen(settings: ScreenSettings | None) -> Iterator[None]:
    """Apply one resolved presentation to every public direct/pipeline stream."""
    logging_config = (
        settings.configuration.logging
        if settings is not None and settings.configuration is not None
        else default_terminal_settings()
    )
    with terminal_presentation(logging_config), _record_screen_streams(settings):
        yield


@contextmanager
def _record_screen_streams(settings: ScreenSettings | None) -> Iterator[None]:
    """Tee process stdout/stderr to one append-only screen log."""
    global _ACTIVE_RECORDERS, _PROGRESS_DISPLAY_CONSOLE
    global _PROGRESS_DISPLAY_STREAM, _PROGRESS_SUMMARY_CONSOLE

    if settings is None:
        with nullcontext():
            yield
        return

    settings.log_file.parent.mkdir(parents=True, exist_ok=True)
    log_fd = os.open(
        settings.log_file,
        os.O_APPEND | os.O_CREAT | os.O_WRONLY,
        0o644,
    )
    streams = (sys.stdout, sys.stderr)
    for stream in streams:
        stream.flush()

    original_fds = (os.dup(1), os.dup(2))
    pipes = (os.pipe(), os.pipe())
    stop_pipes = (os.pipe(), os.pipe())
    log_lock = threading.Lock()
    workers = []
    recording_active = False
    previous_progress_console = _PROGRESS_DISPLAY_CONSOLE
    previous_progress_stream = _PROGRESS_DISPLAY_STREAM
    previous_summary_console = _PROGRESS_SUMMARY_CONSOLE
    try:
        _PROGRESS_DISPLAY_CONSOLE = None
        _PROGRESS_DISPLAY_STREAM = None
        _PROGRESS_SUMMARY_CONSOLE = None
        if settings.show_screen and os.isatty(original_fds[0]):
            _PROGRESS_DISPLAY_STREAM = os.fdopen(
                os.dup(original_fds[0]),
                "w",
                buffering=1,
                encoding=getattr(sys.stdout, "encoding", None) or "utf-8",
                errors="replace",
            )
            _PROGRESS_DISPLAY_CONSOLE = Console(
                file=_PROGRESS_DISPLAY_STREAM,
                force_terminal=True,
                color_system=(None if terminal_settings().terminal_style.color == "never" else "auto"),
                highlight=False,
                markup=False,
            )
            _PROGRESS_SUMMARY_CONSOLE = Console(
                file=_LockedTranscriptStream(log_fd, log_lock),
                force_terminal=False,
                color_system=None,
            )
        stream_descriptors = zip((1, 2), pipes, stop_pipes, original_fds)
        for (
            target_fd,
            (read_fd, write_fd),
            (stop_fd, stop_write_fd),
            display_fd,
        ) in stream_descriptors:
            os.dup2(write_fd, target_fd)
            os.close(write_fd)
            worker = threading.Thread(
                target=copy_context().run,
                args=(
                    _copy_pipe,
                    read_fd,
                    stop_fd,
                    display_fd,
                    log_fd,
                    settings.show_screen,
                    log_lock,
                    _PROGRESS_DISPLAY_CONSOLE,
                    (
                        getattr(_PROGRESS_DISPLAY_STREAM, "encoding", None)
                        or "utf-8"
                    ),
                ),
                daemon=True,
            )
            worker.start()
            workers.append(worker)
        _ACTIVE_RECORDERS += 1
        recording_active = True
        yield
    finally:
        if recording_active:
            _ACTIVE_RECORDERS -= 1
        for stream in streams:
            try:
                stream.flush()
            except Exception:
                pass
        for target_fd, original_fd in zip((1, 2), original_fds):
            os.dup2(original_fd, target_fd)
        for _, stop_write_fd in stop_pipes:
            try:
                _write_all(stop_write_fd, b"\0")
            except BrokenPipeError:
                # The copy worker already observed normal EOF and exited.
                pass
            os.close(stop_write_fd)
        for worker in workers:
            worker.join()
        active_progress_stream = _PROGRESS_DISPLAY_STREAM
        _PROGRESS_DISPLAY_CONSOLE = previous_progress_console
        _PROGRESS_DISPLAY_STREAM = previous_progress_stream
        _PROGRESS_SUMMARY_CONSOLE = previous_summary_console
        if active_progress_stream is not None and (
            active_progress_stream is not previous_progress_stream
        ):
            active_progress_stream.close()
        for original_fd in original_fds:
            os.close(original_fd)
        os.fsync(log_fd)
        os.close(log_fd)


__all__ = [
    "ScreenSettings",
    "record_screen",
    "progress_display_console",
    "progress_summary_console",
    "resolve_screen_settings",
    "screen_recording_active",
]
