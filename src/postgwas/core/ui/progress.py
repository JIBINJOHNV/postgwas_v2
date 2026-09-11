"""Reusable live stage progress for interactive scientific runs."""

from __future__ import annotations

import time
from collections.abc import Callable
from contextlib import contextmanager
from contextvars import ContextVar
from typing import TypeVar

from rich.cells import cell_len
from rich.console import Console
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    ProgressColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.text import Text
from rich.table import Column

from postgwas.core.ui.screen import (
    SYMBOLS,
    normalize_stage_outcome_fields,
    screen_field,
    screen_line,
    style_screen_block,
    terminal_style,
)


_Result = TypeVar("_Result")
_ACTIVE_PROGRESS: ContextVar[tuple["StageProgress", ...]] = ContextVar(
    "postgwas_active_progress", default=(),
)


class _NeutralProgressColumn(ProgressColumn):
    """Retain Rich's count/time algorithms without its independent palette."""

    def __init__(self, column: ProgressColumn):
        super().__init__(table_column=Column(no_wrap=True))
        self.column = column
        # Rich refreshes in its own thread, outside the run ContextVar scope.
        self.style = terminal_style("text")

    def render(self, task):
        return Text(self.column.render(task).plain, style=self.style)


class _OperationColumn(TextColumn):
    """Show the active stage number separately from completed-stage progress."""

    def __init__(self):
        super().__init__(
            "{task.description}", style=terminal_style("analysis"), markup=False,
            table_column=Column(ratio=2, no_wrap=True, overflow="ellipsis"),
        )
        self.prefix = screen_line("analysis", "", indent=6)

    def render(self, task):
        # The numeric progress columns count only validated completed work.
        # Retain the active stage number in a compact form so it is not
        # mistaken for the completed-stage count shown beside the bar, while
        # preserving enough width for the operation name on narrow terminals.
        _counter, separator, operation = task.description.partition(" · ")
        description = operation if separator else task.description
        active_stage = task.fields.get("active_stage")
        active_stage_total = task.fields.get("active_stage_total")
        if active_stage is not None and active_stage_total is not None:
            description = "%s/%s · %s" % (
                active_stage,
                active_stage_total,
                description,
            )
        return Text(self.prefix + description, style=self.style)


class StageOutcome:
    """Mutable completion details yielded by :meth:`StageProgress.step`."""

    __slots__ = ("message", "fields")

    def __init__(self) -> None:
        self.message: str | None = None
        self.fields: tuple[tuple, ...] | None = None

    def outcome(self, message: str, fields=None, **_values) -> None:
        """Attach a concise message and semantic fields to a completed stage."""
        text = str(message).strip()
        if not text:
            raise ValueError("A completed-stage outcome must not be empty")
        self.message = text
        self.fields = normalize_stage_outcome_fields(fields)


class StageProgress:
    """Retain completed stages and aligned, meaning-specific outcome fields."""

    _OUTCOME_FIELD_INDENT = 14
    _OUTCOME_SECTION_FIELD_INDENT = 18

    def __init__(
        self,
        label: str,
        *,
        enabled: bool,
        outcome_label_width: int | None = None,
        console: Console | None = None,
    ):
        self.label = str(label)
        self.enabled = bool(enabled)
        if outcome_label_width is not None and int(outcome_label_width) < 1:
            raise ValueError("outcome_label_width must be a positive integer")
        self.outcome_label_width = (
            None if outcome_label_width is None else int(outcome_label_width)
        )
        self.console = console or Console()
        from postgwas.core.screen_logging import (
            progress_display_console,
            progress_summary_console,
        )

        display_console = progress_display_console()
        self._display_console = display_console
        self._summary_console = (
            progress_summary_console() if display_console is not None else None
        )
        self._live_console = display_console or self.console
        # Rich can redraw a live region only when it owns a real terminal.
        # Captured pipes and durable screen logs must receive deterministic
        # milestones; otherwise every refresh becomes another permanent line.
        self._interactive = bool(
            display_console is not None or self.console.is_terminal
        )
        self._progress: Progress | None = None
        self._task_id = None
        self._step_started: float | None = None
        self._started = False
        self._closed = False
        self._parent: StageProgress | None = None

    @property
    def outcome_separator_column(self) -> int | None:
        """Return the configured terminal-cell column before `` : ``."""
        if self.outcome_label_width is None:
            return None
        return (
            self._OUTCOME_FIELD_INDENT
            + cell_len(SYMBOLS["count"])
            + 2
            + self.outcome_label_width
        )

    def _emit(self, renderable=None, *, soft_wrap: bool = False) -> None:
        """Serialize a milestone with live display and plain transcript copies."""
        if isinstance(renderable, (str, Text)):
            renderable = style_screen_block(renderable)

        def emit(console: Console) -> None:
            if renderable is None:
                console.print()
            else:
                console.print(renderable, soft_wrap=soft_wrap)

        if self._display_console is None:
            emit(self.console)
            return
        emit(self._display_console)
        if self._summary_console is not None:
            emit(self._summary_console)

    def _emit_started(self, renderable) -> None:
        """Record a durable start without leaving a static live-stage copy."""
        if self._summary_console is not None:
            self._summary_console.print(renderable, soft_wrap=True)
        elif not self._interactive:
            self.console.print(renderable, soft_wrap=True)

    def _emit_forced_capture_current(self, description: str) -> None:
        """Preserve one live state when a forced terminal is not a real TTY."""
        if (
            self._display_console is None
            and self._interactive
            and not self.console.file.isatty()
        ):
            self.console.print(
                style_screen_block(screen_line("analysis", description, indent=6)),
                end="\r",
            )

    def print_block(self, block: str | Text, *, console: Console | None = None) -> None:
        """Print one multiline result block without a live progress redraw."""
        active = _ACTIVE_PROGRESS.get()
        current = active[-1] if active else None
        if current is not None:
            current._pause()
        try:
            if console is not None and self._display_console is None:
                console.print(style_screen_block(block), soft_wrap=True)
            else:
                self._emit(block, soft_wrap=True)
        finally:
            if current is not None and not current._closed:
                current._resume()

    def _pause(self) -> None:
        if self._progress is not None:
            self._progress.stop()

    @staticmethod
    def _columns():
        """The same theme for stage-count and native measured progress bars."""
        return (
            _OperationColumn(),
            BarColumn(
                bar_width=None,
                table_column=Column(ratio=1, min_width=1),
                style=terminal_style("text"),
                complete_style=terminal_style("analysis"),
                finished_style=terminal_style("success"),
                pulse_style=terminal_style("analysis"),
            ),
            _NeutralProgressColumn(MofNCompleteColumn()),
            _NeutralProgressColumn(TaskProgressColumn()),
            _NeutralProgressColumn(TimeElapsedColumn()),
        )

    def _resume(self) -> None:
        active = _ACTIVE_PROGRESS.get()
        # Rich's live stack must retain the same parent-to-child order as this
        # logical stack. Restarting an ancestor before its active child closes
        # reverses that order and makes Rich recursively refresh its root Live.
        has_active_child = any(item is self for item in active) and (
            active[-1] is not self
        )
        if self._progress is not None and not self._closed and not has_active_child:
            self._progress.start()

    def _start(self, number: int, total: int, title: str) -> None:
        active = _ACTIVE_PROGRESS.get()
        if active:
            self._parent = active[-1]
            self._parent._pause()
        self._started = True
        self._emit()
        self._emit(screen_line("analysis", self.label, indent=2))
        self._emit_started(
            style_screen_block(
                screen_line("analysis", "%s", indent=6)
                % self._description("Started", number, total, title),
            ),
        )
        if (
            self._display_console is not None
            and not self._display_console.file.isatty()
        ):
            # A terminal-preserving recorder can target a captured descriptor
            # in tests or redirected launchers. Rich cannot retain its first
            # transient redraw there, so emit the truthful active state once.
            self._display_console.print(style_screen_block(
                screen_line("analysis", "%s", indent=6)
                % self._description("Current", number, total, title),
            ))
        self._emit_forced_capture_current(
            self._description("Current", number, total, title),
        )
        if self._interactive:
            self._progress = Progress(
                *self._columns(),
                expand=True,
                console=self._live_console,
                transient=True,
                redirect_stdout=False,
                redirect_stderr=False,
            )
            self._task_id = self._progress.add_task(
                self._description("Current", number, total, title),
                total=total,
                completed=number - 1,
                active_stage=number,
                active_stage_total=total,
            )
            self._progress.start()
            self._progress.refresh()
        _ACTIVE_PROGRESS.set((*active, self))

    @staticmethod
    def _values(number: int, total: int) -> tuple[int, int]:
        number, total = int(number), int(total)
        if total < 1 or number < 1 or number > total:
            raise ValueError(
                "Invalid stage progress: step=%s total=%s" % (number, total)
            )
        return number, total

    @staticmethod
    def _description(status: str, number: int, total: int, title: str) -> str:
        return "%s %d/%d · %s" % (status, number, total, title)

    @staticmethod
    def _duration(seconds: float) -> str:
        if seconds < 60:
            return "%.1f s" % seconds
        return "%d m %02d s" % (int(seconds // 60), int(seconds % 60))

    def start_step(self, number: int, total: int, title: str) -> None:
        if not self.enabled:
            return
        number, total = self._values(number, total)
        title = str(title)
        if not self._started:
            self._start(number, total, title)
        else:
            self._pause()
            self._emit_started(
                style_screen_block(
                    screen_line("analysis", "%s", indent=6)
                    % self._description("Started", number, total, title),
                ),
            )
            self._emit_forced_capture_current(
                self._description("Current", number, total, title),
            )
            if self._progress is not None:
                self._progress.update(
                    self._task_id,
                    description=self._description("Current", number, total, title),
                    total=total,
                    completed=number - 1,
                    active_stage=number,
                    active_stage_total=total,
                )
            self._resume()
        self._step_started = time.monotonic()

    @contextmanager
    def step(self, number: int, total: int, title: str):
        """Run one measurable stage and preserve truthful completion details."""
        outcome = StageOutcome()
        self.start_step(number, total, title)
        try:
            yield outcome
        except BaseException:
            self.fail_step(number, total, title)
            raise
        else:
            self.complete_step(
                number,
                total,
                title,
                outcome=outcome.message,
                outcome_fields=outcome.fields,
            )

    def complete_step(
        self,
        number: int,
        total: int,
        title: str,
        outcome: str | None = None,
        outcome_fields=None,
    ) -> None:
        from postgwas.core.validation_reporting import (
            consolidate_validation_fields, flush_file_validation_display,
        )

        outcome_fields = consolidate_validation_fields(outcome_fields)
        flush_file_validation_display()
        if not self.enabled:
            return
        number, total = self._values(number, total)
        if not self._started:
            self._start(number, total, str(title))
        elapsed = time.monotonic() - (self._step_started or time.monotonic())
        if self._progress is not None:
            self._progress.update(
                self._task_id,
                description=self._description("Completed", number, total, str(title)),
                total=total,
                completed=number,
            )
            self._progress.refresh()
        self._pause()
        self._emit(screen_line("success", "Completed %d/%d · %s · %s", indent=6)
            % (number, total, title, self._duration(elapsed)), soft_wrap=True)
        if outcome_fields:
            self._emit(screen_line("count", "Outcome", indent=10))
            value_fields = [field for field in outcome_fields if len(field) == 3]
            label_width = self.outcome_label_width or max(
                (len(field[1]) for field in value_fields), default=1,
            )
            separator_column = self.outcome_separator_column or (
                self._OUTCOME_FIELD_INDENT
                + cell_len(SYMBOLS["count"])
                + 2
                + label_width
            )
            section_started = False
            for field in outcome_fields:
                kind, label = field[:2]
                if len(field) == 2:
                    if section_started:
                        self._emit("")
                    self._emit(screen_line(
                            kind,
                            label,
                            indent=self._OUTCOME_FIELD_INDENT,
                        ))
                    section_started = True
                    continue
                value = field[2]
                if isinstance(value, int) and not isinstance(value, bool):
                    value = f"{value:,}"
                self._emit(screen_field(
                        kind,
                        label,
                        value,
                        width=self.console.width,
                        indent=(
                            self._OUTCOME_SECTION_FIELD_INDENT
                            if section_started
                            else self._OUTCOME_FIELD_INDENT
                        ),
                        label_width=label_width,
                        separator_column=separator_column,
                        break_long_values=True,
                    ), soft_wrap=True)
        elif outcome:
            self._emit(screen_field("count", "Outcome", str(outcome).strip(),
                             indent=10, width=self._live_console.width), soft_wrap=True)
        if number == total:
            self._emit(screen_line("success", "All %d stages completed" % total, indent=6))
            self.close()
        else:
            self._resume()

    def fail_step(self, number: int, total: int, title: str) -> None:
        from postgwas.core.validation_reporting import flush_file_validation_display

        flush_file_validation_display()
        if not self.enabled:
            return
        number, total = self._values(number, total)
        if not self._started:
            self._start(number, total, str(title))
        elapsed = time.monotonic() - (self._step_started or time.monotonic())
        if self._progress is not None:
            self._progress.update(
                self._task_id,
                description=self._description("Failed", number, total, str(title)),
                total=total,
                completed=number - 1,
            )
            self._progress.refresh()
        self._pause()
        self._emit(screen_line("error", "Failed %d/%d · %s · %s", indent=6)
            % (number, total, title, self._duration(elapsed)), soft_wrap=True)
        # A failed child stage propagates to its enclosing operation. Keep the
        # parent live region paused so the actionable error can be printed
        # without an obsolete outer bar redrawing over it.
        self.close(resume_parent=False)

    def close(self, *, resume_parent: bool = True) -> None:
        if self._closed:
            return
        if self._progress is not None:
            self._progress.stop()
        self._closed = True
        active = _ACTIVE_PROGRESS.get()
        if self in active:
            _ACTIVE_PROGRESS.set(tuple(item for item in active if item is not self))
        if resume_parent and self._parent is not None:
            self._parent._resume()


class PipelineStageController:
    """Advance one validated scientific stage plan across pipeline modules."""

    def __init__(
        self,
        label: str,
        stages,
        *,
        console: Console | None = None,
        outcome_label_width: int | None = None,
    ):
        titles = tuple(str(title) for title in stages)
        if not titles:
            raise ValueError("Pipeline progress requires at least one stage")
        self.stages = titles
        self.progress = StageProgress(
            label,
            enabled=True,
            console=console,
            outcome_label_width=outcome_label_width,
        )
        self.current = 0
        self.completed = 0

    def start(self, number: int) -> None:
        number = int(number)
        if number != self.completed + 1:
            raise ValueError(
                "Pipeline stage %d cannot start after completed stage %d"
                % (number, self.completed)
            )
        self.current = number
        self.progress.start_step(number, len(self.stages), self.stages[number - 1])

    def complete(self, number: int, *, outcome=None, outcome_fields=None) -> None:
        number = int(number)
        if self.current != number:
            raise ValueError(
                "Pipeline stage %d cannot complete while stage %d is active"
                % (number, self.current)
            )
        self.progress.complete_step(
            number,
            len(self.stages),
            self.stages[number - 1],
            outcome=outcome,
            outcome_fields=outcome_fields,
        )
        self.completed = number
        self.current = 0

    def fail_active(self) -> None:
        if self.current:
            self.progress.fail_step(
                self.current,
                len(self.stages),
                self.stages[self.current - 1],
            )
            self.current = 0

    def close(self) -> None:
        self.progress.close()


class MeasuredProgress(StageProgress):
    """Display progress derived from observed external-tool work units.

    Intermediate updates are capped below completion.  Only ``complete`` may
    render 100%, after the caller has validated the command and its output.
    """

    def __init__(self, label: str, *, enabled: bool, console: Console | None = None):
        super().__init__(label, enabled=enabled, console=console)
        self._title = ""
        self._completed = 0
        self._total: int | None = None
        self._last_summary_percentage: int | None = None

    @staticmethod
    def _count_values(
        completed: int,
        total: int | None,
    ) -> tuple[int, int | None]:
        completed = int(completed)
        total = None if total is None else int(total)
        if completed < 0 or (total is not None and total < 1):
            raise ValueError(
                "Invalid measured progress: completed=%s total=%s"
                % (completed, total)
            )
        if total is not None and completed > total:
            raise ValueError(
                "Measured progress exceeds its total: completed=%s total=%s"
                % (completed, total)
            )
        return completed, total

    @staticmethod
    def _count_description(
        status: str,
        completed: int,
        total: int | None,
        title: str,
    ) -> str:
        denominator = "?" if total is None else f"{total:,}"
        return "%s %s/%s · %s" % (
            status,
            f"{completed:,}",
            denominator,
            title,
        )

    def start(self, title: str, *, total: int | None = None) -> None:
        if not self.enabled:
            return
        if self._started:
            raise RuntimeError("Measured progress has already started")
        _, total = self._count_values(0, total)
        active = _ACTIVE_PROGRESS.get()
        if active:
            self._parent = active[-1]
            self._parent._pause()
        self._started = True
        self._step_started = time.monotonic()
        self._title = str(title)
        self._total = total
        self._emit()
        self._emit(screen_line("analysis", self.label, indent=2))
        started = style_screen_block(
            screen_line("analysis", "%s", indent=6)
            % self._count_description("Started", 0, total, self._title),
        )
        self._emit_started(started)
        self._emit_forced_capture_current(
            self._count_description("Current", 0, total, self._title),
        )
        if self._interactive:
            self._progress = Progress(
                *self._columns(),
                expand=True,
                console=self._live_console,
                transient=True,
                redirect_stdout=False,
                redirect_stderr=False,
            )
            self._task_id = self._progress.add_task(
                self._count_description("Current", 0, total, self._title),
                total=total,
                completed=0,
            )
            self._progress.start()
            self._progress.refresh()
        _ACTIVE_PROGRESS.set((*active, self))

    def update(
        self,
        completed: int,
        *,
        total: int | None = None,
        title: str | None = None,
    ) -> None:
        if not self.enabled:
            return
        completed, total = self._count_values(completed, total)
        if not self._started:
            self.start(title or self.label, total=total)
        if total is None:
            total = self._total
        elif self._total is not None and total != self._total:
            raise ValueError(
                "Measured progress total changed from %s to %s"
                % (self._total, total)
            )
        if completed < self._completed:
            raise ValueError(
                "Measured progress decreased from %s to %s"
                % (self._completed, completed)
            )
        self._total = total
        self._completed = completed
        if title is not None:
            self._title = str(title)

        # A complete native file does not prove that the external command and
        # output validation succeeded. Keep the live value below 100% until
        # complete() is called by the validated execution path.
        displayed = completed
        if total is not None and completed == total:
            displayed = max(0, total - 1)
        if self._progress is not None:
            self._progress.update(
                self._task_id,
                description=self._count_description(
                    "Current", displayed, total, self._title,
                ),
                total=total,
                completed=displayed,
            )
            self._progress.refresh()

        if total is None:
            return
        percentage = int(100 * displayed / total)
        if percentage == self._last_summary_percentage:
            return
        self._last_summary_percentage = percentage
        elapsed = time.monotonic() - (self._step_started or time.monotonic())
        self._emit(screen_line("analysis", "Progress %s/%s · %s · %d%% · %s", indent=6)
            % (
                f"{displayed:,}",
                f"{total:,}",
                self._title,
                percentage,
                self._duration(elapsed),
            ))

    def set_phase(self, title: str) -> None:
        """Change an observed phase without pretending that work was measured."""
        if not self.enabled:
            return
        title = str(title)
        if not self._started:
            self.start(title, total=self._total)
            return
        if title == self._title:
            return
        self._title = title
        displayed = self._completed
        if self._total is not None and displayed == self._total:
            displayed = max(0, self._total - 1)
        current = style_screen_block(
            screen_line("analysis", "%s", indent=6)
            % self._count_description(
                "Current", displayed, self._total, self._title,
            ),
        )
        self._emit_started(current)
        if self._progress is not None:
            self._progress.update(
                self._task_id,
                description=self._count_description(
                    "Current", displayed, self._total, self._title,
                ),
                total=self._total,
                completed=displayed,
            )
            self._progress.refresh()

    def complete(self, completed: int, *, total: int, title: str | None = None) -> None:
        if not self.enabled:
            return
        completed, total = self._count_values(completed, total)
        if completed != total:
            raise ValueError(
                "Completed measured progress must equal its total: %s/%s"
                % (completed, total)
            )
        if not self._started:
            self.start(title or self.label, total=total)
        if title is not None:
            self._title = str(title)
        elapsed = time.monotonic() - (self._step_started or time.monotonic())
        if self._progress is not None:
            self._progress.update(
                self._task_id,
                description=self._count_description(
                    "Completed", completed, total, self._title,
                ),
                total=total,
                completed=completed,
            )
            self._progress.refresh()
        self._pause()
        self._emit(screen_line("success", "Completed %s/%s · %s · %s", indent=6)
            % (f"{completed:,}", f"{total:,}", self._title, self._duration(elapsed)))
        self.close()

    def fail(self, *, title: str | None = None) -> None:
        if not self.enabled:
            return
        if not self._started:
            self.start(title or self.label)
        if title is not None:
            self._title = str(title)
        total = self._total
        displayed = self._completed
        if total is not None and displayed == total:
            displayed = max(0, total - 1)
        elapsed = time.monotonic() - (self._step_started or time.monotonic())
        if self._progress is not None:
            self._progress.update(
                self._task_id,
                description=self._count_description(
                    "Failed", displayed, total, self._title,
                ),
                total=total,
                completed=displayed,
            )
            self._progress.refresh()
        self._pause()
        self._emit(screen_line("error", "Failed %s · %s", indent=6)
            % (
                self._count_description("", displayed, total, self._title).strip(),
                self._duration(elapsed),
            ))
        self.close(resume_parent=False)


def print_screen_block(block: str | Text, *, console: Console | None = None) -> None:
    """Print a styled semantic block without corrupting active progress output."""
    active = _ACTIVE_PROGRESS.get()
    if active:
        active[-1].print_block(block, console=console)
        return
    from postgwas.core.screen_logging import progress_display_console, progress_summary_console

    rendered = style_screen_block(block)
    display = progress_display_console()
    if display is not None:
        display.print(rendered, soft_wrap=True)
        summary = progress_summary_console()
        if summary is not None:
            summary.print(rendered, soft_wrap=True)
    else:
        (console or Console()).print(rendered, soft_wrap=True)


def print_screen_message(
    kind: str, message: str, *, stderr: bool = False, console: Console | None = None,
) -> None:
    """Render explicitly classified literal text with shared wrapping and colour."""
    destination = console or Console(stderr=stderr)
    prefix = screen_line(kind, "", indent=6)
    continuation = " " * cell_len(prefix)
    # Rich Text handles terminal-cell widths and never interprets user paths
    # or external diagnostics as markup. Keep each explicit source line.
    lines = []
    for number, line in enumerate(str(message).strip("\n").split("\n")):
        first = prefix if number == 0 else continuation
        available = max(1, destination.width - cell_len(first))
        fragments = Text(line).wrap(
            destination, available, overflow="fold", no_wrap=False,
        )
        lines.extend(
            (first if index == 0 else continuation) + fragment.plain
            for index, fragment in enumerate(fragments)
        )
    print_screen_block("\n".join(lines), console=destination)


def run_with_progress(
    operation: Callable[[], _Result],
    *,
    label: str,
    title: str,
    enabled: bool,
    console: Console | None = None,
) -> _Result:
    """Run one top-level operation with truthful success or failure progress."""
    progress = StageProgress(label, enabled=enabled, console=console)
    progress.start_step(1, 1, title)
    try:
        result = operation()
    except SystemExit as exc:
        if exc.code in (None, 0):
            progress.complete_step(1, 1, title)
        else:
            progress.fail_step(1, 1, title)
        raise
    except BaseException:
        progress.fail_step(1, 1, title)
        raise
    else:
        if isinstance(result, int) and not isinstance(result, bool) and result != 0:
            progress.fail_step(1, 1, title)
        else:
            progress.complete_step(1, 1, title)
        return result
    finally:
        progress.close()


__all__ = [
    "MeasuredProgress",
    "StageProgress",
    "print_screen_block",
    "print_screen_message",
    "run_with_progress",
]
