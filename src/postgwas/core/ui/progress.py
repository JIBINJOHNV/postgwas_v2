"""Reusable live stage progress for interactive scientific runs."""

from __future__ import annotations

import time
from collections.abc import Callable
from contextlib import contextmanager
from contextvars import ContextVar
from typing import TypeVar

from rich.console import Console
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.text import Text

from postgwas.core.ui.screen import screen_field


_Result = TypeVar("_Result")
_ACTIVE_PROGRESS: ContextVar[tuple["StageProgress", ...]] = ContextVar(
    "postgwas_active_progress", default=(),
)


class StageProgress:
    """Retain completed stages and aligned, meaning-specific outcome fields."""

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

    def _emit(self, renderable=None) -> None:
        """Serialize a milestone with live display and plain transcript copies."""
        def emit(console: Console) -> None:
            if renderable is None:
                console.print()
            else:
                console.print(renderable)

        if self._display_console is None:
            emit(self.console)
            return
        emit(self._display_console)
        if self._summary_console is not None:
            emit(self._summary_console)

    def _emit_started(self, renderable) -> None:
        """Record a durable start without leaving a static live-stage copy."""
        if self._summary_console is not None:
            self._summary_console.print(renderable)
        elif not self._interactive:
            self.console.print(renderable)

    def _pause(self) -> None:
        if self._progress is not None:
            self._progress.stop()

    def _resume(self) -> None:
        if self._progress is not None and not self._closed:
            self._progress.start()

    def _start(self, number: int, total: int, title: str) -> None:
        active = _ACTIVE_PROGRESS.get()
        if active:
            self._parent = active[-1]
            self._parent._pause()
        self._started = True
        self._emit()
        self._emit(Text("  🔬  %s" % self.label, style="bold cyan"))
        self._emit_started(
            Text(
                "      🔬  %s"
                % self._description("Started", number, total, title),
                style="cyan",
            ),
        )
        if self._interactive:
            self._progress = Progress(
                TextColumn("      🔬  {task.description}"),
                BarColumn(),
                MofNCompleteColumn(),
                TaskProgressColumn(),
                TimeElapsedColumn(),
                console=self._live_console,
                transient=True,
                redirect_stdout=False,
                redirect_stderr=False,
            )
            self._task_id = self._progress.add_task(
                self._description("Current", number, total, title),
                total=total,
                completed=number - 1,
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
                Text(
                    "      🔬  %s"
                    % self._description("Started", number, total, title),
                    style="cyan",
                ),
            )
            if self._progress is not None:
                self._progress.update(
                    self._task_id,
                    description=self._description("Current", number, total, title),
                    total=total,
                    completed=number - 1,
                )
            self._resume()
        self._step_started = time.monotonic()

    @contextmanager
    def step(self, number: int, total: int, title: str):
        """Run one measurable stage and preserve truthful failure progress."""
        self.start_step(number, total, title)
        try:
            yield
        except BaseException:
            self.fail_step(number, total, title)
            raise
        else:
            self.complete_step(number, total, title)

    def complete_step(
        self,
        number: int,
        total: int,
        title: str,
        outcome: str | None = None,
        outcome_fields=None,
    ) -> None:
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
        self._emit(Text(
            "      ✅  Completed %d/%d · %s · %s"
            % (number, total, title, self._duration(elapsed)),
            style="green",
        ))
        if outcome_fields:
            self._emit(Text("          🧮  Outcome", style="bold cyan"))
            label_width = self.outcome_label_width or max(
                len(label) for _, label, _ in outcome_fields
            )
            for kind, label, value in outcome_fields:
                if isinstance(value, int) and not isinstance(value, bool):
                    value = f"{value:,}"
                self._emit(Text(
                    screen_field(
                        kind,
                        label,
                        value,
                        width=self.console.width,
                        indent=14,
                        label_width=label_width,
                    ),
                    style="cyan",
                ))
        elif outcome:
            self._emit(Text(
                "          🧮  Outcome : %s" % str(outcome).strip(),
                style="cyan",
            ))
        if number == total:
            self._emit(Text(
                "      ✅  All %d stages completed" % total,
                style="green",
            ))
            self.close()
        else:
            self._resume()

    def fail_step(self, number: int, total: int, title: str) -> None:
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
        self._emit(Text(
            "      ❌  Failed %d/%d · %s · %s"
            % (number, total, title, self._duration(elapsed)),
            style="bold red",
        ))
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
        self._emit(Text("  🔬  %s" % self.label, style="bold cyan"))
        started = Text(
            "      🔬  %s"
            % self._count_description("Started", 0, total, self._title),
            style="cyan",
        )
        self._emit_started(started)
        if self._interactive:
            self._progress = Progress(
                TextColumn("      🔬  {task.description}"),
                BarColumn(),
                MofNCompleteColumn(),
                TaskProgressColumn(),
                TimeElapsedColumn(),
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
        self._emit(Text(
            "      🔬  Progress %s/%s · %s · %d%% · %s"
            % (
                f"{displayed:,}",
                f"{total:,}",
                self._title,
                percentage,
                self._duration(elapsed),
            ),
            style="cyan",
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
        current = Text(
            "      🔬  %s"
            % self._count_description(
                "Current", displayed, self._total, self._title,
            ),
            style="cyan",
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
        self._emit(Text(
            "      ✅  Completed %s/%s · %s · %s"
            % (f"{completed:,}", f"{total:,}", self._title, self._duration(elapsed)),
            style="green",
        ))
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
        self._emit(Text(
            "      ❌  Failed %s · %s"
            % (
                self._count_description("", displayed, total, self._title).strip(),
                self._duration(elapsed),
            ),
            style="bold red",
        ))
        self.close(resume_parent=False)


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


__all__ = ["MeasuredProgress", "StageProgress", "run_with_progress"]
