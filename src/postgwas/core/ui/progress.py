"""Reusable live stage progress for interactive scientific runs."""

from __future__ import annotations

import time

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
        self._progress: Progress | None = None
        self._task_id = None
        self._step_started: float | None = None
        self._closed = False

    def _start(self, number: int, total: int, title: str) -> None:
        self.console.print()
        self.console.print(Text("  🔬  %s" % self.label, style="bold cyan"))
        self._progress = Progress(
            TextColumn("      🔬  {task.description}"),
            BarColumn(),
            MofNCompleteColumn(),
            TaskProgressColumn(),
            TimeElapsedColumn(),
            console=self.console,
        )
        self._task_id = self._progress.add_task(
            self._description("Current", number, total, title),
            total=total,
            completed=number - 1,
        )
        self._progress.start()

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
        if self._progress is None:
            self._start(number, total, title)
        else:
            self._progress.update(
                self._task_id,
                description=self._description("Current", number, total, title),
                total=total,
                completed=number - 1,
            )
            self._progress.refresh()
        self._step_started = time.monotonic()

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
        if self._progress is None:
            self._start(number, total, str(title))
        elapsed = time.monotonic() - (self._step_started or time.monotonic())
        self._progress.update(
            self._task_id,
            description=self._description("Completed", number, total, str(title)),
            total=total,
            completed=number,
        )
        self._progress.refresh()
        self.console.print(Text(
            "      ✅  Completed %d/%d · %s · %s"
            % (number, total, title, self._duration(elapsed)),
            style="green",
        ))
        if outcome_fields:
            self.console.print(Text("          🧮  Outcome", style="bold cyan"))
            label_width = self.outcome_label_width or max(
                len(label) for _, label, _ in outcome_fields
            )
            for kind, label, value in outcome_fields:
                if isinstance(value, int) and not isinstance(value, bool):
                    value = f"{value:,}"
                self.console.print(Text(
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
            self.console.print(Text(
                "          🧮  Outcome : %s" % str(outcome).strip(),
                style="cyan",
            ))
        if number == total:
            self._progress.update(
                self._task_id,
                description="All %d stages completed" % total,
                completed=total,
            )
            self.close()

    def fail_step(self, number: int, total: int, title: str) -> None:
        if not self.enabled:
            return
        number, total = self._values(number, total)
        if self._progress is None:
            self._start(number, total, str(title))
        elapsed = time.monotonic() - (self._step_started or time.monotonic())
        self._progress.update(
            self._task_id,
            description=self._description("Failed", number, total, str(title)),
            total=total,
            completed=number - 1,
        )
        self._progress.refresh()
        self.console.print(Text(
            "      ❌  Failed %d/%d · %s · %s"
            % (number, total, title, self._duration(elapsed)),
            style="bold red",
        ))
        self.close()

    def close(self) -> None:
        if self._closed:
            return
        if self._progress is not None:
            self._progress.stop()
        self._closed = True


__all__ = ["StageProgress"]
