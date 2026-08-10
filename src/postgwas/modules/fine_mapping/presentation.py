"""Concise terminal progress and final reporting for fine-mapping."""

from __future__ import annotations

from pathlib import Path

from rich.console import Console

from postgwas.core.ui import StageProgress, screen_field, screen_line


FINE_MAPPING_PROGRESS_STAGES = (
    "Initialize the fine-mapping run",
    "Validate tools and reference resources",
    "Prepare locus boundaries",
    "Prepare locus model inputs",
    "Fit the primary loci",
    "Validate and publish primary results",
    "Summarize the primary analysis",
    "Resolve overlapping loci and finalize results",
)


def _reason_count(reason, count) -> str:
    count = int(count)
    unit = "locus" if count == 1 else "loci"
    normalized = str(reason).replace("_", " ").replace(";", "; ")
    return f"{normalized} ({count} {unit})"


def _output_path(value, root: Path | None) -> str:
    path = Path(value)
    if root is not None:
        try:
            return str(path.resolve().relative_to(root))
        except ValueError:
            pass
    return str(path)


class FineMappingScreen:
    """Render ordered stage outcomes while keeping audit details in the file log."""

    def __init__(
        self,
        engine: str,
        *,
        show_progress: bool,
        label_width: int,
        console: Console | None = None,
    ):
        self.engine = str(engine)
        self.console = console or Console()
        self.progress = StageProgress(
            "Fine-mapping analysis progress",
            enabled=show_progress,
            outcome_label_width=label_width,
            console=self.console,
        )
        self.current_step = 0

    @property
    def total(self) -> int:
        return len(FINE_MAPPING_PROGRESS_STAGES)

    def start(self) -> None:
        if self.current_step:
            return
        self.current_step = 1
        self.progress.start_step(
            self.current_step,
            self.total,
            FINE_MAPPING_PROGRESS_STAGES[self.current_step - 1],
        )

    def complete(self, number: int, outcome_fields=None) -> None:
        number = int(number)
        if number != self.current_step:
            raise ValueError(
                "Fine-mapping screen stages must complete in order: "
                f"expected {self.current_step}, received {number}"
            )
        self.progress.complete_step(
            number,
            self.total,
            FINE_MAPPING_PROGRESS_STAGES[number - 1],
            outcome_fields=outcome_fields,
        )
        if number < self.total:
            self.current_step = number + 1
            self.progress.start_step(
                self.current_step,
                self.total,
                FINE_MAPPING_PROGRESS_STAGES[self.current_step - 1],
            )
        else:
            self.current_step = 0

    def fail(self) -> None:
        if not self.current_step:
            return
        self.progress.fail_step(
            self.current_step,
            self.total,
            FINE_MAPPING_PROGRESS_STAGES[self.current_step - 1],
        )
        self.current_step = 0

    def print_final_summary(self, result: dict, log_file: str | Path) -> None:
        """Print only decision-relevant counts and authoritative output paths."""
        label_width = self.progress.outcome_label_width or 42
        status = str(result.get("status", "unknown"))
        status_kind = (
            "success"
            if status == "success"
            else "warning"
            if status.startswith("completed")
            else "error"
        )
        output_root = (
            Path(result["output_dir"]).resolve()
            if result.get("output_dir")
            else None
        )
        lines = [
            "",
            screen_line("analysis", "Fine-mapping completed", indent=2),
            screen_field(
                "analysis", "Engine", self.engine,
                indent=6, label_width=label_width,
            ),
            screen_field(
                status_kind, "Overall status", status.replace("_", " "),
                indent=6, label_width=label_width,
            ),
            screen_field(
                "count", "Primary loci attempted", result.get("n_attempted", 0),
                indent=6, label_width=label_width,
            ),
            screen_field(
                "success", "Primary loci successful", result.get("n_successful", 0),
                indent=6, label_width=label_width,
            ),
        ]
        failed = int(result.get("n_failed", 0) or 0)
        if failed:
            lines.append(
                screen_field(
                    "warning", "Failed or skipped loci", failed,
                    indent=6, label_width=label_width,
                )
            )
        warnings = int(result.get("n_warnings", 0) or 0)
        lines.append(
            screen_field(
                "warning" if warnings else "info",
                "Loci with scientific warnings",
                warnings,
                indent=6,
                label_width=label_width,
            )
        )
        for reason, count in sorted(
            (result.get("warning_reason_counts") or {}).items()
        ):
            lines.append(
                screen_field(
                    "warning",
                    "Warning reason",
                    _reason_count(reason, count),
                    indent=6,
                    label_width=label_width,
                )
            )
        for reason, count in sorted(
            (result.get("failure_reason_counts") or {}).items()
        ):
            lines.append(
                screen_field(
                    "error",
                    "Failure reason",
                    _reason_count(reason, count),
                    indent=6,
                    label_width=label_width,
                )
            )
        lines.extend(
            [
                screen_field(
                    "analysis", "Connected overlap groups",
                    result.get("n_overlap_groups", 0),
                    indent=6, label_width=label_width,
                ),
                screen_field(
                    "success", "Final credible sets",
                    result.get(
                        "n_final_credible_sets",
                        result.get("n_credible_sets", 0),
                    ),
                    indent=6, label_width=label_width,
                ),
            ]
        )
        if output_root is not None:
            lines.append(
                screen_field(
                    "info",
                    "Output directory",
                    output_root,
                    indent=6,
                    label_width=label_width,
                )
            )
        for kind, label, key in (
            ("success", "Final combined result", "final_combined_credible_sets"),
            ("success", "FLAMES input", "flames_input"),
            ("info", "Locus QC", "locus_status"),
            ("info", "QC summary", "overlap_resolution_summary"),
        ):
            value = result.get(key)
            if value:
                lines.append(
                    screen_field(
                        kind,
                        label,
                        _output_path(value, output_root),
                        indent=6, label_width=label_width,
                    )
                )
        lines.extend(
            [
                screen_field(
                    "info",
                    "Full detailed log",
                    _output_path(log_file, output_root),
                    indent=6, label_width=label_width,
                ),
                "",
            ]
        )
        self.console.print("\n".join(lines))


__all__ = ["FINE_MAPPING_PROGRESS_STAGES", "FineMappingScreen"]
