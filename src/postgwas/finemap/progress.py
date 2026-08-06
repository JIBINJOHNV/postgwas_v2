"""Structured progress and configuration logging shared by both engines."""

import csv
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

from postgwas.finemap.defaults import DEFAULT_PROGRESS_DECIMAL_PLACES


PROGRESS_COLUMNS = (
    "timestamp_utc",
    "scope",
    "stage",
    "status",
    "completed",
    "total",
    "percentage",
    "remaining",
    "message",
)


def progress_values(completed, total, decimal_places=DEFAULT_PROGRESS_DECIMAL_PLACES):
    """Return validated completed, total, percentage, and remaining values."""
    completed = int(completed)
    total = int(total)
    if total < 0 or completed < 0 or completed > total:
        raise ValueError(
            f"Invalid progress values: completed={completed}, total={total}"
        )
    percentage = 100.0 if total == 0 else round(100.0 * completed / total, decimal_places)
    return completed, total, percentage, total - completed


class ProgressRecorder:
    """Write every progress event to both a logger and a machine-readable TSV."""

    def __init__(
        self,
        progress_file,
        logger=None,
        overwrite=True,
        decimal_places=DEFAULT_PROGRESS_DECIMAL_PLACES,
    ):
        self.progress_file = Path(progress_file)
        self.progress_file.parent.mkdir(parents=True, exist_ok=True)
        self.logger = logger or logging.getLogger("postgwas.finemap")
        self.decimal_places = int(decimal_places)
        self.latest = {}
        if overwrite or not self.progress_file.exists():
            with self.progress_file.open("w", newline="", encoding="utf-8") as handle:
                csv.DictWriter(handle, fieldnames=PROGRESS_COLUMNS, delimiter="\t").writeheader()

    def record(self, scope, stage, completed, total, status, message=""):
        completed, total, percentage, remaining = progress_values(
            completed, total, decimal_places=self.decimal_places
        )
        row = {
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "scope": str(scope),
            "stage": str(stage),
            "status": str(status),
            "completed": completed,
            "total": total,
            "percentage": percentage,
            "remaining": remaining,
            "message": str(message),
        }
        with self.progress_file.open("a", newline="", encoding="utf-8") as handle:
            csv.DictWriter(handle, fieldnames=PROGRESS_COLUMNS, delimiter="\t").writerow(row)
        self.latest[row["scope"]] = row
        self.logger.info(
            "[PROGRESS] scope=%s stage=%s status=%s completed=%d/%d "
            "percentage=%s%% remaining=%d message=%s",
            row["scope"], row["stage"], row["status"], completed, total,
            f"{percentage:.{self.decimal_places}f}", remaining, row["message"],
        )
        return row


def write_run_configuration(output_file, configuration):
    """Write a JSON-safe, sorted record of inputs, defaults, and derived values."""
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    def serialise(value):
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, dict):
            return {str(key): serialise(item) for key, item in value.items()}
        if isinstance(value, (list, tuple, set)):
            return [serialise(item) for item in value]
        if value is None or isinstance(value, (str, int, float, bool)):
            return value
        return str(value)

    output_file.write_text(
        json.dumps(serialise(configuration), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return output_file
