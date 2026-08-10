"""Atomic writers for small structured PostGWAS reports."""

from __future__ import annotations

import csv
import os
from pathlib import Path
import tempfile
from typing import Any, Mapping, Sequence

import yaml


def _atomic_text_path(destination: Path):
    destination.parent.mkdir(parents=True, exist_ok=True)
    return tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        dir=destination.parent,
        prefix=".%s." % destination.name,
        delete=False,
    )


def write_yaml_report(value: Mapping[str, Any], path: str | Path) -> Path:
    """Write one ordered YAML document without exposing a partial result."""
    destination = Path(path)
    temporary_path = None
    try:
        with _atomic_text_path(destination) as handle:
            temporary_path = Path(handle.name)
            yaml.safe_dump(
                dict(value), handle, sort_keys=False, allow_unicode=True,
            )
        os.replace(temporary_path, destination)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    return destination


def write_delimited_report(
    records: Sequence[Mapping[str, Any]],
    path: str | Path,
    *,
    fieldnames: Sequence[str],
    delimiter: str,
    null_value: str,
) -> Path:
    """Write homogeneous records atomically using one configured delimiter."""
    destination = Path(path)
    temporary_path = None
    try:
        with _atomic_text_path(destination) as handle:
            temporary_path = Path(handle.name)
            writer = csv.DictWriter(
                handle,
                fieldnames=list(fieldnames),
                delimiter=delimiter,
                extrasaction="ignore",
            )
            writer.writeheader()
            for record in records:
                writer.writerow({
                    name: null_value if record.get(name) is None else record.get(name)
                    for name in fieldnames
                })
        os.replace(temporary_path, destination)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    return destination


__all__ = ["write_delimited_report", "write_yaml_report"]
