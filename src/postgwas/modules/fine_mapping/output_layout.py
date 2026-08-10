"""Resolve the canonical fine-mapping output layout below one run directory."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Mapping

from postgwas.core.paths import configured_output_path


def resolve_output_paths(
    output_directory: str | Path,
    layout: Mapping[str, object],
    dataset_id: str,
) -> dict[str, Path]:
    """Render every configured path once without permitting directory escape."""
    return {
        name: configured_output_path(
            output_directory,
            str(pattern),
            dataset_id=dataset_id,
        )
        for name, pattern in layout.items()
        if name != "cleanup_successful_workers"
    }


def cleanup_successful_workers(
    workers_directory: str | Path,
    *,
    enabled: bool,
) -> bool:
    """Remove validated worker copies when the configured retention policy allows."""
    directory = Path(workers_directory)
    if not enabled or not directory.exists():
        return False
    shutil.rmtree(directory)
    return True


__all__ = ["cleanup_successful_workers", "resolve_output_paths"]
