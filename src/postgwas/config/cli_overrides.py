"""Helpers keeping argparse defaults separate from configuration defaults."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Mapping


GLOBAL_RUNTIME_OVERRIDES = {
    "show_screen": "logging.show_screen",
    "resume": "run.resume",
    "overwrite": "run.overwrite",
}


def get_dotted(value: Any, dotted_path: str) -> Any:
    current = value
    for part in dotted_path.split("."):
        current = getattr(current, part) if hasattr(current, part) else current[part]
    return current


def explicit_overrides(
    namespace: argparse.Namespace, destination_to_path: Mapping[str, str]
) -> dict[str, Any]:
    """Return explicit CLI values in their YAML-compatible representation."""
    paths = dict(GLOBAL_RUNTIME_OVERRIDES)
    paths.update(destination_to_path)
    return {
        dotted_path: str(value) if isinstance(value, Path) else value
        for destination, dotted_path in paths.items()
        if hasattr(namespace, destination)
        for value in (getattr(namespace, destination),)
    }
