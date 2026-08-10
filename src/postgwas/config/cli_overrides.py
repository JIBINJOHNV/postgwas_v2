"""Helpers keeping argparse defaults separate from configuration defaults."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Mapping


def get_dotted(value: Any, dotted_path: str) -> Any:
    current = value
    for part in dotted_path.split("."):
        current = getattr(current, part) if hasattr(current, part) else current[part]
    return current


def help_with_default(description: str, defaults: Any, dotted_path: str) -> str:
    value = get_dotted(defaults, dotted_path)
    return "%s Packaged default: %s." % (description.rstrip(". "), value)


def explicit_overrides(
    namespace: argparse.Namespace, destination_to_path: Mapping[str, str]
) -> dict[str, Any]:
    """Return explicit CLI values in their YAML-compatible representation."""
    return {
        dotted_path: str(value) if isinstance(value, Path) else value
        for destination, dotted_path in destination_to_path.items()
        if hasattr(namespace, destination)
        for value in (getattr(namespace, destination),)
    }
