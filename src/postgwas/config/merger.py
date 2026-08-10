"""Deterministic configuration merging without module-specific special cases."""

from __future__ import annotations

from copy import deepcopy
from typing import Any, Mapping

from postgwas.core.errors import ConfigurationError


def deep_merge(base: Mapping[str, Any], override: Mapping[str, Any]) -> dict[str, Any]:
    """Return a recursive merge where ``override`` always wins."""
    merged = deepcopy(dict(base))
    for key, value in override.items():
        if isinstance(value, Mapping) and isinstance(merged.get(key), Mapping):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = deepcopy(value)
    return merged


def apply_dotted_overrides(
    configuration: Mapping[str, Any], overrides: Mapping[str, Any] | None
) -> dict[str, Any]:
    """Apply explicit CLI overrides expressed as dotted configuration paths."""
    result = deepcopy(dict(configuration))
    for dotted_path, value in (overrides or {}).items():
        if value is None:
            continue
        parts = dotted_path.split(".")
        if not parts or any(not part for part in parts):
            raise ConfigurationError("Invalid CLI configuration path: %r" % dotted_path)
        cursor = result
        for part in parts[:-1]:
            child = cursor.get(part)
            if not isinstance(child, dict):
                child = {}
                cursor[part] = child
            cursor = child
        cursor[parts[-1]] = value
    return result
