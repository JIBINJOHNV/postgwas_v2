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


# Public destinations shared by the report-only QC assessment and physical
# filtering commands.  Each service maps these normalized policy fields onto
# its own canonical YAML location instead of maintaining parallel option maps.
VARIANT_QC_POLICY_CLI_FIELDS = {
    "minimum_neglog10_p": "minimum_neglog10_p",
    "minimum_maf": "maf_min",
    "missing_pvalue_action": "missing_pvalue_action",
    "missing_af_action": "missing_af_action",
    "minimum_info": "info_min",
    "maximum_info": "info_max",
    "missing_info_action": "missing_info_action",
    "maximum_af_difference": "maximum_af_difference",
    "include_indels": "include_indels",
    "remove_palindromic": "remove_palindromic",
    "palindromic_af_lower": "palindromic_lower",
    "palindromic_af_upper": "palindromic_upper",
    "remove_mhc": "remove_mhc",
}


def variant_qc_policy_override_paths(
    *,
    prefix: str = "",
    field_overrides: Mapping[str, str] | None = None,
) -> dict[str, str]:
    """Map shared QC CLI destinations onto one module's canonical fields."""
    fields = dict(VARIANT_QC_POLICY_CLI_FIELDS)
    fields.update(field_overrides or {})
    return {
        destination: "%s%s" % (prefix, field)
        for destination, field in fields.items()
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
