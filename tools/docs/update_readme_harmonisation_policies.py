#!/usr/bin/env python3
"""Synchronize the README harmonisation policy table with canonical YAML."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys
from typing import Any

import yaml


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
README = REPOSITORY_ROOT / "README.md"
POLICIES_YAML = (
    REPOSITORY_ROOT
    / "src"
    / "postgwas"
    / "config"
    / "defaults"
    / "modules"
    / "harmonisation.yaml"
)
START = "<!-- BEGIN GENERATED HARMONISATION POLICY REFERENCE -->"
END = "<!-- END GENERATED HARMONISATION POLICY REFERENCE -->"


def _scalar(value: Any) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, list):
        return "[" + ", ".join(_scalar(item) for item in value) + "]"
    if isinstance(value, dict):
        return "{" + ", ".join(
            f"{key}: {_scalar(item)}" for key, item in value.items()
        ) + "}"
    if isinstance(value, str) and value == "":
        return '""'
    return str(value)


def _code(value: Any) -> str:
    return f"`{_scalar(value).replace('|', '&#124;')}`"


def _members(values: list[Any]) -> str:
    return ", ".join(_code(value) for value in values)


def _allowed(validator: dict[str, Any]) -> str:
    kind = validator["kind"]
    allow_none = bool(validator.get("allow_none"))

    if kind == "flag":
        result = "`true` or `false`"
    elif kind == "enum":
        result = _members(validator["members"])
    elif kind == "num":
        minimum = validator.get("minimum")
        maximum = validator.get("maximum")
        label = "Integer" if validator.get("integer") else "Number"
        if minimum is not None and maximum is not None:
            operator = ">" if validator.get("exclusive_minimum") else "≥"
            result = f"{label} {operator} {_code(minimum)} and ≤ {_code(maximum)}"
        elif minimum is not None:
            operator = ">" if validator.get("exclusive_minimum") else "≥"
            result = f"{label} {operator} {_code(minimum)}"
        elif maximum is not None:
            result = f"{label} ≤ {_code(maximum)}"
        else:
            result = label
    elif kind == "list":
        emptiness = "possibly empty" if validator.get("allow_empty") else "non-empty"
        uniqueness = "unique " if validator.get("unique") else ""
        item = validator["item"]
        if item["kind"] == "enum":
            result = (
                f"{emptiness} {uniqueness}list selected from "
                f"{_members(item['members'])}"
            )
        else:
            result = f"{emptiness} {uniqueness}list of text values"
        if validator.get("nested_alternatives"):
            result += "; nested lists express alternatives"
    elif kind == "chromosome_set":
        result = "Unique chromosome list selected from " + _members(
            validator["members"]
        )
    elif kind == "string_map":
        result = "String map whose values are selected from " + _members(
            validator["value_members"]
        )
    elif kind == "text":
        description = validator.get("description", "text")
        result = f"Text satisfying {description}"
    elif kind == "numeric_pair":
        result = (
            "Two numbers (low < high), each from "
            f"{_code(validator['minimum'])} through {_code(validator['maximum'])}"
        )
    else:  # pragma: no cover - a new validator must be documented deliberately
        raise ValueError(f"Unsupported policy validator kind: {kind}")

    if allow_none:
        result += "; `null` is also allowed"
    return result


def _brief(help_text: str) -> str:
    text = " ".join(str(help_text).split())
    if text.startswith("What it does:"):
        body = text[len("What it does:") :].lstrip()
        action = re.split(
            r"\b(?:When it applies|Values|Example|Note):",
            body,
            maxsplit=1,
        )[0].strip()
        applies_match = re.search(
            r"\bWhen it applies:\s*(.*?)(?=\b(?:Values|Example|Note):|$)",
            body,
        )
        if applies_match:
            action += " Applies: " + applies_match.group(1).strip()
        text = action
    elif ". " in text:
        text = text.split(". ", 1)[0] + "."
    return text.replace("|", "&#124;")


def render_reference() -> str:
    document = yaml.safe_load(POLICIES_YAML.read_text(encoding="utf-8"))
    policies = document["policies"]
    sections: list[str] = [START]

    for group_spec in sorted(document["group_order"], key=lambda group: group["order"]):
        group = group_spec["group"]
        settings = policies[group]
        ordered_keys = group_spec["keys"]
        if set(ordered_keys) != set(settings):
            missing = sorted(set(settings) - set(ordered_keys))
            extra = sorted(set(ordered_keys) - set(settings))
            raise ValueError(
                f"group_order mismatch for {group}: missing={missing}, extra={extra}"
            )

        sections.extend(
            [
                "",
                "<details>",
                "<summary><code>%s</code> — %s (%d %s)</summary>"
                % (
                    group,
                    " ".join(str(group_spec["stage"]).split()),
                    len(ordered_keys),
                    "setting" if len(ordered_keys) == 1 else "settings",
                ),
                "",
                "| Setting | Default | Allowed values | Decision or action |",
                "|---|---|---|---|",
            ]
        )
        for key in ordered_keys:
            spec = settings[key]
            action = _brief(spec["help"])
            recommendation = spec.get("recommendation")
            if (
                recommendation is not None
                and _scalar(recommendation) != _scalar(spec["default"])
            ):
                action += " Registry recommendation: %s." % _code(recommendation)
            sections.append(
                "| `%s.%s` | %s | %s | %s |"
                % (
                    group,
                    key,
                    _code(spec["default"]),
                    _allowed(spec["validator"]),
                    action,
                )
            )
        sections.extend(["", "</details>"])

    sections.extend(["", END])
    return "\n".join(sections)


def synchronized_readme() -> str:
    current = README.read_text(encoding="utf-8")
    if current.count(START) != 1 or current.count(END) != 1:
        raise ValueError("README must contain exactly one policy-reference marker pair")
    before, remainder = current.split(START, 1)
    _old, after = remainder.split(END, 1)
    return before + render_reference() + after


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Synchronize README harmonisation policies with canonical YAML."
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Return non-zero when README is not synchronized; do not write it.",
    )
    args = parser.parse_args()

    current = README.read_text(encoding="utf-8")
    expected = synchronized_readme()
    if args.check:
        if current != expected:
            print("README harmonisation policy reference is out of date.", file=sys.stderr)
            return 1
        print("README harmonisation policy reference is synchronized.")
        return 0

    README.write_text(expected, encoding="utf-8")
    print("Updated README harmonisation policy reference.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
