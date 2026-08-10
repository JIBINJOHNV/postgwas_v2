#!/usr/bin/env python3
"""Validate commands in published Wiki console blocks against the live CLI."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import shlex
import subprocess
from typing import Iterable

import yaml


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = REPOSITORY_ROOT / "docs" / "wiki.yml"
FENCE_PATTERN = re.compile(r"```(?:console|bash|shell)\n(.*?)```", re.DOTALL)
INLINE_COMMAND_PATTERN = re.compile(r"`(postgwas(?:\s+[^`\n]+)?)`")
LONG_OPTION_PATTERN = re.compile(r"(?<![A-Za-z0-9-])(--[A-Za-z][A-Za-z0-9-]*)")
TABLE_COMMAND_PATTERN = re.compile(r"^\s*│\s*([a-z][a-z0-9_-]*)\s*│", re.MULTILINE)
CONFIG_ACTION_PATTERN = re.compile(r"\{([a-z]+(?:,[a-z]+)+)\}")


class CliDocumentationError(ValueError):
    """Raised when a documented command is absent from the installed CLI."""


def _published_sources(manifest: Path) -> tuple[Path, ...]:
    document = yaml.safe_load(manifest.read_text(encoding="utf-8")) or {}
    pages = document.get("wiki", {}).get("pages", [])
    sources = []
    for page in pages:
        source = page.get("source") if isinstance(page, dict) else None
        if not isinstance(source, str):
            raise CliDocumentationError("Every Wiki page must define a source")
        path = (REPOSITORY_ROOT / source).resolve()
        if not path.is_file():
            raise CliDocumentationError(f"Published Wiki source does not exist: {source}")
        sources.append(path)
    return tuple(sources)


def _logical_commands(block: str) -> Iterable[str]:
    buffered = ""
    for raw_line in block.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        continued = line.endswith("\\")
        if continued:
            line = line[:-1].rstrip()
        buffered = f"{buffered} {line}".strip()
        if not continued:
            yield buffered
            buffered = ""
    if buffered:
        yield buffered


def documented_commands(manifest: Path = DEFAULT_MANIFEST) -> tuple[tuple[Path, str], ...]:
    commands = []
    for source in _published_sources(manifest):
        text = source.read_text(encoding="utf-8")
        for block in FENCE_PATTERN.findall(text):
            for command in _logical_commands(block):
                try:
                    tokens = shlex.split(command)
                except ValueError as exc:
                    raise CliDocumentationError(
                        f"Cannot parse command in {source.relative_to(REPOSITORY_ROOT)}: "
                        f"{command}: {exc}"
                    ) from exc
                indexes = [index for index, token in enumerate(tokens) if token == "postgwas"]
                for index in indexes:
                    candidate = tokens[index:]
                    if len(candidate) > 1 and candidate[1] != "activate":
                        commands.append((source, " ".join(candidate)))
        for command in INLINE_COMMAND_PATTERN.findall(text):
            commands.append((source, " ".join(command.split())))
    return tuple(commands)


def _run_help(executable: str, arguments: list[str]) -> str:
    completed = subprocess.run(
        [executable, *arguments],
        cwd=REPOSITORY_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    output = completed.stdout + completed.stderr
    if completed.returncode != 0:
        raise CliDocumentationError(
            f"CLI help failed ({completed.returncode}): "
            f"{executable} {' '.join(arguments)}\n{output.strip()}"
        )
    return output


def _command_inventory(executable: str) -> tuple[set[str], set[str]]:
    global_help = _run_help(executable, ["--help"])
    commands = set(TABLE_COMMAND_PATTERN.findall(global_help))
    if not commands:
        raise CliDocumentationError("Could not read command names from postgwas --help")

    config_help = _run_help(executable, ["config", "--help"])
    match = CONFIG_ACTION_PATTERN.search(config_help)
    if not match:
        raise CliDocumentationError("Could not read actions from postgwas config --help")
    return commands, set(match.group(1).split(","))


def _clean_option(token: str) -> str | None:
    match = re.search(r"--[A-Za-z][A-Za-z0-9-]*", token)
    return match.group(0) if match else None


def _pipeline_help_arguments(tokens: list[str]) -> list[str]:
    if "--modules" not in tokens:
        return ["pipeline", "--help"]
    start = tokens.index("--modules")
    modules = []
    for token in tokens[start + 1 :]:
        if token.startswith("-") or token.startswith("["):
            break
        modules.append(token.strip("[],{}"))
    return ["pipeline", "--modules", *modules, "--help"] if modules else ["pipeline", "--help"]


def _choice_options(help_text: str) -> dict[str, set[str]]:
    """Return argparse choices displayed in each long-option help block."""
    choices: dict[str, set[str]] = {}
    current_options: tuple[str, ...] = ()
    block_lines: list[str] = []

    def finish_block() -> None:
        if not current_options:
            return
        block = " ".join(block_lines)
        match = re.search(r"\{([^{}]+)\}", block_lines[0])
        if not match:
            match = re.search(
                r"(?:choices?|accepted values?)\s*:?\s*\{([^{}]+)\}",
                block,
                re.IGNORECASE,
            )
        if not match:
            return
        values = {value.strip() for value in match.group(1).split(",") if value.strip()}
        for option in current_options:
            choices[option] = values

    for line in help_text.splitlines():
        options = tuple(LONG_OPTION_PATTERN.findall(line))
        if options and line.lstrip().startswith("--"):
            finish_block()
            current_options = options
            block_lines = [line]
        elif current_options:
            block_lines.append(line)
    finish_block()
    return choices


def _documented_choice_errors(
    tokens: list[str], choices: dict[str, set[str]]
) -> list[str]:
    errors = []
    for index, token in enumerate(tokens):
        option = _clean_option(token)
        if option not in choices:
            continue
        values = []
        for value in tokens[index + 1 :]:
            if value.startswith("-") or any(marker in value for marker in "[]{}|"):
                break
            if value.isupper():
                break
            values.append(value)
        invalid = [value for value in values if value not in choices[option]]
        if invalid:
            errors.append(
                f"{option} value(s) {', '.join(invalid)}; accepted: "
                f"{', '.join(sorted(choices[option]))}"
            )
    return errors


def validate_documented_commands(
    *,
    executable: str = "postgwas",
    manifest: Path = DEFAULT_MANIFEST,
    commands_only: bool = False,
) -> int:
    public_commands, config_actions = _command_inventory(executable)
    errors = []
    checked = 0
    help_cache: dict[tuple[str, ...], tuple[set[str], dict[str, set[str]]]] = {}

    for source, command in documented_commands(manifest):
        tokens = shlex.split(command)
        if len(tokens) < 2 or tokens[1] == "--help":
            checked += 1
            continue
        if tokens[1] == "COMMAND":
            continue

        if tokens[1] == "--validate":
            help_arguments = ["--validate", "--help"]
        elif tokens[1] == "config":
            if len(tokens) < 3 or tokens[2].startswith("-"):
                help_arguments = ["config", "--help"]
            elif tokens[2] not in config_actions:
                errors.append(
                    f"{source.relative_to(REPOSITORY_ROOT)}: unsupported config action "
                    f"{tokens[2]!r}: {command}"
                )
                continue
            else:
                help_arguments = ["config", tokens[2]]
                if tokens[2] == "export":
                    for selector in ("--module", "--pipeline"):
                        if selector not in tokens:
                            continue
                        selector_index = tokens.index(selector)
                        help_arguments.append(selector)
                        for value in tokens[selector_index + 1 :]:
                            if value.startswith("-"):
                                break
                            help_arguments.append(value)
                        break
                help_arguments.append("--help")
        elif tokens[1] not in public_commands:
            errors.append(
                f"{source.relative_to(REPOSITORY_ROOT)}: unsupported command "
                f"{tokens[1]!r}: {command}"
            )
            continue
        elif tokens[1] == "pipeline":
            help_arguments = _pipeline_help_arguments(tokens)
        else:
            help_arguments = [tokens[1], "--help"]

        checked += 1
        if commands_only:
            continue
        cache_key = tuple(help_arguments)
        if cache_key not in help_cache:
            help_text = _run_help(executable, help_arguments)
            help_cache[cache_key] = (
                set(LONG_OPTION_PATTERN.findall(help_text)),
                _choice_options(help_text),
            )
        available_options, option_choices = help_cache[cache_key]
        documented_options = {
            option
            for token in tokens[2:]
            if (option := _clean_option(token)) is not None and option != "--help"
        }
        unknown = sorted(documented_options - available_options)
        if unknown:
            errors.append(
                f"{source.relative_to(REPOSITORY_ROOT)}: unsupported option(s) "
                f"{', '.join(unknown)}: {command}"
            )
            continue

        choice_errors = _documented_choice_errors(tokens, option_choices)
        if choice_errors:
            errors.append(
                f"{source.relative_to(REPOSITORY_ROOT)}: unsupported choice in "
                f"{command}: {'; '.join(choice_errors)}"
            )

    if errors:
        raise CliDocumentationError("\n".join(errors))
    return checked


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate commands in published Wiki console blocks against PostGWAS help."
    )
    parser.add_argument("--executable", default="postgwas")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument(
        "--commands-only",
        action="store_true",
        help="Check command names and config actions without importing module CLIs.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        checked = validate_documented_commands(
            executable=args.executable,
            manifest=args.manifest,
            commands_only=args.commands_only,
        )
    except CliDocumentationError as exc:
        raise SystemExit(f"Wiki CLI validation failed:\n{exc}") from exc
    scope = "command names" if args.commands_only else "commands and options"
    print(f"Wiki CLI validation passed: {checked} documented commands ({scope})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
