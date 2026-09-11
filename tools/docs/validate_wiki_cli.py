#!/usr/bin/env python3
"""Validate commands in the README and published Wiki against the live CLI."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import shlex
import subprocess
from typing import Iterable

import yaml

from postgwas.core.ui import cli_option_values


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = REPOSITORY_ROOT / "docs" / "wiki.yml"
ROOT_README = REPOSITORY_ROOT / "README.md"
FENCE_PATTERN = re.compile(r"```(?:console|bash|shell)\n(.*?)```", re.DOTALL)
INLINE_COMMAND_PATTERN = re.compile(r"`(postgwas(?:\s+[^`\n]+)?)`")
LONG_OPTION_PATTERN = re.compile(r"(?<![A-Za-z0-9-])(--[A-Za-z][A-Za-z0-9-]*)")
TABLE_COMMAND_PATTERN = re.compile(r"^\s*│\s*([a-z][a-z0-9_-]*)\s*│", re.MULTILINE)
CONFIG_ACTION_PATTERN = re.compile(r"\{([a-z]+(?:,[a-z]+)+)\}")
# These are help-routing inputs, not analysis defaults. Keep the documented
# selections when asking the live CLI which workflow's options are available.
# GCTA's set source is selected by option presence; help does not read the file.
PIPELINE_CONTEXT_OPTIONS = (
    "--modules", "--clumping-methods", "--finemap-method", "--tools",
    "--method", "--analysis", "--cojo-mode", "--gmt", "--fastbat-set-list",
)
PIPELINE_CONTEXT_FLAGS = (
    "--apply-filter", "--apply-imputation", "--apply-manhattan", "--heritability",
)


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


def documented_commands(
    manifest: Path = DEFAULT_MANIFEST, *, include_inline: bool = True,
) -> tuple[tuple[Path, str], ...]:
    commands = []
    for source in (ROOT_README, *_published_sources(manifest)):
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
                        commands.append((source, shlex.join(candidate)))
        if include_inline:
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
    match = re.match(r"--[A-Za-z][A-Za-z0-9-]*", token)
    return match.group(0) if match else None


def _starts_synopsis(tokens: list[str], index: int) -> bool:
    """Recognise exact shorthand, including the two-token [pipeline options]."""
    return any(
        re.fullmatch(
            r"\[(?:(?:pipeline )?options|overrides|arguments|args)\]",
            " ".join(tokens[index:index + width]),
            re.I,
        )
        for width in (1, 2)
    )


def _pipeline_help_arguments(tokens: list[str]) -> list[str]:
    arguments = ["pipeline"]
    for option in PIPELINE_CONTEXT_OPTIONS:
        values = []
        documented_values = cli_option_values(tokens, option)
        for index, value in enumerate(documented_values):
            if _starts_synopsis(documented_values, index):
                break  # Prose such as [pipeline options], not a CLI value.
            values.append(value)
        if values:
            arguments.extend((option, *values))
        elif option != "--modules" and any(
            token == option or token.startswith(option + "=") for token in tokens
        ):
            raise CliDocumentationError(
                f"Documented pipeline selector {option} has no value"
            )
    arguments.extend(flag for flag in PIPELINE_CONTEXT_FLAGS if flag in tokens)
    # Example configuration paths need not exist. Do not load --run-config or
    # pass analysis inputs: the CLI is inspected, never an analysis executed.
    return [*arguments, "--help"]


def _option_help_blocks(help_text: str) -> Iterable[tuple[tuple[str, ...], str, str]]:
    """Yield long-option aliases, invocation and wrapped help description."""
    current_options: tuple[str, ...] = ()
    block_lines: list[str] = []
    lines = help_text.splitlines()
    option_indent = min(
        (len(line) - len(line.lstrip()) for line in lines
         if line.lstrip().startswith("--")),
        default=0,
    )

    def finish_block() -> tuple[tuple[str, ...], str, str]:
        parts = re.split(r"\s{2,}", block_lines[0].lstrip(), maxsplit=1)
        description = " ".join([*parts[1:], *(line.strip() for line in block_lines[1:])]).strip()
        return current_options, parts[0], description

    for line in lines:
        invocation = re.split(r"\s{2,}", line.lstrip(), maxsplit=1)[0]
        options = tuple(LONG_OPTION_PATTERN.findall(invocation))
        if (
            options and line.lstrip().startswith("--")
            and len(line) - len(line.lstrip()) == option_indent
        ):
            if current_options:
                yield finish_block()
            current_options = options
            block_lines = [line]
        elif current_options:
            block_lines.append(line)
    if current_options:
        yield finish_block()


def _choice_options(help_text: str) -> dict[str, set[str]]:
    """Return argparse choices displayed in each long-option help block."""
    choices: dict[str, set[str]] = {}
    for options, invocation, description in _option_help_blocks(help_text):
        block = f"{invocation} {description}"
        match = re.search(r"\{([^{}]+)\}", invocation)
        if not match:
            match = re.search(
                r"(?:choices?|accepted values?)\s*:?\s*\{([^{}]+)\}",
                block,
                re.IGNORECASE,
            )
        if not match:
            match = re.search(
                r"Available options:\s*(.+?)(?:\.(?:\s|$)|$)",
                block,
            )
        if not match:
            continue
        values = {value.strip() for value in match.group(1).split(",") if value.strip()}
        for option in options:
            choices[option] = values
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
        following = tokens[index + 1 :]
        if "=" in token:
            following = [token.partition("=")[2], *following]
        for value in following:
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


def _usage_required_options(help_text: str) -> set[str]:
    """Intersect unconditional usage options with live Required-labelled help.

    Manual usage can show an optional, defaulted selector outside brackets.
    Check conservative public metadata, not arbitrary conditional help prose.
    Conditional resource requirements are additionally checked for pipelines
    through their registered, method-aware dependency plan.
    """
    match = re.search(r"^usage:\s*(.*?)(?:\n\s*\n|\Z)", help_text, re.I | re.M | re.S)
    if match is None:
        return set()
    mandatory = []
    depth = 0
    for character in match.group(1):
        if character == "[":
            depth += 1
        elif character == "]":
            depth = max(0, depth - 1)
        elif depth == 0:
            mandatory.append(character)
    # Required alternative groups cannot be represented as an all-required set.
    # Leave their value-dependent validation to the module's own preflight.
    text = re.sub(r"\([^()]*\|[^()]*\)", "", "".join(mandatory))
    labelled_required = {
        option
        for options, _, description in _option_help_blocks(help_text)
        if description.startswith("Required:")
        for option in options
    }
    return set(LONG_OPTION_PATTERN.findall(text)) & labelled_required


def _pipeline_required_options(tokens: list[str]) -> set[str]:
    """Reuse the planner and registry without opening study/reference files."""
    from postgwas.config import load_configuration
    from postgwas.config.models.modules.single_cell import single_cell_pipeline_dependencies
    from postgwas.pipeline.planner import build_pipeline_plan, resolve_pipeline_dependency_overrides
    from postgwas.pipeline.registry import PIPELINE_REQUIRED_OPTIONS, REGISTRY

    configuration = load_configuration()
    targets = cli_option_values(tokens, "--modules")
    overrides = resolve_pipeline_dependency_overrides(tokens, configuration)
    selected_tools = cli_option_values(tokens, "--tools")
    if selected_tools:
        overrides["single_cell"] = single_cell_pipeline_dependencies(selected_tools)
    plan = build_pipeline_plan(
        targets,
        dependency_overrides=overrides,
        **{flag[2:].replace("-", "_"): flag in tokens for flag in PIPELINE_CONTEXT_FLAGS},
    )
    required = {item.flag for item in PIPELINE_REQUIRED_OPTIONS}
    for name in plan.steps:
        if name == "formatter" and name not in targets:
            continue  # The analysis formatter's input selections are generated.
        required.update(item.flag for item in REGISTRY.get(name).required_options)
    return required


def _missing_recipe_options(tokens: list[str], help_text: str) -> set[str]:
    """Check CLI-only recipes; a placeholder YAML cannot be resolved safely.

    This is static recipe validation, not an analysis-readiness test. It does
    not certify file contents, alternative inputs or every runtime condition.
    """
    if any(token in tokens for token in ("--help", "-h")):
        return set()
    if cli_option_values(tokens, "--run-config"):
        return set()
    if any(_starts_synopsis(tokens, index) for index in range(len(tokens))):
        raise CliDocumentationError("Executable examples must not abbreviate arguments with [options]")
    required = _usage_required_options(help_text)
    if tokens[1] == "pipeline":
        required.update(_pipeline_required_options(tokens))
    # Require a supplied value, not just a bare option name. All required
    # analysis-input flags represented by these contracts take values.
    aliases = {
        option: options
        for options, _, _ in _option_help_blocks(help_text)
        for option in options
    }
    return {
        option for option in required
        if not any(cli_option_values(tokens, alias) for alias in aliases.get(option, (option,)))
    }


def validate_documented_commands(
    *,
    executable: str = "postgwas",
    manifest: Path = DEFAULT_MANIFEST,
    commands_only: bool = False,
) -> int:
    public_commands, config_actions = _command_inventory(executable)
    errors = []
    checked = 0
    help_cache: dict[tuple[str, ...], tuple[set[str], dict[str, set[str]], str]] = {}
    recipes = set(documented_commands(manifest, include_inline=False))

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
                help_text,
            )
        available_options, option_choices, help_text = help_cache[cache_key]
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
            continue
        if (source, command) in recipes and tokens[1] not in {"config", "resources"}:
            try:
                missing = _missing_recipe_options(tokens, help_text)
            except (ValueError, RuntimeError) as exc:
                errors.append(f"{source.relative_to(REPOSITORY_ROOT)}: {exc}: {command}")
                continue
            if missing:
                errors.append(
                    f"{source.relative_to(REPOSITORY_ROOT)}: incomplete CLI recipe; "
                    f"missing required value(s) {', '.join(sorted(missing))}: {command}"
                )

    if errors:
        raise CliDocumentationError("\n".join(errors))
    return checked


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate commands in README and published Wiki blocks against PostGWAS help."
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
    scope = "command names" if args.commands_only else "commands, options and CLI recipe requirements"
    print(f"Wiki CLI validation passed: {checked} documented commands ({scope})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
