"""Inspect and validate configuration without starting scientific work."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

from postgwas.config.loader import (
    MODULE_NAMES,
    load_configuration,
    load_module_configuration,
    write_resolved_configuration,
)
from postgwas.core.errors import ConfigurationError
from postgwas.config.exporter import (
    EXPORT_STYLES,
    render_module_configuration,
    render_pipeline_configuration,
)
from postgwas.pipeline.registry import REGISTRY
from postgwas.core.ui import AlignedRichHelpFormatter, format_cli_examples


PIPELINE_TARGETS = tuple(
    name for name in REGISTRY.names() if REGISTRY.get(name).pipeline_enabled
)


def _choice_lines(values: tuple[str, ...], size: int = 4) -> str:
    """Lay long choice sets out beneath their option instead of across the screen."""
    rows = [", ".join(values[index : index + size]) for index in range(0, len(values), size)]
    return "\n".join("  " + row for row in rows)


EXPORT_EXAMPLES = format_cli_examples(
    (
        "Export one fully annotated module configuration:",
        "postgwas config export",
        ("--module harmonisation", "--style full", "--output harmonisation.yaml"),
    ),
    (
        "Export a fine-mapping pipeline with all prerequisites:",
        "postgwas config export",
        ("--pipeline finemap", "--style minimal", "--output finemap_pipeline.yaml"),
    ),
    (
        "Export several final analyses with shared prerequisites included once:",
        "postgwas config export",
        ("--pipeline finemap magma", "--style values", "--output analysis_pipeline.yaml"),
    ),
)

CONFIG_EXAMPLES = format_cli_examples(
    (
        "Validate a complete run configuration:",
        "postgwas config validate",
        ("--config run.yaml",),
    ),
    (
        "Show the resolved settings for one module:",
        "postgwas config show",
        ("--config run.yaml", "--module mixer"),
    ),
    (
        "Export a reusable module configuration:",
        "postgwas config export",
        ("--module mixer", "--style full", "--output mixer.yaml"),
    ),
)


def _add_export_arguments(parser: argparse.ArgumentParser) -> None:
    target_group = parser.add_argument_group("Choose what to export")
    target = target_group.add_mutually_exclusive_group(required=True)
    target.add_argument(
        "--module",
        choices=MODULE_NAMES,
        metavar="MODULE",
        help=(
            "Export one standalone module.\n\n"
            "Available modules:\n" + _choice_lines(MODULE_NAMES)
        ),
    )
    target.add_argument(
        "--pipeline",
        nargs="+",
        choices=PIPELINE_TARGETS,
        metavar="MODULE",
        help=(
            "Export one or more final analyses. PostGWAS automatically includes\n"
            "every required preceding module in execution order.\n\n"
            "Available pipeline targets:\n" + _choice_lines(PIPELINE_TARGETS)
        ),
    )

    format_group = parser.add_argument_group("Choose the amount of explanation")
    format_group.add_argument(
        "--style",
        choices=EXPORT_STYLES,
        default="minimal",
        metavar="STYLE",
        help=(
            "full     Detailed scientific explanation for every setting.\n"
            "minimal  Short comments beside settings (default).\n"
            "values   Key-value pairs only; analysis values remain identical."
        ),
    )

    files = parser.add_argument_group("Input and destination")
    files.add_argument(
        "--run-config",
        dest="config",
        metavar="PATH",
        help="Optional existing run configuration whose values should be resolved.",
    )
    files.add_argument(
        "--output",
        metavar="PATH",
        help="Write directly to this file instead of standard output.",
    )


def build_export_parser(prog: str = "postgwas --config") -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        usage="%s (--module MODULE | --pipeline MODULE [MODULE ...]) [options]" % prog,
        description=(
            "Create a reusable YAML configuration from PostGWAS defaults.\n"
            "Long lists, including chromosome sets, remain on one line."
        ),
        formatter_class=AlignedRichHelpFormatter,
        epilog=EXPORT_EXAMPLES,
    )
    _add_export_arguments(parser)
    return parser


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="postgwas config",
        usage="postgwas config {validate,show,export} [options]",
        description="Validate, inspect, or export PostGWAS configuration.",
        formatter_class=AlignedRichHelpFormatter,
        epilog=CONFIG_EXAMPLES,
    )
    subparsers = parser.add_subparsers(dest="action", required=True)

    for action in ("validate", "show"):
        command = subparsers.add_parser(
            action,
            prog="postgwas config %s" % action,
            formatter_class=AlignedRichHelpFormatter,
            description=(
                "Validate a complete or module-specific YAML configuration."
                if action == "validate"
                else "Show the configuration after defaults and supplied values are resolved."
            ),
            epilog=CONFIG_EXAMPLES,
        )
        command.add_argument(
            "--config", metavar="PATH",
            help="YAML configuration to validate or resolve. Omit to use packaged defaults.",
        )
        command.add_argument(
            "--module", choices=MODULE_NAMES, metavar="MODULE",
            help="Restrict validation or display to one module.",
        )
        if action == "show":
            command.add_argument(
                "--output", metavar="PATH",
                help="Write the complete resolved configuration to this file.",
            )

    export = subparsers.add_parser(
        "export",
        prog="postgwas config export",
        help="Write a reusable module or pipeline configuration as YAML.",
        description=(
            "Create a reusable YAML configuration from PostGWAS defaults.\n"
            "Long lists, including chromosome sets, remain on one line."
        ),
        formatter_class=AlignedRichHelpFormatter,
        epilog=EXPORT_EXAMPLES,
    )
    _add_export_arguments(export)
    return parser


def main() -> int:
    shorthand = sys.argv[0].endswith(" --config")
    parser = build_export_parser() if shorthand else build_parser()
    args = parser.parse_args()
    if shorthand:
        args.action = "export"
    try:
        if args.action == "export":
            if args.module:
                rendered = render_module_configuration(
                    args.module,
                    config_file=args.config,
                    style=args.style,
                )
            else:
                rendered = render_pipeline_configuration(
                    args.pipeline,
                    config_file=args.config,
                    style=args.style,
                )
            if args.output:
                destination = Path(args.output).expanduser()
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_text(rendered, encoding="utf-8")
            else:
                print(rendered, end="")
            return 0

        if args.module:
            config = load_module_configuration(
                args.module, args.config
            )
        else:
            config = load_configuration(args.config)
    except ConfigurationError as exc:
        print(str(exc), file=sys.stderr)
        return 2

    if args.action == "validate":
        print("Configuration is valid.")
        return 0

    if args.output:
        if args.module:
            print(
                "--output writes complete resolved configurations; omit --module",
                file=sys.stderr,
            )
            return 2
        destination = write_resolved_configuration(config, args.output)
        print(destination)
    else:
        print(yaml.safe_dump(config.model_dump(mode="json"), sort_keys=False), end="")
    return 0
