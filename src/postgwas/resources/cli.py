"""Installed command-line interface for reproducible scientific resources."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import yaml

from postgwas.core.resource_preparation import ResourcePreparationError
from postgwas.core.ui import (
    AlignedRichHelpFormatter,
    format_cli_examples,
    screen_field,
    screen_line,
)
from postgwas.resources.magma_functional_mapping import packaged_configuration, run


RESOURCE_NAMES = ("magma",)
RESOURCE_CONFIG_PATH, RESOURCE_CONFIGURATION = packaged_configuration()

RESOURCE_EXAMPLES = format_cli_examples(
    (
        "Install the pinned MAGMA functional-mapping bundle:",
        "postgwas resources prepare magma",
        (),
    ),
    (
        "Install it in another user-selected directory:",
        "postgwas resources prepare magma",
        (
            "--output-directory /path/to/postgwas/magma/functional_mapping",
            "--cache-directory /path/to/download-cache",
        ),
    ),
    (
        "Revalidate an installation and rebuild only its metadata:",
        "postgwas resources refresh magma",
        (),
    ),
)


def _add_resource_arguments(
    parser: argparse.ArgumentParser,
    *,
    include_cache: bool,
) -> None:
    parser.add_argument(
        "resource",
        choices=RESOURCE_NAMES,
        metavar="RESOURCE",
        help="Resource bundle to process. Currently available: magma.",
    )
    paths = parser.add_argument_group("Resource locations")
    paths.add_argument(
        "--output-directory",
        type=Path,
        metavar="PATH",
        help=(
            "Installation directory. If omitted, the packaged configuration "
            "uses %s." % RESOURCE_CONFIGURATION["default_output_directory"]
        ),
    )
    if include_cache:
        paths.add_argument(
            "--cache-directory",
            type=Path,
            metavar="PATH",
            help=(
                "Download cache. If omitted, the packaged configuration uses "
                "%s." % RESOURCE_CONFIGURATION["default_cache_directory"]
            ),
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="postgwas resources",
        usage="postgwas resources {prepare,refresh} RESOURCE [options]",
        description=(
            "Install or revalidate pinned scientific reference bundles without "
            "using a repository-specific script path."
        ),
        formatter_class=AlignedRichHelpFormatter,
        epilog=RESOURCE_EXAMPLES,
    )
    actions = parser.add_subparsers(dest="action", required=True)
    prepare_parser = actions.add_parser(
        "prepare",
        prog="postgwas resources prepare",
        usage="postgwas resources prepare RESOURCE [options]",
        help="Download, verify, and atomically install a resource bundle.",
        description=(
            "Download, checksum, validate, and atomically install a pinned "
            "scientific resource bundle."
        ),
        formatter_class=AlignedRichHelpFormatter,
        epilog=RESOURCE_EXAMPLES,
    )
    _add_resource_arguments(prepare_parser, include_cache=True)
    refresh_parser = actions.add_parser(
        "refresh",
        prog="postgwas resources refresh",
        usage="postgwas resources refresh RESOURCE [options]",
        help="Revalidate resources and regenerate configuration metadata.",
        description=(
            "Verify every installed scientific file against its recorded "
            "checksum, then regenerate only the manifest and PostGWAS YAML."
        ),
        formatter_class=AlignedRichHelpFormatter,
        epilog=RESOURCE_EXAMPLES,
    )
    _add_resource_arguments(refresh_parser, include_cache=False)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    print(screen_line("analysis", "MAGMA functional-mapping resources", indent=2))
    print(screen_field("info", "Operation", args.action, indent=4))
    try:
        installed, generated_config = run(
            args.action,
            config_path=RESOURCE_CONFIG_PATH,
            configuration=RESOURCE_CONFIGURATION,
            output_directory=args.output_directory,
            cache_directory=getattr(args, "cache_directory", None),
        )
    except (ResourcePreparationError, OSError, yaml.YAMLError) as exc:
        print(screen_field("error", "Failed", exc, indent=4), file=sys.stderr)
        return 2
    except Exception as exc:
        # The preparation layer finalises its log; the command boundary must not
        # replace that actionable failure with an unformatted Python traceback.
        print(
            screen_field(
                "error",
                "Unexpected failure",
                "%s: %s" % (type(exc).__name__, exc),
                indent=4,
            ),
            file=sys.stderr,
        )
        return 1
    print(screen_field("success", "Resource directory", installed, indent=4))
    print(screen_field("success", "PostGWAS config", generated_config, indent=4))
    return 0


__all__ = ["build_parser", "main"]
