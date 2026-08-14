"""Shared runtime options for every PostGWAS command."""

from __future__ import annotations

import argparse

from rich_argparse import RichHelpFormatter

from postgwas.config import load_configuration
from postgwas.core.ui import help_with_default


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def _positive_float(value: str) -> float:
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than 0")
    return parsed


def get_compute_parser() -> argparse.ArgumentParser:
    """Return reusable runtime flags without assigning argparse defaults."""
    configuration = load_configuration()
    defaults = configuration.execution
    parser = argparse.ArgumentParser(
        add_help=False,
        formatter_class=RichHelpFormatter,
    )
    group = parser.add_argument_group("Performance")
    group.add_argument(
        "--threads",
        type=_positive_int,
        default=argparse.SUPPRESS,
        metavar="N",
        help=help_with_default(
            "Maximum number of tasks PostGWAS may run at the same time. "
            "If omitted, PostGWAS chooses automatically",
            defaults.threads,
        ),
    )
    group.add_argument(
        "--memory-gb",
        type=_positive_float,
        default=argparse.SUPPRESS,
        metavar="GB",
        help=help_with_default(
            "Maximum total memory PostGWAS may use, in gigabytes. "
            "If omitted, PostGWAS uses up to "
            f"{defaults.usable_memory_fraction * 100:g}%% of available physical memory",
            "%g GB" % defaults.memory_gb,
        ),
    )
    group.add_argument(
        "--seed",
        type=int,
        default=argparse.SUPPRESS,
        metavar="N",
        help=help_with_default(
            "Number used to make randomized analyses reproducible. "
            "Use the same value to reproduce a run",
            defaults.random_seed,
        ),
    )
    screen = parser.add_argument_group("Screen output")
    visibility = screen.add_mutually_exclusive_group()
    visibility.add_argument(
        "--show-screen",
        dest="show_screen",
        action="store_true",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Display PostGWAS progress, summaries, warnings, and errors in the terminal",
            configuration.logging.show_screen,
        ),
    )
    visibility.add_argument(
        "--hide-screen",
        dest="show_screen",
        action="store_false",
        default=argparse.SUPPRESS,
        help=(
            "Hide PostGWAS terminal output while continuing to save the same "
            "screen stream to logging.screen_log_file."
        ),
    )
    continuation = parser.add_argument_group("Run continuation")
    continuation.add_argument(
        "--resume",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Reuse existing outputs only after the selected module validates "
            "their completion state and provenance",
            configuration.run.resume,
        ),
    )
    continuation.add_argument(
        "--overwrite",
        action="store_true",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Replace existing outputs owned by the selected run; this takes "
            "precedence over resume",
            configuration.run.overwrite,
        ),
    )
    return parser


def resolve_compute_args(namespace: argparse.Namespace) -> argparse.Namespace:
    """Fill omitted compute values once from the central packaged configuration."""
    if all(hasattr(namespace, name) for name in ("threads", "memory_gb", "seed")):
        return namespace
    defaults = load_configuration().execution
    if not hasattr(namespace, "threads"):
        namespace.threads = defaults.threads
    if not hasattr(namespace, "memory_gb"):
        namespace.memory_gb = defaults.memory_gb
    if not hasattr(namespace, "seed"):
        namespace.seed = defaults.random_seed
    return namespace


def memory_limit(memory_gb: float) -> str:
    """Convert the canonical GB value for tools that require a unit suffix."""
    return f"{memory_gb:g}G"
