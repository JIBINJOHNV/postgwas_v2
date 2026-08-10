#!/usr/bin/env python3
"""Lazy top-level command dispatcher for PostGWAS."""

import sys

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from postgwas._mp_fix import *  # noqa: F401,F403 - multiprocessing bootstrap
from postgwas.core.ui import format_cli_examples
from postgwas.pipeline.registry import REGISTRY, resolve_reference


console = Console()

# Keep the historical misspelling as a compatibility alias without advertising it.
COMMAND_ALIASES = {"pathway_enrichmnet": "pathway_enrichment"}


def print_global_help(prog: str = "postgwas") -> None:
    header = Text("PostGWAS — Unified Toolkit for Post-GWAS Analyses", style="bold cyan")
    console.print(Panel(header, expand=False))

    table = Table(title="Available Modules", title_style="bold magenta", padding=(0, 1))
    table.add_column("Command", style="cyan", no_wrap=True)
    table.add_column("Description", style="green")
    for command, spec in sorted(REGISTRY.commands().items()):
        table.add_row(command, spec.description)
    console.print(table)
    console.print(
        "\n" + format_cli_examples(
            (
                "Prepare raw GWAS summary statistics:",
                "postgwas harmonisation",
                ("--help",),
            ),
            (
                "Validate a harmonised GWAS-VCF:",
                "postgwas --validate",
                ("--help",),
            ),
            (
                "Choose and inspect a multi-module workflow:",
                "postgwas pipeline",
                ("--help",),
            ),
            (
                "Export a complete fine-mapping pipeline configuration:",
                "postgwas config export",
                (
                    "--pipeline finemap",
                    "--style full",
                    "--output finemap_pipeline.yaml",
                ),
            ),
            title="Examples — start here",
        ) + "\n"
    )


def main():
    argv = sys.argv
    prog = argv[0]
    if len(argv) == 1 or argv[1] in ("-h", "--help"):
        print_global_help(prog)
        return 0

    if argv[1] == "--validate":
        from postgwas.modules.harmonisation.concordance.cli import main as validation_main

        sys.argv = [f"{prog} --validate"] + argv[2:]
        return validation_main()

    shorthand_export = argv[1] == "--config"
    command = "config" if shorthand_export else COMMAND_ALIASES.get(argv[1], argv[1])
    commands = REGISTRY.commands()
    if command not in commands:
        console.print(f"[red]{prog}: error: unknown module '{argv[1]}'[/red]\n")
        print_global_help(prog)
        return 2

    spec = commands[command]
    arguments = argv[2:]
    command_prog = f"{prog} --config" if shorthand_export else f"{prog}-{command}"
    sys.argv = [command_prog] + arguments
    try:
        entrypoint = resolve_reference(spec.cli_entrypoint)
    except (ImportError, AttributeError) as exc:
        console.print(
            f"[red]{prog}: module '{command}' could not be loaded: "
            f"{type(exc).__name__}: {exc}[/red]"
        )
        return 1
    result = entrypoint()
    return 0 if result is None else result


if __name__ == "__main__":
    raise SystemExit(main())
