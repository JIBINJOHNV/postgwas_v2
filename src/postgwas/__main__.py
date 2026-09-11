#!/usr/bin/env python3
"""Lazy top-level command dispatcher for PostGWAS."""

import sys
from functools import partial

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from postgwas._mp_fix import *  # noqa: F401,F403 - multiprocessing bootstrap
from postgwas.cli.compute import get_compute_parser
from postgwas.core.screen_logging import record_screen, resolve_screen_settings
from postgwas.core.ui import format_cli_examples, print_screen_message, run_with_progress
from postgwas.pipeline.registry import REGISTRY, resolve_reference


console = Console()

# Keep the historical misspelling as a compatibility alias without advertising it.
COMMAND_ALIASES = {"pathway_enrichmnet": "pathway_enrichment"}


def _run_recorded_command(operation, screen_settings):
    """Keep unhandled runtime errors inside the terminal/log capture boundary."""
    with record_screen(screen_settings):
        try:
            return operation()
        except Exception as exc:
            print_screen_message("error", "%s: %s" % (type(exc).__name__, exc), console=console)
            return 1
        except KeyboardInterrupt:
            print_screen_message("error", "PostGWAS interrupted by user.", console=console)
            return 130


def print_global_help(prog: str = "postgwas") -> None:
    header = Text("PostGWAS — Unified Toolkit for Post-GWAS Analyses", style="bold cyan")
    console.print(Panel(header, expand=False))

    table = Table(title="Available Modules", title_style="bold magenta", padding=(0, 1))
    table.add_column("Command", style="cyan", no_wrap=True)
    table.add_column("Description", style="green")
    for command, spec in sorted(REGISTRY.commands().items()):
        table.add_row(command, spec.description)
    console.print(table)

    shared_actions = get_compute_parser()._actions
    for title, destinations in (
        ("Shared Screen Options", {"show_screen"}),
        ("Shared Run Options", {"resume", "overwrite"}),
    ):
        option_table = Table(
            title=title,
            title_style="bold magenta",
            padding=(0, 1),
        )
        option_table.add_column("Option", style="cyan", no_wrap=True)
        option_table.add_column("Description", style="green")
        for action in shared_actions:
            if action.dest in destinations:
                option_table.add_row(
                    ", ".join(action.option_strings),
                    Text.from_markup(str(action.help or "")),
                )
        console.print(option_table)
    console.print(
        "[dim]Place these options after a direct command or "
        "pipeline selection.[/dim]"
    )
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

        arguments = argv[2:]
        sys.argv = [f"{prog} --validate"] + arguments
        screen_settings = resolve_screen_settings("harmonisation", arguments)

        def invoke_validation():
            if (
                screen_settings is None
                or screen_settings.configuration is None
                or screen_settings.output_directory is None
            ):
                return validation_main()
            from postgwas.core.direct_execution import run_direct_with_checkpoint
            from postgwas.core.validation_reporting import run_with_file_validation

            return run_direct_with_checkpoint(
                partial(
                    run_with_file_validation, validation_main,
                    command="harmonisation_validation",
                    configuration=screen_settings.configuration,
                    output_directory=screen_settings.output_directory,
                ),
                module_name="harmonisation",
                public_command="harmonisation_validation",
                arguments=arguments,
                configuration=screen_settings.configuration,
                output_directory=screen_settings.output_directory,
                console=console,
            )

        def run_validation():
            if screen_settings is None:
                return invoke_validation()
            return run_with_progress(
                invoke_validation,
                label="Harmonisation validation progress",
                title="Validate the harmonised GWAS-VCF",
                enabled=True,
                console=console,
            )
        return _run_recorded_command(run_validation, screen_settings)

    shorthand_export = argv[1] == "--config"
    command = "config" if shorthand_export else COMMAND_ALIASES.get(argv[1], argv[1])
    commands = REGISTRY.commands()
    if command not in commands:
        print_screen_message("error", f"{prog}: error: unknown module '{argv[1]}'", console=console)
        print_global_help(prog)
        return 2

    spec = commands[command]
    arguments = argv[2:]
    command_prog = f"{prog} --config" if shorthand_export else f"{prog}-{command}"
    sys.argv = [command_prog] + arguments
    screen_settings = resolve_screen_settings(command, arguments)

    def invoke_command():
        try:
            entrypoint = resolve_reference(spec.cli_entrypoint)
        except (ImportError, AttributeError) as exc:
            print_screen_message(
                "error", f"{prog}: module '{command}' could not be loaded: "
                f"{type(exc).__name__}: {exc}", console=console,
            )
            return 1
        if (
            spec.direct_checkpoint != "not_applicable"
            and screen_settings is not None
            and screen_settings.configuration is not None
            and screen_settings.output_directory is not None
        ):
            from postgwas.core.validation_reporting import run_with_file_validation

            entrypoint = partial(
                run_with_file_validation, entrypoint, command=command,
                configuration=screen_settings.configuration,
                output_directory=screen_settings.output_directory,
            )
        if (
            spec.direct_checkpoint != "orchestrated"
            or screen_settings is None
            or screen_settings.configuration is None
            or screen_settings.output_directory is None
        ):
            return entrypoint()
        from postgwas.core.direct_execution import run_direct_with_checkpoint

        return run_direct_with_checkpoint(
            entrypoint,
            module_name=spec.name,
            public_command=command,
            arguments=arguments,
            configuration=screen_settings.configuration,
            output_directory=screen_settings.output_directory,
            console=console,
        )

    def run_command():
        # Direct modules may own analysis-stage progress, but the public command
        # also owns final audit publication. Only the pipeline owns that entire
        # outer boundary itself, including its audit and checkpoints.
        if (
            screen_settings is None
            or spec.owns_top_level_progress
        ):
            return invoke_command()
        else:
            return run_with_progress(
                invoke_command,
                label="%s progress" % command.replace("_", " ").title(),
                title=spec.description,
                enabled=True,
                console=console,
            )
    result = _run_recorded_command(run_command, screen_settings)
    return 0 if result is None else result


if __name__ == "__main__":
    raise SystemExit(main())
