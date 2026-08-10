#!/usr/bin/env python3
"""
PostGWAS — Pipeline CLI (v52, Error Deduplication Fix)

Changes:
• FIX: Groups missing arguments to prevent duplicate error messages.
  (e.g., reports '--vcf' once, listing all modules that need it).
"""

import argparse
import sys
from rich.console import Console
from rich.markup import escape
from rich.table import Table
from postgwas.core.errors import (
    ConfigurationError,
    ModuleExecutionError,
    PipelinePlanningError,
)
from postgwas.core.ui import (
    AlignedRichHelpFormatter,
    format_cli_examples,
    help_with_default,
)
from postgwas.config import load_configuration
from postgwas.pipeline.planner import build_pipeline_plan
from postgwas.pipeline.registry import REGISTRY, resolve_reference
from postgwas.cli.compute import resolve_compute_args

# Initialize Rich Console Globally
console = Console()


def _print_error(label, error, *, leading_newline=False):
    """Render an exception as literal text while retaining the label style."""
    prefix = "\n" if leading_newline else ""
    console.print(
        "%s[bold red]%s[/bold red] %s"
        % (prefix, label, escape(str(error)))
    )

# =====================================================================
# IMPORT ORCHESTRATOR
# =====================================================================
from postgwas.pipeline.executor import execute_pipeline

# =====================================================================
# CUSTOM PARSER
# =====================================================================

class HelpOnErrorArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        raise ValueError(message)

def print_full_pipeline_help(modules, parser, *, error=False):
    if error:
        console.print("\n[bold red]The pipeline cannot start yet.[/bold red]")
        console.print("Add the missing values shown above, then run the command again.\n")
    console.print("[bold cyan]Steps that will run[/bold cyan]\n")
    for idx, m in enumerate(modules, 1):
        if m in REGISTRY.names(include_internal=True):
            console.print(f" {idx}) [cyan]{m}[/cyan]")
            console.print(f"      • {REGISTRY.get(m).description}")
        print("")
    console.print("[bold cyan]Options for these steps[/bold cyan]\n")
    parser.print_help()

# =====================================================================
# MAIN
# =====================================================================

def main():

    # 1. Minimal Parse
    mini = argparse.ArgumentParser(add_help=False)
    mini.add_argument("--modules", nargs="*", choices=REGISTRY.names())
    mini.add_argument("--apply-filter", action="store_true")
    mini.add_argument("--apply-imputation", action="store_true")
    mini.add_argument("--apply-manhattan", action="store_true")
    mini.add_argument("--heritability", action="store_true")
    mini.add_argument("-h", "--help", action="store_true")

    a1, _ = mini.parse_known_args()

    # 2. Show Table if Empty
    if not a1.modules and not (
        a1.apply_filter
        or a1.apply_imputation
        or a1.apply_manhattan
        or a1.heritability
    ):
        console.print("\nChoose the final analysis you want PostGWAS to run.\n")
        t = Table(title="PostGWAS Pipeline — Available Modules", header_style="bold cyan", show_lines=True)
        t.add_column("Module", style="magenta", no_wrap=True)
        t.add_column("Description", style="white")
        for name in REGISTRY.names():
            spec = REGISTRY.get(name)
            status = "" if spec.pipeline_enabled else " [unavailable in pipeline]"
            t.add_row(name, spec.description + status)
        console.print(t)
        console.print()
        console.print(
            format_cli_examples(
                (
                    "Inspect every option required for a fine-mapping workflow:",
                    "postgwas pipeline",
                    ("--modules finemap", "--apply-filter", "--help"),
                ),
                (
                    "Export a complete reusable fine-mapping configuration:",
                    "postgwas config export",
                    (
                        "--pipeline finemap",
                        "--style full",
                        "--output finemap_pipeline.yaml",
                    ),
                ),
            )
        )
        sys.exit(0)

    # 3. Validate targets and build the execution plan exactly once.
    raw_modules = list(a1.modules or [])
    try:
        plan = build_pipeline_plan(
            raw_modules,
            apply_filter=a1.apply_filter,
            apply_imputation=a1.apply_imputation,
            apply_manhattan=a1.apply_manhattan,
            heritability=a1.heritability,
        )
    except PipelinePlanningError as exc:
        console.print(f"\n❌ [bold red]Pipeline Planning Error:[/bold red] {exc}")
        sys.exit(2)

    execution_modules = list(plan.steps)
    parser_modules = list(plan.active_modules)

    # 4. Build Full Parser
    parent_parsers = []
    pipeline_supplied_destinations = set()
    seen = set()
    for m in parser_modules:
        pipeline_supplied_destinations.update(
            REGISTRY.get(m).pipeline_supplied_options
        )
        for reference in REGISTRY.get(m).parser_factories:
            if reference not in seen:
                parent = resolve_reference(reference)()
                parent_parsers.append(parent)
                seen.add(reference)

    custom_usage = "postgwas pipeline [options] --modules " + " ".join(
        plan.requested_modules
    )
    selected_modules = " ".join(plan.requested_modules)
    workflow_examples = []
    if len(plan.requested_modules) == 1:
        example_factory = REGISTRY.get(
            plan.requested_modules[0]
        ).pipeline_example_factory
        if example_factory:
            workflow_examples.extend(resolve_reference(example_factory)())
    if not workflow_examples:
        workflow_examples.append(
            (
                "Run this configured workflow:",
                "postgwas pipeline",
                (
                    "--modules %s" % selected_modules,
                    "--vcf study.vcf.gz",
                    "--dataset-id STUDY",
                    "--output-directory results",
                    "--run-config analysis_pipeline.yaml",
                ),
            )
        )
    workflow_examples.append(
        (
            "Export a complete configuration for these analyses:",
            "postgwas config export",
            (
                "--pipeline %s" % selected_modules,
                "--style full",
                "--output analysis_pipeline.yaml",
            ),
        )
    )
    help_sections = [
        format_cli_examples(
            *workflow_examples,
            notes=(
                "The exported configuration includes required preceding modules in execution order.",
            ),
            title="Workflow examples",
        )
    ]
    parser = HelpOnErrorArgumentParser(
        prog="postgwas pipeline",
        formatter_class=AlignedRichHelpFormatter,
        parents=parent_parsers,
        conflict_handler='resolve',
        usage=custom_usage,
        epilog="\n\n".join(help_sections) or None,
    )

    run_config_action = next(
        (action for action in parser._actions if action.dest == "run_config"),
        None,
    )
    if run_config_action is None:
        configuration = parser.add_argument_group("Configuration")
        configuration.add_argument(
            "--run-config",
            default=argparse.SUPPRESS,
            metavar="PATH",
            help=(
                "YAML file containing settings and reference resources for every "
                "selected module."
            ),
        )
    pipeline_controls = parser.add_argument_group("Pipeline continuation")
    pipeline_controls.add_argument(
        "--resume",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Continue after provenance-validated completed steps",
            load_configuration().run.resume,
        ),
    )
    pipeline_controls.add_argument(
        "--overwrite",
        action="store_true",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Rerun steps and replace their existing outputs; this takes precedence "
            "over resume",
            load_configuration().run.overwrite,
        ),
    )
    # -------------------------------------------------------------
    # Context-Sensitive Arguments
    # -------------------------------------------------------------

    if "formatter" not in raw_modules:
        pipeline_supplied_destinations.add("format")

    for action in parser._actions:
        if action.dest in pipeline_supplied_destinations:
            action.required = False
            action.help = argparse.SUPPRESS

    # Re-add CLI Flags
    workflow = parser.add_argument_group("Choose analyses")
    workflow.add_argument(
        "--modules",
        nargs="*",
        metavar="MODULE",
        help="One or more final analyses to run, for example finemap magma pops.",
    )
    workflow.add_argument(
        "--apply-filter",
        action="store_true",
        help="Run quality-control filtering before downstream analyses.",
    )
    workflow.add_argument(
        "--apply-imputation",
        action="store_true",
        help="Impute missing summary statistics before downstream analyses.",
    )
    workflow.add_argument(
        "--apply-manhattan",
        action="store_true",
        help="Create Manhattan and QQ plots.",
    )
    workflow.add_argument(
        "--heritability",
        action="store_true",
        help="Estimate SNP heritability with LDSC.",
    )

    if a1.help:
        print_full_pipeline_help(execution_modules, parser)
        sys.exit(0)

    # 5. Parse & Validate
    try:
        args = parser.parse_args()
    except ValueError as e:
        _print_error("❌ Argument Parsing Error:", e, leading_newline=True)
        print_full_pipeline_help(execution_modules, parser, error=True)
        sys.exit(2)

    if not hasattr(args, "resume"):
        args.resume = load_configuration(
            getattr(args, "run_config", None)
        ).run.resume

    # -----------------------------------------------------------
    # NEW LOGIC: Deduplicate Error Messages
    # -----------------------------------------------------------
    missing_args_map = {} # arg_name -> list of modules

    check_list = parser_modules.copy()
    if "formatter" not in raw_modules and "formatter" in check_list:
        check_list.remove("formatter")

    for m in check_list:
        for option in REGISTRY.get(m).required_options:
            val = getattr(args, option.dest, None)
            if not val:
                if option.flag not in missing_args_map:
                    missing_args_map[option.flag] = []
                missing_args_map[option.flag].append(m)

    if missing_args_map:
        console.print("\n❌ [bold red]Missing Required Arguments:[/bold red]")
        for arg, modules in missing_args_map.items():
            mod_list = ", ".join([f"[cyan]{m}[/cyan]" for m in modules])
            console.print(f"   • Argument [bold red]{arg}[/bold red] is required by module(s): {mod_list}")

        print_full_pipeline_help(execution_modules, parser, error=True)
        sys.exit(2)

    try:
        for module_name in parser_modules:
            preflight = REGISTRY.get(module_name).preflight
            if preflight:
                resolve_reference(preflight)(args)
    except (ConfigurationError, OSError, RuntimeError, ValueError) as exc:
        _print_error("❌ Pipeline preflight failed:", exc, leading_newline=True)
        sys.exit(2)

    # 6. Execute
    try:
        args.modules = parser_modules
        resolve_compute_args(args)
        execute_pipeline(args, plan)

    # ✅ CASE 1: Clean Stop (sys.exit from within the pipeline)
    except SystemExit as e:
        # Just exit with the same code.
        # Do NOT print help. This handles your "Empty File" stop cleanly.
        sys.exit(e.code)

    # ✅ CASE 2: Known Configuration/Usage Errors
    except ModuleExecutionError as e:
        _print_error("Pipeline failed:", e)
        sys.exit(1)
    except (ValueError, OSError, AttributeError, TypeError) as e:
        err_msg = str(e).lower()

        # Only print the Help Menu if the user is missing a tool or argument
        if "required" in err_msg or "executable not found" in err_msg:
             _print_error("❌ Configuration Error:", e, leading_newline=True)
             print_full_pipeline_help(execution_modules, parser)
             sys.exit(2)

        # For other Python errors (bugs), just print the error, no help menu spam
        _print_error("❌ Pipeline Failed:", e, leading_newline=True)
        sys.exit(1)

    # ✅ CASE 3: Unknown Crashes
    except Exception as e:
        _print_error("❌ Critical Unexpected Error:", e, leading_newline=True)
        sys.exit(1)

    console.print("\n🎉 Pipeline complete.\n")

if __name__ == "__main__":
    main()
