#!/usr/bin/env python3
"""
PostGWAS — Pipeline CLI (v52, Error Deduplication Fix)

Changes:
• FIX: Groups missing arguments to prevent duplicate error messages.
  (e.g., reports '--vcf' once, listing all modules that need it).
"""

import argparse
from collections.abc import Collection
import sys

from rich.console import Console
from rich.table import Table

from postgwas.cli.compute import get_compute_parser, resolve_compute_args
from postgwas.config import load_configuration
from postgwas.config.cli_overrides import get_dotted
from postgwas.config.models.modules.single_cell import (
    SINGLE_CELL_TOOLS,
    single_cell_pipeline_dependencies,
)
from postgwas.core.errors import (
    ConfigurationError,
    MissingRequiredArgumentsError,
    ModuleExecutionError,
    PipelinePlanningError,
)
from postgwas.core.input_validation import (
    InputValidationSession,
    current_validation_session,
    record_file_validation,
    validation_scope,
)
from postgwas.core.required_arguments import (
    RequiredArgument,
    require_resolved_arguments,
)
from postgwas.core.paths import resolve_executable
from postgwas.core.validation_reporting import FileValidationDisplay, file_validation_display
from postgwas.core.preflight import (
    PipelinePreflightEvidence,
    capture_preflight_file_identities,
)
from postgwas.core.ui import (
    AlignedRichHelpFormatter,
    format_cli_examples,
    mark_cli_required_help,
    print_screen_block,
    print_screen_message,
    StageProgress,
)
from postgwas.core.ui.screen import (
    SUMMARY_CARD_LABEL_WIDTH,
    screen_field,
    screen_line,
)
from postgwas.pipeline.planner import (
    build_pipeline_plan,
    resolve_pipeline_dependency_overrides,
)
from postgwas.pipeline.registry import PIPELINE_REQUIRED_OPTIONS, REGISTRY, resolve_reference
from postgwas.pipeline.validation_reporting import (
    save_pipeline_validation_report,
    validation_report_path,
)

# Initialize Rich Console Globally
console = Console()


def _print_error(label, error):
    """Use the same literal error rendering as every direct command."""
    print_screen_message(
        "error", "%s %s" % (label.removeprefix("❌ "), error), console=console,
    )


def _resolve_pipeline_genome_build(args, modules, configuration):
    """Resolve one build for every selected module sharing --genome-build."""
    if hasattr(args, "genome_build"):
        return
    configured = {}
    for module_name in modules:
        path = REGISTRY.get(module_name).genome_build_config_path
        if path is None:
            continue
        value = get_dotted(configuration, path)
        value = getattr(value, "value", value)
        if value is not None:
            configured[path] = value
    distinct = set(configured.values())
    if len(distinct) > 1:
        details = ", ".join(
            "%s=%s" % item for item in sorted(configured.items())
        )
        raise ConfigurationError(
            "Selected pipeline modules configure incompatible genome builds: "
            "%s. Set the listed canonical module keys to one build or provide "
            "--genome-build BUILD to override them for this pipeline." % details
        )
    if distinct:
        args.genome_build = distinct.pop()


def _require_pipeline_arguments(
    args,
    modules,
    configuration,
    *,
    explicitly_requested_modules=(),
):
    """Validate resolved pipeline requirements through the shared contract."""
    checked_modules = list(modules)
    if (
        "formatter" not in explicitly_requested_modules
        and "formatter" in checked_modules
    ):
        checked_modules.remove("formatter")

    requirements = []
    seen = set()
    for module_name in checked_modules:
        for option in REGISTRY.get(module_name).required_options:
            identity = (option.flag, option.config_path)
            if identity in seen:
                continue
            seen.add(identity)
            value = getattr(args, option.dest, None)
            if not value and option.config_path is not None:
                value = get_dotted(configuration, option.config_path)
            requirements.append(RequiredArgument(
                option.flag,
                option.config_path,
                value,
            ))
    if "single_cell" in checked_modules:
        from postgwas.modules.single_cell.service import single_cell_required_arguments

        requirements.extend(single_cell_required_arguments(
            configuration, pipeline=True, args=args,
        ))
    require_resolved_arguments(requirements)


def _validate_pipeline_entry_vcf(args, configuration):
    """Validate the common pipeline VCF before module reference resources."""
    from postgwas.core.vcf import (
        FormattingError,
        validate_harmonised_indexed_vcf,
    )

    module = configuration.modules.formatting
    bcftools = resolve_executable(
        configuration.resources.executables.bcftools,
        "bcftools executable",
        error_type=FormattingError,
    )
    dataset_id = (
        getattr(args, "dataset_id", None) or configuration.run.dataset_id
    )
    evidence = validate_harmonised_indexed_vcf(
        args.vcf,
        dataset_id,
        bcftools,
        module,
    )
    evidence["bcftools_identity"] = capture_preflight_file_identities(
        (bcftools,),
        error_type=FormattingError,
        label="bcftools executable",
    )
    indexed = evidence["indexed"]
    record_file_validation(
        indexed.vcf,
        "Input GWAS-VCF",
        checks=(
            "indexed VCF header",
            "required tag declarations",
            "exactly one sample column",
            "PostGWAS provenance",
        ),
        metrics={
            "variants_from_index": indexed.variant_count,
            "declared_genome_build": indexed.genome_build,
            "sample_count": 1,
            "sample": indexed.sample,
            "contigs": list(indexed.contigs),
        },
        message=(
            "Record-level values are checked by the relevant downstream "
            "readers, not by this header/index check."
        ),
    )
    for index_path, _, _ in indexed.index_identities:
        record_file_validation(index_path, "VCF index", checks=("index available for VCF record counting",))
    return evidence


class PipelinePreflightError(RuntimeError):
    """All independent failed module preflights and the successful evidence."""

    def __init__(self, failures, evidence):
        self.failures = tuple(failures)
        self.evidence = dict(evidence)
        super().__init__("Pipeline input validation failed:\n" + "\n".join(
            "- %s: %s" % (name, error) for name, error in self.failures
        ))


def _run_pipeline_preflights(args, modules, *, initial_evidence=None):
    """Check every selected module, aggregating independent input failures."""
    evidence = dict(initial_evidence or {})
    failures = []
    for module_name in dict.fromkeys(modules):
        preflight = REGISTRY.get(module_name).preflight
        if preflight is None:
            continue
        with validation_scope(module_name):
            try:
                # A failed resource check must not publish partially inferred CLI settings.
                module_args = argparse.Namespace(**vars(args))
                validated = resolve_reference(preflight)(
                    module_args,
                    preflight_evidence=evidence,
                )
                if not isinstance(validated, PipelinePreflightEvidence):
                    raise RuntimeError(
                        "Pipeline preflight for %s did not return the shared "
                        "PipelinePreflightEvidence contract" % module_name
                    )
                if validated.module != module_name:
                    raise RuntimeError(
                        "Pipeline preflight registered for %s returned evidence for %s"
                        % (module_name, validated.module)
                    )
                if current_validation_session() is not None:
                    from postgwas.pipeline.resource_validation import record_pipeline_resource_validation

                    record_pipeline_resource_validation(module_name, validated)
                vars(args).update(vars(module_args))
                evidence[module_name] = validated
                for check in validated.deferred_checks:
                    record_file_validation(None, "Generated-input validation", status="deferred", message=check)
            except (ConfigurationError, OSError, RuntimeError, ValueError) as exc:
                failures.append((module_name, exc))
                record_file_validation(None, "Module preflight", status="failed", message=str(exc))
    if failures:
        raise PipelinePreflightError(failures, evidence)
    return evidence


def _prepare_and_execute_pipeline(args, plan, configuration):
    """Keep startup evidence alive through consumers without checkpointing it."""
    # Required-value validation already accepts canonical YAML. Publish those
    # common values before the audit and runners consume the CLI namespace.
    for option in PIPELINE_REQUIRED_OPTIONS:
        if option.config_path and getattr(args, option.dest, None) is None:
            setattr(args, option.dest, get_dotted(configuration, option.config_path))
    report_path = validation_report_path(args, configuration)
    modules = plan.active_modules
    evidence = {}
    with InputValidationSession() as session:
        display = FileValidationDisplay(session, configuration)
        progress = StageProgress(
            "Pipeline input validation", enabled=True,
            outcome_label_width=configuration.logging.terminal_label_width,
        )
        error = None
        try:
            # One active operation: no invented within-file completion percentage.
            with progress.step(1, 1, "Validate pipeline inputs and reference resources"):
                with session.scope("pipeline"):
                    try:
                        entry = _validate_pipeline_entry_vcf(args, configuration)
                    except (ConfigurationError, OSError, RuntimeError, ValueError) as exc:
                        record_file_validation(getattr(args, "vcf", None), "Input GWAS-VCF", status="failed", message=str(exc))
                        for name in modules:
                            with session.scope(name):
                                record_file_validation(None, "Module preflight", status="blocked", message="Entry GWAS-VCF validation failed.")
                        raise
                evidence = {"input_vcf": entry, "current_vcf": entry}
                evidence = _run_pipeline_preflights(args, modules, initial_evidence=evidence)
                save_pipeline_validation_report(session, evidence, modules, report_path)
        except BaseException as exc:
            error = exc
            if isinstance(exc, PipelinePreflightError):
                evidence = exc.evidence
            elif isinstance(exc, (ConfigurationError, OSError, RuntimeError, ValueError)):
                raise PipelinePreflightError((("pipeline", exc),), evidence) from exc
            raise
        finally:
            progress.close()
            try:
                if error is not None:
                    save_pipeline_validation_report(session, evidence, modules, report_path, error=error)
                display.flush(report_path=report_path)
            except Exception as report_error:
                if error is None:
                    raise
                _print_error("Validation audit could not be saved:", report_error)
        _print_pipeline_preflight_summary(entry, evidence, modules)

        def finalize_validation():
            display.flush()
            save_pipeline_validation_report(
                session, evidence, modules, report_path, phase="pipeline_execution",
            )

        try:
            resolve_compute_args(args)
            with file_validation_display(session, configuration, display=display):
                return execute_pipeline(
                    args, plan, configuration, preflight_evidence=evidence,
                    finalize_validation=finalize_validation,
                )
        except BaseException as exc:
            record_file_validation(None, "Pipeline execution", status="failed", message=str(exc))
            display.flush()
            try:
                save_pipeline_validation_report(
                    session, evidence, modules, report_path,
                    error=exc, phase="pipeline_execution",
                )
            except Exception as report_error:
                _print_error("Final validation audit could not be saved:", report_error)
            raise


def _print_pipeline_preflight_summary(input_vcf_evidence, evidence, modules):
    """Show one concise readiness card after ordered pipeline preflight."""
    validated = [
        evidence[name]
        for name in dict.fromkeys(modules)
        if isinstance(evidence.get(name), PipelinePreflightEvidence)
    ]

    def has_resources(item):
        resources = item.resources
        if resources is None:
            return False
        if isinstance(resources, Collection):
            return len(resources) > 0
        return True

    resource_preflights = sum(has_resources(item) for item in validated)
    deferred_checks = sum(len(item.deferred_checks) for item in validated)
    print_screen_block("\n".join((
        "",
        screen_line("analysis", "Pipeline readiness", indent=2),
        screen_field(
            "success", "Stage preflights",
            "%d/%d passed" % (len(validated), len(tuple(dict.fromkeys(modules)))),
            indent=6, label_width=SUMMARY_CARD_LABEL_WIDTH,
        ),
        screen_field(
            "success", "External-resource checks",
            "%d resource-dependent stages passed" % resource_preflights,
            indent=6, label_width=SUMMARY_CARD_LABEL_WIDTH,
        ),
        screen_field(
            "decision", "Generated-input checks",
            "%d deferred to their consuming stages" % deferred_checks,
            indent=6, label_width=SUMMARY_CARD_LABEL_WIDTH,
        ),
        "",
    )))

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


def _add_pipeline_selection_arguments(container):
    """Add the canonical workflow selectors to a parser or argument group."""
    container.add_argument(
        "--modules",
        nargs="*",
        choices=REGISTRY.names(),
        metavar="MODULE",
        help="One or more final analyses to run, for example finemap magma pops.",
    )
    container.add_argument(
        "--apply-filter",
        action="store_true",
        help="Run quality-control filtering before downstream analyses.",
    )
    container.add_argument(
        "--apply-imputation",
        action="store_true",
        help="Impute missing summary statistics before downstream analyses.",
    )
    container.add_argument(
        "--apply-manhattan",
        action="store_true",
        help="Create a Manhattan plot.",
    )
    container.add_argument(
        "--heritability",
        action="store_true",
        help="Estimate SNP heritability with LDSC.",
    )


def print_full_pipeline_help(modules, parser, *, error=False):
    if error:
        print_screen_message(
            "error", "The pipeline cannot start yet.\n"
            "Add the missing values shown above, then run the command again.",
            console=console,
        )
    selection_only = getattr(
        parser, "_postgwas_pipeline_selection_only", False,
    )
    if selection_only:
        console.print("[bold cyan]Choose the fine-mapping workflow[/bold cyan]\n")
    else:
        console.print("[bold cyan]Steps that will run[/bold cyan]\n")
    description_overrides = getattr(
        parser, "_postgwas_pipeline_step_descriptions", {},
    )
    name_overrides = getattr(parser, "_postgwas_pipeline_step_names", {})
    if not selection_only:
        for idx, m in enumerate(modules, 1):
            if m in REGISTRY.names(include_internal=True):
                console.print(
                    f" {idx}) [cyan]{name_overrides.get(m, m)}[/cyan]"
                )
                description = description_overrides.get(
                    m, REGISTRY.get(m).description,
                )
                console.print(f"      • {description}")
            print("")
    heading = (
        "Workflow selection options"
        if selection_only
        else "Options for these steps"
    )
    console.print("[bold cyan]%s[/bold cyan]\n" % heading)
    parser.print_help()

# =====================================================================
# MAIN
# =====================================================================

def main():

    # 1. Minimal Parse
    mini = argparse.ArgumentParser(
        prog="postgwas pipeline",
        add_help=False,
        formatter_class=AlignedRichHelpFormatter,
        parents=[get_compute_parser()],
    )
    _add_pipeline_selection_arguments(mini)
    mini.add_argument(
        "--tools",
        nargs="+",
        choices=SINGLE_CELL_TOOLS,
        metavar="TOOL",
        help="Single-cell integration tools used when selecting single_cell.",
    )
    mini.add_argument(
        "--run-config",
        metavar="PATH",
        help="YAML file containing pipeline and module settings.",
    )
    mini.add_argument(
        "-h", "--help", action="store_true",
        help="Show this help message and exit.",
    )

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
                    (
                        "--modules finemap",
                        "--clumping-methods standard",
                        "--finemap-method susie",
                        "--help",
                    ),
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
        if a1.help:
            console.print("\n[bold cyan]Shared pipeline options[/bold cyan]\n")
            mini.print_help()
        sys.exit(0)

    # 3. Validate targets and build the execution plan exactly once.
    raw_modules = list(a1.modules or [])
    try:
        configuration = load_configuration(a1.run_config)
        dependency_overrides = resolve_pipeline_dependency_overrides(
            tuple(sys.argv[1:]), configuration,
        )
        if "single_cell" in raw_modules:
            tools = a1.tools or configuration.modules.single_cell.tools
            dependency_overrides["single_cell"] = (
                single_cell_pipeline_dependencies(tools)
            )
        plan = build_pipeline_plan(
            raw_modules,
            apply_filter=a1.apply_filter,
            apply_imputation=a1.apply_imputation,
            apply_manhattan=a1.apply_manhattan,
            heritability=a1.heritability,
            dependency_overrides=dependency_overrides or None,
        )
    except (ConfigurationError, PipelinePlanningError) as exc:
        _print_error("Pipeline Planning Error:", exc)
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
    parser = HelpOnErrorArgumentParser(
        prog="postgwas pipeline",
        formatter_class=AlignedRichHelpFormatter,
        parents=parent_parsers,
        conflict_handler='resolve',
        usage=custom_usage,
        epilog=None,
    )
    parser._postgwas_help_arguments = tuple(sys.argv[1:])

    if len(plan.requested_modules) == 1:
        customizer = REGISTRY.get(
            plan.requested_modules[0]
        ).single_target_help_customizer
        if customizer:
            try:
                resolve_reference(customizer)(parser)
            except ConfigurationError as exc:
                _print_error(
                    "❌ Configuration Error:", exc,
                )
                sys.exit(2)

    workflow_examples = list(getattr(
        parser,
        "_postgwas_pipeline_examples",
        workflow_examples,
    ))
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
    parser.epilog = format_cli_examples(
        *workflow_examples,
        notes=(
            "The exported configuration includes required preceding modules in execution order.",
        ),
        title="Workflow examples",
    )

    run_config_action = next(
        (action for action in parser._actions if action.dest == "run_config"),
        None,
    )
    if run_config_action is None:
        configuration_group = parser.add_argument_group("Configuration")
        configuration_group.add_argument(
            "--run-config",
            default=argparse.SUPPRESS,
            metavar="PATH",
            help=(
                "YAML file containing settings and reference resources for every "
                "selected module."
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

    mark_cli_required_help(
        parser,
        (
            option.dest
            for module_name in parser_modules
            for option in REGISTRY.get(module_name).required_options
        ),
    )
    if "single_cell" in parser_modules:
        from postgwas.modules.single_cell.service import single_cell_required_arguments

        mark_cli_required_help(parser, (
            item.option[2:].replace("-", "_")
            for item in single_cell_required_arguments(
                configuration, pipeline=True, args=a1,
            )
        ))

    # Re-add CLI Flags
    workflow = parser.add_argument_group("Choose analyses")
    _add_pipeline_selection_arguments(workflow)

    selection_only = getattr(
        parser, "_postgwas_pipeline_selection_only", False,
    )
    selection_visible = getattr(
        parser, "_postgwas_pipeline_selection_visible_destinations", None,
    )
    if selection_only and selection_visible is not None:
        for action in parser._actions:
            if action.dest not in selection_visible:
                action.required = False
                action.help = argparse.SUPPRESS
    selection_error = getattr(
        parser, "_postgwas_pipeline_selection_error", None,
    )
    if selection_only and not a1.help:
        if selection_error:
            _print_error(
                "❌ Fine-mapping workflow selection error:",
                selection_error,
            )
        else:
            print_screen_message(
                "warning", "Fine-mapping workflow selection is required.\n"
                "Choose both methods below, or provide them through --run-config.",
                console=console,
            )
        print_full_pipeline_help(execution_modules, parser)
        console.print("\n[dim]No analysis was started.[/dim]\n")
        sys.exit(2)

    if a1.help:
        if selection_error:
            _print_error(
                "❌ Fine-mapping workflow selection error:",
                selection_error,
            )
        print_full_pipeline_help(execution_modules, parser)
        sys.exit(2 if selection_error else 0)

    # 5. Parse & Validate
    try:
        args = parser.parse_args()
    except ValueError as e:
        _print_error("❌ Argument Parsing Error:", e)
        print_full_pipeline_help(execution_modules, parser, error=True)
        sys.exit(2)

    try:
        resolved_configuration = load_configuration(
            getattr(args, "run_config", None)
        )
        _resolve_pipeline_genome_build(
            args, parser_modules, resolved_configuration,
        )
    except ConfigurationError as exc:
        _print_error("❌ Configuration Error:", exc)
        sys.exit(2)
    run_defaults = resolved_configuration.run
    if not hasattr(args, "resume"):
        args.resume = run_defaults.resume
    if not hasattr(args, "overwrite"):
        args.overwrite = run_defaults.overwrite

    try:
        _require_pipeline_arguments(
            args,
            parser_modules,
            resolved_configuration,
            explicitly_requested_modules=raw_modules,
        )
    except MissingRequiredArgumentsError as exc:
        _print_error(
            "❌ Missing Required Arguments:", exc,
        )
        print_full_pipeline_help(execution_modules, parser, error=True)
        sys.exit(2)

    args._pipeline_requested_modules = plan.requested_modules
    args.modules = parser_modules
    # 6. Execute
    try:
        _prepare_and_execute_pipeline(args, plan, resolved_configuration)

    except PipelinePreflightError as exc:
        _print_error("❌ Pipeline preflight failed:", exc)
        sys.exit(2)

    # ✅ CASE 1: Clean Stop (sys.exit from within the pipeline)
    except SystemExit as e:
        # Just exit with the same code.
        # Do NOT print help. This handles your "Empty File" stop cleanly.
        sys.exit(e.code)
    except KeyboardInterrupt:
        print_screen_message("error", "Pipeline interrupted by user.", console=console)
        sys.exit(130)

    # ✅ CASE 2: Known Configuration/Usage Errors
    except ModuleExecutionError as e:
        _print_error("Pipeline failed:", e)
        sys.exit(1)
    except (ValueError, OSError, AttributeError, TypeError) as e:
        err_msg = str(e).lower()

        # Only print the Help Menu if the user is missing a tool or argument
        if "required" in err_msg or "executable not found" in err_msg:
             _print_error("❌ Configuration Error:", e)
             print_full_pipeline_help(execution_modules, parser)
             sys.exit(2)

        # For other Python errors (bugs), just print the error, no help menu spam
        _print_error("❌ Pipeline Failed:", e)
        sys.exit(1)

    # ✅ CASE 3: Unknown Crashes
    except Exception as e:
        _print_error("❌ Critical Unexpected Error:", e)
        sys.exit(1)

    print_screen_message("success", "Pipeline complete.", console=console)

if __name__ == "__main__":
    main()
