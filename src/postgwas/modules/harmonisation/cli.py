"""Standalone command-line boundary for GWAS harmonisation."""

from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
import traceback
from pathlib import Path
from typing import Any

from postgwas.config import load_configuration, write_resolved_configuration
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.config.merger import deep_merge
from postgwas.cli.compute import get_compute_parser
from postgwas.core.errors import ConfigurationError
from postgwas.core.paths import configured_output_path
from postgwas.modules.harmonisation.sample_sheet import (
    HarmonisationSampleSheetRow,
    load_harmonisation_sample_sheet,
    to_harmonisation_input,
)
from postgwas.modules.harmonisation.policies import POLICIES, load_policies
from postgwas.core.pipeline_logging import write_log_record as _append_run_log
from postgwas.modules.harmonisation.service import (
    ConfigError,
    PipelineError,
    run_harmonisation_pipeline,
)
from postgwas.modules.harmonisation.concordance.errors import (
    ConcordanceValidationError,
)
from postgwas.core.ui import (
    AlignedRichHelpFormatter,
    format_cli_default,
    format_cli_examples,
    help_with_default,
    screen_field,
    screen_line,
)
from postgwas.modules.harmonisation.vcf_processing import require_binaries


_RESOLVED_RESOURCE_PATHS = (
    "root",
    "executables.bash",
    "executables.bcftools",
    "executables.python",
    "executables.tabix",
    "executables.pigz",
)


CLI_OVERRIDE_PATHS = {
    "resource_directory": "resources.root",
    "output_directory": "run.output_directory",
    "comparison_af_source": "modules.harmonisation.comparison_af.source",
    "comparison_af_column": "modules.harmonisation.comparison_af.column",
    "threads": "execution.threads",
    "memory_gb": "execution.memory_gb",
    "seed": "execution.random_seed",
    "validate": "modules.harmonisation.concordance_validation.enabled",
    "fixed_info": "modules.harmonisation.fixed_info.value",
    "zero_p_se_action": (
        "modules.harmonisation.policies.pvalue.zero_missing_se"
    ),
}


def _unit_interval_float(value: str) -> float:
    """Argparse type for a finite score in the closed unit interval."""
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError("must be a number between 0 and 1") from exc
    if not math.isfinite(number) or not 0 <= number <= 1:
        raise argparse.ArgumentTypeError("must be between 0 and 1 inclusive")
    return number


def get_harmonisation_parser(add_help: bool = False) -> argparse.ArgumentParser:
    """Build a parser whose defaults are displayed but never owned by argparse."""
    defaults = load_configuration()
    comparison = defaults.modules.harmonisation.comparison_af
    policy_defaults = load_policies(defaults.modules.harmonisation.policies)
    zero_p_se_policy = POLICIES["pvalue.zero_missing_se"]
    parser = argparse.ArgumentParser(
        add_help=add_help,
        formatter_class=AlignedRichHelpFormatter,
        parents=[get_compute_parser()],
    )
    files = parser.add_argument_group("Files and folders")
    files.add_argument(
        "--sample-sheet",
        default=argparse.SUPPRESS,
        metavar="PATH",
        help=(
            "REQUIRED. CSV or TSV file listing each GWAS dataset and the names "
            "of its columns. See examples/configs/harmonisation for templates."
        ),
    )
    files.add_argument(
        "--dataset-id",
        default=argparse.SUPPRESS,
        metavar="ID",
        help=(
            "Process only this dataset_id from the sample sheet. Without this option, "
            "all rows are processed."
        ),
    )
    files.add_argument(
        "--run-config",
        default=argparse.SUPPRESS,
        metavar="PATH",
        help=(
            "Optional YAML file with settings for this run. Command-line values "
            "override values in this file."
        ),
    )
    files.add_argument(
        "--resource-directory",
        default=argparse.SUPPRESS,
        metavar="PATH",
        help=(
            "Folder containing the PostGWAS reference data. REQUIRED unless "
            "resources.root is set in --run-config."
        ),
    )
    files.add_argument(
        "--output-directory",
        default=argparse.SUPPRESS,
        metavar="PATH",
        help=help_with_default(
            "Folder where results, logs, and the resolved configuration will be saved",
            defaults.run.output_directory,
        ),
    )
    frequency = parser.add_argument_group("Reference frequency check")
    frequency.add_argument(
        "--comparison-af-source",
        default=argparse.SUPPRESS,
        metavar="SOURCE",
        choices=tuple(comparison.available_sources),
        help=(
            "Reference dataset used to check allele frequencies; it is not used "
            "to fill missing study frequencies.\n"
            "Options: %s.\n"
            "The resource directory must contain one VCF per chromosome:\n"
            "%s\n"
            "%s."
            % (
                ", ".join(comparison.available_sources),
                "\n".join(
                    "  %s: %s" % (source, comparison.resource_examples[source])
                    for source in comparison.available_sources
                ),
                format_cli_default(comparison.source),
            )
        ),
    )
    validation = parser.add_argument_group("Optional accuracy validation")
    validation.add_argument(
        "--validate",
        action="store_true",
        default=argparse.SUPPRESS,
        help=(
            "After harmonisation, compare each input dataset with its automatically "
            "selected same-build merged VCF. Disabled by default."
        ),
    )
    frequency.add_argument(
        "--comparison-af-column",
        default=argparse.SUPPRESS,
        metavar="COLUMN",
        help=help_with_default(
            "Population column to use in that reference, such as EUR, AFR, or EAS",
            defaults.modules.harmonisation.comparison_af.column,
        ),
    )
    quality = parser.add_argument_group("Imputation quality")
    quality.add_argument(
        "--fixed-info",
        type=_unit_interval_float,
        default=argparse.SUPPRESS,
        metavar="VALUE",
        help=(
            "Explicitly assign one INFO value to every variant only for datasets "
            "that provide neither an internal INFO column nor an external INFO "
            "file. Must be between 0 and 1. Internal INFO has first priority and "
            "external INFO has second priority. The assigned value is logged as "
            "user-provided, not measured per variant."
        ),
    )
    statistics = parser.add_argument_group("Effect reconstruction")
    statistics.add_argument(
        "--zero-p-se-action",
        dest="zero_p_se_action",
        choices=zero_p_se_policy.validator.members,
        default=argparse.SUPPRESS,
        metavar="ACTION",
        help=help_with_default(
            "What to do when a raw p-value is exactly zero and SE cannot be "
            "obtained from supplied SE or a usable Z score: 'fail' stops without "
            "discarding the variants, 'approximate' derives SE using the configured "
            "p-value floor and records the approximation, and 'reject' writes the "
            "variants to the rejected-variant output",
            policy_defaults.get("pvalue.zero_missing_se"),
        ),
    )

    return parser


def _module_policy_block(config) -> dict[str, Any]:
    harmonisation = config.modules.harmonisation
    configured = {
        "execution": {
            "total_cpu_budget": config.execution.threads,
        },
        "logging": {"level": config.logging.file_level},
    }
    return deep_merge(configured, harmonisation.policies)


def _harmonisation_executables(config) -> dict[str, str]:
    executables = config.resources.executables
    return {
        "bash": executables.bash,
        "bcftools": executables.bcftools,
        "python": executables.python,
        "tabix": executables.tabix,
    }


def _engine_defaults(config, resolved_executables=None) -> dict[str, Any]:
    """Build the scientific engine configuration from the resolved run config."""
    harmonisation = config.modules.harmonisation
    mapping = harmonisation.external_eaf_mapping

    def reference_mapping(value) -> dict[str, str]:
        return {
            "chr": value.chromosome,
            "pos": value.position,
            "a1": value.effect_allele,
            "a2": value.other_allele,
            "delimiter": value.delimiter,
        }

    policy_block = _module_policy_block(config)
    return {
        "executables": dict(
            resolved_executables
            if resolved_executables is not None
            else _harmonisation_executables(config)
        ),
        "executables_prevalidated": resolved_executables is not None,
        "compression_executable": config.resources.executables.pigz,
        "default_comparison_af_file": harmonisation.comparison_af.source,
        "default_comparison_af_column": harmonisation.comparison_af.column,
        "default_dbsnp": harmonisation.reference.dbsnp_source,
        "default_eaf_colmap": reference_mapping(
            harmonisation.default_eaf_mapping
        ),
        "external_eaf_colmap": reference_mapping(mapping),
        "external_info_colmap": reference_mapping(
            harmonisation.external_info_mapping
        ),
        "build_check_colmap": reference_mapping(
            harmonisation.build_check_mapping
        ),
        "resource_layout": harmonisation.resource_layout.model_dump(),
        "output_layout": dict(harmonisation.output_layout.root),
        "gwas2vcf_input": harmonisation.gwas2vcf_input.model_dump(),
        "vcf_processing": harmonisation.vcf_processing.model_dump(),
        "population_frequency_qc": harmonisation.population_frequency_qc.model_dump(),
        "policies": policy_block,
    }


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def _runtime_path(config, root, setting: str, **values) -> Path:
    pattern = getattr(config.modules.harmonisation.runtime, setting)
    return configured_output_path(root, pattern, **values)


def _top_run_log(config, root) -> Path:
    metadata = _runtime_path(config, root, "top_metadata_directory")
    return configured_output_path(metadata, config.logging.filename)


def _fixed_info_cli_value(args, config) -> float | None:
    """Return the explicit fallback while rejecting run-config activation."""
    if hasattr(args, "fixed_info"):
        return float(args.fixed_info)
    if config.modules.harmonisation.fixed_info.value is not None:
        raise ConfigurationError(
            "fixed_info.value cannot be set in a run configuration. Supply "
            "the explicit command-line option --fixed-info VALUE instead."
        )
    return None


def _validate_info_fallback(
    rows: list[HarmonisationSampleSheetRow], fixed_info: float | None,
) -> None:
    """Fail preflight when any selected dataset has no usable INFO source."""
    missing = [
        row.dataset_id
        for row in rows
        if row.imputation_info_column is None
        and not (row.external_info_file and row.external_info_column)
    ]
    if missing and fixed_info is None:
        raise ConfigurationError(
            "No imputation-quality source is available for dataset(s): %s. "
            "Provide imputation_info_column, provide external_info_file together "
            "with external_info_column, or explicitly supply --fixed-info VALUE."
            % ", ".join(missing)
        )


def _write_run_metadata(config, args, row: HarmonisationSampleSheetRow) -> Path:
    root = Path(config.run.output_directory).expanduser().resolve()
    metadata = _runtime_path(
        config, root, "metadata_directory", dataset_id=row.dataset_id,
    )
    metadata.mkdir(parents=True, exist_ok=True)
    write_resolved_configuration(
        config,
        _runtime_path(config, metadata, "resolved_config_file"),
        modules=("harmonisation",),
        resource_paths=_RESOLVED_RESOURCE_PATHS,
    )
    _write_text(
        _runtime_path(config, metadata, "sample_sheet_row_file"),
        json.dumps(row.model_dump(mode="json"), indent=2, sort_keys=True) + "\n",
    )
    _write_text(
        _runtime_path(config, metadata, "command_file"),
        " ".join(sys.argv) + "\n",
    )
    if hasattr(args, "run_config"):
        shutil.copyfile(
            args.run_config,
            _runtime_path(config, metadata, "supplied_config_file"),
        )
    return metadata


def _run_validated_rows(
    args,
    config,
    rows: list[HarmonisationSampleSheetRow],
    resolved_executables=None,
):
    fixed_info_from_cli = _fixed_info_cli_value(args, config)
    _validate_info_fallback(rows, fixed_info_from_cli)
    resource_root = config.resources.root
    if resource_root is None:
        raise ConfigurationError(
            "resources.root is required; supply --resource-directory or resources.root in --run-config"
        )
    resource_root = Path(resource_root).expanduser().resolve()
    if not resource_root.is_dir():
        raise ConfigurationError("Resource directory does not exist or is not a directory: %s" % resource_root)
    output_root = Path(config.run.output_directory).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    top_log = _top_run_log(config, output_root)
    _append_run_log(top_log, "INFO", "Resolved sample sheet with %d dataset(s)." % len(rows))
    top_metadata = _runtime_path(
        config, output_root, "top_metadata_directory",
    )
    write_resolved_configuration(
        config,
        _runtime_path(config, top_metadata, "resolved_config_file"),
        modules=("harmonisation",),
        resource_paths=_RESOLVED_RESOURCE_PATHS,
    )

    defaults = _engine_defaults(config, resolved_executables=resolved_executables)
    results: dict[str, Any] = {}
    failures = []
    for row in rows:
        metadata = _write_run_metadata(config, args, row)
        dataset_log = metadata / config.logging.filename
        for warning in row.normalisation_warnings:
            _append_run_log(dataset_log, "WARNING", warning)
            print(screen_field(
                "warning", "Sample-sheet value", warning,
                indent=4, label_width=20,
            ))
        engine_input = to_harmonisation_input(
            row,
            resource_directory=resource_root,
            output_directory=output_root,
            output_layout=defaults["output_layout"],
            fixed_info=fixed_info_from_cli,
            fixed_info_column=config.modules.harmonisation.fixed_info.column,
        )
        if fixed_info_from_cli is not None and engine_input["info_source"] != "fixed_cli":
            _append_run_log(
                dataset_log,
                "INFO",
                "--fixed-info %g was not used for dataset %s because %s has priority."
                % (
                    fixed_info_from_cli,
                    row.dataset_id,
                    engine_input["info_source_detail"],
                ),
            )
        elif engine_input["info_source"] == "fixed_cli":
            message = (
                "Dataset %s assigns INFO=%g to every variant from the explicit "
                "--fixed-info option. This is a user-assigned constant, not a "
                "variant-level measured imputation-quality score."
                % (row.dataset_id, fixed_info_from_cli)
            )
            _append_run_log(dataset_log, "WARNING", message)
            print(screen_field(
                "warning", "Fixed INFO", message,
                indent=4, label_width=20,
            ))
        _append_run_log(dataset_log, "INFO", "Starting dataset %s." % row.dataset_id)
        attempts = config.execution.retries + 1
        for attempt in range(1, attempts + 1):
            try:
                result = run_harmonisation_pipeline(
                    sample_column_dict=engine_input,
                    default_cfg=defaults,
                    threads=config.execution.threads,
                )
                if config.modules.harmonisation.concordance_validation.enabled:
                    from postgwas.modules.harmonisation.concordance.service import run_concordance_validation

                    manifest_document = json.loads(
                        Path(result["manifest"]).read_text(encoding="utf-8")
                    )
                    dataset_context = manifest_document.get("dataset") or {}
                    build_context = dataset_context.get("genome_build") or {}
                    input_build = build_context.get("inferred_build")
                    supported_builds = set(
                        defaults["vcf_processing"]["target_builds"]
                    )
                    if input_build not in supported_builds:
                        raise ConcordanceValidationError(
                            "Harmonisation result did not record a resolved input genome build."
                        )
                    validation = run_concordance_validation(
                        row=row,
                        vcf_path=result[input_build],
                        output_root=output_root,
                        settings=config.modules.harmonisation.concordance_validation,
                        policies=load_policies(config.modules.harmonisation.policies),
                        threads=config.execution.threads,
                        bcftools=config.resources.executables.bcftools,
                        output_layout=defaults["output_layout"],
                        vcf_config=defaults["vcf_processing"],
                        file_log_level=config.logging.file_level,
                        screen_log_level=config.logging.console_level,
                        run_manifest=result["manifest"],
                        external_eaf_mapping=config.modules.harmonisation.external_eaf_mapping,
                    )
                    result["concordance_validation"] = validation
                    if validation["status"] == "FAIL":
                        raise ConcordanceValidationError(
                            "Concordance validation failed for dataset %s; see %s"
                            % (row.dataset_id, validation["reports"]["summary"])
                        )
                results[row.dataset_id] = result
                _append_run_log(dataset_log, "INFO", "Dataset completed successfully.")
                _append_run_log(
                    top_log,
                    "INFO",
                    "Dataset %s completed successfully." % row.dataset_id,
                )
                print(screen_line(
                    "success", "Dataset %s completed" % row.dataset_id, indent=4,
                ))
                break
            except KeyboardInterrupt:
                _append_run_log(dataset_log, "WARNING", "Dataset interrupted by the user.")
                _append_run_log(
                    top_log,
                    "WARNING",
                    "Run interrupted by the user while processing dataset %s."
                    % row.dataset_id,
                )
                raise
            except ConfigError as exc:
                _append_run_log(dataset_log, "ERROR", "Configuration failure: %s" % exc)
                _append_run_log(
                    top_log,
                    "ERROR",
                    "Dataset %s failed configuration validation; see %s."
                    % (row.dataset_id, dataset_log),
                )
                failures.append((row.dataset_id, exc))
                break
            except ConcordanceValidationError as exc:
                _append_run_log(
                    dataset_log,
                    "ERROR",
                    "Concordance validation failed after harmonisation: %s" % exc,
                )
                _append_run_log(
                    top_log,
                    "ERROR",
                    "Dataset %s harmonised but failed concordance validation; see %s."
                    % (row.dataset_id, dataset_log),
                )
                failures.append((row.dataset_id, exc))
                break
            except Exception as exc:
                _append_run_log(
                    dataset_log,
                    "ERROR",
                    "Attempt %d/%d failed: %s: %s"
                    % (attempt, attempts, type(exc).__name__, exc),
                )
                if attempt == attempts:
                    _append_run_log(
                        dataset_log,
                        "ERROR",
                        "Dataset failed permanently.\n%s" % traceback.format_exc(),
                    )
                    _append_run_log(
                        top_log,
                        "ERROR",
                        "Dataset %s failed after %d attempt(s); see %s."
                        % (row.dataset_id, attempts, dataset_log),
                    )
                    failures.append((row.dataset_id, exc))
    if failures and not results:
        _append_run_log(
            top_log,
            "ERROR",
            "Run failed: all %d dataset(s) failed." % len(failures),
        )
        raise PipelineError("All harmonisation datasets failed; see %s" % top_log)
    if failures:
        _append_run_log(
            top_log,
            "WARNING",
            "Run completed with %d successful and %d failed dataset(s)."
            % (len(results), len(failures)),
        )
        print(screen_field(
            "warning", "Datasets failed", "%d; see %s" % (len(failures), top_log),
            indent=4, label_width=18,
        ))
    else:
        _append_run_log(top_log, "INFO", "Run completed successfully: %d dataset(s)." % len(results))
    return results


def run_harmonisation(args):
    if not hasattr(args, "sample_sheet"):
        raise ConfigurationError("--sample-sheet is required")

    overrides = explicit_overrides(args, CLI_OVERRIDE_PATHS)
    config = load_configuration(getattr(args, "run_config", None), cli_overrides=overrides)
    preflight_log = _top_run_log(config, config.run.output_directory)
    try:
        rows = load_harmonisation_sample_sheet(args.sample_sheet)
        if hasattr(args, "dataset_id"):
            rows = [row for row in rows if row.dataset_id == args.dataset_id]
            if not rows:
                raise ConfigurationError(
                    "dataset_id %r is not present in sample sheet %s"
                    % (args.dataset_id, args.sample_sheet)
                )
        _validate_info_fallback(rows, _fixed_info_cli_value(args, config))
        resolved_executables = require_binaries(
            _harmonisation_executables(config),
            plugins=(config.modules.harmonisation.vcf_processing.liftover_plugin,),
        )
    except Exception as exc:
        _append_run_log(preflight_log, "ERROR", "Harmonisation preflight failed: %s" % exc)
        raise
    return _run_validated_rows(
        args, config, rows, resolved_executables=resolved_executables,
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        prog="postgwas harmonisation",
        usage="postgwas harmonisation --sample-sheet PATH [options]",
        description=(
            "Prepare one or more GWAS summary-statistics files for PostGWAS. "
            "The command checks the sample sheet, harmonises alleles and genome "
            "coordinates, and writes analysis-ready GWAS-VCF files."
        ),
        epilog=format_cli_examples(
            (
                "Harmonise every dataset in a sample sheet:",
                "postgwas harmonisation",
                (
                    "--sample-sheet studies.csv",
                    "--run-config harmonisation.yaml",
                    "--resource-directory resources",
                    "--output-directory results",
                ),
            ),
            (
                "Harmonise one selected dataset and validate its merged VCF:",
                "postgwas harmonisation",
                (
                    "--sample-sheet studies.csv",
                    "--dataset-id STUDY",
                    "--resource-directory resources",
                    "--output-directory results",
                    "--validate",
                ),
            ),
            (
                "Export a reusable harmonisation configuration:",
                "postgwas config export",
                (
                    "--module harmonisation",
                    "--style full",
                    "--output harmonisation.yaml",
                ),
            ),
            notes=(
                "Sample-sheet templates are available in examples/configs/harmonisation/.",
            ),
        ),
        parents=[get_harmonisation_parser()],
        formatter_class=AlignedRichHelpFormatter,
    )
    if len(sys.argv) == 1:
        parser.print_help()
        return 0
    args = parser.parse_args()
    try:
        result = run_harmonisation(args)
    except KeyboardInterrupt:
        try:
            fallback = load_configuration()
            output = Path(getattr(args, "output_directory", fallback.run.output_directory))
            interruption_log = _top_run_log(fallback, output)
            _append_run_log(
                interruption_log,
                "WARNING",
                "Harmonisation interrupted by the user; the run did not complete.",
            )
        except Exception:
            interruption_log = None
        suffix = "; see %s" % interruption_log if interruption_log else ""
        print(screen_field(
            "warning", "Interrupted", "harmonisation stopped%s" % suffix,
            indent=4, label_width=18,
        ), file=sys.stderr)
        return 130
    except (ConfigurationError, ConfigError) as exc:
        try:
            fallback = load_configuration()
            output = Path(getattr(args, "output_directory", fallback.run.output_directory))
            _append_run_log(
                _top_run_log(fallback, output),
                "ERROR",
                "Configuration error: %s" % exc,
            )
        except Exception:
            pass
        print(screen_field(
            "error", "Configuration", exc,
            indent=4, label_width=18,
        ), file=sys.stderr)
        return 2
    except Exception as exc:
        try:
            fallback = load_configuration()
            output = Path(getattr(args, "output_directory", fallback.run.output_directory))
            failure_log = _top_run_log(fallback, output)
            _append_run_log(
                failure_log,
                "ERROR",
                "Harmonisation failed: %s: %s\n%s"
                % (type(exc).__name__, exc, traceback.format_exc()),
            )
        except Exception:
            failure_log = None
        suffix = (
            "; see %s" % failure_log
            if failure_log and str(failure_log) not in str(exc)
            else ""
        )
        print(screen_field(
            "error", "Harmonisation", "%s%s" % (exc, suffix),
            indent=4, label_width=18,
        ), file=sys.stderr)
        return 1
    print(screen_line(
        "success", "Harmonisation complete — %d dataset(s)" % len(result), indent=4,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
