"""Standalone ``postgwas --validate`` command."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import traceback

from postgwas.cli.compute import get_compute_parser
from postgwas.config import load_configuration
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.core.errors import ConfigurationError
from postgwas.core.ui import (
    AlignedRichHelpFormatter,
    format_cli_examples,
    help_with_default,
    screen_field,
    screen_line,
)
from postgwas.modules.harmonisation.sample_sheet import load_harmonisation_sample_sheet
from postgwas.modules.harmonisation.policies import load_policies

from .errors import ConcordanceValidationError
from .service import run_concordance_validation, write_validation_preflight_failure


CLI_OVERRIDE_PATHS = {
    "output_directory": "run.output_directory",
    "threads": "execution.threads",
    "memory_gb": "execution.memory_gb",
    "seed": "execution.random_seed",
}


def get_validation_parser() -> argparse.ArgumentParser:
    defaults = load_configuration()
    parser = argparse.ArgumentParser(
        prog="postgwas --validate",
        usage="postgwas --validate --sample-sheet PATH --dataset-id ID --vcf PATH [options]",
        description=(
            "Compare one original summary-statistics dataset with its corresponding "
            "same-build harmonised GWAS-VCF.\n\n"
            "PostGWAS checks variant retention, effect estimates, allele frequencies, "
            "standard errors, Z scores, and p-values. If the input has no Z column, "
            "it calculates Z from the available effect statistics. SNP and indel "
            "concordance are reported "
            "separately. Unmatched variants are reported without failing the value audit; "
            "one-to-one unmatched records at the same position receive diagnostic value checks."
        ),
        epilog=format_cli_examples(
            (
                "Validate one harmonised dataset:",
                "postgwas --validate",
                (
                    "--sample-sheet studies.csv",
                    "--dataset-id STUDY",
                    "--vcf STUDY_GRCh37_merged.vcf.gz",
                    "--output-directory results",
                ),
            ),
            notes=(
                "Detailed TSV reports and the complete log are saved under the dataset harmonisation output.",
                "Validation runs only when --validate is explicitly requested.",
            ),
        ),
        parents=[get_compute_parser()],
        formatter_class=AlignedRichHelpFormatter,
    )
    files = parser.add_argument_group("Files and dataset")
    files.add_argument(
        "--sample-sheet", required=True, metavar="PATH",
        help="CSV or TSV sample sheet used for harmonisation.",
    )
    files.add_argument(
        "--dataset-id", required=True, metavar="ID",
        help="Dataset row to validate; matching is exact and case-sensitive.",
    )
    files.add_argument(
        "--vcf", required=True, metavar="PATH",
        help="Unfiltered merged GWAS-VCF in the input summary statistics genome build.",
    )
    files.add_argument(
        "--run-config", default=argparse.SUPPRESS, metavar="PATH",
        help=(
            "Optional YAML configuration controlling concordance tolerances, strand "
            "matching, palindromic variants, and failure thresholds."
        ),
    )
    files.add_argument(
        "--output-directory", default=argparse.SUPPRESS, metavar="PATH",
        help=help_with_default(
            "Root output folder", defaults.run.output_directory,
        ),
    )
    return parser


def run_standalone_validation(args) -> dict:
    overrides = explicit_overrides(args, CLI_OVERRIDE_PATHS)
    overrides["modules.harmonisation.concordance_validation.enabled"] = True
    config = load_configuration(getattr(args, "run_config", None), cli_overrides=overrides)
    args._resolved_output_directory = config.run.output_directory
    args._resolved_output_layout = dict(
        config.modules.harmonisation.output_layout.root
    )
    args._resolved_logging = (
        config.logging.file_level,
        config.logging.console_level,
    )
    rows = load_harmonisation_sample_sheet(args.sample_sheet)
    selected = [row for row in rows if row.dataset_id == args.dataset_id]
    if not selected:
        raise ConfigurationError(
            "dataset_id %r is not present in sample sheet %s"
            % (args.dataset_id, args.sample_sheet)
        )
    if len(selected) != 1:
        raise ConfigurationError(
            "dataset_id %r occurs more than once in the sample sheet" % args.dataset_id
        )
    return run_concordance_validation(
        row=selected[0],
        vcf_path=args.vcf,
        output_root=config.run.output_directory,
        settings=config.modules.harmonisation.concordance_validation,
        policies=load_policies(config.modules.harmonisation.policies),
        threads=config.execution.threads,
        bcftools=config.resources.executables.bcftools,
        output_layout=args._resolved_output_layout,
        vcf_config=config.modules.harmonisation.vcf_processing.model_dump(),
        file_log_level=config.logging.file_level,
        screen_log_level=config.logging.console_level,
        external_eaf_mapping=config.modules.harmonisation.external_eaf_mapping,
    )


def _failure_output_layout(args) -> dict[str, str]:
    configured = getattr(args, "_resolved_output_layout", None)
    if configured is not None:
        return configured
    return dict(load_configuration().modules.harmonisation.output_layout.root)


def _failure_logging(args) -> tuple[str, str]:
    configured = getattr(args, "_resolved_logging", None)
    if configured is not None:
        return configured
    logging = load_configuration().logging
    return logging.file_level, logging.console_level


def main() -> int:
    parser = get_validation_parser()
    if len(sys.argv) == 1:
        parser.print_help()
        return 0
    args = parser.parse_args()
    try:
        result = run_standalone_validation(args)
    except KeyboardInterrupt:
        file_level, screen_level = _failure_logging(args)
        write_validation_preflight_failure(
            getattr(args, "_resolved_output_directory", getattr(args, "output_directory", "results")),
            getattr(args, "dataset_id", "unknown_dataset"),
            "Validation interrupted by the user.",
            _failure_output_layout(args),
            file_level=file_level,
            screen_level=screen_level,
        )
        print(screen_field(
            "warning", "Interrupted", "concordance validation stopped",
            indent=4, label_width=18,
        ), file=sys.stderr)
        return 130
    except (ConfigurationError, ConcordanceValidationError) as exc:
        file_level, screen_level = _failure_logging(args)
        log_path = write_validation_preflight_failure(
            getattr(args, "_resolved_output_directory", getattr(args, "output_directory", "results")),
            getattr(args, "dataset_id", "unknown_dataset"),
            "%s: %s" % (type(exc).__name__, exc),
            _failure_output_layout(args),
            file_level=file_level,
            screen_level=screen_level,
        )
        print(screen_field(
            "error", "Validation", "%s; log: %s" % (exc, log_path),
            indent=4, label_width=18,
        ), file=sys.stderr)
        return 2 if isinstance(exc, ConfigurationError) else 1
    except Exception as exc:
        file_level, screen_level = _failure_logging(args)
        log_path = write_validation_preflight_failure(
            getattr(args, "_resolved_output_directory", getattr(args, "output_directory", "results")),
            getattr(args, "dataset_id", "unknown_dataset"),
            "%s: %s\n%s" % (type(exc).__name__, exc, traceback.format_exc()),
            _failure_output_layout(args),
            file_level=file_level,
            screen_level=screen_level,
        )
        print(screen_field(
            "error", "Validation", "%s; log: %s" % (exc, log_path),
            indent=4, label_width=18,
        ), file=sys.stderr)
        return 1
    if result["status"] == "FAIL":
        print(screen_line("error", "Concordance validation failed", indent=4))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
