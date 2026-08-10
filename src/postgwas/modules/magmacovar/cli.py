"""Command-line interface for MAGMA gene-property analysis."""

from __future__ import annotations

import argparse

from rich.console import Console

from postgwas.cli.common import (
    get_common_magma_covar_parser,
    get_common_out_parser,
    get_magma_binary_parser,
)
from postgwas.cli.compute import get_compute_parser
from postgwas.config import load_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.core.ui import (
    AlignedRichHelpFormatter,
    format_cli_examples,
    help_with_default,
)
from postgwas.modules.magmacovar.errors import MagmaCovarError


def get_magma_covar_gene_results_parser(add_help=False):
    """Add the MAGMA gene-level results consumed only in direct mode."""
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("MAGMA gene-property input")
    group.add_argument(
        "--magma-gene-results-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "MAGMA .genes.raw file produced by a completed gene-association run. "
            "Pipeline mode supplies this automatically."
        ),
    )
    return parser


def _get_direct_controls_parser() -> argparse.ArgumentParser:
    defaults = load_configuration()
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group("Configuration and continuation")
    group.add_argument(
        "--run-config",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help="YAML settings; explicit command-line values override matching keys.",
    )
    group.add_argument(
        "--resume",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Reuse a completed MAGMAcovar result only after input, configuration, "
            "and output fingerprints match",
            defaults.run.resume,
        ),
    )
    group.add_argument(
        "--overwrite",
        action="store_true",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Replace output owned by this configured MAGMAcovar prefix",
            defaults.run.overwrite,
        ),
    )
    return parser


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="postgwas magmacovar",
        usage=(
            "postgwas magmacovar --magma-gene-results-file PATH "
            "--covariates PATH [options]"
        ),
        description=(
            "Run validated MAGMA gene-property analysis using existing MAGMA gene "
            "results and a gene-level covariate table."
        ),
        epilog=format_cli_examples(
            (
                "Run the documented two-sided marginal gene-property tests:",
                "postgwas magmacovar",
                (
                    "--magma-gene-results-file STUDY.genes.raw",
                    "--covariates gene_covariates.tsv",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Create tissue-relevance input as documented by original FLAMES:",
                "postgwas magmacovar",
                (
                    "--magma-gene-results-file STUDY.genes.raw",
                    "--covariates "
                    "reference/GTEx/gtex_v8_ts_avg_log2TPM.txt",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Run the separate FUMA-compatible tissue-specific model:",
                "postgwas magmacovar",
                (
                    "--magma-gene-results-file STUDY.genes.raw",
                    "--covariates "
                    "reference/GTEx/gtex_v8_ts_general_avg_log2TPM.txt",
                    "--covariate-model condition-hide=Average",
                    "--covariate-direction greater",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Export reusable MAGMAcovar settings:",
                "postgwas config export",
                (
                    "--module magmacovar",
                    "--style full",
                    "--output magmacovar.yaml",
                ),
            ),
            (
                "Run entirely from exported settings:",
                "postgwas magmacovar",
                ("--run-config magmacovar.yaml",),
            ),
            notes=(
                "The original FLAMES README uses the 54-specific-tissue GTEx v8 "
                "file and omits --model, which means marginal two-sided tests.",
                "FUMA's FLAMES integration instead consumes its 30-general-tissue "
                "result generated with condition-hide=Average and greater.",
                "The bundled FLAMES example .gsa.out records the FUMA-style model, "
                "so it does not match the literal command printed in its README.",
            ),
        ),
        formatter_class=AlignedRichHelpFormatter,
        parents=[
            get_compute_parser(),
            get_common_out_parser(),
            get_magma_binary_parser(),
            get_common_magma_covar_parser(add_help=False),
            get_magma_covar_gene_results_parser(add_help=False),
            _get_direct_controls_parser(),
        ],
    )
    for action in parser._actions:
        if action.dest in {
            "dataset_id", "output_directory", "threads", "memory_gb", "seed",
        }:
            action.default = argparse.SUPPRESS
            action.type = None
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    from postgwas.modules.magmacovar.service import run_magma_covar_direct

    try:
        run_magma_covar_direct(args)
    except (MagmaCovarError, ConfigurationError, OSError, ValueError) as exc:
        Console(stderr=True).print(
            "\n[bold red]MAGMA gene-property analysis failed.[/bold red] %s\n" % exc
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = ["build_parser", "get_magma_covar_gene_results_parser", "main"]
