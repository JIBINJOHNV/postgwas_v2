"""Command-line interface for MAGMA gene-property analysis."""

from __future__ import annotations

import argparse

from postgwas.cli.common import (
    get_common_magma_covar_parser,
    get_common_out_parser,
)
from postgwas.cli.compute import get_compute_parser
from postgwas.core.errors import ConfigurationError
from postgwas.core.ui import (
    print_screen_message,
    AlignedRichHelpFormatter,
    format_cli_examples,
    mark_cli_required_help,
)
from postgwas.modules.magmacovar.contract import (
    MAGMACOVAR_EXCLUDED_PATHWAY_OPTIONS,
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
            "MAGMA .genes.raw file produced by a completed gene-association run."
        ),
    )
    return parser


def _get_direct_controls_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group("Configuration and continuation")
    group.add_argument(
        "--run-config",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help="YAML settings; explicit command-line values override matching keys.",
    )
    return parser


def get_magmacovar_pipeline_examples():
    """Return complete examples for context-sensitive pipeline help."""
    common_arguments = (
        "--vcf study.vcf.gz",
        "--magma-ld-reference reference/1000G_EUR",
        "--gene-location-file reference/NCBI37.3.gene.loc",
    )
    return (
        (
            "Run marginal gene-property tests, as in the original FLAMES README:",
            "postgwas pipeline",
            (
                "--modules magmacovar",
                *common_arguments,
                "--covariates reference/GTEx/gtex_v8_ts_avg_log2TPM.txt",
                "--dataset-id STUDY",
                "--output-directory results",
            ),
        ),
        (
            "Run FUMA-style tissue-specific gene-property tests:",
            "postgwas pipeline",
            (
                "--modules magmacovar",
                *common_arguments,
                "--covariates "
                "reference/GTEx/gtex_v8_ts_general_avg_log2TPM.txt",
                "--covariate-model condition-hide=Average",
                "--covariate-direction greater",
                "--dataset-id STUDY",
                "--output-directory results",
            ),
        ),
    )


def get_magmacovar_pipeline_parser(add_help=False):
    """Return concise MAGMAcovar options for context-sensitive pipeline help."""
    return get_common_magma_covar_parser(
        add_help=add_help,
        compact_help=True,
    )


def customize_magmacovar_only_pipeline_help(
    parser: argparse.ArgumentParser,
) -> None:
    """Remove pathway-only dependency controls from MAGMAcovar-only help."""
    excluded_destinations = {
        destination
        for destination, _flag in MAGMACOVAR_EXCLUDED_PATHWAY_OPTIONS
    }
    for action in parser._actions:
        if action.dest in excluded_destinations:
            action.required = False
            action.help = argparse.SUPPRESS
        elif action.dest == "mhc_policy" and isinstance(action.help, str):
            action.help = action.help.replace(
                "removes overlapping annotated units from gene and competitive "
                "gene-set analysis",
                "removes overlapping annotated genes before downstream "
                "gene-property analysis",
            )
        elif action.dest == "exclude_chromosomes" and isinstance(action.help, str):
            action.help = action.help.replace(
                "Chromosomes excluded from SNP, gene, and competitive gene-set "
                "analysis",
                "Chromosomes excluded from SNP and gene analysis",
            )
    parser._postgwas_pipeline_step_descriptions = {
        "magma": (
            "Run the MAGMA gene-association analysis required by MAGMAcovar."
        ),
    }


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
    return mark_cli_required_help(
        parser,
        (
            "magma_gene_results_file",
            "covariates",
            "dataset_id",
            "output_directory",
        ),
    )


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    from postgwas.modules.magmacovar.service import run_magma_covar_direct

    try:
        run_magma_covar_direct(args)
    except (MagmaCovarError, ConfigurationError, OSError, ValueError) as exc:
        print_screen_message(
            "error", "MAGMA gene-property analysis failed. %s" % exc, stderr=True,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "build_parser",
    "customize_magmacovar_only_pipeline_help",
    "get_magma_covar_gene_results_parser",
    "get_magmacovar_pipeline_parser",
    "get_magmacovar_pipeline_examples",
    "main",
]
