#!/usr/bin/env python3
"""Command-line interface for the upstream FLAMES integration."""

from __future__ import annotations

import argparse

from rich.console import Console

from postgwas.cli.common import get_common_out_parser, get_flames_common_parser
from postgwas.cli.compute import get_compute_parser
from postgwas.config import load_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.core.ui import (
    AlignedRichHelpFormatter,
    format_cli_examples,
    help_with_default,
)
from postgwas.modules.flames.errors import FlamesError


def get_flames_pipeline_examples():
    """Return a dependency-complete GRCh37 pipeline example for FLAMES."""
    return (
        (
            "Run the complete GRCh37 fine-mapping-to-FLAMES workflow:",
            "postgwas pipeline",
            (
                "--modules flames",
                "--vcf study_GRCh37.vcf.gz",
                "--genome-build GRCh37",
                "--ld-region-dir reference/ld_blocks",
                "--ld-block-populations EUR",
                "--ld-folder reference/pairwise_ld",
                "--population EUR",
                "--magma-ld-reference reference/1000G_EUR",
                "--gene-location-file reference/FUMA/ENSGv102.coding.genes.txt",
                "--covariates "
                "reference/GTEx/gtex_v8_ts_avg_log2TPM.txt",
                "--feature-matrix-prefix "
                "reference/FLAMES/pops_features_full_FUMA_compatible/"
                "features_munged/pops_features",
                "--feature-matrix-chunks 116",
                "--pops-gene-location-file "
                "reference/FLAMES/pops_features_full_FUMA_compatible/"
                "gene_annots.txt",
                "--control-features-file "
                "reference/FLAMES/pops_features_full_FUMA_compatible/"
                "control.features",
                "--pops-genome-build GRCh37",
                "--finemap-method susie",
                "--finemap-ld-reference reference/1000G_EUR",
                "--flames-annotation-directory "
                "reference/FLAMES/Annotation_data",
                "--flames-genome-build GRCh37",
                "--plink plink",
                "--bcftools bcftools",
                "--dataset-id STUDY",
                "--output-directory results",
            ),
        ),
    )


def get_flames_direct_parser(add_help=False):
    """Return direct-mode inputs and run controls without argparse defaults."""
    defaults = load_configuration()
    parser = argparse.ArgumentParser(add_help=add_help)
    inputs = parser.add_argument_group("FLAMES inputs")
    inputs.add_argument(
        "--credible-sets-directory", metavar="PATH", default=argparse.SUPPRESS,
        help=(
            "PostGWAS fine-mapping interchange directory containing the configured "
            "index and credible-set tables."
        ),
    )
    inputs.add_argument(
        "--magma-gene-results-file", metavar="PATH", default=argparse.SUPPRESS,
        help="MAGMA gene-association .genes.out file using compatible Ensembl IDs.",
    )
    inputs.add_argument(
        "--magma-covariate-results-file", metavar="PATH", default=argparse.SUPPRESS,
        help=(
            "MAGMA gene-property .gsa.out used as FLAMES tissue-relevance weights. "
            "The original FLAMES README generates it from the 54-specific-tissue "
            "GTEx v8 table using marginal two-sided MAGMA defaults."
        ),
    )
    inputs.add_argument(
        "--pops-scores-file", metavar="PATH", default=argparse.SUPPRESS,
        help="PoPS .preds file using the same Ensembl gene annotation as MAGMA.",
    )
    controls = parser.add_argument_group("Configuration")
    controls.add_argument(
        "--run-config", metavar="PATH", default=argparse.SUPPRESS,
        help="YAML settings; explicit command-line values override matching keys.",
    )
    controls.add_argument(
        "--resume", action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Reuse outputs only after their completion manifest and inputs validate",
            defaults.run.resume,
        ),
    )
    controls.add_argument(
        "--overwrite", action="store_true", default=argparse.SUPPRESS,
        help="Replace an existing complete FLAMES result after validation.",
    )
    controls.add_argument(
        "--dry-run", action="store_true", default=argparse.SUPPRESS,
        help=(
            "Validate configuration, inputs, resources, and commands without "
            "executing FLAMES."
        ),
    )
    return parser


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="postgwas flames",
        usage="postgwas flames [--run-config PATH] [overrides]",
        description=(
            "Prioritise effector genes with upstream FLAMES using validated "
            "fine-mapping, MAGMA, gene-property, PoPS, and annotation resources."
        ),
        formatter_class=AlignedRichHelpFormatter,
        epilog=format_cli_examples(
            (
                "Run FLAMES from validated existing analysis outputs:",
                "postgwas flames",
                (
                    "--credible-sets-directory "
                    "finemap/downstream_inputs/flames",
                    "--magma-gene-results-file magma/STUDY.genes.out",
                    "--magma-covariate-results-file magma/STUDY.gsa.out",
                    "--pops-scores-file pops/STUDY.preds",
                    "--flames-annotation-directory reference/FLAMES/Annotation_data",
                    "--flames-genome-build GRCh37",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Run FLAMES from a complete exported YAML configuration:",
                "postgwas flames", ("--run-config flames.yaml",),
            ),
            (
                "Export reusable FLAMES settings:", "postgwas config export",
                ("--module flames", "--style full", "--output flames.yaml"),
            ),
            notes=(
                "The bundled published model is calibrated for a 750 kb gene window.",
                "Genome build and Ensembl gene annotation compatibility are "
                "declared, never inferred.",
                "For reproducible production annotation, configure versioned local "
                "VEP and CADD resources instead of live APIs.",
                "Pipeline defaults follow the literal original FLAMES README: the "
                "54-specific-tissue file with marginal two-sided MAGMA tests. "
                "Direct mode cannot infer those settings from an existing .gsa.out "
                "file, so verify its MAGMA log or provenance.",
                "FUMA uses a different 30-general-tissue input with "
                "condition-hide=Average and direction greater; pipeline users can "
                "select that model explicitly with the MAGMACOVAR options.",
            ),
        ),
        parents=[
            get_compute_parser(), get_common_out_parser(),
            get_flames_direct_parser(), get_flames_common_parser(),
        ],
    )
    for action in parser._actions:
        if action.dest in {
            "dataset_id", "output_directory", "threads", "memory_gb", "seed",
        }:
            action.default = argparse.SUPPRESS
            action.type = None
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    from postgwas.modules.flames.service import run_flames_direct

    try:
        run_flames_direct(args)
    except (FlamesError, ConfigurationError, OSError, ValueError) as exc:
        Console(stderr=True).print(
            "\n[bold red]FLAMES analysis failed.[/bold red] %s\n" % exc
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "build_parser",
    "get_flames_direct_parser",
    "get_flames_pipeline_examples",
    "main",
]
