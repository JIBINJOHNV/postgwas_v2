#!/usr/bin/env python3
"""
PostGWAS Fine-Mapping CLI
Supports two engines:
  • SuSiE  (Python/R-based)
  • FINEMAP (C++ binary)
Works in:
  • DIRECT mode   → run fine-mapping on existing locus + sumstats
  • PIPELINE mode → Harmonisation → LD Blocks → QC → Formatter → Fine-mapping
"""

import argparse
import sys
from postgwas.cli.compute import get_compute_parser
from postgwas.core.ui import AlignedRichHelpFormatter, format_cli_examples


# =========================================================
# BACKEND RUNNERS
# =========================================================
from postgwas.modules.fine_mapping.service import run_fine_mapping

# Shared parsers
from postgwas.cli.common import (
    get_common_out_parser,
    get_plink_binary_parser,
)

from postgwas.modules.fine_mapping.arguments import (
    get_finemap_common_parser,
    get_common_susie_arguments,
    get_common_finemap_finemap_arguments,
    get_finemap_susie_inputs_parser,
    get_finemap_finemap_arguments,
)


def get_finemap_pipeline_examples():
    """Return a dependency-complete SuSiE fine-mapping pipeline example."""
    return (
        (
            "Run LD annotation, clumping, formatting, and SuSiE fine-mapping:",
            "postgwas pipeline",
            (
                "--modules finemap",
                "--vcf study_GRCh37.vcf.gz",
                "--genome-build GRCh37",
                "--ld-region-dir reference/ld_blocks",
                "--ld-block-populations EUR",
                "--ld-folder reference/pairwise_ld",
                "--population EUR",
                "--variant-id-type unique",
                "--finemap-method susie",
                "--finemap-ld-reference reference/1000G_EUR",
                "--plink plink",
                "--bcftools bcftools",
                "--dataset-id STUDY",
                "--output-directory results",
            ),
        ),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="postgwas finemap",
        usage="postgwas finemap --finemap-method {susie,finemap} [options]",
        description=(
            "Fine-map selected loci with SuSiE-RSS or FINEMAP.\n"
            "Standalone mode uses formatter-created summary statistics and an existing locus file."
        ),
        epilog=format_cli_examples(
            (
                "Run SuSiE-RSS fine-mapping:",
                "postgwas finemap",
                (
                    "--finemap-method susie",
                    "--susie-input-file formatted/STUDY_susie.tsv.gz",
                    "--locus-file loci.tsv",
                    "--finemap-ld-reference reference/1000G_EUR",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Run FINEMAP:",
                "postgwas finemap",
                (
                    "--finemap-method finemap",
                    "--finemap-in-files formatted/STUDY_finemap.tsv.gz",
                    "--locus-file loci.tsv",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Export the fine-mapping pipeline configuration:",
                "postgwas config export",
                ("--pipeline finemap", "--style full", "--output finemap_pipeline.yaml"),
            ),
        ),
        formatter_class=AlignedRichHelpFormatter,
        parents=[
            get_compute_parser(),
            get_common_out_parser(),
            get_finemap_common_parser(include_genome_build=True),
            get_finemap_susie_inputs_parser(),
            get_common_susie_arguments(),
            get_common_finemap_finemap_arguments(),
            get_finemap_finemap_arguments(),
            get_plink_binary_parser(),
        ],
    )
    configuration = parser.add_argument_group("Configuration")
    configuration.add_argument(
        "--run-config",
        default=argparse.SUPPRESS,
        metavar="PATH",
        help=(
            "YAML file containing fine-mapping settings. Explicit command-line "
            "values override matching YAML values."
        ),
    )
    return parser


def main():
    parser = build_parser()

    # If no arguments provided → show help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    run_fine_mapping(args)


if __name__ == "__main__":
    main()


__all__ = ["build_parser", "get_finemap_pipeline_examples", "main"]
