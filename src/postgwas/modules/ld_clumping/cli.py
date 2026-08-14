#!/usr/bin/env python3
"""Command-line interface for configured LD clumping."""

from __future__ import annotations

import argparse
import sys

from postgwas.cli.common import (
    get_bcftools_binary_parser,
    get_common_out_parser,
    get_inputvcf_parser,
    get_ld_clump_parser,
    get_ld_clumping_genome_build_parser,
    get_ld_clumping_population_parser,
    get_tabix_binary_parser,
)
from postgwas.cli.compute import get_compute_parser
from postgwas.core.errors import ConfigurationError
from postgwas.core.ui import (
    AlignedRichHelpFormatter,
    format_cli_examples,
    mark_cli_required_help,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="postgwas ld_clump",
        usage="postgwas ld_clump --vcf PATH [--ld-folder PATH] [options]",
        description=(
            "Identify association signals using annotated LD blocks, a symmetric "
            "population LD reference, or both configured analyses."
        ),
        epilog=format_cli_examples(
            (
                "Run both configured analyses (the packaged default):",
                "postgwas ld_clump",
                (
                    "--vcf study_ldblock.vcf.gz",
                    "--ld-folder reference/ld",
                    "--genome-build GRCh37",
                    "--population EUR",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Run only annotated-region pruning:",
                "postgwas ld_clump",
                (
                    "--vcf study_ldblock.vcf.gz",
                    "--clumping-methods region",
                    "--genome-build GRCh37",
                    "--population EUR",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Run only standard r² clumping:",
                "postgwas ld_clump",
                (
                    "--vcf study.vcf.gz",
                    "--clumping-methods standard",
                    "--ld-folder reference/ld",
                    "--genome-build GRCh37",
                    "--population EUR",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Export a reusable LD-clumping configuration:",
                "postgwas config export",
                ("--module ld_clumping", "--style full", "--output ld_clumping.yaml"),
            ),
        ),
        formatter_class=AlignedRichHelpFormatter,
        parents=[
            get_compute_parser(),
            get_bcftools_binary_parser(),
            get_tabix_binary_parser(),
            get_inputvcf_parser(),
            get_common_out_parser(),
            get_ld_clumping_genome_build_parser(),
            get_ld_clumping_population_parser(),
            get_ld_clump_parser(add_help=False),
        ],
    )
    parser.add_argument(
        "--run-config",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "YAML file containing LD-clumping settings or a complete PostGWAS "
            "run configuration."
        ),
    )
    mark_cli_required_help(
        parser, ("vcf", "dataset_id", "output_directory"),
    )
    return parser


def main() -> int:
    parser = build_parser()
    if len(sys.argv) == 1:
        parser.print_help()
        return 0
    args = parser.parse_args()
    from postgwas.modules.ld_clumping.service import (
        resolve_ld_clumping_configuration,
        run_ld_clump_direct,
    )

    try:
        configuration = resolve_ld_clumping_configuration(args)
        run_ld_clump_direct(args, configuration=configuration)
    except ConfigurationError as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    main()


__all__ = ["build_parser", "main"]
