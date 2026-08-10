#!/usr/bin/env python3
"""
PostGWAS — LDSC Heritability / Genetic Correlation CLI

Mode:
  • Direct – Run LDSC on existing munged .sumstats.gz
"""
import argparse
import sys
from postgwas.cli.compute import get_compute_parser
from postgwas.core.ui import AlignedRichHelpFormatter, format_cli_examples

# =========================================================
# BACKEND RUNNERS
# =========================================================
from postgwas.modules.ldsc.service import run_ldsc_direct


# Shared parsers
from postgwas.cli.common import (
    get_common_out_parser,
    get_ldsc_common_parser
)

from postgwas.core.execution.runtime import validate_path


# =========================================================
# SHARED LDSC INPUT ARGUMENTS
# =========================================================
def get_ldsc_parser(add_help: bool = False) -> argparse.ArgumentParser:
    """
    Defines input arguments for LDSC direct execution.
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    grp = parser.add_argument_group("Formatted summary statistics")

    grp.add_argument(
        "--ldsc-input",
        metavar="PATH",
        type=validate_path(must_exist=True, must_be_file=True),
        help="REQUIRED. LDSC input TSV created by the PostGWAS formatter.",
    )

    return parser


# =========================================================
# MAIN CLI
# =========================================================
def build_parser() -> argparse.ArgumentParser:
    return argparse.ArgumentParser(
        prog="postgwas heritability",
        usage="postgwas heritability --ldsc-input PATH [options]",
        description="Estimate single-trait SNP heritability from formatter-created LDSC input.",
        epilog=format_cli_examples(
            (
                "Estimate observed-scale SNP heritability:",
                "postgwas heritability",
                (
                    "--ldsc-input formatted/STUDY_ldsc.tsv.gz",
                    "--merge-alleles reference/w_hm3.snplist",
                    "--ref-ld-chr reference/eur_w_ld_chr",
                    "--w-ld-chr reference/eur_w_ld_chr",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Add prevalence values for liability-scale reporting:",
                "postgwas heritability",
                (
                    "--ldsc-input formatted/STUDY_ldsc.tsv.gz",
                    "--ref-ld-chr reference/eur_w_ld_chr",
                    "--w-ld-chr reference/eur_w_ld_chr",
                    "--samp-prev 0.2",
                    "--pop-prev 0.01",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
        ),
        parents=[
            get_compute_parser(),
            get_common_out_parser(),
            get_ldsc_parser(),
            get_ldsc_common_parser(add_help=False),
        ],
        formatter_class=AlignedRichHelpFormatter,
    )


def main():
    parser = build_parser()

    # ---------------------------------------------------
    # EXECUTION
    # ---------------------------------------------------

    # Check if no arguments provided (sys.argv[0] is the script name)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    # Directly dispatch to the direct runner
    run_ldsc_direct(args)


if __name__ == "__main__":
    main()


__all__ = ["build_parser", "get_ldsc_parser", "main"]
