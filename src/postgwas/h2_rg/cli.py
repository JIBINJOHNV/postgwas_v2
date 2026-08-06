#!/usr/bin/env python3
"""
PostGWAS — LDSC Heritability / Genetic Correlation CLI

Mode:
  • Direct – Run LDSC on existing munged .sumstats.gz
"""
import argparse
import sys
from rich_argparse import RichHelpFormatter

# =========================================================
# BACKEND RUNNERS
# =========================================================
from postgwas.h2_rg.workflows import run_ldsc_direct


# Shared parsers
from postgwas.clis.common_cli import (
    get_defaultresourse_parser,
    get_common_out_parser,
    get_ldsc_common_parser
)

from postgwas.utils.main import validate_path


# =========================================================
# SHARED LDSC INPUT ARGUMENTS
# =========================================================
def get_ldsc_parser(add_help: bool = False) -> argparse.ArgumentParser:
    """
    Defines input arguments for LDSC direct execution.
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    grp = parser.add_argument_group("LDSC Inputs")

    grp.add_argument(
        "--ldsc_inut",
        metavar=" ",
        type=validate_path(must_exist=True, must_be_file=True),
        help="Path to the LDSC input tsv file genrated by  postgwas fomatter module",
    )

    return parser


# =========================================================
# MAIN CLI
# =========================================================
def main():
    # If the user runs the script with no arguments, print help and exit
    if len(sys.argv) == 1:
        # We temporarily create a parser just to print the help, 
        # or we can rely on the main parser construction below.
        # To ensure it prints help cleanly before parsing:
        pass 

    parser = argparse.ArgumentParser(
        prog="postgwas-ldsc",
        description="LDSC heritability module.",
        parents=[
            get_defaultresourse_parser(),
            get_common_out_parser(),
            get_ldsc_parser(),
            get_ldsc_common_parser(add_help=False),
        ],
        formatter_class=RichHelpFormatter,
    )

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