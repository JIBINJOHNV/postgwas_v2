#!/usr/bin/env python3
"""
PostGWAS — Formatter (Step 4)

Convert harmonised GWAS VCF files into formats required for:
    - MAGMA
    - FINEMAP
    - LDpred
    - LDSC

This version is DIRECT-ONLY.
No pipeline mode and NO subcommands.
"""

import argparse
import sys
from rich_argparse import RichHelpFormatter

from postgwas.clis.common_cli import (
    get_defaultresourse_parser,
    get_inputvcf_parser,
    get_common_out_parser,
    get_formatter_parser,
    get_bcftools_binary_parser
)

from postgwas.formatter.workflows import run_formatter_direct


# ------------------------------------------------------------
# MAIN CLI (DIRECT ONLY)
# ------------------------------------------------------------
def main(argv=None):

    parser = argparse.ArgumentParser(
        prog="formatter",
        description="Convert harmonised GWAS VCF → MAGMA / FINEMAP / LDpred / LDSC formats.",
        formatter_class=RichHelpFormatter,
        parents=[
            get_defaultresourse_parser(),
            get_inputvcf_parser(),
            get_formatter_parser(),
            get_common_out_parser(),
            get_bcftools_binary_parser()
        ],
    )

    # Selecting the function to run
    parser.set_defaults(func=run_formatter_direct)

    # Parse args
    args = parser.parse_args(argv)

    # Run
    return args.func(args)


if __name__ == "__main__":
    main()
