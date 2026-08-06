#!/usr/bin/env python3
import sys
from rich_argparse import RichHelpFormatter
import argparse

# ---------------------------------------------------------
# IMPORTS
# ---------------------------------------------------------
from postgwas.clis.common_cli import ( get_defaultresourse_parser,
                                    get_common_out_parser,
                                    sumstat_summary_arg_parser,
                                    get_inputvcf_parser)
from postgwas.qc_summary.workflows import run_qc_summary_direct


# =========================================================
# MAIN CLI
# =========================================================
def main():
    parser = argparse.ArgumentParser(
        prog="qc_summary",
        description=(
            "QC Summary Module\n\n"
            "This module performs summary-level QC on VCF files using bcftools. "
        ),
        formatter_class=RichHelpFormatter,
    )
    
    parser = argparse.ArgumentParser(
        prog="qc_summary",
        description="postgwas variant qc summary",
        parents=[
            get_defaultresourse_parser(),
            get_inputvcf_parser(),
            get_common_out_parser(),
            sumstat_summary_arg_parser()
        ],
        formatter_class=RichHelpFormatter,
    )

    # ---------------------------------------------------
    # NO ARGS → HELP
    # ---------------------------------------------------
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    # ---------------------------------------------------
    # EXECUTION (DIRECT ONLY)
    # ---------------------------------------------------
    run_qc_summary_direct(args)


if __name__ == "__main__":
    main()

