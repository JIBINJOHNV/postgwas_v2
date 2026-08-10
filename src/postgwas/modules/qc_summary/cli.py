#!/usr/bin/env python3
import sys
import argparse
from postgwas.cli.compute import get_compute_parser
from postgwas.core.ui import AlignedRichHelpFormatter, format_cli_examples

# ---------------------------------------------------------
# IMPORTS
# ---------------------------------------------------------
from postgwas.cli.common import (
                                    get_common_out_parser,
                                    sumstat_summary_arg_parser,
                                    get_inputvcf_parser)
from postgwas.modules.qc_summary.service import run_qc_summary_direct


# =========================================================
# MAIN CLI
# =========================================================
def build_parser() -> argparse.ArgumentParser:
    return argparse.ArgumentParser(
        prog="postgwas qc",
        usage="postgwas qc --vcf PATH [options]",
        description="Calculate summary-level quality-control metrics for a GWAS-VCF.",
        epilog=format_cli_examples(
            (
                "Summarise one harmonised GWAS-VCF:",
                "postgwas qc",
                (
                    "--vcf study.vcf.gz",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
        ),
        parents=[
            get_compute_parser(),
            get_inputvcf_parser(),
            get_common_out_parser(),
            sumstat_summary_arg_parser()
        ],
        formatter_class=AlignedRichHelpFormatter,
    )


def main():
    parser = build_parser()

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


__all__ = ["build_parser", "main"]
