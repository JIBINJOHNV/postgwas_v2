#!/usr/bin/env python3
import argparse
import sys
from postgwas.cli.compute import get_compute_parser
from postgwas.core.ui import AlignedRichHelpFormatter, format_cli_examples

# ---------------------------------------------------------
# SHARED PARSERS
# ---------------------------------------------------------
from postgwas.cli.common import (
    get_inputvcf_parser,
    get_common_out_parser,
    get_assoc_plot_parser,
    get_genome_build_parser
)

from postgwas.modules.manhattan.service import run_assoc_plot_direct




# =====================================================================
#   MAIN CLI
# =====================================================================
# =========================================================
# MAIN CLI - DIRECT MODE ONLY
# =========================================================
def build_parser() -> argparse.ArgumentParser:
    return argparse.ArgumentParser(
        prog="postgwas manhattan",
        usage="postgwas manhattan --vcf PATH [--png PATH | --pdf PATH] [options]",
        description=(
            "Create Manhattan-style association plots from a harmonised GWAS-VCF, "
            "with optional consequence and allelic-shift annotations."
        ),
        epilog=format_cli_examples(
            (
                "Create a PNG Manhattan plot:",
                "postgwas manhattan",
                (
                    "--vcf study.vcf.gz",
                    "--png STUDY_manhattan.png",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Highlight consequence annotations:",
                "postgwas manhattan",
                ("--vcf study_annotated.vcf.gz", "--png STUDY_csq.png", "--csq"),
            ),
        ),
        formatter_class=AlignedRichHelpFormatter,
        parents=[
            get_compute_parser(),
            get_inputvcf_parser(),
            get_genome_build_parser(),
            get_common_out_parser(),
            get_assoc_plot_parser(),
        ],
    )


def main():
    parser = build_parser()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    # 3. Execution and Dispatch
    args = parser.parse_args()
    # The function is always run_assoc_plot_direct
    run_assoc_plot_direct(args)

if __name__ == "__main__":
    main()


__all__ = ["build_parser", "main"]


