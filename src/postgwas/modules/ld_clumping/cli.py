#!/usr/bin/env python3
import argparse
import sys
from postgwas.cli.compute import get_compute_parser
from postgwas.core.ui import AlignedRichHelpFormatter, format_cli_examples

# ---------------------------------------------------------
# IMPORT SHARED PARSERS
# ---------------------------------------------------------
from postgwas.cli.common import (
    get_inputvcf_parser,
    get_common_out_parser,
    get_bcftools_binary_parser,
    get_ld_clump_parser,
    get_population_parser,

)


from postgwas.modules.ld_clumping.service import (
    run_ld_clump_direct)



def build_parser() -> argparse.ArgumentParser:
    return argparse.ArgumentParser(
        prog="postgwas ld_clump",
        usage="postgwas ld_clump --vcf PATH --ld-folder PATH [options]",
        description=(
            "Identify approximately independent association signals by LD clumping "
            "a harmonised GWAS-VCF against the configured reference population."
        ),
        epilog=format_cli_examples(
            (
                "Run LD clumping:",
                "postgwas ld_clump",
                (
                    "--vcf study.vcf.gz",
                    "--ld-folder reference/ld",
                    "--population EUR",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Export the LD-clumping pipeline configuration:",
                "postgwas config export",
                ("--pipeline ld_clump", "--style full", "--output ld_clump_pipeline.yaml"),
            ),
        ),
        formatter_class=AlignedRichHelpFormatter,
        parents=[
            get_compute_parser(),
            get_bcftools_binary_parser(),
            get_inputvcf_parser(),
            get_common_out_parser(),
            get_population_parser(),
            get_ld_clump_parser(add_help=False)
        ],
    )


def main():
    parser = build_parser()

    # ---------------------------------------------------------
    # EXECUTE
    # ---------------------------------------------------------
    # If no arguments provided → show help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    # Parse all arguments
    args = parser.parse_args()

    # We set the function to be executed explicitly:
    run_ld_clump_direct(args)


if __name__ == "__main__":
    main()


__all__ = ["build_parser", "main"]
