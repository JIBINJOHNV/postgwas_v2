#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
from rich_argparse import RichHelpFormatter

# ---------------------------------------------------------
# IMPORT SHARED PARSERS
# ---------------------------------------------------------
from postgwas.clis.common_cli import (
    get_defaultresourse_parser,
    get_inputvcf_parser,
    get_common_out_parser,
    get_bcftools_binary_parser,
    get_ld_clump_parser

)


from postgwas.ld_clump.workflows import (
    run_ld_clump_direct)



def main():
    # ---------------------------------------------------------
    # Top-level parser (Combined Direct/Pipeline arguments)
    # ---------------------------------------------------------
    parser = argparse.ArgumentParser(
        prog="ld_clump",
        description=(
            "LD Clumping Module (Direct Run Mode)\n\n"
            "Performs LD pruning using precomputed LD blocks or standard r² + window."
        ),
        formatter_class=RichHelpFormatter,
        parents=[
            # These are the combined parsers required for the full functionality:
            get_defaultresourse_parser(),
            get_bcftools_binary_parser(),
            get_inputvcf_parser(),
            get_common_out_parser(),
            get_ld_clump_parser(add_help=False)
        ]
    )

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
    # Ensure necessary imports are present for the script to run
    import argparse
    import sys
    main()
