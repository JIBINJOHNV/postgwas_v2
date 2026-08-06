#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
from datetime import datetime
import subprocess
from rich_argparse import RichHelpFormatter

# ---------------------------------------------------------
# SHARED PARSERS
# ---------------------------------------------------------
from postgwas.clis.common_cli import (
    get_defaultresourse_parser,
    get_inputvcf_parser,
    get_common_out_parser,
    get_assoc_plot_parser,
    get_genomeversion_parser
)

from postgwas.harmonisation.cli import (
    get_harmonisation_parser,
    run_harmonisation
)

from postgwas.annot_ldblock.cli import (
    get_annot_ldblock_parser
    )

# Import the parser we just created
from postgwas.manhattan.workflows import run_assoc_plot_direct




# =====================================================================
#   MAIN CLI
# =====================================================================
# =========================================================
# MAIN CLI - DIRECT MODE ONLY
# =========================================================
def main():
    # 1. Setup the single, top-level parser, including ONLY Direct Mode components
    parser = argparse.ArgumentParser(
        prog="assoc-plot",
        description=(
            "Association Manhattan / CSQ / Allelic-shift plotting module (Direct Mode Only).\n"
            "Plots from an existing VCF file."
        ),
        formatter_class=RichHelpFormatter,
        parents=[
            # Components from the original 'direct' subcommand
            get_defaultresourse_parser(),
            get_inputvcf_parser(),
            get_genomeversion_parser(),
            get_common_out_parser(),
            get_assoc_plot_parser()
            
            # Note: You need to ensure all parser builder functions
            # accept 'add_help=False' to avoid the TypeError seen previously.
        ]
    )
    # 2. Add explicit required arguments here if the parent parsers don't enforce them.
    #    (e.g., if get_inputvcf_parser doesn't use required=True, you'd add it here)
    # If no arguments provided → show help (Uses the existing logic)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    # 3. Execution and Dispatch
    args = parser.parse_args()
    # The function is always run_assoc_plot_direct
    run_assoc_plot_direct(args)

if __name__ == "__main__":
    main()






