#!/usr/bin/env python3
"""
PostGWAS Fine-Mapping CLI
Supports two engines:
  • SuSiE  (Python/R-based)
  • FINEMAP (C++ binary)
Works in:
  • DIRECT mode   → run fine-mapping on existing locus + sumstats
  • PIPELINE mode → Harmonisation → LD Blocks → QC → Formatter → Fine-mapping
"""

import argparse
import sys,os
from pathlib import Path
from rich_argparse import RichHelpFormatter


# =========================================================
# BACKEND RUNNERS
# =========================================================
from postgwas.finemap.susie.workflows import (
    run_parallel_susie,
)

# Shared parsers
from postgwas.clis.common_cli import (
    get_defaultresourse_parser,
    get_common_out_parser,
    get_plink_binary_parser,
)

from postgwas.clis.finemap_cli import (
    get_finemap_common_parser,
    get_common_susie_arguments,
    get_common_finemap_finemap_arguments,
    get_finemap_susie_inputs_parser,
    get_finemap_finemap_arguments,
    )

from postgwas.utils.main import validate_path


# =========================================================
# MAIN CLI - DIRECT MODE ONLY
# =========================================================
def main():

    # 1. Setup the single, top-level parser, including ONLY Direct Mode components
    parser = argparse.ArgumentParser(
        prog="postgwas-finemap",
        description=(
            "SuSiE / FINEMAP fine-mapping module.\n"
            "Runs fine-mapping using existing sumstats and locus file."
        ),
        formatter_class=RichHelpFormatter,
        parents=[
            # --- Common Components ---
            get_defaultresourse_parser(),
            get_common_out_parser(),
            get_finemap_common_parser(), # Defines --finemap_method
            get_finemap_susie_inputs_parser(), # Defines sumstat/locus file inputs
            get_common_susie_arguments(),
            get_common_finemap_finemap_arguments(),
            get_finemap_finemap_arguments(),
            get_plink_binary_parser(),
            
            # --- Direct Mode specific inputs ---
        ],
    )

    # If no arguments provided → show help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    
    # 3. Execution and Dispatch
    args = parser.parse_args()

    # Dispatch logic based only on --finemap_method
    if "susie" in args.finemap_method:
        run_parallel_susie(args)
    elif "finemap" in args.finemap_method :
        from postgwas.finemap.finemap.workflows import (
            run_finemap_pipeline)
        run_finemap_pipeline(args)
    else:
        # This branch should ideally not be reached if choices are used for --finemap_method
        print(f"Error: Unknown fine-mapping method '{args.finemap_method}'.")


if __name__ == "__main__":
    main()
    
