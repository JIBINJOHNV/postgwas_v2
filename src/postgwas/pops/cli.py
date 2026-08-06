#!/usr/bin/env python3
"""
PostGWAS — PoPS Gene Prioritisation

Supports:
  • DIRECT MODE:
        Run PoPS directly using existing MAGMA outputs and PoPS feature matrices.

  • PIPELINE MODE:
        Full workflow:
           Harmonisation → LD Blocks → QC → Formatter → MAGMA → PoPS

This CLI is fully aligned with your PostGWAS architecture.
"""

import argparse,sys
from rich_argparse import RichHelpFormatter
from postgwas.utils.main import validate_path
from postgwas.clis.common_cli import safe_parse_args

# =============================
# Allowed Common CLI Builders
# =============================
from postgwas.clis.common_cli import (
    get_defaultresourse_parser,
    get_common_out_parser,
    get_common_pops_parser,
)



from postgwas.pops.workflows import run_pops_direct
# =====================================================================
#   DIRECT MODE PARSER
# =====================================================================
def get_pops_direct_parser(add_help=False):
    """
    Parser for PoPS DIRECT MODE.
    User must supply all required PoPS inputs manually.
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("PoPS Direct Inputs")

    group.add_argument(
        "--magma_assoc_prefix",
        metavar=" ",
        help=(
            "[bold bright_red]Required[/bold bright_red]: Prefix to MAGMA gene-level output files.\n"
            "This should be the base prefix only — [bold]do NOT include[/bold] "
            "[cyan].genes.out[/cyan] or [cyan].genes.raw[/cyan].\n\n"
            "PoPS will automatically read:\n"
            "  • {prefix}[cyan].genes.out[/cyan]\n"
            "  • {prefix}[cyan].genes.raw[/cyan]\n\n"
            "[Used for direct PoPS mode]"
        ),
    )
    group.add_argument(
        "--pops_verbose",
        action="store_true",
        help=(
            "Enable detailed PoPS logging.\n"
            "When set, PoPS will print additional progress messages, debugging "
            "information, and feature-selection details.\n"
            "[bold green]Recommended[/bold green] when troubleshooting or running "
            "small test datasets."
        ),
    )
    return parser



# =====================================================================
#   MAIN CLI ENTRY
# =====================================================================
def main():
    # --------------------------------------------------------
    # Top-level parser (PoPS Direct Mode Only)
    # --------------------------------------------------------
    parser = argparse.ArgumentParser(
        prog="pops",
        description="PoPS Gene Prioritisation — Direct mode: Run PoPS directly using existing MAGMA outputs.",
        formatter_class=RichHelpFormatter,
        parents=[
            # Parsers required for a PoPS Direct run:
            get_defaultresourse_parser(),
            get_common_out_parser(),
            # PoPS specific inputs/parameters
            get_pops_direct_parser(add_help=False),
            get_common_pops_parser(add_help=False),
        ]
    )

    # --------------------------------------------------------
    # EXECUTION LOGIC
    # --------------------------------------------------------

    # If no arguments provided → show help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    # Parse all arguments
    args = parser.parse_args()

    # The script now assumes a single execution flow (the previous 'direct' function)
    run_pops_direct(args)


if __name__ == "__main__":
    # Ensure necessary imports are present for the script to run
    # import argparse
    # import sys
    # from rich_argparse import RichHelpFormatter
    # ... all necessary 'get_...' and 'run_pops_direct' imports ...

    main()
