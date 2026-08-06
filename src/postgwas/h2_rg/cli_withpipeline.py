#!/usr/bin/env python3
"""
PostGWAS — LDSC Heritability / Genetic Correlation CLI

Modes:
  • direct   – Run LDSC on existing munged .sumstats.gz
  • pipeline – Run munge_sumstats + LDSC (h2 or rg) via Docker / local
  • tools    – Utility helpers (e.g. show docker cmds, check refs, etc.)

Structured to mirror the postgwas-finemap CLI.
"""

import argparse
import sys
from rich_argparse import RichHelpFormatter

# =========================================================
# BACKEND RUNNERS
# =========================================================
from postgwas.h2_rg.workflows import (
    run_ldsc_direct,
    run_ldsc_pipeline,
)

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
    Tools mode:
        Gets the same arguments; backend can decide what to use.
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    grp = parser.add_argument_group("1. LDSC Inputs")

    grp.add_argument(
        "--sumstats",
        metavar=" ",
        type=validate_path(must_exist=True, must_be_file=True),
        help=(
            "Summary statistics file.\n"
            "  • DIRECT mode: LDSC-munged .sumstats.gz file.\n"
            "  • PIPELINE mode: LDSC-input TSV (from vcf_to_ldsc)."
        ),
    )

    return parser


# =========================================================
# MAIN CLI
# =========================================================
def main():

    top = argparse.ArgumentParser(
        prog="postgwas-ldsc",
        description="LDSC heritability / genetic correlation module.",
        formatter_class=RichHelpFormatter,
    )

    sub = top.add_subparsers(dest="mode", help="Choose execution mode")

    # ---------------------------------------------------
    # DIRECT MODE
    # ---------------------------------------------------
    direct = sub.add_parser(
        "direct",
        help="Run LDSC on existing munged .sumstats.gz (no munging step).",
        parents=[
            get_defaultresourse_parser(),
            get_common_out_parser(),
            get_ldsc_parser(),
            get_ldsc_common_parser(add_help=False),
        ],
        formatter_class=RichHelpFormatter,
    )

    direct.set_defaults(func=run_ldsc_direct)

    # ---------------------------------------------------
    # PIPELINE MODE
    # ---------------------------------------------------
    pipeline = sub.add_parser(
        "pipeline",
        help="Run munge_sumstats → LDSC heritability (full LDSC pipeline).",
        parents=[
            get_defaultresourse_parser(),
            get_common_out_parser(),
            get_ldsc_common_parser(add_help=False),
        ],
        formatter_class=RichHelpFormatter,
    )

    pipeline.set_defaults(func=run_ldsc_pipeline)

    # ---------------------------------------------------
    # EXECUTION
    # ---------------------------------------------------
    args = top.parse_args()

    # No subcommand → show top-level help
    if args.mode is None:
        top.print_help()
        sys.exit(0)

    # If user typed only "direct" or "pipeline" etc → show mode-specific help
    if args.mode in ["direct", "pipeline"] and len(sys.argv) == 2:
        top.parse_args([args.mode, "--help"])
        sys.exit(0)

    # Dispatch to backend
    args.func(args)


if __name__ == "__main__":
    main()
