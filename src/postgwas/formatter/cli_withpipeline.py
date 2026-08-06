#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
from rich_argparse import RichHelpFormatter

# ---- Shared CLI imports ----
from postgwas.clis.common_cli import (
    get_defaultresourse_parser,
    get_inputvcf_parser,
    get_genomeversion_parser,
    get_common_out_parser,
)

# ---- Step 1–3 parser imports for pipeline mode ----
from postgwas.harmonisation.cli import get_harmonisation_parser
from postgwas.annot_ldblock.cli import get_annot_ldblock_parser
from postgwas.formatter.workflows import run_formatter_direct, run_formatter_pipeline

# ============================================================
#  FORMATTER PARSER (Step-4 specific options)
# ============================================================
def get_formatter_parser(add_help=False):
    """
    Parser that defines formatter-specific arguments:
    --format   (magma/finemap)
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("Formatter Options")
    group.add_argument(
        "--format",
            choices=["magma", "finemap","ldpred","ldsc"],
            metavar="{magma, finemap, ldpred, ldsc}",
            help=(
                "[bold bright_red]Required[/bold bright_red]: Choose output format for conversion. "
                "[bold]Choices:[/bold] [bright_yellow]{magma, finemap, ldpred, ldsc}[/bright_yellow]."
        )
    )
    return parser


# ============================================================
#  MAIN CLI
# ============================================================
def main():

    top_parser = argparse.ArgumentParser(
        prog="formatter",
        description="Convert GWAS VCF to MAGMA or FINEMAP formats (Direct or Pipeline).",
        formatter_class=RichHelpFormatter
    )
    subparsers = top_parser.add_subparsers(dest="mode", help="Choose execution mode")
    # --------------------------------------------------------
    # DIRECT MODE
    # --------------------------------------------------------
    direct_parser = subparsers.add_parser(
        "direct",
        help="Convert GWAS VCF directly → MAGMA or FINEMAP",
        parents=[
            get_defaultresourse_parser(),
            get_formatter_parser(),
            get_inputvcf_parser(),
            get_common_out_parser(),
        ],
        formatter_class=RichHelpFormatter
    )
    direct_parser.set_defaults(func=run_formatter_direct)
    # --------------------------------------------------------
    # PIPELINE MODE
    # --------------------------------------------------------
    pipeline_parser = subparsers.add_parser(
        "pipeline",
        help="Run Harmonisation → LD Blocks → QC Filter → Format",
        parents=[
            get_defaultresourse_parser(),
            get_harmonisation_parser(),
            get_annot_ldblock_parser(),
            get_formatter_parser(),
            get_common_out_parser(),
        ],
        formatter_class=RichHelpFormatter
    )
    pipeline_parser.set_defaults(func=run_formatter_pipeline)
    # --------------------------------------------------------
    # EXECUTE
    # --------------------------------------------------------
    args = top_parser.parse_args()
    # Auto-help if user does:
    #   formatter direct
    #   formatter pipeline
    if args.mode == "direct" and len(sys.argv) == 2:
        top_parser.parse_args(["direct", "--help"])
        sys.exit(0)
    if args.mode == "pipeline" and len(sys.argv) == 2:
        top_parser.parse_args(["pipeline", "--help"])
        sys.exit(0)
    if hasattr(args, "func"):
        args.func(args)
    else:
        top_parser.print_help()


if __name__ == "__main__":
    main()