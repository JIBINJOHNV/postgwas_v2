import argparse
import sys
from rich_argparse import RichHelpFormatter 
# Assuming utility/logic imports
from postgwas.utils.main import validate_path 
from postgwas.annot_ldblock.annot_ldblock import annotate_ldblocks 
from postgwas.annot_ldblock.workflows import run_annot_ldblock, run_annot_ldblock_pipeline

# --- IMPORTANT IMPORTS FROM STEP 1 (HARMONISATION) ---
from postgwas.harmonisation.cli import get_harmonisation_parser, run_harmonisation 
from postgwas.clis.common_cli import get_defaultresourse_parser,get_inputvcf_parser,get_genomeversion_parser

# ---------------------------------------------------------

# --- 1. BASE PARSER FOR STEP 2 ARGUMENTS ---

def get_annot_ldblock_parser(add_help=False):
    """
    Returns the core ArgumentParser for annot_ldblock arguments, structured using a group. 
    """
    # Create bare parser for inheritance
    parser = argparse.ArgumentParser(add_help=add_help) 
    
    # 💥 DEFINE ARGUMENT GROUP for Step 2 💥
    group = parser.add_argument_group('LD BLOCK ANNOTATION Arguments')

    group.add_argument(
        "--ld_block_population", 
        nargs="+", 
        default=["EUR", "AFR", "EAS"], 
        help=(
            "Populations to annotate (e.g., EUR, AFR). "
            "Files must be named following the pattern: [cyan]Genomeversion_Population_ldetect.bed.gz[/cyan]."
            "[bold green]Default:[/bold green] [cyan]EUR AFR EAS[/cyan]. "
        )
    )
    group.add_argument(
        "--ld-region-dir",
        type=validate_path(must_exist=True, must_be_dir=True),
        metavar='',
        help=(
            "[bold bright_red]Required[/bold bright_red]: Directory containing LD-block BED files. "
            "Each BED file must contain four columns: CHROM, START, END, and Annotation. "
            "The fourth column is used as the LD-block annotation label."
        ),
    )

    
    return parser



# ---------------------------------------------------------
# 3. STANDALONE ENTRY POINT (The "CLI" with Subcommands)
# ---------------------------------------------------------

import argparse
import sys


def main():
    top_parser = argparse.ArgumentParser(
        prog="annot_ldblock",
        description="Run LD Block Annotation in direct or pipeline mode.",
        formatter_class=RichHelpFormatter,
        add_help=False,                     # we add -h manually so it's on top
    )

    # Put -h at the very top of the help
    top_parser.add_argument(
        "-h", "--help",
        action="help",
        default=argparse.SUPPRESS,
        help=(
            "Show main (top-level) help.\n"
            "For detailed help, run:\n"
            "  postgwas annot_ldblock direct --help\n"
            "  postgwas annot_ldblock pipeline --help"
        )
    )

    subparsers = top_parser.add_subparsers(
        dest="mode",
        metavar="{direct,pipeline}",
    )

    # === direct subcommand ===
    direct_parser = subparsers.add_parser(
        "direct",
        help="Run LD block annotation using VCF output from Harmonisation.",
        parents=[
            get_defaultresourse_parser(),
            get_inputvcf_parser(),
            get_genomeversion_parser(),
            get_annot_ldblock_parser(add_help=False),
        ],
        formatter_class=RichHelpFormatter,
    )
    direct_parser.set_defaults(func=run_annot_ldblock)

    # === pipeline subcommand ===
    pipeline_parser = subparsers.add_parser(
        "pipeline",
        help="Run Step 1 (Harmonisation) → Step 2 (LD Block Annotation) of the PostGWAS pipeline.",
        parents=[
            get_defaultresourse_parser(),
            get_harmonisation_parser(),
            get_annot_ldblock_parser(add_help=False),
        ],
        formatter_class=RichHelpFormatter,
    )
    pipeline_parser.set_defaults(func=run_annot_ldblock_pipeline)

    # ------------------------------------------------------------------
    # This is the magic part that makes everything work perfectly
    # ------------------------------------------------------------------
    if len(sys.argv) == 1:
        # python cli.py  → show top-level help
        top_parser.print_help()
        sys.exit(0)

    if len(sys.argv) == 2 and sys.argv[1] in {"direct", "pipeline"}:
        # python cli.py direct   → show help for that subcommand
        top_parser.parse_args([sys.argv[1], "--help"])

    # Normal parsing
    args = top_parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        top_parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()