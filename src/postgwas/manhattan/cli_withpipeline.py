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
)

from postgwas.harmonisation.cli import (
    get_harmonisation_parser,
    run_harmonisation
)

from postgwas.annot_ldblock.cli import (
    get_annot_ldblock_parser,
    run_annot_ldblock_pipeline
)

# Import the parser we just created
from postgwas.manhattan.workflows import run_assoc_plot_direct,run_assoc_plot_pipeline





def get_assoc_plot_parser(add_help=False):
    """
    Parent-safe parser containing common assoc-plot arguments
    shared between DIRECT and PIPELINE modes.
    Matches the original assoc_plot.R defaults:
      • min-lp default = 2
      • spacing default = 20 (R default)
      • cyto-ratio default = 25
      • loglog-pval default = 10
      • width = 7.0
      • fontsize = 12
      • nauto = 22
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    # ------------------------------------------------------------
    # OUTPUTS
    # ------------------------------------------------------------
    output_group = parser.add_argument_group()

    output_group.add_argument(
        "--pdf",
        metavar=" ",
        help="Output PDF file (e.g., result.pdf)."
    )

    output_group.add_argument(
        "--png",
        metavar=" ",
        help="Output PNG file (e.g., result.png)."
    )
    # ------------------------------------------------------------
    # INPUTS
    # ------------------------------------------------------------
    input_group = parser.add_argument_group("Association Plot Input Arguments")

    input_group.add_argument(
        "--genome",
        choices=["GRCh37", "GRCh38"],
        metavar=" ",
        help=(
            "Genome assembly used for chromosome-length annotation.\n"
            "Choices: [cyan]GRCh37[/cyan], [cyan]GRCh38[/cyan]\n\n"
            "[bold yellow]IMPORTANT:[/bold yellow] You must provide either:\n"
            "  • [cyan]--genome[/cyan]\n"
            "        OR\n"
            "  • [cyan]--cytoband[/cyan]\n"
            "Exactly one is required."
        )
    )

    input_group.add_argument(
        "--cytoband",
        metavar=" ",
        help=(
            "Cytoband annotation file (<cytoband.txt.gz>).\n"
            "Use this [cyan]instead[/cyan] of --genome for custom cytobands."
        )
    )

    # input_group.add_argument(
    #     "--tbx",
    #     metavar=" ",
    #     help="Tabix (.tbi) index corresponding to the input VCF."
    # )

    input_group.add_argument(
        "--pheno",
        metavar=" ",
        help="Phenotype name to extract from GWAS-VCF file."
    )
    # ------------------------------------------------------------
    # FLAGS
    # ------------------------------------------------------------
    flag_group = parser.add_argument_group("Flags")

    flag_group.add_argument(
        "--as",
        dest="allelic_shift",
        action="store_true",
        help="Input VCF contains allelic-shift annotation."
    )

    flag_group.add_argument(
        "--csq",
        action="store_true",
        help="Flag coding variants (from CSQ annotation) in red."
    )
    # ------------------------------------------------------------
    # NUMERIC OPTIONS
    # ------------------------------------------------------------
    numeric_group = parser.add_argument_group("Numeric Plot Options")
    numeric_group.add_argument(
        "--nauto",
        type=int,
        default=22,
        metavar=" ",
        help="Number of autosomes. [bold green]Default:[/bold green] [cyan]22[/cyan]"
    )

    numeric_group.add_argument(
        "--min-af",
        type=float,
        default=0.0,
        metavar=" ",
        help="Minimum allele frequency filter. [bold green]Default:[/bold green] [cyan]0.0[/cyan]"
    )

    numeric_group.add_argument(
        "--min-lp",
        type=int,
        default=2,
        metavar=" ",
        help="Minimum −log10(P) threshold. [bold green]Default:[/bold green] [cyan]2[/cyan]"
    )

    numeric_group.add_argument(
        "--loglog-pval",
        type=int,
        default=10,
        metavar=" ",
        help="-log10(P) threshold for switching to log-log scale. "
             "[bold green]Default:[/bold green] [cyan]10[/cyan]"
    )

    numeric_group.add_argument(
        "--cyto-ratio",
        type=int,
        default=25,
        metavar=" ",
        help="Plot height : cytoband height ratio. "
             "[bold green]Default:[/bold green] [cyan]25[/cyan]"
    )

    numeric_group.add_argument(
        "--max-height",
        type=int,
        metavar=" ",
        help="Maximum vertical height of plot (optional)."
    )

    numeric_group.add_argument(
        "--spacing",
        type=int,
        default=20,
        metavar=" ",
        help="Spacing between chromosomes. "
             "[bold green]Default:[/bold green] [cyan]20[/cyan]"
    )
    # ------------------------------------------------------------
    # FIGURE SIZE
    # ------------------------------------------------------------
    fig_group = parser.add_argument_group("Figure Dimensions")

    fig_group.add_argument(
        "--width",
        type=float,
        default=7.0,
        metavar=" ",
        help="Plot width in inches. [bold green]Default:[/bold green] [cyan]7.0[/cyan]"
    )

    fig_group.add_argument(
        "--height",
        type=float,
        metavar=" ",
        help="Plot height in inches (optional)."
    )

    fig_group.add_argument(
        "--fontsize",
        type=int,
        default=12,
        metavar=" ",
        help="Font size. [bold green]Default:[/bold green] [cyan]12[/cyan]"
    )

    return parser




# =====================================================================
#   MAIN CLI
# =====================================================================
def main():
    top = argparse.ArgumentParser(
        prog="assoc-plot",
        description="Association Manhattan / CSQ / Allelic-shift plotting module.",
        formatter_class=RichHelpFormatter
    )

    sub = top.add_subparsers(dest="mode", help="direct or pipeline")

    # ------------------------------------------------------------
    # DIRECT
    # ------------------------------------------------------------
    direct = sub.add_parser(
        "direct",
        help="Plot from an existing VCF.",
        parents=[
            get_defaultresourse_parser(),
            get_inputvcf_parser(),
            get_common_out_parser(),
            get_assoc_plot_parser(add_help=False)
        ],
        formatter_class=RichHelpFormatter
    )
    direct.set_defaults(func=run_assoc_plot_direct)

    # ------------------------------------------------------------
    # PIPELINE
    # ------------------------------------------------------------
    pipeline = sub.add_parser(
        "pipeline",
        help="Harmonisation → Annot-LD → Assoc Plot.",
        parents=[
            get_defaultresourse_parser(),
            get_harmonisation_parser(),
            get_annot_ldblock_parser(),
            get_common_out_parser(),
            get_assoc_plot_parser(add_help=False)
        ],
        formatter_class=RichHelpFormatter
    )
    pipeline.set_defaults(func=run_assoc_plot_pipeline)

    # ============================================================
    # EXECUTE
    # ============================================================
    args = top.parse_args()

    if len(sys.argv) == 1:
        top.print_help()
        sys.exit(0)

    if args.mode == "direct" and len(sys.argv) == 2:
        top.parse_args(["direct", "--help"])
        sys.exit(0)

    if args.mode == "pipeline" and len(sys.argv) == 2:
        top.parse_args(["pipeline", "--help"])
        sys.exit(0)

    if hasattr(args, "func"):
        args.func(args)
    else:
        top.print_help()


if __name__ == "__main__":
    main()






