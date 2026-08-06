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
    get_inputvcf_parser,
    get_genomeversion_parser,
    get_common_out_parser,
    get_common_magma_covar_parser,
    get_common_magma_assoc_parser,
    get_common_sumstat_filter_parser,
    get_common_pops_parser,
    get_common_susie_arguments,
    get_finemap_method_parser,
    get_magma_binary_parser,
    get_plink_binary_parser,
    get_bcftools_binary_parser,
)

# MAGMA Covariate execution
from postgwas.magmacovar.workflows import (
    run_magma_covar_direct,
    run_magma_covar_pipeline,
)


from postgwas.pops.workflows import run_pops_direct, run_pops_pipeline
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

    top = argparse.ArgumentParser(
        prog="pops",
        description="PoPS Gene Prioritisation — Direct mode or Full PostGWAS Pipeline.",
        formatter_class=RichHelpFormatter,
    )

    sub = top.add_subparsers(dest="mode", help="Choose execution mode")

    # --------------------------------------------------------
    # DIRECT MODE
    # --------------------------------------------------------
    direct = sub.add_parser(
        "direct",
        help="Run PoPS directly using existing MAGMA outputs.",
        parents=[
            get_defaultresourse_parser(),
            get_common_out_parser(),
            get_pops_direct_parser(add_help=False),
            get_common_pops_parser(add_help=False),
        ],
        formatter_class=RichHelpFormatter,
    )
    direct.set_defaults(func=run_pops_direct)

    # --------------------------------------------------------
    # PIPELINE MODE
    # --------------------------------------------------------
    pipeline = sub.add_parser(
        "pipeline",
        help="Run full pipeline → Harmonisation → LD → QC → Formatter → MAGMA → PoPS",
        parents=[
            get_defaultresourse_parser(),
            get_genomeversion_parser(),
            get_common_magma_assoc_parser(add_help=False),
            get_common_pops_parser(add_help=False),
            get_common_out_parser(),
        ],
        formatter_class=RichHelpFormatter,
    )
    pipeline.set_defaults(func=run_pops_pipeline)

    # --------------------------------------------------------
    # EXECUTION LOGIC
    # --------------------------------------------------------
    args = top.parse_args()

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
