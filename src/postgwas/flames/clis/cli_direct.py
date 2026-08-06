# postgwas/flames/clis/cli_direct.py

from __future__ import annotations
import argparse
from rich.console import Console
from rich_argparse import ArgumentDefaultsRichHelpFormatter
from postgwas.utils.main import validate_path


# ---------------------------------------------------------
# EXECUTION FUNCTION
# ---------------------------------------------------------
def run_flames_direct(args):
    console = Console()
    console.rule("[bold magenta]FLAMES — Direct Mode[/bold magenta]")
    console.print(vars(args))


# ---------------------------------------------------------
# DIRECT MODE SUBPARSER
# ---------------------------------------------------------
def get_flame_direct_parser(subparsers, parent):
    """
    Add the FLAMES 'direct' subcommand.
    Includes:
      • Inherited global & model arguments
      • Direct-only controls
      • No metavars
      • Rich + default substitution support
    """

    p = subparsers.add_parser(
        "direct",
        help="Run FLAMES using precomputed SuSiE, MAGMA and PoPS outputs.",
        description=(
            "[bold cyan]DIRECT MODE[/bold cyan]\n\n"
            "Use this mode when precomputed data are already available:\n"
            " • SuSiE credible sets\n"
            " • MAGMA gene-level results\n"
            " • PoPS gene scores\n\n"
            "FLAMES integrates these to produce ML-based effector gene prioritisation."
        ),
        formatter_class=ArgumentDefaultsRichHelpFormatter,
        add_help=False,          # custom -h handling
        parents=[parent],        # inherit global + model settings
    )

    # -----------------------------------------------------
    # REQUIRED INPUTS
    # -----------------------------------------------------
    req = p.add_argument_group("Required Input Files")

    req.add_argument(
        "--susie",
        type=validate_path(must_exist=True, must_be_file=True),
        help="SuSiE credible-set summary file.",
        metavar=""
    )

    req.add_argument(
        "--magma",
        type=validate_path(must_exist=True, must_be_file=True),
        help="MAGMA gene-level results file.",
        metavar=""
    )

    req.add_argument(
        "--pops",
        type=validate_path(must_exist=True, must_be_file=True),
        help="PoPS gene-prioritisation file.",
        metavar=""
    )

    # -----------------------------------------------------
    # OPTIONAL SETTINGS
    # -----------------------------------------------------
    opt = p.add_argument_group("Optional Settings")

    opt.add_argument(
        "--label",
        type=str,
        help="Optional label/tag for this analysis run.",
        metavar=""
    )

    opt.add_argument(
        "--dry_run",
        action="store_true",
        help="Validate input files without running FLAMES."
    )

    # Execution handler
    p.set_defaults(func=run_flames_direct)
    return p
