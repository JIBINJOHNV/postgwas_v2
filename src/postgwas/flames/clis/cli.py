#!/usr/bin/env python3

from __future__ import annotations
import argparse
import sys

from rich_argparse import ArgumentDefaultsRichHelpFormatter

from postgwas.flames.clis.cli_banner import print_flames_banner
from postgwas.flames.clis.cli_common import add_common_options
from postgwas.flames.clis.cli_direct import add_direct_mode
from postgwas.flames.clis.cli_pipeline import add_pipeline_mode


# ---------------------------------------------------------
# PARENT PARSER (Global + Model Options)
# ---------------------------------------------------------
def build_parent_parser():
    """
    Builds parent parser containing common global options.
    Used as parent for both subcommands (OPTION B).
    """
    parent = argparse.ArgumentParser(
        add_help=False,
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    add_common_options(parent)
    return parent


# ---------------------------------------------------------
# MAIN PARSER
# ---------------------------------------------------------
def build_parser():
    parent = build_parent_parser()

    parser = argparse.ArgumentParser(
        prog="flames",
        description="FLAMES: Fine-mapped Locus Assessment Model of Effector geneS.",
        formatter_class=ArgumentDefaultsRichHelpFormatter,
        add_help=False,            # We manually handle -h
        parents=[parent],
    )

    # Top-level help flag
    parser.add_argument(
        "-h", "--help",
        action="store_true",
        help="Show help message and exit."
    )

    # Subcommands
    subparsers = parser.add_subparsers(
        title="Subcommands",
        dest="cmd",
        metavar=""                # hide metavar
    )

    # Register each mode
    add_direct_mode(subparsers, parent)
    add_pipeline_mode(subparsers, parent)

    return parser


# ---------------------------------------------------------
# MAIN EXECUTION CONTROL
# ---------------------------------------------------------
def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = build_parser()

    # Parse only the first layer
    args, remaining = parser.parse_known_args(argv)

    # ==========================================================
    # CASE 1 — Top-level help: flames --help
    # ==========================================================
    if args.help and args.cmd is None:
        print_flames_banner()     # show banner only
        return  

    # ==========================================================
    # CASE 2 — Subcommand help, e.g.: flames pipeline --help
    # ==========================================================
    if args.cmd and (args.help or "--help" in remaining or "-h" in remaining):

        # Lookup correct subparser (robust, non-private)
        subparser = None
        for action in parser._actions:
            if isinstance(action, argparse._SubParsersAction):
                subparser = action.choices.get(args.cmd)
                break

        if subparser is not None:
            subparser.print_help()
            return

    # ==========================================================
    # CASE 3 — No subcommand provided → Banner
    # ==========================================================
    if args.cmd is None:
        print_flames_banner()
        return

    # ==========================================================
    # CASE 4 — Normal execution
    # ==========================================================
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
