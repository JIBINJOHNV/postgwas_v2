#!/usr/bin/env python3
"""
PostGWAS — LDSC Heritability CLI

Mode:
  • Direct – Munge formatter output and run single-trait LDSC heritability
"""
import argparse
import sys

from pydantic import ValidationError
from rich.console import Console
from rich.markup import escape

from postgwas.cli.compute import get_compute_parser
from postgwas.core.errors import ConfigurationError
from postgwas.core.ui import AlignedRichHelpFormatter, format_cli_examples

# =========================================================
# BACKEND RUNNERS
# =========================================================
from postgwas.modules.ldsc.ldsc_runner import LDSCError
from postgwas.modules.ldsc.service import run_ldsc_direct


# Shared parsers
from postgwas.cli.common import (
    get_common_out_parser,
    get_ldsc_common_parser,
)

from postgwas.core.execution.runtime import validate_path


# =========================================================
# SHARED LDSC INPUT ARGUMENTS
# =========================================================
def get_ldsc_parser(add_help: bool = False) -> argparse.ArgumentParser:
    """
    Defines input arguments for LDSC direct execution.
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    grp = parser.add_argument_group("Formatted summary statistics")

    grp.add_argument(
        "--ldsc-input",
        metavar="PATH",
        default=argparse.SUPPRESS,
        type=validate_path(
            must_exist=True,
            must_be_file=True,
            must_not_be_empty=True,
        ),
        help="REQUIRED. LDSC input TSV created by the PostGWAS formatter.",
    )

    controls = parser.add_argument_group("Configuration")
    controls.add_argument(
        "--run-config",
        metavar="PATH",
        default=argparse.SUPPRESS,
        type=validate_path(must_exist=True, must_be_file=True),
        help=(
            "YAML file containing LDSC or full-run settings. Explicit command-line "
            "values override matching YAML values."
        ),
    )
    return parser


# =========================================================
# MAIN CLI
# =========================================================
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="postgwas heritability",
        description="Estimate single-trait SNP heritability from formatter-created LDSC input.",
        epilog=format_cli_examples(
            (
                "Estimate observed-scale SNP heritability:",
                "postgwas heritability",
                (
                    "--ldsc-input formatted/STUDY_ldsc.tsv.gz",
                    "--merge-alleles reference/w_hm3.snplist",
                    "--ref-ld-chr reference/eur_w_ld_chr",
                    "--w-ld-chr reference/eur_w_ld_chr",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Add prevalence values for liability-scale reporting:",
                "postgwas heritability",
                (
                    "--ldsc-input formatted/STUDY_ldsc.tsv.gz",
                    "--merge-alleles reference/w_hm3.snplist",
                    "--ref-ld-chr reference/eur_w_ld_chr",
                    "--w-ld-chr reference/eur_w_ld_chr",
                    "--samp-prev 0.2",
                    "--pop-prev 0.01",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
        ),
        parents=[
            get_compute_parser(),
            get_common_out_parser(),
            get_ldsc_parser(),
            get_ldsc_common_parser(add_help=False),
        ],
        formatter_class=AlignedRichHelpFormatter,
    )
    for destination in (
        "ldsc_input",
        "merge_alleles",
        "ref_ld_chr",
        "w_ld_chr",
        "dataset_id",
        "output_directory",
    ):
        action = next(
            item for item in parser._actions if item.dest == destination
        )
        action.required = True
    return parser


def main(argv=None) -> int:
    parser = build_parser()

    # ---------------------------------------------------
    # EXECUTION
    # ---------------------------------------------------

    # Check if no arguments provided (sys.argv[0] is the script name)
    arguments = sys.argv[1:] if argv is None else list(argv)
    if not arguments:
        parser.print_help()
        return 0

    args = parser.parse_args(arguments)

    # Directly dispatch to the direct runner
    try:
        run_ldsc_direct(args)
    except (LDSCError, ConfigurationError, OSError, ValidationError) as exc:
        Console(stderr=True).print(
            "\n[bold red]LDSC heritability stopped.[/bold red]\n"
            "[bold]Reason:[/bold] %s\n"
            "[bold]No new LDSC results were published.[/bold]\n"
            % escape(str(exc))
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = ["build_parser", "get_ldsc_parser", "main"]
