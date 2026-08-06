#!/usr/bin/env python3
"""
PostGWAS — MAGMA Gene-Property (Covariate) Analysis

DIRECT MODE:
    Run covariate analysis using an existing MAGMA .genes.raw file.

PIPELINE MODE:
    Lightweight MAGMA pipeline using available CLI modules.
"""

import argparse
from rich_argparse import RichHelpFormatter
from postgwas.utils.main import validate_path
from postgwas.clis.common_cli import safe_parse_args

# =============================
# Allowed Common CLI Builders
# =============================
from postgwas.clis.common_cli import (
    get_defaultresourse_parser,
    get_common_out_parser,
    get_common_magma_covar_parser,
    get_magma_binary_parser,
    get_plink_binary_parser,
    get_bcftools_binary_parser,
)

# MAGMA Covariate execution
from postgwas.magmacovar.workflows import (
    run_magma_covar_direct,
)


# ------------------------------------------------------------
# Add MAGMA gene_results file argument (with your required name)
# ------------------------------------------------------------
def get_magma_covar_gene_results_parser(add_help=False):
    """
    Adds the required --magama_gene_assoc_raw argument for MAGMA covariate analysis.
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("MAGMA Covariate Inputs")

    group.add_argument(
        "--magama_gene_assoc_raw",
        metavar="",
        type=validate_path(
            must_exist=True,
            must_be_file=True,
            must_not_be_empty=True,
        ),
        required=True,
        help=(
            "MAGMA gene-level results file (.genes.raw). "
            "Must be produced by MAGMA gene-level association."
        ),
    )

    return parser


# ===================================================================
# MAIN ENTRYPOINT
# ===================================================================

def main() -> None:
    # --------------------------------------------------------------
    # Top-level parser (MAGMA Covariate Direct Mode Only)
    # --------------------------------------------------------------
    parser = argparse.ArgumentParser(
        prog="magma-covar",
        description=(
            "[bold]MAGMA Covariate (Gene-Property) Analysis[/bold]\n"
            "Run only the covariate analysis step on existing gene results."
        ),
        formatter_class=RichHelpFormatter,
        parents=[
            # Parsers required for a MAGMA Covariate Direct run:
            get_defaultresourse_parser(),
            get_common_out_parser(),
            get_magma_binary_parser(),
            # Covariate-specific inputs/parameters:
            get_common_magma_covar_parser(add_help=False),
            get_magma_covar_gene_results_parser(add_help=False),
        ],
    )
    # --------------------------------------------------------------
    # Parse + Execute
    # --------------------------------------------------------------
    # The safe_parse_args function should be defined elsewhere in your project
    args = safe_parse_args(parser)
    # The execution logic is simplified for single-mode operation.
    # We call the direct function with the parsed arguments.
    run_magma_covar_direct(args)


if __name__ == "__main__":
    main()




