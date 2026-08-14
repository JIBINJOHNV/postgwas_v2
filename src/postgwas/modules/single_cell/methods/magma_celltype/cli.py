"""Command-line options owned by the MAGMA cell-typing method."""

from __future__ import annotations

import argparse

from postgwas.core.ui import help_with_default


def add_magma_celltype_arguments(
    parser: argparse.ArgumentParser,
    defaults,
) -> None:
    """Add MAGMA cell-typing inputs and settings without CLI defaults."""
    module = defaults.modules.single_cell
    inputs = parser.add_argument_group("MAGMA cell-type inputs")
    inputs.add_argument(
        "--single-cell-covariates",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "Whitespace-delimited MAGMA gene-covariate matrix with gene IDs in "
            "column one, one average-expression column, and one or more cell-type "
            "columns."
        ),
    )
    settings = parser.add_argument_group("MAGMA cell-type settings")
    settings.add_argument(
        "--average-property",
        metavar="COLUMN",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Exact covariate column containing expression averaged across cell types",
            module.magma_celltype.average_property,
        ),
    )
    settings.add_argument(
        "--cell-type-correction",
        nargs="+",
        choices=("bonferroni", "sidak", "holm", "fdr_bh"),
        metavar="METHOD",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Corrections applied across all tested cell types in this dataset",
            " ".join(module.multiple_testing.methods),
        ),
    )
    settings.add_argument(
        "--primary-cell-type-correction",
        choices=("bonferroni", "sidak", "holm", "fdr_bh"),
        metavar="METHOD",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Primary adjusted-p method recorded for downstream reporting",
            module.multiple_testing.primary_method,
        ),
    )
    settings.add_argument(
        "--cell-type-significance-threshold",
        metavar="P",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Adjusted-p reporting threshold; it does not filter result rows",
            module.multiple_testing.significance_threshold,
        ),
    )


__all__ = ["add_magma_celltype_arguments"]
