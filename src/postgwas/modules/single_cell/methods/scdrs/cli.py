"""Command-line options owned by the scDRS method."""

from __future__ import annotations

import argparse

from postgwas.core.ui import help_with_default


def add_scdrs_arguments(parser: argparse.ArgumentParser, defaults) -> None:
    """Add scDRS inputs and settings without introducing CLI-owned defaults."""
    method = defaults.modules.single_cell.scdrs
    inputs = parser.add_argument_group("scDRS inputs")
    inputs.add_argument(
        "--scdrs-h5ad-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "AnnData H5AD file whose adata.X matrix contains the expression "
            "values declared by --scdrs-matrix-state."
        ),
    )
    inputs.add_argument(
        "--scdrs-gene-set-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help="Native scDRS two-column TRAIT/GENESET .gs file.",
    )
    inputs.add_argument(
        "--scdrs-magma-gene-results-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "Headered MAGMA gene result containing configured gene IDs and "
            "ZSTAT values; used when --scdrs-gene-set-source magma is selected."
        ),
    )
    inputs.add_argument(
        "--scdrs-gene-id-map",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "Pinned one-to-one identifier crosswalk from MAGMA genes to the "
            "identifiers in adata.var_names."
        ),
    )
    inputs.add_argument(
        "--scdrs-covariate-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "Optional tab-delimited numeric scDRS covariates keyed by exact "
            "adata.obs_names and containing the configured constant column."
        ),
    )
    inputs.add_argument(
        "--scdrs",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "scDRS executable path or command name",
            defaults.resources.executables.scdrs,
        ),
    )
    settings = parser.add_argument_group("scDRS settings")
    settings.add_argument(
        "--scdrs-h5ad-species",
        choices=("human", "hsapiens", "mouse", "mmusculus"),
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Species of identifiers in the H5AD file",
            method.h5ad_species,
        ),
    )
    settings.add_argument(
        "--scdrs-gene-set-source",
        choices=("file", "magma"),
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Use an exact .gs file or construct it from MAGMA Z statistics",
            method.magma_gene_set.source,
        ),
    )
    settings.add_argument(
        "--scdrs-source-gene-id-type",
        metavar="TYPE",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Declared MAGMA gene-identifier namespace",
            method.magma_gene_set.source_identifier_type,
        ),
    )
    settings.add_argument(
        "--scdrs-target-gene-id-type",
        metavar="TYPE",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Declared adata.var_names identifier namespace",
            method.magma_gene_set.target_identifier_type,
        ),
    )
    settings.add_argument(
        "--scdrs-gene-set-species",
        choices=("human", "hsapiens", "mouse", "mmusculus"),
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Species of identifiers in the scDRS gene sets",
            method.gene_set_species,
        ),
    )
    settings.add_argument(
        "--scdrs-matrix-state",
        choices=("raw_counts", "normalized_log1p"),
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Declared state of adata.X; raw counts are normalized by scDRS",
            method.matrix_state,
        ),
    )
    settings.add_argument(
        "--scdrs-group-analysis",
        nargs="+",
        metavar="OBS_COLUMN",
        default=argparse.SUPPRESS,
        help="One or more categorical adata.obs columns for group analysis.",
    )
    settings.add_argument(
        "--scdrs-correlation-analysis",
        nargs="+",
        metavar="OBS_COLUMN",
        default=argparse.SUPPRESS,
        help="One or more numeric adata.obs columns for correlation analysis.",
    )
    settings.add_argument(
        "--scdrs-gene-analysis",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Correlate gene expression with scDRS disease scores",
            method.downstream.gene_analysis,
        ),
    )
    settings.add_argument(
        "--scdrs-control-gene-sets",
        type=int,
        metavar="N",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Number of matched control gene sets used by scDRS",
            method.control_gene_sets,
        ),
    )


__all__ = ["add_scdrs_arguments"]
