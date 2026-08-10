"""Command-line interface for GWAS and single-cell integration."""

from __future__ import annotations

import argparse

from rich.console import Console

from postgwas.cli.common import (
    get_common_out_parser,
    get_magma_binary_parser,
)
from postgwas.cli.compute import get_compute_parser
from postgwas.config import load_configuration
from postgwas.config.models.modules.single_cell import SINGLE_CELL_TOOLS
from postgwas.core.errors import ConfigurationError
from postgwas.core.ui import (
    AlignedRichHelpFormatter,
    format_cli_examples,
    help_with_default,
)
from postgwas.modules.magmacovar.cli import (
    get_magma_covar_gene_results_parser,
)
from postgwas.modules.magmacovar.errors import MagmaCovarError
from postgwas.modules.single_cell.errors import SingleCellError


def get_single_cell_parser(add_help=False, *, direct_controls=False):
    """Return reusable single-cell options without CLI-owned defaults."""
    defaults = load_configuration()
    module = defaults.modules.single_cell
    parser = argparse.ArgumentParser(add_help=add_help)
    inputs = parser.add_argument_group("Single-cell inputs")
    inputs.add_argument(
        "--single-cell-covariates",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "Whitespace-delimited MAGMA gene-covariate matrix with gene IDs in "
            "column one, one average-expression column, and one or more cell-type "
            "columns. The GWAS pipeline cannot create this biological input."
        ),
    )
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
            "ZSTAT values. The pipeline supplies this from its MAGMA stage."
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
            "adata.obs_names."
        ),
    )
    inputs.add_argument(
        "--scdrs",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "scDRS executable path or command name",
            defaults.resources.executables.scdrs,
            label="Configured default",
        ),
    )
    settings = parser.add_argument_group("Single-cell analysis settings")
    settings.add_argument(
        "--tools",
        nargs="+",
        choices=SINGLE_CELL_TOOLS,
        metavar="TOOL",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "One or more single-cell integration tools to run",
            " ".join(module.tools),
            label="Configured default",
        ),
    )
    settings.add_argument(
        "--average-property",
        metavar="COLUMN",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Exact covariate column containing expression averaged across cell types",
            module.magma_celltype.average_property,
            label="Configured default",
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
            label="Configured default",
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
            label="Configured default",
        ),
    )
    settings.add_argument(
        "--cell-type-significance-threshold",
        metavar="P",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Adjusted-p reporting threshold; it does not filter result rows",
            module.multiple_testing.significance_threshold,
            label="Configured default",
        ),
    )
    settings.add_argument(
        "--scdrs-h5ad-species",
        choices=("human", "hsapiens", "mouse", "mmusculus"),
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Species of identifiers in the H5AD file",
            module.scdrs.h5ad_species,
            label="Configured default",
        ),
    )
    settings.add_argument(
        "--scdrs-gene-set-source",
        choices=("file", "magma"),
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Use an exact .gs file or construct it from MAGMA Z statistics",
            module.scdrs.magma_gene_set.source,
            label="Configured default",
        ),
    )
    settings.add_argument(
        "--scdrs-source-gene-id-type",
        metavar="TYPE",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Declared MAGMA gene-identifier namespace",
            module.scdrs.magma_gene_set.source_identifier_type,
            label="Configured default",
        ),
    )
    settings.add_argument(
        "--scdrs-target-gene-id-type",
        metavar="TYPE",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Declared adata.var_names identifier namespace",
            module.scdrs.magma_gene_set.target_identifier_type,
            label="Configured default",
        ),
    )
    settings.add_argument(
        "--scdrs-gene-set-species",
        choices=("human", "hsapiens", "mouse", "mmusculus"),
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Species of identifiers in the scDRS gene sets",
            module.scdrs.gene_set_species,
            label="Configured default",
        ),
    )
    settings.add_argument(
        "--scdrs-matrix-state",
        choices=("raw_counts", "normalized_log1p"),
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Declared state of adata.X; raw counts are normalized by scDRS",
            module.scdrs.matrix_state,
            label="Configured default",
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
            module.scdrs.downstream.gene_analysis,
            label="Configured default",
        ),
    )
    settings.add_argument(
        "--scdrs-control-gene-sets",
        type=int,
        metavar="N",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Number of matched control gene sets used by scDRS",
            module.scdrs.control_gene_sets,
            label="Configured default",
        ),
    )
    if direct_controls:
        controls = parser.add_argument_group("Configuration and continuation")
        controls.add_argument(
            "--run-config",
            metavar="PATH",
            default=argparse.SUPPRESS,
            help="YAML settings; explicit command-line values override matching keys.",
        )
        controls.add_argument(
            "--resume",
            action=argparse.BooleanOptionalAction,
            default=argparse.SUPPRESS,
            help=help_with_default(
                "Reuse provenance-matched MAGMA and cell-type results",
                defaults.run.resume,
            ),
        )
        controls.add_argument(
            "--overwrite",
            action="store_true",
            default=argparse.SUPPRESS,
            help=help_with_default(
                "Replace outputs owned by this configured single-cell run",
                defaults.run.overwrite,
            ),
        )
    return parser


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="postgwas single_cell",
        usage=(
            "postgwas single_cell --tools TOOL [tool-specific inputs] [options]"
        ),
        description=(
            "Run one or more GWAS/single-cell integration methods. MAGMA cell "
            "typing uses MAGMA gene results and an expression covariate matrix; "
            "scDRS uses an exact H5AD atlas and native .gs gene sets."
        ),
        epilog=format_cli_examples(
            (
                "Run directly from exact MAGMA and cell-type inputs:",
                "postgwas single_cell",
                (
                    "--tools magma_celltype",
                    "--magma-gene-results-file STUDY.genes.raw",
                    "--single-cell-covariates atlas_celltype_average.tsv",
                    "--magma /path/to/magma",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Run native scDRS scores and cell-type group analysis:",
                "postgwas single_cell",
                (
                    "--tools scdrs",
                    "--scdrs-h5ad-file brain_atlas.h5ad",
                    "--scdrs-gene-set-file study.gs",
                    "--scdrs-group-analysis cell_type",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Export reusable settings:",
                "postgwas config export",
                ("--module single_cell", "--style full", "--output single_cell.yaml"),
            ),
        ),
        formatter_class=AlignedRichHelpFormatter,
        parents=[
            get_compute_parser(),
            get_common_out_parser(),
            get_magma_binary_parser(),
            get_magma_covar_gene_results_parser(add_help=False),
            get_single_cell_parser(direct_controls=True),
        ],
    )
    for action in parser._actions:
        if action.dest in {
            "dataset_id", "output_directory", "threads", "memory_gb", "seed",
        }:
            action.default = argparse.SUPPRESS
            action.type = None
    return parser


def get_single_cell_pipeline_examples():
    """Return complete examples for context-sensitive pipeline help."""
    return (
        (
            "Run MAGMA cell typing from a harmonised GWAS-VCF:",
            "postgwas pipeline",
            (
                "--modules single_cell",
                "--tools magma_celltype",
                "--vcf study.vcf.gz",
                "--single-cell-covariates atlas_celltype_average.tsv",
                "--magma-ld-reference reference/g1000_eur",
                "--gene-location-file reference/NCBI37.3.gene.loc",
                "--dataset-id STUDY",
                "--output-directory results",
            ),
        ),
        (
            "Construct a disease gene set from pipeline MAGMA and run scDRS:",
            "postgwas pipeline",
            (
                "--modules single_cell",
                "--tools scdrs",
                "--vcf study.vcf.gz",
                "--scdrs-h5ad-file brain_atlas.h5ad",
                "--scdrs-gene-id-map entrez_to_symbol.tsv",
                "--scdrs-group-analysis cell_type",
                "--magma-ld-reference reference/g1000_eur",
                "--gene-location-file reference/NCBI37.3.gene.loc",
                "--dataset-id STUDY",
                "--output-directory results",
            ),
        ),
    )


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    from postgwas.modules.single_cell.service import run_single_cell_direct

    try:
        result = run_single_cell_direct(args)
    except (
        SingleCellError,
        MagmaCovarError,
        ConfigurationError,
        OSError,
        ValueError,
    ) as exc:
        Console(stderr=True).print(
            "\n[bold red]Single-cell analysis failed.[/bold red] %s\n" % exc
        )
        return 1
    print("\nSingle-cell analysis completed (%d output artifacts).\n" % len(result.artifacts))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "build_parser", "get_single_cell_parser",
    "get_single_cell_pipeline_examples", "main",
]
