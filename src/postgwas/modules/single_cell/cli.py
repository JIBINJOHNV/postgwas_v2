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
from postgwas.modules.single_cell.methods.magma_celltype.cli import (
    add_magma_celltype_arguments,
)
from postgwas.modules.single_cell.methods.ldsc_celltype.cli import (
    add_ldsc_celltype_arguments,
)
from postgwas.modules.single_cell.methods.scdrs.cli import add_scdrs_arguments


def get_single_cell_parser(add_help=False, *, direct_controls=False):
    """Return reusable single-cell options without CLI-owned defaults."""
    defaults = load_configuration()
    module = defaults.modules.single_cell
    parser = argparse.ArgumentParser(add_help=add_help)
    settings = parser.add_argument_group("Single-cell method selection")
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
    add_magma_celltype_arguments(parser, defaults)
    add_scdrs_arguments(parser, defaults)
    add_ldsc_celltype_arguments(parser, defaults)
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
            "scDRS uses an exact H5AD atlas and either native .gs gene sets or "
            "a gene set constructed from MAGMA Z statistics; LDSC cell typing "
            "tests cell-type LD-score annotations with --h2-cts."
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
                "Run LDSC cell-type analysis from munged summary statistics:",
                "postgwas single_cell",
                (
                    "--tools ldsc_celltype",
                    "--ldsc-celltype-sumstats-file STUDY.sumstats.gz",
                    "--ldsc-celltype-ldcts-file brain.ldcts",
                    "--ldsc-celltype-baseline-prefix reference/baselineLD.",
                    "--ldsc-celltype-weights-prefix reference/weights.",
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
        if action.dest == "seed":
            action.help = (
                "%s Native scDRS 1.0.3 does not consume this option and uses "
                "its fixed internal seed 0." % action.help
            )
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
        (
            "Run LDSC cell typing from a harmonised GWAS-VCF:",
            "postgwas pipeline",
            (
                "--modules single_cell",
                "--tools ldsc_celltype",
                "--vcf study.vcf.gz",
                "--ldsc-celltype-ldcts-file brain.ldcts",
                "--ldsc-celltype-baseline-prefix reference/baselineLD.",
                "--ldsc-celltype-weights-prefix reference/weights.",
                "--ldsc-celltype-merge-alleles-file reference/w_hm3.snplist",
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
