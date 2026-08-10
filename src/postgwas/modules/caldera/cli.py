"""Command-line interface for CALDERA."""

from __future__ import annotations

import argparse

from rich.console import Console

from postgwas.cli.common import get_common_out_parser
from postgwas.cli.compute import get_compute_parser
from postgwas.config import load_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.core.ui import AlignedRichHelpFormatter, format_cli_examples, help_with_default
from postgwas.modules.caldera.errors import CalderaError


def get_caldera_pipeline_examples():
    """Return a runnable GRCh37 pipeline through PoPS, fine-mapping and CALDERA."""
    return (
        (
            "Run the complete GRCh37 PoPS, fine-mapping, and CALDERA workflow:",
            "postgwas pipeline",
            (
                "--modules caldera",
                "--vcf study_GRCh37.vcf.gz",
                "--genome-build GRCh37",
                "--ld-region-dir reference/ld_blocks",
                "--ld-block-populations EUR",
                "--ld-folder reference/ld",
                "--population EUR",
                "--magma-ld-reference reference/1000G_EUR",
                "--gene-location-file reference/NCBI37.3.gene.loc",
                "--feature-matrix-prefix reference/pops/features_munged/pops_features",
                "--feature-matrix-chunks 116",
                "--pops-gene-location-file reference/pops/GRCh37_gene_annot.tsv",
                "--pops-genome-build GRCh37",
                "--finemap-method susie",
                "--finemap-ld-reference reference/1000G_EUR",
                "--caldera-genome-build GRCh37",
                "--dataset-id STUDY",
                "--output-directory results",
            ),
        ),
    )


def get_caldera_parser(add_help=False, *, direct_controls=False):
    defaults = load_configuration()
    module = defaults.modules.caldera
    parser = argparse.ArgumentParser(add_help=add_help)
    inputs = parser.add_argument_group("CALDERA inputs")
    inputs.add_argument("--pops-file", metavar="PATH", default=argparse.SUPPRESS)
    inputs.add_argument("--credible-set-file", metavar="PATH", default=argparse.SUPPRESS)
    inputs.add_argument("--caldera-repository", metavar="PATH", default=argparse.SUPPRESS)
    inputs.add_argument("--caldera-adapter-script", metavar="PATH", default=argparse.SUPPRESS)
    inputs.add_argument(
        "--caldera-genome-build", choices=list(defaults.resources.genomes),
        metavar="BUILD", default=argparse.SUPPRESS,
        help=help_with_default(
            "Genome build of credible-set coordinates and CALDERA gene locations",
            module.genome_build.value, label="Configured default",
        ),
    )
    if direct_controls:
        controls = parser.add_argument_group("Configuration")
        controls.add_argument("--run-config", metavar="PATH", default=argparse.SUPPRESS)
        controls.add_argument(
            "--resume", action=argparse.BooleanOptionalAction,
            default=argparse.SUPPRESS,
        )
        controls.add_argument("--overwrite", action="store_true", default=argparse.SUPPRESS)
    return parser


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="postgwas caldera",
        usage="postgwas caldera --pops-file PATH --credible-set-file PATH [options]",
        description="Prioritise causal genes with CALDERA from PoPS and credible sets.",
        formatter_class=AlignedRichHelpFormatter,
        epilog=format_cli_examples(
            (
                "Run CALDERA directly:", "postgwas caldera",
                (
                    "--pops-file pops/STUDY_pops.preds",
                    "--credible-set-file finemap/STUDY_credible_sets_GRCh37.tsv",
                    "--caldera-genome-build GRCh37",
                    "--dataset-id STUDY", "--output-directory results",
                ),
            ),
            *get_caldera_pipeline_examples(),
            notes=(
                "Direct credible-set tables require locus, chr, bp, and pip columns with at least 0.95 cumulative PIP per locus.",
                "Pipeline mode supplies PoPS predictions and credible sets internally and currently supports GRCh37 only.",
            ),
        ),
        parents=[
            get_compute_parser(), get_common_out_parser(),
            get_caldera_parser(direct_controls=True),
        ],
    )
    for action in parser._actions:
        if action.dest in {"dataset_id", "output_directory", "threads", "memory_gb", "seed"}:
            action.default = argparse.SUPPRESS
            action.type = None
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    from postgwas.modules.caldera.service import run_caldera_direct

    try:
        run_caldera_direct(args)
    except (CalderaError, ConfigurationError, OSError, ValueError) as exc:
        Console(stderr=True).print("\n[bold red]CALDERA analysis failed.[/bold red] %s\n" % exc)
        return 1
    return 0


__all__ = [
    "build_parser", "get_caldera_parser", "get_caldera_pipeline_examples", "main",
]
