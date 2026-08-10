"""Command-line interface for K-POPS v1.0.0."""

from __future__ import annotations

import argparse

from rich.console import Console

from postgwas.cli.common import get_common_out_parser
from postgwas.cli.compute import get_compute_parser
from postgwas.config import load_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.core.ui import AlignedRichHelpFormatter, format_cli_examples, help_with_default
from postgwas.modules.kpops.errors import KPopsError


def get_kpops_pipeline_examples():
    """Return a runnable formatter-to-MAGMA-to-K-POPS pipeline example."""
    return (
        (
            "Run MAGMA gene analysis followed by K-POPS prioritisation:",
            "postgwas pipeline",
            (
                "--modules kpops",
                "--vcf study_GRCh37.vcf.gz",
                "--magma-ld-reference reference/1000G_EUR",
                "--gene-location-file reference/NCBI37.3.gene.loc",
                "--kpops-gene-annotation-file reference/kpops/GRCh37_gene_annot.tsv",
                "--kernel-matrix-prefix reference/kpops/kernel_linear",
                "--kpops-genome-build GRCh37",
                "--dataset-id STUDY",
                "--output-directory results",
            ),
        ),
    )


def get_kpops_parser(add_help=False, *, direct_controls=False):
    defaults = load_configuration()
    module = defaults.modules.kpops
    parser = argparse.ArgumentParser(add_help=add_help)
    inputs = parser.add_argument_group("K-POPS inputs")
    inputs.add_argument("--magma-association-prefix", metavar="PREFIX", default=argparse.SUPPRESS)
    inputs.add_argument("--kpops-gene-annotation-file", metavar="PATH", default=argparse.SUPPRESS)
    inputs.add_argument("--kernel-matrix-prefix", metavar="PREFIX", default=argparse.SUPPRESS)
    inputs.add_argument("--kpops-script", metavar="PATH", default=argparse.SUPPRESS)
    inputs.add_argument(
        "--kpops-genome-build", choices=list(defaults.resources.genomes),
        metavar="BUILD", default=argparse.SUPPRESS,
        help=help_with_default(
            "Genome build shared by MAGMA, the K-POPS annotation, and kernel genes",
            module.genome_build, label="Configured default",
        ),
    )
    settings = parser.add_argument_group("K-POPS model settings")
    settings.add_argument(
        "--training-chromosomes", nargs="+", metavar="CHROM",
        default=argparse.SUPPRESS,
        help="Use 'loco', 'all', or an explicit set of annotation chromosome labels.",
    )
    settings.add_argument(
        "--kpops-device", choices=("cpu", "cuda", "mps"),
        default=argparse.SUPPRESS,
        help=help_with_default("PyTorch execution device", module.device, label="Configured default"),
    )
    settings.add_argument(
        "--top-contributor-gene-count", type=int, metavar="N",
        default=argparse.SUPPRESS,
    )
    settings.add_argument("--anchor-genes", nargs="+", metavar="GENE", default=argparse.SUPPRESS)
    settings.add_argument(
        "--anchor-gene-type", choices=("ENSGID", "NAME"),
        default=argparse.SUPPRESS,
    )
    settings.add_argument(
        "--use-magma-covariates", action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
    )
    settings.add_argument(
        "--save-attribution-files", action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
    )
    settings.add_argument(
        "--kpops-verbose", action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
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
        prog="postgwas kpops",
        usage="postgwas kpops --magma-association-prefix PREFIX [options]",
        description="Prioritise genes with K-POPS from validated MAGMA and kernel inputs.",
        formatter_class=AlignedRichHelpFormatter,
        epilog=format_cli_examples(
            (
                "Run K-POPS directly:", "postgwas kpops",
                (
                    "--magma-association-prefix magma/results/STUDY_magma_35up_10down",
                    "--kpops-gene-annotation-file reference/kpops/GRCh37_gene_annot.tsv",
                    "--kernel-matrix-prefix reference/kpops/kernel_linear",
                    "--kpops-genome-build GRCh37",
                    "--dataset-id STUDY", "--output-directory results",
                ),
            ),
            *get_kpops_pipeline_examples(),
            (
                "Export reusable settings:", "postgwas config export",
                ("--module kpops", "--style full", "--output kpops.yaml"),
            ),
            notes=(
                "Direct mode requires PREFIX.genes.out and PREFIX.genes.raw from the same MAGMA run.",
                "Pipeline mode creates MAGMA results, but still requires compatible MAGMA and K-POPS reference resources.",
            ),
        ),
        parents=[
            get_compute_parser(), get_common_out_parser(),
            get_kpops_parser(direct_controls=True),
        ],
    )
    for action in parser._actions:
        if action.dest in {"dataset_id", "output_directory", "threads", "memory_gb", "seed"}:
            action.default = argparse.SUPPRESS
            action.type = None
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    from postgwas.modules.kpops.service import run_kpops_direct

    try:
        run_kpops_direct(args)
    except (KPopsError, ConfigurationError, OSError, ValueError) as exc:
        Console(stderr=True).print("\n[bold red]K-POPS analysis failed.[/bold red] %s\n" % exc)
        return 1
    return 0


__all__ = [
    "build_parser", "get_kpops_parser", "get_kpops_pipeline_examples", "main",
]
