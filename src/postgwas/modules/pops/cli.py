#!/usr/bin/env python3
"""Command-line interface for the upstream PoPS v0.2 integration."""

import argparse
from rich.console import Console

from postgwas.cli.compute import get_compute_parser
from postgwas.config import load_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.core.ui import (
    AlignedRichHelpFormatter,
    format_cli_examples,
    help_with_default,
    mark_cli_required_help,
)
from postgwas.modules.pops.errors import PopsError

# =============================
# Allowed Common CLI Builders
# =============================
from postgwas.cli.common import (
    get_common_out_parser,
    get_common_pops_parser,
)



# =====================================================================
#   DIRECT MODE PARSER
# =====================================================================
def get_pops_direct_parser(add_help=False):
    """Return standalone-only inputs and configuration controls."""
    defaults = load_configuration()
    module = defaults.modules.pops
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group(
        "PoPS input",
        "Direct analysis requires exactly one target source: "
        "--magma-association-prefix or --target-score-file.",
    )
    group.add_argument(
        "--magma-association-prefix",
        metavar="PREFIX",
        default=argparse.SUPPRESS,
        help=(
            "Prefix shared by MAGMA .genes.out and .genes.raw files. Mutually "
            "exclusive with --target-score-file."
        ),
    )
    group.add_argument(
        "--genome-build",
        choices=list(defaults.resources.genomes),
        metavar="BUILD",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Declared genome build shared by MAGMA and the PoPS gene annotation; "
            "PostGWAS does not infer it",
            module.genome_build,
        ),
    )
    controls = parser.add_argument_group("Configuration")
    controls.add_argument(
        "--run-config",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help="YAML settings; explicit command-line values override matching keys.",
    )
    group.add_argument(
        "--verbose",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Enable verbose messages from upstream PoPS",
            module.verbose,
        ),
    )
    return parser


def get_pops_pipeline_parser(add_help=False):
    """Return the conflict-free PoPS genome-build override for pipeline mode."""
    defaults = load_configuration()
    module = defaults.modules.pops
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("PoPS input compatibility")
    group.add_argument(
        "--pops-genome-build",
        dest="pops_genome_build",
        choices=list(defaults.resources.genomes),
        metavar="BUILD",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Genome build shared by MAGMA and PoPS resources",
            module.genome_build,
        ),
    )
    return parser


def get_pops_pipeline_examples():
    """Return a runnable MAGMA-to-PoPS pipeline example."""
    return (
        (
            "Run MAGMA gene analysis followed by PoPS prioritisation:",
            "postgwas pipeline",
            (
                "--modules pops",
                "--vcf study.vcf.gz",
                "--magma-ld-reference reference/1000G_EUR",
                "--gene-location-file reference/FUMA/ENSGv102.coding.genes.txt",
                "--feature-matrix-prefix "
                "reference/FLAMES/pops_features_full_FUMA_compatible/"
                "features_munged/pops_features",
                "--feature-matrix-chunks 116",
                "--pops-gene-location-file "
                "reference/FLAMES/pops_features_full_FUMA_compatible/"
                "gene_annots.txt",
                "--control-features-file "
                "reference/FLAMES/pops_features_full_FUMA_compatible/"
                "control.features",
                "--pops-genome-build GRCh37",
                "--dataset-id STUDY",
                "--output-directory results",
            ),
        ),
        (
            "Run MAGMA and explicitly intersect incompatible PoPS target genes:",
            "postgwas pipeline",
            (
                "--modules pops",
                "--vcf study.vcf.gz",
                "--magma-ld-reference reference/1000G_EUR",
                "--gene-location-file reference/FUMA/ENSGv102.coding.genes.txt",
                "--feature-matrix-prefix "
                "reference/FLAMES/pops_features_full_FUMA_compatible/"
                "features_munged/pops_features",
                "--feature-matrix-chunks 116",
                "--pops-gene-location-file "
                "reference/FLAMES/pops_features_full_FUMA_compatible/"
                "gene_annots.txt",
                "--control-features-file "
                "reference/FLAMES/pops_features_full_FUMA_compatible/"
                "control.features",
                "--pops-genome-build GRCh37",
                "--gene-universe-policy intersect",
                "--dataset-id STUDY",
                "--output-directory results",
            ),
        ),
    )



# =====================================================================
#   MAIN CLI ENTRY
# =====================================================================
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="postgwas pops",
        usage=(
            "postgwas pops (--magma-association-prefix PREFIX | "
            "--target-score-file PATH) [options]"
        ),
        description="Prioritise genes with PoPS using existing MAGMA results and feature matrices.",
        epilog=format_cli_examples(
            (
                "Run PoPS gene prioritisation:",
                "postgwas pops",
                (
                    "--magma-association-prefix magma/STUDY",
                    "--feature-matrix-prefix "
                    "reference/pops/features_munged/pops_features",
                    "--feature-matrix-chunks 2",
                    "--pops-gene-location-file reference/pops/gene_annot.tsv",
                    "--control-features-file reference/pops/control.features",
                    "--genome-build GRCh37",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Explicitly intersect incompatible MAGMA target genes:",
                "postgwas pops",
                (
                    "--magma-association-prefix magma/STUDY",
                    "--feature-matrix-prefix "
                    "reference/pops/features_munged/pops_features",
                    "--feature-matrix-chunks 2",
                    "--pops-gene-location-file reference/pops/gene_annot.tsv",
                    "--gene-universe-policy intersect",
                    "--genome-build GRCh37",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Run PoPS with a custom target-score table:",
                "postgwas pops",
                (
                    "--target-score-file target_scores.tsv",
                    "--feature-matrix-prefix "
                    "reference/pops/features_munged/pops_features",
                    "--feature-matrix-chunks 2",
                    "--pops-gene-location-file reference/pops/gene_annot.tsv",
                    "--genome-build GRCh37",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Export reusable PoPS settings:",
                "postgwas config export",
                ("--module pops", "--style full", "--output pops.yaml"),
            ),
            (
                "Run from exported settings:",
                "postgwas pops",
                ("--run-config pops.yaml",),
            ),
            notes=(
                "PoPS requires MAGMA .genes.out and .genes.raw, a matching "
                "Ensembl gene annotation, and every declared feature-matrix chunk.",
                "PoPS does not use MAGMACOVAR .gsa.out results. "
                "--use-magma-covariates refers to technical covariates read from "
                "MAGMA .genes.raw.",
            ),
        ),
        formatter_class=AlignedRichHelpFormatter,
        parents=[
            get_compute_parser(),
            get_common_out_parser(),
            get_pops_direct_parser(add_help=False),
            get_common_pops_parser(add_help=False),
        ],
    )
    for action in parser._actions:
        if action.dest in {
            "dataset_id", "output_directory", "threads", "memory_gb", "seed",
        }:
            action.default = argparse.SUPPRESS
            action.type = None
    mark_cli_required_help(parser, ("genome_build",))
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    from postgwas.modules.pops.service import run_pops_direct

    try:
        run_pops_direct(args)
    except (PopsError, ConfigurationError, OSError, ValueError) as exc:
        Console(stderr=True).print(
            "\n[bold red]PoPS analysis failed.[/bold red] %s\n" % exc
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "build_parser",
    "get_pops_direct_parser",
    "get_pops_pipeline_examples",
    "get_pops_pipeline_parser",
    "main",
]
