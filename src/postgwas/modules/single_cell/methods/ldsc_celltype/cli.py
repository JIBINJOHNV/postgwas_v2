"""Command-line options for LDSC cell-type-specific analyses."""

from __future__ import annotations

import argparse

from postgwas.core.ui import help_with_default


def add_ldsc_celltype_arguments(parser: argparse.ArgumentParser, defaults) -> None:
    """Add method-owned options without defining a second set of defaults."""
    method = defaults.modules.single_cell.ldsc_celltype
    inputs = parser.add_argument_group("LDSC cell-type inputs")
    inputs.add_argument(
        "--ldsc-celltype-sumstats-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "Pre-munged LDSC .sumstats.gz consumed by direct cell-type analysis."
        ),
    )
    inputs.add_argument(
        "--ldsc-celltype-ldcts-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "Two-column .ldcts manifest: cell-type label and comma-delimited "
            "tested/control LD-score prefixes."
        ),
    )
    inputs.add_argument(
        "--ldsc-celltype-baseline-prefix",
        nargs="+",
        metavar="PREFIX",
        default=argparse.SUPPRESS,
        help="One or more chromosome-aware baseline-model LD-score prefixes.",
    )
    inputs.add_argument(
        "--ldsc-celltype-weights-prefix",
        metavar="PREFIX",
        default=argparse.SUPPRESS,
        help="Chromosome-aware LDSC regression-weights prefix.",
    )
    inputs.add_argument(
        "--ldsc-celltype-merge-alleles-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help="HapMap3 SNP list used when formatter output must be munged.",
    )
    inputs.add_argument(
        "--ldsc",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "LDSC executable path or command name",
            defaults.resources.executables.ldsc,
        ),
    )
    inputs.add_argument(
        "--munge-sumstats",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "LDSC munge_sumstats executable path or command name",
            defaults.resources.executables.munge_sumstats,
        ),
    )
    settings = parser.add_argument_group("LDSC cell-type settings")
    settings.add_argument(
        "--ldsc-celltype-sumstats-source",
        choices=("munged", "formatter"),
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Whether the supplied statistics are already munged or formatter output",
            method.input.source,
        ),
    )
    settings.add_argument(
        "--ldsc-celltype-genome-build",
        choices=("GRCh37", "GRCh38"),
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Declared genome build shared by GWAS and LD-score resources",
            method.genome_build.value,
        ),
    )
    settings.add_argument(
        "--ldsc-celltype-population",
        choices=("EUR", "AFR", "EAS", "SAS", "AMR"),
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Declared ancestry of the LD-score resources",
            method.population.value,
        ),
    )


__all__ = ["add_ldsc_celltype_arguments"]
