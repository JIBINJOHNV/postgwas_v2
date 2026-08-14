"""Command-line interface for MAGMA gene and competitive gene-set analysis."""

from __future__ import annotations

import argparse

from rich.console import Console

from postgwas.cli.common import get_common_out_parser, get_magma_binary_parser
from postgwas.cli.compute import get_compute_parser
from postgwas.config import load_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.core.ui import (
    AlignedRichHelpFormatter,
    format_cli_examples,
    help_with_default,
    mark_cli_required_help,
)
from postgwas.modules.magma.errors import MagmaError


def get_magma_parser(add_help=False, *, direct_controls=False):
    """Return reusable MAGMA options without assigning CLI-owned defaults."""
    defaults = load_configuration()
    module = defaults.modules.magma
    parser = argparse.ArgumentParser(
        add_help=add_help,
        parents=[get_magma_binary_parser()],
    )

    inputs = parser.add_argument_group("MAGMA inputs")
    inputs.add_argument(
        "--snp-location-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help="Formatter-created SNP location table consumed by direct MAGMA analysis.",
    )
    inputs.add_argument(
        "--p-value-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "Formatter-created SNP, raw p-value and per-variant sample-size table. "
            "The SNP identifiers must correspond to --snp-location-file."
        ),
    )
    inputs.add_argument(
        "--magma-ld-reference",
        metavar="PREFIX",
        default=argparse.SUPPRESS,
        help=(
            "LD-reference prefix with configured companion files: %s."
            % ", ".join(module.input.required_reference_extensions)
        ),
    )
    inputs.add_argument(
        "--gene-location-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help="MAGMA gene-location file used to map reference SNPs to genes.",
    )
    inputs.add_argument(
        "--gene-set-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "Optional standard GMT or native MAGMA set-annotation file for "
            "competitive gene-set analysis."
        ),
    )
    settings = parser.add_argument_group("MAGMA analysis settings")
    settings.add_argument(
        "--magma-mapping",
        metavar="NAME",
        nargs="+",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "One or more mapping definitions from modules.magma.mapping.definitions; "
            "each is analysed as a separate hypothesis family. The generated "
            "functional-resource YAML lists every available tissue and cell context",
            " ".join(module.mapping.selected),
        ),
    )
    settings.add_argument(
        "--primary-magma-mapping",
        metavar="NAME",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Mapping designated as the primary MAGMA result. Use positional, "
            "eMAGMA, H-MAGMA, or nMAGMA rather than chromMAGMA when subsequent "
            "analyses require ordinary MAGMA gene units",
            module.mapping.primary,
        ),
    )
    settings.add_argument(
        "--window-upstream",
        metavar="KB",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Kilobases added upstream of each gene",
            module.gene_window_upstream_kb,
        ),
    )
    settings.add_argument(
        "--window-downstream",
        metavar="KB",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Kilobases added downstream of each gene",
            module.gene_window_downstream_kb,
        ),
    )
    settings.add_argument(
        "--gene-model",
        metavar="MODEL",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "MAGMA model supported with SNP p-value input, for example "
            "snp-wise=mean, snp-wise=top or multi=snp-wise",
            module.gene_model,
        ),
    )
    settings.add_argument(
        "--sample-size-column",
        metavar="COLUMN",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Per-variant total-sample-size column in the p-value table",
            module.input.sample_size_column,
        ),
    )
    settings.add_argument(
        "--minimum-snp-overlap",
        metavar="FRACTION",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Minimum fraction of unique formatter variants that must have exact "
            "BIM field-2 matches when --resolve-variants-to-reference is enabled",
            module.snp_harmonisation.minimum_overlap_fraction,
        ),
    )
    settings.add_argument(
        "--minimum-gene-id-overlap",
        metavar="FRACTION",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Minimum gene-ID overlap required between tested MAGMA units and a "
            "gene-set file unless a mapping definition supplies its own validated floor",
            module.gene_sets.minimum_gene_id_overlap_fraction,
        ),
    )
    settings.add_argument(
        "--resolve-variants-to-reference",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Compare formatter-created variant IDs with BIM field 2 and retain "
            "only exact coordinate-and-allele-consistent reference matches. When "
            "disabled, use the formatter files directly without reference filtering",
            module.snp_harmonisation.resolve_variants_to_reference,
        ),
    )
    if direct_controls:
        settings.add_argument(
            "--run-config",
            metavar="PATH",
            default=argparse.SUPPRESS,
            help="YAML settings; explicit command-line values override matching keys.",
        )
    return parser


def build_parser() -> argparse.ArgumentParser:
    defaults = load_configuration()
    formatting = defaults.modules.formatting.exports["magma"].outputs
    location_name = formatting["snp_location"].output_file.format(dataset_id="STUDY")
    p_value_name = formatting["p_values"].output_file.format(dataset_id="STUDY")
    parser = argparse.ArgumentParser(
        prog="postgwas magma",
        usage=(
            "postgwas magma --snp-location-file PATH --p-value-file PATH "
            "--magma-ld-reference PREFIX --gene-location-file PATH "
            "--dataset-id NAME --output-directory PATH [options]"
        ),
        description=(
            "Run MAGMA gene association and optional competitive gene-set analysis "
            "from PostGWAS formatter inputs."
        ),
        epilog=format_cli_examples(
            (
                "Create MAGMA inputs first:",
                "postgwas formatter",
                (
                    "--vcf study.vcf.gz",
                    "--format magma",
                    "--dataset-id STUDY",
                    "--output-directory formatted",
                ),
            ),
            (
                "Run standalone MAGMA gene analysis:",
                "postgwas magma",
                (
                    "--snp-location-file formatted/%s" % location_name,
                    "--p-value-file formatted/%s" % p_value_name,
                    "--magma-ld-reference reference/1000G_EUR",
                    "--gene-location-file reference/NCBI37.3.gene.loc",
                    "--resolve-variants-to-reference",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Export reusable MAGMA settings:",
                "postgwas config export",
                ("--module magma", "--style full", "--output magma.yaml"),
            ),
            (
                "Run entirely from exported settings:",
                "postgwas magma",
                ("--run-config magma.yaml",),
            ),
            (
                "Run configured positional, eMAGMA, H-MAGMA, nMAGMA and chromMAGMA mappings:",
                "postgwas magma",
                (
                    "--run-config resources/magma/functional_mapping/configs/magma_functional_mapping.yaml",
                    "--magma-mapping positional emagma_brain_amygdala h_magma_adult_brain n_magma_cortex chrom_magma_ovarian_h3k27ac",
                    "--primary-magma-mapping positional",
                    "--snp-location-file formatted/%s" % location_name,
                    "--p-value-file formatted/%s" % p_value_name,
                    "--magma-ld-reference reference/g1000_eur",
                    "--gene-location-file reference/NCBI37.3.gene.loc",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
        ),
        formatter_class=AlignedRichHelpFormatter,
        parents=[
            get_compute_parser(),
            get_common_out_parser(),
            get_magma_parser(direct_controls=True),
        ],
    )
    for action in parser._actions:
        if action.dest in {
            "dataset_id", "output_directory", "threads", "memory_gb", "seed",
        }:
            action.default = argparse.SUPPRESS
            action.type = None
    mark_cli_required_help(
        parser,
        (
            "snp_location_file",
            "p_value_file",
            "magma_ld_reference",
            "gene_location_file",
            "dataset_id",
            "output_directory",
        ),
    )
    return parser


def get_magma_pipeline_examples():
    """Return complete MAGMA examples for context-sensitive pipeline help."""
    return (
        (
            "Run conventional positional MAGMA from a GWAS-VCF:",
            "postgwas pipeline",
            (
                "--modules magma",
                "--vcf study.vcf.gz",
                "--dataset-id STUDY",
                "--output-directory results",
                "--magma-ld-reference reference/g1000_eur",
                "--gene-location-file reference/NCBI37.3.gene.loc",
            ),
        ),
        (
            "Run configured positional, eMAGMA, H-MAGMA, nMAGMA and chromMAGMA mappings:",
            "postgwas pipeline",
            (
                "--modules magma",
                "--vcf study.vcf.gz",
                "--dataset-id STUDY",
                "--output-directory results",
                "--run-config resources/magma/functional_mapping/configs/magma_functional_mapping.yaml",
                "--magma-mapping positional emagma_brain_amygdala h_magma_adult_brain n_magma_cortex chrom_magma_ovarian_h3k27ac",
                "--primary-magma-mapping positional",
                "--magma-ld-reference reference/g1000_eur",
                "--gene-location-file reference/NCBI37.3.gene.loc",
            ),
        ),
    )


def main(argv=None):
    args = build_parser().parse_args(argv)
    from postgwas.modules.magma.service import run_magma_direct

    try:
        run_magma_direct(args)
    except (MagmaError, ConfigurationError, OSError, ValueError) as exc:
        Console(stderr=True).print(
            "\n[bold red]MAGMA analysis failed.[/bold red] %s\n" % exc
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "build_parser", "get_magma_parser", "get_magma_pipeline_examples", "main",
]
