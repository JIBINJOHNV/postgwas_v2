"""Command-line interface for MAGMA gene and competitive gene-set analysis."""

from __future__ import annotations

import argparse

from postgwas.cli.common import get_common_out_parser
from postgwas.cli.compute import get_compute_parser, positive_float
from postgwas.config import load_configuration
from postgwas.config.models.modules.magma import (
    MAGMA_ALTERNATE_GENE_ID_DUPLICATE_POLICIES,
    MAGMA_DUPLICATE_POLICIES,
    MAGMA_GENE_SET_MISMATCH_ACTIONS,
    MAGMA_GENE_IDENTIFIER_TYPES,
    MAGMA_MHC_POLICIES,
)
from postgwas.core.errors import ConfigurationError
from postgwas.core.ui import (
    print_screen_message,
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
    bim_identifier_column = module.input.bim_columns.index("variant_id") + 1
    alternate_gene_column = (
        module.input.gene_location_columns.index("alternate_gene_id") + 1
        if "alternate_gene_id" in module.input.gene_location_columns else None
    )
    parser = argparse.ArgumentParser(add_help=add_help)

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
    positional = parser.add_argument_group(
        "MAGMA positional gene-location declaration"
    )
    positional_definition = module.mapping.definitions["positional"]
    positional.add_argument(
        "--magma-positional-gene-id-type",
        choices=MAGMA_GENE_IDENTIFIER_TYPES,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Identifier system stored in column 1 of the positional MAGMA "
            "gene-location file",
            positional_definition.gene_id_type,
        ),
    )
    positional.add_argument(
        "--magma-positional-source-name",
        metavar="TEXT",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Name of the source that supplied the positional gene-location file",
            positional_definition.source_name,
        ),
    )
    positional.add_argument(
        "--magma-positional-source-version",
        metavar="TEXT",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Version of the positional gene-location source",
            positional_definition.source_version,
        ),
    )
    positional.add_argument(
        "--magma-positional-source-url",
        metavar="URL",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Provenance URL for the positional gene-location source",
            positional_definition.source_url,
        ),
    )
    positional.add_argument(
        "--magma-positional-context",
        metavar="TEXT",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Scientific context represented by the positional gene-location file",
            positional_definition.context,
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
            "BIM field-%d matches when --resolve-variants-to-reference is enabled"
            % bim_identifier_column,
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
        "--gene-set-identifier-mismatch",
        choices=MAGMA_GENE_SET_MISMATCH_ACTIONS,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Action when gene-set identifiers are incompatible with the selected "
            "gene mapping: skip preserves gene-association analysis and omits "
            "competitive pathway analysis; error stops the complete MAGMA run",
            module.gene_sets.identifier_mismatch_action,
        ),
    )
    settings.add_argument(
        "--alternate-gene-id-duplicate-policy",
        choices=MAGMA_ALTERNATE_GENE_ID_DUPLICATE_POLICIES,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Handling when gene-location column %s is selected for positional "
            "MAGMA but one alternate identifier labels multiple intervals: "
            "longest_interval retains the longest interval and error stops "
            "pathway-compatible reference preparation"
            % (alternate_gene_column or "alternate"),
            module.gene_sets.alternate_id_duplicate_policy,
        ),
    )
    settings.add_argument(
        "--gene-location-alternate-id-type",
        choices=MAGMA_GENE_IDENTIFIER_TYPES,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Identifier system stored in the optional alternate gene-location "
            "column and used by a derived positional MAGMA reference when it "
            "matches the pathway file",
            module.input.alternate_gene_id_type,
        ),
    )
    settings.add_argument(
        "--duplicate-policy",
        choices=MAGMA_DUPLICATE_POLICIES,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Handling for repeated SNP IDs after optional LD-reference filtering: "
            "err stops, lowest_p retains one row with the lowest valid p-value, "
            "and remove excludes every row in each duplicated-ID group",
            module.snp_harmonisation.duplicate_policy,
        ),
    )
    settings.add_argument(
        "--resolve-variants-to-reference",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Compare formatter-created variant IDs with the configured BIM "
            "variant-ID field and retain only exact identifier matches. When "
            "disabled, use the formatter files directly without reference filtering",
            module.snp_harmonisation.resolve_variants_to_reference,
        ),
    )
    settings.add_argument(
        "--magma-memory-per-worker-gb",
        type=positive_float,
        default=argparse.SUPPRESS,
        metavar="GB",
        help=help_with_default(
            "Memory budget assigned to each concurrent MAGMA gene-association "
            "worker. PostGWAS limits the worker count using this value and the "
            "total --memory-gb budget; it does not impose a memory limit on MAGMA",
            "%g GB" % module.batching.memory_per_process_gb,
        ),
    )
    exclusions = parser.add_argument_group("MAGMA analysis scope")
    exclusions.add_argument(
        "--mhc-policy",
        choices=MAGMA_MHC_POLICIES,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "MHC handling: include keeps all MHC SNPs and genes; exclude_snps "
            "removes MHC SNPs before annotation and gene testing; exclude_genes "
            "removes overlapping annotated units from gene and competitive "
            "gene-set analysis; exclude_both applies both exclusions",
            module.mhc.policy,
        ),
    )
    genome_region = defaults.resources.genomes[module.genome_build.value].regions["mhc"]
    exclusions.add_argument(
        "--mhc-chrom",
        metavar="CHROM",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Override the MHC chromosome. Supply all three MHC override options "
            "together",
            genome_region.chromosome,
        ),
    )
    exclusions.add_argument(
        "--mhc-start",
        metavar="POSITION",
        type=int,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Override the one-based inclusive MHC start position. Supply all "
            "three MHC override options together",
            genome_region.start,
        ),
    )
    exclusions.add_argument(
        "--mhc-end",
        metavar="POSITION",
        type=int,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Override the one-based inclusive MHC end position. Supply all three "
            "MHC override options together",
            genome_region.end,
        ),
    )
    exclusions.add_argument(
        "--exclude-chromosomes",
        metavar="CHROM",
        nargs="+",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Chromosomes excluded from SNP, gene, and competitive gene-set "
            "analysis. X is retained unless explicitly listed",
            " ".join(module.chromosomes.exclude),
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
        print_screen_message(
            "error", "MAGMA analysis failed. %s" % exc, stderr=True,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "build_parser", "get_magma_parser", "get_magma_pipeline_examples", "main",
]
