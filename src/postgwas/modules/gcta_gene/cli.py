"""Command-line interface for GCTA fastBAT and mBAT-combo gene tests."""

from __future__ import annotations

import argparse
from typing import get_args

from rich.console import Console

from postgwas.cli.common import get_common_out_parser
from postgwas.cli.compute import get_compute_parser
from postgwas.config import load_configuration
from postgwas.config.models.modules.gcta_gene import (
    GctaChromosomeLabelPolicy,
    GctaDuplicateGenePolicy,
    GctaEmptyPathwayPolicy,
    GctaGeneMethod,
    GctaOversizedSetPolicy,
    GctaUnmappedGenePolicy,
)
from postgwas.core.errors import ConfigurationError
from postgwas.core.execution.runtime import validate_path
from postgwas.core.ui import (
    AlignedRichHelpFormatter,
    format_cli_examples,
    help_with_default,
)
from postgwas.modules.gcta_gene.errors import GctaGeneError


def _with_configured_default(description: str, value) -> str:
    return help_with_default(
        description,
        value,
        label="Configured YAML default",
    )


def get_gcta_gene_parser(add_help=False, *, direct_controls=False):
    defaults = load_configuration()
    module = defaults.modules.gcta_gene
    parser = argparse.ArgumentParser(add_help=add_help)
    inputs = parser.add_argument_group("GCTA association inputs")
    if direct_controls:
        inputs.add_argument(
            "--vcf",
            type=validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True),
            metavar="PATH",
            default=argparse.SUPPRESS,
            help=(
                "Harmonised GWAS-VCF used to reconcile formatter variant IDs to "
                "the exact PLINK BIM IDs. Pipeline mode supplies this automatically"
            ),
        )
        inputs.add_argument(
            "--gcta-input-file",
            type=validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True),
            metavar="PATH",
            default=argparse.SUPPRESS,
            help=_with_configured_default(
                "Formatter-created input for the selected method. In pipeline mode "
                "this is supplied automatically",
                module.input_file,
            ),
        )
    inputs.add_argument(
        "--gcta-reference-prefix",
        metavar="PREFIX",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Prefix shared by the PLINK BED, BIM, and FAM LD-reference files",
            module.reference.prefix,
        ),
    )
    inputs.add_argument(
        "--gene-list",
        type=validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True),
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Four-column chromosome, start, end, gene-ID file used by GCTA",
            module.gene_annotation.file,
        ),
    )
    inputs.add_argument(
        "--fastbat-set-list",
        type=validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True),
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Prepared custom fastBAT set file containing a set ID, SNP IDs, and "
            "END. For fastbat_set, use this or --gmt, not both",
            module.set_annotation.file,
        ),
    )
    inputs.add_argument(
        "--gmt",
        type=validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True),
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "GMT pathway file converted to a fastBAT set list using --gene-list "
            "and the reference BIM. Valid only for fastbat_set and mutually "
            "exclusive with --fastbat-set-list",
            module.set_annotation.gmt_file,
        ),
    )
    compatibility = parser.add_argument_group("Scientific compatibility")
    compatibility.add_argument(
        "--genome-build",
        choices=list(defaults.resources.genomes),
        metavar="BUILD",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Shared genome build of the harmonised GWAS, PLINK LD reference, "
            "and gene list when used; all coordinate-bearing inputs must match",
            module.genome_build,
        ),
    )
    compatibility.add_argument(
        "--gcta-reference-population",
        choices=list(defaults.resources.populations),
        metavar="CODE",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Population represented by the LD reference; no ancestry is inferred",
            module.reference.population,
        ),
    )
    settings = parser.add_argument_group("GCTA association settings")
    settings.add_argument(
        "--method",
        dest="gcta_gene_method",
        choices=get_args(GctaGeneMethod),
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Run fastBAT by genes, fixed segments, or custom SNP sets, or run "
            "gene-based mBAT-combo",
            module.method,
        ),
    )
    settings.add_argument(
        "--gene-window-kb", type=int, metavar="KB", default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Window added to both sides of each gene boundary, including GMT conversion",
            module.gene_window_kb,
        ),
    )
    settings.add_argument(
        "--gmt-chromosome-label-policy",
        choices=get_args(GctaChromosomeLabelPolicy),
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "How GMT conversion reconciles gene-list and BIM chromosome labels",
            module.set_annotation.conversion.chromosome_label_policy,
        ),
    )
    settings.add_argument(
        "--gmt-duplicate-gene-policy",
        choices=get_args(GctaDuplicateGenePolicy),
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "How repeated gene symbols within a GMT pathway are handled",
            module.set_annotation.conversion.duplicate_gene_policy,
        ),
    )
    settings.add_argument(
        "--gmt-unmapped-gene-policy",
        choices=get_args(GctaUnmappedGenePolicy),
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Whether GMT genes absent from the coordinate file fail or are reported",
            module.set_annotation.conversion.unmapped_gene_policy,
        ),
    )
    settings.add_argument(
        "--gmt-empty-pathway-policy",
        choices=get_args(GctaEmptyPathwayPolicy),
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Whether pathways with no usable variants after BIM/GWAS "
            "intersection fail or are omitted",
            module.set_annotation.conversion.empty_pathway_policy,
        ),
    )
    settings.add_argument(
        "--fastbat-oversized-set-policy",
        choices=get_args(GctaOversizedSetPolicy),
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "How sets above GCTA's 20,000-variant hard limit are handled",
            module.set_annotation.oversized_set_policy,
        ),
    )
    settings.add_argument(
        "--fastbat-segment-size-kb", type=int, metavar="KB",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Fixed segment size for the fastbat_segment method",
            module.segment_size_kb,
        ),
    )
    settings.add_argument(
        "--gcta-coordinate-fallback",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Resolve nonmatching summary-statistic IDs by an unambiguous "
            "chromosome, position, and allele-pair match to the PLINK BIM",
            module.variant_harmonisation.coordinate_fallback,
        ),
    )
    settings.add_argument(
        "--gcta-minimum-reference-overlap",
        type=float,
        metavar="FRACTION",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Minimum fraction of unique input variants that must resolve to "
            "exact PLINK BIM IDs",
            module.variant_harmonisation.minimum_overlap_fraction,
        ),
    )
    settings.add_argument(
        "--gcta-chromosome-label-policy",
        choices=get_args(GctaChromosomeLabelPolicy),
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "How variant-ID reconciliation compares chromosome labels in the "
            "harmonised VCF and PLINK BIM",
            module.variant_harmonisation.chromosome_label_policy,
        ),
    )
    settings.add_argument(
        "--gcta-allow-strand-complement",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Allow complementary single-nucleotide allele pairs when only "
            "direct ID matching is available",
            module.variant_harmonisation.allow_strand_complement,
        ),
    )
    settings.add_argument(
        "--gcta-reference-maf-min", type=float, metavar="FRACTION",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Minimum MAF retained from the LD reference",
            module.reference_maf_min,
        ),
    )
    settings.add_argument(
        "--fastbat-ld-cutoff", type=float, metavar="R2", default=argparse.SUPPRESS,
        help=_with_configured_default(
            "LD r-squared pruning cutoff used by fastBAT",
            module.fastbat_ld_cutoff,
        ),
    )
    settings.add_argument(
        "--mbat-svd-gamma", type=float, metavar="FRACTION",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Minimum LD variance retained by mBAT principal components",
            module.mbat_svd_gamma,
        ),
    )
    settings.add_argument(
        "--frequency-difference-max", type=float, metavar="FRACTION",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Maximum A1-frequency difference between GWAS and LD reference",
            module.frequency_difference_max,
        ),
    )
    settings.add_argument(
        "--print-component-p-values",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "For mBAT-combo, include mBAT and fastBAT component p-values",
            module.print_component_p_values,
        ),
    )
    settings.add_argument(
        "--write-snpset",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Write the SNP sets included in the selected fastBAT or mBAT analysis",
            module.write_snpset,
        ),
    )
    settings.add_argument(
        "--gcta-top-results", type=int, metavar="N", default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Number of lowest-p-value associations shown in the terminal summary",
            module.reporting.top_result_count,
        ),
    )
    settings.add_argument(
        "--gcta-nominal-alpha",
        type=float,
        metavar="ALPHA",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Nominal p-value threshold reported for individual associations",
            module.reporting.nominal_alpha,
        ),
    )
    settings.add_argument(
        "--gcta-reporting-alpha",
        type=float,
        metavar="ALPHA",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Family-wise alpha used only to report the Bonferroni threshold; "
            "results are never filtered by this value",
            module.reporting.familywise_alpha,
        ),
    )
    settings.add_argument(
        "--gcta-fdr-alpha",
        type=float,
        metavar="ALPHA",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Benjamini-Hochberg false-discovery-rate threshold",
            module.reporting.fdr_alpha,
        ),
    )
    settings.add_argument(
        "--gcta-p-value-digits",
        type=int,
        metavar="N",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Significant digits used to display association p-values",
            module.reporting.p_value_significant_digits,
        ),
    )
    settings.add_argument(
        "--gcta", metavar="PATH", default=argparse.SUPPRESS,
        help=_with_configured_default(
            "GCTA executable path or command name",
            defaults.resources.executables.gcta,
        ),
    )
    if direct_controls:
        settings.add_argument(
            "--bcftools", metavar="PATH", default=argparse.SUPPRESS,
            help=_with_configured_default(
                "bcftools executable used to read variant coordinates from --vcf",
                defaults.resources.executables.bcftools,
            ),
        )
        settings.add_argument(
            "--run-config", metavar="PATH", default=argparse.SUPPRESS,
            help="YAML settings; explicit command-line values override matching keys.",
        )
        settings.add_argument(
            "--resume",
            action=argparse.BooleanOptionalAction,
            default=argparse.SUPPRESS,
            help=_with_configured_default(
                "Reuse an existing validated primary result",
                defaults.run.resume,
            ),
        )
        settings.add_argument(
            "--overwrite", action="store_true", default=argparse.SUPPRESS,
            help=_with_configured_default(
                "Replace completed or isolated partial outputs",
                defaults.run.overwrite,
            ),
        )
        settings.add_argument(
            "--dry-run", action="store_true", default=argparse.SUPPRESS,
            help="Validate inputs and record the exact command without running analysis.",
        )
    return parser


def get_gcta_gene_pipeline_examples():
    """Return scientifically distinct GCTA examples for pipeline help."""
    return (
        (
            "Run gene-based mBAT-combo:",
            "postgwas pipeline",
            (
                "--modules gcta_gene",
                "--method mbat_combo",
                "--vcf study.vcf.gz",
                "--gcta-reference-prefix reference/1000G_EUR",
                "--gene-list genes_grch37.txt",
                "--genome-build GRCh37",
                "--gcta-reference-population EUR",
                "--dataset-id STUDY",
                "--output-directory results",
            ),
        ),
        (
            "Convert a GMT pathway file and run custom-set fastBAT:",
            "postgwas pipeline",
            (
                "--modules gcta_gene",
                "--method fastbat_set",
                "--vcf study.vcf.gz",
                "--gcta-reference-prefix reference/1000G_EUR",
                "--gene-list genes_grch37.txt",
                "--gmt pathways.gmt",
                "--genome-build GRCh37",
                "--gcta-reference-population EUR",
                "--dataset-id STUDY",
                "--output-directory results",
            ),
        ),
    )


def build_parser():
    defaults = load_configuration()
    parser = argparse.ArgumentParser(
        prog="postgwas gcta_gene",
        usage="postgwas gcta_gene --gcta-input-file PATH --method METHOD [options]",
        description=(
            "Run GCTA fastBAT gene, segment, or custom-set association, or "
            "gene-based mBAT-combo, using a "
            "genome-build- and population-matched PLINK LD reference."
        ),
        epilog=format_cli_examples(
            (
                "Run standalone fastBAT:",
                "postgwas gcta_gene",
                (
                    "--method fastbat_gene",
                    "--gcta-input-file formatted/STUDY_gcta.ma",
                    "--gcta-reference-prefix reference/1000G_EUR",
                    "--gene-list genes_grch37.txt",
                    "--genome-build GRCh37",
                    "--gcta-reference-population EUR",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Export reusable module settings:",
                "postgwas config export",
                ("--module gcta_gene", "--style full", "--output gcta_gene.yaml"),
            ),
            (
                "Convert a GMT file and run custom-set fastBAT:",
                "postgwas gcta_gene",
                (
                    "--method fastbat_set",
                    "--gcta-input-file formatted/STUDY_gcta.ma",
                    "--gcta-reference-prefix reference/1000G_EUR",
                    "--gene-list genes_grch37.txt",
                    "--gmt pathways.gmt",
                    "--genome-build GRCh37",
                    "--gcta-reference-population EUR",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
        ),
        formatter_class=AlignedRichHelpFormatter,
        parents=[
            get_compute_parser(),
            get_common_out_parser(),
            get_gcta_gene_parser(direct_controls=True),
        ],
    )
    for action in parser._actions:
        if action.dest == "dataset_id":
            action.help = _with_configured_default(
                "Short, unique name used in output filenames",
                defaults.run.dataset_id,
            )
            action.default = argparse.SUPPRESS
        elif action.dest == "output_directory":
            action.help = _with_configured_default(
                "Folder where PostGWAS saves results and logs",
                defaults.run.output_directory,
            )
            action.default = argparse.SUPPRESS
    return parser


def main(argv=None):
    from postgwas.modules.gcta_gene.service import run_gcta_gene_direct

    args = build_parser().parse_args(argv)
    try:
        run_gcta_gene_direct(args)
    except (GctaGeneError, ConfigurationError, OSError, ValueError) as exc:
        Console(stderr=True).print(
            "\n[bold red]GCTA gene analysis failed.[/bold red] %s\n" % exc
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "build_parser",
    "get_gcta_gene_parser",
    "get_gcta_gene_pipeline_examples",
    "main",
]
