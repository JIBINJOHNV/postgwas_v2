import argparse
from rich_argparse import RichHelpFormatter
from rich_argparse import RawTextRichHelpFormatter

from postgwas.core.execution.runtime import validate_alphanumeric, validate_path, validate_prefix_files


# get_inputvcf_parser()
# get_genome_build_parser()
# get_common_out_parser()
# get_flames_common_parser(add_help=False)
# get_common_imputation_parser


def add_variant_qc_policy_arguments(
    parser: argparse.ArgumentParser,
    *,
    defaults,
    reference_af_column,
    maximum_af_difference,
    group_title: str,
):
    """Add the policy options shared by QC assessment and VCF filtering."""
    from postgwas.core.ui import help_with_default

    group = parser.add_argument_group(group_title)
    group.add_argument(
        "--minimum-neglog10-p",
        dest="minimum_neglog10_p",
        type=float,
        metavar="VALUE",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Minimum -log10(P) accepted by the QC policy. For P <= 5e-8, "
            "use 7.30103",
            defaults.minimum_neglog10_p,
        ),
    )
    group.add_argument(
        "--minimum-maf",
        dest="minimum_maf",
        type=float,
        metavar="VALUE",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Minimum minor-allele frequency accepted by the QC policy",
            defaults.maf_min,
        ),
    )
    group.add_argument(
        "--reference-af-column",
        dest="reference_af_column",
        metavar="TAG",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Reference-population allele-frequency INFO tag",
            reference_af_column,
        ),
    )
    group.add_argument(
        "--maximum-af-difference",
        dest="maximum_af_difference",
        type=float,
        metavar="VALUE",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Maximum accepted absolute difference between study and reference "
            "allele frequencies",
            maximum_af_difference,
        ),
    )
    group.add_argument(
        "--minimum-info",
        dest="minimum_info",
        type=float,
        metavar="VALUE",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Minimum INFO or MaCH Rsq score accepted by the QC policy",
            defaults.info_min,
        ),
    )
    group.add_argument(
        "--maximum-info",
        dest="maximum_info",
        type=float,
        metavar="VALUE",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Maximum INFO or MaCH Rsq score accepted by the QC policy",
            defaults.info_max,
        ),
    )
    for option, destination, label, value in (
        (
            "--missing-pvalue-action",
            "missing_pvalue_action",
            "Action when the association P value is missing",
            defaults.missing_pvalue_action,
        ),
        (
            "--missing-af-action",
            "missing_af_action",
            "Action when a required allele frequency is missing",
            defaults.missing_af_action,
        ),
        (
            "--missing-info-action",
            "missing_info_action",
            "Action when imputation quality is missing",
            defaults.missing_info_action,
        ),
    ):
        group.add_argument(
            option,
            dest=destination,
            choices=("keep", "remove"),
            metavar="ACTION",
            default=argparse.SUPPRESS,
            help=help_with_default(label, value),
        )

    group.add_argument(
        "--include-indels",
        action="store_true",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Allow indels and other non-SNP variants to pass this policy",
            (
                "indels are included"
                if defaults.include_indels
                else "indels are not included"
            ),
        ),
    )
    group.add_argument(
        "--remove-palindromic",
        action="store_true",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Fail frequency-ambiguous A/T and C/G SNPs",
            (
                "palindromic variants are removed"
                if defaults.remove_palindromic
                else "palindromic variants are not removed"
            ),
        ),
    )
    group.add_argument(
        "--palindromic-af-lower",
        type=float,
        default=argparse.SUPPRESS,
        metavar="VALUE",
        help=help_with_default(
            "Lower allele-frequency ambiguity bound used with "
            "--remove-palindromic",
            defaults.palindromic_lower,
        ),
    )
    group.add_argument(
        "--palindromic-af-upper",
        type=float,
        default=argparse.SUPPRESS,
        metavar="VALUE",
        help=help_with_default(
            "Upper allele-frequency ambiguity bound used with "
            "--remove-palindromic",
            defaults.palindromic_upper,
        ),
    )
    group.add_argument(
        "--remove-mhc",
        action="store_true",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Fail variants in the configured build-specific MHC region",
            (
                "MHC variants are removed"
                if defaults.remove_mhc
                else "MHC variants are not removed"
            ),
        ),
    )
    return group


# ==================================================================
# ARGUMENT PARSER BUILDER (BASE + DIRECT ONLY)
# ==================================================================
def sumstat_summary_arg_parser() -> argparse.ArgumentParser:
    from postgwas.config import load_configuration
    from postgwas.core.ui import help_with_default

    defaults = load_configuration().modules.qc_summary
    parser = argparse.ArgumentParser(
        add_help=False,
        conflict_handler="resolve",
    )
    group = add_variant_qc_policy_arguments(
        parser,
        defaults=defaults.rules,
        reference_af_column=defaults.reference_af_column,
        maximum_af_difference=defaults.rules.maximum_af_difference,
        group_title="GWAS-VCF QC policy",
    )
    group.add_argument(
        "--sample-size-reference-quantile",
        type=float,
        metavar="VALUE",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Raw usable-Neff quantile used as the low-Neff reference",
            defaults.rules.sample_size_reference_quantile,
        ),
    )
    group.add_argument(
        "--sample-size-minimum-fraction",
        type=float,
        metavar="VALUE",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Fraction of the raw Neff reference used as the low-Neff threshold",
            defaults.rules.sample_size_minimum_fraction_of_reference,
        ),
    )
    return parser




# --- 1. BASE PARSER FOR STEP 2 ARGUMENTS ---

def get_annot_ldblock_parser(add_help=False):
    """
    Returns the core ArgumentParser for annot_ldblock arguments, structured using a group.
    """
    from postgwas.config import load_configuration
    from postgwas.core.ui import help_with_default, style_cli_requirement

    configuration = load_configuration()
    defaults = configuration.modules.ld_annotation
    populations = [population.value for population in defaults.populations]
    bed_filename_template = defaults.bed_filename_template
    info_field_template = defaults.info_field_template
    displayed_bed_pattern = bed_filename_template.replace(
        "{genome_build}", "<BUILD>"
    ).replace("{population}", "<POPULATION>")
    displayed_info_field = info_field_template.replace(
        "{population}", "<POP>"
    )
    example_bed_filename = bed_filename_template.format(
        genome_build=next(iter(configuration.resources.genomes)),
        population=populations[0],
    )
    example_bed_row = (
        "1<TAB>16103<TAB>2047216<TAB>"
        f"{populations[0]}-1_16103_2047216"
    )
    # Create bare parser for inheritance
    parser = argparse.ArgumentParser(add_help=add_help)

    # 💥 DEFINE ARGUMENT GROUP for Step 2 💥
    group = parser.add_argument_group("LD-block annotation")

    group.add_argument(
        "--ld-block-populations",
        dest="ld_block_populations",
        nargs="+",
        choices=list(configuration.resources.populations),
        metavar="POPULATION",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Populations to annotate (e.g., EUR AFR). "
            "One exactly named BED file is required for each selected population",
            " ".join(populations),
        ),
    )
    group.add_argument(
        "--ld-region-dir",
        default=argparse.SUPPRESS,
        metavar="PATH",
        help=style_cli_requirement(
            "Directory containing LD-block BED files.\n\n"
            "[bold cyan]File naming:[/bold cyan]\n"
            "PostGWAS reads BUILD from the input VCF and expects each selected "
            f"population file to be named exactly {displayed_bed_pattern}; for "
            f"example, {example_bed_filename}. The spelling and .bed.gz suffix "
            "must be exact.\n\n"
            "[bold cyan]Required BED columns (tab-separated):[/bold cyan]\n"
            "1. CHROM — exact input-VCF contig name, such as 1 or chr1.\n"
            "2. START — zero-based block start; included.\n"
            "3. END — zero-based block end; excluded.\n"
            f"4. BLOCK_LABEL — nonempty value written to {displayed_info_field}. For "
            "region clumping, it must end with _<START>_<END>, using the same "
            "values as columns 2 and 3.\n"
            "Additional trailing columns are allowed but ignored.\n"
            f"Example row: {example_bed_row}\n\n"
            "[bold bright_yellow]Warning:[/bold bright_yellow] "
            "[bright_yellow]Every BED file must have been generated for the same "
            "genome build declared inside the input GWAS-VCF; PostGWAS cannot "
            "infer a BED file's build from its coordinates.[/bright_yellow]",
            required=True,
        ),
    )


    return parser




# ------------------------------------------------------------
# Formatter-specific parser
# ------------------------------------------------------------
def get_formatter_parser(*, direct_controls=False):
    from typing import get_args

    from postgwas.config import load_configuration
    from postgwas.config.models.modules.formatting import FormattingDuplicatePolicy
    from postgwas.core.variant_identifiers import VariantIdentifierType
    from postgwas.core.ui import format_cli_choices, help_with_default
    from postgwas.modules.formatting.contracts import FORMAT_CONTRACTS

    configuration = load_configuration()
    defaults = configuration.modules.formatting
    available = list(defaults.format_order)
    parser = argparse.ArgumentParser(
        add_help=False,
        formatter_class=RawTextRichHelpFormatter,
    )
    group = parser.add_argument_group("Output formats")

    group.add_argument(
        "--format",
        choices=available,
        nargs="+",
        default=argparse.SUPPRESS,
        metavar="FORMAT",
        help=(
            "One or more downstream inputs to create.\n%s.\n"
            "Format details:\n%s\n"
            "If omitted, modules.formatting.formats is read from --run-config."
            % (
                format_cli_choices(available),
                "\n".join(
                    "  %-16s %s" % (target, FORMAT_CONTRACTS[target].name)
                    for target in available
                ),
            )
        ),
    )
    identifiers = parser.add_argument_group("Variant identifiers")
    identifiers.add_argument(
        "--variant-id-type",
        choices=get_args(VariantIdentifierType),
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Identifier convention written to every explicitly selected format. "
            "Use 'rsid' to extract an rsID from the VCF ID field or 'unique' to "
            "construct the configured chromosome-position-REF-ALT identifier. "
            "Per-target YAML settings can use different conventions",
            defaults.variant_identifiers.default_type,
        ),
    )
    identifiers.add_argument(
        "--duplicate-id-policy",
        choices=get_args(FormattingDuplicatePolicy),
        default=argparse.SUPPRESS,
        metavar="POLICY",
        help=help_with_default(
            "Action for conflicting duplicate IDs after collapsing exact repeated "
            "records and applying supported reference matching (currently LDSC "
            "--merge-alleles). exclude_all "
            "removes every ambiguous record; error stops; most_significant, "
            "highest_maf, and highest_info retain only a unique best-ranked row. "
            "A missing ranking value or tied best rank excludes the whole group. "
            "Per-target YAML settings can use different policies",
            defaults.variant_identifiers.default_duplicate_policy,
        ),
    )
    if direct_controls:
        group.add_argument(
            "--run-config",
            default=argparse.SUPPRESS,
            metavar="PATH",
            help=(
                "YAML file containing run and module settings. Command-line values "
                "override the matching YAML values."
            ),
        )

    return parser




def get_ld_clump_parser(add_help=False):
    """Return YAML-backed LD-clumping options without argparse defaults."""
    from typing import get_args

    from postgwas.config import load_configuration
    from postgwas.config.models.modules.ld_clumping import LDClumpingMethod
    from postgwas.core.ui import help_with_default
    from postgwas.modules.gcta_cojo.cli import get_gcta_cojo_selection_parser

    defaults = load_configuration().modules.ld_clumping
    parser = argparse.ArgumentParser(
        add_help=add_help,
        parents=[get_gcta_cojo_selection_parser(
            include_genome_build=False,
            include_reference_population=False,
            reference_requirement="when cojo-slct is selected",
        )],
    )
    general = parser.add_argument_group("LD-clumping analysis")
    general.add_argument(
        "--clumping-methods",
        nargs="+",
        choices=get_args(LDClumpingMethod),
        metavar="METHOD",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Analyses to run: region selects one association per annotated LD "
            "block; standard performs two-stage r² clumping against the tabix "
            "LD reference; cojo-slct runs GCTA stepwise selection and groups "
            "the selected signals into physical loci",
            " ".join(defaults.methods),
        ),
    )
    general.add_argument(
        "--lead-p", metavar="P", type=float, default=argparse.SUPPRESS,
        help=help_with_default(
            "P-value threshold for lead SNPs", defaults.lead_pvalue,
        ),
    )
    general.add_argument(
        "--candidate-p", metavar="P", type=float, default=argparse.SUPPRESS,
        help=help_with_default(
            "Maximum P-value for GWAS variants included as candidate SNPs",
            defaults.candidate_pvalue,
        ),
    )
    general.add_argument(
        "--remove-mhc",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Exclude the build-specific MHC interval before each selected analysis",
            defaults.remove_mhc,
        ),
    )

    standard = parser.add_argument_group("Standard r² clumping")
    standard.add_argument(
        "--ld-folder",
        metavar="PATH",
        default=argparse.SUPPRESS,
        # Deliberately unvalidated here. The directory is checked once, against
        # the resolved configuration, so the same message is reported whether
        # the value came from the command line or from YAML.
        help=(
            "Required when standard clumping is selected: directory containing "
            "the configured LD manifest, variant inventories, and orientation-"
            "specific tabix-indexed files."
        ),
    )
    standard.add_argument(
        "--r2-clump", metavar="R2", type=float, default=argparse.SUPPRESS,
        help=help_with_default(
            "Minimum r² for defining independent significant SNPs",
            defaults.clump_r2,
        ),
    )
    standard.add_argument(
        "--r2-lead", metavar="R2", type=float, default=argparse.SUPPRESS,
        help=help_with_default(
            "Minimum r² for the second clumping step that defines lead SNPs",
            defaults.lead_r2,
        ),
    )
    standard.add_argument(
        "--ld-window-kb", dest="ld_window_kb", metavar="KB", type=int,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Maximum distance from an index SNP for assigning LD partners",
            defaults.window_kb,
        ),
    )
    standard.add_argument(
        "--merge-dist", metavar="BP", type=int, default=argparse.SUPPRESS,
        help=help_with_default(
            "Maximum distance between independent-significant-SNP LD blocks "
            "that are merged into one genomic locus",
            defaults.merge_distance_bp,
        ),
    )
    standard.add_argument(
        "--missing-index-action",
        choices=("error", "warning_skip"),
        metavar="ACTION",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Action when a significant index SNP is absent from, allele-"
            "incompatible with, or below the configured MAF threshold of the "
            "LD-reference inventory. warning_skip excludes and audits the index; "
            "error stops the analysis. A verified index without LD pairs is "
            "always retained as self-only",
            defaults.missing_index_action,
        ),
    )
    standard.add_argument(
        "--missing-chromosome-action",
        choices=("error", "warning"),
        metavar="ACTION",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Action when a chromosome containing significant variants has no "
            "required LD, reverse-index, inventory, or tabix resource. warning "
            "skips that chromosome and records the excluded variant count; "
            "error stops the analysis",
            defaults.missing_chromosome_action,
        ),
    )
    cojo = parser.add_argument_group("GCTA-COJO physical locus definition")
    cojo.add_argument(
        "--cojo-merge-dist", type=int, metavar="BP",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Maximum consecutive distance between GCTA-selected signals "
            "grouped into one physical locus after model fitting",
            defaults.cojo.merge_distance_bp,
        ),
    )
    return parser





def get_inputvcf_parser():
    """
    Common parser for shared input arguments across PostGWAS modules:
    - VCF input
    - Genome version
    - Sample/phenotype ID
    - Output folder
    Safe for inheritance across multiple pipeline steps.
    """
    parser = argparse.ArgumentParser(
        add_help=False,
        formatter_class=RichHelpFormatter
    )

    input_group = parser.add_argument_group("Input file")

    # -------- VCF INPUT --------
    input_group.add_argument(
        "--vcf",
        metavar="PATH",
        type=validate_path(
            must_exist=True,
            must_be_file=True,
            must_not_be_empty=True,
            allowed_suffixes=[".vcf", ".vcf.gz"]
        ),
        help="REQUIRED. Harmonised GWAS-VCF file to analyse."
    )

    return parser




def get_genome_build_parser(
    *,
    available_builds=None,
    default_build=None,
    suppress_default=False,
    display_default=True,
):
    """
    Common parser for shared input arguments across PostGWAS modules:
    - VCF input
    - Genome version
    - Sample/phenotype ID
    - Output folder
    Safe for inheritance across multiple pipeline steps.
    """
    if available_builds is None:
        from postgwas.config import load_configuration

        available_builds = list(load_configuration().resources.genomes)
    from postgwas.core.ui import help_with_default

    available_builds = list(available_builds)
    if not available_builds:
        raise ValueError("At least one configured genome build is required")
    resolved_default = default_build
    if resolved_default is None and (not suppress_default or display_default):
        resolved_default = available_builds[0]
    if resolved_default is not None and resolved_default not in available_builds:
        raise ValueError("The default genome build must be one of available_builds")
    parser = argparse.ArgumentParser(
        add_help=False,
        formatter_class=RichHelpFormatter
    )

    input_group = parser.add_argument_group("Genome coordinates")

    # -------- GENOME VERSION --------
    build_help = (
        "Genome build used by the input VCF and all reference files. "
        "Do not mix resources from different builds"
    )
    input_group.add_argument(
        "--genome-build",
        default=argparse.SUPPRESS if suppress_default else resolved_default,
        choices=available_builds,
        metavar="BUILD",
        help=(
            help_with_default(build_help, resolved_default)
            if display_default
            else build_help + ". If omitted, the selected modules' canonical "
            "run configuration supplies the value."
        ),
    )
    return parser


def get_pipeline_genome_build_parser(add_help=False):
    """Return one default-free genome-build override shared by pipeline modules."""
    from postgwas.config import load_configuration

    configuration = load_configuration()
    builds = list(configuration.resources.genomes)
    return get_genome_build_parser(
        available_builds=builds,
        suppress_default=True,
        display_default=False,
    )


def _get_population_parser(
    module_name: str,
    add_help=False,
    *,
    suppress_default: bool = False,
):
    """Return a module-specific, YAML-backed reference-population option."""
    from postgwas.config import load_configuration
    from postgwas.core.ui import help_with_default

    configuration = load_configuration()
    configured = getattr(configuration.modules, module_name).population.value
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("Reference population")
    group.add_argument(
        "--population",
        choices=list(configuration.resources.populations),
        default=argparse.SUPPRESS if suppress_default else configured,
        metavar="CODE",
        help=help_with_default(
            "Population represented by the reference data",
            configured,
        ),
    )
    return parser


def get_imputation_population_parser(add_help=False):
    return _get_population_parser(
        "imputation", add_help=add_help, suppress_default=True,
    )


def get_ld_clumping_population_parser(add_help=False):
    return _get_population_parser(
        "ld_clumping", add_help=add_help, suppress_default=True,
    )


def get_ld_clumping_genome_build_parser(add_help=False):
    """Return the LD-clumping build override backed by canonical YAML."""
    from postgwas.config import load_configuration

    configuration = load_configuration()
    return get_genome_build_parser(
        available_builds=list(configuration.resources.genomes),
        default_build=configuration.modules.ld_clumping.genome_build.value,
        suppress_default=True,
    )

def get_common_out_parser():
    """
    Common parser for shared OUTPUT arguments across PostGWAS modules.

    NOTE:
    -----
    - No required=True
    - No filesystem calls (Path.cwd) during parser creation
    - Defaults resolved AFTER parsing
    """
    parser = argparse.ArgumentParser(
        add_help=False,
        formatter_class=RichHelpFormatter,
    )

    group = parser.add_argument_group("Output")

    # -------- SAMPLE ID --------
    group.add_argument(
        "--dataset-id",
        dest="dataset_id",
        metavar="NAME",
        type=validate_alphanumeric,
        help=(
            "REQUIRED. Short, unique name for the dataset, such as PGC3_SCZ_EUR. "
            "It is used in output file names."
        ),
    )

    # -------- OUTPUT DIRECTORY --------
    group.add_argument(
        "--output-directory",
        dest="output_directory",
        metavar="PATH",
        default=None,  # 🔑 SAFE
        type=validate_path(
            must_exist=False,
            must_be_dir=True,
            create_if_missing=True,
        ),
        help="REQUIRED. Folder where PostGWAS will save results and logs.",
    )

    return parser


def get_common_magma_covar_parser(add_help=False, *, compact_help=False):
    """Return MAGMAcovar options with optional concise pipeline guidance."""
    from typing import get_args

    from rich.markup import escape

    from postgwas.config import load_configuration
    from postgwas.config.models.modules.magmacovar import (
        MagmaCovarDirection,
        MagmaCovarMissingGenes,
        MagmaCovarMissingValues,
    )
    from postgwas.core.ui import (
        format_cli_default,
        format_cli_help_block,
        help_with_default,
    )

    defaults = load_configuration().modules.magmacovar
    accepted_modifiers = ", ".join(defaults.model_modifiers)
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("MAGMA gene-property settings")
    group.add_argument(
        "--covariates",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "Gene-level covariate table with a header, gene IDs in column one, "
            "numeric property columns, and missing values encoded as NA."
        ),
    )

    group.add_argument(
        "--covariate-model",
        metavar="MODEL",
        nargs="+",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Choose how MAGMA tests the properties in the covariate table. "
            "Accepted modifiers: %s. Omit this option to test every property "
            "separately" % accepted_modifiers,
            " ".join(defaults.model) if defaults.model else "none",
        ),
    )

    group.add_argument(
        "--covariate-direction",
        choices=get_args(MagmaCovarDirection),
        metavar="DIRECTION",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Testing direction passed through MAGMA's direction-covar model "
            "modifier. Greater and smaller are one-sided tests. Two-sided is "
            "MAGMA's default for any gene property. The original FLAMES README "
            "relies on that default, while tissue-specific analysis uses greater",
            defaults.direction,
        ),
    )
    group.add_argument(
        "--minimum-genes",
        type=int,
        metavar="N",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Minimum overlapping genes required before execution and in every "
            "published COVAR result",
            defaults.minimum_genes,
        ),
    )
    group.add_argument(
        "--covariate-missing-values",
        choices=get_args(MagmaCovarMissingValues),
        metavar="ACTION",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "How MAGMA handles NA values after a property passes max-miss. "
            "Drop excludes genes having any NA in the loaded table; median "
            "and mean replace each NA using that property's observed values",
            defaults.input.missing_values,
        ),
    )
    group.add_argument(
        "--covariate-max-miss",
        type=float,
        metavar="FRACTION",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Largest permitted missing fraction for every property. Values "
            "must be between 0 and MAGMA's maximum 0.2; equality is accepted",
            defaults.input.maximum_missing_fraction,
        ),
    )
    group.add_argument(
        "--covariate-missing-genes",
        choices=get_args(MagmaCovarMissingGenes),
        metavar="ACTION",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "How MAGMA handles eligible .genes.raw genes absent from the "
            "covariate table. Drop excludes them; fill counts them as missing "
            "for every property before applying max-miss",
            defaults.input.missing_genes,
        ),
    )

    if compact_help:
        return parser

    accepted_directions = ", ".join(get_args(MagmaCovarDirection))
    direction_lines = [
        "[bold bright_yellow]Exact MAGMA direction modifier[/bold bright_yellow]",
        "  PostGWAS emits [cyan]direction-covar=<DIRECTION>[/cyan] for gene "
        "properties, not the blanket direction= or gene-set-only direction-sets=.",
    ]
    direction_lines.append(
        "  Accepted values: [cyan]%s[/cyan]. MAGMA uses [cyan]smaller[/cyan], "
        "not less." % escape(accepted_directions)
    )
    safety_lines = [
        "[bold cyan]Property missingness guard[/bold cyan]",
        "MAGMA removes a gene property when its missing fraction is greater "
        "than max-miss. PostGWAS checks every property first and stops with "
        "the property name, so a removed property cannot look like a test that "
        "was never requested.",
        "",
        "[bold bright_yellow]Failure rule[/bold bright_yellow]",
        "  A property fails only when missing fraction > max-miss; equality is "
        "accepted.",
        "",
        "[bold bright_yellow]Missingness denominator[/bold bright_yellow]",
        "[cyan]missing-genes=fill[/cyan]",
        "  All eligible genes from the .genes.raw file. Genes absent from the "
        "covariate table are filled and count as missing.",
        "[cyan]missing-genes=drop[/cyan]",
        "  Only gene IDs overlapping the .genes.raw file and covariate table. "
        "Absent genes are excluded.",
        "",
        *direction_lines,
        "",
        format_cli_default(
            "max-miss=%.12g, missing-genes=%s, missing-values=%s"
            % (
                defaults.input.maximum_missing_fraction,
                defaults.input.missing_genes,
                defaults.input.missing_values,
            )
        ),
    ]
    parser.add_argument_group(
        "Pre-MAGMA safety checks",
        format_cli_help_block("\n".join(safety_lines)),
    )

    choice_lines = ["Choose the model that matches your research question:"]
    for use_case in defaults.model_use_cases.values():
        selected_model = (
            "--covariate-model " + " ".join(use_case.model)
            if use_case.model else "omit --covariate-model"
        )
        choice_lines.extend((
            "",
            "[bold cyan]%s[/bold cyan]" % escape(use_case.label),
            "  [bold]Question:[/bold] %s" % escape(use_case.question),
            "  [bold bright_yellow]Model:[/bold bright_yellow] "
            "[cyan]%s[/cyan]" % escape(selected_model),
            "  [bold bright_yellow]Direction:[/bold bright_yellow] "
            "[cyan]--covariate-direction %s[/cyan]"
            % escape(use_case.direction),
        ))
    parser.add_argument_group(
        "How to choose a MAGMA model",
        format_cli_help_block("\n".join(choice_lines)),
    )

    modifier_lines = [
        "Use these only for a planned conditional or multivariable analysis:",
        "",
        "[bold bright_yellow]Values used below[/bold bright_yellow]",
    ]
    for placeholder, description in defaults.model_placeholders.items():
        modifier_lines.extend((
            "[cyan]<%s>[/cyan]" % escape(placeholder),
            "  %s" % escape(description),
        ))
    modifier_lines.extend((
        "",
        "Example property names are illustrative; use exact headers from your table.",
        "",
        "Interaction and gene-set-only modifiers are unavailable.",
        "This interface does not provide --set-annot.",
    ))
    parser.add_argument_group(
        "Advanced MAGMA model modifiers",
        format_cli_help_block("\n".join(modifier_lines)),
    )

    for modifier, specification in defaults.model_modifiers.items():
        detail_lines = [
            "[bold]Purpose:[/bold] %s" % escape(specification.description),
            "",
            "[bold bright_yellow]Accepted forms[/bold bright_yellow]",
        ]
        detail_lines.extend(
            "  [cyan]%s[/cyan]" % escape(syntax)
            for syntax in specification.syntax
        )
        detail_lines.extend((
            "",
            "[bold magenta]Example command[/bold magenta]",
            "  [cyan]%s[/cyan]" % escape(specification.example),
            "",
            "[bold]Published GWAS example[/bold]",
            "  %s" % escape(specification.published_example),
        ))
        parser.add_argument_group(
            "MAGMA model modifier: %s" % modifier,
            format_cli_help_block("\n".join(detail_lines)),
        )
    return parser




def get_common_magma_assoc_parser(add_help=False):
    """
    Parser containing only MAGMA association-related arguments.
    This is parent-safe and can be inserted into direct/pipeline mode.
    """
    from postgwas.config import load_configuration
    from postgwas.core.ui import help_with_default

    defaults = load_configuration().modules.magma
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("MAGMA settings")

    # Required inputs (but NOT marked required here → handled by subparser)
    group.add_argument(
        "--magma-ld-reference",
        dest="magma_ld_reference",
        type=lambda p: validate_prefix_files(
            p, defaults.input.required_reference_extensions,
        ),
        metavar="PREFIX",
        default=argparse.SUPPRESS,
        help=(
        "Prefix of the PLINK reference panel for MAGMA.\n"
        "Configured companion files: "
        f"{', '.join(defaults.input.required_reference_extensions)}.\n"
        "Make sure the summary statistics and the LD reference "
        "use the same genome build."
        ) ,
    )

    group.add_argument(
        "--gene-location-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        type=validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True),
        help=help_with_default(
            "[bold bright_red]Required[/bold bright_red]: MAGMA gene location file (.loc) matching the genome build. "
            "Must contain columns in the following order: Gene ID, Chromosome, Gene Start, Gene End. "
            "Gene IDs must match those used in pathway files and gene-covariate files",
            defaults.input.gene_location_file,
        ),
    )

    group.add_argument(
        "--gene-set-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        type=validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True),
        help=help_with_default(
            "Pathway file (.gmt). MSigDB or custom GMT supported",
            defaults.input.gene_set_file,
        ),
    )

    group.add_argument(
        "--window-upstream",
        type=int,
        default=argparse.SUPPRESS,
        metavar="KB",
        help=help_with_default(
            "Upstream gene window in kb",
            defaults.gene_window_upstream_kb,
        ),
    )

    group.add_argument(
        "--window-downstream",
        type=int,
        default=argparse.SUPPRESS,
        metavar="KB",
        help=help_with_default(
            "Downstream gene window in kb",
            defaults.gene_window_downstream_kb,
        ),
    )

    group.add_argument(
        "--gene-model",
        default=argparse.SUPPRESS,
        metavar="MODEL",
        help=help_with_default(
            "MAGMA model supported with SNP p-value input",
            defaults.gene_model,
        ),
    )

    group.add_argument(
        "--sample-size-column",
        dest="sample_size_column",
        default=argparse.SUPPRESS,
        metavar="COLUMN",
        help=help_with_default(
            "Column name for sample sizes in *.pval input",
            defaults.input.sample_size_column,
        ),
    )

    return parser


def get_plink_binary_parser(add_help=False):
    """
    Parser containing only MAGMA association-related arguments.
    This is parent-safe and can be inserted into direct/pipeline mode.
    """
    from postgwas.config import load_configuration
    from postgwas.core.ui import help_with_default

    configured = load_configuration().resources.executables.plink
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("PLINK binary")
    # Required inputs (but NOT marked required here → handled by subparser)
    group.add_argument(
        "--plink",
        default=argparse.SUPPRESS,
        metavar="PATH",
        help=help_with_default(
            "Full path or command name for the PLINK executable used for LD",
            configured,
        ),
    )
    return parser


def get_tabix_binary_parser(add_help=False):
    """Return the shared tabix executable option."""
    from postgwas.config import load_configuration
    from postgwas.core.ui import help_with_default

    configured = load_configuration().resources.executables.tabix
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("External software")
    group.add_argument(
        "--tabix",
        default=argparse.SUPPRESS,
        metavar="PATH",
        help=help_with_default(
            "Full path or command name for tabix",
            configured,
        ),
    )
    return parser


def get_common_pops_parser(add_help=False):
    """Return reusable PoPS options whose defaults come only from YAML."""
    from typing import get_args

    from postgwas.config import load_configuration
    from postgwas.config.models.modules.pops import (
        PopsGeneUniversePolicy,
        PopsMethod,
    )
    from postgwas.core.ui import help_with_default

    module = load_configuration().modules.pops
    parser = argparse.ArgumentParser(add_help=add_help)

    grp_feats = parser.add_argument_group("PoPS input")
    grp_feats.add_argument(
        "--feature-matrix-prefix",
        metavar="PREFIX",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Required standard-PoPS prefix shared by its matrix, column, and row "
            "files",
            module.feature_matrix_prefix,
        ),
    )
    grp_feats.add_argument(
        "--feature-matrix-chunks",
        type=int,
        metavar="N",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Number of consecutively numbered feature-matrix chunks; this must "
            "match the files below --feature-matrix-prefix",
            module.feature_matrix_chunks,
        ),
    )
    grp_feats.add_argument(
        "--pops-gene-location-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Required standard-PoPS gene annotation containing the configured "
            "Ensembl gene, chromosome, and TSS columns",
            module.gene_location_file,
        ),
    )
    grp_feats.add_argument(
        "--gene-universe-policy",
        choices=get_args(PopsGeneUniversePolicy),
        metavar="POLICY",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "How MAGMA genes with Z-scores absent from PoPS resources are handled: "
            "strict fails before fitting; intersect explicitly derives aligned "
            "compatible and excluded MAGMA files with a full audit report",
            module.gene_universe_policy,
        ),
    )
    grp_cov = parser.add_argument_group("PoPS covariates")
    grp_cov.add_argument(
        "--use-magma-covariates",
        action="store_true",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Project out technical covariates read from MAGMA .genes.raw, such as "
            "gene size and density, before fitting. This does not consume a "
            "MAGMACOVAR .gsa.out file or use --covariate-model",
            module.use_magma_covariates,
        ),
    )
    grp_cov.add_argument(
        "--ignore-magma-covariates",
        action="store_false",
        dest="use_magma_covariates",
        default=argparse.SUPPRESS,
        help="Do not project out MAGMA covariates.",
    )
    grp_cov.add_argument(
        "--use-magma-error-covariance",
        action="store_true",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Use the MAGMA gene-error covariance during fitting",
            module.use_magma_error_covariance,
        ),
    )
    grp_cov.add_argument(
        "--ignore-magma-error-covariance",
        action="store_false",
        dest="use_magma_error_covariance",
        default=argparse.SUPPRESS,
        help="Do not use the MAGMA gene-error covariance.",
    )
    grp_cov.add_argument(
        "--target-score-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "Custom gene-score table containing the configured gene-ID and "
            "score columns; mutually exclusive with --magma-association-prefix. "
            "The option name is retained for compatibility."
        ),
    )
    grp_cov.add_argument(
        "--target-covariates-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help="Optional covariate table aligned to the custom gene scores.",
    )
    grp_cov.add_argument(
        "--target-error-covariance-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "Optional SciPy NPZ or NumPy NPY covariance for custom gene scores."
        ),
    )
    grp_cov.add_argument(
        "--covariate-projection-chromosomes",
        nargs="+",
        metavar="CHROM",
        default=argparse.SUPPRESS,
        help="Chromosomes used to estimate covariate effects; omit to use all annotated chromosomes.",
    )
    grp_cov.add_argument(
        "--remove-hla-during-covariate-projection",
        action="store_true",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Apply upstream PoPS HLA exclusion during covariate projection",
            module.remove_hla_during_covariate_projection,
        ),
    )
    grp_cov.add_argument(
        "--keep-hla-during-covariate-projection",
        action="store_false",
        dest="remove_hla_during_covariate_projection",
        default=argparse.SUPPRESS,
        help="Keep HLA-region genes during covariate projection.",
    )

    grp_sel = parser.add_argument_group("PoPS feature selection")
    grp_sel.add_argument(
        "--feature-subset-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help="Optional file containing the only feature names eligible for selection.",
    )
    grp_sel.add_argument(
        "--control-features-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "Optional file listing features retained in the model; the official "
            "PoPS example and FUMA FLAMES workflow supply this file."
        ),
    )
    grp_sel.add_argument(
        "--feature-selection-chromosomes",
        nargs="+",
        metavar="CHROM",
        default=argparse.SUPPRESS,
        help="Chromosomes used for marginal feature selection; omit to use all annotated chromosomes.",
    )
    grp_sel.add_argument(
        "--feature-selection-p-cutoff",
        type=float,
        metavar="P",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Strict marginal-association p-value cutoff for feature selection",
            module.feature_selection_p_cutoff,
        ),
    )
    grp_sel.add_argument(
        "--maximum-selected-features",
        type=int,
        metavar="N",
        default=argparse.SUPPRESS,
        help="Optional positive maximum excluding configured control features.",
    )
    grp_sel.add_argument(
        "--forward-selected-features",
        type=int,
        metavar="N",
        default=argparse.SUPPRESS,
        help="Optional positive number selected by upstream forward stepwise selection.",
    )
    grp_sel.add_argument(
        "--remove-hla-during-feature-selection",
        action="store_true",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Apply upstream PoPS HLA exclusion during feature selection",
            module.remove_hla_during_feature_selection,
        ),
    )
    grp_sel.add_argument(
        "--keep-hla-during-feature-selection",
        action="store_false",
        dest="remove_hla_during_feature_selection",
        default=argparse.SUPPRESS,
        help="Keep HLA-region genes during feature selection.",
    )

    grp_train = parser.add_argument_group("PoPS training")
    grp_train.add_argument(
        "--training-chromosomes",
        nargs="+",
        metavar="CHROM",
        default=argparse.SUPPRESS,
        help="Chromosomes used to fit coefficients; omit to use all annotated chromosomes.",
    )
    grp_train.add_argument(
        "--remove-hla-during-training",
        action="store_true",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Apply upstream PoPS HLA exclusion during coefficient fitting",
            module.remove_hla_during_training,
        ),
    )
    grp_train.add_argument(
        "--keep-hla-during-training",
        action="store_false",
        dest="remove_hla_during_training",
        default=argparse.SUPPRESS,
        help="Keep HLA-region genes during coefficient fitting.",
    )

    grp_misc = parser.add_argument_group("PoPS model settings")
    grp_misc.add_argument(
        "--method",
        choices=get_args(PopsMethod),
        metavar="METHOD",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Upstream PoPS coefficient model: ridge, lasso or linreg",
            module.method,
        ),
    )
    grp_misc.add_argument(
        "--save-matrix-files",
        action="store_true",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Write upstream PoPS training and prediction matrices",
            module.save_matrix_files,
        ),
    )
    grp_misc.add_argument(
        "--no-save-matrix-files",
        action="store_false",
        dest="save_matrix_files",
        default=argparse.SUPPRESS,
        help="Do not write upstream PoPS matrix files.",
    )
    return parser


from rich_argparse import RawTextRichHelpFormatter


def get_common_sumstat_filter_parser(add_help=False):
    from postgwas.config import load_configuration
    from postgwas.core.ui import help_with_default

    module = load_configuration().modules.filtering
    mhc_defaults = {
        field: ", ".join(
            "%s=%s" % (build.value, getattr(region, field))
            for build, region in module.mhc_regions.items()
        )
        for field in ("chromosome", "start", "end")
    }
    parser = argparse.ArgumentParser(
        formatter_class=RawTextRichHelpFormatter,
        add_help=add_help,
        conflict_handler="resolve",
    )

    group = add_variant_qc_policy_arguments(
        parser,
        defaults=module,
        reference_af_column=module.reference_population_tag,
        maximum_af_difference=module.frequency_difference_max,
        group_title="Quality-control filters",
    )

    group.add_argument(
        "--mhc-chrom",
        default=argparse.SUPPRESS,
        metavar="CHROM",
        help=help_with_default(
            "Override the MHC chromosome for the genome build inferred from the "
            "input VCF; otherwise use mhc_regions",
            mhc_defaults["chromosome"],
        ),
    )

    group.add_argument(
        "--mhc-start",
        type=int,
        default=argparse.SUPPRESS,
        metavar="POSITION",
        help=help_with_default(
            "Override the one-based inclusive MHC start for the genome build "
            "inferred from the input VCF; otherwise use mhc_regions",
            mhc_defaults["start"],
        ),
    )

    group.add_argument(
        "--mhc-end",
        type=int,
        default=argparse.SUPPRESS,
        metavar="POSITION",
        help=help_with_default(
            "Override the one-based inclusive MHC end for the genome build "
            "inferred from the input VCF; otherwise use mhc_regions",
            mhc_defaults["end"],
        ),
    )

    group.add_argument(
        "--write-soft-filter-vcf",
        action="store_true",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Also retain every input variant in an indexed audit VCF and add "
            "configured FILTER IDs for each active removal rule it fails",
            module.write_soft_filter_vcf,
        ),
    )

    return parser




# ============================================================
# COMMON FLAMES PARSER (shared by DIRECT + PIPELINE)
# ============================================================
def get_flames_common_parser(add_help=False):
    """Return YAML-backed FLAMES options shared by direct and pipeline modes."""
    from postgwas.config import load_configuration
    from postgwas.core.ui import help_with_conditional_requirement, help_with_default

    defaults = load_configuration()
    module = defaults.modules.flames
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("FLAMES analysis settings")
    group.add_argument(
        "--flames-annotation-directory",
        dest="annotation_resource_directory",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "Directory containing the pinned upstream FLAMES annotation bundle. "
            "PostGWAS validates its configured structure before execution."
        ),
    )
    group.add_argument(
        "--flames-genome-build",
        choices=list(defaults.resources.genomes),
        metavar="BUILD",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Genome build shared by credible-set coordinates and upstream gene evidence; "
            "PostGWAS does not infer it",
            module.genome_build.value,
        ),
    )
    group.add_argument(
        "--flames-model-directory",
        dest="flames_model_directory",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "Directory containing the configured upstream model and feature manifest. "
            "If omitted, the packaged model resource declared in YAML is used."
        ),
    )
    group.add_argument(
        "--flames-vep-mode",
        choices=("api", "local"),
        metavar="MODE",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Use the upstream Ensembl VEP API or a configured local VEP installation",
            module.vep_mode,
        ),
    )
    group.add_argument(
        "--vep-command",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=help_with_conditional_requirement(
            "VEP executable",
            "when --flames-vep-mode local is selected",
        ),
    )
    group.add_argument(
        "--vep-cache",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=help_with_conditional_requirement(
            "Build-compatible VEP cache directory",
            "when --flames-vep-mode local is selected",
        ),
    )
    group.add_argument(
        "--vep-cache-genome-build",
        metavar="BUILD",
        choices=list(defaults.resources.genomes),
        default=argparse.SUPPRESS,
        help=help_with_conditional_requirement(
            "Genome build explicitly declared for the local VEP cache",
            "when --flames-vep-mode local is selected",
        ),
    )
    group.add_argument(
        "--flames-cadd-mode",
        choices=("api", "local"),
        metavar="MODE",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Use the upstream CADD API or a configured local tabix-indexed CADD file",
            module.cadd_mode,
        ),
    )
    group.add_argument(
        "--cadd-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=help_with_conditional_requirement(
            "Genome-build-matched tabix-indexed CADD file",
            "when --flames-cadd-mode local is selected",
        ),
    )
    group.add_argument(
        "--cadd-genome-build",
        metavar="BUILD",
        choices=list(defaults.resources.genomes),
        default=argparse.SUPPRESS,
        help=help_with_conditional_requirement(
            "Genome build explicitly declared for the local CADD file",
            "when --flames-cadd-mode local is selected",
        ),
    )
    return parser



def get_common_imputation_parser(add_help: bool = False) -> argparse.ArgumentParser:
    from postgwas.config import load_configuration
    from postgwas.core.ui import help_with_default

    module = load_configuration().modules.imputation
    pred_ld = module.engines.pred_ld
    parser = argparse.ArgumentParser(add_help=add_help)
    grp = parser.add_argument_group("Imputation settings")

    grp.add_argument(
        "--imputation-engine",
        dest="imputation_engine",
        metavar="ENGINE",
        type=str,
        choices=[module.engine],
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Choose the configured imputation tool",
            module.engine,
        ),
    )

    grp.add_argument(
        "--imputation-ld-reference",
        dest="imputation_ld_reference",
        metavar="PATH",
        type=validate_path(must_exist=True, must_be_dir=True),
        default=argparse.SUPPRESS,
        help=(
            "[bold bright_red]Required[/bold bright_red]: Directory containing the LD "
            "reference panel for PRED-LD. "
            "Ensure that the [cyan]genome build[/cyan] of the LD reference "
            "matches the genome version of your summary statistics "
            "(e.g., both GRCh38)."
        ),
    )

    grp.add_argument(
        "--imputation-r2-threshold",
        dest="imputation_r2_threshold",
        metavar="R2",
        type=float,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "LD r² tagging threshold", pred_ld.minimum_r2,
        ),
    )

    grp.add_argument(
        "--imputation-minimum-maf",
        dest="imputation_minimum_maf",
        metavar="MAF",
        type=float,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Minor allele frequency threshold", pred_ld.minimum_maf,
        ),
    )

    grp.add_argument(
        "--ref",
        metavar="NAME",
        type=str,
        default=argparse.SUPPRESS,
        choices=[pred_ld.mode],
        help=help_with_default(
            "Reference panel name label", pred_ld.mode,
        ),
    )

    grp.add_argument(
        "--correlation-method",
        dest="corr_method",
        metavar="METHOD",
        type=str,
        default=argparse.SUPPRESS,
        choices=["pearson", "spearman"],
        help=help_with_default(
            "Correlation method used in QC/post-processing",
            pred_ld.correlation_method,
        ),

    )

    grp.add_argument(
        "--resource-directory",
        metavar="PATH",
        default=argparse.SUPPRESS,
        type=validate_path(must_exist=True, must_be_dir=True),
        help="[bold bright_red]Required[/bold bright_red]: PostGWAS resource directory.",
    )

    return parser



def get_ldsc_merge_alleles_parser(
    add_help: bool = False,
) -> argparse.ArgumentParser:
    """Expose one shared HapMap3 option to formatter and LDSC commands."""
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("LDSC allele matching")
    _add_ldsc_merge_alleles_argument(
        group,
        help_text=(
            "HapMap3 SNP-and-allele table (for example w_hm3.snplist). "
            "The direct formatter uses it optionally to select one rsID-and-allele "
            "compatible LDSC record. Without it, the selected duplicate-ID "
            "policy handles duplicated rsIDs. Heritability workflows require it "
            "and reuse the same file during munge_sumstats."
        ),
    )
    return parser


def _add_ldsc_merge_alleles_argument(
    group,
    *,
    help_text: str,
) -> None:
    """Add the shared LDSC allele-reference option to a context-specific group."""
    group.add_argument(
        "--merge-alleles",
        metavar="PATH",
        default=argparse.SUPPRESS,
        type=validate_path(
            must_exist=True,
            must_be_file=True,
            must_not_be_empty=True,
        ),
        help=help_text,
    )


def get_ldsc_common_parser(
    add_help: bool = False,
    *,
    pipeline_mode: bool = False,
) -> argparse.ArgumentParser:
    """
    Arguments shared by LDSC DIRECT and PIPELINE modes.

    DIRECT mode:
        --ldsc-input    → formatter-created LDSC input table
        --merge-alleles → HapMap3 SNP+allele list (w_hm3.snplist)
        --ref-ld-chr    → eur_w_ld_chr/ (directory)
        --w-ld-chr      → eur_w_ld_chr/ (directory)
        --samp-prev / --pop-prev

    PIPELINE mode:
        formatter output is supplied to the service by the pipeline context
        --merge-alleles → HapMap3 SNP+allele list (w_hm3.snplist)
        --info-min / --maf-min used during munge_sumstats
        (ref/w ld-chr + prevs also used)
    """
    from postgwas.config import load_configuration
    from postgwas.config.models.common import GenomeBuild, Population
    from postgwas.core.ui import help_with_default

    defaults = load_configuration()
    module = defaults.modules.ldsc
    if pipeline_mode:
        sample_prevalence_help = (
            "Optional sample-prevalence override for liability-scale h². When "
            "population prevalence is set and no CLI or YAML sample prevalence "
            "is supplied, PostGWAS uses the formatter-returned value; it does "
            "not recalculate it"
        )
        population_prevalence_help = (
            "Population prevalence that activates liability-scale h². Sample "
            "prevalence comes from CLI or YAML, or from the formatter result "
            "when neither supplies it"
        )
    else:
        sample_prevalence_help = (
            "Sample prevalence for liability-scale h². Direct mode requires "
            "it from --samp-prev or YAML together with population prevalence; "
            "PostGWAS does not calculate or infer it"
        )
        population_prevalence_help = (
            "Population prevalence for liability-scale h². Direct mode also "
            "requires sample prevalence from --samp-prev or YAML; omit population "
            "prevalence to run observed-scale h² only"
        )
    parser = argparse.ArgumentParser(add_help=add_help)

    reference_files = parser.add_argument_group("LDSC reference files")
    _add_ldsc_merge_alleles_argument(
        reference_files,
        help_text=(
            "HapMap3 SNP-and-allele table (for example w_hm3.snplist) used to "
            "retain allele-compatible variants and passed to munge_sumstats."
        ),
    )
    reference_files.add_argument(
        "--ref-ld-chr",
        metavar="PATH",
        default=argparse.SUPPRESS,
        type=validate_path(must_exist=True, must_be_dir=True),
        help="Directory containing chromosome-split reference LD scores.",
    )
    reference_files.add_argument(
        "--w-ld-chr",
        metavar="PATH",
        default=argparse.SUPPRESS,
        type=validate_path(must_exist=True, must_be_dir=True),
        help=(
            "Directory containing chromosome-split LDSC regression weights; "
            "this is often the same directory as --ref-ld-chr."
        ),
    )

    metadata = parser.add_argument_group("LDSC reference metadata")

    metadata.add_argument(
        "--ldsc-population",
        choices=[population.value for population in Population],
        default=argparse.SUPPRESS,
        metavar="CODE",
        help=help_with_default(
            "Population represented by the supplied LDSC reference files",
            module.population,
        ),
    )
    metadata.add_argument(
        "--ldsc-genome-build",
        choices=[build.value for build in GenomeBuild],
        default=argparse.SUPPRESS,
        metavar="BUILD",
        help=help_with_default(
            "Declared genome build of the supplied LDSC reference release",
            module.genome_build,
        ),
    )

    software = parser.add_argument_group("LDSC executables")
    software.add_argument(
        "--ldsc-executable",
        metavar="PATH_OR_COMMAND",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "CBIIT ldsc.py executable path or command name",
            defaults.resources.executables.ldsc,
        ),
    )
    software.add_argument(
        "--munge-sumstats-executable",
        metavar="PATH_OR_COMMAND",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "CBIIT munge_sumstats.py executable path or command name",
            defaults.resources.executables.munge_sumstats,
        ),
    )

    prev = parser.add_argument_group("Liability-scale settings")

    prev.add_argument(
        "--samp-prev",
        metavar="VALUE",
        type=float,
        default=argparse.SUPPRESS,
        help=help_with_default(
            sample_prevalence_help,
            module.sample_prevalence,
        ),
    )

    prev.add_argument(
        "--pop-prev",
        metavar="VALUE",
        type=float,
        default=argparse.SUPPRESS,
        help=help_with_default(
            population_prevalence_help,
            module.population_prevalence,
        ),
    )
    if pipeline_mode:
        prev.add_argument(
            "--samp-prev-warning-threshold",
            dest="ldsc_sample_prevalence_warning_threshold",
            metavar="VALUE",
            type=float,
            default=argparse.SUPPRESS,
            help=help_with_default(
                "Absolute difference above which an explicit sample prevalence "
                "and the formatter-derived value produce a warning; values are "
                "expressed as proportions",
                module.sample_prevalence_comparison.warning_absolute_difference,
            ),
        )

    regression = parser.add_argument_group("LDSC heritability settings")
    regression.add_argument(
        "--intercept-h2",
        dest="ldsc_intercept",
        metavar="VALUE",
        type=float,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Constrain the single-trait LDSC intercept to this value",
            module.intercept,
        ),
    )
    regression.add_argument(
        "--two-step",
        dest="ldsc_two_step",
        metavar="CHI_SQUARE",
        type=float,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Two-step estimator chi-square cutoff; when unset, LDSC applies "
            "its data-dependent single-annotation behaviour",
            module.two_step,
        ),
    )
    regression.add_argument(
        "--chisq-max",
        dest="ldsc_chisq_max",
        metavar="VALUE",
        type=float,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Maximum chi-square statistic retained by LDSC",
            module.chisq_max,
        ),
    )
    regression.add_argument(
        "--n-blocks",
        dest="ldsc_n_blocks",
        metavar="N",
        type=int,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Number of block-jackknife blocks",
            module.n_blocks,
        ),
    )
    not_m_5_50_default = not module.use_m_5_50
    regression.add_argument(
        "--not-M-5-50",
        dest="ldsc_use_m_5_50",
        action="store_false",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Use .l2.M files instead of .l2.M_5_50 files, matching the "
            "original LDSC flag; this switch is %s by default"
            % ("used" if not_m_5_50_default else "not used"),
            not_m_5_50_default,
        ),
    )
    regression.add_argument(
        "--print-cov",
        dest="ldsc_print_covariance",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Write the covariance matrix of LDSC estimates",
            module.print_covariance,
        ),
    )
    regression.add_argument(
        "--print-delete-vals",
        dest="ldsc_print_delete_values",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Write block-jackknife delete values",
            module.print_delete_values,
        ),
    )

    mung = parser.add_argument_group("Munge-sumstats settings")

    mung.add_argument(
        "--info-min", "--ldsc-minimum-info",
        dest="ldsc_minimum_info",
        metavar="VALUE",
        type=float,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Minimum INFO score during munge_sumstats",
            module.minimum_info,
        ),
    )

    mung.add_argument(
        "--maf-min", "--ldsc-minimum-maf",
        dest="ldsc_minimum_maf",
        metavar="VALUE",
        type=float,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Minimum MAF during munge_sumstats",
            module.minimum_maf,
        ),
    )
    mung.add_argument(
        "--n-min",
        dest="ldsc_minimum_n",
        metavar="N",
        type=float,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Minimum per-variant sample size; when unset, upstream LDSC uses "
            "the 90th percentile divided by 1.5",
            module.minimum_n,
        ),
    )
    mung.add_argument(
        "--chunksize",
        dest="ldsc_chunksize",
        metavar="N",
        type=int,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Rows read per munge_sumstats chunk",
            module.chunksize,
        ),
    )
    mung.add_argument(
        "--keep-maf",
        dest="ldsc_keep_maf",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Retain the frequency column in munged statistics",
            module.keep_maf,
        ),
    )

    return parser


def get_ldsc_pipeline_parser(
    add_help: bool = False,
) -> argparse.ArgumentParser:
    """Return LDSC options with pipeline-specific prevalence guidance."""
    return get_ldsc_common_parser(add_help=add_help, pipeline_mode=True)



def get_assoc_plot_parser(add_help=False):
    """
    Parent-safe parser containing common assoc-plot arguments
    shared between DIRECT and PIPELINE modes. Every displayed
    default is loaded from canonical ``modules.manhattan`` YAML.
    """
    from postgwas.config import load_configuration
    from postgwas.core.ui import help_with_default

    defaults = load_configuration().modules.manhattan
    parser = argparse.ArgumentParser(add_help=add_help)
    # ------------------------------------------------------------
    # OUTPUTS
    # ------------------------------------------------------------
    output_group = parser.add_argument_group("Plot output").add_mutually_exclusive_group()

    output_group.add_argument(
        "--pdf",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help="Output PDF file (e.g., result.pdf)."
    )

    output_group.add_argument(
        "--png",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help="Output PNG file (e.g., result.png)."
    )
    # ------------------------------------------------------------
    # INPUTS
    # ------------------------------------------------------------
    input_group = parser.add_argument_group("Plot input")

    input_group.add_argument(
        "--pheno",
        metavar="NAME",
        default=argparse.SUPPRESS,
        help="Phenotype name to extract from GWAS-VCF file."
    )
    # ------------------------------------------------------------
    # FLAGS
    # ------------------------------------------------------------
    flag_group = parser.add_argument_group("Flags")

    flag_group.add_argument(
        "--allelic-shift",
        action="store_true",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Input VCF contains allelic-shift annotation",
            defaults.allelic_shift,
        ),
    )

    flag_group.add_argument(
        "--csq",
        action="store_true",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Flag coding variants from CSQ annotation in red",
            defaults.flag_coding,
        ),
    )
    # ------------------------------------------------------------
    # NUMERIC OPTIONS
    # ------------------------------------------------------------
    numeric_group = parser.add_argument_group("Plot thresholds")
    numeric_group.add_argument(
        "--nauto",
        type=int,
        default=argparse.SUPPRESS,
        metavar="N",
        help=help_with_default(
            "Number of autosomes", defaults.autosome_count,
        ),
    )

    numeric_group.add_argument(
        "--min-af",
        type=float,
        default=argparse.SUPPRESS,
        metavar="AF",
        help=help_with_default(
            "Minimum allele frequency filter", defaults.minimum_af,
        ),
    )

    numeric_group.add_argument(
        "--min-lp",
        type=float,
        default=argparse.SUPPRESS,
        metavar="LP",
        help=help_with_default(
            "Minimum −log10(P) threshold", defaults.minimum_neglog10_p,
        ),
    )

    numeric_group.add_argument(
        "--loglog-pval",
        type=float,
        default=argparse.SUPPRESS,
        metavar="LP",
        help=help_with_default(
            "−log10(P) threshold for switching to log-log scale",
            defaults.loglog_pvalue,
        ),
    )

    numeric_group.add_argument(
        "--cyto-ratio",
        type=float,
        default=argparse.SUPPRESS,
        metavar="RATIO",
        help=help_with_default(
            "Plot-height to cytoband-height ratio", defaults.cytoband_ratio,
        ),
    )

    numeric_group.add_argument(
        "--max-height",
        type=float,
        default=argparse.SUPPRESS,
        metavar="VALUE",
        help=help_with_default(
            "Maximum raw -log10(P), before the log-log display transformation", defaults.maximum_height,
        ),
    )

    numeric_group.add_argument(
        "--spacing",
        type=int,
        default=argparse.SUPPRESS,
        metavar="VALUE",
        help=help_with_default(
            "Spacing between chromosomes", defaults.chromosome_spacing,
        ),
    )
    # ------------------------------------------------------------
    # FIGURE SIZE
    # ------------------------------------------------------------
    fig_group = parser.add_argument_group("Figure dimensions")

    fig_group.add_argument(
        "--width",
        type=float,
        default=argparse.SUPPRESS,
        metavar="INCHES",
        help=help_with_default(
            "Plot width in inches", defaults.width,
        ),
    )

    fig_group.add_argument(
        "--height",
        type=float,
        default=argparse.SUPPRESS,
        metavar="INCHES",
        help=help_with_default(
            "Plot height in inches", defaults.height,
        ),
    )

    fig_group.add_argument(
        "--fontsize",
        type=int,
        default=argparse.SUPPRESS,
        metavar="POINTS",
        help=help_with_default(
            "Font size", defaults.font_size,
        ),
    )

    return parser
