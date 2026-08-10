import argparse
from rich_argparse import RichHelpFormatter
from rich_argparse import RawTextRichHelpFormatter

from postgwas.core.execution.runtime import validate_alphanumeric, validate_path, validate_prefix_files


# get_inputvcf_parser()
# get_genome_build_parser()
# get_common_out_parser()
# get_magma_binary_parser(add_help=False)
# get_flames_common_parser(add_help=False)
# get_common_imputation_parser


# ==================================================================
# ARGUMENT PARSER BUILDER (BASE + DIRECT ONLY)
# ==================================================================
def sumstat_summary_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    # 💥 DEFINE ARGUMENT GROUP for Step 2 💥
    group = parser.add_argument_group("GWAS-VCF QC settings")

    group.add_argument(
        "--reference-af-column",
        dest="reference_af_column",
        metavar="TAG",
        default="EUR",
        help=(
            "Name of the INFO/<AF> tag used as external reference AF.\n"
            "Examples: EUR, AFR, EAS.\n"
            "Default: EUR"
        ),
    )

    group.add_argument(
        "--maximum-af-difference",
        dest="maximum_af_difference",
        type=float,
        metavar="VALUE",
        default=0.2,
        help=(
            "Maximum allowed absolute difference between FORMAT/AF and "
            "INFO/<external_af_name>.\n"
            "Default: 0.2"
        ),
    )

    return parser




# --- 1. BASE PARSER FOR STEP 2 ARGUMENTS ---

def get_annot_ldblock_parser(add_help=False):
    """
    Returns the core ArgumentParser for annot_ldblock arguments, structured using a group.
    """
    # Create bare parser for inheritance
    parser = argparse.ArgumentParser(add_help=add_help)

    # 💥 DEFINE ARGUMENT GROUP for Step 2 💥
    group = parser.add_argument_group("LD-block annotation")

    group.add_argument(
        "--ld-block-populations",
        dest="ld_block_populations",
        nargs="+",
        metavar="POPULATION",
        default=["EUR", "AFR", "EAS"],
        help=(
            "Populations to annotate (e.g., EUR AFR). "
            "Files must be named following the pattern: [cyan]Genomeversion_Population_ldetect.bed.gz ;Eg. GRCh37_EUR_ldetect.bed.gz [/cyan]."
            "[bold green]Default:[/bold green] [cyan]EUR AFR EAS[/cyan] "
        )
    )
    group.add_argument(
        "--ld-region-dir",
        type=validate_path(must_exist=True, must_be_dir=True,dir_must_have_files=True),
        metavar="PATH",
        help=(
            "[bold bright_red]Required[/bold bright_red]: Directory containing LD-block BED files. "
            "Each BED file must contain four columns: CHROM, START, END, and Annotation. "
            "The fourth column is used as the LD-block annotation label."
            "The Genome build of bed fle should be same as input sumstat Genome build"
        ),
    )


    return parser




# ------------------------------------------------------------
# Formatter-specific parser
# ------------------------------------------------------------
def get_formatter_parser(*, direct_controls=False):
    from typing import get_args

    from postgwas.config import load_configuration
    from postgwas.core.variant_identifiers import VariantIdentifierType
    from postgwas.core.ui import help_with_default
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
            "One or more downstream inputs to create. Available values:\n%s\n"
            "If omitted, modules.formatting.formats is read from --run-config."
            % "\n".join(
                "  %-16s %s" % (target, FORMAT_CONTRACTS[target].name)
                for target in available
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
            label="Configured default",
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
        group.add_argument(
            "--resume",
            action=argparse.BooleanOptionalAction,
            default=argparse.SUPPRESS,
            help=help_with_default(
                "Reuse existing formatter outputs only after their provenance and "
                "content have been validated",
                configuration.run.resume,
            ),
        )
        group.add_argument(
            "--overwrite",
            action="store_true",
            default=argparse.SUPPRESS,
            help=help_with_default(
                "Rerun formatting and replace existing outputs; this takes "
                "precedence over resume",
                configuration.run.overwrite,
            ),
        )

    return parser




def get_ld_clump_parser(add_help=False):
    """
    Returns a parent-safe parser containing arguments
    that are shared between DIRECT and PIPELINE modes.

    COMMON arguments:
        --ld-mode
        --population
        --r2-cutoff
        --window-kb

    (VCF input, outdir, resource folder, bcftools, harmonisation,
     annot_ldblock etc. are added separately in main() via their own
     get_*_parser functions.)
    """
    parser = argparse.ArgumentParser(add_help=add_help)

    # -------------------------------
    # General LD clumping options
    # -------------------------------
    general = parser.add_argument_group("LD-clumping settings")

    # general.add_argument(
    #     "--ld-mode",
    #     default="by_regions",
    #     choices=["by_regions", "standard"],
    #     metavar=" ",
    #     help=(
    #         "LD clumping strategy.\n"
    #         "Choices: [cyan]by_regions[/cyan], [cyan]standard[/cyan]\n"
    #         f"[bold green]Default:[/bold green] [cyan]by_regions[/cyan]"
    #     )
    # )

    # -------------------------------
    # Standard-mode arguments
    # -------------------------------
    standard = parser.add_argument_group(
        "PLINK-style clumping"
    )

    # Required Inputs
    standard.add_argument("--ld-folder", metavar="PATH",
                     type=validate_path(must_exist=True,must_not_be_empty=True),
                     help="Folder containing EUR_chr*.ld.gz files")

    # Thresholds
    standard.add_argument("--lead-p", metavar="P", type=float, default=5e-8,
                          help="P-value threshold for Lead SNPs [bold green]Default:[/bold green] [cyan]5e-8[/cyan]")
    standard.add_argument("--r2-clump", metavar="R2", type=float, default=0.6,
                          help="The minimum r2 for defining independent significant SNPs. [bold green]Default:[/bold green] [cyan]0.6[/cyan]")
    standard.add_argument("--r2-lead", metavar="R2", type=float, default=0.1,
                          help="The minimum r2 for defining lead SNPs, which is used for the second clumping (clumping of the independent significant SNPs). [bold green]Default:[/bold green] [cyan]0.1[/cyan]")

    # System & Output
    standard.add_argument("--merge-dist", metavar="BP", type=int, default=250000,
                          help="The maximum distance between LD blocks of independent significant SNPs to merge into a single genomic locus. [bold green]Default:[/bold green] [cyan]250000[/cyan]")


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
    available_builds = list(available_builds)
    if not available_builds:
        raise ValueError("At least one configured genome build is required")
    resolved_default = default_build or available_builds[0]
    if resolved_default not in available_builds:
        raise ValueError("The default genome build must be one of available_builds")
    parser = argparse.ArgumentParser(
        add_help=False,
        formatter_class=RichHelpFormatter
    )

    input_group = parser.add_argument_group("Genome coordinates")

    # -------- GENOME VERSION --------
    input_group.add_argument(
        "--genome-build",
        default=argparse.SUPPRESS if suppress_default else resolved_default,
        choices=available_builds,
        metavar="BUILD",
        help=(
            "Genome build used by the input VCF and all reference files "
            "(default: %s). Do not mix resources from different builds."
            % resolved_default
        )
    )
    return parser


def get_population_parser(add_help=False):
    """Return the shared reference-population option."""
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("Reference population")
    group.add_argument(
        "--population",
        choices=["EUR", "AFR", "EAS", "SAS", "AMR"],
        default="EUR",
        metavar="CODE",
        help="Population represented by the reference data (default: EUR).",
    )
    return parser

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


def get_common_magma_covar_parser(add_help=False):
    from textwrap import wrap
    from typing import get_args

    from postgwas.config import load_configuration
    from postgwas.config.models.modules.magmacovar import MagmaCovarDirection
    from postgwas.core.ui import help_with_default

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
            "separately. See 'How to choose a MAGMA model' and 'Advanced MAGMA "
            "model modifiers' below" % accepted_modifiers,
            " ".join(defaults.model) if defaults.model else "none",
            label="Configured default",
        ),
    )

    group.add_argument(
        "--covariate-direction",
        choices=get_args(MagmaCovarDirection),
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Testing direction passed through MAGMA's direction-covar model "
            "modifier. Greater and smaller are one-sided tests. Two-sided is "
            "MAGMA's default for any gene property. The original FLAMES README "
            "relies on that default, while tissue-specific analysis uses greater",
            defaults.direction,
            label="Configured default",
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
            label="Configured default",
        ),
    )

    choice_lines = ["Choose the model that matches your scientific question:"]
    for use_case in defaults.model_use_cases.values():
        selected_model = (
            "--covariate-model " + " ".join(use_case.model)
            if use_case.model else "omit --covariate-model"
        )
        choice_lines.extend((
            "",
            "  %s" % use_case.label,
            "    Question: %s" % use_case.question,
            "    Use: %s" % selected_model,
            "    Direction: --covariate-direction %s" % use_case.direction,
        ))
    parser.add_argument_group(
        "How to choose a MAGMA model",
        "\n".join(choice_lines),
    )

    def append_help_detail(lines, label, value):
        prefix = "  %s: " % label
        continuation = " " * len(prefix)
        chunks = wrap(
            value,
            width=76 - len(prefix),
            break_long_words=False,
            break_on_hyphens=False,
        )
        lines.append(prefix + chunks[0])
        lines.extend(continuation + chunk for chunk in chunks[1:])

    modifier_lines = [
        "Use these only for a planned conditional or multivariable analysis:",
        "",
        "Values used below:",
    ]
    for placeholder, description in defaults.model_placeholders.items():
        append_help_detail(modifier_lines, "<%s>" % placeholder, description)
    modifier_lines.extend((
        "",
        "Example property names are illustrative; use exact headers from your table.",
    ))
    for specification in defaults.model_modifiers.values():
        modifier_lines.append("")
        append_help_detail(modifier_lines, "Purpose", specification.description)
        append_help_detail(
            modifier_lines,
            "Published GWAS example",
            specification.published_example,
        )
        append_help_detail(modifier_lines, "Example command", specification.example)
        modifier_lines.append("  Accepted forms:")
        modifier_lines.extend("    %s" % syntax for syntax in specification.syntax)
    modifier_lines.extend((
        "",
        "Interaction and gene-set-only modifiers are unavailable.",
        "This interface does not provide --set-annot.",
    ))
    parser.add_argument_group(
        "Advanced MAGMA model modifiers",
        "\n".join(modifier_lines),
    )
    return parser




def get_common_magma_assoc_parser(add_help=False):
    """
    Parser containing only MAGMA association-related arguments.
    This is parent-safe and can be inserted into direct/pipeline mode.
    """
    from postgwas.config import load_configuration
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
        help=(
            "[bold bright_red]Required[/bold bright_red]: MAGMA gene location file (.loc) matching the genome build. "
            "Must contain columns in the following order: Gene ID, Chromosome, Gene Start, Gene End. "
            "Gene IDs must match those used in pathway files and gene-covariate files. "
            "[bold green]Default:[/bold green] [cyan]None[/cyan]"
        ),
    )

    group.add_argument(
        "--gene-set-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        type=validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True),
        help=(
            "Pathway file (.gmt). MSigDB or custom GMT supported."
            f"\n[bold green]Default:[/bold green] [cyan]None[/cyan]"
        ),
    )

    group.add_argument(
        "--window-upstream",
        type=int,
        default=argparse.SUPPRESS,
        metavar="KB",
        help=(
            "Upstream gene window in kb. "
            f"[bold green]Configured default:[/bold green] [cyan]{defaults.gene_window_upstream_kb}[/cyan]"
        ),
    )

    group.add_argument(
        "--window-downstream",
        type=int,
        default=argparse.SUPPRESS,
        metavar="KB",
        help=(
            "Downstream gene window in kb. "
            f"[bold green]Configured default:[/bold green] [cyan]{defaults.gene_window_downstream_kb}[/cyan]"
        ),
    )

    group.add_argument(
        "--gene-model",
        default=argparse.SUPPRESS,
        metavar="MODEL",
        help=(
            "MAGMA model supported with SNP p-value input. "
            f"[bold green]Configured default:[/bold green] [cyan]{defaults.gene_model}[/cyan]"
        ),
    )

    group.add_argument(
        "--sample-size-column",
        dest="sample_size_column",
        default=argparse.SUPPRESS,
        metavar="COLUMN",
        help=(
            "Column name for sample sizes in *.pval input. "
            f"[bold green]Configured default:[/bold green] [cyan]{defaults.input.sample_size_column}[/cyan]"
        ),
    )

    return parser


def get_magma_binary_parser(add_help=False):
    """
    Parser containing only MAGMA association-related arguments.
    This is parent-safe and can be inserted into direct/pipeline mode.
    """
    from postgwas.config import load_configuration

    configured = load_configuration().resources.executables.magma
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("MAGMA binary")

    # Required inputs (but NOT marked required here → handled by subparser)
    group.add_argument(
        "--magma",
        default=argparse.SUPPRESS,
        metavar="PATH",
        help=(
            "Path to the MAGMA binary. "
            "If not provided, assumes 'magma' is available in your system PATH. "
            f"[bold green]Configured default:[/bold green] [cyan]{configured}[/cyan]"
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
            label="Configured default",
        ),
    )
    return parser




def get_bcftools_binary_parser(add_help=False):
    """Return the shared bcftools executable option."""
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("External software")

    # Required inputs (but NOT marked required here → handled by subparser)
    group.add_argument(
        "--bcftools",
        default="bcftools",
        metavar="PATH",
        help="Location of bcftools. Omit this when bcftools is available from your terminal.",
    )
    return parser





def get_common_pops_parser(add_help=False):
    """Return reusable PoPS options whose defaults come only from YAML."""
    from typing import get_args

    from postgwas.config import load_configuration
    from postgwas.config.models.modules.pops import PopsMethod
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
            label="Configured default",
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
            label="Configured default",
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
            label="Configured default",
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
            label="Configured default",
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
            label="Configured default",
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
        help="Custom target-score table used only when no MAGMA prefix is supplied.",
    )
    grp_cov.add_argument(
        "--target-covariates-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help="Optional covariate table aligned to the custom target scores.",
    )
    grp_cov.add_argument(
        "--target-error-covariance-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help="Optional SciPy NPZ or NumPy NPY covariance for custom target scores.",
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
            label="Configured default",
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
            label="Configured default",
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
            label="Configured default",
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
            label="Configured default",
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
            label="Configured default",
        ),
    )
    grp_misc.add_argument(
        "--save-matrix-files",
        action="store_true",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Write upstream PoPS training and prediction matrices",
            module.save_matrix_files,
            label="Configured default",
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
    parser = argparse.ArgumentParser(formatter_class=RawTextRichHelpFormatter,add_help=add_help)

    group = parser.add_argument_group("Quality-control filters")

    # =====================================================
    # NUMERIC FILTERS
    # =====================================================
    group.add_argument(
        "--minimum-neglog10-p",
        dest="minimum_neglog10_p",
        type=float,
        metavar="VALUE",
        default=None,
        help=(
            "Keep variants with -log10(P) at or above this value. "
            "For genome-wide significance (P <= 5e-8), use 7.3. "
            "Omit to keep all P-values."
        )
    )

    group.add_argument(
        "--minimum-maf",
        dest="minimum_maf",
        type=float,
        metavar="VALUE",
        default=0.002,
        help="Remove variants below this minor-allele frequency (default: 0.002)."
    )

    group.add_argument(
        "--reference-af-column",
        dest="reference_af_column",
        default="EUR",
        metavar="TAG",
        help="Reference-population frequency tag in the VCF, such as EUR, AFR, or EAS (default: EUR)."
    )

    group.add_argument(
        "--maximum-af-difference",
        dest="maximum_af_difference",
        type=float,
        metavar="VALUE",
        default=0.2,
        help="Remove variants whose study and reference allele frequencies differ by more than this value (default: 0.2)."
    )

    group.add_argument(
        "--minimum-info",
        dest="minimum_info",
        type=float,
        metavar="VALUE",
        default=0.3,
        help="Minimum INFO score (default: 0.3)."
    )

    group.add_argument(
        "--maximum-info",
        dest="maximum_info",
        type=float,
        metavar="VALUE",
        default=None,
        help="Remove variants above this INFO score. Omit to set no maximum."
    )

    group.add_argument(
        "--missing-info-action",
        dest="missing_info_action",
        type=str,
        metavar="ACTION",
        default="remove",
        choices=["keep", "remove"],
        help="What to do when INFO is missing: remove or keep (default: remove)."
    )


    # =====================================================
    # FLAGS
    # =====================================================
    group.add_argument(
        "--include-indels",
        action="store_true",
        help="Keep insertion/deletion variants. By default, only SNPs are kept."
    )

    group.add_argument(
        "--remove-palindromic",
        action="store_true",
        help=(
            "Remove ambiguous A/T and C/G SNPs when their allele frequency falls "
            "between the two palindromic AF bounds."
        )
    )



    # =====================================================
    # PALINDROMIC AF BOUNDS
    # =====================================================
    group.add_argument(
        "--palindromic-af-lower",
        type=float,
        default=0.4,
        metavar="VALUE",
        help="Lower ambiguity bound used with --remove-palindromic (default: 0.4)."
    )

    group.add_argument(
        "--palindromic-af-upper",
        type=float,
        default=0.6,
        metavar="VALUE",
        help="Upper ambiguity bound used with --remove-palindromic (default: 0.6)."
    )

    # =====================================================
    # MHC REGION
    # =====================================================
    group.add_argument(
        "--remove-mhc",
        action="store_true",
        default=False,
        help="Remove variants in the MHC region. By default, the region is retained."
    )

    group.add_argument(
        "--mhc-chrom",
        default="6",
        metavar="CHROM",
        help="Chromosome containing the MHC region (default: 6)."
    )

    group.add_argument(
        "--mhc-start",
        type=int,
        default=25_000_000,
        metavar="POSITION",
        help="Start coordinate of the MHC region (default: 25,000,000)."
    )

    group.add_argument(
        "--mhc-end",
        type=int,
        default=34_000_000,
        metavar="POSITION",
        help="End coordinate of the MHC region (default: 34,000,000)."
    )

    return parser




# ============================================================
# COMMON FLAMES PARSER (shared by DIRECT + PIPELINE)
# ============================================================
def get_flames_common_parser(add_help=False):
    """Return YAML-backed FLAMES options shared by direct and pipeline modes."""
    from postgwas.config import load_configuration
    from postgwas.core.ui import help_with_default

    defaults = load_configuration()
    module = defaults.modules.flames
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("FLAMES scientific settings")
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
            label="Configured default",
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
            label="Configured default",
        ),
    )
    group.add_argument(
        "--vep-command",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help="VEP executable required when --flames-vep-mode local is selected.",
    )
    group.add_argument(
        "--vep-cache",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help="Build-compatible VEP cache directory required for local VEP.",
    )
    group.add_argument(
        "--flames-cadd-mode",
        choices=("api", "local"),
        metavar="MODE",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Use the upstream CADD API or a configured local tabix-indexed CADD file",
            module.cadd_mode,
            label="Configured default",
        ),
    )
    group.add_argument(
        "--cadd-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help="Genome-build-matched CADD file required for local CADD mode.",
    )
    group.add_argument(
        "--tabix",
        metavar="COMMAND",
        default=argparse.SUPPRESS,
        help="Tabix executable override for local CADD mode.",
    )
    return parser



def get_common_imputation_parser(add_help: bool = False) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=add_help)
    grp = parser.add_argument_group("Imputation settings")

    grp.add_argument(
        "--imputation-engine",
        dest="imputation_engine",
        metavar="ENGINE",
        type=str,
        choices=["pred_ld"],
        default="pred_ld",
        help=(
            "Choose imputation tool. "
            "[bold]Available options:[/bold] [cyan]pred_ld[/cyan] "
            "[bold green]Default:[/bold green] [cyan]pred_ld[/cyan]"
        ),
    )

    grp.add_argument(
        "--imputation-ld-reference",
        dest="imputation_ld_reference",
        metavar="PATH",
        type=validate_path(must_exist=True, must_be_dir=True),
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
        default=0.8,
        help=("LD r² tagging threshold."
              "[bold green]Default:[/bold green] [cyan]0.8[/cyan]")
    )

    grp.add_argument(
        "--imputation-minimum-maf",
        dest="imputation_minimum_maf",
        metavar="MAF",
        type=float,
        default=0.001,
        help=("Minor allele frequency threshold."
              "[bold green]Default:[/bold green] [cyan]0.001[/cyan]")
    )

    grp.add_argument(
        "--ref",
        metavar="NAME",
        type=str,
        default="TOP_LD",
        help=("Reference panel name label."
              "[bold green]Default:[/bold green] [cyan]TOP_LD[/cyan]")
    )

    grp.add_argument(
        "--correlation-method",
        metavar="METHOD",
        type=str,
        default="pearson",
        choices=["pearson", "spearman"],
        help=("Correlation method used in QC/post-processing."
              "[bold green]Default:[/bold green] [cyan]pearson[/cyan]")

    )

    grp.add_argument(
        "--resource-directory",
        metavar="PATH",
        default=None,
        type=validate_path(must_exist=True, must_be_dir=True),
        help="[bold bright_red]Required[/bold bright_red]: PostGWAS resource directory.",
    )

    return parser



def get_ldsc_common_parser(add_help: bool = False) -> argparse.ArgumentParser:
    """
    Arguments shared by LDSC DIRECT and PIPELINE modes.

    DIRECT mode:
        --sumstats      → munged .sumstats.gz (LDSC-ready)
        --ref-ld-chr    → eur_w_ld_chr/ (directory)
        --w-ld-chr      → eur_w_ld_chr/ (directory)
        --samp-prev / --pop-prev

    PIPELINE mode:
        --sumstats      → raw LDSC-input TSV (your vcf_to_ldsc output)
        --merge-alleles → HapMap3 SNP+allele list (w_hm3.snplist)
        --info-min / --maf-min used during munge_sumstats
        (ref/w ld-chr + prevs also used)

    Tools mode:
        Gets the same arguments; backend can decide what to use.
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    grp = parser.add_argument_group("Allele matching")

    grp.add_argument(
        "--merge-alleles",
        metavar="PATH",
        type=validate_path(must_exist=True, must_be_file=True),
        help=(
            "[bold bright_red]Required[/bold bright_red]: HapMap3 SNP list with alleles (w_hm3.snplist). Used during munge_sumstats.\n"
            "Optional in DIRECT mode; recommended in PIPELINE mode."
        ),
    )

    ref = parser.add_argument_group("LDSC reference")

    ref.add_argument(
        "--ref-ld-chr",
        metavar="PATH",
        type=validate_path(must_exist=True, must_be_dir=True),
        help="[bold bright_red]Required[/bold bright_red]: Directory containing reference LD scores (e.g., eur_w_ld_chr/).",
    )

    ref.add_argument(
        "--w-ld-chr",
        metavar="PATH",
        type=validate_path(must_exist=True, must_be_dir=True),
        help="[bold bright_red]Required[/bold bright_red]: Directory containing LDSC regression weights (often same as --ref-ld-chr).",
    )

    prev = parser.add_argument_group("Liability-scale reporting")

    prev.add_argument(
        "--samp-prev",
        metavar="VALUE",
        type=float,
        default=None,
        help=("Sample prevalence (cases / total) used for liability h²."
              "[bold green]Default:[/bold green] [cyan]None[/cyan]" )

    )

    prev.add_argument(
        "--pop-prev",
        metavar="VALUE",
        type=float,
        default=None,
        help=("Population prevalence (e.g., 0.01)."
              "[bold green]Default:[/bold green] [cyan]None[/cyan]")
    )

    mung = parser.add_argument_group("LDSC formatting filters")

    mung.add_argument(
        "--ldsc-minimum-info",
        dest="ldsc_minimum_info",
        metavar="VALUE",
        type=float,
        default=0.9,
        help=("Minimum INFO score during munge_sumstats. "
              "[bold green]Default:[/bold green] [cyan]0.9[/cyan]" )
    )

    mung.add_argument(
        "--ldsc-minimum-maf",
        dest="ldsc_minimum_maf",
        metavar="VALUE",
        type=float,
        default=0.01,
        help=("Minimum MAF during munge_sumstats (default: 0.01). [bold green]Default:[/bold green] [cyan]0.01[/cyan]" )
    )

    # docker_grp = parser.add_argument_group("LDSC Docker Execution ")

    # docker_grp.add_argument(
    #     "--docker-image",
    #     metavar=" ",
    #     type=str,
    #     default="jibinjv/ldsc:1.0.1",
    #     help=("Docker image used to run LDSC. [bold green]Default:[/bold green] [cyan]0.7[/cyan]"),
    # )

    # docker_grp.add_argument(
    #     "--platform",
    #     metavar=" ",
    #     type=str,
    #     default="linux/amd64",
    #     help="Docker platform (default: linux/amd64, for Apple Silicon + x86 image). [bold green]Default:[/bold green] [cyan]linux/amd64[/cyan]",
    # )

    return parser



def get_assoc_plot_parser(add_help=False):
    """
    Parent-safe parser containing common assoc-plot arguments
    shared between DIRECT and PIPELINE modes.
    Matches the original assoc_plot.R defaults:
      • min-lp default = 2
      • spacing default = 20 (R default)
      • cyto-ratio default = 25
      • loglog-pval default = 10
      • width = 7.0
      • fontsize = 12
      • nauto = 22
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    # ------------------------------------------------------------
    # OUTPUTS
    # ------------------------------------------------------------
    output_group = parser.add_argument_group("Plot output")

    output_group.add_argument(
        "--pdf",
        metavar="PATH",
        help="Output PDF file (e.g., result.pdf)."
    )

    output_group.add_argument(
        "--png",
        metavar="PATH",
        help="Output PNG file (e.g., result.png)."
    )
    # ------------------------------------------------------------
    # INPUTS
    # ------------------------------------------------------------
    input_group = parser.add_argument_group("Plot input")

    input_group.add_argument(
        "--pheno",
        metavar="NAME",
        help="Phenotype name to extract from GWAS-VCF file."
    )
    # ------------------------------------------------------------
    # FLAGS
    # ------------------------------------------------------------
    flag_group = parser.add_argument_group("Flags")

    flag_group.add_argument(
        "--allelic-shift",
        action="store_true",
        help="Input VCF contains allelic-shift annotation."
    )

    flag_group.add_argument(
        "--csq",
        action="store_true",
        help="Flag coding variants (from CSQ annotation) in red."
    )
    # ------------------------------------------------------------
    # NUMERIC OPTIONS
    # ------------------------------------------------------------
    numeric_group = parser.add_argument_group("Plot thresholds")
    numeric_group.add_argument(
        "--nauto",
        type=int,
        default=22,
        metavar="N",
        help="Number of autosomes. [bold green]Default:[/bold green] [cyan]22[/cyan]"
    )

    numeric_group.add_argument(
        "--min-af",
        type=float,
        default=0.0,
        metavar="AF",
        help="Minimum allele frequency filter. [bold green]Default:[/bold green] [cyan]0.0[/cyan]"
    )

    numeric_group.add_argument(
        "--min-lp",
        type=int,
        default=2,
        metavar="LP",
        help="Minimum −log10(P) threshold. [bold green]Default:[/bold green] [cyan]2[/cyan]"
    )

    numeric_group.add_argument(
        "--loglog-pval",
        type=int,
        default=10,
        metavar="LP",
        help="-log10(P) threshold for switching to log-log scale. "
             "[bold green]Default:[/bold green] [cyan]10[/cyan]"
    )

    numeric_group.add_argument(
        "--cyto-ratio",
        type=int,
        default=25,
        metavar="RATIO",
        help="Plot height : cytoband height ratio. "
             "[bold green]Default:[/bold green] [cyan]25[/cyan]"
    )

    numeric_group.add_argument(
        "--max-height",
        type=int,
        metavar="VALUE",
        help="Maximum vertical height of plot (optional)."
    )

    numeric_group.add_argument(
        "--spacing",
        type=int,
        default=20,
        metavar="VALUE",
        help="Spacing between chromosomes. "
             "[bold green]Default:[/bold green] [cyan]20[/cyan]"
    )
    # ------------------------------------------------------------
    # FIGURE SIZE
    # ------------------------------------------------------------
    fig_group = parser.add_argument_group("Figure dimensions")

    fig_group.add_argument(
        "--width",
        type=float,
        default=7.0,
        metavar="INCHES",
        help="Plot width in inches. [bold green]Default:[/bold green] [cyan]7.0[/cyan]"
    )

    fig_group.add_argument(
        "--height",
        type=float,
        metavar="INCHES",
        help="Plot height in inches (optional)."
    )

    fig_group.add_argument(
        "--fontsize",
        type=int,
        default=12,
        metavar="POINTS",
        help="Font size. [bold green]Default:[/bold green] [cyan]12[/cyan]"
    )

    return parser
