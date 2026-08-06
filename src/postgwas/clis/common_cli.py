import argparse
import multiprocessing
import sys
from rich_argparse import RichHelpFormatter 
from rich_argparse import RawTextRichHelpFormatter
from postgwas.utils.main import validate_path
import textwrap
from pathlib import Path

from postgwas.utils.main import validate_path,detect_total_memory_gb,validate_alphanumeric,validate_prefix_files


# get_defaultresourse_parser()
# get_inputvcf_parser()
# get_genomeversion_parser()
# get_common_out_parser()
# get_magma_binary_parser(add_help=False)
# get_flames_common_parser(add_help=False)
# get_common_imputation_parser


TOTAL_CORES = multiprocessing.cpu_count()
AUTO_THREADS = max(1, TOTAL_CORES - 2)
TOTAL_MEM_GB = detect_total_memory_gb()
AUTO_MEM_GB = max(1, TOTAL_MEM_GB - 2)         # keep ~2GB for system
AUTO_MAX_MEM = f"{AUTO_MEM_GB}G"               # string format for CLI



# ==================================================================
# ARGUMENT PARSER BUILDER (BASE + DIRECT ONLY)
# ==================================================================
def sumstat_summary_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    # 💥 DEFINE ARGUMENT GROUP for Step 2 💥
    group = parser.add_argument_group('sumstat qc summary Arguments')

    group.add_argument(
        "--qc_summary_external_af_name",
        metavar=" ",
        default="EUR",
        help=(
            "Name of the INFO/<AF> tag used as external reference AF.\n"
            "Examples: EUR, AFR, EAS.\n"
            "Default: EUR"
        ),
    )

    group.add_argument(
        "--qc_summary_allelefreq-diff-cutoff",
        type=float,
        metavar=" ",
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
    group = parser.add_argument_group('LD BLOCK ANNOTATION Arguments')

    group.add_argument(
        "--ld_block_population", 
        nargs="+", 
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
        metavar='',
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
def get_formatter_parser():
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group("Formatter Options")

    group.add_argument(
        "--format",
        choices=["magma", "finemap_susie","finemap_finemap", "ldpred", "ldsc"],
        nargs="+",
        required=True,
        metavar=" ",
        help=(
            "[bold bright_red]Required[/bold bright_red]: Choose one or more "
            "output formats for conversion. choices are magma, finemap_susie,finemap_finemap,ldpred, ldsc.\n"
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
    general = parser.add_argument_group("Common LD Clumping Arguments")

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

    general.add_argument(
        "--ld_clump_population",
        default="EUR",
        choices=["EUR", "AFR", "EAS"],
        metavar=" ",
        help=(
            "Population LD blocks (by_regions mode only).\n"
            f"[bold green]Default:[/bold green] [cyan]EUR[/cyan]"
        )
    )

    # -------------------------------
    # Standard-mode arguments
    # -------------------------------
    standard = parser.add_argument_group(
        "Standard LD Clumping Arguments (PLINK-style)"
    )

    # Required Inputs
    standard.add_argument("--ld-folder", metavar="",
                     type=validate_path(must_exist=True,must_not_be_empty=True), 
                     help="Folder containing EUR_chr*.ld.gz files")
    
    # Thresholds
    standard.add_argument("--lead-p", metavar="", type=float, default=5e-8, 
                          help="P-value threshold for Lead SNPs [bold green]Default:[/bold green] [cyan]5e-8[/cyan]")
    standard.add_argument("--r2-clump", metavar="", type=float, default=0.6, 
                          help="The minimum r2 for defining independent significant SNPs. [bold green]Default:[/bold green] [cyan]0.6[/cyan]")
    standard.add_argument("--r2-lead", metavar="", type=float, default=0.1, 
                          help="The minimum r2 for defining lead SNPs, which is used for the second clumping (clumping of the independent significant SNPs). [bold green]Default:[/bold green] [cyan]0.1[/cyan]")

    # System & Output
    standard.add_argument("--merge-dist", metavar="", type=int, default=250000, 
                          help="The maximum distance between LD blocks of independent significant SNPs to merge into a single genomic locus. [bold green]Default:[/bold green] [cyan]250000[/cyan]")


    return parser





def safe_parse_args(top_parser):
    """
    A wrapper around argparse.parse_args() that prevents premature exit and
    reroutes invalid input to the appropriate subcommand help.
    """
    try:
        return top_parser.parse_args()

    except SystemExit:
        # argparse already printed an error message
        argv = sys.argv

        # If user typed: cli.py direct ...
        if len(argv) > 1 and argv[1] == "direct":
            top_parser.parse_args(["direct", "--help"])
            sys.exit(1)

        # If user typed: cli.py pipeline ...
        if len(argv) > 1 and argv[1] == "pipeline":
            top_parser.parse_args(["pipeline", "--help"])
            sys.exit(1)

        # If no mode provided → show top-level help
        top_parser.print_help()
        sys.exit(1)




def get_defaultresourse_parser():
    """
    Common parser providing system resource arguments (threads + memory).
    Safe for inheritance by other parsers.
    """
    parser = argparse.ArgumentParser(
        add_help=False,
        formatter_class=RichHelpFormatter
    )

    # Create a dedicated group to avoid showing under "Optional Arguments"
    resource_group = parser.add_argument_group("SYSTEM RESOURCES")

    # -------- THREADS ----------
    resource_group.add_argument(
        "--nthreads",
        type=int,
        metavar='',
        default=AUTO_THREADS,
        help=f"Threads to use [bold green]Default:[/bold green] [cyan]{AUTO_THREADS}[/cyan]"
    )

    # -------- MEMORY ----------
    resource_group.add_argument(
        "--max-mem",
        metavar='',
        default=AUTO_MAX_MEM,
        help=(
            f"Maximum memory allowed. "
            "Formats accepted: 4G, 800M, 1200M."
            f"[bold green]Default:[/bold green] [cyan]{AUTO_MAX_MEM}[/cyan]"
        )
    )

    resource_group.add_argument(
        "--seed",
        type=int,
        default=10,
        metavar="",
        help="Random seed. [bold green]Default:[/bold green] [cyan]10[/cyan]",
        
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

    input_group = parser.add_argument_group("INPUT sumstat Arguments")

    # -------- VCF INPUT --------
    input_group.add_argument(
        "--vcf",
        metavar=" ",
        type=validate_path(
            must_exist=True,
            must_be_file=True,
            must_not_be_empty=True,
            allowed_suffixes=[".vcf", ".vcf.gz"]
        ),
        help=(
            "[bold bright_red]Required[/bold bright_red]: Harmonised GWAS summary statistics "
            "VCF file generated by the [cyan]PostGWAS harmonisation module[/cyan].\n"
        )
    )

    return parser




def get_genomeversion_parser():
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

    input_group = parser.add_argument_group("Genome Version Arguments")

    # -------- GENOME VERSION --------
    input_group.add_argument(
        "--genome-version",
        default="GRCh37",
        choices=["GRCh37", "GRCh38"],
        metavar="",
        help=(
            "[bold bright_red]Required[/bold bright_red]: Genome build of the GWAS summary statistics "
            "([cyan]GRCh37[/cyan] or [cyan]GRCh38[/cyan]).\n"
            "This value [yellow]must match[/yellow] all genomic coordinate–based reference files used "
            "in downstream steps, including:\n"
            " • LD reference panels\n"
            " • LD block BED files\n"
            " • Gene location (.loc) files\n"
            " • Annotation/functional reference datasets (MAGMA, PoPS)\n\n"
            "[bold green]Default:[/bold green]: [cyan]GRCh37[/cyan]."
        )
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

    group = parser.add_argument_group("OUTPUT Arguments")

    # -------- SAMPLE ID --------
    group.add_argument(
        "--sample_id",
        metavar="",
        type=validate_alphanumeric,
        help=(
            "[bold bright_red]Required[/bold bright_red]: Unique identifier for this GWAS dataset.\n"
            "Use the phenotype/trait name (e.g., [cyan]PGC3_SCZ_eur[/cyan], [cyan]UKB_BMI[/cyan]).\n"
            "This value is used as the prefix for ALL generated output files."
        ),
    )

    # -------- OUTPUT DIRECTORY --------
    group.add_argument(
        "--outdir",
        metavar="",
        default=None,  # 🔑 SAFE
        type=validate_path(
            must_exist=False,
            must_be_dir=True,
            create_if_missing=True,
        ),
        help=(
            "[bold bright_red]Required[/bold bright_red]: Output directory where all results "
            "and intermediate files will be stored.\n"
            "[bold green]Default:[/bold green] current working directory"
        ),
    )

    return parser


def get_common_magma_covar_parser(add_help=False):
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("MAGMA covar Arguments")
    group.add_argument(
        "--covariates",
        metavar="",
        type=validate_path(must_exist=True,must_be_file=True,must_not_be_empty=True),
        help="[bold bright_red]Required[/bold bright_red]:  Covariate matrix file (TSV)"
    )
    
    group.add_argument(
        "--covar_model",
        metavar=" ",
        default=None,
        help=(
            "Specify a single MAGMA covariate model.\n"
            "Supported formats:\n"
            "  condition=<v1,v2,...>\n"
            "  condition-hide=<v1,v2,...>\n"
            "  condition-residualize=<v1,v2,...>\n"
            "  condition-interaction=<v1,v2,...>\n"
            "  interaction-each=<v1,v2,...>\n"
            "  interaction-all=<v1,v2,...>\n"
            "  analyse\n"
            "  joint\n"
            "  interaction\n\n"
            "Examples:\n"
            "  --covar_model condition-hide=Average\n"
            "  --covar_model condition=Brain_Cortex,Brain_Cerebellum\n\n"
            "NOTE for FLAMES users:\n"
            "  FLAMES typically DOES NOT recommend using a custom --covar_model\n"
            "  unless you know exactly what you're doing."
        ),
    )

    group.add_argument(
        "--covar_direction",
        metavar=" ",
        default=None,
        choices=["two-sided", "greater", "less"],
        help=(
            "Specify test direction.\n"
            "MAGMA does not take a separate --direction flag; instead, the\n"
            "direction token is passed as an additional argument after the\n"
            "model specification.\n\n"
            "Example CLI:\n"
            "--covar_model condition-hide=Average --covar_direction greater\n"
            "Internally MAGMA will see:\n"
            "--model condition-hide=Average direction=greater"
        ),
    )
    return parser



    
def get_common_magma_assoc_parser(add_help=False):
    """
    Parser containing only MAGMA association-related arguments.
    This is parent-safe and can be inserted into direct/pipeline mode.
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("MAGMA Association Arguments")

    # Required inputs (but NOT marked required here → handled by subparser)
    group.add_argument(
        "--ld_ref",
        type=lambda p: validate_prefix_files(p, [".bed", ".bim", ".fam"]),
        metavar="",
        help=(
        "[bold bright_red]Required[/bold bright_red]: Prefix of PLINK reference panel for MAGMA "
        "(e.g., 1000G_EUR).\n"
        "Should point to files: PREFIX.bed / PREFIX.bim / PREFIX.fam.\n"
        "Make sure the summary statistics and the LD reference "
        "use the same genome build."
        "\n[bold green]Default:[/bold green] [cyan]None[/cyan]"
        ) ,
    )

    group.add_argument(
        "--gene_loc_file",
        metavar="",
        type=validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True),
        help=(
            "[bold bright_red]Required[/bold bright_red]: MAGMA gene location file (.loc) matching the genome build. "
            "Must contain columns in the following order: Gene ID, Chromosome, Gene Start, Gene End. "
            "Gene IDs must match those used in pathway files and gene-covariate files. "
            "[bold green]Default:[/bold green] [cyan]None[/cyan]"
        ),
    )
    
    group.add_argument(
        "--geneset_file",
        metavar="",
        type=validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True),
        help=(
            "Pathway file (.gmt). MSigDB or custom GMT supported."
            f"\n[bold green]Default:[/bold green] [cyan]None[/cyan]"
        ),
    )

    group.add_argument(
        "--window_upstream",
        type=int,
        default=35,
        metavar="",
        help=(
            "Upstream gene window in kb. "
            f"[bold green]Default:[/bold green] [cyan]35[/cyan]"
        ),
    )

    group.add_argument(
        "--window_downstream",
        type=int,
        default=10,
        metavar="",
        help=(
            "Downstream gene window in kb. "
            f"[bold green]Default:[/bold green] [cyan]10[/cyan]"
        ),
    )

    group.add_argument(
        "--gene_model",
        default="snp-wise=mean",
        choices=["snp-wise=mean", "snp-wise=top", "snp-wise=all"],
        metavar="",
        help=(
            "MAGMA gene model. "
            f"[bold green]Default:[/bold green] [cyan]snp-wise=mean[/cyan]"
        ),
    )

    group.add_argument(
        "--n_sample_col",
        default="N_COL",
        metavar="",
        help=(
            "Column name for sample sizes in *.pval input. "
            f"[bold green]Default:[/bold green] [cyan]N_COL[/cyan]"
        ),
    )

    return parser


def get_magma_binary_parser(add_help=False):
    """
    Parser containing only MAGMA association-related arguments.
    This is parent-safe and can be inserted into direct/pipeline mode.
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("MAGMA binary")

    # Required inputs (but NOT marked required here → handled by subparser)
    group.add_argument(
        "--magma",
        default="magma",
        metavar="",
        help=(
            "Path to the MAGMA binary. "
            "If not provided, assumes 'magma' is available in your system PATH. "
            f"[bold green]Default:[/bold green] [cyan]magma[/cyan]"
        ),
    )
    return parser



def get_plink_binary_parser(add_help=False):
    """
    Parser containing only MAGMA association-related arguments.
    This is parent-safe and can be inserted into direct/pipeline mode.
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("PLINK binary")
    # Required inputs (but NOT marked required here → handled by subparser)
    group.add_argument(
        "--plink",
        default="plink",
        metavar="",
        help=(
        "Full path to the PLINK executable. "
        "PLINK is required for LD computation used in SuSiE fine-mapping. "
        f"[bold green]Default:[/bold green] [cyan]None[/cyan]"),
    )
    return parser




def get_bcftools_binary_parser(add_help=False):
    """
    Parser containing only MAGMA association-related arguments.
    This is parent-safe and can be inserted into direct/pipeline mode.
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    group = parser.add_argument_group("BCFTOOLS binary")

    # Required inputs (but NOT marked required here → handled by subparser)
    group.add_argument(
        "--bcftools",
        default="bcftools",
        metavar="",
        help=(
            "Path to the bcftools binary. "
            "If not provided, assumes 'bcftools' is available in your system PATH. "
            f"[bold green]Default:[/bold green] [cyan]bcftools[/cyan]"
        ),
    )
    return parser





def get_common_pops_parser(add_help=False):
    """
    Single combined PoPS parser with all argument groups.
    No direct mode, no pipeline mode — all PoPS arguments in one parser.
    Parent-safe for PostGWAS/FLAMES.
    """
    parser = argparse.ArgumentParser(add_help=add_help)

    # ------------------------------------------------------------
    # 1. PoPS Required Inputs
    # ------------------------------------------------------------
    grp_feats = parser.add_argument_group("PoPS Inputs Arguments")

    grp_feats.add_argument(
        "--feature_mat_prefix",
        metavar="",
        help=(
            "[bold bright_red]Required[/bold bright_red]: Prefix to split feature matrix files "
            "(.mat.*.npy, .cols.*.txt, .rows.txt). "
            f"[bold green]Default:[/bold green] [cyan]None[/cyan]"
        ),
    )

    grp_feats.add_argument(
        "--num_feature_chunks",
        type=int,
        metavar="",
        default=2,
        help=(
            "Number of feature matrix chunks. "
            f"[bold green]Default:[/bold green] [cyan]2[/cyan]"
        ),
    )

    grp_feats.add_argument(
        "--pops_gene_loc_file",
        metavar="",
        help=(
            "[bold bright_red]Required[/bold bright_red]: PoPS gene annotation TSV "
            "(must contain ENSGID, CHR, TSS). "
            f"[bold green]Default:[/bold green] [cyan]None[/cyan]"
        ),
    )

    # ------------------------------------------------------------
    # 2. PoPS Covariates
    # ------------------------------------------------------------
    grp_cov = parser.add_argument_group("PoPS Covariates Arguments")

    grp_cov.add_argument(
        "--use_magma_covariates",
        help=(
            "Project out MAGMA covariates before fitting. "
            f"[bold green]Default:[/bold green] [cyan]True[/cyan]"
        ),
        action="store_true",
    )
    grp_cov.add_argument(
        "--ignore_magma_covariates",
        help=(
            "Ignore MAGMA covariates. "
            f"[bold green]Default:[/bold green] [cyan]False[/cyan]"
        ),
        action="store_false",
        dest="use_magma_covariates",
    )
    parser.set_defaults(use_magma_covariates=True)

    grp_cov.add_argument(
        "--use_magma_error_cov",
        help=(
            "Use MAGMA error covariance. "
            f"[bold green]Default:[/bold green] [cyan]True[/cyan]"
        ),
        action="store_true",
    )
    grp_cov.add_argument(
        "--ignore_magma_error_cov",
        help=(
            "Ignore MAGMA error covariance. "
            f"[bold green]Default:[/bold green] [cyan]False[/cyan]"
        ),
        action="store_false",
        dest="use_magma_error_cov",
    )
    parser.set_defaults(use_magma_error_cov=True)

    grp_cov.add_argument(
        "--y_path",
        metavar="",
        help=(
            "Custom target score (must contain ENSGID and Score). "
            f"[bold green]Default:[/bold green] [cyan]None[/cyan]"
        ),
    )

    grp_cov.add_argument(
        "--y_covariates_path",
        metavar="",
        help=(
            "Optional covariates for --y_path. "
            f"[bold green]Default:[/bold green] [cyan]None[/cyan]"
        ),
    )

    grp_cov.add_argument(
        "--y_error_cov_path",
        metavar="",
        help=(
            "Optional error covariance for --y_path (npz or npy). "
            f"[bold green]Default:[/bold green] [cyan]None[/cyan]"
        ),
    )

    grp_cov.add_argument(
        "--project_out_covariates_chromosomes",
        nargs="*",
        metavar="",
        help=(
            "Chromosomes included when projecting out covariates. "
            f"[bold green]Default:[/bold green] [cyan][][/cyan]"
        ),
    )

    grp_cov.add_argument(
        "--project_out_covariates_remove_hla",
        action="store_true",
        help=(
            "Remove HLA genes during covariate projection. "
            f"[bold green]Default:[/bold green] [cyan]True[/cyan]"
        ),
    )
    grp_cov.add_argument(
        "--project_out_covariates_keep_hla",
        action="store_false",
        dest="project_out_covariates_remove_hla",
        help=(
            "Keep HLA genes during covariate projection. "
            f"[bold green]Default:[/bold green] [cyan]False[/cyan]"
        ),
    )
    parser.set_defaults(project_out_covariates_remove_hla=True)

    # ------------------------------------------------------------
    # 3. PoPS Feature Selection
    # ------------------------------------------------------------
    grp_sel = parser.add_argument_group("PoPS Feature Selection")

    grp_sel.add_argument(
        "--subset_features_path",
        metavar="",
        help="Optional feature subset list. "
             f"[bold green]Default:[/bold green] [cyan]None[/cyan]"
    )

    grp_sel.add_argument(
        "--control_features_path",
        metavar="",
        help="Features always included. "
             f"[bold green]Default:[/bold green] [cyan]None[/cyan]"
    )

    grp_sel.add_argument(
        "--feature_selection_chromosomes",
        nargs="*",
        metavar="",
        help="Chromosomes used in feature selection. "
             f"[bold green]Default:[/bold green] [cyan][][/cyan]"
    )

    grp_sel.add_argument(
        "--feature_selection_p_cutoff",
        type=float,
        default=0.05,
        metavar="",
        help=(
            "P-value cutoff for feature selection. "
            f"[bold green]Default:[/bold green] [cyan]0.05[/cyan]"
        ),
    )

    grp_sel.add_argument(
        "--feature_selection_max_num",
        type=int,
        metavar="",
        help="Maximum number of selected features. "
             f"[bold green]Default:[/bold green] [cyan]None[/cyan]"
    )

    grp_sel.add_argument(
        "--feature_selection_fss_num_features",
        type=int,
        metavar="",
        help=("Number of FSS-selected features. "
             f"[bold green]Default:[/bold green] [cyan]None[/cyan]"),
        )

    grp_sel.add_argument(
        "--feature_selection_remove_hla",
        action="store_true",
        help=(
            "Remove HLA genes during feature selection. "
            f"[bold green]Default:[/bold green] [cyan]True[/cyan]"
        ),
    )
    grp_sel.add_argument(
        "--feature_selection_keep_hla",
        action="store_false",
        dest="feature_selection_remove_hla",
        help=(
            "Keep HLA genes during feature selection. "
            f"[bold green]Default:[/bold green] [cyan]False[/cyan]"
        ),
    )
    parser.set_defaults(feature_selection_remove_hla=True)

    # ------------------------------------------------------------
    # 4. PoPS Training
    # ------------------------------------------------------------
    grp_train = parser.add_argument_group("PoPS Training Arguments")

    grp_train.add_argument(
        "--training_chromosomes",
        nargs="*",
        metavar="",
        help=(
            "Chromosomes used when computing coefficients. "
            f"[bold green]Default:[/bold green] [cyan][][/cyan]"
        ),
    )

    grp_train.add_argument(
        "--training_remove_hla",
        action="store_true",
        help=(
            "Remove HLA genes during training. "
            f"[bold green]Default:[/bold green] [cyan]True[/cyan]"
        ),
    )
    grp_train.add_argument(
        "--training_keep_hla",
        action="store_false",
        dest="training_remove_hla",
        help=(
            "Keep HLA genes during training. "
            f"[bold green]Default:[/bold green] [cyan]False[/cyan]"
        ),
    )
    parser.set_defaults(training_remove_hla=True)

    # ------------------------------------------------------------
    # 5. Runtime Settings
    # ------------------------------------------------------------
    grp_misc = parser.add_argument_group("PoPS Runtime Settings")

    grp_misc.add_argument(
        "--method",
        default="ridge",
        metavar="",
        help=(
            "Regularization method: ridge, lasso, or linreg. "
            f"[bold green]Default:[/bold green] [cyan]ridge[/cyan]"
        ),
    )

    grp_misc.add_argument(
        "--save_matrix_files",
        action="store_true",
        help=(
            "Save matrices used during PoPS. "
            f"[bold green]Default:[/bold green] [cyan]False[/cyan]"
        ),
    )
    grp_misc.add_argument(
        "--no_save_matrix_files",
        action="store_false",
        dest="save_matrix_files",
        help=(
            "Do not save matrices. "
            f"[bold green]Default:[/bold green] [cyan]False[/cyan]"
        ),
    )
    parser.set_defaults(save_matrix_files=False)

    return parser


from rich_argparse import RawTextRichHelpFormatter

def get_common_sumstat_filter_parser(add_help=False):
    parser = argparse.ArgumentParser(formatter_class=RawTextRichHelpFormatter,add_help=add_help)
    
    group = parser.add_argument_group("SUMSTAT QC Filtering Arguments")

    # =====================================================
    # NUMERIC FILTERS
    # =====================================================
    group.add_argument(
        "--pval-cutoff",
        type=float,
        metavar=" ",
        default=None,
        help=(
            "Filter variants (Retain) by [bold]−log10(P)[/bold] ≥ the provided cutoff.\n"
            "Example: To retain variants with P <=5e-8, use [cyan]7.3[/cyan]\n"
            "[bold green]Default:[/bold green] [cyan]None[/cyan] (no filtering applied). "
            "Make sure the LP field is present in the VCF FORMAT column."
        )
    )
    
    group.add_argument(
        "--maf-cutoff",
        type=float,
        metavar=" ",
        default=0.002,
        help=(
            " Minimum Minor Allele Frequency (MAF) required to retain a variant .\n"
            "[bold green]Default:[/bold green] [cyan]0.001[/cyan] (no MAF filtering). "
            "Make sure the AF field is present in the VCF FORMAT column."
        )
    )

    group.add_argument(
        "--external-af-name",
        default="EUR",
        metavar=" ",
        help=(
            "Name of external allele-frequency tag inside the VCF header "
            "(e.g., EUR, AFR, EAS).\n"
            "[bold green]Default:[/bold green] [cyan]EUR[/cyan] . "
            "Make sure the mentioned external-af-name should be present in the VCF INFO column."
        )
    )
    
    group.add_argument(
        "--allelefreq-diff-cutoff",
        type=float,
        metavar=" ",
        default=0.2,
        help=(
            "Maximum allowed difference between dataset allele frequency and "
            "external reference AF (e.g., EUR/EAS/AFR AF tags).\n"
            "[bold green]Default:[/bold green] [cyan]0.2[/cyan]"
        )
    ) 
    
    group.add_argument(
        "--info-cutoff",
        type=float,
        metavar=" ",
        default=None,
        help=(
            "Primary minimum INFO (imputation quality) threshold.\n"
            "Variants with INFO < cutoff will be removed.\n\n"
            "Priority:\n"
            "  • Overrides --info-min if both are provided\n"
            "  • Can be combined with --info-max\n\n"
            "Typical values: 0.6, 0.8\n"
            "Default: None"
        )
    )

    group.add_argument(
        "--info-min",
        type=float,
        metavar=" ",
        default=0.3,
        help=(
            "Minimum INFO threshold (used only if --info-cutoff is NOT provided).\n"
            "Variants with INFO < info-min will be removed.\n\n"
            "Ignored when --info-cutoff is set.\n"
            "Default: 0.3"
        )
    )

    group.add_argument(
        "--info-max",
        type=float,
        metavar=" ",
        default=None,
        help=(
            "Maximum INFO threshold.\n"
            "Variants with INFO > info-max will be removed.\n\n"
            "Can be used together with --info-cutoff or --info-min.\n"
            "Default: None"
        )
    )

    group.add_argument(
        "--info-missing",
        type=str,
        metavar=" ",
        default="remove",
        choices=["keep", "remove"],
        help=(
            "How to handle variants with missing INFO values.\n\n"
            "Options:\n"
            "  • remove → exclude variants with missing INFO (recommended)\n"
            "  • keep   → retain variants even if INFO is missing\n\n"
            "Default: remove"
        )
    )


    # =====================================================
    # FLAGS
    # =====================================================
    group.add_argument(
        "--include-indels",
        action="store_true",
        help=(
            "Include INDEL variants in the QC output.\n\n"
            "If this flag is provided, INDELs will be retained.\n"
            "If omitted, INDELs will be excluded (default behaviour).\n\n"
            "[bold green]Default:[/bold green] [cyan]False (INDELs excluded)[/cyan]"
        )
    )
    
    group.add_argument(
        "--exclude-palindromic",
        action="store_true",
        help=textwrap.dedent(
            """
            Control whether ambiguous palindromic SNPs (A/T, C/G) are removed.

            By default, all palindromic SNPs are kept. When this flag is used,
            palindromic SNPs whose allele frequency falls within the
            ambiguity interval defined by [cyan]--palindromic-af-lower[/cyan]
            and [cyan]--palindromic-af-upper[/cyan] will be removed.

            Examples:
            Keep all palindromic SNPs (default):
                (omit this flag)

            Remove ambiguous palindromic SNPs:
                --exclude-palindromic

            [bold green]Default:[/bold green] [cyan]Keep all palindromic SNPs[/cyan]
            """
        )
    )



    # =====================================================
    # PALINDROMIC AF BOUNDS
    # =====================================================
    group.add_argument(
        "--palindromic-af-lower",
        type=float,
        default=0.4,
        metavar=" ",
        help=(
            "Lower AF bound used to detect ambiguous palindromic SNPs when the flag \n"
            "[cyan]--exclude-palindromic[/cyan] is used.\n"
            "[bold green]Default:[/bold green] [cyan]0.4[/cyan]"
        )
    )

    group.add_argument(
        "--palindromic-af-upper",
        type=float,
        default=0.6,
        metavar=" ",
        help=(
            "Upper AF bound used to detect ambiguous palindromic SNPs when the flag \n"
            "[cyan]--exclude-palindromic[/cyan].\n"
            "[bold green]Default:[/bold green] [cyan]0.6[/cyan]"
        )
    )

    # =====================================================
    # MHC REGION
    # =====================================================
    group.add_argument(
        "--remove-mhc",
        action="store_true",
        default=False,
        help=(
            "If you want to Remove MHC Region use --remove-mhc \n"
            "[bold green]Default:[/bold green] [cyan]MHC Region will not be removed[/cyan]"
        )
    )

    group.add_argument(
        "--mhc-chrom",
        default="6",
        metavar=" ",
        help=(
            "Chromosome of the MHC region.\n"
            "[bold green]Default:[/bold green] [cyan]6[/cyan]"
        )
    )
        
    group.add_argument(
        "--mhc-start",
        type=int,
        default=25_000_000,
        metavar=" ",
        help=(
            "Start coordinate of the MHC region.\n"
            "[bold green]Default:[/bold green] [cyan]25,000,000[/cyan]"
        )
    )

    group.add_argument(
        "--mhc-end",
        type=int,
        default=34_000_000,
        metavar=" ",
        help=(
            "End coordinate of the MHC region.\n"
            "[bold green]Default:[/bold green] [cyan]34,000,000[/cyan]"
        )
    )

    return parser




# ============================================================
# COMMON FLAMES PARSER (shared by DIRECT + PIPELINE)
# ============================================================
def get_flames_common_parser(add_help=False):
    """
    Arguments shared by BOTH:
      • FLAMES Direct
      • FLAMES Pipeline
    """
    parser = argparse.ArgumentParser(add_help=add_help)

    # ────────────────────────────────
    # FLAMES Model Options
    # ────────────────────────────────
    model_group = parser.add_argument_group("FLAMES Model Options")

    model_group.add_argument(
        "--flames_annot_dir",
        metavar="",
        type=validate_path(must_exist=True, must_be_dir=True),
        help="[bold bright_red]Required[/bold bright_red]: Directory where FLAMES annotation files will be written.",
    )
    
    model_group.add_argument(
        "--distance",
        type=int,
        default=750000,
        help="Gene distance (bp) for locus-to-gene mapping.",
        metavar=""
    )

    model_group.add_argument(
        "--weight",
        metavar="",
        type=float,
        default=0.725,
        help="Weight for XGB model (default 0.725). Remaining weight applied to PoPS."
    )

    model_group.add_argument(
        "--modelpath",
        metavar="",
        help="Directory containing pretrained FLAMES model.",
        required=False
    )

    model_group.add_argument(
        "--true_positives",
        metavar="",
        type=validate_path(must_exist=True, must_be_dir=True),
        help="File describing true-positive genes per locus.",
        default=False
    )

    model_group.add_argument("--cmd_vep", metavar="", default=False,
                             help="Path to the VEP executable.")
    model_group.add_argument("--vep_cache", metavar="", default=False,
                             help="Path to the VEP cache.")
    model_group.add_argument("--tabix", metavar="", default=False,
                             help="Path to tabix executable.")
    model_group.add_argument(
        "--CADD_file",
        metavar="",
        default=False,
        help="Local CADD file (must match genome build)."
    )

    return parser



def get_common_imputation_parser(add_help: bool = False) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=add_help)
    grp = parser.add_argument_group("Imputation Inputs Arguments")

    grp.add_argument(
        "--imputation_tool",
        metavar=" ",
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
        "--ref_ld",
        metavar=" ",
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
        "--r2threshold",
        metavar=" ",
        type=float,
        default=0.8,
        help=("LD r² tagging threshold."
              "[bold green]Default:[/bold green] [cyan]0.8[/cyan]")
    )

    grp.add_argument(
        "--maf",
        metavar=" ",
        type=float,
        default=0.001,
        help=("Minor allele frequency threshold."
              "[bold green]Default:[/bold green] [cyan]0.001[/cyan]")
    )

    grp.add_argument(
        "--population",
        metavar=" ",
        type=str,
        default="EUR",
        help=("Population code for LD reference."
             "[bold green]Default:[/bold green] [cyan]EUR[/cyan]")
    )

    grp.add_argument(
        "--ref",
        metavar=" ",
        type=str,
        default="TOP_LD",
        help=("Reference panel name label."
              "[bold green]Default:[/bold green] [cyan]TOP_LD[/cyan]") 
    )

    grp.add_argument(
        "--corr_method",
        metavar=" ",
        type=str,
        default="pearson",
        choices=["pearson", "spearman"],
        help=("Correlation method used in QC/post-processing."
              "[bold green]Default:[/bold green] [cyan]pearson[/cyan]") 
              
    )

    grp.add_argument(
        "--gwas2vcf_resource",
        metavar=" ",
        default=None,
        type=validate_path(must_exist=True, must_be_dir=True),
        help="[bold bright_red]Required[/bold bright_red]: GWAS2VCF resource folder.",
    )
    
    grp.add_argument(
        "--gwas2vcf_default_config",
        metavar=" ",
        type=validate_path(must_exist=True, must_be_file=True),
        help="[bold bright_red]Required[/bold bright_red]: GWAS2VCF default config file"
              
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
    grp = parser.add_argument_group("LDSC Inputs")

    grp.add_argument(
        "--heritability_tool",
        metavar=" ",
        type=str,
        choices=["ldsc"],
        default="ldsc",
        help="Choose imputation engine (default: ldsc). [bold green]Default:[/bold green] [cyan]ldsc[/cyan]",
    )
    
    grp.add_argument(
        "--merge-alleles",
        metavar=" ",
        type=validate_path(must_exist=True, must_be_file=True),
        help=(
            "[bold bright_red]Required[/bold bright_red]: HapMap3 SNP list with alleles (w_hm3.snplist). Used during munge_sumstats.\n"
            "Optional in DIRECT mode; recommended in PIPELINE mode."
        ),
    )

    ref = parser.add_argument_group("LDSC Reference ")

    ref.add_argument(
        "--ref-ld-chr",
        metavar=" ",
        type=validate_path(must_exist=True, must_be_dir=True),
        help="[bold bright_red]Required[/bold bright_red]: Directory containing reference LD scores (e.g., eur_w_ld_chr/).",
    )

    ref.add_argument(
        "--w-ld-chr",
        metavar=" ",
        type=validate_path(must_exist=True, must_be_dir=True),
        help="[bold bright_red]Required[/bold bright_red]: Directory containing LDSC regression weights (often same as --ref-ld-chr).",
    )

    prev = parser.add_argument_group("LDSC Liability-Scale Parameters")

    prev.add_argument(
        "--samp-prev",
        metavar=" ",
        type=float,
        default=None,
        help=("Sample prevalence (cases / total) used for liability h²."
              "[bold green]Default:[/bold green] [cyan]None[/cyan]" )
        
    )

    prev.add_argument(
        "--pop-prev",
        metavar=" ",
        type=float,
        default=None,
        help=("Population prevalence (e.g., 0.01)."
              "[bold green]Default:[/bold green] [cyan]None[/cyan]")
    )

    mung = parser.add_argument_group("LDSC Munge Settings")

    mung.add_argument(
        "--heritability_info-min",
        metavar=" ",
        type=float,
        default=0.9,
        help=("Minimum INFO score during munge_sumstats (default: 0.9)."
              "[bold green]Default:[/bold green] [cyan]0.7[/cyan]" )
    )

    mung.add_argument(
        "--heritability_maf-min",
        metavar=" ",
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
    output_group = parser.add_argument_group()

    output_group.add_argument(
        "--pdf",
        metavar=" ",
        help="Output PDF file (e.g., result.pdf)."
    )

    output_group.add_argument(
        "--png",
        metavar=" ",
        help="Output PNG file (e.g., result.png)."
    )
    # ------------------------------------------------------------
    # INPUTS
    # ------------------------------------------------------------
    input_group = parser.add_argument_group("Association Plot Input Arguments")

    input_group.add_argument(
        "--pheno",
        metavar=" ",
        help="Phenotype name to extract from GWAS-VCF file."
    )
    # ------------------------------------------------------------
    # FLAGS
    # ------------------------------------------------------------
    flag_group = parser.add_argument_group("Flags")

    flag_group.add_argument(
        "--as",
        dest="allelic_shift",
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
    numeric_group = parser.add_argument_group("Numeric Plot Options")
    numeric_group.add_argument(
        "--nauto",
        type=int,
        default=22,
        metavar=" ",
        help="Number of autosomes. [bold green]Default:[/bold green] [cyan]22[/cyan]"
    )

    numeric_group.add_argument(
        "--min-af",
        type=float,
        default=0.0,
        metavar=" ",
        help="Minimum allele frequency filter. [bold green]Default:[/bold green] [cyan]0.0[/cyan]"
    )

    numeric_group.add_argument(
        "--min-lp",
        type=int,
        default=2,
        metavar=" ",
        help="Minimum −log10(P) threshold. [bold green]Default:[/bold green] [cyan]2[/cyan]"
    )

    numeric_group.add_argument(
        "--loglog-pval",
        type=int,
        default=10,
        metavar=" ",
        help="-log10(P) threshold for switching to log-log scale. "
             "[bold green]Default:[/bold green] [cyan]10[/cyan]"
    )

    numeric_group.add_argument(
        "--cyto-ratio",
        type=int,
        default=25,
        metavar=" ",
        help="Plot height : cytoband height ratio. "
             "[bold green]Default:[/bold green] [cyan]25[/cyan]"
    )

    numeric_group.add_argument(
        "--max-height",
        type=int,
        metavar=" ",
        help="Maximum vertical height of plot (optional)."
    )

    numeric_group.add_argument(
        "--spacing",
        type=int,
        default=20,
        metavar=" ",
        help="Spacing between chromosomes. "
             "[bold green]Default:[/bold green] [cyan]20[/cyan]"
    )
    # ------------------------------------------------------------
    # FIGURE SIZE
    # ------------------------------------------------------------
    fig_group = parser.add_argument_group("Figure Dimensions")

    fig_group.add_argument(
        "--width",
        type=float,
        default=7.0,
        metavar=" ",
        help="Plot width in inches. [bold green]Default:[/bold green] [cyan]7.0[/cyan]"
    )

    fig_group.add_argument(
        "--height",
        type=float,
        metavar=" ",
        help="Plot height in inches (optional)."
    )

    fig_group.add_argument(
        "--fontsize",
        type=int,
        default=12,
        metavar=" ",
        help="Font size. [bold green]Default:[/bold green] [cyan]12[/cyan]"
    )

    return parser


