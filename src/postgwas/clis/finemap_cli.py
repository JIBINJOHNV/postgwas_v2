#!/usr/bin/env python3
import argparse
from textwrap import dedent
from rich_argparse import RichHelpFormatter

from postgwas.utils.main import validate_path,detect_total_memory_gb,validate_alphanumeric,validate_prefix_files

class HelpOnErrorArgumentParser(argparse.ArgumentParser):
    """Suppress default argparse error spam and show rich help instead."""
    def error(self, message):
        raise ValueError(message)


def get_finemap_common_parser(add_help=False):
    """
    Create an argument parser for selecting the fine-mapping method.
    This helper is designed to be used as a parent parser in other CLIs.
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    grp = parser.add_argument_group("Fine-map common arguments")

    # ------------------------------------------------------------------
    # Method selection
    # ------------------------------------------------------------------
    grp.add_argument(
        "--finemap_method",
        metavar=" ",
        nargs="+",
        choices=["susie", "finemap"],
        default="susie",
        help=(
            "Fine-mapping method to use. "
            "[bold]Choices:[/bold] [bright_yellow]{susie, finemap}[/bright_yellow]. "
            "[bold green]Default:[/bold green] [cyan]susie[/cyan]"
        ),
    )

    grp.add_argument(
            "--locus_file",
            metavar="FILE",
            type=validate_path(must_exist=True, must_be_file=True),
            help=(
                "[bold bright_red]Required[/bold bright_red]: Path to the file defining genomic regions. "
                "\n- If [bold]locus_type[/bold] is 'range': Must contain [cyan]CHROM, START, END[/cyan] columns. "
                "\n- If [bold]locus_type[/bold] is 'point': Must contain [cyan]CHROM, POS[/cyan] columns."
            ),
        )
    grp.add_argument(
            "--locus_type",
            choices=["range", "point"],
            default="range",
            help=(
                "Determines how boundaries are constructed (Default: [bold]range[/bold]):"
                "\n[bold]range[/bold]: Uses START/END columns. [italic]--window_kb[/italic] acts as an optional flank. "
                "\n[bold]point[/bold]: Uses the POS column. [italic]--window_kb[/italic] is used to create the window."
            ),
        )
    grp.add_argument(
            "--window_kb",
            metavar="INT",
            type=int,
            default=500,
            help=(
                "The distance (in Kilobases) to extend the locus boundaries."
                "\n- [bold]point[/bold] mode: Creates a symmetric window ([italic]POS ± window_kb[/italic]). "
                "\n- [bold]range[/bold] mode: Adds a flank to the existing coordinates ([italic]START - window_kb[/italic] and [italic]END + window_kb[/italic])."
                "\nSet to 0 if you want to use the exact coordinates in 'range' mode. Default: 500."
            ),
        )

    # ------------------------------------------------------------------
    # Locus selection / tuning
    # ------------------------------------------------------------------
    grp.add_argument(
        "--lp_threshold",
        metavar=" ",
        type=float,
        default=7.3,
        help=(
            "Minus log10(P) threshold used to include a locus for fine-mapping. "
            "The locus file must contain an LP column representing −log10(P) values. "
            "[bold green]Default:[/bold green] [cyan]7.3[/cyan]"
        ),
    )

    grp.add_argument(
        "--min_ram_per_worker_gb",
        type=int,
        metavar=" ",
        default=14,
        help=(
            "Minimum RAM (in GB) reserved per worker when running fine-mapping "
            "in parallel. Used for auto-detecting optimal worker count. "
            "[bold green]Default:[/bold green] [cyan]14GB[/cyan]"
        ),
    )

    # ------------------------------------------------------------------
    # MHC handling (mutually exclusive, clean default)
    # ------------------------------------------------------------------
    grp.add_argument(
        "--finemap_skip_mhc",
        action="store_true",
        help=(
            "Skip the extended MHC region (chr6:25–35Mb). "
            "[bold green]Default behavior[/bold green]."
        ),
    )

    grp.add_argument(
        "--finemap_include_mhc",
        action="store_true",
        help=(
            "Include the extended MHC region (chr6:25–35Mb). "
            "Overrides the default skip behavior."
        ), 
    ) 

    grp.add_argument(
        "--finemap_mhc_chrom",
        type=int,
        default=6,
        metavar=" ",
        help=(
            "Chromosome containing the MHC region to skip during fine-mapping "
            "(used unless --finemap_include_mhc is set). "
            "[bold green]Default:[/bold green] [cyan]6[/cyan]"
        ),
    )

    grp.add_argument(
        "--finemap_mhc_start",
        type=int,
        default=25000000,
        metavar=" ",
        help=(
            "Start coordinate for the MHC region to skip "
            "(used unless --finemap_include_mhc is set). "
            "[bold green]Default:[/bold green] [cyan]25000000[/cyan]"
        ),
    )

    grp.add_argument(
        "--finemap_mhc_end",
        type=int,
        default=35000000,
        metavar=" ",
        help=(
            "End coordinate for the MHC region to skip "
            "(used unless --finemap_include_mhc is set). "
            "[bold green]Default:[/bold green] [cyan]35000000[/cyan]"
        ),
    )

    # ------------------------------------------------------------------
    # LD reference
    # ------------------------------------------------------------------
    grp.add_argument(
        "--finemap_ld_ref",
        metavar=" ",
        help=(
            "[bold bright_red]Required[/bold bright_red]: Prefix of the PLINK LD "
            "reference panel (e.g., 1000G EUR). Should correspond to files: "
            "PREFIX.bed, PREFIX.bim, PREFIX.fam."
        ),
    )

    # Default behavior: skip MHC
    parser.set_defaults(finemap_skip_mhc=True)

    return parser


# =========================================================
# SHARED DIRECT-MODE INPUT ARGUMENTS
# =========================================================
def get_finemap_susie_inputs_parser(add_help=False):
    """
    SuSiE/FINEMAP input files that exist ONLY in direct mode.
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    grp = parser.add_argument_group("SuSiE Specific Inputs")
    grp.add_argument(
        "--susie_input_file",
        metavar=" ",
        type=validate_path(must_exist=True, must_be_file=True),
        help=(
            "[bold bright_red]Required[/bold bright_red]: "
            "Fine-mapping–ready summary statistics file. "
            "Should be generated by the postgwas formatter module."
        )
    )

    return parser



def get_common_susie_arguments(add_help=False):
    """
    Add SuSiE tuning arguments. 
    Arguments are added directly to the 'susie' group object to force 
    them to appear under the SuSiE heading in the help menu.
    """
    parser = argparse.ArgumentParser(
        add_help=add_help,
        formatter_class=RichHelpFormatter,
    )

    # All arguments added to this group WILL stay under this heading
    susie = parser.add_argument_group("SuSiE Specific Arguments")

    susie.add_argument(
        "--L",
        type=int,
        default=10,
        metavar="",
        help=(
            "Maximum number of SuSiE credible sets per locus. "
            "Increasing this may increase runtime and memory usage. "
            "[bold green]Default:[/bold green] [cyan]10[/cyan]"
        ),
    )

    # ------------------------------------------------------------------
    # Resource / timeout controls
    # ------------------------------------------------------------------
    susie.add_argument(
        "--timeout_ld_seconds",
        type=int,
        default=180,
        metavar=" ",
        help=(
            "Maximum time (in seconds) allowed for PLINK LD-matrix computation "
            "per locus. Execution aborts for loci exceeding this limit. "
            "[bold green]Default:[/bold green] [cyan]180[/cyan]"
        ),
    )

    susie.add_argument(
        "--timeout_susie_seconds",
        type=int,
        default=180,
        metavar=" ",
        help=(
            "Maximum time (in seconds) allowed for SuSiE model fitting per locus. "
            "Loci exceeding the limit are skipped with a warning. "
            "[bold green]Default:[/bold green] [cyan]180[/cyan]"
        ),
    )

    return parser


def get_finemap_finemap_arguments(add_help=False):
    """
    FINEMAP input files that exist ONLY in direct mode.
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    grp = parser.add_argument_group("FINEMAP Specific Inputs")

    # ===============================================================
    # Input / execution
    # ===============================================================
    grp.add_argument(
        "--finemap-in-files",
        metavar=" ",
        type=validate_path(must_exist=True, must_be_file=True),
        help=(
            "FINEMAP input file created by postgwas formtter module with finemap_finemap option "
        )
    )

    return parser


def get_common_finemap_finemap_arguments(add_help=False):
    """
    Add FINEMAP tuning arguments with 'sss' (Shotgun Stochastic Search) set as the default mode.
    Includes stylized default values and mode selection clarification.
    """
    try:
        formatter = RichHelpFormatter
    except NameError:
        formatter = argparse.ArgumentDefaultsHelpFormatter

    parser = argparse.ArgumentParser(
        add_help=add_help,
        formatter_class=formatter,
    )

    finemap = parser.add_argument_group("FINEMAP Specific Arguments")

    # ===============================================================
    # Subprograms (Modes)
    # ===============================================================
    finemap.add_argument(
        "--sss",
        action="store_true",
        help="Fine-mapping with [bold]shotgun stochastic search[/bold]. "
            "[bold yellow]Use either [cyan]--sss[/cyan] or [cyan]--cond[/cyan] (not both). If both used then --sss will be consider[/bold yellow] "
            "[bold green]Default:[/bold green] [cyan]--sss[/cyan]"
    )

    finemap.add_argument(
        "--cond",
        action="store_true",
        help=( "Fine-mapping with [bold]stepwise conditional search[/bold]. "
            "[bold yellow]Use either [cyan]--sss[/cyan] or [cyan]--cond[/cyan] (not both).[/bold yellow]" )
    )

    # ===============================================================
    # Shotgun stochastic search (SSS)
    # ===============================================================
    finemap.add_argument(
        "--n-iter",
        metavar=" ",
        type=int,
        default=100000,
        help=("Maximum number of SSS iterations. "
              "[bold green]Default:[/bold green] [cyan]100000[/cyan]" )
    )

    finemap.add_argument(
        "--n-conv-sss",
        type=int,
        metavar=" ",
        default=1000,
        help=("Iterations used to check SSS convergence. "
              "[bold green]Default:[/bold green] [cyan]1000[/cyan]")
    )

    finemap.add_argument(
        "--prob-conv-sss-tol",
        metavar=" ",
        type=float,
        default=0.001,
        help=("Option to set the tolerance at which the added probability mass (over --n-conv-sss iterations) is considered small enough to terminate the shotgun stochastic search. "
              "[bold green]Default:[/bold green] [cyan]0.001[/cyan]" )
    )

    finemap.add_argument(
        "--n-configs-top",
        metavar=" ",
        type=int,
        default=50000,
        help=("Number of top causal configurations to save. "
              "[bold green]Default:[/bold green] [cyan]50000[/cyan]" )
    )

    # ===============================================================
    # Causal priors
    # ===============================================================
    finemap.add_argument(
        "--n-causal-snps",
        type=int,
        metavar=" ", 
        default=5,
        help=("Maximum number of causal SNPs. "
              "[bold green]Default:[/bold green] [cyan]5[/cyan]"  )
    )

    finemap.add_argument(
        "--prior-k",
        action="store_true",
        help=("Use prior probabilities for number of causal SNPs from K files. "
              )
    )

    # finemap.add_argument(
    #     "--prior-k0",
    #     type=float,
    #     metavar=" ",
    #     help=("Prior probability that there are zero causal SNPs. ")
    # )

    # finemap.add_argument(
    #     "--prior-snps",
    #     action="store_true",
    #     help="Use per-SNP causal prior probabilities from Z file"
    # )

    # finemap.add_argument(
    #     "--prior-std",
    #     metavar=" ",
    #     default="0.05",
    #     help=( "Comma-separated prior SDs of effect sizes. "
    #           "[bold green]Default:[/bold green] [cyan]0.05[/cyan]" )
    # )

    # ===============================================================
    # LD / correlation handling
    # ===============================================================
    finemap.add_argument(
        "--corr-config",
        metavar=" ",
        type=float,
        default=0.95,
        help=("Option to set the posterior probability of a causal configuration to zero if it includes a pair of SNPs with absolute correlation above this threshold "
              "[bold green]Default:[/bold green] [cyan]0.95[/cyan]" )
    )

    finemap.add_argument(
        "--collinear-tol",
        metavar=" ",
        type=float,
        default=0.95,
        help=("Collinearity tolerance for conditional search. "
              "[bold green]Default:[/bold green] [cyan]0.95[/cyan]" )
    )

    finemap.add_argument(
        "--force-n-samples",
        action="store_true",
        help=("Allow LD computed using different sample size than GWAS "
              ) 
    )

    # ===============================================================
    # Filtering / credibility
    # ===============================================================
    finemap.add_argument(
        "--pvalue-snps",
        metavar=" ",
        type=float,
        default=1.0,
        help=("Option to set a p-value threshold at which SNPs are included "
              "[bold green]Default:[/bold green] [cyan]1.0[/cyan]"  )
    )

    finemap.add_argument(
        "--cond-pvalue",
        metavar=" ",
        type=float,
        default=5e-8,
        help=("Genome-wide significance threshold (conditional search). "
              "[bold green]Default:[/bold green] [cyan]5e-8[/cyan]" )
    )

    finemap.add_argument(
        "--prob-cred-set",
        metavar=" ",
        type=float,
        default=0.95,
        help=( "Option to set the probability at which the credible interval includes a causal SNP. "
              "[bold green]Default:[/bold green] [cyan]0.95[/cyan]"  )
    )

    finemap.add_argument(
        "--std-effects",
        action="store_true",
        help="Option to print mean and standard deviation of the posterior effect size distribution for standardized dosages"
    )

    return parser
