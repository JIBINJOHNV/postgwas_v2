#!/usr/bin/env python3
import argparse
from typing import get_args

from postgwas.config import load_configuration
from postgwas.config.models.modules.fine_mapping import FineMappingEngine
from postgwas.core.ui import (
    AlignedRichHelpFormatter,
    format_cli_default,
    help_with_default,
    style_cli_requirement,
)
from postgwas.core.execution.runtime import validate_path


def get_finemap_common_parser(add_help=False, *, include_genome_build=False):
    """
    Create an argument parser for selecting the fine-mapping method.
    This helper is designed to be used as a parent parser in other CLIs.
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    fine_mapping_defaults = load_configuration().modules.fine_mapping
    configured_engine = fine_mapping_defaults.engine
    grp = parser.add_argument_group("Fine-mapping settings")

    # ------------------------------------------------------------------
    # Method selection
    # ------------------------------------------------------------------
    grp.add_argument(
        "--finemap-method",
        metavar="METHOD",
        choices=get_args(FineMappingEngine),
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Fine-mapping method to use",
            configured_engine,
        ),
    )

    grp.add_argument(
        "--locus-file",
        metavar="FILE",
        type=validate_path(must_exist=True, must_be_file=True),
        default=argparse.SUPPRESS,
        help=style_cli_requirement(
            "Path to a predefined "
            "locus file. "
            "\n- If [bold]locus_type[/bold] is 'range': Must contain "
            "[cyan]CHR, START, END, LP[/cyan] columns. Use "
            "[italic]--window-kb 0[/italic] to preserve START/END as the exact "
            "fine-mapping boundaries. "
            "\n- If [bold]locus_type[/bold] is 'point': Must contain "
            "[cyan]CHR, POS, LP[/cyan] columns.",
            required=True,
        ),
    )
    grp.add_argument(
        "--locus-type",
        choices=["range", "point"],
        default=argparse.SUPPRESS,
        help=(
            "Determines how boundaries are constructed (%s):"
            "\n[bold]range[/bold]: Uses START/END columns. [italic]--window-kb[/italic] acts as an optional flank. "
            "\n[bold]point[/bold]: Uses the POS column. [italic]--window-kb[/italic] is used to create the window."
        ) % format_cli_default(
            fine_mapping_defaults.locus_type,
        ),
    )
    grp.add_argument(
        "--window-kb",
        metavar="INT",
        type=int,
        default=argparse.SUPPRESS,
        help=(
            "The flank, in kilobases on each side, used to construct the "
            "fine-mapping boundary."
            "\n- [bold]point[/bold] mode: Creates a symmetric window ([italic]POS ± window_kb[/italic]). "
            "\n- [bold]range[/bold] mode: Adds a flank to the existing coordinates ([italic]START - window_kb[/italic] and [italic]END + window_kb[/italic])."
            "\nSet to 0 if you want to use the exact coordinates in 'range' mode. %s."
        ) % format_cli_default(
            fine_mapping_defaults.locus_window_kb,
        ),
    )

    # ------------------------------------------------------------------
    # Locus selection / tuning
    # ------------------------------------------------------------------
    grp.add_argument(
        "--lp-threshold",
        metavar="LP",
        type=float,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Minus log10(P) threshold used to include a locus for fine-mapping. "
            "The locus file must contain an LP column representing −log10(P) values",
            fine_mapping_defaults.locus_lp_threshold,
        ),
    )

    grp.add_argument(
        "--minimum-memory-per-worker-gb",
        type=float,
        metavar="GB",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Minimum RAM (in GB) reserved per worker when running fine-mapping "
            "in parallel. Used for auto-detecting optimal worker count",
            fine_mapping_defaults.memory_per_worker_gb,
        ),
    )

    grp.add_argument(
        "--maximum-variants-per-locus",
        type=int,
        metavar="COUNT",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Maximum harmonized variants permitted in one primary or joint "
            "fine-mapping locus before dense LD construction. Increase only "
            "when sufficient memory is available",
            fine_mapping_defaults.ld_resource_guard.maximum_variants_per_locus,
        ),
    )

    # ------------------------------------------------------------------
    # MHC handling (mutually exclusive, clean default)
    # ------------------------------------------------------------------
    mhc_mode = grp.add_mutually_exclusive_group()
    mhc_mode.add_argument(
        "--finemap-skip-mhc",
        action="store_true",
        dest="finemap_skip_mhc",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Skip the configured extended MHC region",
            fine_mapping_defaults.skip_mhc,
        ),
    )

    mhc_mode.add_argument(
        "--finemap-include-mhc",
        action="store_false",
        dest="finemap_skip_mhc",
        default=argparse.SUPPRESS,
        help=(
            "Include loci in the configured extended MHC region for this run."
        ),
    )

    grp.add_argument(
        "--finemap-mhc-chromosome",
        type=int,
        default=argparse.SUPPRESS,
        metavar="CHROM",
        help=help_with_default(
            "Chromosome containing the MHC region",
            fine_mapping_defaults.mhc_chromosome,
        ),
    )

    grp.add_argument(
        "--finemap-mhc-start",
        type=int,
        default=argparse.SUPPRESS,
        metavar="BP",
        help=help_with_default(
            "Start coordinate for the MHC region",
            fine_mapping_defaults.mhc_start,
        ),
    )

    grp.add_argument(
        "--finemap-mhc-end",
        type=int,
        default=argparse.SUPPRESS,
        metavar="BP",
        help=help_with_default(
            "End coordinate for the MHC region",
            fine_mapping_defaults.mhc_end,
        ),
    )

    if include_genome_build:
        grp.add_argument(
            "--genome-build",
            choices=["GRCh37", "GRCh38"],
            default=argparse.SUPPRESS,
            metavar="BUILD",
            help=help_with_default(
                "Genome build used by locus coordinates, summary statistics, and LD",
                fine_mapping_defaults.genome_build.value,
            ),
        )

    # ------------------------------------------------------------------
    # LD reference
    # ------------------------------------------------------------------
    grp.add_argument(
        "--finemap-ld-reference",
        metavar="PREFIX",
        default=argparse.SUPPRESS,
        help=style_cli_requirement(
            "Prefix of the PLINK LD "
            "reference panel (e.g., 1000G EUR). Should correspond to files: "
            "PREFIX.bed, PREFIX.bim, PREFIX.fam.",
            required=True,
        ),
    )

    return parser


# =========================================================
# SHARED DIRECT-MODE INPUT ARGUMENTS
# =========================================================
def get_finemap_susie_inputs_parser(add_help=False):
    """
    SuSiE/FINEMAP input files that exist ONLY in direct mode.
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    grp = parser.add_argument_group("SuSiE input")
    grp.add_argument(
        "--susie-input-file",
        metavar="PATH",
        type=validate_path(must_exist=True, must_be_file=True),
        default=argparse.SUPPRESS,
        help=style_cli_requirement(
            "Fine-mapping–ready summary statistics file. Should be generated "
            "by the postgwas formatter module.",
            required=True,
        ),
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
        formatter_class=AlignedRichHelpFormatter,
    )

    defaults = load_configuration().modules.fine_mapping.engines.susie

    # All arguments added to this group WILL stay under this heading
    susie = parser.add_argument_group("SuSiE settings")

    susie.add_argument(
        "--L",
        type=int,
        default=argparse.SUPPRESS,
        metavar="N",
        help=help_with_default(
            "Maximum number of SuSiE credible sets per locus. "
            "Increasing this may increase runtime and memory usage",
            defaults.max_causal_components,
        ),
    )

    susie.add_argument(
        "--minimum-purity",
        type=float,
        default=argparse.SUPPRESS,
        metavar="CORRELATION",
        help=help_with_default(
            "Minimum absolute within-set correlation required for a SuSiE "
            "credible set",
            defaults.minimum_purity,
        ),
    )

    # ------------------------------------------------------------------
    # Resource / timeout controls
    # ------------------------------------------------------------------
    susie.add_argument(
        "--ld-timeout-seconds",
        type=int,
        default=argparse.SUPPRESS,
        metavar="SECONDS",
        help=help_with_default(
            "Maximum time (in seconds) allowed for PLINK LD-matrix computation "
            "per locus. Execution aborts for loci exceeding this limit",
            defaults.execution.ld_timeout_seconds,
        ),
    )

    susie.add_argument(
        "--susie-timeout-seconds",
        type=int,
        default=argparse.SUPPRESS,
        metavar="SECONDS",
        help=help_with_default(
            "Maximum time (in seconds) allowed for SuSiE model fitting per locus. "
            "Loci exceeding the limit are skipped with a warning",
            defaults.execution.susie_timeout_seconds,
        ),
    )

    susie.add_argument(
        "--susie-main-max-iter",
        type=int,
        default=argparse.SUPPRESS,
        metavar="N",
        help=help_with_default(
            "Maximum iterations for the primary SuSiE fit",
            defaults.fitting.main_max_iter,
        ),
    )
    susie.add_argument(
        "--susie-recovery-max-iter",
        type=int,
        default=argparse.SUPPRESS,
        metavar="N",
        help=help_with_default(
            "Maximum iterations for a recovery fit",
            defaults.recovery.max_iter,
        ),
    )
    susie.add_argument(
        "--susie-repaired-max-iter",
        type=int,
        default=argparse.SUPPRESS,
        metavar="N",
        help=help_with_default(
            "Maximum iterations for a fit using repaired and revalidated LD",
            defaults.recovery.repaired_max_iter,
        ),
    )
    susie.add_argument(
        "--susie-reduced-l-max",
        type=int,
        default=argparse.SUPPRESS,
        metavar="N",
        help=help_with_default(
            "Maximum causal components used by reduced-L recovery",
            defaults.recovery.reduced_l_max,
        ),
    )
    susie.add_argument(
        "--susie-ld-eigenvalue-tolerance",
        type=float,
        default=argparse.SUPPRESS,
        metavar="VALUE",
        help=help_with_default(
            "Negative-eigenvalue tolerance used to classify LD as non-PSD",
            defaults.ld_validation.eigenvalue_tolerance,
        ),
    )
    susie.add_argument(
        "--susie-ld-mismatch-warning-threshold",
        type=float,
        default=argparse.SUPPRESS,
        metavar="LAMBDA",
        help=help_with_default(
            "SuSiE-RSS LD/z mismatch lambda that produces a warning",
            defaults.ld_validation.mismatch_warning_threshold,
        ),
    )
    susie.add_argument(
        "--susie-ld-mismatch-failure-threshold",
        type=float,
        default=argparse.SUPPRESS,
        metavar="LAMBDA",
        help=help_with_default(
            "SuSiE-RSS LD/z mismatch lambda that invalidates a locus",
            defaults.ld_validation.mismatch_failure_threshold,
        ),
    )
    susie.add_argument(
        "--susie-ld-repair-maximum-change",
        type=float,
        default=argparse.SUPPRESS,
        metavar="CORRELATION",
        help=help_with_default(
            "Largest absolute correlation change permitted during LD repair",
            defaults.ld_validation.repair_maximum_change,
        ),
    )
    susie.add_argument(
        "--susie-plink-timeout-seconds",
        type=int,
        default=argparse.SUPPRESS,
        metavar="SECONDS",
        help=help_with_default(
            "Hard timeout for each PLINK subprocess invoked by the R worker",
            defaults.execution.plink_timeout_seconds,
        ),
    )
    susie.add_argument(
        "--susie-memory-used-threshold-percent",
        type=float,
        default=argparse.SUPPRESS,
        metavar="PERCENT",
        help=help_with_default(
            "System memory-use percentage that activates the R memory throttle",
            defaults.memory.used_threshold_percent,
        ),
    )
    susie.add_argument(
        "--susie-memory-maximum-wait-seconds",
        type=int,
        default=argparse.SUPPRESS,
        metavar="SECONDS",
        help=help_with_default(
            "Maximum time an R worker waits for memory pressure to subside",
            defaults.memory.maximum_wait_seconds,
        ),
    )

    return parser


def get_finemap_finemap_arguments(add_help=False):
    """
    FINEMAP input files that exist ONLY in direct mode.
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    grp = parser.add_argument_group("FINEMAP input")

    # ===============================================================
    # Input / execution
    # ===============================================================
    grp.add_argument(
        "--finemap-in-files",
        metavar="PATH",
        type=validate_path(must_exist=True, must_be_file=True),
        default=argparse.SUPPRESS,
        help=style_cli_requirement(
            "FINEMAP input file created by `postgwas formatter --format "
            "finemap`.",
            required=True,
        ),
    )

    return parser


def get_common_finemap_finemap_arguments(add_help=False):
    """
    Add FINEMAP tuning arguments with 'sss' (Shotgun Stochastic Search) set as the default mode.
    Includes stylized default values and mode selection clarification.
    """
    parser = argparse.ArgumentParser(
        add_help=add_help,
        formatter_class=AlignedRichHelpFormatter,
    )

    finemap = parser.add_argument_group("FINEMAP settings")
    defaults = load_configuration().modules.fine_mapping.engines.finemap

    # ===============================================================
    # Subprograms (Modes)
    # ===============================================================
    algorithm = finemap.add_mutually_exclusive_group()
    algorithm.add_argument(
        "--sss",
        action="store_true",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Fine-mapping with [bold]shotgun stochastic search[/bold]. "
            "Use either [cyan]--sss[/cyan] or [cyan]--cond[/cyan], not both",
            defaults.algorithm,
        ),
    )

    algorithm.add_argument(
        "--cond",
        action="store_true",
        default=argparse.SUPPRESS,
        help=(
            "Fine-mapping with [bold]stepwise conditional search[/bold]. "
            "Use either [cyan]--sss[/cyan] or [cyan]--cond[/cyan], not both."
        ),
    )

    # ===============================================================
    # External-process supervision
    # ===============================================================
    finemap.add_argument(
        "--ldstore-timeout-seconds",
        type=int,
        metavar="SECONDS",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Maximum total runtime for LDstore processing within one locus",
            defaults.ldstore_timeout_seconds,
        ),
    )

    finemap.add_argument(
        "--finemap-timeout-seconds",
        type=int,
        metavar="SECONDS",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Maximum FINEMAP model-fitting runtime for one locus",
            defaults.finemap_timeout_seconds,
        ),
    )

    finemap.add_argument(
        "--finemap-termination-grace-seconds",
        type=int,
        metavar="SECONDS",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Grace period after timeout before the external process tree is "
            "force-killed",
            defaults.termination_grace_seconds,
        ),
    )

    # ===============================================================
    # Shotgun stochastic search (SSS)
    # ===============================================================
    finemap.add_argument(
        "--n-iter",
        metavar="N",
        type=int,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Maximum number of SSS iterations",
            defaults.n_iter,
        ),
    )

    finemap.add_argument(
        "--n-conv-sss",
        type=int,
        metavar="N",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Consecutive iterations for which added posterior mass must stay "
            "below --prob-conv-sss-tol",
            defaults.n_conv_sss,
        ),
    )

    finemap.add_argument(
        "--prob-conv-sss-tol",
        metavar="VALUE",
        type=float,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Added-posterior-mass tolerance used to terminate SSS",
            defaults.prob_conv_sss_tol,
        ),
    )

    finemap.add_argument(
        "--n-configs-top",
        metavar="N",
        type=int,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Number of top causal configurations to save",
            defaults.n_configs_top,
        ),
    )

    # ===============================================================
    # Causal priors
    # ===============================================================
    finemap.add_argument(
        "--n-causal-snps",
        type=int,
        metavar="N",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Maximum, rather than exact, number of causal SNPs per locus",
            defaults.n_causal_snps,
        ),
    )

    finemap.add_argument(
        "--prior-k",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=(
            "Use prior probabilities for the number of causal SNPs from FINEMAP "
            "K files. This requires K-file input support."
        ),
    )

    finemap.add_argument(
        "--prior-std",
        type=float,
        metavar="VALUE",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Prior standard deviation for causal effect sizes",
            defaults.prior_std,
        ),
    )

    # ===============================================================
    # LD / correlation handling
    # ===============================================================
    finemap.add_argument(
        "--corr-config",
        metavar="VALUE",
        type=float,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Set a causal configuration's posterior probability to zero when it "
            "contains a SNP pair with absolute correlation above this threshold",
            defaults.corr_config,
        ),
    )

    finemap.add_argument(
        "--force-n-samples",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help="Allow BCOR LD and GWAS sample sizes to differ.",
    )

    # ===============================================================
    # Filtering / credibility
    # ===============================================================
    finemap.add_argument(
        "--pvalue-snps",
        metavar="P",
        type=float,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Marginal p-value threshold for including SNPs",
            defaults.pvalue_snps,
        ),
    )

    finemap.add_argument(
        "--cond-pvalue",
        metavar="P",
        type=float,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Significance threshold for conditional search",
            defaults.cond_pvalue,
        ),
    )

    finemap.add_argument(
        "--prob-cred-set",
        metavar="PROBABILITY",
        type=float,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Posterior coverage target for each credible set",
            load_configuration().modules.fine_mapping.credible_set_coverage,
        ),
    )

    finemap.add_argument(
        "--std-effects",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=(
            "Print posterior effect-size means and standard deviations for "
            "standardized dosages."
        ),
    )

    return parser
