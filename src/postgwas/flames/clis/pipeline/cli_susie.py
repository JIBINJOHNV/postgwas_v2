# postgwas/flames/clis/pipeline/cli_susie.py

from __future__ import annotations
from postgwas.utils.main import (
    validate_path,
    validate_plink_prefix,
    validate_alphanumeric,
)


def add_susie_arguments(parser):
    """
    Add SuSiE fine-mapping CLI arguments.
    Follows FLAMES OPTION-B style:
      - NO required=True (validated later)
      - Rich-friendly formatting
      - clean metavars
    """

    susie = parser.add_argument_group("3. SuSiE Fine-Mapping")

    # -------------------------------------------------------------
    # Core SuSiE configuration
    # -------------------------------------------------------------
    susie.add_argument(
        "--susie_L",
        type=int,
        default=10,
        help="Maximum number of credible sets per locus.",
        metavar=""
    )

    susie.add_argument(
        "--susie_min_pip",
        type=float,
        default=0.01,
        help="Minimum posterior inclusion probability (PIP).",
        metavar=""
    )

    susie.add_argument(
        "--skip_mhc",
        action="store_true",
        help="Skip the extended MHC region (chr6:25–34 Mb)."
    )

    # -------------------------------------------------------------
    # Required pipeline inputs (OPTION-B: validated later in code)
    # -------------------------------------------------------------
    susie.add_argument(
        "--locus_file",
        type=validate_path(must_exist=True, must_be_file=True),
        help="SuSiE locus definition file. [bold magenta]Required[/bold magenta]",
        metavar=""
    )


    susie.add_argument(
        "--plink",
        type=validate_path(must_exist=True, must_be_file=True),
        help="Path to PLINK executable. [bold magenta]Required[/bold magenta]",
        metavar=""
    )

    susie.add_argument(
        "--SUSIE_Analysis_folder",
        type=str,
        help="Output folder to store SuSiE results. [bold magenta]Required[/bold magenta]",
        metavar=""
    )

    # -------------------------------------------------------------
    # Optional SuSiE engine tuning
    # -------------------------------------------------------------
    susie.add_argument(
        "--lp_threshold",
        type=float,
        default=7.3,
        help="Log-P threshold to define fine-mapping loci.",
        metavar=""
    )

    susie.add_argument(
        "--L",
        type=int,
        default=10,
        help="Upper bound for number of SuSiE single-effect components.",
        metavar=""
    )

    susie.add_argument(
        "--workers",
        default="auto",
        help="Number of parallel workers (auto detects optimal).",
        metavar=""
    )

    susie.add_argument(
        "--min_ram_per_worker_gb",
        type=float,
        default=4,
        help="Minimum RAM per worker (GB).",
        metavar=""
    )

    susie.add_argument(
        "--timeout_ld_seconds",
        type=int,
        default=180,
        help="Timeout for LD matrix computation (seconds).",
        metavar=""
    )

    susie.add_argument(
        "--timeout_susie_seconds",
        type=int,
        default=180,
        help="Timeout for SuSiE fitting per locus (seconds).",
        metavar=""
    )

    susie.add_argument(
        "--mhc_start",
        type=int,
        default=25_000_000,
        help="Start coordinate of MHC region.",
        metavar=""
    )

    susie.add_argument(
        "--mhc_end",
        type=int,
        default=35_000_000,
        help="End coordinate of MHC region.",
        metavar=""
    )


    return parser
