# postgwas/flames/clis/cli_common.py

from __future__ import annotations
import argparse
from pathlib import Path

from rich_argparse import ArgumentDefaultsRichHelpFormatter
from postgwas.utils.main import auto_detect_workers, validate_path, validate_alphanumeric


def add_common_options(parser: argparse.ArgumentParser):
    """
    Add FLAMES global and model-level CLI options.
    Clean formatting, no metavars, Rich styling.
    """

    parser.formatter_class = ArgumentDefaultsRichHelpFormatter

    threads_default = auto_detect_workers()

    # ────────────────────────────────
    # Global Settings
    # ────────────────────────────────
    global_group = parser.add_argument_group("Global Settings")

    global_group.add_argument(
        "--out_dir",
        type=str,
        default=str(Path.cwd() / "flames_analysis"),
        help="Output directory for FLAMES results.",
        metavar=""
    )

    global_group.add_argument(
        "--sample_id",
        type=validate_alphanumeric,
        help="Sample identifier / GWAS phenotype name. ([bold magenta]Required[/bold magenta])",
        metavar=""
    )

    global_group.add_argument(
        "--threads",
        type=int,
        default=threads_default,
        help="Number of worker threads.",
        metavar=""
    )

    global_group.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing outputs."
    )

    global_group.add_argument(
        "--log_file",
        type=str,
        help="Optional log file to save run details.",
        metavar=""
    )

    # ────────────────────────────────
    # FLAMES Model Options
    # ────────────────────────────────
    model_group = parser.add_argument_group("FLAMES Model Options")

    model_group.add_argument(
        "--distance",
        type=int,
        default=750000,
        help="Gene distance (bp) for locus-to-gene mapping.",
        metavar=""
    )

    model_group.add_argument(
        "--weight",
        help="XGB weight to use, default is 0.725 XGBoost, 0.275 PoPS",
        required=False,
         metavar="",
        default=0.725,
    )
    
    model_group.add_argument(
        "--modelpath",
        type=validate_path(must_exist=True, must_be_dir=True),
        default=str(Path(__file__).resolve().parent.parent / "model"),
        help="Directory containing pretrained FLAMES model.",
        metavar=""
    )

    model_group.add_argument(
        "--annotation_dir",
        type=validate_path(must_exist=True, must_be_dir=True),
        help="Directory containing the flame annotation files",
        metavar="",
    )
    
    model_group.add_argument(
        "--true_positives",
        type=validate_path(must_exist=True, must_be_dir=True),
        help="Path to the file that describes the true positive genes in each locus",
        metavar="",
        default = False
    )

    model_group.add_argument( "--cmd_vep", help="path to the vep executable", metavar="", default=False )
    
    model_group.add_argument( "--vep_cache", help="path to the vep cache",  metavar="", default=False  )
    

    model_group.add_argument(
       "--tabix", help="path to tabix if using local CADD scores",  metavar="", default=False )
    
    model_group.add_argument(
        "--CADD_file", metavar="",
        help="path to CADD scores file if using local CADD scores, must match build of inputted credible variants",
        default=False )
    
    
    return parser
