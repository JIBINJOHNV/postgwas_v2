import argparse
import sys
from rich_argparse import RichHelpFormatter 
from pathlib import Path
import multiprocessing

# --- Imports ---
from postgwas.utils.main import validate_path
from postgwas.harmonisation.cli import get_harmonisation_parser, run_harmonisation
from postgwas.annot_ldblock.cli import get_annot_ldblock_parser, run_annot_ldblock_pipeline
from postgwas.sumstat_filter.sumstat_filter import filter_gwas_vcf_bcftools

from postgwas.sumstat_filter.workflows import run_sumstat_filter_pipeline, run_sumstat_filter_direct
from postgwas.clis.common_cli import (
    get_defaultresourse_parser,
    get_inputvcf_parser,
    get_common_out_parser,
    get_common_sumstat_filter_parser
)


# ---------------------------------------------------------
# Main CLI
# ---------------------------------------------------------
def main():
    top_parser = argparse.ArgumentParser(
        prog="sumstat_filter",
        description=(
            "sumstat_filter Module:\n\n"
            "This tool removes low-quality variants from GWAS summary statistics "
            "based on user-defined thresholds and filtering rules. "
            "Quality control steps include:\n"
            "  • Filtering variants with p-values\n"
            "  • Removing variants with low INFO scores\n"
            "  • Filtering by MAF / allele-frequency thresholds\n"
            "  • Removing highly palindromic SNPs (optional)\n"
            "  • Excluding variants with large AF mismatches to an external reference\n"
            "  • Optional removal of variants in the MHC region\n\n"
            "All QC thresholds are controlled by command-line flags. "
            "The resulting filtered summary statistics are written to the user’s "
            "output directory for use in downstream modules such as MAGMA, PoPS, "
            "fine-mapping, and FLAMES."
        ),
        formatter_class=RichHelpFormatter
    )

    subparsers = top_parser.add_subparsers(dest="mode", help="direct or pipeline")

    # ---------------- DIRECT ----------------
    direct_parser = subparsers.add_parser(
        "direct",
        help="Run QC only on an existing harmonised VCF.",
        parents=[
            get_defaultresourse_parser(),
            get_inputvcf_parser(),
            get_common_out_parser(),
            get_common_sumstat_filter_parser(add_help=False)
        ],
        formatter_class=RichHelpFormatter
    )
    direct_parser.set_defaults(func=run_sumstat_filter_direct)

    # ---------------- PIPELINE ----------------
    pipeline_parser = subparsers.add_parser(
        "pipeline",
        help="Run Harmonisation → Annot LDBlock → QC (Steps 1–3).",
        parents=[
            get_defaultresourse_parser(),
            get_harmonisation_parser(),
            get_annot_ldblock_parser(),
            get_common_sumstat_filter_parser(add_help=False)
        ],
        formatter_class=RichHelpFormatter
    )
    pipeline_parser.set_defaults(func=run_sumstat_filter_pipeline)

    # ---------------- EXECUTE ----------------
    args = top_parser.parse_args()

    # If a subcommand is chosen but no arguments are provided → show help instead of running
    if args.mode == "direct" and len(sys.argv) == 2:
        top_parser.parse_args(["direct", "--help"])
        sys.exit(0)

    if args.mode == "pipeline" and len(sys.argv) == 2:
        top_parser.parse_args(["pipeline", "--help"])
        sys.exit(0)
    
    if hasattr(args, "func"):
        args.func(args)
    else:
        top_parser.print_help()


if __name__ == "__main__":
    main()
