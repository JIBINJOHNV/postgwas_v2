import argparse
import sys
from rich_argparse import RichHelpFormatter 
from pathlib import Path
import multiprocessing

# --- Imports ---
from postgwas.utils.main import validate_path
from postgwas.harmonisation.cli import get_harmonisation_parser, run_harmonisation
from postgwas.annot_ldblock.cli import get_annot_ldblock_parser
from postgwas.sumstat_filter.sumstat_filter import filter_gwas_vcf_bcftools

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

    # Top-level direct-only CLI
    parser = argparse.ArgumentParser(
        prog="sumstat_filter",
        description=(
            "GWAS Summary Statistics QC Module\n\n"
            "This tool removes low-quality variants from harmonised GWAS VCF files. "
            "QC operations include:\n"
            "  • P-value filtering\n"
            "  • INFO score thresholds\n"
            "  • MAF / allele-frequency thresholds\n"
            "  • Optional palindromic variant removal\n"
            "  • Optional allele-frequency mismatch filtering\n"
            "  • Optional exclusion of the MHC region\n\n"
            "Produces cleaned summary statistics for downstream analyses "
            "(MAGMA, SuSiE, PoPS, FLAMES, etc.)."
        ),
        formatter_class=RichHelpFormatter,
        add_help=False,
        parents=[
            get_defaultresourse_parser(),
            get_inputvcf_parser(),
            get_common_out_parser(),
            get_common_sumstat_filter_parser(add_help=False)
        ],
    )

    # Custom -h at the top
    parser.add_argument(
        "-h", "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show help for sumstat_filter (direct QC mode)."
    )

    # ---------------------------------------------------
    # If no arguments → show help
    # ---------------------------------------------------
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    # Parse arguments
    args = parser.parse_args()

    # ---------------------------------------------------
    # Validate required inputs (direct mode)
    # ---------------------------------------------------
    missing = []

    # VCF is mandatory
    if not hasattr(args, "vcf") or args.vcf is None:
        missing.append("--vcf")

    # Output directory is mandatory
    if not hasattr(args, "outdir") or args.outdir is None:
        missing.append("--outdir")

    # Show help if required arguments missing
    if missing:
        print(f"\n❗ Missing required arguments: {', '.join(missing)}\n")
        parser.print_help()
        sys.exit(1)
    # ---------------------------------------------------
    # Run direct QC workflow
    # ---------------------------------------------------
    from postgwas.sumstat_filter.workflows import run_sumstat_filter_direct
    run_sumstat_filter_direct(args)


if __name__ == "__main__":
    main()
