import argparse
import sys
from rich_argparse import RichHelpFormatter 
# Assuming utility/logic imports
from postgwas.utils.main import validate_path 
from postgwas.annot_ldblock.workflows import run_annot_ldblock

# --- IMPORTANT IMPORTS FROM STEP 1 (HARMONISATION) ---
from postgwas.clis.common_cli import (
        get_defaultresourse_parser,
        get_inputvcf_parser,
        get_genomeversion_parser,
        get_annot_ldblock_parser,
        get_common_out_parser)



def main():
    # Top-level parser (direct-only CLI)
    parser = argparse.ArgumentParser(
        prog="annot_ldblock",
        description="Run LD Block Annotation using VCF output from harmonisation.",
        formatter_class=RichHelpFormatter,
        add_help=False,
        parents=[
            get_defaultresourse_parser(),
            get_inputvcf_parser(),
            get_genomeversion_parser(),
            get_annot_ldblock_parser(add_help=False),
            get_common_out_parser()
        ],
    )

    # Custom help
    parser.add_argument(
        "-h", "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show help for annot_ldblock"
    )

    # -----------------------------
    # If no args → show help
    # -----------------------------
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    # Parse arguments
    args = parser.parse_args()

    # -----------------------------
    # Validate required arguments
    # -----------------------------
    missing = []
    if not hasattr(args, "vcf") or args.vcf is None:
        missing.append("--vcf")
    if not hasattr(args, "ld_region_dir") or args.ld_region_dir is None:
        missing.append("--ld-region-dir")
    if missing:
        print("\n❗ Missing required arguments: " + ", ".join(missing) + "\n")
        parser.print_help()
        sys.exit(1)

    # -----------------------------
    # Run the workflow
    # -----------------------------
    run_annot_ldblock(args)


if __name__ == "__main__":
    main()
