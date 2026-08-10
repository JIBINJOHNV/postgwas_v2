import argparse
import sys
from postgwas.modules.ld_annotation.service import run_annot_ldblock
from postgwas.cli.compute import get_compute_parser
from postgwas.core.ui import AlignedRichHelpFormatter, format_cli_examples

# --- IMPORTANT IMPORTS FROM STEP 1 (HARMONISATION) ---
from postgwas.cli.common import (
        get_inputvcf_parser,
        get_genome_build_parser,
        get_annot_ldblock_parser,
        get_common_out_parser)



def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="postgwas annot_ldblock",
        usage="postgwas annot_ldblock --vcf PATH --ld-region-dir PATH [options]",
        description="Annotate a harmonised GWAS-VCF with population-specific LD blocks.",
        epilog=format_cli_examples(
            (
                "Annotate one harmonised GWAS-VCF:",
                "postgwas annot_ldblock",
                (
                    "--vcf study.vcf.gz",
                    "--genome-build GRCh37",
                    "--ld-region-dir reference/ld_blocks",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
        ),
        formatter_class=AlignedRichHelpFormatter,
        add_help=False,
        parents=[
            get_compute_parser(),
            get_inputvcf_parser(),
            get_genome_build_parser(),
            get_annot_ldblock_parser(add_help=False),
            get_common_out_parser()
        ],
    )

    # Custom help
    parser.add_argument(
        "-h", "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit.",
    )
    return parser


def main():
    parser = build_parser()

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


__all__ = ["build_parser", "main"]
