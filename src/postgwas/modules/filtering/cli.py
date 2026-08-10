import argparse
import sys

from postgwas.config import load_module_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.cli.compute import get_compute_parser
from postgwas.core.ui import AlignedRichHelpFormatter, format_cli_examples

from postgwas.cli.common import (
    get_inputvcf_parser,
    get_common_out_parser,
    get_common_sumstat_filter_parser
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="postgwas sumstat_filter",
        usage="postgwas sumstat_filter --vcf PATH [options]",
        description=(
            "Apply configured quality-control filters to a harmonised GWAS-VCF.\n"
            "Filtering can assess P-values, INFO, allele frequency, variant type, "
            "strand ambiguity, reference-frequency concordance, and the MHC region."
        ),
        epilog=format_cli_examples(
            (
                "Filter one harmonised GWAS-VCF:",
                "postgwas sumstat_filter",
                (
                    "--vcf study.vcf.gz",
                    "--dataset-id STUDY",
                    "--output-directory results",
                    "--run-config filtering.yaml",
                ),
            ),
            (
                "Export a reusable filtering configuration:",
                "postgwas config export",
                ("--module filtering", "--style full", "--output filtering.yaml"),
            ),
        ),
        formatter_class=AlignedRichHelpFormatter,
        add_help=False,
        parents=[
            get_compute_parser(),
            get_inputvcf_parser(),
            get_common_out_parser(),
            get_common_sumstat_filter_parser(add_help=False)
        ],
    )

    parser.add_argument(
        "-h", "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit.",
    )
    parser.add_argument(
        "--run-config",
        metavar="PATH",
        help="YAML file containing filtering settings or a complete PostGWAS run configuration.",
    )
    return parser


def main():
    parser = build_parser()

    # ---------------------------------------------------
    # If no arguments → show help
    # ---------------------------------------------------
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    # Parse arguments
    args = parser.parse_args()

    # Only explicitly supplied CLI flags override YAML. Historical argparse
    # defaults are intentionally ignored by assigning the resolved values.
    supplied = set(sys.argv[1:])
    cli_map = {
        "--vcf": ("inputs.vcf", args.vcf),
        "--dataset-id": ("inputs.dataset_id", getattr(args, "dataset_id", None)),
        "--output-directory": ("output_directory", args.output_directory),
        "--minimum-neglog10-p": ("minimum_neglog10_p", args.minimum_neglog10_p),
        "--minimum-maf": ("maf_min", args.minimum_maf),
        "--reference-af-column": ("reference_population_tag", args.reference_af_column),
        "--maximum-af-difference": ("frequency_difference_max", args.maximum_af_difference),
        "--minimum-info": ("info_min", args.minimum_info),
        "--maximum-info": ("info_max", args.maximum_info),
        "--missing-info-action": ("missing_info_action", args.missing_info_action),
        "--include-indels": ("include_indels", args.include_indels),
        "--remove-palindromic": ("remove_palindromic", args.remove_palindromic),
        "--palindromic-af-lower": ("palindromic_lower", args.palindromic_af_lower),
        "--palindromic-af-upper": ("palindromic_upper", args.palindromic_af_upper),
        "--remove-mhc": ("remove_mhc", args.remove_mhc),
        "--mhc-chrom": ("mhc.chromosome", args.mhc_chrom),
        "--mhc-start": ("mhc.start", args.mhc_start),
        "--mhc-end": ("mhc.end", args.mhc_end),
    }
    overrides = {path: value for flag, (path, value) in cli_map.items() if flag in supplied}
    try:
        config = load_module_configuration(
            "filtering", args.run_config, cli_overrides=overrides
        )
    except ConfigurationError as exc:
        parser.error(str(exc))

    args.vcf = args.vcf if "--vcf" in supplied else config.inputs.vcf
    args.dataset_id = (
        args.dataset_id
        if "--dataset-id" in supplied
        else config.inputs.dataset_id
    )
    args.output_directory = (
        args.output_directory
        if "--output-directory" in supplied
        else config.output_directory
    )
    args.minimum_neglog10_p = config.minimum_neglog10_p
    args.minimum_maf = config.maf_min
    args.reference_af_column = config.reference_population_tag
    args.maximum_af_difference = config.frequency_difference_max
    args.minimum_info = config.info_min
    args.maximum_info = config.info_max
    args.missing_info_action = config.missing_info_action
    args.include_indels = config.include_indels
    args.remove_palindromic = config.remove_palindromic
    args.palindromic_af_lower = config.palindromic_lower
    args.palindromic_af_upper = config.palindromic_upper
    args.remove_mhc = config.remove_mhc
    args.mhc_chrom = config.mhc.chromosome
    args.mhc_start = config.mhc.start
    args.mhc_end = config.mhc.end
    args.module_config = config

    # ---------------------------------------------------
    # Validate required inputs (direct mode)
    # ---------------------------------------------------
    missing = []

    # VCF is mandatory
    if not hasattr(args, "vcf") or args.vcf is None:
        missing.append("--vcf")

    # Output directory is mandatory
    if not hasattr(args, "output_directory") or args.output_directory is None:
        missing.append("--output-directory")

    if not hasattr(args, "dataset_id") or args.dataset_id is None:
        missing.append("--dataset-id")

    # Show help if required arguments missing
    if missing:
        print(f"\n❗ Missing required arguments: {', '.join(missing)}\n")
        parser.print_help()
        sys.exit(1)
    # ---------------------------------------------------
    # Run direct QC workflow
    # ---------------------------------------------------
    from postgwas.modules.filtering.service import run_sumstat_filter_direct
    run_sumstat_filter_direct(args)


if __name__ == "__main__":
    main()
