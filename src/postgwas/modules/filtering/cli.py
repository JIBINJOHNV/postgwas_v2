"""Command-line interface for configured GWAS-VCF filtering."""

import argparse
import sys

from postgwas.cli.common import (
    get_bcftools_binary_parser,
    get_common_out_parser,
    get_common_sumstat_filter_parser,
    get_inputvcf_parser,
)
from postgwas.cli.compute import get_compute_parser
from postgwas.core.errors import ConfigurationError
from postgwas.core.ui import AlignedRichHelpFormatter, format_cli_examples


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
            get_common_sumstat_filter_parser(add_help=False),
            get_bcftools_binary_parser(add_help=False),
        ],
    )
    parser.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit.",
    )
    parser.add_argument(
        "--run-config",
        metavar="PATH",
        help=(
            "YAML file containing filtering settings or a complete PostGWAS "
            "run configuration."
        ),
    )
    return parser


def main() -> int:
    parser = build_parser()
    if len(sys.argv) == 1:
        parser.print_help()
        return 0

    args = parser.parse_args()
    from postgwas.modules.filtering.service import (
        resolve_filtering_configuration,
        run_sumstat_filter_direct,
    )

    try:
        configuration = resolve_filtering_configuration(args)
    except ConfigurationError as exc:
        parser.error(str(exc))

    module = configuration.modules.filtering
    missing = []
    if module.inputs.vcf is None:
        missing.append("--vcf")
    if module.output_directory is None and configuration.run.output_directory is None:
        missing.append("--output-directory")
    if module.inputs.dataset_id is None and not configuration.run.dataset_id:
        missing.append("--dataset-id")
    if missing:
        parser.error("Missing required arguments: %s" % ", ".join(missing))

    run_sumstat_filter_direct(args, configuration=configuration)
    return 0


if __name__ == "__main__":
    main()
