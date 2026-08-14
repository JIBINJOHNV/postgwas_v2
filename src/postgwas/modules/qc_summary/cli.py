#!/usr/bin/env python3
import sys
import argparse
from postgwas.cli.compute import get_compute_parser
from postgwas.core.ui import AlignedRichHelpFormatter, format_cli_examples

# ---------------------------------------------------------
# IMPORTS
# ---------------------------------------------------------
from postgwas.cli.common import (
    get_bcftools_binary_parser,
    get_common_out_parser,
    get_genome_build_parser,
    get_inputvcf_parser,
    sumstat_summary_arg_parser,
)
from postgwas.config import load_configuration


# =========================================================
# MAIN CLI
# =========================================================
def build_parser() -> argparse.ArgumentParser:
    configuration = load_configuration()
    parser = argparse.ArgumentParser(
        prog="postgwas qc",
        usage="postgwas qc --vcf PATH [options]",
        description="Calculate summary-level quality-control metrics for a GWAS-VCF.",
        epilog=format_cli_examples(
            (
                "Summarise one harmonised GWAS-VCF:",
                "postgwas qc",
                (
                    "--vcf study.vcf.gz",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
        ),
        parents=[
            get_compute_parser(),
            get_inputvcf_parser(),
            get_common_out_parser(),
            get_genome_build_parser(
                available_builds=list(configuration.resources.genomes),
                default_build=configuration.modules.qc_summary.target_build.value,
                suppress_default=True,
            ),
            sumstat_summary_arg_parser(),
            get_bcftools_binary_parser(),
        ],
        formatter_class=AlignedRichHelpFormatter,
    )
    parser.add_argument(
        "--run-config",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "YAML file containing QC-summary settings or a complete PostGWAS "
            "run configuration."
        ),
    )
    return parser


def main() -> int:
    parser = build_parser()

    # ---------------------------------------------------
    # NO ARGS → HELP
    # ---------------------------------------------------
    if len(sys.argv) == 1:
        parser.print_help()
        return 0

    args = parser.parse_args()

    from postgwas.core.errors import ConfigurationError
    from postgwas.modules.qc_summary.service import (
        resolve_qc_summary_configuration,
        run_qc_summary_direct,
    )

    try:
        configuration = resolve_qc_summary_configuration(args)
    except ConfigurationError as exc:
        parser.error(str(exc))

    module = configuration.modules.qc_summary
    missing = []
    if module.inputs.vcf is None:
        missing.append("--vcf")
    if module.inputs.dataset_id is None and not configuration.run.dataset_id:
        missing.append("--dataset-id")
    if module.output_directory is None and configuration.run.output_directory is None:
        missing.append("--output-directory")
    if missing:
        parser.error("Missing required arguments: %s" % ", ".join(missing))

    # ---------------------------------------------------
    # EXECUTION (DIRECT ONLY)
    # ---------------------------------------------------
    run_qc_summary_direct(args, configuration=configuration)
    return 0


if __name__ == "__main__":
    main()


__all__ = ["build_parser", "main"]
