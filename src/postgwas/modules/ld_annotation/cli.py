import argparse
import sys

from postgwas.cli.common import (
    get_annot_ldblock_parser,
    get_common_out_parser,
    get_inputvcf_parser,
)
from postgwas.cli.compute import get_compute_parser
from postgwas.core.errors import ConfigurationError
from postgwas.core.ui import (
    AlignedRichHelpFormatter,
    format_cli_examples,
    mark_cli_required_help,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="postgwas annot_ldblock",
        usage="postgwas annot_ldblock --vcf PATH --ld-region-dir PATH [options]",
        description=(
            "Annotate a harmonised GWAS-VCF with population-specific LD blocks."
        ),
        epilog=format_cli_examples(
            (
                "Annotate one harmonised GWAS-VCF:",
                "postgwas annot_ldblock",
                (
                    "--vcf study.vcf.gz",
                    "--ld-region-dir reference/ld_blocks",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Run entirely from a reusable configuration:",
                "postgwas annot_ldblock",
                ("--run-config ld_annotation.yaml",),
            ),
            (
                "Export a reusable LD-annotation configuration:",
                "postgwas config export",
                (
                    "--module ld_annotation",
                    "--style full",
                    "--output ld_annotation.yaml",
                ),
            ),
        ),
        formatter_class=AlignedRichHelpFormatter,
        add_help=False,
        parents=[
            get_compute_parser(),
            get_inputvcf_parser(),
            get_annot_ldblock_parser(add_help=False),
            get_common_out_parser(),
        ],
    )

    # Custom help
    parser.add_argument(
        "-h", "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit.",
    )
    parser.add_argument(
        "--run-config",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "YAML file containing LD-annotation settings or a complete "
            "PostGWAS run configuration."
        ),
    )
    mark_cli_required_help(
        parser,
        ("vcf", "ld_region_dir", "dataset_id", "output_directory"),
    )
    return parser


def main() -> int:
    parser = build_parser()

    # -----------------------------
    # If no args → show help
    # -----------------------------
    if len(sys.argv) == 1:
        parser.print_help()
        return 0

    args = parser.parse_args()
    from postgwas.modules.ld_annotation.service import (
        resolve_ld_annotation_configuration,
        run_annot_ldblock,
    )

    try:
        configuration = resolve_ld_annotation_configuration(args)
        run_annot_ldblock(args, configuration=configuration)
    except ConfigurationError as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    main()


__all__ = ["build_parser", "main"]
