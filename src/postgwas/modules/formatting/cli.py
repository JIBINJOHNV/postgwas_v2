"""Command-line interface for downstream format export."""

from __future__ import annotations

import argparse

from rich.console import Console
from postgwas.cli.common import (
    get_bcftools_binary_parser,
    get_common_out_parser,
    get_formatter_parser,
    get_inputvcf_parser,
)
from postgwas.cli.compute import get_compute_parser
from postgwas.core.errors import ConfigurationError
from postgwas.core.ui import AlignedRichHelpFormatter, format_cli_examples
from postgwas.modules.formatting.service import run_formatter_direct
from postgwas.modules.formatting.table import FormattingError


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="postgwas formatter",
        usage="postgwas formatter --vcf PATH [--format FORMAT [FORMAT ...]] [options]",
        description=(
            "Create scientifically validated MAGMA, SuSiE, FINEMAP, PRED-LD, "
            "LDSC, or MiXeR inputs from one harmonised GWAS-VCF. The VCF is read once "
            "even when several formats are requested."
        ),
        epilog=format_cli_examples(
            (
                "Create MAGMA and LDSC inputs:",
                "postgwas formatter",
                (
                    "--vcf study.vcf.gz",
                    "--dataset-id STUDY",
                    "--output-directory results",
                    "--format magma ldsc",
                ),
            ),
            (
                "Read output formats and parameters from YAML:",
                "postgwas formatter",
                (
                    "--vcf study.vcf.gz",
                    "--dataset-id STUDY",
                    "--output-directory results",
                    "--run-config formatting.yaml",
                ),
            ),
            (
                "Create coordinate-position-REF-ALT unique IDs:",
                "postgwas formatter",
                (
                    "--vcf study.vcf.gz",
                    "--dataset-id STUDY",
                    "--output-directory results",
                    "--format magma",
                    "--variant-id-type unique",
                ),
            ),
            (
                "Export a reusable formatter configuration:",
                "postgwas config export",
                ("--module formatting", "--style minimal", "--output formatting.yaml"),
            ),
        ),
        formatter_class=AlignedRichHelpFormatter,
        parents=[
            get_compute_parser(),
            get_inputvcf_parser(),
            get_formatter_parser(direct_controls=True),
            get_common_out_parser(),
            get_bcftools_binary_parser(),
        ],
    )
    for action in parser._actions:
        if action.dest == "vcf":
            action.required = True
        elif action.dest in {"bcftools", "dataset_id", "output_directory"}:
            action.default = argparse.SUPPRESS
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        run_formatter_direct(args)
    except (FormattingError, ConfigurationError, OSError, ValueError) as exc:
        Console(stderr=True).print(
            "\n[bold red]Formatter failed.[/bold red] %s\n" % exc
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = ["build_parser", "main"]
