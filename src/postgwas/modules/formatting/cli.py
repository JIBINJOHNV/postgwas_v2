"""Command-line interface for downstream format export."""

from __future__ import annotations

import argparse
import sys
from typing import get_args

from rich.console import Console

from postgwas.cli.common import (
    get_bcftools_binary_parser,
    get_common_out_parser,
    get_formatter_parser,
    get_inputvcf_parser,
    get_ldsc_merge_alleles_parser,
)
from postgwas.cli.compute import get_compute_parser
from postgwas.config.models.modules.formatting import FormattingCustomField
from postgwas.core.errors import ConfigurationError
from postgwas.core.ui import AlignedRichHelpFormatter, format_cli_examples
from postgwas.modules.formatting.service import run_formatter_direct
from postgwas.modules.formatting.table import FormattingError


class _CustomColumnAction(argparse.Action):
    """Collect custom columns in the exact order supplied on the command line."""

    def __init__(self, option_strings, dest, *, field, **kwargs):
        self.field = field
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        columns = dict(getattr(namespace, self.dest, {}))
        if self.field in columns:
            raise argparse.ArgumentError(
                self,
                "%s may be supplied only once" % option_string,
            )
        columns[self.field] = values
        setattr(namespace, self.dest, columns)


def _add_custom_output_arguments(parser: argparse.ArgumentParser) -> None:
    """Add the optional custom-table interface without CLI-owned defaults."""
    group = parser.add_argument_group(
        "custom output",
        "Write an additional table with user-named columns. Every requested "
        "field is required per row; invalid rows are excluded and counted.",
    )
    group.add_argument(
        "--custom-output",
        dest="custom_output_file",
        metavar="FILE",
        default=argparse.SUPPRESS,
        help="Relative output filename for the additional custom table.",
    )
    field_help = {
        "id": "Variant ID selected by --variant-id-type (required).",
        "chr": "Chromosome.",
        "pos": "One-based position.",
        "ref": "Other/non-effect allele.",
        "alt": "Effect allele.",
        "beta": "Effect estimate.",
        "se": "Standard error.",
        "z": "Z statistic.",
        "lp": "Input -log10(P).",
        "p": "Raw P derived from -log10(P).",
        "eaf": "Effect-allele frequency.",
        "maf": "Minor-allele frequency derived as min(EAF, 1-EAF).",
        "n": "Total sample size.",
        "neff": "Effective sample size.",
        "n_case": "Case sample size.",
        "n_control": "Control sample size.",
        "info": "Imputation INFO score.",
    }
    for field in get_args(FormattingCustomField):
        group.add_argument(
            "--%s" % field.replace("_", "-"),
            dest="custom_columns",
            metavar="NAME",
            action=_CustomColumnAction,
            field=field,
            default=argparse.SUPPRESS,
            help=field_help[field],
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="postgwas formatter",
        usage=(
            "postgwas formatter --vcf PATH --output-directory PATH "
            "[--format FORMAT [FORMAT ...]] [--custom-output FILE] [options]"
        ),
        description=(
            "Create validated MAGMA, SuSiE, FINEMAP, PRED-LD, "
            "LDSC, or MiXeR inputs, plus an optional custom table, from one "
            "harmonised GWAS-VCF. The VCF is read once even when several "
            "outputs are requested."
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
                "Export a reusable LDSC-only formatter configuration:",
                "postgwas config export",
                (
                    "--module formatting",
                    "--format ldsc",
                    "--style minimal",
                    "--output formatting.yaml",
                ),
            ),
            (
                "Create an additional custom table without editing YAML:",
                "postgwas formatter",
                (
                    "--vcf study.vcf.gz",
                    "--dataset-id STUDY",
                    "--output-directory results",
                    "--custom-output study_custom.tsv",
                    "--id SNP --chr CHR --pos BP --alt A1 --ref A2 --p P",
                ),
            ),
        ),
        formatter_class=AlignedRichHelpFormatter,
        parents=[
            get_compute_parser(),
            get_inputvcf_parser(),
            get_formatter_parser(direct_controls=True),
            get_common_out_parser(),
            get_bcftools_binary_parser(),
            get_ldsc_merge_alleles_parser(),
        ],
    )
    _add_custom_output_arguments(parser)
    for action in parser._actions:
        if action.dest in {"vcf", "output_directory"}:
            action.required = True
        if action.dest in {"bcftools", "dataset_id", "output_directory"}:
            action.default = argparse.SUPPRESS
    return parser


def main(argv=None):
    parser = build_parser()
    arguments = sys.argv[1:] if argv is None else list(argv)
    if not arguments:
        parser.print_help()
        return 0
    args = parser.parse_args(arguments)
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
