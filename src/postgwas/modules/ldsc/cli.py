#!/usr/bin/env python3
"""
PostGWAS — LDSC Heritability CLI

Mode:
  • Direct – Munge formatter output and run single-trait LDSC heritability
"""
import argparse
import sys

from pydantic import ValidationError

from postgwas.cli.compute import get_compute_parser
from postgwas.config import load_configuration
from postgwas.core.errors import ConfigurationError
from postgwas.core.ui import (
    AlignedRichHelpFormatter, format_cli_examples, print_screen_message,
)

# =========================================================
# BACKEND RUNNERS
# =========================================================
from postgwas.modules.ldsc.ldsc_runner import LDSCError
from postgwas.modules.ldsc.service import run_ldsc_direct


# Shared parsers
from postgwas.cli.common import (
    get_common_out_parser,
    get_ldsc_common_parser,
)

from postgwas.core.execution.runtime import validate_path


# =========================================================
# SHARED LDSC INPUT ARGUMENTS
# =========================================================
def get_ldsc_parser(add_help: bool = False) -> argparse.ArgumentParser:
    """
    Defines input arguments for LDSC direct execution.
    """
    parser = argparse.ArgumentParser(add_help=add_help)
    grp = parser.add_argument_group("LDSC input")

    grp.add_argument(
        "--ldsc-input",
        metavar="PATH",
        default=argparse.SUPPRESS,
        type=validate_path(
            must_exist=True,
            must_be_file=True,
            must_not_be_empty=True,
        ),
        help=(
            "Formatter-created summary-statistics table consumed by direct "
            "LDSC analysis."
        ),
    )

    controls = parser.add_argument_group("Configuration")
    controls.add_argument(
        "--run-config",
        metavar="PATH",
        default=argparse.SUPPRESS,
        type=validate_path(must_exist=True, must_be_file=True),
        help=(
            "YAML file containing LDSC or full-run settings. Explicit command-line "
            "values override matching YAML values."
        ),
    )
    return parser


# =========================================================
# MAIN CLI
# =========================================================
def build_parser() -> argparse.ArgumentParser:
    defaults = load_configuration()
    ldsc_export = defaults.modules.formatting.exports["ldsc"]
    if ldsc_export.output_file is None:
        raise ConfigurationError(
            "modules.formatting.exports.ldsc.output_file must be configured"
        )
    ldsc_input_name = ldsc_export.output_file.format(dataset_id="STUDY")
    parser = argparse.ArgumentParser(
        prog="postgwas heritability",
        usage=(
            "postgwas heritability --ldsc-input PATH --merge-alleles PATH "
            "--ref-ld-chr PATH --w-ld-chr PATH --dataset-id NAME "
            "--output-directory PATH [options]"
        ),
        description=(
            "Run single-trait LDSC to estimate observed-scale SNP heritability "
            "and, when prevalence values are available, liability-scale "
            "heritability."
        ),
        epilog=format_cli_examples(
            (
                "Create the LDSC input first:",
                "postgwas formatter",
                (
                    "--vcf study.vcf.gz",
                    "--format ldsc",
                    "--merge-alleles reference/w_hm3.snplist",
                    "--dataset-id STUDY",
                    "--output-directory formatted",
                ),
            ),
            (
                "Run observed-scale LDSC heritability:",
                "postgwas heritability",
                (
                    "--ldsc-input formatted/%s" % ldsc_input_name,
                    "--merge-alleles reference/w_hm3.snplist",
                    "--ref-ld-chr reference/eur_w_ld_chr",
                    "--w-ld-chr reference/eur_w_ld_chr",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Run observed- and liability-scale LDSC heritability:",
                "postgwas heritability",
                (
                    "--ldsc-input formatted/%s" % ldsc_input_name,
                    "--merge-alleles reference/w_hm3.snplist",
                    "--ref-ld-chr reference/eur_w_ld_chr",
                    "--w-ld-chr reference/eur_w_ld_chr",
                    "--samp-prev 0.2",
                    "--pop-prev 0.01",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Export reusable LDSC settings:",
                "postgwas config export",
                ("--module ldsc", "--style full", "--output ldsc.yaml"),
            ),
            (
                "Run with exported LDSC settings and explicit input resources:",
                "postgwas heritability",
                (
                    "--run-config ldsc.yaml",
                    "--ldsc-input formatted/%s" % ldsc_input_name,
                    "--merge-alleles reference/w_hm3.snplist",
                    "--ref-ld-chr reference/eur_w_ld_chr",
                    "--w-ld-chr reference/eur_w_ld_chr",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
        ),
        parents=[
            get_compute_parser(),
            get_common_out_parser(),
            get_ldsc_parser(),
            get_ldsc_common_parser(add_help=False),
        ],
        formatter_class=AlignedRichHelpFormatter,
    )
    for destination in (
        "ldsc_input",
        "merge_alleles",
        "ref_ld_chr",
        "w_ld_chr",
        "dataset_id",
        "output_directory",
    ):
        action = next(
            item for item in parser._actions if item.dest == destination
        )
        action.required = True
    return parser


def get_ldsc_pipeline_examples():
    """Return complete LDSC examples for context-sensitive pipeline help."""
    required = (
        "--modules heritability",
        "--vcf study.vcf.gz",
        "--merge-alleles reference/w_hm3.snplist",
        "--ref-ld-chr reference/eur_w_ld_chr",
        "--w-ld-chr reference/eur_w_ld_chr",
        "--dataset-id STUDY",
        "--output-directory results",
    )
    return (
        (
            "Run observed-scale LDSC heritability from a GWAS-VCF:",
            "postgwas pipeline",
            required,
        ),
        (
            "Add liability-scale heritability using the GWAS-VCF case fraction:",
            "postgwas pipeline",
            required + ("--pop-prev 0.01",),
        ),
        (
            "Override the GWAS-VCF case fraction for liability-scale heritability:",
            "postgwas pipeline",
            required + ("--samp-prev 0.2", "--pop-prev 0.01"),
        ),
    )


def main(argv=None) -> int:
    parser = build_parser()

    # ---------------------------------------------------
    # EXECUTION
    # ---------------------------------------------------

    # Check if no arguments provided (sys.argv[0] is the script name)
    arguments = sys.argv[1:] if argv is None else list(argv)
    if not arguments:
        parser.print_help()
        return 0

    args = parser.parse_args(arguments)

    # Directly dispatch to the direct runner
    try:
        run_ldsc_direct(args)
    except (LDSCError, ConfigurationError, OSError, ValidationError) as exc:
        print_screen_message(
            "error", "LDSC heritability stopped.\n"
            "Reason: %s\nNo new LDSC results were published." % exc,
            stderr=True,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "build_parser", "get_ldsc_parser", "get_ldsc_pipeline_examples", "main",
]
