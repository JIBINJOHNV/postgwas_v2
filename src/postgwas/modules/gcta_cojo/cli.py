"""Direct GCTA-input CLI and pipeline-facing options for GCTA-COJO."""

from __future__ import annotations

import argparse
from typing import get_args

from postgwas.cli.common import get_common_out_parser
from postgwas.cli.compute import get_compute_parser
from postgwas.config import load_configuration
from postgwas.config.models.modules.gcta_cojo import GctaCojoMode
from postgwas.core.errors import ConfigurationError
from postgwas.core.execution.runtime import validate_path
from postgwas.core.ui import (
    AlignedRichHelpFormatter,
    format_cli_examples,
    help_with_conditional_requirement,
    help_with_default,
    mark_cli_required_help,
)
from postgwas.modules.gcta_cojo.errors import GctaCojoError


def _default(description: str, value) -> str:
    return help_with_default(description, value)


def get_gcta_cojo_parser(add_help=False, *, direct_controls=False):
    defaults = load_configuration()
    module = defaults.modules.gcta_cojo
    parser = argparse.ArgumentParser(add_help=add_help)
    inputs = parser.add_argument_group("GCTA-COJO inputs")
    if direct_controls:
        inputs.add_argument(
            "--cojo-file",
            dest="gcta_cojo_input_file",
            type=validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True),
            metavar="PATH",
            default=argparse.SUPPRESS,
            help=_default(
                "Existing GCTA summary-statistics .ma file with columns SNP, A1, "
                "A2, freq, b, se, p, and N",
                module.input_file,
            ),
        )
    else:
        inputs.add_argument(
            "--gcta-cojo-input-file",
            type=validate_path(
                must_exist=True,
                must_be_file=True,
                must_not_be_empty=True,
            ),
            metavar="PATH",
            default=argparse.SUPPRESS,
            help=argparse.SUPPRESS,
        )
    inputs.add_argument(
        "--cojo-reference-prefix",
        metavar="PREFIX",
        default=argparse.SUPPRESS,
        help=_default(
            "Prefix shared by the ancestry-matched PLINK BED, BIM, and FAM "
            "LD-reference files",
            module.reference.prefix,
        ),
    )
    inputs.add_argument(
        "--condition-snps",
        type=validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True),
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=_default(
            help_with_conditional_requirement(
                "Known SNPs to adjust for, with one BIM SNP ID per line. GCTA "
                "tests the other included SNPs conditional on these variants. "
                "Supplying this option selects cond mode when --cojo-mode is "
                "omitted; any explicit mode must be cond",
                "for cond mode",
            ),
            module.inputs.condition_snps,
        ),
    )
    inputs.add_argument(
        "--joint-snps",
        type=validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True),
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=_default(
            help_with_conditional_requirement(
                "SNPs whose effects are estimated together in one multiple-SNP "
                "model, with one BIM SNP ID per line. Supplying this option "
                "selects joint mode when --cojo-mode is omitted; any explicit "
                "mode must be joint",
                "for joint mode",
            ),
            module.inputs.joint_snps,
        ),
    )
    inputs.add_argument(
        "--cojo-extract",
        type=validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True),
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=_default(
            "Only test or select SNPs listed in this file, with one BIM SNP "
            "ID per line. GCTA still reads the full .ma file",
            module.inputs.extract_snps,
        ),
    )
    inputs.add_argument(
        "--cojo-exclude",
        type=validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True),
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=_default(
            "Do not test or select SNPs listed in this file, with one BIM SNP ID "
            "per line",
            module.inputs.exclude_snps,
        ),
    )

    compatibility = parser.add_argument_group("Input and reference compatibility")
    compatibility.add_argument(
        "--genome-build",
        choices=list(defaults.resources.genomes),
        metavar="BUILD",
        default=argparse.SUPPRESS,
        help=_default(
            "Shared genome build of the GWAS summary statistics and PLINK LD reference",
            module.genome_build,
        ),
    )
    compatibility.add_argument(
        "--cojo-reference-population",
        choices=list(defaults.resources.populations),
        metavar="CODE",
        default=argparse.SUPPRESS,
        help=_default(
            "Population represented by the LD reference; ancestry is never inferred",
            module.reference.population,
        ),
    )
    compatibility.add_argument(
        "--cojo-minimum-reference-overlap",
        type=float,
        metavar="FRACTION",
        default=argparse.SUPPRESS,
        help=_default(
            "Minimum fraction of GCTA .ma variants that must exactly match BIM IDs",
            module.input_validation.minimum_reference_overlap_fraction,
        ),
    )

    settings = parser.add_argument_group("GCTA-COJO analysis options")
    settings.add_argument(
        "--cojo-mode",
        choices=get_args(GctaCojoMode),
        default=argparse.SUPPRESS,
        help=_default(
            "Analysis mode: stepwise selection, fixed-count top SNPs, joint "
            "estimation, or conditional association",
            module.mode,
        ),
    )
    settings.add_argument(
        "--cojo-p",
        type=float,
        metavar="P",
        default=argparse.SUPPRESS,
        help=_default(
            "Genome-wide significance threshold used by cojo-slct",
            module.analysis.significance_threshold,
        ),
    )
    settings.add_argument(
        "--cojo-top-snps",
        type=int,
        metavar="N",
        default=argparse.SUPPRESS,
        help=_default(
            "Fixed number of independent signals selected by GCTA",
            module.analysis.top_snp_count,
        ),
    )
    settings.add_argument(
        "--cojo-window-kb",
        type=int,
        metavar="KB",
        default=argparse.SUPPRESS,
        help=_default(
            "Distance beyond which GCTA assumes complete linkage equilibrium",
            module.analysis.window_kb,
        ),
    )
    settings.add_argument(
        "--cojo-collinear",
        type=float,
        metavar="R2",
        default=argparse.SUPPRESS,
        help=_default(
            "Multiple-regression R-squared collinearity cutoff",
            module.analysis.collinearity_cutoff,
        ),
    )
    settings.add_argument(
        "--cojo-diff-freq",
        type=float,
        metavar="FRACTION",
        default=argparse.SUPPRESS,
        help=_default(
            "Maximum A1-frequency difference between GWAS and LD reference",
            module.analysis.frequency_difference_max,
        ),
    )
    settings.add_argument(
        "--cojo-maf",
        type=float,
        metavar="FRACTION",
        default=argparse.SUPPRESS,
        help=_default(
            "Minimum MAF retained from the LD reference",
            module.analysis.reference_maf_min,
        ),
    )
    settings.add_argument(
        "--cojo-gc",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=_default(
            "Apply GCTA genomic control to summary-statistic p-values",
            module.analysis.genomic_control,
        ),
    )
    settings.add_argument(
        "--cojo-gc-lambda",
        type=float,
        metavar="LAMBDA",
        default=argparse.SUPPRESS,
        help=_default(
            "Explicit genomic-inflation factor; supplying it enables --cojo-gc",
            module.analysis.genomic_control_lambda,
        ),
    )
    settings.add_argument(
        "--cojo-chromosome",
        metavar="CHR",
        default=argparse.SUPPRESS,
        help=_default(
            "Optional chromosome passed to GCTA --chr",
            module.analysis.chromosome,
        ),
    )
    settings.add_argument(
        "--cojo-top-results",
        type=int,
        metavar="N",
        default=argparse.SUPPRESS,
        help=_default(
            "Number of lowest adjusted-p findings displayed at completion",
            module.reporting.top_result_count,
        ),
    )
    settings.add_argument(
        "--cojo-finding-threshold",
        type=float,
        metavar="P",
        default=argparse.SUPPRESS,
        help=_default(
            "Adjusted-p threshold used only for the terminal findings summary",
            module.reporting.finding_threshold,
        ),
    )
    settings.add_argument(
        "--cojo-p-value-digits",
        type=int,
        metavar="N",
        default=argparse.SUPPRESS,
        help=_default(
            "Significant digits used for p-values in the terminal summary",
            module.reporting.p_value_significant_digits,
        ),
    )
    settings.add_argument(
        "--gcta",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=_default(
            "GCTA executable path or command name",
            defaults.resources.executables.gcta,
        ),
    )
    if direct_controls:
        settings.add_argument(
            "--run-config",
            metavar="PATH",
            default=argparse.SUPPRESS,
            help="YAML settings; explicit CLI values override matching keys.",
        )
        settings.add_argument(
            "--dry-run",
            action="store_true",
            default=argparse.SUPPRESS,
            help="Validate existing COJO inputs and record the command without running GCTA.",
        )
    return parser


def get_gcta_cojo_pipeline_examples():
    return (
        (
            "Select independent signals when no conditioning SNPs are supplied:",
            "postgwas pipeline",
            (
                "--modules gcta_cojo",
                "--vcf study.vcf.gz",
                "--cojo-reference-prefix reference/EUR_ld",
                "--genome-build GRCh37",
                "--cojo-reference-population EUR",
                "--dataset-id STUDY",
                "--output-directory results",
            ),
        ),
        (
            "Run conditional association on an explicit SNP list:",
            "postgwas pipeline",
            (
                "--modules gcta_cojo",
                "--vcf study.vcf.gz",
                "--condition-snps lead_snps.txt",
                "--cojo-reference-prefix reference/EUR_ld",
                "--genome-build GRCh37",
                "--cojo-reference-population EUR",
                "--dataset-id STUDY",
                "--output-directory results",
            ),
        ),
    )


def build_parser():
    defaults = load_configuration()
    parser = argparse.ArgumentParser(
        prog="postgwas gcta_cojo",
        usage=(
            "postgwas gcta_cojo --cojo-file PATH "
            "--cojo-reference-prefix PREFIX [options]"
        ),
        description=(
            "Run GCTA-COJO stepwise, fixed-count, joint, or conditional analysis "
            "from an existing GCTA .ma file using an ancestry- and build-matched "
            "PLINK LD reference."
        ),
        epilog=format_cli_examples(
            (
                "Select independent signals (default when no SNP list is supplied):",
                "postgwas gcta_cojo",
                (
                    "--cojo-file study.ma",
                    "--cojo-reference-prefix reference/EUR_ld",
                    "--genome-build GRCh37",
                    "--cojo-reference-population EUR",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Condition on known lead variants:",
                "postgwas gcta_cojo",
                (
                    "--cojo-file study.ma",
                    "--condition-snps lead_snps.txt",
                    "--cojo-reference-prefix reference/EUR_ld",
                    "--genome-build GRCh37",
                    "--cojo-reference-population EUR",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Export reusable module settings:",
                "postgwas config export",
                ("--module gcta_cojo", "--style full", "--output gcta_cojo.yaml"),
            ),
        ),
        formatter_class=AlignedRichHelpFormatter,
        parents=[
            get_compute_parser(),
            get_common_out_parser(),
            get_gcta_cojo_parser(direct_controls=True),
        ],
    )
    for action in parser._actions:
        if action.dest == "dataset_id":
            action.help = _default(
                "Short, unique name used in output filenames",
                defaults.run.dataset_id,
            )
            action.default = argparse.SUPPRESS
        elif action.dest == "output_directory":
            action.help = _default(
                "Folder where PostGWAS saves COJO results, metadata, and logs",
                defaults.run.output_directory,
            )
            action.default = argparse.SUPPRESS
    mark_cli_required_help(
        parser,
        (
            "gcta_cojo_input_file",
            "cojo_reference_prefix",
            "genome_build",
            "cojo_reference_population",
        ),
    )
    return parser


def main(argv=None):
    from postgwas.modules.gcta_cojo.service import run_gcta_cojo_direct

    args = build_parser().parse_args(argv)
    try:
        run_gcta_cojo_direct(args)
    except (GctaCojoError, ConfigurationError, OSError, ValueError):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "build_parser", "get_gcta_cojo_parser", "get_gcta_cojo_pipeline_examples", "main",
]
