"""Command-line interface for CALDERA."""

from __future__ import annotations

import argparse

from rich.console import Console

from postgwas.cli.common import get_common_out_parser
from postgwas.cli.compute import get_compute_parser
from postgwas.config import load_configuration
from postgwas.core.errors import ConfigurationError, MissingRequiredArgumentsError
from postgwas.core.ui import (
    print_screen_message,
    AlignedRichHelpFormatter,
    format_cli_examples,
    help_with_default,
    mark_cli_required_help,
    move_cli_help_actions,
)
from postgwas.modules.caldera.errors import CalderaError


_CALDERA_OUTPUT_DESTINATIONS = (
    "dataset_id",
    "output_directory",
    "save_matrix_files",
)
_CALDERA_INPUT_FILE_DESTINATIONS = (
    "vcf",
    "pops_file",
    "credible_set_file",
    "magma_association_prefix",
    "target_score_file",
    "target_covariates_file",
    "target_error_covariance_file",
)
_CALDERA_PIPELINE_GENERATED_INPUT_DESTINATIONS = (
    "target_score_file",
    "target_covariates_file",
    "target_error_covariance_file",
)
_CALDERA_REFERENCE_FILE_DESTINATIONS = (
    "magma_ld_reference",
    "gene_location_file",
    "gene_set_file",
    "feature_matrix_prefix",
    "pops_gene_location_file",
    "feature_subset_file",
    "control_features_file",
    "ld_region_dir",
    "cojo_reference_prefix",
    "ld_folder",
    "finemap_ld_reference",
    "caldera_repository",
)
_CALDERA_REFERENCE_SETTING_DESTINATIONS = (
    "magma_positional_gene_id_type",
    "magma_positional_source_name",
    "magma_positional_source_version",
    "magma_positional_source_url",
    "magma_positional_context",
    "caldera_genome_build",
)
_POPS_RESOURCE_SETTING_DESTINATIONS = (
    "feature_matrix_chunks",
    "gene_universe_policy",
    "pops_genome_build",
)
_CALDERA_SOFTWARE_DESTINATIONS = ("caldera_adapter_script",)
_EXTERNAL_SOFTWARE_DESTINATIONS = ("tabix", "gcta")
_CALDERA_HELP_GROUP_ORDER = (
    "Output",
    "Input files",
    "Reference files",
    "CALDERA reference settings",
    "Variant identifiers",
    "MAGMA analysis settings",
    "MAGMA analysis scope",
    "PoPS resource settings",
    "PoPS covariate settings",
    "PoPS feature-selection settings",
    "PoPS training settings",
    "PoPS model settings",
    "LD-block annotation",
    "Genome coordinates",
    "Reference population",
    "Input and reference compatibility",
    "GCTA-COJO stepwise selection",
    "LD-clumping analysis",
    "Standard r² clumping",
    "GCTA-COJO physical locus definition",
    "Fine-mapping settings",
    "SuSiE settings",
    "FINEMAP settings",
    "External software",
    "CALDERA software",
    "Performance",
    "Screen output",
    "Run continuation",
    "Configuration",
    "Choose analyses",
)


def organize_caldera_help(parser: argparse.ArgumentParser) -> None:
    """Place outputs, analysis inputs, and reference resources before settings."""
    original_groups = list(parser._action_groups)
    pipeline_mode = any(action.dest == "vcf" for action in parser._actions)
    if pipeline_mode:
        for action in parser._actions:
            if action.dest in _CALDERA_PIPELINE_GENERATED_INPUT_DESTINATIONS:
                action.help = argparse.SUPPRESS

    output = next(group for group in original_groups if group.title == "Output")
    inputs = parser.add_argument_group(
        "Input files",
        "Primary study data consumed directly by the selected execution mode.",
    )
    references = parser.add_argument_group(
        "Reference files",
        "External genomic resources and installed CALDERA model data.",
    )
    reference_settings = parser.add_argument_group(
        "CALDERA reference settings",
        "Declarations used to select and validate build-matched resources.",
    )
    pops_resource_settings = parser.add_argument_group("PoPS resource settings")
    software = parser.add_argument_group("CALDERA software")

    move_cli_help_actions(parser, output, _CALDERA_OUTPUT_DESTINATIONS)
    move_cli_help_actions(parser, inputs, _CALDERA_INPUT_FILE_DESTINATIONS)
    move_cli_help_actions(
        parser,
        references,
        _CALDERA_REFERENCE_FILE_DESTINATIONS,
    )
    move_cli_help_actions(
        parser,
        reference_settings,
        _CALDERA_REFERENCE_SETTING_DESTINATIONS,
    )
    move_cli_help_actions(
        parser,
        pops_resource_settings,
        _POPS_RESOURCE_SETTING_DESTINATIONS,
    )
    move_cli_help_actions(parser, software, _CALDERA_SOFTWARE_DESTINATIONS)
    external_software = next(
        (
            group
            for group in original_groups
            if group.title == "External software"
        ),
        None,
    )
    if external_software is not None:
        move_cli_help_actions(
            parser,
            external_software,
            _EXTERNAL_SOFTWARE_DESTINATIONS,
        )

    for group in parser._action_groups:
        if group.title == "PoPS covariates":
            group.title = "PoPS covariate settings"
        elif group.title == "PoPS feature selection":
            group.title = "PoPS feature-selection settings"
        elif group.title == "PoPS training":
            group.title = "PoPS training settings"

    leading_groups = [
        group
        for group in original_groups
        if group.title in {"positional arguments", "options"}
    ]
    ordered_groups = [
        group
        for title in _CALDERA_HELP_GROUP_ORDER
        for group in parser._action_groups
        if group.title == title
    ]
    remaining_groups = [
        group
        for group in parser._action_groups
        if group not in leading_groups and group not in ordered_groups
        and not (
            group.title in {
                "Input file",
                "MAGMA inputs",
                "PoPS input",
                "PoPS input compatibility",
                "GCTA-COJO inputs",
                "CALDERA inputs",
            }
            and not any(
                action.help is not argparse.SUPPRESS
                for action in group._group_actions
            )
        )
    ]
    parser._action_groups[:] = leading_groups + ordered_groups + remaining_groups


def get_caldera_pipeline_examples():
    """Return a runnable GRCh37 pipeline through PoPS, fine-mapping and CALDERA."""
    return (
        (
            "Run the complete GRCh37 PoPS, fine-mapping, and CALDERA workflow:",
            "postgwas pipeline",
            (
                "--modules caldera",
                "--vcf study_GRCh37.vcf.gz",
                "--genome-build GRCh37",
                "--ld-region-dir reference/ld_blocks",
                "--ld-block-populations EUR",
                "--ld-folder reference/ld",
                "--population EUR",
                "--magma-ld-reference reference/1000G_EUR",
                "--gene-location-file reference/NCBI37.3.gene.loc",
                "--feature-matrix-prefix reference/pops/features_munged/pops_features",
                "--feature-matrix-chunks 116",
                "--pops-gene-location-file reference/pops/GRCh37_gene_annot.tsv",
                "--pops-genome-build GRCh37",
                "--finemap-method susie",
                "--finemap-ld-reference reference/1000G_EUR",
                "--caldera-genome-build GRCh37",
                "--dataset-id STUDY",
                "--output-directory results",
            ),
        ),
    )


def get_caldera_parser(add_help=False, *, direct_controls=False):
    defaults = load_configuration()
    module = defaults.modules.caldera
    parser = argparse.ArgumentParser(add_help=add_help)
    inputs = parser.add_argument_group("CALDERA inputs")
    inputs.add_argument(
        "--pops-file", metavar="PATH", default=argparse.SUPPRESS,
        help="PoPS .preds file containing the gene-prioritisation scores.",
    )
    inputs.add_argument(
        "--credible-set-file", metavar="PATH", default=argparse.SUPPRESS,
        help="Annotated credible-set table consumed by direct CALDERA analysis.",
    )
    inputs.add_argument(
        "--caldera-repository", metavar="PATH", default=argparse.SUPPRESS,
        help=(
            "Optional CALDERA repository override. When omitted, PostGWAS uses "
            "the installed repository location configured in YAML."
        ),
    )
    inputs.add_argument(
        "--caldera-adapter-script", metavar="PATH", default=argparse.SUPPRESS,
        help=(
            "Optional PostGWAS CALDERA R adapter override. When omitted, the "
            "packaged adapter is used."
        ),
    )
    inputs.add_argument(
        "--caldera-genome-build", choices=list(defaults.resources.genomes),
        metavar="BUILD", default=argparse.SUPPRESS,
        help=help_with_default(
            "Genome build of credible-set coordinates and CALDERA gene locations",
            module.genome_build.value,
        ),
    )
    if direct_controls:
        controls = parser.add_argument_group("Configuration")
        controls.add_argument("--run-config", metavar="PATH", default=argparse.SUPPRESS)
    return parser


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="postgwas caldera",
        usage="postgwas caldera --pops-file PATH --credible-set-file PATH [options]",
        description="Prioritise causal genes with CALDERA from PoPS and credible sets.",
        formatter_class=AlignedRichHelpFormatter,
        epilog=format_cli_examples(
            (
                "Run CALDERA directly:", "postgwas caldera",
                (
                    "--pops-file pops/STUDY_pops.preds",
                    "--credible-set-file finemap/STUDY_credible_sets_GRCh37.tsv",
                    "--caldera-genome-build GRCh37",
                    "--dataset-id STUDY", "--output-directory results",
                ),
            ),
            notes=(
                "Direct credible-set tables require locus, chr, bp, and pip columns with at least 0.95 cumulative PIP per locus.",
            ),
        ),
        parents=[
            get_compute_parser(), get_common_out_parser(),
            get_caldera_parser(direct_controls=True),
        ],
    )
    for action in parser._actions:
        if action.dest in {"dataset_id", "output_directory", "threads", "memory_gb", "seed"}:
            action.default = argparse.SUPPRESS
            action.type = None
    mark_cli_required_help(parser, ("pops_file", "credible_set_file"))
    organize_caldera_help(parser)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    from postgwas.modules.caldera.service import run_caldera_direct

    try:
        run_caldera_direct(args)
    except MissingRequiredArgumentsError as exc:
        console = Console(stderr=True)
        print_screen_message("error", "CALDERA analysis failed. %s" % exc, console=console)
        parser.print_help(file=console.file)
        return 1
    except (CalderaError, ConfigurationError, OSError, ValueError) as exc:
        print_screen_message("error", "CALDERA analysis failed. %s" % exc, stderr=True)
        return 1
    return 0


__all__ = [
    "build_parser", "get_caldera_parser", "get_caldera_pipeline_examples", "main",
    "organize_caldera_help",
]
