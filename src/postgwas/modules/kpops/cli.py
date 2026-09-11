"""Command-line interface for K-POPS v1.0.0."""

from __future__ import annotations

import argparse
from typing import get_args

from rich.console import Console

from postgwas.cli.common import get_common_out_parser
from postgwas.cli.compute import get_compute_parser
from postgwas.config import load_configuration
from postgwas.config.models.modules.kpops import KPopsGeneUniversePolicy
from postgwas.core.errors import ConfigurationError, MissingRequiredArgumentsError
from postgwas.core.ui import (
    print_screen_message,
    AlignedRichHelpFormatter,
    format_cli_examples,
    help_with_default,
    mark_cli_required_help,
    move_cli_help_actions,
)
from postgwas.modules.kpops.errors import KPopsError


_KPOPS_OUTPUT_DESTINATIONS = (
    "dataset_id",
    "output_directory",
    "save_attribution_files",
)
_KPOPS_INPUT_FILE_DESTINATIONS = ("vcf", "magma_association_prefix")
_KPOPS_REFERENCE_FILE_DESTINATIONS = (
    "magma_ld_reference",
    "gene_location_file",
    "gene_set_file",
    "kpops_gene_annotation_file",
    "kernel_matrix_prefix",
)
_KPOPS_REFERENCE_SETTING_DESTINATIONS = (
    "kpops_genome_build",
    "gene_universe_policy",
    "kpops_gene_universe_policy",
)
_KPOPS_SOFTWARE_DESTINATIONS = ("kpops_script",)
_KPOPS_HELP_GROUP_ORDER = (
    "Output",
    "Input files",
    "Reference files",
    "K-POPS reference settings",
    "Variant identifiers",
    "MAGMA analysis settings",
    "MAGMA analysis scope",
    "K-POPS training and explanation settings",
    "K-POPS software",
    "Performance",
    "Screen output",
    "Run continuation",
    "Configuration",
    "Choose analyses",
)


def organize_kpops_help(parser: argparse.ArgumentParser) -> None:
    """Place outputs, analysis inputs, and reference resources before settings."""
    original_groups = list(parser._action_groups)
    output = next(group for group in original_groups if group.title == "Output")
    inputs = parser.add_argument_group(
        "Input files",
        "Primary study data consumed directly by the selected execution mode.",
    )
    references = parser.add_argument_group(
        "Reference files",
        "MAGMA and K-POPS resources defining LD, genes, and the kernel.",
    )
    reference_settings = parser.add_argument_group(
        "K-POPS reference settings",
        "Declarations controlling compatibility across K-POPS resources.",
    )
    software = parser.add_argument_group("K-POPS software")

    move_cli_help_actions(parser, output, _KPOPS_OUTPUT_DESTINATIONS)
    move_cli_help_actions(parser, inputs, _KPOPS_INPUT_FILE_DESTINATIONS)
    move_cli_help_actions(
        parser,
        references,
        _KPOPS_REFERENCE_FILE_DESTINATIONS,
    )
    move_cli_help_actions(
        parser,
        reference_settings,
        _KPOPS_REFERENCE_SETTING_DESTINATIONS,
    )
    move_cli_help_actions(parser, software, _KPOPS_SOFTWARE_DESTINATIONS)

    leading_groups = [
        group
        for group in original_groups
        if group.title in {"positional arguments", "options"}
    ]
    ordered_groups = [
        group
        for title in _KPOPS_HELP_GROUP_ORDER
        for group in parser._action_groups
        if group.title == title
    ]
    remaining_groups = [
        group
        for group in parser._action_groups
        if group not in leading_groups and group not in ordered_groups
        and not (
            group.title in {"Input file", "MAGMA inputs", "K-POPS inputs"}
            and not any(
                action.help is not argparse.SUPPRESS
                for action in group._group_actions
            )
        )
    ]
    parser._action_groups[:] = leading_groups + ordered_groups + remaining_groups


def get_kpops_pipeline_examples():
    """Return a runnable formatter-to-MAGMA-to-K-POPS pipeline example."""
    return (
        (
            "Run MAGMA gene analysis followed by K-POPS prioritisation:",
            "postgwas pipeline",
            (
                "--modules kpops",
                "--vcf study_GRCh37.vcf.gz",
                (
                    "--magma-ld-reference /path/to/postgwas-resources/magma/"
                    "functional_mapping/base/ld_reference/g1000_eur/g1000_eur"
                ),
                (
                    "--gene-location-file reference/PoPS_GRCh37_strand_aware.loc"
                ),
                (
                    "--kpops-gene-annotation-file /path/to/postgwas-resources/"
                    "pops/GRCh37_gene_annot_jun10.txt"
                ),
                (
                    "--kernel-matrix-prefix /path/to/postgwas-resources/kpops/"
                    "kernels/GRCh37/pops_features_standardized_linear"
                ),
                "--kpops-genome-build GRCh37",
                "--dataset-id STUDY",
                "--output-directory results",
            ),
        ),
    )


def get_kpops_direct_examples():
    """Return direct-run examples, including the upstream ApoB model settings."""
    return (
        (
            "Run K-POPS directly:",
            "postgwas kpops",
            (
                (
                    "--magma-association-prefix magma/02_intermediates/"
                    "positional/native_outputs/STUDY_magma_35up_10down"
                ),
                (
                    "--kpops-gene-annotation-file /path/to/postgwas-resources/"
                    "pops/GRCh37_gene_annot_jun10.txt"
                ),
                (
                    "--kernel-matrix-prefix /path/to/postgwas-resources/kpops/"
                    "kernels/GRCh37/pops_features_standardized_linear"
                ),
                "--kpops-genome-build GRCh37",
                "--dataset-id STUDY",
                "--output-directory results",
            ),
        ),
        (
            "Reproduce the upstream v1.0.0 ApoB command settings:",
            "postgwas kpops",
            (
                "--magma-association-prefix /path/to/ApoB",
                (
                    "--kpops-gene-annotation-file /path/to/postgwas-resources/"
                    "pops/GRCh37_gene_annot_jun10.txt"
                ),
                "--kernel-matrix-prefix /path/to/ApoB/kernel_linear",
                "--kpops-genome-build GRCh37",
                "--gene-universe-policy strict",
                "--use-magma-covariates",
                "--training-chromosomes loco",
                "--kpops-device cuda",
                "--top-contributor-gene-count 5",
                (
                    "--anchor-genes PCSK9 ANGPTL3 APOB ABCG5 ALB NPC1L1 GIGYF1 "
                    "JAK2 A1CF PDE3B APOC3 CETP ASGR1 LDLR ZNF234 ZNF229 "
                    "NECTIN2 RRBP1"
                ),
                "--anchor-gene-type NAME",
                "--seed 42",
                "--kpops-verbose",
                "--dataset-id ApoB",
                "--output-directory results/ApoB",
            ),
        ),
    )


def get_kpops_pipeline_parser():
    """Keep K-POPS model overrides independent of PoPS in joint pipelines."""
    return get_kpops_parser(pipeline=True)


def get_kpops_parser(add_help=False, *, direct_controls=False, pipeline=False):
    defaults = load_configuration()
    module = defaults.modules.kpops
    schema = module.input_schema
    parser = argparse.ArgumentParser(add_help=add_help)
    inputs = parser.add_argument_group("K-POPS inputs")
    inputs.add_argument(
        "--magma-association-prefix",
        metavar="PREFIX",
        default=argparse.SUPPRESS,
        help=(
            "MAGMA output prefix before the %s and %s suffixes. Both files must "
            "come from the same gene-association run and have identical gene order."
            % (schema.magma_genes_out_suffix, schema.magma_genes_raw_suffix)
        ),
    )
    inputs.add_argument(
        "--kpops-gene-annotation-file",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "Gene annotation for the declared build. It must contain the configured "
            "%s, %s, %s, and %s columns; gene IDs must be unique."
            % (
                schema.gene_id_column,
                schema.gene_name_column,
                schema.chromosome_column,
                schema.tss_column,
            )
        ),
    )
    inputs.add_argument(
        "--kernel-matrix-prefix",
        metavar="PREFIX",
        default=argparse.SUPPRESS,
        help=(
            "Kernel prefix before the %s matrix and %s gene-order suffixes. The "
            "float32 matrix must be square, finite, symmetric, and ordered exactly "
            "like the gene file."
            % (schema.kernel_matrix_suffix, schema.kernel_genes_suffix)
        ),
    )
    inputs.add_argument(
        "--kpops-gene-universe-policy" if pipeline else "--gene-universe-policy",
        choices=get_args(KPopsGeneUniversePolicy),
        metavar="POLICY",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "How MAGMA genes absent from the K-POPS annotation or kernel are "
            "handled: intersect retains the validated shared genes and writes "
            "an audit; strict stops before fitting",
            module.gene_universe_policy,
        ),
    )
    inputs.add_argument(
        "--kpops-script",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Executable name or path for the pinned K-POPS program. Installation "
            "normally provides it; override only for a nonstandard installation",
            module.script_path,
        ),
    )
    inputs.add_argument(
        "--kpops-genome-build", choices=list(defaults.resources.genomes),
        metavar="BUILD", default=argparse.SUPPRESS,
        help=help_with_default(
            "Genome build shared by MAGMA coordinates, the gene annotation, and "
            "kernel genes. PostGWAS validates this declaration but does not infer "
            "or lift coordinates",
            module.genome_build,
        ),
    )
    settings = parser.add_argument_group("K-POPS training and explanation settings")
    settings.add_argument(
        "--kpops-training-chromosomes" if pipeline else "--training-chromosomes",
        nargs="+", metavar="CHROM",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Training design: 'loco' fits one leave-one-chromosome-out model per "
            "chromosome; 'all' trains and predicts with all chromosomes; explicit "
            "annotation chromosome labels train once on those chromosomes and "
            "predict the remaining chromosomes",
            " ".join(module.training_chromosomes),
        ),
    )
    settings.add_argument(
        "--kpops-device", choices=("cpu", "cuda", "mps"),
        default=argparse.SUPPRESS,
        help=help_with_default(
            "PyTorch compute device: cpu for general use, cuda for a supported "
            "NVIDIA GPU, or mps for supported Apple Silicon",
            module.device,
        ),
    )
    settings.add_argument(
        "--top-contributor-gene-count", type=int, metavar="N",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Number of highest positive kernel-contribution genes reported for "
            "each prediction; this changes explanation columns, not the K-POPS score",
            module.top_contributor_gene_count,
        ),
    )
    settings.add_argument(
        "--anchor-genes",
        nargs="+",
        metavar="GENE",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Known genes used to calculate each prediction's Anchor_Score. Supply "
            "space-separated identifiers in the format selected by "
            "--anchor-gene-type",
            "none" if not module.anchor_genes else " ".join(module.anchor_genes),
        ),
    )
    settings.add_argument(
        "--anchor-gene-type", choices=("ENSGID", "NAME"),
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Identifier format used by --anchor-genes: ENSGID for Ensembl gene IDs "
            "or NAME for gene symbols. Ambiguous annotation symbols are rejected",
            module.anchor_gene_type,
        ),
    )
    settings.add_argument(
        "--kpops-use-magma-covariates" if pipeline else "--use-magma-covariates",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Project MAGMA gene-size, gene-density, and inverse-MAC covariates "
            "derived from the .genes.raw file out of target Z statistics before "
            "kernel fitting; use %s to fit unadjusted targets"
            % ("--no-kpops-use-magma-covariates" if pipeline else "--no-use-magma-covariates"),
            module.use_magma_covariates,
        ),
    )
    settings.add_argument(
        "--save-attribution-files", action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Write the dense gene-by-gene contribution matrix and its row and "
            "column gene-order files. Storage grows quadratically with gene count",
            module.save_attribution_files,
        ),
    )
    settings.add_argument(
        "--kpops-verbose", action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Include verbose messages from the upstream K-POPS program in the "
            "captured command log",
            module.verbose,
        ),
    )
    if direct_controls:
        controls = parser.add_argument_group("Configuration")
        controls.add_argument(
            "--run-config",
            metavar="PATH",
            default=argparse.SUPPRESS,
            help=(
                "YAML file containing K-POPS and shared run settings. Explicit "
                "command-line values override matching YAML keys."
            ),
        )
    return parser


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="postgwas kpops",
        usage="postgwas kpops --magma-association-prefix PREFIX [options]",
        description="Prioritise genes with K-POPS from validated MAGMA and kernel inputs.",
        formatter_class=AlignedRichHelpFormatter,
        epilog=format_cli_examples(
            *get_kpops_direct_examples(),
            (
                "Export reusable settings:", "postgwas config export",
                ("--module kpops", "--style full", "--output kpops.yaml"),
            ),
            notes=(
                "Direct mode requires PREFIX.genes.out and PREFIX.genes.raw from the same MAGMA run.",
                (
                    "The ApoB example reproduces upstream command settings; the "
                    "upstream repository does not bundle ApoB MAGMA outputs or "
                    "its prepared test kernel, so supply those two prefixes."
                ),
                "The ApoB CUDA setting requires a supported NVIDIA GPU; use --kpops-device cpu when CUDA is unavailable.",
            ),
        ),
        parents=[
            get_compute_parser(), get_common_out_parser(),
            get_kpops_parser(direct_controls=True),
        ],
    )
    for action in parser._actions:
        if action.dest in {"dataset_id", "output_directory", "threads", "memory_gb", "seed"}:
            action.default = argparse.SUPPRESS
            action.type = None
    mark_cli_required_help(
        parser,
        (
            "magma_association_prefix",
            "kpops_gene_annotation_file",
            "kernel_matrix_prefix",
            "kpops_genome_build",
            "dataset_id",
            "output_directory",
        ),
    )
    organize_kpops_help(parser)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    from postgwas.modules.kpops.service import run_kpops_direct

    try:
        run_kpops_direct(args)
    except MissingRequiredArgumentsError as exc:
        console = Console(stderr=True)
        print_screen_message("error", "K-POPS analysis failed. %s" % exc, console=console)
        parser.print_help(file=console.file)
        return 1
    except (KPopsError, ConfigurationError, OSError, ValueError) as exc:
        print_screen_message("error", "K-POPS analysis failed. %s" % exc, stderr=True)
        return 1
    return 0


__all__ = [
    "build_parser",
    "get_kpops_direct_examples",
    "get_kpops_parser",
    "get_kpops_pipeline_parser",
    "get_kpops_pipeline_examples",
    "main",
    "organize_kpops_help",
]
