"""Command-line interface for official single-trait MiXeR analyses."""

from __future__ import annotations

import argparse
from typing import get_args

from rich.console import Console

from postgwas.config import load_configuration
from postgwas.cli.common import get_common_out_parser, get_genome_build_parser
from postgwas.cli.compute import get_compute_parser
from postgwas.core.errors import ConfigurationError
from postgwas.core.execution.runtime import validate_path
from postgwas.core.ui import (
    AlignedRichHelpFormatter,
    format_cli_examples,
    help_with_conditional_requirement,
    help_with_default,
)
from postgwas.config.models.modules.mixer import MixerAnalysis, MixerExecutionBackend
from postgwas.modules.mixer.service import MixerError, run_mixer_direct


def get_mixer_parser(add_help=False, *, direct_controls=False):
    defaults = load_configuration()
    mixer_defaults = defaults.modules.mixer
    chromosome_placeholder = mixer_defaults.workflow.chromosome_placeholder
    mixer_input_pattern = defaults.modules.formatting.exports["mixer"].output_file
    container_defaults = defaults.resources.containers.mixer
    parser = argparse.ArgumentParser(add_help=add_help)
    inputs = parser.add_argument_group("MiXeR inputs")
    inputs.add_argument(
        "--mixer-input-file",
        type=validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True),
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "MiXeR input created by postgwas formatter. Configured filename "
            "pattern: %s." % mixer_input_pattern
        ),
    )
    inputs.add_argument(
        "--bim-file-pattern",
        metavar="PATTERN",
        default=argparse.SUPPRESS,
        help=help_with_conditional_requirement(
            "Per-chromosome BIM path containing %s, for example reference/chr%s.bim."
            % (chromosome_placeholder, chromosome_placeholder),
            "for every MiXeR analysis",
        ),
    )
    inputs.add_argument(
        "--ld-file-pattern",
        metavar="PATTERN",
        default=argparse.SUPPRESS,
        help=(
            "Per-chromosome MiXeR LD path containing %s, for example "
            "reference/chr%s.ld."
            % (chromosome_placeholder, chromosome_placeholder)
        ),
    )
    gsa = parser.add_argument_group("GSA-MiXeR reference inputs")
    gsa.add_argument(
        "--gsa-annotation-file-pattern",
        metavar="PATTERN",
        default=argparse.SUPPRESS,
        help=help_with_conditional_requirement(
            "Per-chromosome SNP annotation path containing %s"
            % chromosome_placeholder,
            "for gsa/all",
        ),
    )
    gsa.add_argument(
        "--gsa-loadlib-file-pattern",
        metavar="PATTERN",
        default=argparse.SUPPRESS,
        help=(
            "Optional precomputed GSA load-library path containing %s. When "
            "omitted, GSA-MiXeR uses --ld-file-pattern."
            % chromosome_placeholder
        ),
    )
    for option, label in (
        ("baseline", "baseline annotation"),
        ("model", "full-model gene annotation"),
        ("test", "gene-set test annotation"),
    ):
        gsa.add_argument(
            "--gsa-%s-go-file" % option,
            type=validate_path(
                must_exist=True, must_be_file=True, must_not_be_empty=True,
            ),
            metavar="PATH",
            default=argparse.SUPPRESS,
            help=help_with_conditional_requirement(
                "GSA-MiXeR %s table" % label,
                "for gsa/all",
            ),
        )
    settings = parser.add_argument_group("MiXeR execution")
    settings.add_argument(
        "--analysis",
        choices=get_args(MixerAnalysis),
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Single-trait analysis to run: architecture only (univariate), "
            "gene-set enrichment only (gsa), or both (all)",
            mixer_defaults.analysis,
        ),
    )
    settings.add_argument(
        "--mixer-backend",
        choices=get_args(MixerExecutionBackend),
        default=argparse.SUPPRESS,
        help=help_with_default(
            "How MiXeR is started. 'auto' uses a local installation when available, "
            "otherwise the configured GSA-MiXeR container",
            mixer_defaults.execution_backend,
        ),
    )
    settings.add_argument(
        "--mixer", metavar="PATH", default=argparse.SUPPRESS,
        help=help_with_default(
            "Path to the official mixer.py script",
            defaults.resources.executables.mixer,
        ),
    )
    settings.add_argument(
        "--mixer-figures",
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Path to mixer_figures.py, used for official univariate QQ and "
            "power diagnostics",
            defaults.resources.executables.mixer_figures,
        ),
    )
    settings.add_argument(
        "--mixer-container-image",
        metavar="IMAGE",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "GSA-MiXeR container image",
            container_defaults.image,
        ),
    )
    settings.add_argument(
        "--mixer-container-runtime",
        metavar="COMMAND",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Container command",
            container_defaults.runtime,
        ),
    )
    settings.add_argument(
        "--mixer-container-platform",
        metavar="PLATFORM",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Container platform",
            container_defaults.platform,
        ),
    )
    if direct_controls:
        settings.add_argument(
            "--run-config", metavar="PATH", default=argparse.SUPPRESS,
            help="YAML settings for MiXeR, references, software paths, and execution.",
        )
        settings.add_argument(
            "--dry-run", action="store_true", default=argparse.SUPPRESS,
            help="Validate argument structure and log commands without executing MiXeR.",
        )
    return parser


def build_parser():
    defaults = load_configuration()
    mixer = defaults.modules.mixer
    formatter_pattern = defaults.modules.formatting.exports["mixer"].output_file
    example_input = formatter_pattern.format(dataset_id="STUDY")
    chromosome_placeholder = mixer.workflow.chromosome_placeholder
    parser = argparse.ArgumentParser(
        prog="postgwas mixer",
        usage=(
            "postgwas mixer --mixer-input-file PATH\n"
            "                      [--analysis {univariate,gsa,all}] [options]"
        ),
        description=(
            "Run official single-trait MiXeR architecture analysis, GSA-MiXeR "
            "gene-set enrichment, or both.\n"
            "Both analyses use one PostGWAS formatter input. "
            "No multi-trait analysis is exposed."
        ),
        epilog=format_cli_examples(
            (
                "Create the MiXeR input first:",
                "postgwas formatter",
                (
                    "--vcf study.vcf.gz",
                    "--dataset-id STUDY",
                    "--output-directory formatted",
                    "--format mixer",
                ),
            ),
            (
                "Run single-trait architecture analysis:",
                "postgwas mixer",
                (
                    "--analysis univariate",
                    "--mixer-input-file formatted/%s" % example_input,
                    "--dataset-id STUDY",
                    "--output-directory results",
                    "--bim-file-pattern 'reference/chr%s.bim'" % chromosome_placeholder,
                    "--ld-file-pattern 'reference/chr%s.ld'" % chromosome_placeholder,
                ),
            ),
            (
                "Run single-trait GSA-MiXeR:",
                "postgwas mixer",
                (
                    "--analysis gsa",
                    "--mixer-input-file formatted/%s" % example_input,
                    "--dataset-id STUDY",
                    "--output-directory results",
                    "--bim-file-pattern 'reference/chr%s.bim'" % chromosome_placeholder,
                    "--ld-file-pattern 'reference/chr%s.ld'" % chromosome_placeholder,
                    "--gsa-annotation-file-pattern 'reference/chr%s.annot.gz'"
                    % chromosome_placeholder,
                    "--gsa-baseline-go-file baseline.tsv",
                    "--gsa-model-go-file genes.tsv",
                    "--gsa-test-go-file gene_sets.tsv",
                ),
            ),
            (
                "Export a complete reloadable MiXeR configuration:",
                "postgwas config export",
                ("--module mixer", "--style full", "--output mixer.yaml"),
            ),
            notes=(
                "Set analysis: all in mixer.yaml to run normal MiXeR and GSA-MiXeR together.",
            ),
        ),
        formatter_class=AlignedRichHelpFormatter,
        parents=[
            get_compute_parser(),
            get_genome_build_parser(
                available_builds=list(defaults.resources.genomes),
                default_build=str(mixer.genome_build),
                suppress_default=True,
            ),
            get_common_out_parser(),
            get_mixer_parser(direct_controls=True),
        ],
    )
    for action in parser._actions:
        if action.dest == "mixer_input_file":
            action.required = True
        elif action.dest in {"genome_build", "dataset_id", "output_directory"}:
            # Omission must preserve a run-config value instead of silently
            # replacing it with a parser-owned value.
            action.default = argparse.SUPPRESS
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        run_mixer_direct(args)
    except (MixerError, ConfigurationError, OSError, ValueError) as exc:
        Console(stderr=True).print("\n[bold red]MiXeR failed.[/bold red] %s\n" % exc)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = ["build_parser", "get_mixer_parser", "main"]
