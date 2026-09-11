#!/usr/bin/env python3
"""
PostGWAS Fine-Mapping CLI
Supports two engines:
  • SuSiE  (Python/R-based)
  • FINEMAP (C++ binary)
Works in:
  • DIRECT mode   → run fine-mapping on existing locus + sumstats
  • PIPELINE mode → Harmonisation → LD Blocks → QC → Formatter → Fine-mapping
"""

import argparse
import sys
from typing import get_args

from postgwas.cli.compute import get_compute_parser
from postgwas.config import load_configuration
from postgwas.config.models.modules.fine_mapping import FineMappingEngine
from postgwas.config.models.modules.ld_clumping import LDClumpingMethod
from postgwas.core.ui import (
    AlignedRichHelpFormatter,
    cli_option_value,
    cli_option_values,
    format_cli_examples,
    format_cli_help_block,
    help_with_conditional_requirement,
    help_with_default,
)


# =========================================================
# BACKEND RUNNERS
# =========================================================
from postgwas.modules.fine_mapping.service import run_fine_mapping

# Shared parsers
from postgwas.cli.common import (
    get_common_out_parser,
    get_plink_binary_parser,
)

from postgwas.modules.fine_mapping.arguments import (
    get_finemap_common_parser,
    get_common_susie_arguments,
    get_common_finemap_finemap_arguments,
    get_finemap_susie_inputs_parser,
    get_finemap_finemap_arguments,
)


_FINEMAP_ENGINE_GROUPS = frozenset({"SuSiE settings", "FINEMAP settings"})
_FINEMAP_UNGROUPED_DESTINATIONS = frozenset({"sss", "cond"})
_COJO_GROUPS = frozenset({
    "GCTA-COJO inputs",
    "Input and reference compatibility",
    "GCTA-COJO stepwise selection",
    "GCTA-COJO physical locus definition",
})
_SELECTION_VISIBLE_DESTINATIONS = frozenset({
    "help",
    "threads",
    "memory_gb",
    "seed",
    "show_screen",
    "resume",
    "overwrite",
    "clumping_methods",
    "finemap_method",
    "run_config",
    "modules",
})
_FINEMAP_HELP_GROUP_ORDER = (
    "Input file",
    "Output",
    "LD-block annotation",
    "Genome coordinates",
    "Reference population",
    "External software",
    "LD-clumping analysis",
    "Standard r² clumping",
    "GCTA-COJO inputs",
    "Input and reference compatibility",
    "GCTA-COJO stepwise selection",
    "GCTA-COJO physical locus definition",
    "Variant identifiers",
    "Fine-mapping settings",
    "SuSiE settings",
    "FINEMAP settings",
    "Performance",
    "Screen output",
    "Run continuation",
    "Configuration",
    "Choose analyses",
)


def _engine_label(engine: str) -> str:
    return "SuSiE-RSS" if engine == "susie" else "FINEMAP"


def _clumping_label(methods) -> str:
    labels = {
        "region": "annotated-region",
        "standard": "standard r²",
        "cojo-slct": "GCTA-COJO",
    }
    return " + ".join(labels[method] for method in methods)


def _clumping_step_label(methods) -> str:
    if tuple(methods) == ("standard",):
        return "Standard LD clumping (r²)"
    return "LD clumping (%s)" % _clumping_label(methods)


def _suppress_groups(parser: argparse.ArgumentParser, titles) -> None:
    for group in parser._action_groups:
        if group.title in titles:
            for action in group._group_actions:
                action.required = False
                action.help = argparse.SUPPRESS


def _order_help_groups(parser: argparse.ArgumentParser, *, selection_only=False):
    original = list(parser._action_groups)
    leading = [
        group
        for group in original
        if group.title in {"positional arguments", "options"}
    ]
    if selection_only:
        preferred = (
            "LD-clumping analysis",
            "Fine-mapping settings",
            "Performance",
            "Screen output",
            "Run continuation",
            "Configuration",
            "Choose analyses",
        )
    else:
        preferred = _FINEMAP_HELP_GROUP_ORDER
    ordered = [
        group
        for title in preferred
        for group in parser._action_groups
        if group.title == title
    ]
    remaining = [
        group
        for group in parser._action_groups
        if group not in leading and group not in ordered
        and any(
            action.help is not argparse.SUPPRESS
            for action in group._group_actions
        )
    ]
    parser._action_groups[:] = leading + ordered + remaining


def _finemap_selection_examples():
    return (
        (
            "Show every option for standard clumping with SuSiE-RSS:",
            "postgwas pipeline",
            (
                "--modules finemap",
                "--clumping-methods standard",
                "--finemap-method susie",
                "--help",
            ),
        ),
        (
            "Show every option for standard clumping with FINEMAP:",
            "postgwas pipeline",
            (
                "--modules finemap",
                "--clumping-methods standard",
                "--finemap-method finemap",
                "--help",
            ),
        ),
    )


def get_finemap_pipeline_examples(
    *, clumping_methods=None, engine=None,
):
    """Return one dependency-complete example for the selected workflow."""
    defaults = load_configuration()
    methods = tuple(clumping_methods or defaults.modules.ld_clumping.methods)
    selected_engine = engine or defaults.modules.fine_mapping.engine
    arguments = [
        "--modules finemap",
        "--clumping-methods %s" % " ".join(methods),
        "--finemap-method %s" % selected_engine,
        "--vcf study_GRCh37.vcf.gz",
        "--genome-build GRCh37",
        "--population EUR",
    ]
    if "region" in methods:
        arguments.extend((
            "--ld-region-dir reference/ld_blocks",
            "--ld-block-populations EUR",
        ))
    if "standard" in methods:
        arguments.append("--ld-folder reference/pairwise_ld")
    if "cojo-slct" in methods:
        arguments.append("--cojo-reference-prefix reference/1000G_EUR")
    arguments.extend((
        "--variant-id-type unique",
        "--finemap-ld-reference reference/1000G_EUR",
        "--plink %s" % ("plink" if selected_engine == "susie" else "plink2"),
        "--dataset-id STUDY",
        "--output-directory results",
    ))
    return (
        (
            "Run %s clumping with %s fine-mapping:"
            % (_clumping_label(methods), _engine_label(selected_engine)),
            "postgwas pipeline",
            tuple(arguments),
        ),
    )


def organize_finemap_pipeline_help(parser: argparse.ArgumentParser) -> None:
    """Route fine-mapping pipeline help through explicit method selection."""
    arguments = tuple(getattr(parser, "_postgwas_help_arguments", ()))
    run_config = cli_option_value(arguments, "--run-config")
    configuration = load_configuration(run_config)
    configured_clumping = tuple(configuration.modules.ld_clumping.methods)
    configured_engine = configuration.modules.fine_mapping.engine
    explicit_clumping = cli_option_values(arguments, "--clumping-methods")
    explicit_engine = cli_option_value(arguments, "--finemap-method")
    selected_clumping = explicit_clumping or configured_clumping
    selected_engine = explicit_engine or configured_engine

    actions = {action.dest: action for action in parser._actions}
    clumping_action = actions["clumping_methods"]
    engine_action = actions["finemap_method"]
    clumping_description = (
        "Choose one or more analyses. standard defines the loci passed to "
        "fine-mapping; region and cojo-slct add parallel locus summaries"
    )
    engine_description = "Choose the Bayesian fine-mapping engine"
    clumping_action.help = help_with_default(
        clumping_description, " ".join(configured_clumping),
    )
    engine_action.help = help_with_default(
        engine_description, configured_engine,
    )
    plink_action = actions.get("plink")
    if plink_action is not None:
        plink_action.help = help_with_default(
            "Full path or command name for the PLINK executable used for LD",
            (
                configuration.resources.executables.plink
                if selected_engine == "susie"
                else configuration.resources.executables.plink2
            ),
        )

    valid_clumping = set(get_args(LDClumpingMethod))
    valid_engines = set(get_args(FineMappingEngine))
    selection_error = None
    unknown_clumping = sorted(set(explicit_clumping) - valid_clumping)
    if unknown_clumping:
        selection_error = "Unknown LD-clumping method(s): %s." % ", ".join(
            unknown_clumping
        )
    elif explicit_engine is not None and explicit_engine not in valid_engines:
        selection_error = "Unknown fine-mapping method: %s." % explicit_engine

    selection_complete = bool(run_config) or (
        bool(explicit_clumping) and explicit_engine is not None
    )
    if selection_complete and "standard" not in selected_clumping:
        selection_error = (
            "Fine-mapping requires standard LD clumping because the standard "
            "genomic-risk-locus table is the validated locus source. Include "
            "standard in --clumping-methods or in "
            "modules.ld_clumping.methods."
        )

    if not selection_complete or selection_error is not None:
        clumping_action.help = help_with_conditional_requirement(
            clumping_action.help,
            "unless --run-config supplies modules.ld_clumping.methods",
        )
        engine_action.help = help_with_conditional_requirement(
            engine_action.help,
            "unless --run-config supplies modules.fine_mapping.engine",
        )
        for action in parser._actions:
            if action.dest not in _SELECTION_VISIBLE_DESTINATIONS:
                action.required = False
                action.help = argparse.SUPPRESS
        parser.usage = (
            "postgwas pipeline --modules finemap "
            "--clumping-methods METHOD [METHOD ...] "
            "--finemap-method METHOD [--help]"
        )
        parser.description = format_cli_help_block(
            "Fine-mapping supports multiple workflow combinations.\n"
            "Select the LD-clumping analyses and fine-mapping engine before "
            "viewing the complete workflow options or starting execution."
        )
        parser._postgwas_pipeline_examples = _finemap_selection_examples()
        parser._postgwas_pipeline_selection_only = True
        parser._postgwas_pipeline_selection_error = selection_error
        parser._postgwas_pipeline_selection_visible_destinations = (
            _SELECTION_VISIBLE_DESTINATIONS
        )
        _order_help_groups(parser, selection_only=True)
        return

    if "cojo-slct" not in selected_clumping:
        _suppress_groups(parser, _COJO_GROUPS)
    _suppress_groups(
        parser,
        _FINEMAP_ENGINE_GROUPS - {
            "SuSiE settings" if selected_engine == "susie" else "FINEMAP settings"
        },
    )
    clumping_action.help = argparse.SUPPRESS
    engine_action.help = argparse.SUPPRESS
    if selected_engine == "susie":
        for destination in _FINEMAP_UNGROUPED_DESTINATIONS:
            action = actions.get(destination)
            if action is not None:
                action.required = False
                action.help = argparse.SUPPRESS
    parser.usage = (
        "postgwas pipeline [options] --modules finemap "
        "--clumping-methods %s --finemap-method %s"
        % (" ".join(selected_clumping), selected_engine)
    )
    parser.description = format_cli_help_block(
        "Selected fine-mapping workflow\n"
        "  LD-clumping methods : %s\n"
        "  Fine-mapping method : %s"
        % (_clumping_label(selected_clumping), _engine_label(selected_engine))
    )
    parser._postgwas_pipeline_examples = get_finemap_pipeline_examples(
        clumping_methods=selected_clumping,
        engine=selected_engine,
    )
    parser._postgwas_pipeline_step_names = {
        "formatter": "%s input preparation" % _engine_label(selected_engine),
        "ld_clump": _clumping_step_label(selected_clumping),
        "finemap": "%s fine-mapping" % _engine_label(selected_engine),
    }
    parser._postgwas_pipeline_step_descriptions = {
        "formatter": (
            "Validate the GWAS-VCF and create the %s summary-statistics table."
            % _engine_label(selected_engine)
        ),
        "ld_clump": (
            "Define the loci passed to fine-mapping with standard r² clumping."
        ),
        "finemap": "Run %s and validate its credible sets."
        % _engine_label(selected_engine),
    }
    _order_help_groups(parser)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="postgwas finemap",
        usage="postgwas finemap --finemap-method {susie,finemap} [options]",
        description=(
            "Fine-map selected loci with SuSiE-RSS or FINEMAP.\n"
            "Standalone mode uses formatter-created summary statistics and an existing locus file."
        ),
        epilog=format_cli_examples(
            (
                "Run SuSiE-RSS fine-mapping:",
                "postgwas finemap",
                (
                    "--finemap-method susie",
                    "--susie-input-file formatted/STUDY_susie.tsv.gz",
                    "--locus-file loci.tsv",
                    "--locus-type range",
                    "--window-kb 0",
                    "--finemap-ld-reference reference/1000G_EUR",
                    "--plink plink",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Run FINEMAP:",
                "postgwas finemap",
                (
                    "--finemap-method finemap",
                    "--finemap-in-files formatted/STUDY_finemap.tsv.gz",
                    "--locus-file loci.tsv",
                    "--locus-type range",
                    "--window-kb 0",
                    "--finemap-ld-reference reference/1000G_EUR",
                    "--plink plink2",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Export the fine-mapping pipeline configuration:",
                "postgwas config export",
                ("--pipeline finemap", "--style full", "--output finemap_pipeline.yaml"),
            ),
        ),
        formatter_class=AlignedRichHelpFormatter,
        parents=[
            get_compute_parser(),
            get_common_out_parser(),
            get_finemap_common_parser(include_genome_build=True),
            get_finemap_susie_inputs_parser(),
            get_common_susie_arguments(),
            get_common_finemap_finemap_arguments(),
            get_finemap_finemap_arguments(),
            get_plink_binary_parser(),
        ],
    )
    configuration = parser.add_argument_group("Configuration")
    configuration.add_argument(
        "--run-config",
        default=argparse.SUPPRESS,
        metavar="PATH",
        help=(
            "YAML file containing fine-mapping settings. Explicit command-line "
            "values override matching YAML values."
        ),
    )
    return parser


def main():
    parser = build_parser()

    # If no arguments provided → show help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    run_fine_mapping(args)


if __name__ == "__main__":
    main()


__all__ = [
    "build_parser",
    "get_finemap_pipeline_examples",
    "main",
    "organize_finemap_pipeline_help",
]
