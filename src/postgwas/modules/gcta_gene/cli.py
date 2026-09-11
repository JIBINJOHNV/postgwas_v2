"""Command-line interface for GCTA fastBAT and mBAT-combo gene tests."""

from __future__ import annotations

import argparse
from typing import get_args

from rich.console import Console

from postgwas.cli.common import get_common_out_parser
from postgwas.cli.compute import get_compute_parser
from postgwas.config import load_configuration
from postgwas.config.models.modules.gcta_gene import (
    GctaChromosomeLabelPolicy,
    GctaDuplicateGenePolicy,
    GctaEmptyPathwayPolicy,
    GctaGeneMethod,
    GctaOversizedSetPolicy,
    GctaUnmappedGenePolicy,
)
from postgwas.config.models.common import MHC_POLICIES
from postgwas.core.errors import ConfigurationError, MissingRequiredArgumentsError
from postgwas.core.execution.runtime import validate_path
from postgwas.core.ui import (
    print_screen_message,
    AlignedRichHelpFormatter,
    cli_option_value,
    format_cli_examples,
    help_with_conditional_requirement,
    help_with_default,
    mark_cli_required_help,
    move_cli_help_actions,
)
from postgwas.modules.gcta_gene.errors import GctaGeneError


_GCTA_OUTPUT_DESTINATIONS = ("dataset_id", "output_directory")
_GCTA_INPUT_FILE_DESTINATIONS = (
    "vcf",
    "gcta_input_file",
    "fastbat_set_list",
    "gmt",
)
_GCTA_REFERENCE_FILE_DESTINATIONS = (
    "gcta_reference_prefix",
    "gene_list",
)
_GCTA_REFERENCE_SETTING_DESTINATIONS = (
    "genome_build",
    "gcta_reference_population",
)
_GCTA_HELP_GROUP_ORDER = (
    "Output",
    "Input files",
    "Reference files",
    "GCTA reference settings",
    "Variant identifiers",
    "GCTA analysis settings",
    "GCTA analysis scope",
    "GCTA variant matching",
    "GCTA pathway conversion",
    "GCTA result reporting",
    "Performance",
    "Screen output",
    "Run continuation",
    "Configuration",
    "Choose analyses",
)
_GCTA_METHOD_DESCRIPTIONS = {
    "fastbat_gene": (
        "Test SNPs assigned to each gene and its configured window with "
        "LD-aware fastBAT."
    ),
    "fastbat_segment": (
        "Test consecutive fixed-size genomic segments with LD-aware fastBAT."
    ),
    "fastbat_set": (
        "Test custom SNP sets supplied through --fastbat-set-list or converted "
        "from --gmt pathways."
    ),
    "mbat_combo": (
        "Test each gene by combining signed mBAT and unsigned fastBAT evidence."
    ),
}
_GCTA_METHOD_SCOPED_DESTINATIONS = {
    "gene_list",
    "fastbat_set_list",
    "gmt",
    "gene_window_kb",
    "fastbat_segment_size_kb",
    "mbat_svd_gamma",
    "print_component_p_values",
    "gmt_chromosome_label_policy",
    "gmt_duplicate_gene_policy",
    "gmt_unmapped_gene_policy",
    "gcta_minimum_gene_id_overlap",
    "gmt_empty_pathway_policy",
    "fastbat_oversized_set_policy",
}
_GCTA_COMMON_ANALYSIS_DESTINATIONS = {
    "gcta_gene_method",
    "gcta_reference_maf_min",
    "fastbat_ld_cutoff",
    "frequency_difference_max",
    "write_snpset",
}
_GCTA_ANALYSIS_SCOPE_DESTINATIONS = {
    "mhc_policy",
    "mhc_chrom",
    "mhc_start",
    "mhc_end",
    "exclude_chromosomes",
}
_GCTA_METHOD_VISIBLE_DESTINATIONS = {
    "fastbat_gene": _GCTA_COMMON_ANALYSIS_DESTINATIONS | {
        "gene_list",
        "gene_window_kb",
    },
    "fastbat_segment": _GCTA_COMMON_ANALYSIS_DESTINATIONS | {
        "fastbat_segment_size_kb",
    },
    "mbat_combo": _GCTA_COMMON_ANALYSIS_DESTINATIONS | {
        "gene_list",
        "gene_window_kb",
        "mbat_svd_gamma",
        "print_component_p_values",
    },
}


def _argument_has_option(arguments: tuple[str, ...], option: str) -> bool:
    return any(
        argument == option or argument.startswith(option + "=")
        for argument in arguments
    )


def _method_help(
    default_method: str,
    selected_method: str | None = None,
    set_source: str | None = None,
) -> str:
    if selected_method is None:
        descriptions = "\n".join(
            "  %s: %s" % (method, _GCTA_METHOD_DESCRIPTIONS[method])
            for method in get_args(GctaGeneMethod)
        )
        introduction = "Choose the tested unit and GCTA procedure:"
    else:
        description = _GCTA_METHOD_DESCRIPTIONS[selected_method]
        if selected_method == "fastbat_set" and set_source == "prepared":
            description = (
                "Test custom SNP sets supplied through a prepared fastBAT "
                "set list."
            )
        elif selected_method == "fastbat_set" and set_source == "gmt":
            description = (
                "Convert GMT pathways to exact reference-BIM variant sets and "
                "test them with fastBAT."
            )
        descriptions = "  %s: %s" % (
            selected_method,
            description,
        )
        introduction = "This help page is filtered for the selected procedure:"
    return _with_configured_default(
        "%s\n%s" % (introduction, descriptions),
        default_method,
    )


def _set_source_from_arguments(arguments: tuple[str, ...]) -> str | None:
    has_prepared = _argument_has_option(arguments, "--fastbat-set-list")
    has_gmt = _argument_has_option(arguments, "--gmt")
    if has_prepared == has_gmt:
        return None
    return "prepared" if has_prepared else "gmt"


def _customize_gcta_gene_pipeline_help(
    parser: argparse.ArgumentParser,
) -> None:
    """Filter GCTA-only pipeline help to an explicitly selected procedure."""
    arguments = tuple(getattr(parser, "_postgwas_help_arguments", ()))
    method = cli_option_value(arguments, "--method")
    if method not in get_args(GctaGeneMethod):
        return

    module = load_configuration().modules.gcta_gene
    set_source = _set_source_from_arguments(arguments)
    visible = set(_GCTA_METHOD_VISIBLE_DESTINATIONS.get(method, ()))
    visible.update(_GCTA_ANALYSIS_SCOPE_DESTINATIONS)
    if method == "fastbat_set":
        visible.update(_GCTA_COMMON_ANALYSIS_DESTINATIONS)
        visible.update({"gmt_empty_pathway_policy", "fastbat_oversized_set_policy"})
        if set_source in {None, "prepared"}:
            visible.add("fastbat_set_list")
        if set_source in {None, "gmt"}:
            visible.update({
                "gmt",
                "gene_list",
                "gene_window_kb",
                "gmt_chromosome_label_policy",
                "gmt_duplicate_gene_policy",
                "gmt_unmapped_gene_policy",
                "gcta_minimum_gene_id_overlap",
            })

    actions = {action.dest: action for action in parser._actions}
    for destination in _GCTA_METHOD_SCOPED_DESTINATIONS:
        if destination not in visible and destination in actions:
            actions[destination].required = False
            actions[destination].help = argparse.SUPPRESS

    groups = {group.title: group for group in parser._action_groups}
    if method == "fastbat_set":
        groups["Input files"].description = (
            "Harmonised GWAS summary statistics and the selected custom-set "
            "input."
        )
    else:
        groups["Input files"].description = (
            "Harmonised GWAS summary statistics supplied to the pipeline."
        )
    groups["Reference files"].description = (
        "PLINK LD and gene-coordinate resources required by this procedure."
        if "gene_list" in visible
        else "PLINK LD resources required by this procedure."
    )

    method_action = actions["gcta_gene_method"]
    method_action.help = _method_help(
        module.method,
        selected_method=method,
        set_source=set_source,
    )

    gene_list_action = actions.get("gene_list")
    if gene_list_action is not None and gene_list_action.help is not argparse.SUPPRESS:
        if method == "fastbat_set":
            description = (
                "Gene-coordinate file used to convert GMT genes into "
                "BIM-compatible variant sets."
            )
        else:
            description = (
                "Gene-coordinate file with configured columns used to define "
                "the tested genes: %s."
                % ", ".join(module.gene_annotation.columns)
            )
        gene_list_action.help = description
        if method in {"fastbat_gene", "mbat_combo"} or set_source == "gmt":
            mark_cli_required_help(parser, ("gene_list",))
        elif method == "fastbat_set":
            gene_list_action.help = help_with_conditional_requirement(
                description,
                "when --gmt is selected",
            )

    if method == "fastbat_set":
        set_list_action = actions["fastbat_set_list"]
        gmt_action = actions["gmt"]
        if set_source == "prepared":
            set_list_action.help = (
                "Prepared custom fastBAT set file containing a set ID, SNP "
                "IDs, and END."
            )
            mark_cli_required_help(parser, ("fastbat_set_list",))
        elif set_source == "gmt":
            gmt_action.help = (
                "GMT pathway file converted to a fastBAT set list using "
                "--gene-list and exact reference-BIM identifiers."
            )
            mark_cli_required_help(parser, ("gmt",))
        else:
            set_list_action.help = help_with_conditional_requirement(
                "Prepared custom fastBAT set file containing a set ID, SNP "
                "IDs, and END.",
                "with --method fastbat_set unless --gmt is supplied",
            )
            gmt_action.help = help_with_conditional_requirement(
                "GMT pathway file converted using --gene-list and exact "
                "reference-BIM identifiers.",
                "with --method fastbat_set unless --fastbat-set-list is supplied",
            )

    parser._postgwas_pipeline_examples = get_gcta_gene_pipeline_examples(
        method=method,
        set_source=set_source,
    )
    method_step_descriptions = {
        "fastbat_gene": "Run GCTA gene-based fastBAT association analysis.",
        "fastbat_segment": "Run GCTA fixed-segment fastBAT association analysis.",
        "fastbat_set": (
            "Convert GMT pathways and run GCTA custom-set fastBAT analysis."
            if set_source == "gmt"
            else "Run GCTA custom-set fastBAT association analysis."
        ),
        "mbat_combo": "Run gene-based GCTA mBAT-combo association analysis.",
    }
    parser._postgwas_pipeline_step_descriptions = {
        "gcta_gene": method_step_descriptions[method],
    }


def organize_gcta_gene_help(parser: argparse.ArgumentParser) -> None:
    """Place GCTA outputs, analysis inputs, and references before settings."""
    _customize_gcta_gene_pipeline_help(parser)
    original_groups = list(parser._action_groups)
    output = next(group for group in original_groups if group.title == "Output")
    inputs = next(
        group for group in original_groups if group.title == "Input files"
    )
    references = next(
        group for group in original_groups if group.title == "Reference files"
    )
    reference_settings = next(
        group
        for group in original_groups
        if group.title == "GCTA reference settings"
    )

    move_cli_help_actions(parser, output, _GCTA_OUTPUT_DESTINATIONS)
    move_cli_help_actions(parser, inputs, _GCTA_INPUT_FILE_DESTINATIONS)
    move_cli_help_actions(
        parser,
        references,
        _GCTA_REFERENCE_FILE_DESTINATIONS,
    )
    move_cli_help_actions(
        parser,
        reference_settings,
        _GCTA_REFERENCE_SETTING_DESTINATIONS,
    )

    leading_groups = [
        group
        for group in original_groups
        if group.title in {"positional arguments", "options"}
    ]
    ordered_groups = [
        group
        for title in _GCTA_HELP_GROUP_ORDER
        for group in parser._action_groups
        if group.title == title
    ]
    remaining_groups = [
        group
        for group in parser._action_groups
        if group not in leading_groups and group not in ordered_groups
        and not (
            group.title in {"Input file", "GCTA inputs"}
            and not any(
                action.help is not argparse.SUPPRESS
                for action in group._group_actions
            )
        )
    ]
    parser._action_groups[:] = (
        leading_groups + ordered_groups + remaining_groups
    )


def _with_configured_default(description: str, value) -> str:
    return help_with_default(
        description,
        value,
    )


def get_gcta_gene_parser(add_help=False, *, direct_controls=False):
    """Return reusable GCTA options without assigning CLI-owned defaults."""
    defaults = load_configuration()
    module = defaults.modules.gcta_gene
    parser = argparse.ArgumentParser(add_help=add_help)
    inputs = parser.add_argument_group(
        "Input files",
        "Primary GWAS summary statistics and optional custom-set or pathway data.",
    )
    references = parser.add_argument_group(
        "Reference files",
        "PLINK LD and gene-coordinate resources used by the selected GCTA method.",
    )
    reference_settings = parser.add_argument_group(
        "GCTA reference settings",
        "Declarations used to validate compatibility across coordinate-bearing inputs.",
    )
    if direct_controls:
        inputs.add_argument(
            "--gcta-input-file",
            type=validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True),
            metavar="PATH",
            default=argparse.SUPPRESS,
            help=(
                "Formatter-created GCTA summary-statistics input for the "
                "selected method."
            ),
        )
    references.add_argument(
        "--gcta-reference-prefix",
        metavar="PREFIX",
        default=argparse.SUPPRESS,
        help=(
            "LD-reference prefix with configured companion files: %s."
            % ", ".join(module.reference.required_extensions)
        ),
    )
    reference_settings.add_argument(
        "--genome-build",
        choices=list(defaults.resources.genomes),
        metavar="BUILD",
        default=argparse.SUPPRESS,
        help=(
            "Shared genome build of the harmonised GWAS, PLINK LD reference, "
            "and gene list when used; all coordinate-bearing inputs must match."
        ),
    )
    reference_settings.add_argument(
        "--gcta-reference-population",
        choices=list(defaults.resources.populations),
        metavar="CODE",
        default=argparse.SUPPRESS,
        help=(
            "Population represented by the LD reference; no ancestry is inferred."
        ),
    )
    references.add_argument(
        "--gene-list",
        type=validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True),
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=help_with_conditional_requirement(
            "Gene-coordinate file with configured columns: %s."
            % ", ".join(module.gene_annotation.columns),
            "for fastbat_gene, mbat_combo, or fastbat_set with --gmt",
        ),
    )
    inputs.add_argument(
        "--fastbat-set-list",
        type=validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True),
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "Optional prepared custom fastBAT set file containing a set ID, "
            "SNP IDs, and END. Use this or --gmt, not both."
        ),
    )
    inputs.add_argument(
        "--gmt",
        type=validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True),
        metavar="PATH",
        default=argparse.SUPPRESS,
        help=(
            "Optional GMT pathway file converted to a fastBAT set list using "
            "--gene-list and the reference BIM; mutually exclusive with "
            "--fastbat-set-list."
        ),
    )
    settings = parser.add_argument_group("GCTA analysis settings")
    exclusions = parser.add_argument_group("GCTA analysis scope")
    matching = parser.add_argument_group("GCTA variant matching")
    pathway = parser.add_argument_group("GCTA pathway conversion")
    reporting = parser.add_argument_group("GCTA result reporting")
    settings.add_argument(
        "--method",
        dest="gcta_gene_method",
        choices=get_args(GctaGeneMethod),
        metavar="METHOD",
        default=argparse.SUPPRESS,
        help=_method_help(module.method),
    )
    exclusions.add_argument(
        "--mhc-policy",
        choices=MHC_POLICIES,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "MHC handling: include keeps all MHC SNPs and coordinate-defined "
            "genes; exclude_snps removes MHC SNPs before testing; "
            "exclude_genes removes overlapping genes from gene tests and "
            "GMT-derived pathway sets; exclude_both applies both exclusions. "
            "Gene-only exclusion is unavailable for fixed segments and prepared "
            "SNP-set files because they contain no gene identities",
            module.mhc.policy,
        ),
    )
    mhc_defaults = {
        str(build): genome.regions.get("mhc")
        for build, genome in defaults.resources.genomes.items()
    }

    def mhc_resource_default(field: str) -> str:
        return "; ".join(
            "%s=%s" % (build, getattr(region, field))
            for build, region in mhc_defaults.items()
            if region is not None
        )

    exclusions.add_argument(
        "--mhc-chrom",
        metavar="CHROM",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Override the MHC chromosome. Supply all three MHC override options "
            "together; when omitted, use resources.genomes.<BUILD>.regions.mhc",
            mhc_resource_default("chromosome"),
        ),
    )
    exclusions.add_argument(
        "--mhc-start",
        metavar="POSITION",
        type=int,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Override the one-based inclusive MHC start position. Supply all "
            "three MHC override options together; when omitted, use "
            "resources.genomes.<BUILD>.regions.mhc",
            mhc_resource_default("start"),
        ),
    )
    exclusions.add_argument(
        "--mhc-end",
        metavar="POSITION",
        type=int,
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Override the one-based inclusive MHC end position. Supply all three "
            "MHC override options together; when omitted, use "
            "resources.genomes.<BUILD>.regions.mhc",
            mhc_resource_default("end"),
        ),
    )
    exclusions.add_argument(
        "--exclude-chromosomes",
        metavar="CHROM",
        nargs="+",
        default=argparse.SUPPRESS,
        help=help_with_default(
            "Chromosomes excluded from LD-reference variants, gene tests, and "
            "GMT-derived pathway sets. X is retained unless explicitly listed",
            " ".join(module.chromosomes.exclude),
        ),
    )
    settings.add_argument(
        "--gene-window-kb", type=int, metavar="KB", default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Window added to both sides of each gene boundary",
            module.gene_window_kb,
        ),
    )
    pathway.add_argument(
        "--gmt-chromosome-label-policy",
        choices=get_args(GctaChromosomeLabelPolicy),
        metavar="POLICY",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "How GMT conversion reconciles gene-list and BIM chromosome labels",
            module.set_annotation.conversion.chromosome_label_policy,
        ),
    )
    pathway.add_argument(
        "--gmt-duplicate-gene-policy",
        choices=get_args(GctaDuplicateGenePolicy),
        metavar="POLICY",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "How repeated gene symbols within a GMT pathway are handled",
            module.set_annotation.conversion.duplicate_gene_policy,
        ),
    )
    pathway.add_argument(
        "--gmt-unmapped-gene-policy",
        choices=get_args(GctaUnmappedGenePolicy),
        metavar="POLICY",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Whether GMT genes absent from the coordinate file fail or are reported",
            module.set_annotation.conversion.unmapped_gene_policy,
        ),
    )
    pathway.add_argument(
        "--minimum-gene-id-overlap",
        dest="gcta_minimum_gene_id_overlap",
        type=float,
        metavar="FRACTION",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Minimum fraction of gene-coordinate identifiers that must occur "
            "exactly among the unique GMT genes. GMT-side mappability is "
            "reported but is not thresholded",
            module.set_annotation.conversion.minimum_gene_id_overlap_fraction,
        ),
    )
    pathway.add_argument(
        "--gmt-empty-pathway-policy",
        choices=get_args(GctaEmptyPathwayPolicy),
        metavar="POLICY",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Whether pathways with no usable variants after BIM/GWAS "
            "intersection fail or are omitted",
            module.set_annotation.conversion.empty_pathway_policy,
        ),
    )
    pathway.add_argument(
        "--fastbat-oversized-set-policy",
        choices=get_args(GctaOversizedSetPolicy),
        metavar="POLICY",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "How sets above the configured GCTA limit of %d variants are handled"
            % module.set_annotation.maximum_set_variants,
            module.set_annotation.oversized_set_policy,
        ),
    )
    settings.add_argument(
        "--fastbat-segment-size-kb", type=int, metavar="KB",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Fixed segment size",
            module.segment_size_kb,
        ),
    )
    if not direct_controls:
        matching.add_argument(
            "--gcta-minimum-reference-overlap",
            type=float,
            metavar="FRACTION",
            default=argparse.SUPPRESS,
            help=_with_configured_default(
                "Minimum fraction of unique pipeline variants that must "
                "resolve to exact PLINK BIM IDs",
                module.variant_harmonisation.minimum_overlap_fraction,
            ),
        )
    matching.add_argument(
        "--gcta-chromosome-label-policy",
        choices=get_args(GctaChromosomeLabelPolicy),
        metavar="POLICY",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "How PLINK-reference and gene-coordinate chromosome labels are "
            "normalized before compatibility checks",
            module.variant_harmonisation.chromosome_label_policy,
        ),
    )
    matching.add_argument(
        "--gcta-allow-strand-complement",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Allow complementary single-nucleotide allele pairs when only "
            "direct ID matching is available",
            module.variant_harmonisation.allow_strand_complement,
        ),
    )
    settings.add_argument(
        "--gcta-reference-maf-min", type=float, metavar="FRACTION",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Minimum MAF retained from the LD reference",
            module.reference_maf_min,
        ),
    )
    settings.add_argument(
        "--fastbat-ld-cutoff", type=float, metavar="R2", default=argparse.SUPPRESS,
        help=_with_configured_default(
            "LD r-squared pruning cutoff used by fastBAT",
            module.fastbat_ld_cutoff,
        ),
    )
    settings.add_argument(
        "--mbat-svd-gamma", type=float, metavar="FRACTION",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Minimum LD variance retained by mBAT principal components",
            module.mbat_svd_gamma,
        ),
    )
    settings.add_argument(
        "--frequency-difference-max", type=float, metavar="FRACTION",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Maximum A1-frequency difference between GWAS and LD reference",
            module.frequency_difference_max,
        ),
    )
    settings.add_argument(
        "--print-component-p-values",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Include mBAT and fastBAT component p-values",
            module.print_component_p_values,
        ),
    )
    settings.add_argument(
        "--write-snpset",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Write the SNP sets included in the selected fastBAT or mBAT analysis",
            module.write_snpset,
        ),
    )
    reporting.add_argument(
        "--gcta-top-results", type=int, metavar="N", default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Number of lowest-p-value associations shown in the terminal summary",
            module.reporting.top_result_count,
        ),
    )
    reporting.add_argument(
        "--gcta-nominal-alpha",
        type=float,
        metavar="ALPHA",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Nominal p-value threshold reported for individual associations",
            module.reporting.nominal_alpha,
        ),
    )
    reporting.add_argument(
        "--gcta-reporting-alpha",
        type=float,
        metavar="ALPHA",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Family-wise alpha used only to report the Bonferroni threshold; "
            "results are never filtered by this value",
            module.reporting.familywise_alpha,
        ),
    )
    reporting.add_argument(
        "--gcta-fdr-alpha",
        type=float,
        metavar="ALPHA",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Benjamini-Hochberg false-discovery-rate threshold",
            module.reporting.fdr_alpha,
        ),
    )
    reporting.add_argument(
        "--gcta-p-value-digits",
        type=int,
        metavar="N",
        default=argparse.SUPPRESS,
        help=_with_configured_default(
            "Significant digits used to display association p-values",
            module.reporting.p_value_significant_digits,
        ),
    )
    if direct_controls:
        configuration = parser.add_argument_group("Configuration")
        configuration.add_argument(
            "--run-config", metavar="PATH", default=argparse.SUPPRESS,
            help="YAML settings; explicit command-line values override matching keys.",
        )
        configuration.add_argument(
            "--dry-run", action="store_true", default=argparse.SUPPRESS,
            help="Validate inputs and record the exact command without running analysis.",
        )
    return parser


def get_gcta_gene_pipeline_examples(
    method: str | None = None,
    set_source: str | None = None,
):
    """Return distinct GCTA analysis examples for pipeline help."""
    examples = (
        (
            "Run gene-based fastBAT:",
            "postgwas pipeline",
            (
                "--modules gcta_gene",
                "--vcf study.vcf.gz",
                "--dataset-id STUDY",
                "--output-directory results",
                "--method fastbat_gene",
                "--gcta-reference-prefix reference/1000G_EUR",
                "--gene-list genes_grch37.txt",
                "--genome-build GRCh37",
                "--gcta-reference-population EUR",
            ),
        ),
        (
            "Run fixed-segment fastBAT:",
            "postgwas pipeline",
            (
                "--modules gcta_gene",
                "--vcf study.vcf.gz",
                "--dataset-id STUDY",
                "--output-directory results",
                "--method fastbat_segment",
                "--gcta-reference-prefix reference/1000G_EUR",
                "--genome-build GRCh37",
                "--gcta-reference-population EUR",
            ),
        ),
        (
            "Run custom-set fastBAT from a prepared set list:",
            "postgwas pipeline",
            (
                "--modules gcta_gene",
                "--vcf study.vcf.gz",
                "--dataset-id STUDY",
                "--output-directory results",
                "--method fastbat_set",
                "--gcta-reference-prefix reference/1000G_EUR",
                "--fastbat-set-list pathways.fastbat.set",
                "--genome-build GRCh37",
                "--gcta-reference-population EUR",
            ),
        ),
        (
            "Convert a GMT pathway file and run custom-set fastBAT:",
            "postgwas pipeline",
            (
                "--modules gcta_gene",
                "--vcf study.vcf.gz",
                "--dataset-id STUDY",
                "--output-directory results",
                "--method fastbat_set",
                "--gcta-reference-prefix reference/1000G_EUR",
                "--gene-list genes_grch37.txt",
                "--gmt pathways.gmt",
                "--genome-build GRCh37",
                "--gcta-reference-population EUR",
            ),
        ),
        (
            "Run gene-based mBAT-combo:",
            "postgwas pipeline",
            (
                "--modules gcta_gene",
                "--vcf study.vcf.gz",
                "--dataset-id STUDY",
                "--output-directory results",
                "--method mbat_combo",
                "--gcta-reference-prefix reference/1000G_EUR",
                "--gene-list genes_grch37.txt",
                "--genome-build GRCh37",
                "--gcta-reference-population EUR",
            ),
        ),
    )
    if method is None:
        return examples
    selected = tuple(
        example
        for example in examples
        if "--method %s" % method in example[2]
    )
    if method != "fastbat_set" or set_source is None:
        return selected
    required_source = "--gmt" if set_source == "gmt" else "--fastbat-set-list"
    return tuple(
        example
        for example in selected
        if any(argument.startswith(required_source + " ") for argument in example[2])
    )


def build_parser() -> argparse.ArgumentParser:
    defaults = load_configuration()
    formatting = defaults.modules.formatting.exports["gcta_gene"].outputs
    input_name = formatting["summary_statistics"].output_file.format(
        dataset_id="STUDY",
    )
    parser = argparse.ArgumentParser(
        prog="postgwas gcta_gene",
        usage=(
            "postgwas gcta_gene --gcta-input-file PATH "
            "--gcta-reference-prefix PREFIX --genome-build BUILD "
            "--gcta-reference-population CODE --method METHOD [options]"
        ),
        description=(
            "Run GCTA fastBAT gene, segment, or custom-set association, or "
            "gene-based mBAT-combo, using a "
            "genome-build- and population-matched PLINK LD reference."
        ),
        epilog=format_cli_examples(
            (
                "Run gene-based fastBAT:",
                "postgwas gcta_gene",
                (
                    "--method fastbat_gene",
                    "--gcta-input-file formatted/%s" % input_name,
                    "--gcta-reference-prefix reference/1000G_EUR",
                    "--gene-list genes_grch37.txt",
                    "--genome-build GRCh37",
                    "--gcta-reference-population EUR",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Run fixed-segment fastBAT:",
                "postgwas gcta_gene",
                (
                    "--method fastbat_segment",
                    "--gcta-input-file formatted/%s" % input_name,
                    "--gcta-reference-prefix reference/1000G_EUR",
                    "--genome-build GRCh37",
                    "--gcta-reference-population EUR",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Run custom-set fastBAT from a prepared set list:",
                "postgwas gcta_gene",
                (
                    "--method fastbat_set",
                    "--gcta-input-file formatted/%s" % input_name,
                    "--gcta-reference-prefix reference/1000G_EUR",
                    "--fastbat-set-list pathways.fastbat.set",
                    "--genome-build GRCh37",
                    "--gcta-reference-population EUR",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Convert a GMT file and run custom-set fastBAT:",
                "postgwas gcta_gene",
                (
                    "--method fastbat_set",
                    "--gcta-input-file formatted/%s" % input_name,
                    "--gcta-reference-prefix reference/1000G_EUR",
                    "--gene-list genes_grch37.txt",
                    "--gmt pathways.gmt",
                    "--genome-build GRCh37",
                    "--gcta-reference-population EUR",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Run gene-based mBAT-combo:",
                "postgwas gcta_gene",
                (
                    "--method mbat_combo",
                    "--gcta-input-file formatted/%s" % input_name,
                    "--gcta-reference-prefix reference/1000G_EUR",
                    "--gene-list genes_grch37.txt",
                    "--genome-build GRCh37",
                    "--gcta-reference-population EUR",
                    "--dataset-id STUDY",
                    "--output-directory results",
                ),
            ),
            (
                "Export reusable module settings:",
                "postgwas config export",
                ("--module gcta_gene", "--style full", "--output gcta_gene.yaml"),
            ),
            (
                "Run entirely from exported settings:",
                "postgwas gcta_gene",
                ("--run-config gcta_gene.yaml",),
            ),
        ),
        formatter_class=AlignedRichHelpFormatter,
        parents=[
            get_compute_parser(),
            get_common_out_parser(),
            get_gcta_gene_parser(direct_controls=True),
        ],
    )
    for action in parser._actions:
        if action.dest == "dataset_id":
            action.help = _with_configured_default(
                "Short, unique name used in output filenames",
                defaults.run.dataset_id,
            )
            action.default = argparse.SUPPRESS
        elif action.dest == "output_directory":
            action.help = _with_configured_default(
                "Folder where PostGWAS saves results and logs",
                defaults.run.output_directory,
            )
            action.default = argparse.SUPPRESS
    mark_cli_required_help(
        parser,
        (
            "gcta_input_file",
            "gcta_reference_prefix",
            "genome_build",
            "gcta_reference_population",
        ),
    )
    organize_gcta_gene_help(parser)
    return parser


def main(argv=None):
    from postgwas.modules.gcta_gene.service import run_gcta_gene_direct

    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        run_gcta_gene_direct(args)
    except MissingRequiredArgumentsError as exc:
        console = Console(stderr=True)
        print_screen_message(
            "error", "GCTA gene analysis failed. %s" % exc, console=console,
        )
        parser.print_help(file=console.file)
        return 1
    except (GctaGeneError, ConfigurationError, OSError, ValueError) as exc:
        print_screen_message(
            "error", "GCTA gene analysis failed. %s" % exc, stderr=True,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "build_parser",
    "get_gcta_gene_parser",
    "get_gcta_gene_pipeline_examples",
    "main",
    "organize_gcta_gene_help",
]
