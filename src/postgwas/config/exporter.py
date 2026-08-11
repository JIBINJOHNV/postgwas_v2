"""Render reproducible user configuration without duplicating defaults."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from postgwas.config.loader import (
    canonical_module_name,
    load_configuration,
    load_module_configuration,
)
from postgwas.config.models.modules.single_cell import (
    single_cell_formatter_targets,
    single_cell_pipeline_dependencies,
    single_cell_supporting_configurations,
)
from postgwas.core.errors import ConfigurationError
from postgwas.pipeline.planner import build_pipeline_plan
from postgwas.modules.formatting.contracts import required_formats


EXPORT_STYLES = ("full", "minimal", "values")


class _CompactDumper(yaml.SafeDumper):
    """Keep sequences such as chromosome sets on one readable YAML line."""


def _represent_compact_sequence(dumper, value):
    return dumper.represent_sequence("tag:yaml.org,2002:seq", value, flow_style=True)


_CompactDumper.add_representer(list, _represent_compact_sequence)
_CompactDumper.add_representer(tuple, _represent_compact_sequence)


def _dump_values(value: dict[str, Any]) -> str:
    return yaml.dump(
        value,
        Dumper=_CompactDumper,
        sort_keys=False,
        default_flow_style=False,
        allow_unicode=True,
        width=4096,
    )


def _indent(text: str, spaces: int) -> str:
    prefix = " " * spaces
    return "\n".join(prefix + line if line else "" for line in text.splitlines())


def _render_harmonisation(module_config, style: str) -> str:
    from postgwas.modules.harmonisation.policies import as_yaml_template, load_policies

    values = module_config.model_dump(mode="json")
    overrides = values.pop("policies", {})
    policies = load_policies(overrides)
    policy_text = as_yaml_template(
        policies=policies,
        style=style,
        include_header=False,
        include_required_inputs=False,
    )
    return _dump_values(values).rstrip() + "\n" + policy_text


_FORMATTING_COMMENTS = {
    "enabled": (
        "Whether the formatter is enabled when this file is used in a pipeline.",
    ),
    "formats": (
        "Downstream inputs to create. CLI --format overrides this list.",
    ),
    "format_order": (
        "Stable order used when several configured downstream formats are required.",
    ),
    "module_formats": (
        "Pipeline module -> formatter artifacts required before that module runs.",
    ),
    "overwrite": (
        "Replace existing formatter outputs. False protects completed results.",
    ),
    "chromosomes": (
        "Chromosome files written for PRED-LD; order is retained in the run.",
    ),
    "minimum_p_value": (
        "Smallest raw p-value written when 10^(-LP) would underflow.",
        "Every bounded value is counted in the formatter audit log.",
    ),
    "vcf_fields": (
        "Canonical formatter-table column -> bcftools query expression.",
        "This is the single source of truth for structural and FORMAT-field extraction.",
    ),
    "numeric_columns": ("Extracted canonical columns parsed as numeric values.",),
    "canonical_columns": ("Semantic roles used to validate and identify variants.",),
    "chromosome_labels": (
        "Normalization applied to chromosome labels before any target is exported.",
    ),
    "chromosome_labels.prefix_pattern": (
        "Configured leading chromosome prefix removed from VCF chromosome labels.",
    ),
    "chromosome_labels.aliases": (
        "Optional exact chromosome-label replacements applied after prefix removal.",
    ),
    "variant_identifiers": (
        "General variant-ID policy applied independently to each formatter target.",
        "Pipeline modules may set one target type after inspecting their reference.",
    ),
    "variant_identifiers.default_type": (
        "Identifier type used unless CLI or a target-specific YAML value overrides it.",
    ),
    "variant_identifiers.target_types": (
        "Optional formatter target -> identifier type overrides.",
    ),
    "variant_identifiers.rsid_pattern": (
        "Pattern used to identify rsIDs in an external reference.",
    ),
    "variant_identifiers.rsid_extraction_pattern": (
        "One-capture-group pattern used to extract an rsID from the VCF ID field.",
    ),
    "variant_identifiers.unique_id_template": (
        "Template used to construct an ID from VCF chromosome, position, REF, and ALT.",
    ),
    "study_design": (
        "Columns and target formats used to infer binary versus quantitative traits.",
    ),
    "study_design.required_formats": (
        "Formatter targets that require trait-specific sample-size interpretation.",
        "The validated scientific contract requires exactly LDSC and MiXeR.",
    ),
    "runtime": ("Configured formatter logs, metadata paths, and I/O settings.",),
    "vcf_include_expression": ("Optional bcftools expression applied during extraction.",),
    "exports": (
        "Per-tool canonical-column -> output-column mappings, validation, transformations, and paths.",
    ),
    "mixer": ("MiXeR-specific validity and recommended pre-analysis QC checks.",),
    "mixer.minimum_info": (
        "Minimum imputation INFO; null disables this export-time check.",
    ),
    "mixer.minimum_sample_size_fraction": (
        "Keep N at or above this fraction of the median N (MiXeR recommends 0.5).",
    ),
    "mixer.snps_only": ("Keep single-nucleotide variants only for MiXeR reference matching.",),
    "mixer.allowed_alleles": ("Allele symbols accepted when SNP-only checking is active.",),
}


def _render_commented(module_config, style: str, comments_by_path) -> str:
    """Render ordered values with optional comments from a small path registry."""
    values = (
        module_config
        if isinstance(module_config, dict)
        else module_config.model_dump(mode="json")
    )
    text = _dump_values(values)
    if style == "values":
        return text
    stack: list[tuple[int, str]] = []
    rendered = []
    for line in text.splitlines():
        indentation = len(line) - len(line.lstrip(" "))
        key = line.strip().split(":", 1)[0]
        while stack and stack[-1][0] >= indentation:
            stack.pop()
        path = ".".join([item[1] for item in stack] + [key])
        comments = comments_by_path.get(path, ())
        if style == "minimal":
            comments = comments[:1]
        rendered.extend(" " * indentation + "# " + comment for comment in comments)
        rendered.append(line)
        if line.rstrip().endswith(":"):
            stack.append((indentation, key))
    return "\n".join(rendered) + "\n"


def _render_formatting(module_config, style: str) -> str:
    return _render_commented(module_config, style, _FORMATTING_COMMENTS)


_MIXER_COMMENTS = {
    "enabled": ("Enable single-trait MiXeR when this module is part of a pipeline.",),
    "analysis": (
        "Run univariate architecture, GSA-MiXeR gene-set enrichment, or both.",
        "Every choice analyzes one trait; multi-trait commands are rejected.",
    ),
    "execution_backend": (
        "Use native MiXeR, the configured GSA-MiXeR container, or choose automatically.",
    ),
    "genome_build": (
        "Genome build shared by the formatted coordinates and MiXeR reference.",
        "MiXeR does not lift coordinates; a build mismatch invalidates reference matching.",
    ),
    "bim_file_pattern": (
        "Per-chromosome BIM path using MiXeR's literal @ placeholder.",
    ),
    "ld_file_pattern": (
        "Per-chromosome MiXeR LD path using the same @ placeholder.",
    ),
    "chromosomes": (
        "Ordered consecutive autosomes included in reference validation and --chr2use.",
    ),
    "fit_arguments": (
        "Additional official fit1 options. Protected and trait-2 options are rejected.",
    ),
    "test_arguments": (
        "Additional official test1 options. Protected and trait-2 options are rejected.",
    ),
    "gsa": (
        "Reference inputs, scientific options, schemas, and reporting policy for single-trait GSA-MiXeR.",
    ),
    "gsa.annotation_file_pattern": (
        "Per-chromosome SNP annotation using the configured chromosome placeholder.",
    ),
    "gsa.loadlib_file_pattern": (
        "Optional precomputed load library; null makes GSA-MiXeR use the configured LD files.",
    ),
    "gsa.baseline_go_file": (
        "Baseline GO-format annotation used by plsa --gsa-base.",
    ),
    "gsa.model_go_file": (
        "Gene-level GO-format annotation used to fit the full GSA model.",
    ),
    "gsa.test_go_file": (
        "GO-format gene sets tested for enrichment in the full GSA model.",
    ),
    "gsa.exclude_ranges": (
        "Named upstream exclusion ranges passed explicitly to every GSA model.",
    ),
    "gsa.hardprune_maf": ("Minimum MAF used by GSA-MiXeR hard pruning.",),
    "gsa.hardprune_r2": ("Maximum LD r² used by GSA-MiXeR hard pruning.",),
    "gsa.result_columns": (
        "Semantic result fields mapped to official go_test_enrich.csv columns.",
    ),
    "gsa.evidence_aic_threshold": (
        "Minimum loglike_aic required for model-support reporting; no p-values are inferred.",
    ),
    "workflow": (
        "Official command names, literal placeholders, and UTC run-ID format.",
    ),
    "output_layout": (
        "Relative raw-result, summary, figure, metadata, and per-run log paths.",
    ),
    "reporting": (
        "Compact result files, official diagnostic figures, and quality actions.",
    ),
    "reporting.generate_figures": (
        "Run official mixer_figures after test1 to create QQ and power diagnostics.",
    ),
    "reporting.quality": (
        "Configured error, warning, or ignore actions for scientifically material checks.",
    ),
}


def _render_mixer(module_config, style: str) -> str:
    return _render_commented(module_config, style, _MIXER_COMMENTS)


_MAGMA_COMMENTS = {
    "enabled": ("Enable MAGMA when this module is selected in a pipeline.",),
    "genome_build": (
        "Genome build declared for provenance; file-level validation is deferred.",
        "Summary statistics, BIM, and gene locations must nevertheless use this build.",
    ),
    "population": (
        "LD-reference population declared for provenance; validation is deferred.",
    ),
    "gene_window_upstream_kb": ("Kilobases added upstream during SNP-to-gene annotation.",),
    "gene_window_downstream_kb": ("Kilobases added downstream during SNP-to-gene annotation.",),
    "gene_model": ("MAGMA gene model supported with SNP p-value input.",),
    "input": (
        "Module inputs and the canonical roles of every table column.",
        "Null paths must be supplied by pipeline artifacts, another YAML layer, or CLI.",
    ),
    "input.sample_size_column": (
        "Per-variant total sample size passed to MAGMA with ncol.",
    ),
    "input.chromosome_prefix_pattern": (
        "Configured leading prefix removed before chromosome labels are compared.",
    ),
    "input.chromosome_aliases": (
        "Configured exact aliases applied after prefix removal and uppercasing.",
    ),
    "input.invalid_chromosome_labels": (
        "Chromosome labels that make a MAGMA input or BIM record invalid.",
    ),
    "input.invalid_allele_labels": (
        "Allele labels that make a MAGMA input or BIM record invalid.",
    ),
    "input.required_reference_extensions": (
        "Companion files required for the configured PLINK reference prefix.",
    ),
    "input.bim_extension": ("Companion file containing the PLINK variant map.",),
    "input.bim_columns": ("Ordered semantic roles of the six PLINK BIM columns.",),
    "input.output_table_delimiter": (
        "Single-character delimiter written to prepared MAGMA input tables.",
    ),
    "snp_harmonisation": (
        "Optional exact reference intersection, overlap protection, and duplicate policy.",
    ),
    "snp_harmonisation.resolve_variants_to_reference": (
        "When true, retain only exact BIM field-2 IDs with matching coordinates and alleles.",
        "When false, use validated formatter tables directly without PostGWAS pre-filtering.",
    ),
    "snp_harmonisation.minimum_overlap_fraction": (
        "Minimum exact BIM-ID overlap required only when reference intersection is enabled.",
    ),
    "snp_harmonisation.duplicate_policy": (
        "Retain the lowest valid p-value when several rows resolve to one reference SNP.",
    ),
    "gene_sets": (
        "Gene-set input format and gene-ID compatibility required before testing.",
    ),
    "gene_sets.input_format": (
        "Accept standard GMT, native MAGMA set annotations, or detect either format.",
    ),
    "annotation_validation": (
        "Structural and exact BIM-ID checks applied to supplied functional annotations.",
    ),
    "annotation_validation.minimum_bim_variant_overlap_fraction": (
        "Minimum fraction of unique annotation variants required in BIM field 2.",
    ),
    "annotation_validation.invalid_variant_identifiers": (
        "Configured annotation placeholders excluded from assignment and overlap counts.",
    ),
    "mapping": (
        "Independent positional or functional SNP-to-unit analyses and their provenance.",
        "Results from different definitions remain separate hypothesis families.",
    ),
    "mapping.selected": (
        "Ordered mapping definitions to run; CLI --magma-mapping overrides this list.",
    ),
    "mapping.primary": (
        "Mapping used for top-level outputs; downstream consumers require a calibrated result.",
    ),
    "mapping.definitions": (
        "Named method, resource paths, identifiers, context, and statistic interpretation.",
    ),
    "chrom_magma_mapping": (
        "Configured regulatory-element schemas and published lowest-element-p assignment.",
    ),
    "chrom_magma_mapping.minimum_element_mapping_fraction": (
        "Minimum tested-element coverage required in the element-to-gene mapping.",
    ),
    "chrom_magma_mapping.gene_assignment_policy": (
        "Select the lowest-p linked element per gene; this is a ranking statistic.",
    ),
    "batching": ("Bound MAGMA batches by configured compute resources and gene count.",),
    "multiple_testing": (
        "Gene-level methods plus global and named-family gene-set corrections.",
    ),
    "multiple_testing.reporting_significance_threshold": (
        "Threshold used only to count and report nominal and adjusted significant results.",
    ),
    "multiple_testing.reporting_method_labels": (
        "User-facing labels for configured correction methods in stage outcomes.",
    ),
    "multiple_testing.gene_methods": (
        "Corrections applied across all valid gene p-values in one MAGMA result.",
    ),
    "result_schema": (
        "MAGMA result-column roles and PostGWAS report names and serialization.",
    ),
    "result_schema.report_gene_set_description_column": (
        "Output column containing the GMT description; null for native MAGMA input.",
    ),
    "result_schema.report_input_genes_column": (
        "Output column containing the gene IDs supplied for each tested set.",
    ),
    "minimum_magma_version": ("Oldest MAGMA version accepted by the runner.",),
    "version_arguments": ("Arguments used only to obtain MAGMA's version output.",),
    "version_pattern": ("Pattern whose first capture group is the MAGMA version.",),
    "output_layout": (
        "Relative log, staging, provenance, raw-result, and final-result paths.",
        "All analysis files are staged before successful outputs are published.",
    ),
}


def _render_magma(module_config, style: str) -> str:
    return _render_commented(module_config, style, _MAGMA_COMMENTS)


def render_module_configuration(
    module: str,
    *,
    config_file: str | Path | None = None,
    style: str = "minimal",
) -> str:
    """Render one complete module configuration as reloadable YAML."""
    if style not in EXPORT_STYLES:
        raise ConfigurationError("Unknown export style: %s" % style)
    name = canonical_module_name(module)
    config = load_module_configuration(name, config_file)
    if name == "harmonisation":
        return _render_harmonisation(config, style)
    if name == "formatting":
        return _render_formatting(config, style)
    if name == "mixer":
        return _render_mixer(config, style)
    if name == "magma":
        return _render_magma(config, style)
    return _dump_values(config.model_dump(mode="json"))


def render_pipeline_configuration(
    targets: list[str] | tuple[str, ...],
    *,
    config_file: str | Path | None = None,
    style: str = "minimal",
) -> str:
    """Render executable steps plus supporting method configuration."""
    if style not in EXPORT_STYLES:
        raise ConfigurationError("Unknown export style: %s" % style)
    config = load_configuration(config_file)
    single_cell_tools = (
        config.modules.single_cell.tools if "single_cell" in targets else []
    )
    dependency_overrides = (
        {
            "single_cell": single_cell_pipeline_dependencies(
                single_cell_tools
            )
        }
        if single_cell_tools else None
    )
    plan = build_pipeline_plan(
        targets, dependency_overrides=dependency_overrides,
    )
    selected = list(
        dict.fromkeys(canonical_module_name(step) for step in plan.steps)
    )
    supporting = single_cell_supporting_configurations(single_cell_tools)
    rendered_modules = []
    for name in selected:
        if name == "single_cell":
            rendered_modules.extend(supporting)
        rendered_modules.append(name)
    rendered_modules = list(dict.fromkeys(rendered_modules))

    values = config.model_dump(mode="json")
    all_modules = values.pop("modules")
    values["pipeline"]["modules"] = selected
    if "formatting" in selected:
        requested_formats = list(all_modules["formatting"].get("formats", ()))
        requested_formats.extend(
            single_cell_formatter_targets(single_cell_tools)
        )
        all_modules["formatting"]["formats"] = required_formats(
            config.modules.formatting,
            plan.active_modules,
            requested_formats,
        )
    output = _dump_values(values).rstrip() + "\nmodules:\n"
    for name in rendered_modules:
        if name in selected:
            all_modules[name]["enabled"] = True
        output += "  %s:\n" % name
        if name == "harmonisation":
            rendered = _render_harmonisation(getattr(config.modules, name), style)
        elif name == "formatting":
            rendered = _render_formatting(all_modules[name], style)
        elif name == "mixer":
            rendered = _render_mixer(all_modules[name], style)
        elif name == "magma":
            rendered = _render_magma(all_modules[name], style)
        else:
            rendered = _dump_values(all_modules[name])
        output += _indent(rendered.rstrip(), 4) + "\n"
    return output
