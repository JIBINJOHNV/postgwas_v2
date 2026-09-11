"""MAGMA workflow reports and concise terminal result summaries."""

from __future__ import annotations

import csv
from html import escape
import json
import os
from pathlib import Path
from typing import Any, Mapping, Sequence

import polars as pl

from postgwas.core.io.reports import write_delimited_report, write_html_report
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.core.values import format_number
from postgwas.core.vcf import vcf_query_field_label


MAGMA_PIPELINE_STAGES = (
    "Validate the input GWAS-VCF",
    "Validate the PLINK LD-reference files and determine the BIM variant-ID format",
    "Validate the gene-location file",
    "Validate the original pathway GMT file",
    "Compare pathway identifiers with gene-location columns and select the "
    "MAGMA gene ID",
    "Prepare variant-level inputs for MAGMA",
    "Create the MAGMA SNP-to-gene annotation",
    "Calculate MAGMA gene-association statistics",
    "Annotate gene results and adjust gene p-values",
    "Run competitive pathway analysis with the original pathway file",
    "Annotate pathway results and adjust pathway p-values",
    "Validate and publish all outputs",
)

MAGMA_PIPELINE_STAGE_KEYS = (
    "vcf",
    "ld_reference",
    "gene_location",
    "pathway_file",
    "pathway_identifiers",
    "variant_inputs",
    "gene_annotation",
    "gene_analysis",
    "gene_results",
    "pathway_analysis",
    "pathway_results",
    "publish",
)

MAGMA_GENE_ONLY_STAGE_KEYS = (
    "vcf",
    "ld_reference",
    "gene_location",
    "variant_inputs",
    "gene_annotation",
    "gene_analysis",
    "gene_results",
    "publish",
)

MAGMA_STAGE_TITLES = dict(zip(MAGMA_PIPELINE_STAGE_KEYS, MAGMA_PIPELINE_STAGES))


def _numbered_stages(keys: Sequence[str]) -> dict[str, int]:
    return {key: number for number, key in enumerate(keys, 1)}


def magma_pipeline_progress_plan(args) -> dict[str, Any] | None:
    """Return the canonical stage plan for a MAGMA-only pipeline."""
    from postgwas.modules.magma.service import _resolved_configuration

    module = _resolved_configuration(args).modules.magma
    if (
        len(module.mapping.selected) != 1
        or module.mapping.definitions[module.mapping.selected[0]].method != "positional"
    ):
        return None
    return {
        "kind": "magma",
        "label": "MAGMA pipeline execution progress",
        "stages": MAGMA_PIPELINE_STAGES,
        "magma_stage_numbers": _numbered_stages(MAGMA_PIPELINE_STAGE_KEYS),
        "modules": {"formatter": (1, 6), "magma": (6, 12)},
        "deferred_completion_modules": {"formatter": 6},
    }


def magma_pipeline_stage_numbers(args) -> Mapping[str, int] | None:
    """Return the active detailed MAGMA stage mapping, if one is attached."""
    plan = getattr(args, "_pipeline_progress_plan", None)
    if isinstance(plan, Mapping):
        stage_numbers = plan.get("magma_stage_numbers")
        if isinstance(stage_numbers, Mapping):
            return stage_numbers
    if getattr(args, "_pipeline_stage_progress", None) is not None:
        # Preserve programmatic callers that attach the original MAGMA
        # controller without routing through the pipeline executor.
        return _numbered_stages(MAGMA_PIPELINE_STAGE_KEYS)
    return None


def _number(value: Any) -> str:
    return "not available" if value is None else f"{int(value):,}"


def _percent(numerator: Any, denominator: Any) -> str:
    if numerator is None or denominator in {None, 0}:
        return "not available"
    return "%.2f%%" % (100.0 * int(numerator) / int(denominator))


def _record(
    step: int,
    stage: str,
    status: str,
    summary: str,
    details: Sequence[str] = (),
    outputs: Sequence[str | Path] = (),
) -> dict[str, Any]:
    return {
        "step": int(step),
        "stage": str(stage),
        "status": str(status),
        "summary": str(summary),
        "details": " | ".join(str(value) for value in details if value),
        "output": "; ".join(str(value) for value in outputs if value),
    }


def _pathway_plan(preflight, primary: str) -> Mapping[str, Any]:
    return preflight.gene_set_plans[primary]


def magma_variant_input_outcome_fields(
    configuration,
    formatter_result: Mapping[str, Any],
    variant_id_observation: Mapping[str, Any],
    qc: Mapping[str, Any],
    output_paths: Mapping[str, str | Path],
) -> list[tuple]:
    """Return one configuration-derived MAGMA input-preparation account."""
    formatting = configuration.modules.formatting
    magma = configuration.modules.magma
    canonical = formatting.canonical_columns
    fields = formatting.vcf_fields.root
    export = formatting.exports["magma"]
    location = export.outputs["snp_location"]
    association = export.outputs["p_values"]
    identifier = formatting.variant_identifiers.unique_id_template.format(
        chromosome=canonical.chromosome,
        position=canonical.position,
        reference_allele=canonical.reference_allele,
        alternate_allele=canonical.alternate_allele,
    )
    identifier_type = variant_id_observation.get(
        "variant_id_type", formatter_result.get("variant_id_type", "configured")
    )
    identifier_source = (
        " + ".join((
            canonical.chromosome,
            canonical.position,
            canonical.reference_allele,
            canonical.alternate_allele,
        )) + " → " + identifier
        if identifier_type == "unique" else
        "%s → %s" % (
            vcf_query_field_label(fields[canonical.variant_id]),
            canonical.resolved_variant_id,
        )
    )
    p_source = next(
        source for source, saved in association.columns.items()
        if saved == magma.input.p_value_column
    )
    n_source = next(
        source for source, saved in association.columns.items()
        if saved == magma.input.sample_size_column
    )
    p_transformation = association.transformations.get(p_source)
    p_processing = (
        "raw P = 10⁻ᴸᴾ"
        if p_transformation == "negative_log10_to_raw_p" else
        "copied" if p_transformation is None else
        "configured transformation: %s" % p_transformation
    )
    n_transformation = association.transformations.get(n_source)
    n_processing = (
        "per-variant total N"
        if n_transformation is None else
        "configured transformation: %s" % n_transformation
    )
    chromosome_source = next(
        source for source, saved in location.columns.items()
        if saved == magma.input.chromosome_column
    )
    position_source = next(
        source for source, saved in location.columns.items()
        if saved == magma.input.position_column
    )
    matched = int(qc["reference_unique_id_matches"])
    total = int(qc["input_unique_variants"])
    absent = int(qc["not_in_reference_rows"])
    intersection = bool(qc["reference_intersection_enabled"])
    excluded_chromosomes = list(qc["excluded_chromosomes"])
    region = qc.get("mhc_region")
    mhc_region = (
        "%s:%s–%s" % (region["chromosome"], region["start"], region["end"])
        if region is not None else "not configured"
    )
    duplicate_policy = {
        "lowest_p": "retain lowest-P record",
        "remove": "remove every duplicated-ID group",
        "err": "stop when duplicated IDs are present",
    }[qc["duplicate_policy"]]
    location_path = Path(output_paths["snp_loc_file"])
    association_path = Path(output_paths["pval_file"])
    output_directory = (
        str(location_path.parent)
        if location_path.parent == association_path.parent else
        "separate configured output locations"
    )
    source_fields = [
        ("analysis", "Extract and transform GWAS-VCF fields"),
        ("count", "Variants extracted", formatter_result["rows_in"]),
        ("genetic", magma.input.variant_id_column, identifier_source),
        (
            "info", magma.input.chromosome_column,
            "%s → normalised chromosome"
            % vcf_query_field_label(fields[chromosome_source]),
        ),
        (
            "info", magma.input.position_column,
            "%s → positive integer"
            % vcf_query_field_label(fields[position_source]),
        ),
        (
            "info", magma.input.p_value_column,
            "%s → %s"
            % (vcf_query_field_label(fields[p_source]), p_processing),
        ),
        (
            "info", magma.input.sample_size_column,
            "%s → %s"
            % (vcf_query_field_label(fields[n_source]), n_processing),
        ),
    ] if formatter_result else [
        ("analysis", "Read and validate supplied MAGMA input files"),
        ("count", "Variants read", total),
        (
            "info", "SNP-location columns",
            ", ".join((
                magma.input.variant_id_column,
                magma.input.chromosome_column,
                magma.input.position_column,
            )),
        ),
        (
            "info", "Association columns",
            ", ".join((
                magma.input.variant_id_column,
                magma.input.p_value_column,
                magma.input.sample_size_column,
            )),
        ),
    ]
    return source_fields + [
        ("analysis", "Compare variants with the LD reference"),
        ("count", "GWAS variants", total),
        (
            "success", "Present in BIM",
            "%s (%s)" % (_number(matched), _percent(matched, total)),
        ),
        (
            "loss" if absent else "success",
            "Absent from BIM",
            "%s (%s)" % (_number(absent), _percent(absent, total)),
        ),
        (
            "info", "Minimum required overlap",
            "%.2f%%" % (magma.snp_harmonisation.minimum_overlap_fraction * 100),
        ),
        ("success", "Identifier compatibility", "passed"),
        ("info", "BIM intersection", "applied" if intersection else "not applied"),
        (
            "info", "Unmatched variants",
            "excluded from MAGMA input" if intersection else "retained in MAGMA input",
        ),
        ("analysis", "Apply analysis-scope policies"),
        (
            "info", "Chromosome policy",
            "exclude %s" % " and ".join(excluded_chromosomes)
            if excluded_chromosomes else "include all chromosomes",
        ),
        (
            "loss" if qc["excluded_chromosome_rows"] else "success",
            "Variants excluded by chromosome", qc["excluded_chromosome_rows"],
        ),
        (
            "info", "MHC SNP policy",
            "exclude" if magma.mhc.excludes_snps else "include",
        ),
        ("info", "MHC region", mhc_region),
        (
            "loss" if qc["excluded_mhc_rows"] else "success",
            "Variants excluded from MHC", qc["excluded_mhc_rows"],
        ),
        ("info", "Duplicate-ID policy", duplicate_policy),
        (
            "loss" if qc["duplicate_groups_detected"] else "success",
            "Duplicate groups detected", qc["duplicate_groups_detected"],
        ),
        (
            "loss" if qc["duplicate_rows_removed"] else "success",
            "Duplicate rows removed", qc["duplicate_rows_removed"],
        ),
        ("success", "Variants retained", qc["retained_rows"]),
        ("analysis", "Write final MAGMA input files"),
        ("info", "Output directory", output_directory),
        ("info", "SNP-location file", location_path.name),
        (
            "info", "Saved columns",
            ", ".join((
                magma.input.variant_id_column,
                magma.input.chromosome_column,
                magma.input.position_column,
            )),
        ),
        ("success", "SNP-location records written", qc["retained_rows"]),
        ("info", "Association file", association_path.name),
        (
            "info", "Saved columns",
            ", ".join((
                magma.input.variant_id_column,
                magma.input.p_value_column,
                magma.input.sample_size_column,
            )),
        ),
        ("success", "Association records written", qc["retained_rows"]),
        ("info", "Retention after analysis", "retained"),
        ("info", "Purpose", "reproducibility, inspection and validated resume"),
    ]


def build_magma_pipeline_summary(
    dataset_id: str,
    configuration,
    preflight,
    result: Mapping[str, Any],
    *,
    formatter_result: Mapping[str, Any] | None = None,
    variant_id_observation: Mapping[str, Any] | None = None,
    input_vcf: str | Path | None = None,
    include_pathway_stages: bool = True,
) -> list[dict[str, Any]]:
    """Build the ordered pipeline account exclusively from validated metrics."""
    module = configuration.modules.magma
    formatter = dict(formatter_result or {})
    observation = dict(variant_id_observation or {})
    primary = result["primary_mapping"]
    mapping = result["mapping_analyses"][primary]
    plan = _pathway_plan(preflight, primary)
    validation = plan.get("validation") or {}
    resolution = validation.get("identifier_resolution") or {}
    qc = result["variant_preparation"]["qc"]
    input_metadata = formatter.get("input_vcf_metadata") or {}
    gene_significance = mapping.get("gene_significance") or {}
    pathway_significance = mapping.get("pathway_significance") or {}
    gene_scope = mapping.get("gene_scope") or {}
    pathway_location = (
        (mapping.get("annotation_validation") or {}).get(
            "pathway_compatible_gene_location"
        ) or {}
    )
    gene_adjusted = gene_significance.get("adjusted_significant") or {}
    primary_gene_count = resolution.get(
        "location_primary_unique_ids", gene_scope.get("input_units"),
    )
    alternate_gene_count = resolution.get("location_alternate_unique_ids")
    pipeline_mode = bool(formatter)
    bim_identifier_column = module.input.bim_columns.index("variant_id") + 1

    pathway_requested = plan["status"] != "not_requested"
    pathway_ready = plan["status"] == "ready"
    pathway_completed = mapping["gene_set_analysis"]["status"] == "completed"
    identifier_strategy = resolution.get("identifier_source")
    observed_variant_id_type = observation.get("variant_id_type", "the configured")
    variant_id_label = (
        "rsID"
        if observed_variant_id_type == "rsid" else
        "chromosome-position-allele ID"
        if observed_variant_id_type == "unique" else
        str(observed_variant_id_type)
    )
    if identifier_strategy == "alternate_gene_id":
        identifier_decision = (
            "a derived gene-location file uses alternate %s IDs as primary IDs; "
            "pathway identifiers remain unchanged"
            % module.input.alternate_gene_id_type
        )
    elif identifier_strategy == "primary_gene_id":
        identifier_decision = "primary gene identifiers used directly"
    else:
        identifier_decision = "identifier conversion not required"

    prepared_outputs = [
        result["variant_preparation"].get("snp_loc_file"),
        result["variant_preparation"].get("pval_file"),
    ]
    variant_input_fields = magma_variant_input_outcome_fields(
        configuration,
        formatter,
        observation,
        qc,
        {
            "snp_loc_file": prepared_outputs[0],
            "pval_file": prepared_outputs[1],
        },
    )
    records = [
        _record(
            1, "Validate the input GWAS-VCF",
            "completed" if pipeline_mode else "not requested",
            (
                "GWAS-VCF structural validation passed; %s variant records were read."
                % _number(formatter["rows_in"])
                if pipeline_mode else
                "Direct-module mode received prepared MAGMA tables; no GWAS-VCF was supplied."
            ),
            (
                "dataset=%s" % dataset_id,
                "total_variants=%s" % _number(formatter.get("rows_in")),
                "genome_build=%s" % input_metadata.get(
                    "genome_build", module.genome_build.value,
                ),
                "embedded_sample=%s" % input_metadata.get(
                    "postgwas_dataset_id", "not available",
                ),
                "structural_validation=passed",
            ),
            (input_vcf,) if pipeline_mode else (),
        ),
        _record(
            2,
            "Validate the PLINK LD-reference files and determine the BIM variant-ID format",
            "completed",
            "%s BIM variants were validated; BIM column %s uses %s identifiers."
            % (
                _number(observation.get("variants", qc["reference_variant_count"])),
                bim_identifier_column,
                variant_id_label,
            ),
            (
                "required_PLINK_files=present and non-empty",
                "BIM_structural_validation=passed",
                "declared_build=%s" % module.genome_build.value,
                "declared_population=%s" % module.population.value,
                "build_and_population_provenance=configuration declarations; "
                "not independently verifiable from PLINK file contents",
            ),
            (preflight.ld_reference_prefix,),
        ),
        _record(
            3, "Validate the gene-location file", "completed",
            "%s primary gene records were validated%s."
            % (
                _number(primary_gene_count),
                (
                    "; %s alternate identifiers are available"
                    % _number(alternate_gene_count)
                    if alternate_gene_count is not None else ""
                ),
            ),
            (
                "rows=%s" % _number(
                    resolution.get("location_rows", gene_scope.get("input_units"))
                ),
                "biological_context=%s" % (
                    module.mapping.definitions[primary].context or "not provided"
                ),
                "primary_column_1_unique_IDs=%s" % _number(primary_gene_count),
                "alternate_unique_IDs=%s" % _number(alternate_gene_count),
                "duplicated_alternate_ids=%s" % _number(
                    resolution.get("location_ambiguous_alternate_ids")
                ),
                "file_structure_validation=passed",
            ),
            (plan.get("gene_reference_file"),),
        ),
        _record(
            4, "Validate the original pathway GMT file",
            "completed" if pathway_requested else "not requested",
            (
                "%s pathways containing %s unique source gene identifiers were validated."
                % (
                    _number((plan.get("input_metadata") or {}).get("gene_sets")),
                    _number(resolution.get("input_unique_ids")),
                )
                if pathway_requested else "No pathway file was requested."
            ),
            (
                "format=%s" % (plan.get("input_metadata") or {}).get(
                    "detected_format", "not available"
                ),
                "pathways=%s" % _number(
                    (plan.get("input_metadata") or {}).get("gene_sets")
                ),
                "unique_gene_identifiers=%s" % _number(
                    resolution.get("input_unique_ids")
                ),
                "pathway_identifiers_modified=false",
                "file_structure_validation=passed",
            ) if pathway_requested else (),
            (plan.get("pathway_file"),) if pathway_requested else (),
        ),
        _record(
            5,
            "Compare pathway identifiers with gene-location columns and select the "
            "MAGMA gene ID",
            "completed" if pathway_ready else (
                "skipped" if pathway_requested else "not requested"
            ),
            (
                "%s; final overlap %s/%s (%s)."
                % (
                    identifier_decision,
                    _number(validation.get("overlapping_unique_ids")),
                    _number(validation.get("comparison_unique_ids")),
                    _percent(
                        validation.get("overlapping_unique_ids"),
                        validation.get("comparison_unique_ids"),
                    ),
                )
                if pathway_requested else "Pathway identifier alignment was not required."
            ),
            (
                "primary_ID_overlap=%s/%s (%s)" % (
                    _number(resolution.get("direct_primary_matches")),
                    _number(resolution.get("primary_comparison_unique_ids")),
                    _percent(
                        resolution.get("direct_primary_matches"),
                        resolution.get("primary_comparison_unique_ids"),
                    ),
                ),
                "alternate_ID_overlap=%s/%s (%s)" % (
                    _number(resolution.get("alternate_candidate_matches")),
                    _number(resolution.get("location_alternate_unique_ids")),
                    _percent(
                        resolution.get("alternate_candidate_matches"),
                        resolution.get("location_alternate_unique_ids"),
                    ),
                ),
                "ambiguous_source_ids=%s" % _number(
                    resolution.get("ambiguous_input_alternate_ids")
                ),
                "unmatched_source_ids=%s" % _number(
                    resolution.get("unmatched_input_ids")
                ),
                "required_overlap=%.2f%%" % (
                    float(plan.get("minimum_overlap", 0.0)) * 100
                ),
                "pathway_analysis=%s" % (
                    "available" if pathway_ready else "unavailable; gene analysis continues"
                ),
                "pathway_identifiers=unchanged",
            ) if pathway_requested else (),
        ),
        _record(
            6, "Prepare variant-level inputs for MAGMA", "completed",
            "%s variants were retained and written to both final MAGMA input files."
            % _number(qc["retained_rows"]),
            tuple(
                "%s=%s" % (field[1].replace(" ", "_"), field[2])
                for field in variant_input_fields if len(field) == 3
            ),
            (
                result["variant_preparation"].get("excluded_variants"),
                *prepared_outputs,
            ),
        ),
        _record(
            7, "Create the MAGMA SNP-to-gene annotation", "completed",
            "%s annotation units were retained for %s."
            % (
                _number(gene_scope.get("retained_units")),
                mapping["display_name"],
            ),
            (
                "mapping_method=%s" % mapping["mapping_method"],
                "gene_identifier_source=%s" % mapping.get(
                    "effective_gene_id_source", "mapping_definition"
                ),
                "window=%s kb upstream,%s kb downstream" % (
                    module.gene_window_upstream_kb,
                    module.gene_window_downstream_kb,
                ),
                "excluded_units=%s" % _number(gene_scope.get("excluded_units")),
                (
                    "alternate_ID_duplicate_policy=%s"
                    % pathway_location["duplicate_policy"]
                    if pathway_location else ""
                ),
                (
                    "ambiguous_alternate_IDs=%s"
                    % _number(pathway_location["ambiguous_alternate_ids"])
                    if pathway_location else ""
                ),
                (
                    "missing_alternate_ID_rows=%s"
                    % _number(pathway_location["missing_alternate_id_rows"])
                    if pathway_location else ""
                ),
            ),
            (
                mapping.get("pathway_compatible_gene_location"),
                mapping.get("gene_annotation"),
            ),
        ),
        _record(
            8, "Calculate MAGMA gene-association statistics", "completed",
            "%s annotation units were analysed in %s batch(es) using %s worker(s)."
            % (
                _number(mapping["batching"]["annotated_units"]),
                _number(mapping["batching"]["batches"]),
                _number(mapping["batching"]["workers"]),
            ),
            ("gene_model=%s" % module.gene_model,),
            (mapping.get("magma_genes_raw"), mapping.get("magma_genes_out")),
        ),
        _record(
            9, "Annotate gene results and adjust gene p-values", "completed",
            "%s genes were tested; %s were nominally significant."
            % (
                _number(gene_significance.get("tested")),
                _number(gene_significance.get("nominal_significant")),
            ),
            (
                "Bonferroni_significant=%s" % _number(
                    gene_adjusted.get("bonferroni")
                ),
                "BH-FDR_significant=%s" % _number(
                    gene_adjusted.get("fdr_bh")
                ),
            ),
            (mapping.get("magma_gene_results"),),
        ),
        _record(
            10, "Run competitive pathway analysis with the original pathway file",
            "completed" if pathway_completed else (
                "skipped" if pathway_requested else "not requested"
            ),
            (
                "The original pathway file containing %s pathways was passed "
                "directly to MAGMA without modification."
                % _number((plan.get("input_metadata") or {}).get("gene_sets"))
                if pathway_completed else
                "Pathway analysis was not performed because pathway and "
                "gene-location identifiers were incompatible."
                if pathway_requested else "Competitive pathway analysis was not requested."
            ),
            (
                "pathway_file_modified=false",
                "tested_genes=%s" % _number(
                    (mapping.get("tested_gene_coverage") or {}).get(
                        "reference_unique_ids"
                    )
                ),
                "pathway_genes_represented=%s" % _number(
                    (mapping.get("tested_gene_coverage") or {}).get(
                        "overlapping_unique_ids"
                    )
                ),
            ) if pathway_completed else (),
            (
                plan.get("pathway_file"), mapping.get("magma_gene_sets_raw"),
            ) if pathway_completed else (),
        ),
        _record(
            11, "Annotate pathway results and adjust pathway p-values",
            "completed" if pathway_completed else (
                "skipped" if pathway_requested else "not requested"
            ),
            (
                "%s pathways were tested; %s were significant by the "
                "configured primary correction."
                % (
                    _number(pathway_significance.get("tested")),
                    _number(pathway_significance.get("primary_significant")),
                )
                if pathway_completed else
                "No pathway result table was produced."
            ),
            (
                "primary_correction=%s" % module.multiple_testing.primary_method,
                "nominally_significant=%s" % _number(
                    pathway_significance.get("nominal_significant")
                ),
            ) if pathway_completed else (),
            (mapping.get("magma_pathway"),) if pathway_completed else (),
        ),
        _record(
            12, "Validate and publish all outputs", "completed",
            (
                "Gene and pathway analyses completed and all declared outputs passed validation."
                if pathway_completed else
                "Gene analysis completed; pathway analysis was not requested or was unavailable."
            ),
            ("selected_mappings=%s" % len(result["mapping_analyses"]),),
            (result.get("magma_mapping_comparison"),),
        ),
    ]
    if include_pathway_stages:
        return records
    pathway_steps = {4, 5, 10, 11}
    gene_only_records = [
        record for record in records if record["step"] not in pathway_steps
    ]
    for step, record in enumerate(gene_only_records, 1):
        record["step"] = step
    gene_only_records[-1]["summary"] = (
        "Gene-association analysis completed and all declared outputs passed "
        "validation."
    )
    return gene_only_records


def render_magma_results_screen(
    result: Mapping[str, Any],
    configuration,
    *,
    dataset_id: str,
    output_directory: str | Path,
    log_file: str | Path,
    label_width: int,
    resumed: bool = False,
    include_pathway_results: bool = True,
) -> str:
    """Render concise scientific findings after the detailed live stages."""
    module = configuration.modules.magma
    threshold = module.multiple_testing.reporting_significance_threshold
    method_labels = module.multiple_testing.reporting_method_labels
    analyses = result["mapping_analyses"]
    common_methods = [
        method for method in module.multiple_testing.gene_methods
        if method in module.multiple_testing.global_methods
    ]
    primary_method = module.multiple_testing.primary_method
    pathway_methods = list(dict.fromkeys((*common_methods, primary_method)))
    screen_summary = module.screen_summary

    lines = [
        "",
        screen_line("analysis", "MAGMA analysis results", indent=2),
        screen_field(
            "info", "Dataset", dataset_id, indent=6, label_width=label_width,
        ),
        screen_field(
            "analysis", "Reporting threshold", "%g" % threshold,
            indent=6, label_width=label_width,
        ),
    ]
    if resumed:
        lines.append(screen_field(
            "success", "Checkpoint reuse", "validated completed results",
            indent=6, label_width=label_width,
        ))

    multiple_analyses = len(analyses) > 1
    for name, analysis in analyses.items():
        if multiple_analyses:
            lines.extend(("", screen_line(
                "analysis", analysis.get("display_name", name), indent=6,
            )))
        section_indent = 10 if multiple_analyses else 6
        field_indent = section_indent + 4
        gene_significance = analysis.get("gene_significance") or {}
        gene_adjusted = gene_significance.get("adjusted_significant") or {}
        lines.extend((
            "",
            screen_line("genetic", "Gene association", indent=section_indent),
            screen_field(
                "count", "Genes tested", _number(gene_significance.get("tested")),
                indent=field_indent, label_width=label_width,
            ),
            screen_field(
                "info", "Nominally significant · P ≤ %g" % threshold,
                _number(gene_significance.get("nominal_significant")),
                indent=field_indent, label_width=label_width,
            ),
        ))
        for method in module.multiple_testing.gene_methods:
            lines.append(screen_field(
                "success",
                "%s significant · adjusted P ≤ %g"
                % (method_labels[method], threshold),
                _number(gene_adjusted.get(method)),
                indent=field_indent,
                label_width=label_width,
            ))
        gene_name_columns = (
            (module.result_schema.gene_id_column,)
            if analysis.get("gene_id_type") == "symbol" else
            (
                module.result_schema.gene_reference_alternate_id_column,
                module.result_schema.gene_id_column,
            )
        )
        gene_rows = _top_associations(
            analysis.get("magma_gene_results"),
            name_columns=gene_name_columns,
            p_value_column=module.result_schema.gene_p_value_column,
            methods=module.multiple_testing.gene_methods,
            correction_column_pattern=(
                module.result_schema.global_correction_column_pattern
            ),
            delimiter=module.result_schema.report_delimiter,
            null_value=module.result_schema.report_null_value,
            limit=screen_summary.top_gene_rows,
        )
        lines.extend(_ranked_association_lines(
            gene_rows,
            heading="Top %s gene associations · ranked by unadjusted P"
            % screen_summary.top_gene_rows,
            method_labels=method_labels,
            p_value_significant_digits=(
                screen_summary.p_value_significant_digits
            ),
            indent=field_indent,
        ))

        pathway_status = (analysis.get("gene_set_analysis") or {}).get("status")
        if include_pathway_results:
            lines.extend(("", screen_line(
                "analysis", "Competitive pathway analysis", indent=section_indent,
            )))
        if include_pathway_results and pathway_status == "completed":
            coverage = analysis.get("tested_gene_coverage") or {}
            pathway_significance = analysis.get("pathway_significance") or {}
            pathway_adjusted = (
                pathway_significance.get("adjusted_significant") or {}
            )
            supplied = coverage.get("total_gene_sets")
            represented = coverage.get(
                "gene_sets_with_at_least_one_reference_gene"
            )
            without_tested_genes = (
                int(supplied) - int(represented)
                if supplied is not None and represented is not None else None
            )
            singleton_pathways = coverage.get(
                "gene_sets_with_exactly_one_reference_gene"
            )
            tested_pathways = pathway_significance.get("tested")
            residual_pathways = (
                int(represented) - int(tested_pathways)
                if singleton_pathways is None
                and represented is not None
                and tested_pathways is not None else None
            )
            lines.extend((
                screen_field(
                    "count", "Pathways supplied", _number(supplied),
                    indent=field_indent, label_width=label_width,
                ),
                screen_field(
                    "loss" if without_tested_genes else "success",
                    "Pathways without a tested gene", _number(without_tested_genes),
                    indent=field_indent, label_width=label_width,
                ),
                *(
                    (screen_field(
                        "loss" if singleton_pathways else "success",
                        "Pathways with only one tested gene",
                        _number(singleton_pathways),
                        indent=field_indent, label_width=label_width,
                    ),)
                    if singleton_pathways is not None else ()
                ),
                *(
                    (screen_field(
                        "loss" if residual_pathways else "success",
                        "Pathways with tested genes but no MAGMA result",
                        _number(residual_pathways),
                        indent=field_indent, label_width=label_width,
                    ),)
                    if residual_pathways is not None else ()
                ),
                screen_field(
                    "count", "Pathways tested",
                    _number(tested_pathways),
                    indent=field_indent, label_width=label_width,
                ),
                screen_field(
                    "info", "Nominally significant · P ≤ %g" % threshold,
                    _number(pathway_significance.get("nominal_significant")),
                    indent=field_indent, label_width=label_width,
                ),
            ))
            for method in pathway_methods:
                lines.append(screen_field(
                    "success",
                    "%s significant · adjusted P ≤ %g"
                    % (method_labels[method], threshold),
                    _number(pathway_adjusted.get(method)),
                    indent=field_indent,
                    label_width=label_width,
                ))
            pathway_rows = _top_associations(
                analysis.get("magma_pathway"),
                name_columns=(
                    module.result_schema.gene_set_full_name_column,
                    module.result_schema.gene_set_name_column,
                ),
                p_value_column=module.result_schema.gene_set_p_value_column,
                methods=pathway_methods,
                correction_column_pattern=(
                    module.result_schema.global_correction_column_pattern
                ),
                delimiter=module.result_schema.report_delimiter,
                null_value=module.result_schema.report_null_value,
                limit=screen_summary.top_pathway_rows,
            )
            lines.extend(_ranked_association_lines(
                pathway_rows,
                heading="Top %s pathways · ranked by unadjusted P"
                % screen_summary.top_pathway_rows,
                method_labels=method_labels,
                p_value_significant_digits=(
                    screen_summary.p_value_significant_digits
                ),
                indent=field_indent,
            ))
        elif include_pathway_results:
            reason = (analysis.get("gene_set_analysis") or {}).get("reason")
            if pathway_status == "skipped" and reason:
                reason = (
                    "not performed because gene identifiers were incompatible: %s"
                    % reason
                )
            lines.append(screen_field(
                "warning" if pathway_status == "skipped" else "info",
                "Pathway analysis",
                reason or str(pathway_status or "not requested").replace("_", " "),
                indent=field_indent,
                label_width=label_width,
            ))

    lines.extend((
        "",
        screen_line("info", "Result files", indent=6),
    ))
    for name, analysis in analyses.items():
        suffix = (
            " · %s" % analysis.get("display_name", name)
            if multiple_analyses else ""
        )
        lines.append(screen_field(
            "genetic", "Gene results%s" % suffix,
            Path(analysis["magma_gene_results"]).name,
            indent=10, label_width=label_width,
        ))
        if include_pathway_results and analysis.get("magma_pathway"):
            lines.append(screen_field(
                "analysis", "Pathway results%s" % suffix,
                Path(analysis["magma_pathway"]).name,
                indent=10, label_width=label_width,
            ))
    lines.extend((
        screen_field(
            "info", "Mapping comparison",
            Path(result["magma_mapping_comparison"]).name,
            indent=10, label_width=label_width,
        ),
        screen_field(
            "success", "HTML report",
            Path(result["magma_pipeline_summary_html"]).name,
            indent=10, label_width=label_width,
        ),
        screen_field(
            "success", "Analysis summary CSV",
            Path(result["magma_pipeline_summary_csv"]).name,
            indent=10, label_width=label_width,
        ),
        screen_field(
            "info", "Full analysis log", Path(log_file).name,
            indent=10, label_width=label_width,
        ),
        screen_field(
            "info", "Results folder", Path(output_directory).name,
            indent=10, label_width=label_width,
        ),
        "",
        screen_line("success", "MAGMA analysis completed successfully", indent=2),
        "",
    ))
    return "\n".join(lines)


def _top_associations(
    source: str | Path | None,
    *,
    name_columns: Sequence[str],
    p_value_column: str,
    methods: Sequence[str],
    correction_column_pattern: str,
    delimiter: str,
    null_value: str,
    limit: int,
) -> list[dict[str, Any]]:
    """Read the strongest validated associations from a corrected result table."""
    correction_columns = [
        (method, correction_column_pattern.format(method=method))
        for method in methods
    ]
    rows, columns = _result_rows(
        source,
        (*name_columns, p_value_column, *(
            column for _, column in correction_columns
        )),
        p_value_column,
        delimiter,
        null_value,
        limit=limit,
    )
    names = [column for column in name_columns if column in columns]
    if not rows or not names:
        return []
    associations = []
    for values in rows:
        row = dict(zip(columns, values))
        label = next(
            (
                str(row[column]) for column in names
                if row.get(column) not in {None, "", null_value}
            ),
            "not available",
        )
        associations.append({
            "name": label,
            "p_value": row[p_value_column],
            "adjusted": {
                method: row[column]
                for method, column in correction_columns if column in row
            },
        })
    return associations


def _ranked_association_lines(
    rows: Sequence[Mapping[str, Any]],
    *,
    heading: str,
    method_labels: Mapping[str, str],
    p_value_significant_digits: int,
    indent: int,
) -> list[str]:
    """Render a compact ranked list using the shared aligned screen fields."""
    if not rows:
        return []
    lines = ["", screen_line("analysis", heading, indent=indent)]
    rank_width = len(str(len(rows)))
    p_value_pattern = "%%.%dg" % p_value_significant_digits
    for rank, row in enumerate(rows, 1):
        values = [
            str(row["name"]),
            "P = %s"
            % format_number(row["p_value"], p_value_pattern, missing="NA"),
        ]
        values.extend(
            "%s = %s"
            % (
                method_labels[method],
                format_number(value, p_value_pattern, missing="NA"),
            )
            for method, value in row["adjusted"].items()
        )
        lines.append(screen_field(
            "info", str(rank), " · ".join(values),
            indent=indent + 4, label_width=rank_width,
        ))
    return lines


def _configured_records(records, schema) -> list[dict[str, Any]]:
    return [
        {
            schema.step_column: record["step"],
            schema.stage_column: record["stage"],
            schema.status_column: record["status"],
            schema.summary_column: record["summary"],
            schema.details_column: record["details"],
            schema.output_column: record["output"],
        }
        for record in records
    ]


def write_magma_pipeline_csv(records, destination, schema) -> Path:
    fieldnames = [
        schema.step_column,
        schema.stage_column,
        schema.status_column,
        schema.summary_column,
        schema.details_column,
        schema.output_column,
    ]
    return write_delimited_report(
        _configured_records(records, schema),
        destination,
        fieldnames=fieldnames,
        delimiter=schema.csv_delimiter,
        null_value=schema.null_value,
    )


def read_magma_pipeline_csv(source, schema) -> list[dict[str, Any]]:
    with Path(source).open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter=schema.csv_delimiter))
    return [
        {
            "step": int(row[schema.step_column]),
            "stage": row[schema.stage_column],
            "status": row[schema.status_column],
            "summary": row[schema.summary_column],
            "details": row[schema.details_column],
            "output": row[schema.output_column],
        }
        for row in rows
    ]


def _report_value(value: Any) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        return "NA" if value != value else f"{value:.6g}"
    return str(value)


def _report_table(rows: Sequence[Mapping[str, Any]], columns: Sequence[str]) -> str:
    if not rows or not columns:
        return '<p class="muted">No result rows are available.</p>'
    header = "".join("<th>%s</th>" % escape(column) for column in columns)
    body = "".join(
        "<tr>%s</tr>" % "".join(
            "<td>%s</td>" % escape(_report_value(row.get(column)))
            for column in columns
        )
        for row in rows
    )
    return (
        '<div class="table-wrap"><table><thead><tr>%s</tr></thead>'
        '<tbody>%s</tbody></table></div>' % (header, body)
    )


def _result_rows(
    source: str | Path | None,
    configured_columns: Sequence[str],
    p_value_column: str,
    delimiter: str,
    null_value: str,
    *,
    limit: int | None = None,
) -> tuple[list[list[str]], list[str]]:
    if source is None or not Path(source).is_file():
        return [], []
    scan = pl.scan_csv(
        source,
        separator=delimiter,
        null_values=[null_value],
        truncate_ragged_lines=False,
    )
    available = scan.columns
    columns = [column for column in configured_columns if column in available]
    if p_value_column not in available:
        return [], columns
    if p_value_column not in columns:
        columns.append(p_value_column)
    ranked = (
        scan.select(columns)
        .with_columns(
            pl.col(p_value_column).cast(pl.Float64, strict=False).alias("__rank")
        )
        .sort("__rank", nulls_last=True)
        .drop("__rank")
    )
    if limit is not None:
        ranked = ranked.filter(
            pl.col(p_value_column).cast(pl.Float64, strict=False).is_finite()
        ).head(limit)
    frame = ranked.collect()
    return [
        [_report_value(value) for value in row]
        for row in frame.iter_rows()
    ], columns


def _paginated_result_table(
    rows: Sequence[Sequence[str]],
    columns: Sequence[str],
    *,
    table_id: str,
    page_size: int,
) -> str:
    if not rows or not columns:
        return '<p class="muted">No result rows are available.</p>'
    header = "".join(
        '<th><button type="button" data-column="%s">%s</button></th>'
        % (index, escape(column))
        for index, column in enumerate(columns)
    )
    payload = json.dumps(rows, ensure_ascii=False, separators=(",", ":")).replace(
        "</", "<\\/"
    )
    return """
<div class="paginated-results" id="%s" data-page-size="%s" data-source="%s-data">
  <div class="table-controls">
    <label>Search all rows <input type="search" class="result-search" autocomplete="off"></label>
    <div class="page-controls">
      <button type="button" class="previous-page">Previous</button>
      <span class="page-status" aria-live="polite"></span>
      <button type="button" class="next-page">Next</button>
    </div>
  </div>
  <div class="table-wrap"><table><thead><tr>%s</tr></thead><tbody></tbody></table></div>
  <noscript><p class="notice">Enable JavaScript to browse this table. The complete TSV result is available from the link below.</p></noscript>
</div>
<script type="application/json" id="%s-data">%s</script>""" % (
        escape(table_id, quote=True), page_size, escape(table_id, quote=True),
        header, escape(table_id, quote=True), payload,
    )


def _metric(label: str, value: Any, tone: str = "") -> str:
    return (
        '<div class="metric %s"><span>%s</span><strong>%s</strong></div>'
        % (escape(tone), escape(label), escape(_report_value(value)))
    )


def _file_link(path: str | Path | None, report_path: Path) -> str:
    if not path:
        return "not produced"
    source = Path(path)
    href = os.path.relpath(source, report_path.parent)
    return '<a href="%s">%s</a>' % (
        escape(href, quote=True), escape(source.name),
    )


def write_magma_pipeline_html(
    records: Sequence[Mapping[str, Any]],
    destination,
    *,
    dataset_id: str,
    schema,
    result: Mapping[str, Any],
    configuration,
    include_pathway_results: bool = True,
) -> Path:
    """Write the scientific MAGMA report from validated analysis outputs."""
    destination = Path(destination)
    module = configuration.modules.magma
    report_config = module.html_report
    qc = result["variant_preparation"]["qc"]
    analyses = result["mapping_analyses"]
    primary = analyses[result["primary_mapping"]]
    primary_genes = primary.get("gene_significance") or {}
    primary_pathways = primary.get("pathway_significance") or {}
    gene_adjusted = primary_genes.get("adjusted_significant") or {}
    threshold = module.multiple_testing.reporting_significance_threshold
    primary_method_label = module.multiple_testing.reporting_method_labels[
        module.multiple_testing.primary_method
    ]

    overview_metrics = [
        _metric("GWAS variants", qc.get("input_rows")),
        _metric(
            "Variants represented in LD reference",
            "%s (%s)" % (
                _number(qc.get("reference_unique_id_matches")),
                _percent(
                    qc.get("reference_unique_id_matches"),
                    qc.get("input_unique_variants"),
                ),
            ),
        ),
        _metric("Genes tested", primary_genes.get("tested")),
        _metric("Genes significant by BH-FDR", gene_adjusted.get("fdr_bh"), "good"),
    ]
    if include_pathway_results:
        overview_metrics.extend((
            _metric("Pathways tested", primary_pathways.get("tested")),
            _metric(
                "Pathways significant by %s" % primary_method_label,
                primary_pathways.get("primary_significant"),
                "good",
            ),
        ))
    overview = "".join(overview_metrics)

    scope = primary.get("gene_scope") or {}
    coverage = primary.get("tested_gene_coverage") or {}
    qc_rows = [
        {"Assessment": "GWAS variant records", "Result": qc.get("input_rows")},
        {
            "Assessment": "Unique GWAS variants",
            "Result": qc.get("input_unique_variants"),
        },
        {
            "Assessment": "Variants represented in PLINK BIM",
            "Result": "%s (%s)" % (
                _number(qc.get("reference_unique_id_matches")),
                _percent(
                    qc.get("reference_unique_id_matches"),
                    qc.get("input_unique_variants"),
                ),
            ),
        },
        {
            "Assessment": "Variants absent from PLINK BIM",
            "Result": "%s (%s)" % (
                _number(qc.get("not_in_reference_rows")),
                _percent(
                    qc.get("not_in_reference_rows"),
                    qc.get("input_unique_variants"),
                ),
            ),
        },
        {
            "Assessment": "BIM intersection",
            "Result": "applied" if qc.get("reference_intersection_enabled") else "not applied",
        },
        {"Assessment": "Final variant-level input rows", "Result": qc.get("retained_rows")},
        {"Assessment": "Duplicate rows removed", "Result": qc.get("duplicate_rows_removed")},
        {
            "Assessment": "Chromosome-policy exclusions",
            "Result": qc.get("excluded_chromosome_rows"),
        },
        {"Assessment": "MHC variant exclusions", "Result": qc.get("excluded_mhc_rows")},
        {
            "Assessment": "Annotated genes before scope exclusions",
            "Result": scope.get("input_units"),
        },
        {"Assessment": "Annotated genes retained", "Result": scope.get("retained_units")},
        {"Assessment": "Annotated genes excluded", "Result": scope.get("excluded_units")},
        {"Assessment": "Genes with valid MAGMA results", "Result": primary_genes.get("tested")},
        {
            "Assessment": "Annotated genes without valid MAGMA results",
            "Result": (
                int(scope["retained_units"]) - int(primary_genes["tested"])
                if scope.get("retained_units") is not None
                and primary_genes.get("tested") is not None else None
            ),
        },
    ]
    if include_pathway_results:
        qc_rows.append({
            "Assessment": "Tested genes represented in pathway file",
            "Result": "%s/%s (%s)" % (
                _number(coverage.get("overlapping_unique_ids")),
                _number(coverage.get("reference_unique_ids")),
                _percent(
                    coverage.get("overlapping_unique_ids"),
                    coverage.get("reference_unique_ids"),
                ),
            ) if coverage else "not assessed",
        })

    analysis_sections = []
    for analysis_index, (name, analysis) in enumerate(analyses.items(), start=1):
        gene_stats = analysis.get("gene_significance") or {}
        pathway_stats = analysis.get("pathway_significance") or {}
        gene_rows, gene_columns = _result_rows(
            analysis.get("magma_gene_results"),
            report_config.gene_columns,
            module.result_schema.gene_p_value_column,
            module.result_schema.report_delimiter,
            module.result_schema.report_null_value,
        )
        pathway_rows, pathway_columns = ([], [])
        if include_pathway_results:
            pathway_rows, pathway_columns = _result_rows(
                analysis.get("magma_gene_sets_corrected"),
                report_config.pathway_columns,
                module.result_schema.gene_set_p_value_column,
                module.result_schema.report_delimiter,
                module.result_schema.report_null_value,
            )
        gene_counts = gene_stats.get("adjusted_significant") or {}
        pathway_counts = pathway_stats.get("adjusted_significant") or {}
        pathway_status = (analysis.get("gene_set_analysis") or {}).get("status")
        pathway_metrics = ""
        pathway_section = ""
        pathway_download = ""
        if include_pathway_results:
            pathway_metrics = "".join((
                _metric("Pathways tested", pathway_stats.get("tested")),
                _metric(
                    "Nominal pathway associations",
                    pathway_stats.get("nominal_significant"),
                ),
                _metric(
                    "Bonferroni-significant pathways",
                    pathway_counts.get("bonferroni"),
                    "good",
                ),
                _metric(
                    "BH-FDR-significant pathways",
                    pathway_counts.get("fdr_bh"),
                    "good",
                ),
            ))
            pathway_section = """
<h3>Competitive pathway results</h3>
%s
%s""" % (
                (
                    '<p class="muted">All %s pathway results are included and initially '
                    'ranked by the unadjusted competitive-test p-value. A positive MAGMA '
                    'BETA indicates stronger association among genes in the pathway '
                    'relative to other tested genes.</p>' % _number(len(pathway_rows))
                    if pathway_status == "completed" else
                    '<p class="notice">Pathway analysis was not completed: %s</p>'
                    % escape(str(
                        (analysis.get("gene_set_analysis") or {}).get("reason")
                        or pathway_status
                    ))
                ),
                _paginated_result_table(
                    pathway_rows,
                    pathway_columns,
                    table_id="pathway-results-%s" % analysis_index,
                    page_size=report_config.page_size,
                ),
            )
            pathway_download = "; pathway results %s" % _file_link(
                analysis.get("magma_pathway"), destination,
            )
        analysis_sections.append("""
<section><h2>%s</h2>
<p class="lead">%s</p>
<div class="metrics">%s</div>
<h3>Most strongly associated genes</h3>
<p class="muted">All %s gene results are included and initially ranked by the unadjusted MAGMA gene p-value. Search, sort, or move between pages without changing the result set.</p>
%s
%s
<div class="downloads"><strong>Complete result tables:</strong> gene results %s%s.</div>
</section>""" % (
            escape(analysis.get("display_name", name)),
            escape(analysis.get("result_statistic_interpretation", "MAGMA association result")),
            "".join((
                _metric("Genes tested", gene_stats.get("tested")),
                _metric("Nominal gene associations", gene_stats.get("nominal_significant")),
                _metric("Bonferroni-significant genes", gene_counts.get("bonferroni"), "good"),
                _metric("BH-FDR-significant genes", gene_counts.get("fdr_bh"), "good"),
            )) + pathway_metrics,
            _number(len(gene_rows)),
            _paginated_result_table(
                gene_rows,
                gene_columns,
                table_id="gene-results-%s" % analysis_index,
                page_size=report_config.page_size,
            ),
            pathway_section,
            _file_link(analysis.get("magma_gene_results"), destination),
            pathway_download,
        ))

    workflow_rows = "".join(
        "<tr><td>%s</td><td>%s</td><td><span class='status %s'>%s</span></td><td>%s</td></tr>"
        % (
            escape(str(record["step"])),
            escape(str(record["stage"])),
            escape(str(record["status"]).replace(" ", "-")),
            escape(str(record["status"])),
            escape(str(record["summary"])),
        )
        for record in records
    )
    parameters = [
        {"Parameter": "Genome build", "Value": module.genome_build.value},
        {"Parameter": "LD-reference population", "Value": module.population.value},
        {"Parameter": "Gene association model", "Value": module.gene_model},
        {"Parameter": "Gene window upstream (kb)", "Value": module.gene_window_upstream_kb},
        {"Parameter": "Gene window downstream (kb)", "Value": module.gene_window_downstream_kb},
        {"Parameter": "MHC policy", "Value": module.mhc.policy},
        {
            "Parameter": "Excluded chromosomes",
            "Value": ", ".join(module.chromosomes.exclude) or "none",
        },
        {"Parameter": "Reporting threshold", "Value": "p ≤ %s" % threshold},
        {"Parameter": "MAGMA version", "Value": result.get("magma_version")},
    ]
    if include_pathway_results:
        parameters.insert(-1, {
            "Parameter": "Primary pathway correction",
            "Value": primary_method_label,
        })
    report_name = (
        "MAGMA gene and pathway association"
        if include_pathway_results else "MAGMA gene association"
    )
    interpretation = (
        "Gene results test aggregate SNP association within configured gene "
        "regions while accounting for linkage disequilibrium. Competitive "
        "pathway results test whether genes in a pathway are more strongly "
        "associated than other tested genes. Association does not establish "
        "causality."
        if include_pathway_results else
        "Gene results test aggregate SNP association within configured gene "
        "regions while accounting for linkage disequilibrium. Association does "
        "not establish causality."
    )
    document = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>%s report · %s</title>
<style>
:root{--ink:#172033;--muted:#607084;--line:#dbe2ea;--blue:#075985;--pale:#eef6fb;--green:#08783e;--amber:#9a6700}
*{box-sizing:border-box}body{font-family:Inter,ui-sans-serif,system-ui,sans-serif;margin:0;color:var(--ink);background:#f3f6fa;line-height:1.45}
main{max-width:1500px;margin:auto;padding:2rem}header,section{background:white;padding:1.5rem 1.7rem;margin-bottom:1.2rem;border:1px solid var(--line);border-radius:14px;box-shadow:0 2px 10px #1f293708}
h1{margin:.1rem 0}.subtitle,.muted{color:var(--muted)}.lead{max-width:95ch}.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:.8rem;margin:1rem 0}
.metric{border-left:4px solid var(--blue);background:var(--pale);padding:.8rem 1rem;border-radius:7px}.metric span{display:block;color:var(--muted);font-size:.82rem}.metric strong{font-size:1.2rem}.metric.good{border-color:var(--green)}
.table-wrap{overflow:auto;border:1px solid var(--line);border-radius:8px;margin:.8rem 0 1.2rem}table{border-collapse:collapse;width:100%%;font-size:.88rem;white-space:nowrap}th,td{border-bottom:1px solid var(--line);padding:.55rem .7rem;text-align:left}th{background:var(--pale);position:sticky;top:0}.status{font-weight:700}.completed{color:var(--green)}.skipped{color:var(--amber)}.not-requested{color:var(--muted)}
.table-controls{display:flex;align-items:center;justify-content:space-between;gap:1rem;flex-wrap:wrap;margin:.8rem 0}.table-controls label{font-weight:700}.table-controls input{margin-left:.45rem;min-width:18rem;padding:.48rem .6rem;border:1px solid var(--line);border-radius:6px}.page-controls{display:flex;align-items:center;gap:.65rem}.page-controls button,th button{font:inherit}.page-controls button{padding:.42rem .7rem;border:1px solid var(--line);border-radius:6px;background:white;color:var(--blue);cursor:pointer}.page-controls button:disabled{color:var(--muted);cursor:not-allowed}th button{border:0;background:transparent;color:inherit;font-weight:700;padding:0;cursor:pointer}.page-status{min-width:15rem;text-align:center;color:var(--muted)}
.notice{border-left:4px solid var(--amber);background:#fff8e6;padding:.8rem 1rem}.downloads{background:#f8fafc;padding:.8rem 1rem;border-radius:8px}a{color:var(--blue)}details summary{cursor:pointer;font-weight:700}
@media(max-width:700px){main{padding:.7rem}header,section{padding:1rem}table{font-size:.8rem}.table-controls input{min-width:12rem}.page-status{min-width:auto}}
</style></head><body><main>
<header><p class="subtitle">PostGWAS scientific results report</p><h1>%s</h1><p><strong>Dataset:</strong> %s</p>
<p class="lead">%s</p></header>
<section><h2>Results at a glance</h2><div class="metrics">%s</div>
<p class="muted">Counts use the configured reporting threshold of p ≤ %s. Adjusted counts are based on the named multiple-testing procedure, not the unadjusted p-value.</p></section>
<section><h2>Quality control and analysis coverage</h2>%s</section>
%s
<section><h2>Analysis parameters and provenance</h2>%s
<p class="muted">The PLINK genome build and population are configuration declarations; they cannot be independently established from BED/BIM/FAM contents alone.</p></section>
<section><details><summary>Complete workflow and validation summary</summary><div class="table-wrap"><table><thead><tr><th>%s</th><th>%s</th><th>%s</th><th>%s</th></tr></thead><tbody>%s</tbody></table></div></details></section>
</main>
<script>
document.querySelectorAll('.paginated-results').forEach(function(container) {
  const rows = JSON.parse(document.getElementById(container.dataset.source).textContent);
  const pageSize = Number(container.dataset.pageSize);
  const body = container.querySelector('tbody');
  const search = container.querySelector('.result-search');
  const previous = container.querySelector('.previous-page');
  const next = container.querySelector('.next-page');
  const status = container.querySelector('.page-status');
  let filtered = rows.slice();
  let page = 0;
  let sortColumn = null;
  let ascending = true;

  function compare(left, right) {
    const leftNumber = Number(String(left).replaceAll(',', ''));
    const rightNumber = Number(String(right).replaceAll(',', ''));
    if (Number.isFinite(leftNumber) && Number.isFinite(rightNumber)) {
      return leftNumber - rightNumber;
    }
    return String(left).localeCompare(String(right), undefined, {numeric: true});
  }

  function render() {
    const pageCount = Math.max(1, Math.ceil(filtered.length / pageSize));
    page = Math.min(page, pageCount - 1);
    const start = page * pageSize;
    const end = Math.min(start + pageSize, filtered.length);
    body.replaceChildren();
    filtered.slice(start, end).forEach(function(row) {
      const tr = document.createElement('tr');
      row.forEach(function(value) {
        const td = document.createElement('td');
        td.textContent = value;
        tr.appendChild(td);
      });
      body.appendChild(tr);
    });
    status.textContent = filtered.length
      ? `${start + 1}–${end} of ${filtered.length.toLocaleString()} · page ${page + 1} of ${pageCount}`
      : '0 matching rows';
    previous.disabled = page === 0;
    next.disabled = page + 1 >= pageCount;
  }

  search.addEventListener('input', function() {
    const term = search.value.trim().toLocaleLowerCase();
    filtered = term
      ? rows.filter(function(row) {
          return row.some(function(value) {
            return String(value).toLocaleLowerCase().includes(term);
          });
        })
      : rows.slice();
    if (sortColumn !== null) {
      filtered.sort(function(left, right) {
        const order = compare(left[sortColumn], right[sortColumn]);
        return ascending ? order : -order;
      });
    }
    page = 0;
    render();
  });
  previous.addEventListener('click', function() { page -= 1; render(); });
  next.addEventListener('click', function() { page += 1; render(); });
  container.querySelectorAll('th button').forEach(function(button) {
    button.addEventListener('click', function() {
      const column = Number(button.dataset.column);
      ascending = sortColumn === column ? !ascending : true;
      sortColumn = column;
      filtered.sort(function(left, right) {
        const order = compare(left[column], right[column]);
        return ascending ? order : -order;
      });
      page = 0;
      render();
    });
  });
  render();
});
</script>
</body></html>""" % (
        escape(report_name), escape(dataset_id), escape(report_name),
        escape(dataset_id), escape(interpretation), overview,
        escape(str(threshold)),
        _report_table(qc_rows, ("Assessment", "Result")),
        "".join(analysis_sections),
        _report_table(parameters, ("Parameter", "Value")),
        escape(schema.step_column), escape(schema.stage_column),
        escape(schema.status_column), escape(schema.summary_column), workflow_rows,
    )
    return write_html_report(document, destination)


__all__ = [
    "MAGMA_GENE_ONLY_STAGE_KEYS",
    "MAGMA_PIPELINE_STAGES",
    "MAGMA_PIPELINE_STAGE_KEYS",
    "MAGMA_STAGE_TITLES",
    "build_magma_pipeline_summary",
    "magma_variant_input_outcome_fields",
    "magma_pipeline_progress_plan",
    "magma_pipeline_stage_numbers",
    "read_magma_pipeline_csv",
    "render_magma_results_screen",
    "write_magma_pipeline_csv",
    "write_magma_pipeline_html",
]
