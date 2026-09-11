"""Validated execution service for GCTA fastBAT and mBAT-combo."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
import re
import shlex
import shutil
import sys
import tempfile
from collections.abc import Mapping
from contextlib import contextmanager
from types import SimpleNamespace
from typing import Any

import polars as pl
import yaml

from postgwas.config import (
    load_configuration,
    load_run_configuration_for_module,
    write_resolved_configuration,
)
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.core.contracts import Artifact, ModuleResult
from postgwas.core.errors import ConfigurationError
from postgwas.core.gene_coordinates import read_gene_coordinates
from postgwas.core.genomic_scope import resolve_genomic_analysis_scope
from postgwas.core.io.reports import write_yaml_report
from postgwas.core.paths import configured_output_path, resolve_executable
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.core.plink import validate_plink_bundle_dimensions, validate_plink_files
from postgwas.core.preflight import (
    PipelinePreflightEvidence,
    PreflightFileIdentity,
    capture_preflight_file_identities,
    pipeline_preflight_evidence,
    require_pipeline_input_vcf,
    require_unchanged_preflight_files,
)
from postgwas.core.required_arguments import (
    RequiredAlternative,
    RequiredArgument,
    require_resolved_arguments,
)
from postgwas.core.resource_preparation import ResourcePreparationError, sha256
from postgwas.core.snp_sets import open_fastbat_set_memberships, validate_fastbat_set_list
from postgwas.core.ui import StageProgress
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.modules.gcta_gene.adapters import (
    build_gcta_command,
    require_supported_gcta,
    run_gcta_command,
)
from postgwas.modules.gcta_gene.errors import GctaGeneError
from postgwas.modules.formatting.reference_identifiers import (
    BimIdentifierRequirement,
    configure_reference_variant_identifiers,
)
from postgwas.modules.gcta_gene.pathway_sets import (
    PathwayGenePreflight,
    prepare_resource,
    scan_gmt_pathways,
    validate_gene_coordinate_reference,
    validate_pathway_gene_compatibility,
)
from postgwas.modules.gcta_gene.progress import (
    GctaResultProgress,
    gcta_native_log_path,
)
from postgwas.modules.gcta_gene.reporting import (
    build_gcta_scientific_summary,
    gcta_gene_reference_outcome_fields,
    gcta_gmt_preparation_outcome_fields,
    gcta_ld_reference_outcome_fields,
    gcta_ma_input_outcome_fields,
    gcta_set_source_outcome_fields,
    record_gcta_scientific_summary,
    render_gcta_scientific_summary,
    write_gcta_html_report,
)
from postgwas.modules.gcta_gene.results import (
    add_multiple_testing_results,
    multiple_testing_configuration,
    normalize_gcta_results,
    validate_raw_gcta_results,
)
from postgwas.modules.gcta_gene.stages import (
    complete_pipeline_stage,
    configure_pipeline_stage_callbacks,
    pipeline_stage_number,
    start_pipeline_stage,
)
_DNA_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_COMPLETION_SCHEMA_VERSION = 3
_LEGACY_COMPLETION_SCHEMA_VERSIONS = {2}


@dataclass(frozen=True)
class GctaInputPreflight:
    """Validated formatter output and BIM intersection reused by GCTA."""

    source: Path
    input_metrics: dict
    identifiers: set[str]
    reference_metrics: dict
    reference_ids: set[str]
    excluded_reference_ids: tuple[str, ...]
    reference_prefix: Path


@dataclass(frozen=True)
class GctaGenePipelineResources:
    """External GCTA resources validated before formatter input creation."""

    configuration: Any
    reference_prefix: Path
    reference_paths: tuple[Path, ...]
    executable: str
    version: str
    identifier_observation: Mapping[str, Any]
    analysis_scope: Mapping[str, Any]
    gene_list: Path | None
    gene_annotation_metrics: Mapping[str, Any] | None
    gmt: Path | None
    pathway_gene_preflight: PathwayGenePreflight | None
    set_list: Path | None
    set_source_metrics: Mapping[str, Any] | None
    file_identities: tuple[PreflightFileIdentity, ...]


def validate_gcta_reference_files(module) -> tuple[Path, tuple[Path, ...]]:
    """Validate the configured PLINK prefix and every required companion file."""
    if module.reference.prefix is None:
        raise GctaGeneError("--gcta-reference-prefix is required.")
    prefix = Path(module.reference.prefix).expanduser().resolve()
    paths = validate_plink_files(
        prefix, module.reference.required_extensions, error_type=GctaGeneError,
    )
    return prefix, tuple(paths.values())


def validate_pipeline_gcta_input(
    input_file: str | Path,
    configuration,
    *,
    pipeline_resources: GctaGenePipelineResources | None = None,
) -> GctaInputPreflight:
    """Validate one formatter-created ``.ma`` against BIM without rewriting it."""
    module = configuration.modules.gcta_gene
    if pipeline_resources is None:
        _validate_build_and_population(module, configuration)
    source = _required_path(str(input_file), "GCTA summary-statistics input")
    if pipeline_resources is not None:
        require_unchanged_preflight_files(
            pipeline_resources.file_identities,
            error_type=GctaGeneError,
            label="GCTA gene resource",
        )
        reference_prefix = pipeline_resources.reference_prefix
    else:
        reference_prefix, _ = validate_gcta_reference_files(module)
    input_metrics, _frame, _identifier_column, identifiers, alleles = (
        _validate_formatted_input(source, configuration, module.method)
    )
    analysis_scope = (
        dict(pipeline_resources.analysis_scope)
        if pipeline_resources is not None
        else _analysis_scope(configuration)
    )
    reference_metrics, reference_ids, excluded_reference_ids = _validate_reference(
        reference_prefix,
        module,
        identifiers,
        alleles,
        analysis_scope=analysis_scope,
        retain_variant_ids=module.method == "fastbat_set",
        direct_input_mode=False,
    )
    input_metrics.update({
        "analysis_input_variants": input_metrics["variants"],
        "variant_ids_replaced": 0,
        "unresolved_variants_removed": 0,
        "direct_input_unmodified": False,
        "formatter_input_unmodified": True,
    })
    reference_metrics["variant_id_match_policy"] = (
        "validate_pipeline_formatter_bim_compatibility_without_rewriting"
    )
    return GctaInputPreflight(
        source=source,
        input_metrics=input_metrics,
        identifiers=identifiers,
        reference_metrics=reference_metrics,
        reference_ids=reference_ids,
        excluded_reference_ids=excluded_reference_ids,
        reference_prefix=reference_prefix,
    )


def _resolved_configuration(args):
    has_set_list = hasattr(args, "fastbat_set_list")
    has_gmt = hasattr(args, "gmt")
    if has_set_list and has_gmt:
        raise ConfigurationError(
            "--fastbat-set-list and --gmt are mutually exclusive."
        )
    module_overrides = explicit_overrides(args, {
        "gcta_gene_method": "method",
        "gcta_input_file": "input_file",
        "genome_build": "genome_build",
        "gcta_reference_prefix": "reference.prefix",
        "gcta_reference_population": "reference.population",
        "gene_list": "gene_annotation.file",
        "fastbat_set_list": "set_annotation.file",
        "gmt": "set_annotation.gmt_file",
        "gmt_chromosome_label_policy": (
            "set_annotation.conversion.chromosome_label_policy"
        ),
        "gmt_duplicate_gene_policy": "set_annotation.conversion.duplicate_gene_policy",
        "gmt_unmapped_gene_policy": "set_annotation.conversion.unmapped_gene_policy",
        "gcta_minimum_gene_id_overlap": (
            "set_annotation.conversion.minimum_gene_id_overlap_fraction"
        ),
        "gmt_empty_pathway_policy": "set_annotation.conversion.empty_pathway_policy",
        "fastbat_oversized_set_policy": "set_annotation.oversized_set_policy",
        "gene_window_kb": "gene_window_kb",
        "fastbat_segment_size_kb": "segment_size_kb",
        "gcta_reference_maf_min": "reference_maf_min",
        "fastbat_ld_cutoff": "fastbat_ld_cutoff",
        "mbat_svd_gamma": "mbat_svd_gamma",
        "frequency_difference_max": "frequency_difference_max",
        "print_component_p_values": "print_component_p_values",
        "write_snpset": "write_snpset",
        "gcta_top_results": "reporting.top_result_count",
        "gcta_nominal_alpha": "reporting.nominal_alpha",
        "gcta_reporting_alpha": "reporting.familywise_alpha",
        "gcta_fdr_alpha": "reporting.fdr_alpha",
        "gcta_p_value_digits": "reporting.p_value_significant_digits",
        "gcta_chromosome_label_policy": (
            "variant_harmonisation.chromosome_label_policy"
        ),
        "gcta_minimum_reference_overlap": (
            "variant_harmonisation.minimum_overlap_fraction"
        ),
        "gcta_allow_strand_complement": (
            "variant_harmonisation.allow_strand_complement"
        ),
        "mhc_policy": "mhc.policy",
        "mhc_chrom": "mhc.region_override.chromosome",
        "mhc_start": "mhc.region_override.start",
        "mhc_end": "mhc.region_override.end",
        "exclude_chromosomes": "chromosomes.exclude",
    })
    set_list_override = module_overrides.pop("set_annotation.file", None)
    gmt_override = module_overrides.pop("set_annotation.gmt_file", None)
    global_overrides = explicit_overrides(args, {
        "dataset_id": "run.dataset_id",
        "output_directory": "run.output_directory",
        "threads": "execution.threads",
        "memory_gb": "execution.memory_gb",
        "seed": "execution.random_seed",
        "gcta": "resources.executables.gcta",
        "bcftools": "resources.executables.bcftools",
        "resume": "run.resume",
        "overwrite": "run.overwrite",
    })
    configuration = load_run_configuration_for_module(
        "gcta_gene",
        getattr(args, "run_config", None),
        module_overrides=module_overrides,
        global_overrides=global_overrides,
    )
    if has_set_list or has_gmt:
        module = configuration.modules.gcta_gene
        values = module.set_annotation.model_dump(mode="python")
        values["file"] = set_list_override if has_set_list else None
        values["gmt_file"] = gmt_override if has_gmt else None
        module.set_annotation = type(module.set_annotation).model_validate(values)
    return configuration


def _require_gcta_gene_arguments(
    configuration,
    input_file,
    *,
    include_generated_input: bool = True,
) -> None:
    """Validate all requirements for the resolved GCTA method together."""
    module = configuration.modules.gcta_gene
    requirements = [
        RequiredArgument(
            "--gcta-reference-prefix",
            "modules.gcta_gene.reference.prefix",
            module.reference.prefix,
        ),
        RequiredArgument(
            "--genome-build",
            "modules.gcta_gene.genome_build",
            module.genome_build,
        ),
        RequiredArgument(
            "--gcta-reference-population",
            "modules.gcta_gene.reference.population",
            module.reference.population,
        ),
    ]
    if include_generated_input:
        requirements.insert(0, RequiredArgument(
            "--gcta-input-file", "modules.gcta_gene.input_file", input_file,
        ))
    if module.method in {"fastbat_gene", "mbat_combo"} or (
        module.method == "fastbat_set"
        and module.set_annotation.gmt_file is not None
    ):
        requirements.append(RequiredArgument(
            "--gene-list",
            "modules.gcta_gene.gene_annotation.file",
            module.gene_annotation.file,
        ))
    alternatives = ()
    if module.method == "fastbat_set":
        alternatives = (
            RequiredAlternative((
                RequiredArgument(
                    "--fastbat-set-list",
                    "modules.gcta_gene.set_annotation.file",
                    module.set_annotation.file,
                ),
                RequiredArgument(
                    "--gmt",
                    "modules.gcta_gene.set_annotation.gmt_file",
                    module.set_annotation.gmt_file,
                ),
            )),
        )
    require_resolved_arguments(requirements, alternatives=alternatives)


def _required_path(value: str | None, label: str) -> Path:
    if value is None:
        raise GctaGeneError("%s is required." % label)
    path = Path(value).expanduser().resolve()
    if not path.is_file() or path.stat().st_size <= 0:
        raise GctaGeneError("%s does not exist or is empty: %s" % (label, path))
    return path


def _validate_build_and_population(module, configuration) -> None:
    if module.genome_build is None:
        raise GctaGeneError(
            "Declare --genome-build. It applies to the harmonised GWAS, PLINK LD "
            "reference, and gene list, which must all use the same genome build."
        )
    if module.genome_build not in configuration.resources.genomes:
        raise GctaGeneError("Unknown configured genome build: %s" % module.genome_build)
    population = module.reference.population
    if population is None:
        raise GctaGeneError(
            "Declare --gcta-reference-population; ancestry matching cannot be inferred."
        )
    if population not in configuration.resources.populations:
        raise GctaGeneError("Unknown configured reference population: %s" % population)


def _validate_formatted_input(
    path: Path, configuration, method: str,
) -> tuple[dict, pl.DataFrame, str, set[str], dict[str, tuple[str, str]]]:
    formatting = configuration.modules.formatting
    schema = formatting.exports["gcta_gene"].outputs["summary_statistics"]
    try:
        frame = pl.read_csv(
            path,
            separator=formatting.runtime.table_delimiter,
            null_values=formatting.runtime.input_null_values,
            infer_schema_length=configuration.modules.gcta_gene.results.infer_schema_length,
        )
    except (OSError, pl.exceptions.PolarsError) as exc:
        raise GctaGeneError("Cannot read GCTA input %s: %s" % (path, exc)) from exc
    required = [schema.columns[source] for source in schema.validation.required_columns]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise GctaGeneError(
            "GCTA %s input is missing required columns: %s"
            % (method, ", ".join(missing))
        )
    if frame.is_empty():
        raise GctaGeneError("GCTA input contains no variants: %s" % path)
    numeric_sources = set(formatting.numeric_columns)
    numeric_columns = [
        destination for source, destination in schema.columns.items()
        if source in numeric_sources
    ]
    checked = frame.with_columns(
        *[pl.col(column).cast(pl.Float64, strict=False) for column in numeric_columns]
    )
    invalid = pl.any_horizontal(pl.col(required).is_null())
    for column in numeric_columns:
        if column in required:
            invalid |= pl.col(column).is_null() | ~pl.col(column).is_finite()
    for source in schema.validation.positive_columns:
        invalid |= pl.col(schema.columns[source]) <= 0
    for source in schema.validation.nonnegative_columns:
        invalid |= pl.col(schema.columns[source]) < 0
    for source in schema.validation.open_unit_interval_columns:
        column = schema.columns[source]
        invalid |= (pl.col(column) <= 0) | (pl.col(column) >= 1)
    p_sources = [
        source for source, transform in schema.transformations.items()
        if transform == "negative_log10_to_raw_p"
    ]
    if len(p_sources) != 1:
        raise GctaGeneError(
            "GCTA formatter schema must define one raw-p transformation."
        )
    p_column = schema.columns[p_sources[0]]
    invalid |= pl.col(p_column) > 1
    invalid_rows = checked.filter(invalid.fill_null(True)).height
    if invalid_rows:
        raise GctaGeneError(
            "GCTA %s input contains %d rows with missing, non-finite, or "
            "out-of-range scientific values." % (method, invalid_rows)
        )
    canonical = formatting.canonical_columns
    identifier_column = schema.columns[canonical.resolved_variant_id]
    identifiers = [str(value).strip() for value in checked[identifier_column].to_list()]
    if any(not value for value in identifiers):
        raise GctaGeneError("GCTA input contains an empty variant identifier.")
    if len(identifiers) != len(set(identifiers)):
        raise GctaGeneError("GCTA input contains duplicate variant identifiers.")
    alleles: dict[str, tuple[str, str]] = {}
    a1_column = schema.columns[canonical.alternate_allele]
    a2_column = schema.columns[canonical.reference_allele]
    for variant, a1, a2 in checked.select(
        identifier_column, a1_column, a2_column,
    ).iter_rows():
        first, second = str(a1).strip().upper(), str(a2).strip().upper()
        if not first or not second or first == second:
            raise GctaGeneError(
                "GCTA .ma input contains invalid alleles for variant %s." % variant
            )
        alleles[str(variant)] = (first, second)
    return (
        {"variants": checked.height, "columns": checked.columns},
        frame,
        identifier_column,
        set(identifiers),
        alleles,
    )


def _normalise_chromosome(value: str, policy: str) -> str:
    chromosome = str(value).strip()
    if policy == "strip_chr_prefix" and chromosome.lower().startswith("chr"):
        chromosome = chromosome[3:]
    return chromosome


def _analysis_scope(configuration, chromosome_policy: str | None = None) -> dict:
    """Resolve the declared build-specific GCTA analysis scope once."""
    module = configuration.modules.gcta_gene
    policy = chromosome_policy or module.variant_harmonisation.chromosome_label_policy
    return resolve_genomic_analysis_scope(
        genome_build=module.genome_build,
        genomes=configuration.resources.genomes,
        mhc=module.mhc,
        chromosomes=module.chromosomes,
        normalize_chromosome=lambda value: _normalise_chromosome(
            str(value), policy,
        ),
        analysis_name="GCTA fastBAT/mBAT-combo",
        error_type=GctaGeneError,
    )


def _allele_key(first: str, second: str) -> tuple[str, str]:
    alleles = tuple(sorted((str(first).strip().upper(), str(second).strip().upper())))
    if not all(alleles) or alleles[0] == alleles[1]:
        raise GctaGeneError("Variant matching requires two nonempty, distinct alleles.")
    return alleles


def _validate_reference(
    prefix: Path,
    module,
    identifiers: set[str],
    alleles: dict[str, tuple[str, str]],
    *,
    analysis_scope: dict,
    retain_variant_ids: bool,
    direct_input_mode: bool,
) -> tuple[dict, set[str], tuple[str, ...]]:
    validate_plink_files(
        prefix, module.reference.required_extensions, error_type=GctaGeneError,
    )
    bim_extensions = [
        suffix for suffix in module.reference.required_extensions
        if suffix.lower() == ".bim"
    ]
    if len(bim_extensions) != 1:
        raise GctaGeneError(
            "reference.required_extensions must contain exactly one .bim suffix."
        )
    bim_path = Path(str(prefix) + bim_extensions[0])
    column_count = len(module.reference.bim_columns)
    column_index = {
        role: index for index, role in enumerate(module.reference.bim_columns)
    }
    reference_ids = set()
    reference_chromosomes = set()
    compatible_alleles = 0
    incompatible_alleles = 0
    exact_matches: set[str] = set()
    analysis_exact_matches: set[str] = set()
    excluded_reference_ids: list[str] = []
    excluded_reference_chromosome_counts = {
        chromosome: 0
        for chromosome in analysis_scope["exclude_chromosomes"]
    }
    excluded_reference_mhc_variants = 0
    excluded_input_chromosome_variants = 0
    excluded_input_mhc_variants = 0
    policy = module.variant_harmonisation
    delimiter = re.compile(module.reference.table_delimiter_pattern)
    with bim_path.open("r", encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            fields = delimiter.split(raw.strip())
            if len(fields) != column_count:
                raise GctaGeneError(
                    "PLINK BIM line %d has %d fields; expected %d."
                    % (number, len(fields), column_count)
                )
            chromosome = fields[column_index["chromosome"]]
            variant = fields[column_index["variant_id"]]
            if not variant or variant in reference_ids:
                raise GctaGeneError(
                    "PLINK BIM contains an empty or duplicate variant ID at line %d."
                    % number
                )
            reference_ids.add(variant)
            normalized_chromosome = _normalise_chromosome(
                chromosome, policy.chromosome_label_policy,
            )
            reference_chromosomes.add(normalized_chromosome)
            try:
                position = int(fields[column_index["position"]])
            except ValueError as exc:
                raise GctaGeneError(
                    "PLINK BIM has an invalid position at line %d." % number
                ) from exc
            if position < 1:
                raise GctaGeneError(
                    "PLINK BIM has a non-positive position at line %d." % number
                )
            exclusion_reason = None
            if normalized_chromosome in excluded_reference_chromosome_counts:
                exclusion_reason = "chromosome"
                excluded_reference_chromosome_counts[normalized_chromosome] += 1
            elif analysis_scope["exclude_mhc_snps"]:
                mhc = analysis_scope["mhc_region"]
                if mhc is None:
                    raise GctaGeneError("Resolved GCTA MHC region is missing")
                if (
                    normalized_chromosome == mhc["chromosome"]
                    and mhc["start"] <= position <= mhc["end"]
                ):
                    exclusion_reason = "mhc"
                    excluded_reference_mhc_variants += 1
            if exclusion_reason is not None:
                excluded_reference_ids.append(variant)
            reference_alleles = _allele_key(
                fields[column_index["allele1"]],
                fields[column_index["allele2"]],
            )
            if variant in identifiers:
                exact_matches.add(variant)
                if exclusion_reason is None:
                    analysis_exact_matches.add(variant)
                elif exclusion_reason == "chromosome":
                    excluded_input_chromosome_variants += 1
                else:
                    excluded_input_mhc_variants += 1
                summary = set(alleles[variant])
                reference = set(reference_alleles)
                complemented = {
                    allele.translate(_DNA_COMPLEMENT) for allele in summary
                }
                complement_allowed = (
                    policy.allow_strand_complement
                    and all(len(allele) == 1 and allele in "ACGT" for allele in summary)
                )
                if summary == reference or (
                    complement_allowed and complemented == reference
                ):
                    compatible_alleles += 1
                else:
                    incompatible_alleles += 1
    overlap = len(exact_matches)
    if overlap == 0:
        if direct_input_mode:
            raise GctaGeneError(
                "No direct GCTA summary-statistic SNP IDs occur exactly in "
                "PLINK BIM column 2. Direct mode does not rewrite identifiers. "
                "Provide --gcta-input-file containing BIM-compatible IDs."
            )
        raise GctaGeneError(
            "No GCTA input variant identifiers occur in the PLINK BIM reference."
        )
    overlap_fraction = overlap / len(identifiers)
    if (
        not direct_input_mode
        and overlap_fraction < policy.minimum_overlap_fraction
    ):
        raise GctaGeneError(
            "Only %d/%d unique GCTA input variants (%.2f%%) match exact PLINK "
            "BIM IDs; the configured minimum is %.2f%%. Review the selected "
            "LD reference and formatter identifier configuration."
            % (
                overlap,
                len(identifiers),
                overlap_fraction * 100,
                policy.minimum_overlap_fraction * 100,
            )
        )
    if incompatible_alleles:
        raise GctaGeneError(
            "%d overlapping GCTA .ma variants have allele pairs that cannot be "
            "matched to the PLINK BIM reference." % incompatible_alleles
        )
    if not analysis_exact_matches:
        raise GctaGeneError(
            "The configured chromosome and MHC policies exclude every "
            "GWAS/BIM-shared variant; no GCTA analysis input remains."
        )
    metrics = {
        "reference_variants": len(reference_ids),
        "overlapping_variants": overlap,
        "overlap_fraction": overlap_fraction,
        "input_variants_absent_from_reference": len(identifiers) - overlap,
        "reference_variants_absent_from_input": len(reference_ids) - overlap,
        "exact_id_matches": overlap,
        "compatible_allele_pairs": compatible_alleles,
        "incompatible_allele_pairs": incompatible_alleles,
        "analyzable_variants": len(analysis_exact_matches),
        "excluded_reference_variants": len(excluded_reference_ids),
        "excluded_reference_chromosome_variants": sum(
            excluded_reference_chromosome_counts.values()
        ),
        "excluded_reference_chromosome_counts": (
            excluded_reference_chromosome_counts
        ),
        "excluded_reference_mhc_variants": excluded_reference_mhc_variants,
        "excluded_input_variants": (
            excluded_input_chromosome_variants + excluded_input_mhc_variants
        ),
        "excluded_input_chromosome_variants": (
            excluded_input_chromosome_variants
        ),
        "excluded_input_mhc_variants": excluded_input_mhc_variants,
        "analysis_scope": analysis_scope,
        "reference_chromosomes": sorted(reference_chromosomes),
        "variant_id_match_policy": (
            "report_exact_bim_overlap_without_rewriting"
            if direct_input_mode
            else "validate_pipeline_formatter_bim_compatibility_without_rewriting"
        ),
    }
    return (
        metrics,
        analysis_exact_matches if retain_variant_ids else set(),
        tuple(excluded_reference_ids),
    )


def _reference_bim_path(prefix: Path, module) -> Path:
    extensions = [
        suffix for suffix in module.reference.required_extensions
        if suffix.lower() == ".bim"
    ]
    if len(extensions) != 1:
        raise GctaGeneError(
            "reference.required_extensions must contain exactly one .bim suffix."
        )
    return Path(str(prefix) + extensions[0])


def _write_deterministic_scope_file(
    path: Path,
    text: str,
    *,
    overwrite: bool,
) -> Path:
    """Atomically publish one deterministic, PostGWAS-owned scope input."""
    if path.exists():
        try:
            if path.read_text(encoding="utf-8") == text:
                return path
        except (OSError, UnicodeError) as exc:
            raise GctaGeneError("Cannot validate analysis-scope file %s: %s" % (path, exc)) from exc
        if not overwrite:
            raise GctaGeneError(
                "Analysis-scope file exists but does not match the current "
                "inputs and policies: %s. Use --overwrite after review." % path
            )
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name = None
    try:
        with tempfile.NamedTemporaryFile(
            "w", encoding="utf-8", dir=path.parent, delete=False,
        ) as handle:
            temporary_name = handle.name
            handle.write(text)
        Path(temporary_name).replace(path)
        temporary_name = None
    finally:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)
    return path


def _prepare_variant_exclusion_list(
    output: Path,
    dataset_id: str,
    module,
    excluded_reference_ids: tuple[str, ...],
    *,
    overwrite: bool,
) -> Path | None:
    """Write exact BIM IDs for GCTA's documented ``--exclude`` option."""
    if not excluded_reference_ids:
        return None
    destination = configured_output_path(
        output,
        module.output_layout.variant_exclusion_list,
        error_type=GctaGeneError,
        dataset_id=dataset_id,
        method=module.method,
    )
    return _write_deterministic_scope_file(
        destination,
        "".join("%s\n" % variant for variant in excluded_reference_ids),
        overwrite=overwrite,
    )


def _prepare_scoped_gene_list(
    source: Path,
    output: Path,
    dataset_id: str,
    module,
    analysis_scope: dict,
    *,
    overwrite: bool,
) -> tuple[Path, Path | None, dict, set[str]]:
    """Remove excluded chromosome/MHC-overlapping tested gene intervals."""
    column_index = {
        role: index for index, role in enumerate(module.gene_annotation.columns)
    }
    window_bp = module.gene_window_kb * 1000
    retained_lines: list[str] = []
    excluded_lines: list[str] = []
    excluded_gene_ids: set[str] = set()
    excluded_chromosome_genes = 0
    excluded_mhc_genes = 0
    retained_gene_chromosomes: set[str] = set()
    mhc = analysis_scope["mhc_region"]
    with source.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip():
                continue
            fields = raw.split()
            chromosome = _normalise_chromosome(
                fields[column_index["chromosome"]],
                (
                    module.set_annotation.conversion.chromosome_label_policy
                    if module.method == "fastbat_set"
                    else module.variant_harmonisation.chromosome_label_policy
                ),
            )
            start = int(fields[column_index["start"]])
            end = int(fields[column_index["end"]])
            gene = fields[column_index["gene"]]
            reason = None
            if chromosome in analysis_scope["exclude_chromosomes"]:
                reason = "chromosome"
                excluded_chromosome_genes += 1
            elif analysis_scope["exclude_mhc_genes"]:
                if mhc is None:
                    raise GctaGeneError("Resolved GCTA MHC region is missing")
                tested_start = max(1, start - window_bp)
                tested_end = end + window_bp
                if (
                    chromosome == mhc["chromosome"]
                    and tested_start <= mhc["end"]
                    and tested_end >= mhc["start"]
                ):
                    reason = "mhc"
                    excluded_mhc_genes += 1
            line = raw.rstrip("\r\n")
            if reason is None:
                retained_lines.append(line)
                retained_gene_chromosomes.add(chromosome)
            else:
                excluded_lines.append(line)
                excluded_gene_ids.add(gene)
    metrics = {
        "input_genes": len(retained_lines) + len(excluded_lines),
        "retained_genes": len(retained_lines),
        "excluded_genes": len(excluded_lines),
        "excluded_chromosome_genes": excluded_chromosome_genes,
        "excluded_mhc_genes": excluded_mhc_genes,
        "retained_gene_chromosomes": sorted(retained_gene_chromosomes),
        "mhc_gene_interval_definition": "gene_coordinates_plus_configured_window",
    }
    if not excluded_lines:
        return source, None, metrics, excluded_gene_ids
    if not retained_lines:
        raise GctaGeneError(
            "The configured chromosome and MHC policies exclude every gene "
            "in the GCTA coordinate file."
        )
    scoped = configured_output_path(
        output,
        module.output_layout.scoped_gene_list,
        error_type=GctaGeneError,
        dataset_id=dataset_id,
        method=module.method,
    )
    excluded = configured_output_path(
        output,
        module.output_layout.excluded_gene_list,
        error_type=GctaGeneError,
        dataset_id=dataset_id,
        method=module.method,
    )
    _write_deterministic_scope_file(
        scoped,
        "\n".join(retained_lines) + "\n",
        overwrite=overwrite,
    )
    _write_deterministic_scope_file(
        excluded,
        "\n".join(excluded_lines) + "\n",
        overwrite=overwrite,
    )
    return scoped, excluded, metrics, excluded_gene_ids


def _prepared_resource_is_current(
    directory: Path,
    module,
    gmt: Path,
    gene_list: Path,
    bim: Path,
    analysis_source: Path,
    analyzable_variant_count: int,
    analysis_scope: dict,
) -> Path:
    conversion = module.set_annotation.conversion
    names = conversion.output_names
    manifest_path = directory / names.manifest
    try:
        manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise GctaGeneError(
            "Cannot validate the existing GMT conversion resource %s: %s. "
            "Use --overwrite after review." % (directory, exc)
        ) from exc
    expected_sources = {
        "gmt": gmt,
        "gene_coordinates": gene_list,
        "plink_bim": bim,
    }
    expected_policies = {
        "allowed_chromosomes": sorted(conversion.allowed_chromosomes),
        "chromosome_label_policy": conversion.chromosome_label_policy,
        "duplicate_gene_policy": conversion.duplicate_gene_policy,
        "unmapped_gene_policy": conversion.unmapped_gene_policy,
        "minimum_gene_id_overlap_fraction": (
            conversion.minimum_gene_id_overlap_fraction
        ),
        "empty_pathway_policy": conversion.empty_pathway_policy,
        "oversized_set_policy": module.set_annotation.oversized_set_policy,
        "maximum_set_variants": module.set_annotation.maximum_set_variants,
        "maximum_set_variants_applied": True,
        "variant_universe": "gwas_bim_intersection",
        "analysis_scope": analysis_scope,
        "audit": conversion.audit.model_dump(mode="json"),
        "disk": conversion.disk.model_dump(mode="json"),
    }
    expected_outputs = {
        "pathway_mapping": names.pathway_mapping,
        "pathway_gene_mapping": (
            names.pathway_gene_mapping
            if conversion.audit.level in {"normalized", "expanded"} else None
        ),
        "gene_variant_mapping": (
            names.gene_variant_mapping
            if conversion.audit.level in {"normalized", "expanded"} else None
        ),
        "expanded_mapping": (
            names.expanded_mapping
            if conversion.audit.level == "expanded" else None
        ),
        "unmapped_genes": names.unmapped_genes,
        "manifest": names.manifest,
        "readme": names.readme,
        "checksums": names.checksums,
    }
    try:
        current = (
            manifest["schema_version"] == "gcta_fastbat_pathway_resource.v6"
            and manifest["resource"]["genome_build"] == module.genome_build
            and manifest["resource"]["gene_window_kb"] == module.gene_window_kb
            and all(
                manifest["sources"][label]["path"] == str(path)
                and manifest["sources"][label]["sha256"] == sha256(path)
                for label, path in expected_sources.items()
            )
            and all(
                manifest["policies"][key] == value
                for key, value in expected_policies.items()
            )
            and manifest["sources"]["analysis_variants"] == {
                "path": str(analysis_source),
                "sha256": sha256(analysis_source),
                "analyzable_unique_variants": analyzable_variant_count,
            }
            and manifest["outputs"] == expected_outputs
        )
        set_path = directory / manifest["resource"]["file"]
        checksum_lines = (
            directory / names.checksums
        ).read_text(encoding="utf-8").splitlines()
        recorded_checksums = {
            filename: digest
            for line in checksum_lines
            for digest, filename in [line.split(maxsplit=1)]
        }
        tracked_outputs = [
            value for key, value in expected_outputs.items()
            if value is not None and key != "checksums"
        ]
        current = (
            current
            and set_path.name == names.set_list
            and set_path.is_file()
            and manifest["resource"]["sha256"] == sha256(set_path)
            and all(
                (directory / name).is_file()
                and recorded_checksums.get(name) == sha256(directory / name)
                for name in tracked_outputs
            )
        )
    except (KeyError, TypeError, ValueError, OSError):
        current = False
        set_path = directory / names.set_list
    if not current:
        raise GctaGeneError(
            "Existing GMT conversion resource does not match the current inputs, "
            "genome build, window, policies, or checksums: %s. Use --overwrite "
            "after review." % directory
        )
    return set_path


def _prepare_gmt_set(
    output: Path,
    dataset_id: str,
    reference_prefix: Path,
    module,
    configuration,
    logger: PipelineLogger,
    bim_variant_total: int,
    analyzable_variant_ids: set[str],
    analysis_variant_source: Path,
    *,
    pathway_gene_preflight=None,
    excluded_gene_ids: set[str] | None = None,
    analysis_scope: dict | None = None,
    pipeline_args=None,
) -> tuple[Path, dict, Path]:
    gmt = _required_path(module.set_annotation.gmt_file, "GMT pathway file")
    gene_list = _required_path(
        module.gene_annotation.file,
        "GCTA gene list required for GMT gene-to-variant conversion",
    )
    bim = _reference_bim_path(reference_prefix, module)
    conversion = module.set_annotation.conversion
    names = conversion.output_names
    prepared_directory = configured_output_path(
        output,
        module.output_layout.prepared_set_directory,
        error_type=GctaGeneError,
        dataset_id=dataset_id,
        method=module.method,
    )
    if prepared_directory.exists():
        if configuration.run.overwrite:
            shutil.rmtree(prepared_directory)
        elif configuration.run.resume:
            set_path = _prepared_resource_is_current(
                prepared_directory,
                module,
                gmt,
                gene_list,
                bim,
                analysis_variant_source,
                len(analyzable_variant_ids),
                analysis_scope or {},
            )
            logger.record(
                "SKIP", "gmt_to_fastbat_set", reason="validated_resume",
                output=str(set_path),
            )
            manifest = yaml.safe_load(
                (prepared_directory / names.manifest).read_text(encoding="utf-8")
            )
            if pipeline_args is not None:
                from postgwas.modules.gcta_gene.stages import (
                    complete_pipeline_stage,
                    start_pipeline_stage,
                )

                validation = manifest["validation"]
                resumed_stages = (
                    (
                        "map_variants",
                        (
                            (
                                "count", "PLINK BIM variants scanned",
                                validation["bim_variants"],
                            ),
                            (
                                "count", "Analyzable BIM variants cached",
                                validation["analyzable_bim_variants_cached"],
                            ),
                            (
                                "success", "Variants mapped to matched genes",
                                validation[
                                    "analyzable_bim_variants_mapped_to_requested_genes"
                                ],
                            ),
                            (
                                "genetic", "Shared chromosomes",
                                ", ".join(validation["shared_chromosomes"]),
                            ),
                        ),
                    ),
                    (
                        "candidate_memberships",
                        (
                            (
                                "count", "Input pathways",
                                validation["input_pathways"],
                            ),
                            (
                                "success", "Candidate non-empty pathways",
                                validation["input_pathways"]
                                - validation["pathways_omitted_empty"],
                            ),
                            (
                                "warning"
                                if validation["pathways_omitted_empty"]
                                else "success",
                                "Pathways without analyzable variants",
                                validation["pathways_omitted_empty"],
                            ),
                            (
                                "success",
                                "Pathway genes represented in coordinate reference",
                                validation["genes_with_coordinates"],
                            ),
                            (
                                "warning"
                                if validation["genes_excluded_by_analysis_scope"]
                                else "success",
                                "Matched genes excluded by analysis scope",
                                validation["genes_excluded_by_analysis_scope"],
                            ),
                            (
                                "success",
                                "Matched genes eligible for analysis",
                                validation["genes_eligible_for_analysis"],
                            ),
                            (
                                "warning"
                                if validation["unmapped_unique_genes"]
                                else "success",
                                "Unmatched pathway genes",
                                validation["unmapped_unique_genes"],
                            ),
                            (
                                "success",
                                "Pathways with complete gene-ID mapping",
                                validation[
                                    "pathways_with_complete_gene_id_mapping"
                                ],
                            ),
                            (
                                "warning"
                                if validation[
                                    "pathways_with_partial_gene_id_mapping"
                                ] else "success",
                                "Pathways with partial gene-ID mapping",
                                validation[
                                    "pathways_with_partial_gene_id_mapping"
                                ],
                            ),
                            (
                                "warning"
                                if validation["pathways_without_gene_id_mapping"]
                                else "success",
                                "Pathways without gene-ID mapping",
                                validation["pathways_without_gene_id_mapping"],
                            ),
                        ),
                    ),
                    (
                        "set_policies",
                        (
                            (
                                "count", "Input pathways",
                                validation["input_pathways"],
                            ),
                            (
                                "success", "Final fastBAT sets written",
                                validation["pathways_written"],
                            ),
                            (
                                "warning"
                                if validation["pathways_omitted_empty"]
                                else "success",
                                "Empty pathways omitted",
                                validation["pathways_omitted_empty"],
                            ),
                            (
                                "warning"
                                if validation["pathways_omitted_oversized"]
                                else "success",
                                "Oversized pathways omitted",
                                validation["pathways_omitted_oversized"],
                            ),
                            (
                                "count", "Unique variants written",
                                validation["unique_variants_written"],
                            ),
                        ),
                    ),
                )
                for key, outcome_fields in resumed_stages:
                    start_pipeline_stage(pipeline_args, key, logger)
                    complete_pipeline_stage(
                        pipeline_args,
                        key,
                        outcome="Reused checksum-validated pathway resource",
                        outcome_fields=outcome_fields,
                        logger=logger,
                    )
            return set_path, manifest, prepared_directory
        else:
            raise GctaGeneError(
                "Prepared GMT set resource already exists: %s. Use --resume to "
                "reuse the validated resource or --overwrite to rebuild it."
                % prepared_directory
            )
    arguments = SimpleNamespace(
        gmt=gmt,
        gene_list=gene_list,
        bim=bim,
        output_directory=prepared_directory,
        output_name=names.set_list,
        pathway_mapping_name=names.pathway_mapping,
        pathway_gene_mapping_name=names.pathway_gene_mapping,
        gene_variant_mapping_name=names.gene_variant_mapping,
        expanded_mapping_name=names.expanded_mapping,
        unmapped_genes_name=names.unmapped_genes,
        manifest_name=names.manifest,
        readme_name=names.readme,
        checksums_name=names.checksums,
        resource_name="%s GCTA fastBAT pathway sets" % dataset_id,
        genome_build=module.genome_build,
        gene_window_kb=module.gene_window_kb,
        allowed_chromosomes=conversion.allowed_chromosomes,
        chromosome_label_policy=conversion.chromosome_label_policy,
        duplicate_gene_policy=conversion.duplicate_gene_policy,
        unmapped_gene_policy=conversion.unmapped_gene_policy,
        minimum_gene_id_overlap_fraction=(
            conversion.minimum_gene_id_overlap_fraction
        ),
        empty_pathway_policy=conversion.empty_pathway_policy,
        gene_columns=module.gene_annotation.columns,
        generation_command=shlex.join(sys.argv),
    )
    try:
        conversion_progress = StageProgress(
            "GMT-to-fastBAT preparation",
            enabled=(
                configuration.logging.show_progress
                and pipeline_args is None
            ),
        )
        manifest = prepare_resource(
            arguments,
            pathway_gene_preflight=pathway_gene_preflight,
            stage_progress=conversion_progress,
            pipeline_args=pipeline_args,
            pipeline_logger=logger,
            measured_progress_enabled=configuration.logging.show_progress,
            bim_variant_total=bim_variant_total,
            progress_refresh_seconds=(
                configuration.logging.progress_refresh_seconds
            ),
            mapping_workers=configuration.execution.threads,
            mapping_memory_gb=configuration.execution.memory_gb,
            worker_memory_multiplier=(
                conversion.parallelism.worker_memory_multiplier
            ),
            minimum_gene_id_overlap_fraction=(
                conversion.minimum_gene_id_overlap_fraction
            ),
            bim_ids_prevalidated=True,
            analyzable_variant_ids=analyzable_variant_ids,
            analysis_variant_source=analysis_variant_source,
            maximum_set_variants=module.set_annotation.maximum_set_variants,
            oversized_set_policy=module.set_annotation.oversized_set_policy,
            audit_level=conversion.audit.level,
            audit_format=conversion.audit.format,
            audit_compression=conversion.audit.compression,
            audit_batch_rows=conversion.audit.batch_rows,
            minimum_free_disk_gb=conversion.disk.minimum_free_gb,
            disk_estimation_safety_factor=(
                conversion.disk.estimation_safety_factor
            ),
            excluded_gene_ids=excluded_gene_ids,
            analysis_scope=analysis_scope,
        )
    except (OSError, ResourcePreparationError, yaml.YAMLError) as exc:
        raise GctaGeneError("GMT-to-fastBAT conversion failed: %s" % exc) from exc
    set_path = prepared_directory / names.set_list
    logger.record(
        "COMPLETED",
        "gmt_to_fastbat_set",
        input_pathways=manifest["validation"]["input_pathways"],
        pathways_written=manifest["validation"]["pathways_written"],
        unique_variants_written=manifest["validation"]["unique_variants_written"],
        pathways_omitted_empty=manifest["validation"]["pathways_omitted_empty"],
        pathways_omitted_oversized=(
            manifest["validation"]["pathways_omitted_oversized"]
        ),
        gmt_gene_mappability_fraction=manifest["validation"][
            "gmt_gene_mappability_fraction"
        ],
        coordinate_reference_coverage_fraction=manifest["validation"][
            "coordinate_reference_coverage_fraction"
        ],
        coordinate_reference_genes_absent_from_gmt=manifest["validation"][
            "coordinate_reference_genes_absent_from_gmt"
        ],
        pathways_with_complete_gene_id_mapping=manifest["validation"][
            "pathways_with_complete_gene_id_mapping"
        ],
        pathways_with_partial_gene_id_mapping=manifest["validation"][
            "pathways_with_partial_gene_id_mapping"
        ],
        pathways_without_gene_id_mapping=manifest["validation"][
            "pathways_without_gene_id_mapping"
        ],
        audit_level=manifest["policies"]["audit"]["level"],
        estimated_output_bytes=manifest["disk_preflight"][
            "estimated_output_bytes"
        ],
        mapping_workers=manifest["generation"]["mapping_workers"],
        mapping_algorithm=manifest["generation"]["mapping_algorithm"],
        output=str(set_path),
    )
    return set_path, manifest, prepared_directory


def _validate_gene_list(
    path: Path,
    module,
    reference_chromosomes: set[str] | None = None,
) -> dict:
    rows = read_gene_coordinates(
        path, column_roles=module.gene_annotation.columns,
        error_type=GctaGeneError,
    )
    chromosomes = {
        _normalise_chromosome(row.chromosome, module.variant_harmonisation.chromosome_label_policy)
        for row in rows
    }
    overlap = (
        chromosomes & reference_chromosomes
        if reference_chromosomes is not None else set()
    )
    if reference_chromosomes is not None and not overlap:
        raise GctaGeneError(
            "Gene list and LD reference do not share any chromosome labels."
        )
    return {
        "genes": len(rows),
        "gene_chromosomes": sorted(chromosomes),
        "shared_chromosomes": sorted(overlap),
    }


def preflight_gcta_gene_pipeline(
    args,
    *,
    preflight_evidence=None,
) -> PipelinePreflightEvidence:
    """Validate every external GCTA gene resource before VCF formatting."""
    entry_vcf = require_pipeline_input_vcf(preflight_evidence)
    configuration = _resolved_configuration(args)
    module = configuration.modules.gcta_gene
    _require_gcta_gene_arguments(
        configuration,
        None,
        include_generated_input=False,
    )
    observed_build = str(entry_vcf["harmonised"]["genome_build"])
    if str(module.genome_build) != observed_build:
        raise GctaGeneError(
            "The harmonised GWAS-VCF declares genome build %s, but GCTA gene "
            "analysis resolves to %s. The GWAS, PLINK LD reference, and gene "
            "coordinates must use the same build."
            % (observed_build, module.genome_build)
        )
    _validate_build_and_population(module, configuration)
    analysis_scope = _analysis_scope(configuration)
    reference_prefix, reference_paths = validate_gcta_reference_files(module)
    executable = resolve_executable(
        configuration.resources.executables.gcta,
        "GCTA executable",
        error_type=GctaGeneError,
    )
    version = require_supported_gcta(
        executable,
        module,
        None,
        configuration.execution.timeout_seconds,
    )
    bim_path = _reference_bim_path(reference_prefix, module)
    configure_reference_variant_identifiers(
        args,
        configuration.modules.formatting,
        [BimIdentifierRequirement(
            consumer="GCTA %s analysis" % module.method,
            formatter_target="gcta_gene",
            bim_file=bim_path,
            column_roles=module.reference.bim_columns,
            delimiter_pattern=module.reference.table_delimiter_pattern,
        )],
    )
    identifier_observation = dict(args.variant_id_observations["gcta_gene"])
    validate_plink_bundle_dimensions(
        {path.suffix: path for path in reference_paths},
        variants=identifier_observation["variants"], error_type=GctaGeneError,
    )

    gene_list = None
    gene_annotation_metrics = None
    gmt = None
    pathway_preflight = None
    set_list = None
    set_source_metrics = None
    if module.method in {"fastbat_gene", "mbat_combo"}:
        if module.set_annotation.gmt_file is not None:
            raise GctaGeneError(
                "--gmt is valid only with --method fastbat_set; the analysis "
                "method is never inferred from an input filename."
            )
        gene_list = _required_path(
            module.gene_annotation.file,
            "GCTA gene-coordinate file",
        )
        gene_annotation_metrics = _validate_gene_list(gene_list, module)
    elif module.method == "fastbat_set":
        if module.set_annotation.gmt_file is not None:
            conversion = module.set_annotation.conversion
            gene_list = _required_path(
                module.gene_annotation.file,
                "GCTA gene-coordinate file",
            )
            gmt = _required_path(
                module.set_annotation.gmt_file,
                "GMT pathway file",
            )
            try:
                all_gene_coordinates, gene_chromosomes = (
                    validate_gene_coordinate_reference(
                        gene_list,
                        set(conversion.allowed_chromosomes),
                        conversion.chromosome_label_policy,
                        module.gene_window_kb * 1000,
                        module.gene_annotation.columns,
                    )
                )
                pathway_source = scan_gmt_pathways(
                    gmt,
                    conversion.duplicate_gene_policy,
                )
                pathway_preflight = validate_pathway_gene_compatibility(
                    gmt=gmt,
                    gene_list=gene_list,
                    allowed_chromosomes=set(conversion.allowed_chromosomes),
                    chromosome_policy=conversion.chromosome_label_policy,
                    window_bp=module.gene_window_kb * 1000,
                    column_roles=module.gene_annotation.columns,
                    duplicate_gene_policy=conversion.duplicate_gene_policy,
                    unmapped_gene_policy=conversion.unmapped_gene_policy,
                    minimum_gene_id_overlap_fraction=(
                        conversion.minimum_gene_id_overlap_fraction
                    ),
                    all_gene_coordinates=all_gene_coordinates,
                    pathway_scan=pathway_source,
                )
            except ResourcePreparationError as exc:
                raise GctaGeneError(str(exc)) from exc
            gene_annotation_metrics = {
                "genes": len(all_gene_coordinates),
                "gene_chromosomes": sorted(gene_chromosomes),
            }
        else:
            set_list = _required_path(
                module.set_annotation.file,
                "GCTA fastBAT set list",
            )
            set_source_metrics = validate_fastbat_set_list(set_list, error_type=GctaGeneError)
    elif module.set_annotation.gmt_file is not None:
        raise GctaGeneError("--gmt is valid only with --method fastbat_set.")

    if (
        module.mhc.policy == "exclude_genes"
        and (
            module.method == "fastbat_segment"
            or (
                module.method == "fastbat_set"
                and module.set_annotation.gmt_file is None
            )
        )
    ):
        raise GctaGeneError(
            "--mhc-policy exclude_genes requires a coordinate-defined gene "
            "input. Fixed segments and prepared SNP-set files contain no gene "
            "identities; use include, exclude_snps, or exclude_both."
        )

    resource_paths = [*reference_paths, executable]
    resource_paths.extend(
        path for path in (gene_list, gmt, set_list) if path is not None
    )
    resources = GctaGenePipelineResources(
        configuration=configuration,
        reference_prefix=reference_prefix,
        reference_paths=reference_paths,
        executable=executable,
        version=version,
        identifier_observation=identifier_observation,
        analysis_scope=analysis_scope,
        gene_list=gene_list,
        gene_annotation_metrics=gene_annotation_metrics,
        gmt=gmt,
        pathway_gene_preflight=pathway_preflight,
        set_list=set_list,
        set_source_metrics=set_source_metrics,
        file_identities=capture_preflight_file_identities(
            resource_paths,
            error_type=GctaGeneError,
            label="GCTA gene resource",
        ),
    )
    return pipeline_preflight_evidence(
        "gcta_gene",
        preflight_evidence,
        resources=resources,
        deferred_checks=(
            "Validate the formatter-created GCTA .ma table.",
            "Validate exact SNP-ID and allele-pair compatibility with the BIM.",
            "Apply analysis-scope and set-membership policies to analyzable variants.",
        ),
    )


def _gcta_gene_pipeline_execution_configuration(args, resources):
    """Retain validated settings while applying the orchestrator stage path."""
    configuration = resources.configuration
    run = configuration.run.model_copy(update={
        "output_directory": Path(args.output_directory).expanduser().resolve(),
    })
    return configuration.model_copy(update={"run": run}, deep=True)


def _prepare_analysis_set_list(
    path: Path,
    destination: Path,
    identifiers: set[str],
    reference_ids: set[str],
    empty_set_policy: str,
    oversized_set_policy: str,
    maximum_set_variants: int,
) -> tuple[Path, dict]:
    """Validate and atomically write sets restricted to analyzable variants.

    GCTA first restricts custom-set IDs to the LD reference and then to the
    summary-statistics input.  Removing unavailable IDs here is scientifically
    equivalent for retained sets and lets the configured empty-set policy
    handle sets that would otherwise reach GCTA with zero usable variants.
    """
    sets = set()
    omitted_sets: list[str] = []
    oversized_sets: list[tuple[str, int]] = []
    matched_variants = 0
    analysis_variants = 0
    current_matched: list[str] = []
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_name = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=destination.parent,
            delete=False,
        ) as output:
            temporary_name = output.name
            with open_fastbat_set_memberships(path, error_type=GctaGeneError) as (source_metrics, memberships):
                for current_set, value in memberships:
                    if value is None:
                        if len(current_matched) > maximum_set_variants:
                            oversized_sets.append(
                                (current_set, len(current_matched))
                            )
                        elif current_matched:
                            output.write(current_set + "\n")
                            output.writelines(
                                variant + "\n" for variant in current_matched
                            )
                            output.write("END\n\n")
                            sets.add(current_set)
                            analysis_variants += len(current_matched)
                        else:
                            omitted_sets.append(current_set)
                        current_matched = []
                        continue
                    if value in identifiers and value in reference_ids:
                        current_matched.append(value)
                        matched_variants += 1
        if omitted_sets and empty_set_policy == "error":
            raise GctaGeneError(
                "%d fastBAT sets have no variants shared by the GWAS input and "
                "PLINK LD reference; examples: %s"
                % (len(omitted_sets), ", ".join(omitted_sets[:10]))
            )
        if oversized_sets and oversized_set_policy == "error":
            raise GctaGeneError(
                "%d fastBAT sets exceed GCTA's %s-variant hard limit; "
                "examples: %s"
                % (
                    len(oversized_sets),
                    format(maximum_set_variants, ","),
                    ", ".join(
                        "%s (%d)" % item for item in oversized_sets[:10]
                    ),
                )
            )
        if not sets:
            raise GctaGeneError(
                "No fastBAT sets remain after GWAS/reference intersection and "
                "GCTA set-size validation."
            )
        Path(temporary_name).replace(destination)
        temporary_name = None
    finally:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)
    return destination, {
        "sets": len(sets),
        "input_sets": source_metrics["input_sets"],
        "omitted_empty_sets": len(omitted_sets),
        "omitted_empty_set_examples": omitted_sets[:10],
        "omitted_oversized_sets": len(oversized_sets),
        "omitted_oversized_set_examples": [
            "%s (%d)" % item for item in oversized_sets[:10]
        ],
        "maximum_set_variants": maximum_set_variants,
        "requested_set_variants": source_metrics["requested_set_variants"],
        "unique_requested_set_variants": source_metrics["unique_requested_set_variants"],
        "matched_set_variants": matched_variants,
        "analysis_set_variants": analysis_variants,
        "unmatched_set_variants": source_metrics["requested_set_variants"] - matched_variants,
        "analysis_set_list": str(destination),
    }


def _method_result_path(output: Path, dataset_id: str, module) -> Path:
    pattern = module.output_layout.primary_results[module.method]
    return configured_output_path(output, pattern, error_type=GctaGeneError, dataset_id=dataset_id)


def _completion_manifest_path(
    output: Path, dataset_id: str, module,
) -> Path:
    return configured_output_path(
        output,
        module.output_layout.completion_manifest,
        error_type=GctaGeneError,
        dataset_id=dataset_id,
        method=module.method,
    )


def _gcta_configuration_digest(configuration) -> str:
    module = configuration.modules.gcta_gene.model_dump(mode="json")
    # Reporting and converter worker scheduling do not determine the GCTA
    # scientific result. Excluding them preserves validated analysis reuse when
    # no GCTA input or command changes.
    module.pop("reporting", None)
    module.pop("html_report", None)
    module["output_layout"].pop("html_report", None)
    conversion = module["set_annotation"]["conversion"]
    conversion.pop("parallelism", None)
    conversion.pop("audit", None)
    conversion.pop("disk", None)
    module["results"].pop("normalized_schema_version", None)
    for schema in module["results"]["schemas"].values():
        schema.pop("reportable_columns", None)
    payload = {
        "module": module,
        "threads": configuration.execution.threads,
        "gcta": str(configuration.resources.executables.gcta),
    }
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _gcta_report_configuration_digest(configuration) -> str:
    module = configuration.modules.gcta_gene
    payload = {
        "reporting": module.reporting.model_dump(mode="json"),
        "html_report": module.html_report.model_dump(mode="json"),
        "output_path": module.output_layout.html_report,
    }
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _completion_fingerprint(path: Path) -> dict:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size <= 0:
        raise GctaGeneError(
            "GCTA completion input or output is missing or empty: %s" % resolved
        )
    return {
        "path": str(resolved),
        "size": resolved.stat().st_size,
        "sha256": sha256(resolved),
    }


def _completion_files(
    source: Path,
    source_annotation_file: Path | None,
    annotation_file: Path | None,
    variant_exclusion_file: Path | None,
    reference_prefix: Path,
    module,
) -> dict[str, Path]:
    files = {"summary_statistics": source}
    if source_annotation_file is not None:
        files["source_annotation"] = source_annotation_file
    if annotation_file is not None:
        files["analysis_annotation"] = annotation_file
    if variant_exclusion_file is not None:
        files["variant_exclusion_list"] = variant_exclusion_file
    for suffix in module.reference.required_extensions:
        files["reference%s" % suffix] = Path(str(reference_prefix) + suffix)
    return files


def _write_gcta_completion_manifest(
    path: Path,
    *,
    dataset_id: str,
    version: str,
    configuration,
    source: Path,
    source_annotation_file: Path | None,
    annotation_file: Path | None,
    variant_exclusion_file: Path | None,
    reference_prefix: Path,
    final_result: Path,
    normalized: Path,
    html_report: Path,
    result_metrics: dict,
) -> None:
    module = configuration.modules.gcta_gene
    inputs = _completion_files(
        source,
        source_annotation_file,
        annotation_file,
        variant_exclusion_file,
        reference_prefix,
        module,
    )
    write_yaml_report(
        {
            "schema_version": _COMPLETION_SCHEMA_VERSION,
            "status": "COMPLETED",
            "dataset_id": dataset_id,
            "method": module.method,
            "gcta_version": version,
            "configuration_sha256": _gcta_configuration_digest(configuration),
            "report_configuration_sha256": (
                _gcta_report_configuration_digest(configuration)
            ),
            "inputs": {
                name: _completion_fingerprint(file_path)
                for name, file_path in inputs.items()
            },
            "outputs": {
                "raw_result": _completion_fingerprint(final_result),
                "normalized_result": _completion_fingerprint(normalized),
                "html_report": _completion_fingerprint(html_report),
            },
            "result_metrics": result_metrics,
        },
        path,
    )


def _validate_gcta_completion_manifest(
    path: Path,
    *,
    dataset_id: str,
    version: str,
    configuration,
    source: Path,
    source_annotation_file: Path | None,
    annotation_file: Path | None,
    variant_exclusion_file: Path | None,
    reference_prefix: Path,
    final_result: Path,
    normalized: Path,
    html_report: Path,
) -> tuple[dict, bool]:
    if not path.is_file():
        raise GctaGeneError(
            "Completed GCTA output has no provenance manifest: %s. Use "
            "--overwrite after reviewing the existing result." % path
        )
    try:
        manifest = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise GctaGeneError(
            "Cannot read GCTA completion manifest %s: %s" % (path, exc)
        ) from exc
    schema_version = manifest.get("schema_version")
    supported_versions = {
        _COMPLETION_SCHEMA_VERSION, *_LEGACY_COMPLETION_SCHEMA_VERSIONS,
    }
    expected_scalars = {
        "status": "COMPLETED",
        "dataset_id": dataset_id,
        "method": configuration.modules.gcta_gene.method,
        "gcta_version": version,
        "configuration_sha256": _gcta_configuration_digest(configuration),
    }
    mismatched = [
        key for key, value in expected_scalars.items()
        if manifest.get(key) != value
    ]
    if schema_version not in supported_versions:
        mismatched.insert(0, "schema_version")
    if mismatched:
        raise GctaGeneError(
            "GCTA completion manifest does not match the current run (%s). Use "
            "--overwrite after review." % ", ".join(mismatched)
        )
    current_inputs = {
        name: _completion_fingerprint(file_path)
        for name, file_path in _completion_files(
            source,
            source_annotation_file,
            annotation_file,
            variant_exclusion_file,
            reference_prefix,
            configuration.modules.gcta_gene,
        ).items()
    }
    current_scientific_outputs = {
        "raw_result": _completion_fingerprint(final_result),
        "normalized_result": _completion_fingerprint(normalized),
    }
    if manifest.get("inputs") != current_inputs:
        raise GctaGeneError(
            "GCTA completion inputs changed since the recorded run. Use "
            "--overwrite after review."
        )
    recorded_outputs = manifest.get("outputs")
    if not isinstance(recorded_outputs, dict) or {
        name: recorded_outputs.get(name)
        for name in current_scientific_outputs
    } != current_scientific_outputs:
        raise GctaGeneError(
            "GCTA completion outputs changed since the recorded run. Use "
            "--overwrite after review."
        )
    report_current = False
    if schema_version == _COMPLETION_SCHEMA_VERSION:
        report_record = recorded_outputs.get("html_report")
        if not isinstance(report_record, dict) or not report_record.get("path"):
            raise GctaGeneError(
                "GCTA completion manifest has no HTML-report fingerprint. Use "
                "--overwrite after review."
            )
        try:
            recorded_report_fingerprint = _completion_fingerprint(
                Path(report_record["path"])
            )
        except GctaGeneError as exc:
            raise GctaGeneError(
                "GCTA HTML report changed or is missing since the recorded run. "
                "Use --overwrite after review."
            ) from exc
        if recorded_report_fingerprint != report_record:
            raise GctaGeneError(
                "GCTA HTML report changed since the recorded run. Use "
                "--overwrite after review."
            )
        recorded_report_path = Path(report_record["path"]).resolve()
        if (
            recorded_report_path != html_report.resolve()
            and (html_report.exists() or html_report.is_symlink())
        ):
            raise GctaGeneError(
                "The newly configured GCTA HTML-report path already exists but "
                "is not tracked by the completion manifest: %s. Move the file "
                "or use --overwrite after review." % html_report
            )
        report_current = (
            manifest.get("report_configuration_sha256")
            == _gcta_report_configuration_digest(configuration)
            and recorded_report_path == html_report.resolve()
        )
    elif html_report.exists() or html_report.is_symlink():
        raise GctaGeneError(
            "The legacy GCTA completion manifest does not track the existing "
            "HTML-report path: %s. Move the file or use --overwrite after "
            "review." % html_report
        )
    metrics = manifest.get("result_metrics")
    if not isinstance(metrics, dict):
        raise GctaGeneError("GCTA completion manifest has no result metrics.")
    return metrics, report_current


def _staged_result_path(final_prefix: Path, final_result: Path, staged_prefix: Path) -> Path:
    if final_result.parent != final_prefix.parent or not final_result.name.startswith(final_prefix.name):
        raise GctaGeneError(
            "Configured primary GCTA result must share the configured output prefix."
        )
    return Path(str(staged_prefix) + final_result.name[len(final_prefix.name):])


def _publish_staged_outputs(staged_prefix: Path, final_prefix: Path, overwrite: bool) -> list[Path]:
    published = []
    for source in sorted(staged_prefix.parent.glob(staged_prefix.name + "*")):
        if not source.is_file():
            continue
        suffix = source.name[len(staged_prefix.name):]
        destination = final_prefix.parent / (final_prefix.name + suffix)
        if destination.exists() and not overwrite:
            raise GctaGeneError("Output already exists: %s" % destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        source.replace(destination)
        published.append(destination)
    return published


def run_gcta_gene(
    input_file: str | Path,
    output_directory: str | Path,
    dataset_id: str,
    configuration,
    logger: PipelineLogger,
    *,
    formatter_variant_id_type: str | None = None,
    input_preflight: GctaInputPreflight | None = None,
    pathway_gene_preflight=None,
    pipeline_resources: GctaGenePipelineResources | None = None,
    pipeline_args=None,
    dry_run: bool = False,
) -> ModuleResult:
    module = configuration.modules.gcta_gene
    pipeline_mode = pipeline_args is not None
    pipeline_stage_mode = bool(
        pipeline_mode
        and getattr(pipeline_args, "_pipeline_progress_plan", {}).get("kind")
        == "gcta_gene"
    )
    if pipeline_resources is not None:
        require_unchanged_preflight_files(
            pipeline_resources.file_identities,
            error_type=GctaGeneError,
            label="GCTA gene resource",
        )
        analysis_scope = dict(pipeline_resources.analysis_scope)
        pathway_gene_preflight = pipeline_resources.pathway_gene_preflight
    else:
        analysis_scope = _analysis_scope(configuration)
    direct_validation_stages = 2
    if pipeline_resources is None and module.method in {"fastbat_gene", "mbat_combo"}:
        direct_validation_stages += 1
    elif pipeline_resources is None and module.method == "fastbat_set":
        direct_validation_stages += 1
    direct_result_stages = 1 if dry_run else 4
    stage_total = direct_validation_stages + direct_result_stages
    stage_number = 0

    @contextmanager
    def stage(title: str, function_name: str, *, pipeline_key: str | None = None):
        nonlocal stage_number
        if pipeline_key is not None and pipeline_stage_mode:
            number = pipeline_stage_number(pipeline_args, pipeline_key)
            if number is None:
                raise GctaGeneError(
                    "The resolved GCTA pipeline plan has no %s stage."
                    % pipeline_key
                )
            active_title = pipeline_args._pipeline_progress_plan["stages"][number - 1]
            active_total = len(pipeline_args._pipeline_progress_plan["stages"])
            start_pipeline_stage(pipeline_args, pipeline_key, logger)
        else:
            stage_number += 1
            number = stage_number
            active_title = title
            active_total = stage_total
        with logger.step(
            number, active_total, active_title, function_name,
        ) as step_context:
            yield step_context
        if pipeline_key is not None and pipeline_stage_mode:
            outcome_values = step_context.extra.get("outcome", {})
            complete_pipeline_stage(
                pipeline_args,
                pipeline_key,
                outcome=outcome_values.get("message"),
                outcome_fields=step_context.outcome_fields,
                logger=logger,
            )

    if pipeline_resources is None:
        _validate_build_and_population(module, configuration)
    output = Path(output_directory).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    layout = module.output_layout
    annotation_file = None
    source_annotation_file = None
    scoped_gene_list = None
    excluded_gene_list = None
    gene_scope_metrics = {}
    excluded_gene_ids: set[str] = set()
    variant_exclusion_file = None
    prepared_manifest = None
    prepared_directory = None
    analysis_set_path = None
    if module.method in {"fastbat_gene", "mbat_combo"}:
        if module.set_annotation.gmt_file is not None:
            raise GctaGeneError(
                "--gmt is valid only with --method fastbat_set; the analysis "
                "method is never inferred from an input filename."
            )
        annotation_file = (
            pipeline_resources.gene_list
            if pipeline_resources is not None
            else _required_path(module.gene_annotation.file, "GCTA gene list")
        )
        if annotation_file is None:
            raise GctaGeneError(
                "GCTA gene-coordinate preflight evidence is incomplete."
            )
        source_annotation_file = annotation_file
    elif module.method == "fastbat_set":
        if (
            module.set_annotation.file is None
            and module.set_annotation.gmt_file is None
        ):
            raise GctaGeneError(
                "fastbat_set requires exactly one of --gmt or --fastbat-set-list."
            )
    elif module.set_annotation.gmt_file is not None:
        raise GctaGeneError("--gmt is valid only with --method fastbat_set.")
    if (
        module.mhc.policy == "exclude_genes"
        and (
            module.method == "fastbat_segment"
            or (
                module.method == "fastbat_set"
                and module.set_annotation.gmt_file is None
            )
        )
    ):
        raise GctaGeneError(
            "--mhc-policy exclude_genes requires a coordinate-defined gene "
            "input. Fixed segments and prepared SNP-set files contain no gene "
            "identities; use include, exclude_snps, or exclude_both."
        )

    if pipeline_mode:
        validated_input = input_preflight or validate_pipeline_gcta_input(
            input_file,
            configuration,
            pipeline_resources=pipeline_resources,
        )
        source = validated_input.source
        input_metrics = dict(validated_input.input_metrics)
        input_metrics["formatter_variant_id_type"] = formatter_variant_id_type
        identifiers = set(validated_input.identifiers)
        reference_metrics = dict(validated_input.reference_metrics)
        reference_ids = set(validated_input.reference_ids)
        excluded_reference_ids = tuple(validated_input.excluded_reference_ids)
        reference_prefix = validated_input.reference_prefix
    else:
        with stage(
            "Validate the GCTA summary-statistics input",
            "validate_gcta_input",
        ) as input_stage:
            source = _required_path(
                str(input_file), "GCTA summary-statistics input",
            )
            reference_prefix, reference_paths = validate_gcta_reference_files(
                module,
            )
            (
                input_metrics,
                input_frame,
                _identifier_column,
                identifiers,
                alleles,
            ) = _validate_formatted_input(source, configuration, module.method)
            input_stage.outcome(
                "Validated the original GCTA summary-statistics input.",
                fields=gcta_ma_input_outcome_fields(
                    configuration,
                    source,
                    input_metrics,
                ),
                input_file=str(source),
                input_format="gcta_ma",
                input_variants=input_metrics["variants"],
                input_columns=input_metrics["columns"],
                required_columns_validated=True,
                scientific_ranges_validated=True,
                unique_variant_identifiers_validated=True,
                allele_pairs_validated=True,
                file_structure_validated=True,
            )
        with stage(
            "Validate the PLINK LD reference and compare it with the GWAS input",
            "validate_gcta_reference",
        ) as reference_stage:
            reference_metrics, reference_ids, excluded_reference_ids = _validate_reference(
                reference_prefix,
                module,
                identifiers,
                alleles,
                analysis_scope=analysis_scope,
                retain_variant_ids=module.method == "fastbat_set",
                direct_input_mode=True,
            )
            absent_from_reference = reference_metrics[
                "input_variants_absent_from_reference"
            ]
            reference_stage.outcome(
                "Validated the PLINK files and compared exact IDs and alleles "
                "without rewriting the direct GCTA input.",
                fields=gcta_ld_reference_outcome_fields(
                    module,
                    reference_prefix,
                    reference_paths,
                    reference_metrics=reference_metrics,
                    variant_id_type=formatter_variant_id_type,
                    summary_variants=len(identifiers),
                    enforce_minimum_overlap=False,
                ),
                reference_prefix=str(reference_prefix),
                reference_files=[str(path) for path in reference_paths],
                declared_genome_build=module.genome_build,
                declared_population=module.reference.population,
                reference_chromosomes=reference_metrics[
                    "reference_chromosomes"
                ],
                input_unique_variant_ids=len(identifiers),
                reference_unique_variant_ids=reference_metrics[
                    "reference_variants"
                ],
                exact_ids_shared=reference_metrics["overlapping_variants"],
                input_ids_absent_from_reference=absent_from_reference,
                reference_ids_absent_from_input=reference_metrics[
                    "reference_variants_absent_from_input"
                ],
                compatible_allele_pairs=reference_metrics[
                    "compatible_allele_pairs"
                ],
                incompatible_allele_pairs=reference_metrics[
                    "incompatible_allele_pairs"
                ],
                bim_structure_validated=True,
                input_rewritten=False,
            )
        input_metrics.update({
            "analysis_input_variants": input_frame.height,
            "variant_ids_replaced": 0,
            "unresolved_variants_removed": 0,
            "direct_input_unmodified": True,
            "formatter_input_unmodified": False,
        })
    annotation_metrics = dict(
        (pipeline_resources.gene_annotation_metrics or {})
        if pipeline_resources is not None
        else {}
    )
    if module.method in {"fastbat_gene", "mbat_combo"}:
        if not annotation_metrics:
            stage_context = stage(
                "Validate gene-coordinate annotation",
                "validate_gcta_gene_list",
            )
            with stage_context as annotation_stage:
                annotation_metrics = _validate_gene_list(
                    annotation_file, module,
                    set(reference_metrics["reference_chromosomes"]),
                )
                annotation_stage.outcome(
                    "Validated the GCTA gene-coordinate annotation.",
                    fields=gcta_gene_reference_outcome_fields(
                        module,
                        annotation_file,
                        annotation_metrics,
                    ),
                    gene_coordinate_file=str(annotation_file),
                    unique_gene_identifiers=annotation_metrics["genes"],
                    gene_chromosomes=annotation_metrics[
                        "gene_chromosomes"
                    ],
                    shared_reference_chromosomes=annotation_metrics[
                        "shared_chromosomes"
                    ],
                    declared_genome_build=module.genome_build,
                    coordinate_structure_validated=True,
                )
        elif not (
            set(annotation_metrics["gene_chromosomes"])
            & set(reference_metrics["reference_chromosomes"])
        ):
            raise GctaGeneError(
                "Gene list and LD reference do not share any chromosome labels."
            )
        else:
            annotation_metrics["shared_chromosomes"] = sorted(
                set(annotation_metrics["gene_chromosomes"])
                & set(reference_metrics["reference_chromosomes"])
            )
        (
            annotation_file,
            excluded_gene_list,
            gene_scope_metrics,
            excluded_gene_ids,
        ) = _prepare_scoped_gene_list(
            source_annotation_file,
            output,
            dataset_id,
            module,
            analysis_scope,
            overwrite=configuration.run.overwrite,
        )
        scoped_gene_list = (
            annotation_file if annotation_file != source_annotation_file else None
        )
        annotation_metrics.update(gene_scope_metrics)
    elif module.method == "fastbat_set":
        if module.set_annotation.gmt_file is not None:
            source_gene_list = (
                pipeline_resources.gene_list
                if pipeline_resources is not None
                else _required_path(
                    module.gene_annotation.file,
                    "GCTA gene list required for GMT gene-to-variant conversion",
                )
            )
            gene_analysis_scope = _analysis_scope(
                configuration,
                module.set_annotation.conversion.chromosome_label_policy,
            )
            (
                scoped_gene_list,
                excluded_gene_list,
                gene_scope_metrics,
                excluded_gene_ids,
            ) = _prepare_scoped_gene_list(
                source_gene_list,
                output,
                dataset_id,
                module,
                gene_analysis_scope,
                overwrite=configuration.run.overwrite,
            )
            if scoped_gene_list == source_gene_list:
                scoped_gene_list = None
            if pipeline_mode:
                (
                    annotation_file,
                    prepared_manifest,
                    prepared_directory,
                ) = _prepare_gmt_set(
                    output,
                    dataset_id,
                    reference_prefix,
                    module,
                    configuration,
                    logger,
                    reference_metrics["reference_variants"],
                    reference_ids,
                    source,
                    pathway_gene_preflight=pathway_gene_preflight,
                    excluded_gene_ids=excluded_gene_ids,
                    analysis_scope=gene_analysis_scope,
                    pipeline_args=pipeline_args,
                )
            else:
                with stage(
                    "Validate GMT and gene resources and prepare fastBAT sets",
                    "prepare_fastbat_set_source",
                ) as set_source_stage:
                    (
                        annotation_file,
                        prepared_manifest,
                        prepared_directory,
                    ) = _prepare_gmt_set(
                        output,
                        dataset_id,
                        reference_prefix,
                        module,
                        configuration,
                        logger,
                        reference_metrics["reference_variants"],
                        reference_ids,
                        source,
                        excluded_gene_ids=excluded_gene_ids,
                        analysis_scope=gene_analysis_scope,
                    )
                    set_source_stage.outcome(
                        "Validated the pathway and gene resources and published "
                        "the BIM-compatible fastBAT set list.",
                        fields=gcta_gmt_preparation_outcome_fields(
                            module,
                            module.set_annotation.gmt_file,
                            module.gene_annotation.file,
                            annotation_file,
                            prepared_manifest,
                        ),
                        gmt_file=str(module.set_annotation.gmt_file),
                        gene_coordinate_file=str(module.gene_annotation.file),
                        input_pathways=prepared_manifest["validation"][
                            "input_pathways"
                        ],
                        input_unique_genes=prepared_manifest["validation"][
                            "input_unique_genes"
                        ],
                        gene_coordinate_reference_genes=prepared_manifest[
                            "validation"
                        ]["gene_coordinate_reference_genes"],
                        matched_gene_identifiers=prepared_manifest["validation"][
                            "genes_with_coordinates"
                        ],
                        gmt_gene_mappability_fraction=prepared_manifest[
                            "validation"
                        ]["gmt_gene_mappability_fraction"],
                        coordinate_reference_coverage_fraction=prepared_manifest[
                            "validation"
                        ]["coordinate_reference_coverage_fraction"],
                        coordinate_reference_genes_absent_from_gmt=(
                            prepared_manifest["validation"][
                                "coordinate_reference_genes_absent_from_gmt"
                            ]
                        ),
                        pathways_with_complete_gene_id_mapping=(
                            prepared_manifest["validation"][
                                "pathways_with_complete_gene_id_mapping"
                            ]
                        ),
                        pathways_with_partial_gene_id_mapping=(
                            prepared_manifest["validation"][
                                "pathways_with_partial_gene_id_mapping"
                            ]
                        ),
                        pathways_without_gene_id_mapping=prepared_manifest[
                            "validation"
                        ]["pathways_without_gene_id_mapping"],
                        pathway_structure_validated=True,
                        gene_coordinate_structure_validated=True,
                        gene_id_compatibility_validated=True,
                        final_fastbat_sets=prepared_manifest["validation"][
                            "pathways_written"
                        ],
                        final_set_file=str(annotation_file),
                        published_set_resource_validated=True,
                    )
            validation = prepared_manifest["validation"]
            annotation_metrics = {
                "sets": validation["pathways_written"],
                "input_sets": validation["input_pathways"],
                "omitted_empty_sets": validation["pathways_omitted_empty"],
                "omitted_oversized_sets": validation[
                    "pathways_omitted_oversized"
                ],
                "requested_set_variants": validation[
                    "total_pathway_variant_memberships"
                ],
                "matched_set_variants": validation[
                    "total_pathway_variant_memberships"
                ],
                "analysis_set_variants": validation[
                    "total_pathway_variant_memberships"
                ],
                "unmatched_set_variants": 0,
                "maximum_set_variants": module.set_annotation.maximum_set_variants,
                "analysis_set_list": str(annotation_file),
                "input_unique_genes": validation["input_unique_genes"],
                "gene_coordinate_reference_genes": validation[
                    "gene_coordinate_reference_genes"
                ],
                "genes_with_coordinates": validation["genes_with_coordinates"],
                "genes_eligible_for_analysis": validation[
                    "genes_eligible_for_analysis"
                ],
                "genes_excluded_by_analysis_scope": validation[
                    "genes_excluded_by_analysis_scope"
                ],
                "unmapped_unique_genes": validation["unmapped_unique_genes"],
                "gmt_gene_mappability_fraction": validation[
                    "gmt_gene_mappability_fraction"
                ],
                "coordinate_reference_coverage_fraction": validation[
                    "coordinate_reference_coverage_fraction"
                ],
                "coordinate_reference_genes_absent_from_gmt": validation[
                    "coordinate_reference_genes_absent_from_gmt"
                ],
                "pathways_with_complete_gene_id_mapping": validation[
                    "pathways_with_complete_gene_id_mapping"
                ],
                "pathways_with_partial_gene_id_mapping": validation[
                    "pathways_with_partial_gene_id_mapping"
                ],
                "pathways_without_gene_id_mapping": validation[
                    "pathways_without_gene_id_mapping"
                ],
                **gene_scope_metrics,
            }
            analysis_set_path = annotation_file
        else:
            annotation_file = (
                pipeline_resources.set_list
                if pipeline_resources is not None
                else _required_path(
                    module.set_annotation.file, "GCTA fastBAT set list",
                )
            )
            source_set_path = annotation_file
            with stage(
                "Validate and prepare the fastBAT set-list input",
                "prepare_fastbat_analysis_sets",
                pipeline_key="set_policies" if pipeline_mode else None,
            ) as set_stage:
                analysis_set_path = configured_output_path(
                    output,
                    layout.analysis_set_list,
                    error_type=GctaGeneError,
                    dataset_id=dataset_id,
                )
                annotation_file, annotation_metrics = _prepare_analysis_set_list(
                    annotation_file,
                    analysis_set_path,
                    identifiers,
                    reference_ids,
                    module.set_annotation.conversion.empty_pathway_policy,
                    module.set_annotation.oversized_set_policy,
                    module.set_annotation.maximum_set_variants,
                )
                source_validation_fields = (
                    ()
                    if pipeline_mode
                    else tuple(gcta_set_source_outcome_fields(
                        source_set_path,
                        annotation_metrics,
                    ))
                )
                set_stage.outcome(
                    "Validated the supplied set list and wrote final fastBAT "
                    "sets from analyzable variants.",
                    fields=(
                        *source_validation_fields,
                        ("analysis", "GWAS/BIM set-membership compatibility"),
                        (
                            "success",
                            "Variant memberships shared by GWAS and BIM",
                            annotation_metrics["matched_set_variants"],
                        ),
                        (
                            "warning"
                            if annotation_metrics["unmatched_set_variants"]
                            else "success",
                            "Unavailable variant memberships",
                            annotation_metrics["unmatched_set_variants"],
                        ),
                        ("analysis", "Final fastBAT set list"),
                        ("success", "Final sets", annotation_metrics["sets"]),
                        (
                            "warning"
                            if annotation_metrics["omitted_empty_sets"]
                            else "success",
                            "Empty sets omitted",
                            annotation_metrics["omitted_empty_sets"],
                        ),
                        (
                            "warning"
                            if annotation_metrics["omitted_oversized_sets"]
                            else "success",
                            "Oversized sets omitted",
                            annotation_metrics["omitted_oversized_sets"],
                        ),
                        (
                            "count", "Retained variant memberships",
                            annotation_metrics["analysis_set_variants"],
                        ),
                        (
                            "success", "Final set-list file",
                            str(analysis_set_path),
                        ),
                    ),
                    input_set_file=str(source_set_path),
                    input_sets=annotation_metrics["input_sets"],
                    requested_variant_memberships=annotation_metrics[
                        "requested_set_variants"
                    ],
                    unique_requested_variants=annotation_metrics[
                        "unique_requested_set_variants"
                    ],
                    matched_variant_memberships=annotation_metrics[
                        "matched_set_variants"
                    ],
                    unavailable_variant_memberships=annotation_metrics[
                        "unmatched_set_variants"
                    ],
                    retained_sets=annotation_metrics["sets"],
                    omitted_empty_sets=annotation_metrics[
                        "omitted_empty_sets"
                    ],
                    omitted_oversized_sets=annotation_metrics[
                        "omitted_oversized_sets"
                    ],
                    set_block_structure_validated=True,
                    output=str(analysis_set_path),
                )
    variant_exclusion_file = _prepare_variant_exclusion_list(
        output,
        dataset_id,
        module,
        excluded_reference_ids,
        overwrite=configuration.run.overwrite,
    )
    del reference_ids
    logger.record("OBSERVED", "gcta_input", **input_metrics)
    logger.record("OBSERVED", "gcta_reference", **reference_metrics)
    if annotation_metrics:
        logger.record("OBSERVED", "gcta_annotation", **annotation_metrics)
    logger.record(
        "TRANSFORM",
        "gcta_analysis_scope",
        **analysis_scope,
        excluded_reference_variants=reference_metrics[
            "excluded_reference_variants"
        ],
        excluded_gwas_bim_variants=reference_metrics[
            "excluded_input_variants"
        ],
        retained_gwas_bim_variants=reference_metrics["analyzable_variants"],
        variant_exclusion_file=(
            str(variant_exclusion_file) if variant_exclusion_file else None
        ),
        scoped_gene_list=(str(scoped_gene_list) if scoped_gene_list else None),
        excluded_gene_list=(
            str(excluded_gene_list) if excluded_gene_list else None
        ),
        gene_scope=gene_scope_metrics or None,
    )
    if analysis_set_path is not None:
        logger.record(
            "TRANSFORM",
            "fastbat_set_intersection",
            input_sets=annotation_metrics["input_sets"],
            retained_sets=annotation_metrics["sets"],
            omitted_empty_sets=annotation_metrics["omitted_empty_sets"],
            requested_variant_memberships=(
                annotation_metrics["requested_set_variants"]
            ),
            retained_variant_memberships=(
                annotation_metrics["analysis_set_variants"]
            ),
            policy=module.set_annotation.conversion.empty_pathway_policy,
            oversized_set_policy=module.set_annotation.oversized_set_policy,
            maximum_set_variants=annotation_metrics["maximum_set_variants"],
            output=str(analysis_set_path),
        )

    final_prefix = configured_output_path(
        output, layout.output_prefix, error_type=GctaGeneError, dataset_id=dataset_id,
    )
    final_result = _method_result_path(output, dataset_id, module)
    normalized = configured_output_path(
        output, layout.normalized_result, error_type=GctaGeneError,
        dataset_id=dataset_id, method=module.method,
    )
    completion_manifest = _completion_manifest_path(
        output, dataset_id, module,
    )
    html_report = configured_output_path(
        output, layout.html_report, error_type=GctaGeneError,
        dataset_id=dataset_id, method=module.method,
    )
    final_prefix.parent.mkdir(parents=True, exist_ok=True)
    resumed = False
    report_current = False
    normalized_backup = None
    pending_publication = None
    if (
        final_result.exists()
        and configuration.run.resume
        and not configuration.run.overwrite
    ):
        with stage(
            "Validate GCTA and reuse the completed analysis",
            "resume_gcta_results",
            pipeline_key="run" if pipeline_mode else None,
        ) as run_stage:
            if pipeline_resources is not None:
                executable = pipeline_resources.executable
                version = pipeline_resources.version
            else:
                executable = resolve_executable(
                    configuration.resources.executables.gcta,
                    "GCTA executable",
                    error_type=GctaGeneError,
                )
                version = require_supported_gcta(
                    executable, module, logger,
                    configuration.execution.timeout_seconds,
                )
            result_metrics, report_current = _validate_gcta_completion_manifest(
                completion_manifest,
                dataset_id=dataset_id,
                version=version,
                configuration=configuration,
                source=source,
                source_annotation_file=source_annotation_file,
                annotation_file=annotation_file,
                variant_exclusion_file=variant_exclusion_file,
                reference_prefix=reference_prefix,
                final_result=final_result,
                normalized=normalized,
                html_report=html_report,
            )
            logger.record(
                "SKIP", "gcta_execution", reason="validated_resume",
                output=str(final_result), manifest=str(completion_manifest),
            )
            run_stage.outcome(
                "Reused the checksum-validated GCTA execution.",
                fields=(
                    ("analysis", "Method", module.method),
                    ("info", "GCTA version", version),
                    ("success", "Native GCTA execution", "validated resume"),
                    ("success", "Raw result", final_result.name),
                ),
            )
        with stage(
            "Validate the raw GCTA results",
            "validate_raw_gcta_results",
            pipeline_key="raw_results" if pipeline_mode else None,
        ) as raw_stage:
            validated_results = validate_raw_gcta_results(final_result, module)
            raw_stage.outcome(
                "Validated the raw GCTA result schema and scientific values.",
                fields=(
                    ("count", "%s tested" % validated_results.unit_label.capitalize(), validated_results.frame.height),
                    ("info", "P-value column", validated_results.p_column),
                    ("success", "Required result columns", "validated"),
                    ("success", "Unit identifiers", "non-empty and unique"),
                    ("success", "P-value range", "0 to 1"),
                ),
            )
        with stage(
            "Validate multiple-testing corrections",
            "validate_gcta_multiple_testing",
            pipeline_key="corrections" if pipeline_mode else None,
        ) as correction_stage:
            recorded_correction = result_metrics.get("multiple_testing", {}).get(
                "configuration"
            )
            current_correction = multiple_testing_configuration(module)
            if recorded_correction != current_correction:
                backup_handle = tempfile.NamedTemporaryFile(
                    dir=normalized.parent,
                    prefix=".%s.pre_correction." % normalized.name,
                    delete=False,
                )
                backup = Path(backup_handle.name)
                backup_handle.close()
                shutil.copy2(normalized, backup)
                try:
                    result_metrics = normalize_gcta_results(
                        final_result,
                        normalized,
                        module,
                    )
                    report_current = False
                except BaseException:
                    backup.replace(normalized)
                    raise
                normalized_backup = backup
                logger.record(
                    "TRANSFORM",
                    "gcta_multiple_testing_corrections",
                    source=str(final_result),
                    output=str(normalized),
                    **result_metrics["multiple_testing"],
                )
            correction_metrics = result_metrics["multiple_testing"]
            correction_stage.outcome(
                "Validated nominal, Bonferroni and Benjamini-Hochberg results.",
                fields=(
                    ("count", "Correction family size", correction_metrics["family_size"]),
                    ("count", "Nominally significant", correction_metrics["nominal_significant"]),
                    ("count", "Bonferroni significant", correction_metrics["bonferroni_significant"]),
                    ("count", "Benjamini-Hochberg FDR significant", correction_metrics["fdr_bh_significant"]),
                ),
            )
            resumed = True
    else:
        existing = [
            path for path in (final_result, normalized, html_report)
            if path.exists()
        ]
        if existing and not configuration.run.overwrite:
            raise GctaGeneError(
                "Output already exists: %s. Use --resume or --overwrite."
                % ", ".join(str(path) for path in existing)
            )
        staging = configured_output_path(
            output, layout.staging_directory, error_type=GctaGeneError,
            dataset_id=dataset_id, method=module.method,
        )
        if staging.exists():
            if not configuration.run.overwrite:
                raise GctaGeneError(
                    "Incomplete staging output exists: %s. Use --overwrite after review."
                    % staging
                )
            shutil.rmtree(staging)
        staging.mkdir(parents=True)
        staged_prefix = staging / final_prefix.name
        staged_result = _staged_result_path(final_prefix, final_result, staged_prefix)
        with stage(
            "Run GCTA %s" % module.method,
            "run_gcta_command",
            pipeline_key="run" if pipeline_mode else None,
        ) as step:
            if pipeline_resources is not None:
                executable = pipeline_resources.executable
                version = pipeline_resources.version
            else:
                executable = resolve_executable(
                    configuration.resources.executables.gcta,
                    "GCTA executable",
                    error_type=GctaGeneError,
                )
                version = require_supported_gcta(
                    executable, module, logger,
                    configuration.execution.timeout_seconds,
                )
            command = build_gcta_command(
                executable,
                source,
                reference_prefix,
                annotation_file,
                staged_prefix,
                configuration.execution.threads,
                module,
                exclude_variants_file=variant_exclusion_file,
            )
            step.input(
                "summary_statistics",
                path=str(source),
                rows=input_metrics["analysis_input_variants"],
            )
            measured_progress = None
            if not dry_run:
                result_schema = module.results.schemas[module.method]
                measured_progress = GctaResultProgress(
                    method=module.method,
                    unit_label=result_schema.unit_label,
                    native_log=gcta_native_log_path(staged_prefix),
                    result_file=staged_result,
                    logger=logger,
                    enabled=configuration.logging.show_progress,
                    total_hint=(
                        annotation_metrics.get("sets")
                        if module.method == "fastbat_set" else None
                    ),
                )
                measured_progress.start()
            try:
                run_gcta_command(
                    command,
                    staged_result,
                    configuration,
                    logger,
                    dry_run=dry_run,
                    progress_callback=(
                        measured_progress.refresh
                        if measured_progress is not None
                        and measured_progress.enabled else None
                    ),
                )
            except BaseException:
                if measured_progress is not None:
                    measured_progress.fail()
                raise
            if dry_run:
                step.output("validated_command", output_prefix=str(final_prefix))
                step.outcome(
                    "Validated the exact GCTA command without executing it.",
                    fields=(
                        ("analysis", "Method", module.method),
                        ("info", "GCTA version", version),
                        ("success", "Command validation", "passed"),
                    ),
                )
                shutil.rmtree(staging)
                dry_artifacts = {}
                if prepared_directory is not None:
                    dry_artifacts["prepared_set_resource"] = Artifact(
                        "gcta_fastbat_set_resource",
                        prepared_directory,
                        {
                            "dataset_id": dataset_id,
                            "method": module.method,
                            "genome_build": module.genome_build,
                            "source": "gmt_conversion",
                        },
                    )
                if analysis_set_path is not None:
                    dry_artifacts["analysis_set_list"] = Artifact(
                        "gcta_fastbat_analysis_set",
                        analysis_set_path,
                        {
                            "dataset_id": dataset_id,
                            "method": module.method,
                            "genome_build": module.genome_build,
                            "source": "gwas_bim_intersection",
                        },
                    )
                if variant_exclusion_file is not None:
                    dry_artifacts["variant_exclusion_list"] = Artifact(
                        "gcta_variant_exclusion_list",
                        variant_exclusion_file,
                        {
                            "dataset_id": dataset_id,
                            "method": module.method,
                            "genome_build": module.genome_build,
                            "source": "configured_analysis_scope",
                        },
                    )
                if scoped_gene_list is not None:
                    dry_artifacts["scoped_gene_list"] = Artifact(
                        "gcta_scoped_gene_list",
                        scoped_gene_list,
                        {
                            "dataset_id": dataset_id,
                            "method": module.method,
                            "genome_build": module.genome_build,
                            "source": "configured_analysis_scope",
                        },
                    )
                if excluded_gene_list is not None:
                    dry_artifacts["excluded_gene_list"] = Artifact(
                        "gcta_excluded_gene_list",
                        excluded_gene_list,
                        {
                            "dataset_id": dataset_id,
                            "method": module.method,
                            "genome_build": module.genome_build,
                            "source": "configured_analysis_scope",
                        },
                    )
                return ModuleResult(
                    "gcta_gene",
                    artifacts=dry_artifacts,
                    metrics={
                        "method": module.method,
                        "dry_run": True,
                        "gcta_version": version,
                        **input_metrics,
                        **reference_metrics,
                        **annotation_metrics,
                    },
                )
            step.outcome(
                "GCTA completed and produced the expected native result file.",
                fields=(
                    ("analysis", "Method", module.method),
                    ("info", "GCTA version", version),
                    ("count", "Input variants", input_metrics["analysis_input_variants"]),
                    (
                        "count",
                        "Variants eligible after analysis-scope exclusions",
                        reference_metrics["analyzable_variants"],
                    ),
                    ("success", "Native result created", staged_result.name),
                ),
            )
        with stage(
            "Validate the raw GCTA results",
            "validate_raw_gcta_results",
            pipeline_key="raw_results" if pipeline_mode else None,
        ) as raw_stage:
            try:
                validated_results = validate_raw_gcta_results(
                    staged_result, module,
                )
            except BaseException:
                if measured_progress is not None:
                    measured_progress.fail()
                raise
            if measured_progress is not None:
                measured_progress.complete(validated_results.frame.height)
            raw_stage.outcome(
                "Validated the raw GCTA result schema and scientific values.",
                fields=(
                    ("count", "%s tested" % validated_results.unit_label.capitalize(), validated_results.frame.height),
                    ("info", "P-value column", validated_results.p_column),
                    ("success", "Required result columns", "validated"),
                    ("success", "Unit identifiers", "non-empty and unique"),
                    ("success", "P-value range", "0 to 1"),
                ),
            )
        stage_normalized = staging / "normalized" / normalized.name
        with stage(
            "Add and validate multiple-testing corrections",
            "add_multiple_testing_results",
            pipeline_key="corrections" if pipeline_mode else None,
        ) as correction_stage:
            result_metrics = add_multiple_testing_results(
                validated_results, stage_normalized, module,
            )
            correction_metrics = result_metrics["multiple_testing"]
            correction_stage.outcome(
                "Added nominal, Bonferroni and Benjamini-Hochberg results.",
                fields=(
                    ("count", "Correction family size", correction_metrics["family_size"]),
                    ("count", "Nominally significant", correction_metrics["nominal_significant"]),
                    ("count", "Bonferroni significant", correction_metrics["bonferroni_significant"]),
                    ("count", "Benjamini-Hochberg FDR significant", correction_metrics["fdr_bh_significant"]),
                    ("success", "Normalized result", stage_normalized.name),
                ),
            )

        pending_publication = (staging, staged_prefix, stage_normalized)

    with stage(
        "Publish validated GCTA outputs and build the scientific summary",
        "publish_gcta_results",
        pipeline_key="publish" if pipeline_mode else None,
    ) as publish_stage:
        published = []
        if pending_publication is not None:
            staging, staged_prefix, stage_normalized = pending_publication
            published = _publish_staged_outputs(
                staged_prefix, final_prefix, configuration.run.overwrite,
            )
            normalized.parent.mkdir(parents=True, exist_ok=True)
            if normalized.exists() and not configuration.run.overwrite:
                raise GctaGeneError("Output already exists: %s" % normalized)
            stage_normalized.replace(normalized)
            result_metrics["normalized_result"] = str(normalized)
            publish_stage.set_rows(result_metrics["tested_units"], removed=0)
            publish_stage.output("gcta_results", files=[str(path) for path in published])
            publish_stage.output("normalized_results", path=str(normalized))
        if not final_result.is_file() or not normalized.is_file():
            raise GctaGeneError(
                "Validated GCTA output publication is incomplete."
            )
        metrics = {
            "method": module.method,
            "dry_run": False,
            "resumed": resumed,
            "gcta_version": version,
            **input_metrics,
            **reference_metrics,
            **annotation_metrics,
            **result_metrics,
        }
        metrics["scientific_summary"] = build_gcta_scientific_summary(
            normalized,
            module,
            metrics,
        )
        report_outputs = {
            "Complete normalized results": normalized,
            "Original GCTA results": final_result,
            "Completion manifest": completion_manifest,
            "Canonical GCTA log": configured_output_path(
                output, layout.log_file, error_type=GctaGeneError,
                dataset_id=dataset_id, method=module.method,
            ),
        }
        frequency_qc = configured_output_path(
            output, layout.frequency_qc_result, error_type=GctaGeneError,
            dataset_id=dataset_id,
        )
        if frequency_qc.is_file():
            report_outputs["GCTA frequency-QC result"] = frequency_qc
        snpset = configured_output_path(
            output, layout.mbat_snpset_result, error_type=GctaGeneError,
            dataset_id=dataset_id,
        )
        if snpset.is_file():
            report_outputs["GCTA mBAT SNP-set result"] = snpset
        report_backup = None
        report_existed = html_report.is_file()
        if resumed and not report_current and report_existed:
            report_backup_handle = tempfile.NamedTemporaryFile(
                dir=html_report.parent,
                prefix=".%s.pre_refresh." % html_report.name,
                delete=False,
            )
            report_backup = Path(report_backup_handle.name)
            report_backup_handle.close()
            shutil.copy2(html_report, report_backup)
        try:
            if not report_current:
                write_gcta_html_report(
                    normalized,
                    html_report,
                    dataset_id=dataset_id,
                    module_config=module,
                    summary=metrics["scientific_summary"],
                    gcta_version=version,
                    output_files=report_outputs,
                )
                logger.record(
                    "OUTPUT", "gcta_html_report", path=str(html_report),
                    rows=result_metrics["tested_units"],
                    page_size=module.html_report.page_size,
                    method=module.method,
                )
            else:
                logger.record(
                    "SKIP", "gcta_html_report", reason="validated_resume",
                    path=str(html_report),
                )
            _write_gcta_completion_manifest(
                completion_manifest,
                dataset_id=dataset_id,
                version=version,
                configuration=configuration,
                source=source,
                source_annotation_file=source_annotation_file,
                annotation_file=annotation_file,
                variant_exclusion_file=variant_exclusion_file,
                reference_prefix=reference_prefix,
                final_result=final_result,
                normalized=normalized,
                html_report=html_report,
                result_metrics=result_metrics,
            )
        except BaseException:
            if normalized_backup is not None:
                normalized_backup.replace(normalized)
            if report_backup is not None:
                report_backup.replace(html_report)
            elif resumed and not report_current and not report_existed:
                html_report.unlink(missing_ok=True)
            raise
        finally:
            if normalized_backup is not None:
                normalized_backup.unlink(missing_ok=True)
            if report_backup is not None:
                report_backup.unlink(missing_ok=True)
        if pending_publication is not None:
            shutil.rmtree(staging)
        publish_stage.output("html_report", path=str(html_report))
        publish_stage.outcome(
            "Validated and published all outputs and built the scientific summary.",
            fields=(
                (
                    "success", "Published native GCTA files",
                    len(published) if pending_publication is not None else "validated resume",
                ),
                ("success", "Normalized results", normalized.name),
                ("success", "HTML report", html_report.name),
                ("success", "Completion manifest", completion_manifest.name),
                (
                    "count", "%s tested" % result_metrics["unit_label"].capitalize(),
                    result_metrics["tested_units"],
                ),
                (
                    "count", "Nominally significant",
                    result_metrics["multiple_testing"]["nominal_significant"],
                ),
                (
                    "count", "Bonferroni significant",
                    result_metrics["multiple_testing"]["bonferroni_significant"],
                ),
                (
                    "count", "Benjamini-Hochberg FDR significant",
                    result_metrics["multiple_testing"]["fdr_bh_significant"],
                ),
            ),
        )

    metadata = {
        "dataset_id": dataset_id,
        "method": module.method,
        "genome_build": module.genome_build,
        "reference_population": module.reference.population,
        "effect_allele": "A1=ALT",
        "schema_version": module.results.normalized_schema_version,
        "gcta_version": version,
        "mhc_policy": analysis_scope["mhc_policy"],
        "excluded_chromosomes": analysis_scope["exclude_chromosomes"],
    }
    artifacts = {
        "raw_results": Artifact("gcta_gene_raw", final_result, metadata),
        "normalized_results": Artifact("gcta_gene_results", normalized, metadata),
        "html_report": Artifact("gcta_gene_html_report", html_report, metadata),
        "completion_manifest": Artifact(
            "gcta_gene_completion_manifest", completion_manifest, metadata,
        ),
    }
    if prepared_directory is not None:
        artifacts["prepared_set_resource"] = Artifact(
            "gcta_fastbat_set_resource",
            prepared_directory,
            {
                **metadata,
                "source": "gmt_conversion",
                "pathways_written": prepared_manifest["validation"]["pathways_written"],
            },
        )
    if analysis_set_path is not None:
        artifacts["analysis_set_list"] = Artifact(
            "gcta_fastbat_analysis_set",
            analysis_set_path,
            {
                **metadata,
                "source": "gwas_bim_intersection",
                "sets_retained": annotation_metrics["sets"],
                "sets_omitted": annotation_metrics["omitted_empty_sets"],
            },
        )
    if variant_exclusion_file is not None:
        artifacts["variant_exclusion_list"] = Artifact(
            "gcta_variant_exclusion_list", variant_exclusion_file, metadata,
        )
    if scoped_gene_list is not None:
        artifacts["scoped_gene_list"] = Artifact(
            "gcta_scoped_gene_list", scoped_gene_list, metadata,
        )
    if excluded_gene_list is not None:
        artifacts["excluded_gene_list"] = Artifact(
            "gcta_excluded_gene_list", excluded_gene_list, metadata,
        )
    frequency_qc = configured_output_path(
        output, layout.frequency_qc_result, error_type=GctaGeneError,
        dataset_id=dataset_id,
    )
    if frequency_qc.is_file():
        artifacts["frequency_qc"] = Artifact("gcta_frequency_qc", frequency_qc, metadata)
    snpset = configured_output_path(
        output, layout.mbat_snpset_result, error_type=GctaGeneError,
        dataset_id=dataset_id,
    )
    if snpset.is_file():
        artifacts["snp_set"] = Artifact("gcta_gene_snpset", snpset, metadata)
    if not pipeline_mode and stage_number != stage_total:
        raise GctaGeneError(
            "Internal GCTA progress plan mismatch: completed "
            "%d of %d planned stages." % (stage_number, stage_total)
        )
    return ModuleResult(
        "gcta_gene",
        artifacts=artifacts,
        metrics=metrics,
    )


def run_gcta_gene_direct(
    args,
    ctx=None,
    *,
    pipeline_resources: GctaGenePipelineResources | None = None,
) -> ModuleResult:
    """Resolve configuration once, execute the selected test, and finalize logging."""
    try:
        configuration = (
            _gcta_gene_pipeline_execution_configuration(args, pipeline_resources)
            if pipeline_resources is not None
            else _resolved_configuration(args)
        )
    except BaseException as exc:
        fallback = load_configuration()
        output = Path(
            getattr(args, "output_directory", None) or fallback.run.output_directory
        ).expanduser().resolve()
        dataset = str(getattr(args, "dataset_id", None) or fallback.run.dataset_id)
        method = str(
            getattr(args, "gcta_gene_method", None)
            or fallback.modules.gcta_gene.method
        )
        log_path = configured_output_path(
            output,
            fallback.modules.gcta_gene.output_layout.log_file,
            error_type=GctaGeneError,
            dataset_id=dataset,
            method=method,
        )
        write_log_record(
            log_path,
            "ERROR",
            "GCTA gene configuration failed: %s: %s" % (type(exc).__name__, exc),
            sample_id=dataset,
            file_level=fallback.logging.file_level,
            screen_level=fallback.logging.console_level,
        )
        raise
    module = configuration.modules.gcta_gene
    if pipeline_resources is not None:
        require_unchanged_preflight_files(
            pipeline_resources.file_identities,
            error_type=GctaGeneError,
            label="GCTA gene resource",
        )
    output = Path(
        getattr(args, "output_directory", None) or configuration.run.output_directory
    ).expanduser().resolve()
    dataset = str(
        getattr(args, "dataset_id", None) or configuration.run.dataset_id
    ).strip()
    log_path = configured_output_path(
        output,
        module.output_layout.log_file,
        error_type=GctaGeneError,
        dataset_id=dataset,
        method=module.method,
    )
    logger = PipelineLogger(
        dataset,
        "run",
        str(log_path.parent),
        level=configuration.logging.file_level,
        screen_level=configuration.logging.console_level,
        log_path=str(log_path),
        stage_progress=(
            None
            if getattr(args, "_pipeline_stage_progress", None) is not None
            else StageProgress(
                "GCTA %s stages" % module.method.replace("_", " "),
                enabled=configuration.logging.show_progress,
            )
        ),
    )
    if ctx is not None:
        configure_pipeline_stage_callbacks(args, configuration, log_path)
    try:
        formatter_output = (
            ctx.get("formatter", {}).get("gcta_gene", {})
            if ctx is not None else {}
        )
        input_file = getattr(args, "gcta_input_file", None) or module.input_file
        input_key = "summary_statistics_input_file"
        if input_file is None and ctx is not None:
            input_file = formatter_output.get(input_key)
        pipeline_mode = ctx is not None
        supplied_vcf = getattr(args, "vcf", None)
        if not pipeline_mode and supplied_vcf is not None:
            raise GctaGeneError(
                "Direct gcta_gene does not accept --vcf because direct input "
                "mode compares exact SNP IDs without coordinate-based "
                "rewriting. Remove --vcf; summary-statistic IDs absent from "
                "the BIM will be reported and ignored by GCTA. Use postgwas "
                "pipeline --modules gcta_gene --vcf PATH when starting from "
                "a harmonised GWAS-VCF; the formatter will select the "
                "BIM-compatible identifier format before creating the .ma."
            )
        formatter_variant_id_type = (
            formatter_output.get("variant_id_type")
            if pipeline_mode else None
        )
        if pipeline_mode and not formatter_variant_id_type:
            raise GctaGeneError(
                "Pipeline formatter output is missing variant_id_type provenance; "
                "the GCTA stage cannot validate the BIM-compatible identifier "
                "contract safely. Re-run the formatter stage."
            )
        _require_gcta_gene_arguments(configuration, input_file)
        resource_paths = ["executables.gcta"]
        if module.genome_build is not None:
            resource_paths.append("genomes.%s" % module.genome_build)
        if module.reference.population is not None:
            resource_paths.append(
                "populations.%s" % module.reference.population
            )
        write_resolved_configuration(
            configuration,
            configured_output_path(
                output,
                module.output_layout.resolved_config_file,
                error_type=GctaGeneError,
                dataset_id=dataset,
                method=module.method,
            ),
            modules=("formatting", "gcta_gene"),
            resource_paths=resource_paths,
        )
        logger.record(
            "PARAM",
            "gcta_gene",
            method=module.method,
            genome_build=module.genome_build,
            reference_population=module.reference.population,
            gene_window_kb=(
                module.gene_window_kb
                if module.method in {"fastbat_gene", "mbat_combo"}
                or module.set_annotation.gmt_file is not None else None
            ),
            set_source=(
                "gmt_conversion"
                if module.set_annotation.gmt_file is not None
                else "prepared_set_list"
                if module.set_annotation.file is not None
                else None
            ),
            oversized_set_policy=(
                module.set_annotation.oversized_set_policy
                if module.method == "fastbat_set" else None
            ),
            gmt_conversion_policies=(
                module.set_annotation.conversion.model_dump(mode="json")
                if module.set_annotation.gmt_file is not None else None
            ),
            segment_size_kb=(
                module.segment_size_kb
                if module.method == "fastbat_segment" else None
            ),
            reference_maf_min=module.reference_maf_min,
            fastbat_ld_cutoff=module.fastbat_ld_cutoff,
            mbat_svd_gamma=(module.mbat_svd_gamma if module.method == "mbat_combo" else None),
            frequency_difference_max=module.frequency_difference_max,
            print_component_p_values=(
                module.print_component_p_values
                if module.method == "mbat_combo" else None
            ),
            write_snpset=module.write_snpset,
            reporting=module.reporting.model_dump(mode="json"),
            html_report=module.html_report.model_dump(mode="json"),
            variant_harmonisation=module.variant_harmonisation.model_dump(
                mode="json"
            ),
            harmonised_vcf=(str(supplied_vcf) if pipeline_mode else None),
            formatter_variant_id_type=formatter_variant_id_type,
            variant_id_match_policy=(
                "validate_pipeline_formatter_bim_compatibility_without_rewriting"
                if pipeline_mode else "report_exact_bim_overlap_without_rewriting"
            ),
            threads=configuration.execution.threads,
        )
        if pipeline_resources is not None:
            logger.record(
                "PASS",
                "gcta_gene_pipeline_resource_preflight",
                reference_prefix=str(pipeline_resources.reference_prefix),
                reference_files=[
                    str(path) for path in pipeline_resources.reference_paths
                ],
                gcta_executable=pipeline_resources.executable,
                gcta_version=pipeline_resources.version,
                reference_identifier_observation=dict(
                    pipeline_resources.identifier_observation
                ),
                gene_annotation_metrics=(
                    None
                    if pipeline_resources.gene_annotation_metrics is None
                    else dict(pipeline_resources.gene_annotation_metrics)
                ),
                set_source_metrics=(
                    None
                    if pipeline_resources.set_source_metrics is None
                    else dict(pipeline_resources.set_source_metrics)
                ),
                gmt_gene_mappability_fraction=(
                    None
                    if pipeline_resources.pathway_gene_preflight is None
                    else pipeline_resources.pathway_gene_preflight
                    .gmt_gene_mappability_fraction
                ),
                coordinate_reference_coverage_fraction=(
                    None
                    if pipeline_resources.pathway_gene_preflight is None
                    else pipeline_resources.pathway_gene_preflight
                    .coordinate_coverage_fraction
                ),
                coordinate_reference_genes_absent_from_gmt=(
                    None
                    if pipeline_resources.pathway_gene_preflight is None
                    else len(
                        pipeline_resources.pathway_gene_preflight
                        .coordinate_genes_absent_from_gmt
                    )
                ),
            )
        input_preflight = (
            ctx.validation("gcta_gene_input")
            if pipeline_mode and hasattr(ctx, "validation") else None
        )
        result = run_gcta_gene(
            input_file,
            output,
            dataset,
            configuration,
            logger,
            formatter_variant_id_type=formatter_variant_id_type,
            input_preflight=input_preflight,
            pipeline_resources=pipeline_resources,
            pipeline_args=args if pipeline_mode else None,
            dry_run=bool(getattr(args, "dry_run", False)),
        )
        if ctx is not None:
            ctx.publish(result)
        if not result.metrics.get("dry_run"):
            summary = result.metrics["scientific_summary"]
            record_gcta_scientific_summary(logger, summary)
        if result.metrics.get("resumed"):
            print(screen_line(
                "success",
                "GCTA outputs validated; continuing from completed step",
                indent=2,
            ))
        status = "VALIDATED" if result.metrics.get("dry_run") else "COMPLETED"
        logger.record(
            "STATUS", "gcta_gene", status=status, method=module.method,
            resumed=bool(result.metrics.get("resumed")),
            scientific_warnings=(
                len(result.metrics["scientific_summary"]["warnings"])
                if not result.metrics.get("dry_run") else None
            ),
        )
        if not result.metrics.get("dry_run"):
            print(render_gcta_scientific_summary(
                result.metrics["scientific_summary"],
                dataset,
                module,
                result.artifacts,
                log_path,
                configuration.logging.terminal_label_width,
            ))
            return result
        lines = [
            "",
            screen_line(
                "analysis",
                "GCTA %s %s" % (
                    module.method,
                    "command validated" if result.metrics.get("dry_run") else "completed",
                ),
                indent=2,
            ),
            screen_field(
                "info", "Dataset", dataset, indent=6,
                label_width=configuration.logging.terminal_label_width,
            ),
            screen_field(
                "info", "Full log", log_path, indent=6,
                label_width=configuration.logging.terminal_label_width,
            ),
        ]
        lines.append("")
        print("\n".join(lines))
        return result
    except BaseException as exc:
        if not logger.summary()["failed"]:
            logger.error("GCTA gene analysis failed: %s: %s" % (type(exc).__name__, exc))
        raise
    finally:
        logger.close()


__all__ = [
    "GctaGenePipelineResources",
    "GctaInputPreflight",
    "preflight_gcta_gene_pipeline",
    "run_gcta_gene",
    "run_gcta_gene_direct",
    "validate_gcta_reference_files",
    "validate_pipeline_gcta_input",
]
