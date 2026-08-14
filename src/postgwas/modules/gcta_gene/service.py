"""Validated execution service for GCTA fastBAT and mBAT-combo."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import shlex
import shutil
import sys
import tempfile
from types import SimpleNamespace

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
from postgwas.core.io.tables import write_dataframe_table
from postgwas.core.io.reports import write_yaml_report
from postgwas.core.paths import configured_output_path, resolve_executable
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.core.required_arguments import (
    RequiredAlternative,
    RequiredArgument,
    require_resolved_arguments,
)
from postgwas.core.resource_preparation import ResourcePreparationError, sha256
from postgwas.core.ui import StageProgress
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.core.vcf import extract_vcf_table
from postgwas.modules.gcta_gene.adapters import (
    build_gcta_command,
    require_supported_gcta,
    run_gcta_command,
)
from postgwas.modules.gcta_gene.errors import GctaGeneError
from postgwas.modules.gcta_gene.pathway_sets import prepare_resource
from postgwas.modules.gcta_gene.progress import (
    GctaResultProgress,
    gcta_native_log_path,
)
from postgwas.modules.gcta_gene.reporting import (
    build_gcta_scientific_summary,
    record_gcta_scientific_summary,
    render_gcta_scientific_summary,
)
from postgwas.modules.gcta_gene.results import (
    multiple_testing_configuration,
    normalize_gcta_results,
)


_DNA_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_COMPLETION_SCHEMA_VERSION = 1


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
        "gcta_coordinate_fallback": "variant_harmonisation.coordinate_fallback",
        "gcta_chromosome_label_policy": (
            "variant_harmonisation.chromosome_label_policy"
        ),
        "gcta_minimum_reference_overlap": (
            "variant_harmonisation.minimum_overlap_fraction"
        ),
        "gcta_allow_strand_complement": (
            "variant_harmonisation.allow_strand_complement"
        ),
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


def _require_gcta_gene_arguments(configuration, input_file) -> None:
    """Validate all requirements for the resolved GCTA method together."""
    module = configuration.modules.gcta_gene
    requirements = [
        RequiredArgument(
            "--gcta-input-file", "modules.gcta_gene.input_file", input_file,
        ),
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


def _allele_key(first: str, second: str) -> tuple[str, str]:
    alleles = tuple(sorted((str(first).strip().upper(), str(second).strip().upper())))
    if not all(alleles) or alleles[0] == alleles[1]:
        raise GctaGeneError("Variant matching requires two nonempty, distinct alleles.")
    return alleles


def _variant_locations_from_vcf(
    vcf_path: str | Path,
    identifiers: set[str],
    summary_alleles: dict[str, tuple[str, str]],
    output: Path,
    dataset_id: str,
    configuration,
    logger: PipelineLogger,
) -> tuple[dict[str, tuple[str, int, tuple[str, str]]], dict]:
    """Read the configured VCF projection needed for allele-aware BIM-ID matching."""
    vcf = _required_path(str(vcf_path), "Harmonised GWAS-VCF")
    formatting = configuration.modules.formatting
    policy = configuration.modules.gcta_gene.variant_harmonisation
    canonical = formatting.canonical_columns
    source_columns = [
        canonical.chromosome,
        canonical.position,
        canonical.variant_id,
        canonical.reference_allele,
        canonical.alternate_allele,
    ]
    projection = {
        column: formatting.vcf_fields.root[column] for column in source_columns
    }
    bcftools = resolve_executable(
        configuration.resources.executables.bcftools,
        "bcftools executable",
        error_type=GctaGeneError,
    )
    handle = tempfile.NamedTemporaryFile(
        mode="w",
        suffix=formatting.runtime.temporary_table_suffix,
        prefix=formatting.runtime.temporary_table_prefix.format(
            dataset_id=dataset_id,
        ),
        dir=output,
        delete=False,
    )
    work_table = Path(handle.name)
    handle.close()
    try:
        extract_vcf_table(
            vcf,
            work_table,
            dataset_id,
            projection,
            bcftools,
            delimiter=formatting.runtime.table_delimiter,
            io_buffer_bytes=formatting.runtime.io_buffer_bytes,
            include_expression=formatting.vcf_include_expression,
            allow_undefined_tags=True,
            logger=logger,
            error_type=GctaGeneError,
            purpose="Extracting GWAS variant coordinates for GCTA ID reconciliation",
        )
        try:
            frame = pl.read_csv(
                work_table,
                separator=formatting.runtime.table_delimiter,
                null_values=formatting.runtime.input_null_values,
                schema={column: pl.String for column in projection},
                low_memory=True,
            )
        except (OSError, pl.exceptions.PolarsError) as exc:
            raise GctaGeneError(
                "Cannot read the GCTA variant-reconciliation VCF table: %s" % exc
            ) from exc
    finally:
        work_table.unlink(missing_ok=True)

    resolved = canonical.resolved_variant_id
    frame = frame.with_columns(
        pl.when(
            pl.col(canonical.variant_id).is_null()
            | (pl.col(canonical.variant_id).str.strip_chars() == "")
        )
        .then(pl.concat_str([
            canonical.chromosome,
            canonical.position,
            canonical.reference_allele,
            canonical.alternate_allele,
        ], separator="_"))
        .otherwise(pl.col(canonical.variant_id).str.strip_chars())
        .alias(resolved),
        pl.col(canonical.position).cast(pl.Int64, strict=False),
        pl.col(canonical.reference_allele).str.strip_chars().str.to_uppercase(),
        pl.col(canonical.alternate_allele).str.strip_chars().str.to_uppercase(),
    ).filter(pl.col(resolved).is_in(list(identifiers)))
    duplicate_ids = frame.select(pl.col(resolved).is_duplicated().sum()).item()
    if duplicate_ids:
        raise GctaGeneError(
            "Harmonised GWAS-VCF assigns multiple records to %d formatted "
            "variant identifiers." % int(duplicate_ids)
        )
    found = set(str(value) for value in frame[resolved].to_list())
    missing = sorted(identifiers - found)
    if missing:
        raise GctaGeneError(
            "The supplied harmonised GWAS-VCF does not contain %d formatted "
            "GCTA variant IDs; examples: %s"
            % (len(missing), ", ".join(missing[:10]))
        )
    locations = {}
    for variant, chromosome, position, reference, alternate in frame.select(
        resolved,
        canonical.chromosome,
        canonical.position,
        canonical.reference_allele,
        canonical.alternate_allele,
    ).iter_rows():
        if position is None or int(position) < 1:
            raise GctaGeneError(
                "Harmonised GWAS-VCF has an invalid position for variant %s."
                % variant
            )
        key = (
            _normalise_chromosome(
                chromosome, policy.chromosome_label_policy,
            ),
            int(position),
            _allele_key(reference, alternate),
        )
        locations[str(variant)] = key
        if str(variant) in summary_alleles:
            if _allele_key(*summary_alleles[str(variant)]) != key[2]:
                raise GctaGeneError(
                    "GCTA .ma alleles disagree with the harmonised GWAS-VCF "
                    "for variant %s." % variant
                )
    return locations, {
        "variant_coordinate_source": str(vcf),
        "variant_coordinates": len(locations),
    }


def _validate_reference(
    prefix: Path,
    module,
    identifiers: set[str],
    alleles: dict[str, tuple[str, str]],
    variant_locations: dict[str, tuple[str, int, tuple[str, str]]] | None,
    *,
    retain_variant_ids: bool,
    direct_input_mode: bool,
) -> tuple[dict, set[str], dict[str, str]]:
    if direct_input_mode and variant_locations is not None:
        raise GctaGeneError(
            "Direct exact-ID validation cannot use coordinate reconciliation."
        )
    required_paths = [Path(str(prefix) + suffix) for suffix in module.reference.required_extensions]
    missing = [str(path) for path in required_paths if not path.is_file() or path.stat().st_size <= 0]
    if missing:
        raise GctaGeneError(
            "PLINK LD reference files are missing or empty: %s" % ", ".join(missing)
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
    direct_matches: dict[str, str] = {}
    coordinate_candidates: dict[str, str] = {}
    ambiguous_coordinate_sources: set[str] = set()
    direct_location_mismatches: list[str] = []
    policy = module.variant_harmonisation
    unique_input_locations = {}
    duplicate_input_locations = set()
    if variant_locations is not None:
        for source, key in variant_locations.items():
            previous = unique_input_locations.get(key)
            if previous is None and key not in duplicate_input_locations:
                unique_input_locations[key] = source
            elif previous != source:
                unique_input_locations.pop(key, None)
                duplicate_input_locations.add(key)
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
            reference_alleles = _allele_key(
                fields[column_index["allele1"]],
                fields[column_index["allele2"]],
            )
            reference_key = (
                normalized_chromosome, position, reference_alleles,
            )
            if variant in identifiers:
                if (
                    variant_locations is not None
                    and variant_locations[variant] != reference_key
                ):
                    direct_location_mismatches.append(variant)
                else:
                    direct_matches[variant] = variant
            if (
                variant_locations is not None
                and policy.coordinate_fallback
                and reference_key in unique_input_locations
            ):
                source = unique_input_locations[reference_key]
                previous = coordinate_candidates.get(source)
                if previous is None and source not in ambiguous_coordinate_sources:
                    coordinate_candidates[source] = variant
                elif previous != variant:
                    coordinate_candidates.pop(source, None)
                    ambiguous_coordinate_sources.add(source)
            if variant_locations is None and variant in alleles:
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
    if direct_location_mismatches:
        raise GctaGeneError(
            "%d direct ID matches have different GWAS-VCF and PLINK BIM "
            "coordinates or allele pairs; examples: %s"
            % (
                len(direct_location_mismatches),
                ", ".join(direct_location_mismatches[:10]),
            )
        )
    mapping = dict(direct_matches)
    target_owners = {target: source for source, target in mapping.items()}
    target_conflicts = 0
    for source, target in coordinate_candidates.items():
        if source in mapping or source in ambiguous_coordinate_sources:
            continue
        owner = target_owners.get(target)
        if owner is not None and owner != source:
            target_conflicts += 1
            continue
        mapping[source] = target
        target_owners[target] = source
    overlap = len(mapping)
    if overlap == 0:
        if direct_input_mode:
            raise GctaGeneError(
                "No direct GCTA summary-statistic SNP IDs occur exactly in "
                "PLINK BIM column 2. Direct mode does not rewrite identifiers. "
                "Provide --gcta-input-file containing BIM-compatible IDs, or "
                "use postgwas pipeline --modules gcta_gene --vcf PATH for "
                "coordinate-and-allele reconciliation."
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
            "Only %d/%d unique GCTA input variants (%.2f%%) resolved to exact "
            "PLINK BIM IDs; the configured minimum is %.2f%%. Supply the "
            "harmonised GWAS-VCF for coordinate-and-allele reconciliation or "
            "review the LD reference."
            % (
                overlap,
                len(identifiers),
                overlap_fraction * 100,
                policy.minimum_overlap_fraction * 100,
            )
        )
    if variant_locations is None and alleles and compatible_alleles == 0:
        raise GctaGeneError(
            "No overlapping GCTA .ma variants have compatible allele pairs in "
            "the PLINK BIM reference."
        )
    if incompatible_alleles:
        raise GctaGeneError(
            "%d overlapping GCTA .ma variants have allele pairs that cannot be "
            "matched to the PLINK BIM reference." % incompatible_alleles
        )
    if variant_locations is not None and alleles:
        compatible_alleles = overlap
    metrics = {
        "reference_variants": len(reference_ids),
        "overlapping_variants": overlap,
        "overlap_fraction": overlap_fraction,
        "input_variants_absent_from_reference": len(identifiers) - overlap,
        "reference_variants_absent_from_input": len(reference_ids) - overlap,
        "direct_id_matches": len(direct_matches),
        "coordinate_and_allele_matches": overlap - len(direct_matches),
        "ambiguous_input_coordinate_keys": len(duplicate_input_locations),
        "ambiguous_reference_coordinate_matches": len(
            ambiguous_coordinate_sources
        ) + target_conflicts,
        "compatible_allele_pairs": compatible_alleles,
        "incompatible_allele_pairs": incompatible_alleles,
        "reference_chromosomes": sorted(reference_chromosomes),
        "variant_id_match_policy": (
            "report_exact_bim_overlap_without_rewriting"
            if direct_input_mode else "allow_pipeline_reconciliation"
        ),
    }
    return metrics, set(mapping.values()) if retain_variant_ids else set(), mapping


def _harmonise_formatted_input(
    frame: pl.DataFrame,
    identifier_column: str,
    mapping: dict[str, str],
    source: Path,
    output: Path,
    dataset_id: str,
    configuration,
) -> tuple[Path, dict, Path | None]:
    """Write only reference-resolved rows and replace IDs with exact BIM IDs."""
    changed_identifiers = sum(source_id != target for source_id, target in mapping.items())
    unresolved = frame.height - len(mapping)
    if changed_identifiers == 0 and unresolved == 0:
        return source, {
            "harmonised_input_variants": frame.height,
            "variant_ids_replaced": 0,
            "unresolved_variants_removed": 0,
        }, None
    harmonised = frame.filter(
        pl.col(identifier_column).is_in(list(mapping))
    ).with_columns(
        pl.col(identifier_column).replace(mapping).alias(identifier_column)
    )
    duplicate_ids = harmonised.select(
        pl.col(identifier_column).is_duplicated().sum()
    ).item()
    if duplicate_ids:
        raise GctaGeneError(
            "Variant-ID reconciliation produced %d duplicate PLINK BIM IDs; "
            "coordinate fallback must be one-to-one." % int(duplicate_ids)
        )
    destination = configured_output_path(
        output,
        configuration.modules.gcta_gene.output_layout.harmonised_input,
        error_type=GctaGeneError,
        dataset_id=dataset_id,
        method=configuration.modules.gcta_gene.method,
    )
    resumed = False
    if destination.exists() and not configuration.run.overwrite:
        if not configuration.run.resume:
            raise GctaGeneError(
                "Harmonised GCTA input already exists: %s. Use --resume to "
                "validate and reuse it or --overwrite to replace it." % destination
            )
        destination.parent.mkdir(parents=True, exist_ok=True)
        handle = tempfile.NamedTemporaryFile(
            mode="w", dir=destination.parent, delete=False,
        )
        candidate = Path(handle.name)
        handle.close()
        try:
            write_dataframe_table(
                harmonised,
                candidate,
                overwrite=True,
                runtime=configuration.modules.formatting.runtime,
                error_type=GctaGeneError,
            )
            if sha256(candidate) != sha256(destination):
                raise GctaGeneError(
                    "Existing harmonised GCTA input does not match the current "
                    "summary statistics, VCF, reference, and policies: %s. Use "
                    "--overwrite after review." % destination
                )
            resumed = True
        finally:
            candidate.unlink(missing_ok=True)
    else:
        write_dataframe_table(
            harmonised,
            destination,
            overwrite=configuration.run.overwrite,
            runtime=configuration.modules.formatting.runtime,
            error_type=GctaGeneError,
        )
    return destination, {
        "harmonised_input_variants": harmonised.height,
        "variant_ids_replaced": changed_identifiers,
        "unresolved_variants_removed": unresolved,
        "harmonised_input_resumed": resumed,
    }, destination


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


def _prepared_resource_is_current(
    directory: Path,
    module,
    gmt: Path,
    gene_list: Path,
    bim: Path,
    analysis_source: Path,
    analyzable_variant_count: int,
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
        "empty_pathway_policy": conversion.empty_pathway_policy,
        "oversized_set_policy": module.set_annotation.oversized_set_policy,
        "maximum_set_variants": module.set_annotation.maximum_set_variants,
        "maximum_set_variants_applied": True,
        "variant_universe": "gwas_bim_intersection",
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
            manifest["schema_version"] == "gcta_fastbat_pathway_resource.v2"
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
            )
            logger.record(
                "SKIP", "gmt_to_fastbat_set", reason="validated_resume",
                output=str(set_path),
            )
            manifest = yaml.safe_load(
                (prepared_directory / names.manifest).read_text(encoding="utf-8")
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
        empty_pathway_policy=conversion.empty_pathway_policy,
        generation_command=shlex.join(sys.argv),
    )
    try:
        conversion_progress = StageProgress(
            "GMT-to-fastBAT preparation",
            enabled=configuration.logging.show_progress,
        )
        manifest = prepare_resource(
            arguments,
            stage_progress=conversion_progress,
            bim_variant_total=bim_variant_total,
            progress_refresh_seconds=(
                configuration.logging.progress_refresh_seconds
            ),
            mapping_workers=configuration.execution.threads,
            mapping_memory_gb=configuration.execution.memory_gb,
            worker_memory_multiplier=(
                conversion.parallelism.worker_memory_multiplier
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
        audit_level=manifest["policies"]["audit"]["level"],
        estimated_output_bytes=manifest["disk_preflight"][
            "estimated_output_bytes"
        ],
        mapping_workers=manifest["generation"]["mapping_workers"],
        mapping_algorithm=manifest["generation"]["mapping_algorithm"],
        output=str(set_path),
    )
    return set_path, manifest, prepared_directory


def _validate_gene_list(path: Path, module, reference_chromosomes: set[str]) -> dict:
    expected = len(module.gene_annotation.columns)
    column_index = {
        role: index for index, role in enumerate(module.gene_annotation.columns)
    }
    genes = set()
    chromosomes = set()
    with path.open("r", encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.split()
            if len(fields) != expected:
                raise GctaGeneError(
                    "Gene-list line %d has %d fields; expected %d."
                    % (number, len(fields), expected)
                )
            chromosome = fields[column_index["chromosome"]]
            start_text = fields[column_index["start"]]
            end_text = fields[column_index["end"]]
            gene = fields[column_index["gene"]]
            try:
                start, end = int(start_text), int(end_text)
            except ValueError as exc:
                raise GctaGeneError(
                    "Gene-list line %d has non-integer coordinates." % number
                ) from exc
            if start < 1 or end < start or not gene:
                raise GctaGeneError(
                    "Gene-list line %d has invalid coordinates or gene ID." % number
                )
            if gene in genes:
                raise GctaGeneError("Gene list contains duplicate gene ID: %s" % gene)
            genes.add(gene)
            chromosomes.add(chromosome.removeprefix("chr"))
    if not genes:
        raise GctaGeneError("Gene list contains no genes: %s" % path)
    overlap = chromosomes & reference_chromosomes
    if not overlap:
        raise GctaGeneError(
            "Gene list and LD reference do not share any chromosome labels."
        )
    return {
        "genes": len(genes),
        "gene_chromosomes": sorted(chromosomes),
        "shared_chromosomes": sorted(overlap),
    }


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
    seen_sets = set()
    input_sets = 0
    omitted_sets: list[str] = []
    oversized_sets: list[tuple[str, int]] = []
    requested_variants = 0
    matched_variants = 0
    analysis_variants = 0
    current_set = None
    current_variants: set[str] = set()
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
            with path.open("r", encoding="utf-8") as handle:
                for number, raw in enumerate(handle, 1):
                    value = raw.strip()
                    if not value:
                        continue
                    if current_set is None:
                        if value.upper() == "END":
                            raise GctaGeneError(
                                "fastBAT set list has END without a set ID at line "
                                "%d." % number
                            )
                        if len(value.split()) != 1:
                            raise GctaGeneError(
                                "fastBAT set ID at line %d must not contain "
                                "whitespace." % number
                            )
                        if value in seen_sets:
                            raise GctaGeneError(
                                "fastBAT set list repeats set ID: %s" % value
                            )
                        seen_sets.add(value)
                        current_set = value
                        current_variants = set()
                        current_matched = []
                        continue
                    if value.upper() == "END":
                        input_sets += 1
                        if not current_variants:
                            raise GctaGeneError(
                                "fastBAT set %s contains no variant IDs."
                                % current_set
                            )
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
                        current_set = None
                        continue
                    if len(value.split()) != 1:
                        raise GctaGeneError(
                            "fastBAT set-list line %d must contain one variant ID."
                            % number
                        )
                    if value in current_variants:
                        raise GctaGeneError(
                            "fastBAT set %s repeats variant ID %s."
                            % (current_set, value)
                        )
                    current_variants.add(value)
                    requested_variants += 1
                    if value in identifiers and value in reference_ids:
                        current_matched.append(value)
                        matched_variants += 1
        if current_set is not None:
            raise GctaGeneError(
                "fastBAT set %s is missing its terminating END line." % current_set
            )
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
        "input_sets": input_sets,
        "omitted_empty_sets": len(omitted_sets),
        "omitted_empty_set_examples": omitted_sets[:10],
        "omitted_oversized_sets": len(oversized_sets),
        "omitted_oversized_set_examples": [
            "%s (%d)" % item for item in oversized_sets[:10]
        ],
        "maximum_set_variants": maximum_set_variants,
        "requested_set_variants": requested_variants,
        "matched_set_variants": matched_variants,
        "analysis_set_variants": analysis_variants,
        "unmatched_set_variants": requested_variants - matched_variants,
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
    conversion = module["set_annotation"]["conversion"]
    conversion.pop("parallelism", None)
    conversion.pop("audit", None)
    conversion.pop("disk", None)
    module["results"].pop("normalized_schema_version", None)
    payload = {
        "module": module,
        "threads": configuration.execution.threads,
        "gcta": str(configuration.resources.executables.gcta),
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
    annotation_file: Path | None,
    reference_prefix: Path,
    module,
) -> dict[str, Path]:
    files = {"summary_statistics": source}
    if annotation_file is not None:
        files["annotation"] = annotation_file
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
    annotation_file: Path | None,
    reference_prefix: Path,
    final_result: Path,
    normalized: Path,
    result_metrics: dict,
) -> None:
    module = configuration.modules.gcta_gene
    inputs = _completion_files(
        source, annotation_file, reference_prefix, module,
    )
    write_yaml_report(
        {
            "schema_version": _COMPLETION_SCHEMA_VERSION,
            "status": "COMPLETED",
            "dataset_id": dataset_id,
            "method": module.method,
            "gcta_version": version,
            "configuration_sha256": _gcta_configuration_digest(configuration),
            "inputs": {
                name: _completion_fingerprint(file_path)
                for name, file_path in inputs.items()
            },
            "outputs": {
                "raw_result": _completion_fingerprint(final_result),
                "normalized_result": _completion_fingerprint(normalized),
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
    annotation_file: Path | None,
    reference_prefix: Path,
    final_result: Path,
    normalized: Path,
) -> dict:
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
    expected_scalars = {
        "schema_version": _COMPLETION_SCHEMA_VERSION,
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
    if mismatched:
        raise GctaGeneError(
            "GCTA completion manifest does not match the current run (%s). Use "
            "--overwrite after review." % ", ".join(mismatched)
        )
    current_inputs = {
        name: _completion_fingerprint(file_path)
        for name, file_path in _completion_files(
            source,
            annotation_file,
            reference_prefix,
            configuration.modules.gcta_gene,
        ).items()
    }
    current_outputs = {
        "raw_result": _completion_fingerprint(final_result),
        "normalized_result": _completion_fingerprint(normalized),
    }
    if manifest.get("inputs") != current_inputs:
        raise GctaGeneError(
            "GCTA completion inputs changed since the recorded run. Use "
            "--overwrite after review."
        )
    if manifest.get("outputs") != current_outputs:
        raise GctaGeneError(
            "GCTA completion outputs changed since the recorded run. Use "
            "--overwrite after review."
        )
    metrics = manifest.get("result_metrics")
    if not isinstance(metrics, dict):
        raise GctaGeneError("GCTA completion manifest has no result metrics.")
    return metrics


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
    harmonised_vcf: str | Path | None = None,
    reconcile_variant_ids: bool = False,
    dry_run: bool = False,
) -> ModuleResult:
    if harmonised_vcf is not None and not reconcile_variant_ids:
        raise GctaGeneError(
            "Direct GCTA input cannot use a harmonised VCF for ID rewriting."
        )
    module = configuration.modules.gcta_gene
    stage_total = (
        5
        + int(reconcile_variant_ids)
        + int(harmonised_vcf is not None)
    )
    if module.method in {"fastbat_gene", "mbat_combo"}:
        stage_total += 1
    elif module.method == "fastbat_set":
        stage_total += 2
    if dry_run:
        stage_total -= 1
    stage_number = 0

    def stage(title: str, function_name: str):
        nonlocal stage_number
        stage_number += 1
        return logger.step(
            stage_number, stage_total, title, function_name,
        )

    with stage(
        "Validate method and formatted GWAS input",
        "validate_gcta_input",
    ):
        _validate_build_and_population(module, configuration)
        source = _required_path(str(input_file), "GCTA summary-statistics input")
        output = Path(output_directory).expanduser().resolve()
        output.mkdir(parents=True, exist_ok=True)
        layout = module.output_layout
        annotation_file = None
        prepared_manifest = None
        prepared_directory = None
        analysis_set_path = None
        harmonised_input_path = None
        if module.method in {"fastbat_gene", "mbat_combo"}:
            if module.set_annotation.gmt_file is not None:
                raise GctaGeneError(
                    "--gmt is valid only with --method fastbat_set; the analysis "
                    "method is never inferred from an input filename."
                )
            annotation_file = _required_path(
                module.gene_annotation.file, "GCTA gene list",
            )
        elif module.method == "fastbat_set":
            if (
                module.set_annotation.file is None
                and module.set_annotation.gmt_file is None
            ):
                raise GctaGeneError(
                    "fastbat_set requires exactly one of --gmt or "
                    "--fastbat-set-list."
                )
        elif module.set_annotation.gmt_file is not None:
            raise GctaGeneError("--gmt is valid only with --method fastbat_set.")
        if module.reference.prefix is None:
            raise GctaGeneError("--gcta-reference-prefix is required.")
        reference_prefix = Path(module.reference.prefix).expanduser().resolve()
        (
            input_metrics,
            input_frame,
            identifier_column,
            identifiers,
            alleles,
        ) = _validate_formatted_input(
            source, configuration, module.method,
        )
    variant_locations = None
    location_metrics = {}
    if harmonised_vcf is not None:
        with stage(
            "Resolve GWAS variant coordinates",
            "variant_locations_from_vcf",
        ):
            variant_locations, location_metrics = _variant_locations_from_vcf(
                harmonised_vcf,
                identifiers,
                alleles,
                output,
                dataset_id,
                configuration,
                logger,
            )
    reference_stage_title = (
        "Validate and match the PLINK LD reference"
        if reconcile_variant_ids
        else "Compare GWAS and PLINK BIM variant IDs"
    )
    with stage(
        reference_stage_title,
        "validate_gcta_reference",
    ) as reference_stage:
        reference_metrics, reference_ids, identifier_mapping = _validate_reference(
            reference_prefix,
            module,
            identifiers,
            alleles,
            variant_locations,
            retain_variant_ids=module.method == "fastbat_set",
            direct_input_mode=not reconcile_variant_ids,
        )
        if not reconcile_variant_ids:
            absent_from_reference = reference_metrics[
                "input_variants_absent_from_reference"
            ]
            reference_stage.outcome(
                "Compared exact IDs without rewriting the direct GCTA input.",
                fields=(
                    ("count", "Summary-statistic unique IDs", len(identifiers)),
                    (
                        "count",
                        "PLINK BIM unique IDs",
                        reference_metrics["reference_variants"],
                    ),
                    (
                        "success",
                        "Exact IDs shared",
                        reference_metrics["overlapping_variants"],
                    ),
                    (
                        "warning" if absent_from_reference else "success",
                        "Summary IDs absent from BIM",
                        (
                            "%s; will not be used by GCTA"
                            % format(absent_from_reference, ",")
                            if absent_from_reference else "0"
                        ),
                    ),
                    (
                        "count",
                        "BIM IDs absent from summary",
                        reference_metrics["reference_variants_absent_from_input"],
                    ),
                ),
                input_unique_variant_ids=len(identifiers),
                reference_unique_variant_ids=reference_metrics[
                    "reference_variants"
                ],
                exact_ids_shared=reference_metrics["overlapping_variants"],
                input_ids_absent_from_reference=absent_from_reference,
                reference_ids_absent_from_input=reference_metrics[
                    "reference_variants_absent_from_input"
                ],
                input_rewritten=False,
            )
    if reconcile_variant_ids:
        with stage(
            "Harmonise GWAS variants to reference IDs",
            "harmonise_gcta_input",
        ):
            (
                source,
                harmonisation_metrics,
                harmonised_input_path,
            ) = _harmonise_formatted_input(
                input_frame,
                identifier_column,
                identifier_mapping,
                source,
                output,
                dataset_id,
                configuration,
            )
        identifiers = set(identifier_mapping.values())
    else:
        harmonisation_metrics = {
            "harmonised_input_variants": input_frame.height,
            "variant_ids_replaced": 0,
            "unresolved_variants_removed": 0,
            "harmonised_input_resumed": False,
            "direct_input_unmodified": True,
        }
    input_metrics.update(location_metrics)
    input_metrics.update(harmonisation_metrics)
    annotation_metrics = {}
    if module.method in {"fastbat_gene", "mbat_combo"}:
        with stage(
            "Validate gene-coordinate annotation",
            "validate_gcta_gene_list",
        ):
            annotation_metrics = _validate_gene_list(
                annotation_file, module,
                set(reference_metrics["reference_chromosomes"]),
            )
    elif module.method == "fastbat_set":
        with stage(
            "Prepare the fastBAT set source",
            "prepare_fastbat_set_source",
        ):
            if module.set_annotation.gmt_file is not None:
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
                )
            else:
                annotation_file = _required_path(
                    module.set_annotation.file, "GCTA fastBAT set list",
                )
        with stage(
            "Restrict fastBAT sets to analyzable variants",
            "prepare_fastbat_analysis_sets",
        ):
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
    del reference_ids
    logger.record("OBSERVED", "gcta_input", **input_metrics)
    logger.record("OBSERVED", "gcta_reference", **reference_metrics)
    if harmonised_input_path is not None:
        logger.record(
            "TRANSFORM",
            "gcta_variant_id_harmonisation",
            source=str(input_file),
            output=str(harmonised_input_path),
            direct_id_matches=reference_metrics["direct_id_matches"],
            coordinate_and_allele_matches=(
                reference_metrics["coordinate_and_allele_matches"]
            ),
            unresolved_variants_removed=(
                harmonisation_metrics["unresolved_variants_removed"]
            ),
            resumed=harmonisation_metrics["harmonised_input_resumed"],
            minimum_overlap_fraction=(
                module.variant_harmonisation.minimum_overlap_fraction
            ),
        )
    if annotation_metrics:
        logger.record("OBSERVED", "gcta_annotation", **annotation_metrics)
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

    with stage(
        "Validate the GCTA executable and version",
        "validate_gcta_software",
    ):
        executable = resolve_executable(
            configuration.resources.executables.gcta,
            "GCTA executable",
            error_type=GctaGeneError,
        )
        version = require_supported_gcta(
            executable, module, logger, configuration.execution.timeout_seconds,
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
    final_prefix.parent.mkdir(parents=True, exist_ok=True)
    resumed = False
    if (
        final_result.exists()
        and configuration.run.resume
        and not configuration.run.overwrite
    ):
        with stage(
            "Validate and reuse completed GCTA results",
            "resume_gcta_results",
        ):
            result_metrics = _validate_gcta_completion_manifest(
                completion_manifest,
                dataset_id=dataset_id,
                version=version,
                configuration=configuration,
                source=source,
                annotation_file=annotation_file,
                reference_prefix=reference_prefix,
                final_result=final_result,
                normalized=normalized,
            )
            logger.record(
                "SKIP", "gcta_execution", reason="validated_resume",
                output=str(final_result), manifest=str(completion_manifest),
            )
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
                    _write_gcta_completion_manifest(
                        completion_manifest,
                        dataset_id=dataset_id,
                        version=version,
                        configuration=configuration,
                        source=source,
                        annotation_file=annotation_file,
                        reference_prefix=reference_prefix,
                        final_result=final_result,
                        normalized=normalized,
                        result_metrics=result_metrics,
                    )
                except BaseException:
                    backup.replace(normalized)
                    raise
                finally:
                    backup.unlink(missing_ok=True)
                logger.record(
                    "TRANSFORM",
                    "gcta_multiple_testing_corrections",
                    source=str(final_result),
                    output=str(normalized),
                    **result_metrics["multiple_testing"],
                )
            resumed = True
    else:
        existing = [path for path in (final_result, normalized) if path.exists()]
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
        command = build_gcta_command(
            executable,
            source,
            reference_prefix,
            annotation_file,
            staged_prefix,
            configuration.execution.threads,
            module,
        )
        with stage(
            "Run GCTA %s and validate results" % module.method,
            "run_gcta_command",
        ) as step:
            step.input(
                "summary_statistics",
                path=str(source),
                rows=input_metrics["harmonised_input_variants"],
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
                if harmonised_input_path is not None:
                    dry_artifacts["harmonised_input"] = Artifact(
                        "gcta_reference_harmonised_input",
                        harmonised_input_path,
                        {
                            "dataset_id": dataset_id,
                            "method": module.method,
                            "genome_build": module.genome_build,
                            "source": "exact_bim_variant_ids",
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
            stage_normalized = staging / "normalized" / normalized.name
            try:
                result_metrics = normalize_gcta_results(
                    staged_result, stage_normalized, module,
                )
            except BaseException:
                if measured_progress is not None:
                    measured_progress.fail()
                raise
            if measured_progress is not None:
                measured_progress.complete(result_metrics["tested_units"])
            published = _publish_staged_outputs(
                staged_prefix, final_prefix, configuration.run.overwrite,
            )
            normalized.parent.mkdir(parents=True, exist_ok=True)
            if normalized.exists() and not configuration.run.overwrite:
                raise GctaGeneError("Output already exists: %s" % normalized)
            stage_normalized.replace(normalized)
            result_metrics["normalized_result"] = str(normalized)
            step.set_rows(result_metrics["tested_units"], removed=0)
            step.output("gcta_results", files=[str(path) for path in published])
            step.output("normalized_results", path=str(normalized))
            shutil.rmtree(staging)
            _write_gcta_completion_manifest(
                completion_manifest,
                dataset_id=dataset_id,
                version=version,
                configuration=configuration,
                source=source,
                annotation_file=annotation_file,
                reference_prefix=reference_prefix,
                final_result=final_result,
                normalized=normalized,
                result_metrics=result_metrics,
            )

    metadata = {
        "dataset_id": dataset_id,
        "method": module.method,
        "genome_build": module.genome_build,
        "reference_population": module.reference.population,
        "effect_allele": "A1=ALT",
        "schema_version": module.results.normalized_schema_version,
        "gcta_version": version,
    }
    artifacts = {
        "raw_results": Artifact("gcta_gene_raw", final_result, metadata),
        "normalized_results": Artifact("gcta_gene_results", normalized, metadata),
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
    if harmonised_input_path is not None:
        artifacts["harmonised_input"] = Artifact(
            "gcta_reference_harmonised_input",
            harmonised_input_path,
            {
                **metadata,
                "source": "exact_bim_variant_ids",
                "variant_ids_replaced": harmonisation_metrics[
                    "variant_ids_replaced"
                ],
                "unresolved_variants_removed": harmonisation_metrics[
                    "unresolved_variants_removed"
                ],
            },
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
    if stage_number != stage_total - 1:
        raise GctaGeneError(
            "Internal GCTA progress plan mismatch before summary: completed "
            "%d of %d planned stages." % (stage_number, stage_total)
        )
    with stage(
        "Build the scientific result summary",
        "build_gcta_scientific_summary",
    ):
        metrics["scientific_summary"] = build_gcta_scientific_summary(
            normalized,
            module,
            metrics,
        )
    return ModuleResult(
        "gcta_gene",
        artifacts=artifacts,
        metrics=metrics,
    )


def run_gcta_gene_direct(args, ctx=None) -> ModuleResult:
    """Resolve configuration once, execute the selected test, and finalize logging."""
    try:
        configuration = _resolved_configuration(args)
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
        stage_progress=StageProgress(
            "GCTA %s stages" % module.method.replace("_", " "),
            enabled=configuration.logging.show_progress,
        ),
    )
    try:
        input_file = getattr(args, "gcta_input_file", None) or module.input_file
        input_key = "summary_statistics_input_file"
        if input_file is None and ctx is not None:
            input_file = ctx.get("formatter", {}).get("gcta_gene", {}).get(
                input_key
            )
        pipeline_mode = ctx is not None
        supplied_vcf = getattr(args, "vcf", None)
        if not pipeline_mode and supplied_vcf is not None:
            raise GctaGeneError(
                "Direct gcta_gene does not accept --vcf because direct input "
                "mode compares exact SNP IDs without coordinate-based "
                "rewriting. Remove --vcf; summary-statistic IDs absent from "
                "the BIM will be reported and ignored by GCTA. Use postgwas "
                "pipeline --modules gcta_gene --vcf PATH when coordinate-and-"
                "allele reconciliation is required."
            )
        harmonised_vcf = supplied_vcf if pipeline_mode else None
        _require_gcta_gene_arguments(configuration, input_file)
        resource_paths = ["executables.gcta"]
        if harmonised_vcf is not None:
            resource_paths.append("executables.bcftools")
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
            frequency_difference_max=(
                module.frequency_difference_max if module.method == "mbat_combo" else None
            ),
            print_component_p_values=(
                module.print_component_p_values
                if module.method == "mbat_combo" else None
            ),
            write_snpset=module.write_snpset,
            reporting=module.reporting.model_dump(mode="json"),
            variant_harmonisation=module.variant_harmonisation.model_dump(
                mode="json"
            ),
            harmonised_vcf=(
                str(harmonised_vcf) if harmonised_vcf is not None else None
            ),
            variant_id_match_policy=(
                "allow_pipeline_reconciliation"
                if pipeline_mode else "report_exact_bim_overlap_without_rewriting"
            ),
            threads=configuration.execution.threads,
        )
        result = run_gcta_gene(
            input_file,
            output,
            dataset,
            configuration,
            logger,
            harmonised_vcf=harmonised_vcf,
            reconcile_variant_ids=pipeline_mode,
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


__all__ = ["run_gcta_gene", "run_gcta_gene_direct"]
