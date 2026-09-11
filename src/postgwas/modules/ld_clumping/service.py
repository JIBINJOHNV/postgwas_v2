"""Validated service boundary for LD-clumping analyses."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal
from uuid import uuid4

from pydantic import Field, field_validator

from postgwas.config import (
    load_run_configuration_for_module,
    write_resolved_configuration,
)
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.config.models.common import GenomeBuild, Population, StrictModel
from postgwas.core.paths import (
    configured_output_path,
    require_nonempty_file,
    resolve_executable,
)
from postgwas.core.contracts import RunContext
from postgwas.core.pipeline_logging import PipelineLogger
from postgwas.core.plink import validate_plink_bundle_dimensions
from postgwas.core.reference_resources import read_yaml_manifest
from postgwas.core.preflight import (
    PipelinePreflightEvidence,
    pipeline_preflight_evidence,
    require_pipeline_input_vcf,
)
from postgwas.core.required_arguments import (
    RequiredArgument,
    require_resolved_arguments,
)
from postgwas.core.vcf import (
    IndexedVcfValidation,
    required_vcf_query_tags,
    validate_indexed_vcf,
)
from postgwas.core.ui import StageProgress
from postgwas.modules.ld_clumping.common import LDClumpingError
from postgwas.modules.ld_clumping.ld_prune_cojo import (
    group_cojo_selected_signals,
    prepare_cojo_exclusion_list,
)
from postgwas.modules.ld_clumping.ld_prune_region import ld_clump_by_regions
from postgwas.modules.ld_clumping.ld_prune_standard import ld_clump_standard
from postgwas.modules.ld_clumping.reporting import (
    build_ld_clumping_summary,
    format_cojo_output_scope,
    format_cojo_reference_match,
    render_ld_clumping_summary,
    write_ld_clumping_html_report,
    write_ld_clumping_summary_csv,
)
from postgwas.modules.formatting.reference_identifiers import (
    BimIdentifierRequirement,
    configure_reference_variant_identifiers,
)
from postgwas.modules.gcta_cojo.parallel import (
    chromosome_parallel_slct_requested,
)


class LDReferenceManifest(StrictModel):
    """Scientific and indexing contract written by LD reference preparation."""

    format_version: int = Field(ge=1)
    genome_build: GenomeBuild
    populations: list[Population]
    orientation: Literal[
        "symmetric_first_endpoint", "upper_triangle_dual_index"
    ]
    file_pattern: str
    reverse_file_pattern: str
    variant_inventory_pattern: str
    columns: list[str]
    variant_inventory_columns: list[str]
    window_kb: int = Field(gt=0)
    minimum_r2: float = Field(ge=0, le=1, allow_inf_nan=False)
    minimum_maf: float = Field(ge=0, le=0.5, allow_inf_nan=False)
    allele_order_preserved: bool
    plink_version: str

    @field_validator("populations")
    @classmethod
    def unique_populations(cls, values: list[Population]) -> list[Population]:
        if not values or len(values) != len(set(values)):
            raise ValueError("populations must contain one or more unique values")
        return values

    @field_validator("plink_version")
    @classmethod
    def plink_19_only(cls, value: str) -> str:
        value = value.strip()
        if not value.startswith("PLINK v1.9"):
            raise ValueError("must identify the PLINK 1.9 executable used")
        return value


@dataclass(frozen=True)
class LDClumpingInputValidation:
    """GWAS-VCF evidence established before reference validation."""

    configuration: Any
    vcf: Path
    output_directory: Path
    dataset_id: str
    bcftools: str
    genome_build: str
    sample: str
    variant_count: int
    contigs: tuple[str, ...]
    required_fields: tuple[str, ...]


@dataclass(frozen=True)
class LDClumpingPreflight:
    configuration: Any
    vcf: Path
    output_directory: Path
    dataset_id: str
    bcftools: str
    tabix: str | None
    reference_directory: Path | None
    reference_manifest: LDReferenceManifest | None
    gcta: str | None = None
    cojo_reference_prefix: Path | None = None
    cojo_reference_files: dict[str, Path] | None = None
    cojo_reference_validation: Any | None = None
    input_genome_build: str | None = None
    input_sample: str | None = None
    input_variant_count: int | None = None
    input_contigs: tuple[str, ...] = ()
    input_required_fields: tuple[str, ...] = ()
    cojo_reference_identifier_observation: dict[str, Any] | None = None


def _owned_analysis_paths(output_directory: Path, module, dataset_id: str) -> list[Path]:
    """Resolve only module-owned scientific and method-log paths."""
    values = {"dataset_id": dataset_id, "population": module.population.value}
    paths = []
    for name, pattern in module.output_layout.model_dump().items():
        if name in module.output_layout.shared_cojo_directory_fields:
            # These are shared containers. Subordinate formatter/GCTA services
            # own and validate their dataset-scoped artifacts within them.
            continue
        if not name.startswith(("region_", "standard_", "cojo_")) and name not in {
            "summary_csv",
            "html_report",
        }:
            continue
        paths.append(
            configured_output_path(
                output_directory,
                pattern,
                error_type=LDClumpingError,
                **values,
            )
        )
    return paths


def _report_paths(output_directory: Path, module, dataset_id: str) -> dict[str, Path]:
    values = {"dataset_id": dataset_id, "population": module.population.value}
    return {
        name: configured_output_path(
            output_directory,
            getattr(module.output_layout, name),
            error_type=LDClumpingError,
            **values,
        )
        for name in ("summary_csv", "html_report")
    }


def _reported_output_paths(outputs: dict[str, Any]) -> dict[str, str | None]:
    """Flatten declared method artifacts for the CSV and HTML output index."""
    paths = {
        "canonical_log": outputs.get("canonical_log"),
        "resolved_configuration": outputs.get("resolved_configuration"),
        "summary_csv": outputs.get("summary_csv"),
        "html_report": outputs.get("html_report"),
    }
    for method, result_key in (
        ("region", "ld_clump_region"),
        ("standard", "ld_clump_standard"),
        ("cojo", "ld_clump_cojo"),
    ):
        result = outputs.get(result_key) or {}
        for name, path in (result.get("output_files") or {}).items():
            paths["%s_%s" % (method, name)] = path
    return paths


def _path_identity(path: Path):
    try:
        stat = path.stat()
    except FileNotFoundError:
        return None
    return (stat.st_size, stat.st_mtime_ns)


def _quarantine_changed_outputs(
    paths: list[Path],
    before: dict[Path, tuple[int, int] | None],
    logger: PipelineLogger,
    *,
    validated_failure_outputs: tuple[str, ...] = (),
) -> None:
    """Quarantine incomplete outputs while retaining validated failure audits."""
    preserved = {
        Path(value).expanduser().resolve() for value in validated_failure_outputs
    }
    for path in paths:
        if _path_identity(path) == before[path] or not path.exists():
            continue
        if path.resolve() in preserved:
            logger.record(
                "OUTPUT",
                "validated_failure_audit_preserved",
                path=str(path),
            )
            continue
        partial = path.with_name("%s.partial.%s" % (path.name, uuid4().hex))
        try:
            path.replace(partial)
        except OSError as exc:
            logger.error(
                "Cannot quarantine incomplete output %s: %s" % (path, exc)
            )
        else:
            logger.record(
                "ACTION",
                "quarantine_incomplete_output",
                original=str(path),
                partial=str(partial),
            )


MODULE_CLI_OVERRIDES = {
    "vcf": "inputs.vcf",
    "dataset_id": "inputs.dataset_id",
    "output_directory": "output_directory",
    "clumping_methods": "methods",
    "genome_build": "genome_build",
    "population": "population",
    "lead_p": "lead_pvalue",
    "candidate_p": "candidate_pvalue",
    "r2_clump": "clump_r2",
    "r2_lead": "lead_r2",
    "ld_window_kb": "window_kb",
    "merge_dist": "merge_distance_bp",
    "cojo_merge_dist": "cojo.merge_distance_bp",
    "missing_index_action": "missing_index_action",
    "missing_chromosome_action": "missing_chromosome_action",
    "remove_mhc": "remove_mhc",
    "ld_folder": "reference.directory",
    "cojo_reference_prefix": "modules.gcta_cojo.reference.prefix",
    "cojo_p": "modules.gcta_cojo.analysis.significance_threshold",
    "cojo_window_kb": "modules.gcta_cojo.analysis.window_kb",
    "cojo_collinear": "modules.gcta_cojo.analysis.collinearity_cutoff",
    "cojo_diff_freq": "modules.gcta_cojo.analysis.frequency_difference_max",
    "cojo_maf": "modules.gcta_cojo.analysis.reference_maf_min",
    "cojo_minimum_reference_overlap": (
        "modules.gcta_cojo.input_validation.minimum_reference_overlap_fraction"
    ),
    "cojo_gc": "modules.gcta_cojo.analysis.genomic_control",
    "cojo_gc_lambda": "modules.gcta_cojo.analysis.genomic_control_lambda",
    "cojo_chromosome": "modules.gcta_cojo.analysis.chromosome",
}

GLOBAL_CLI_OVERRIDES = {
    "dataset_id": "run.dataset_id",
    "output_directory": "run.output_directory",
    "threads": "execution.threads",
    "memory_gb": "execution.memory_gb",
    "seed": "execution.random_seed",
    "resume": "run.resume",
    "overwrite": "run.overwrite",
    "bcftools": "resources.executables.bcftools",
    "tabix": "resources.executables.tabix",
    "gcta": "resources.executables.gcta",
}


def resolve_ld_clumping_configuration(args: argparse.Namespace):
    """Resolve explicit CLI overrides over module or full-run YAML once."""
    module_overrides = {
        path: value
        for path, value in explicit_overrides(args, MODULE_CLI_OVERRIDES).items()
        if value is not None
    }
    global_overrides = {
        path: value
        for path, value in explicit_overrides(args, GLOBAL_CLI_OVERRIDES).items()
        if value is not None
    }
    if "modules.gcta_cojo.analysis.genomic_control_lambda" in module_overrides:
        module_overrides["modules.gcta_cojo.analysis.genomic_control"] = True
    configuration = load_run_configuration_for_module(
        "ld_clumping",
        getattr(args, "run_config", None),
        module_overrides=module_overrides,
        global_overrides=global_overrides,
    )
    module = configuration.modules.ld_clumping
    if "cojo-slct" in module.methods:
        cojo = configuration.modules.gcta_cojo
        cojo_reference = cojo.reference.model_copy(update={
            "population": module.population.value,
        })
        configuration.modules.gcta_cojo = cojo.model_copy(update={
            "genome_build": module.genome_build.value,
            "reference": cojo_reference,
        })
    return configuration


def _read_reference_manifest(
    directory: Path,
    filename: str,
) -> LDReferenceManifest:
    return read_yaml_manifest(
        directory / filename, LDReferenceManifest, "LD reference manifest",
        error_type=LDClumpingError,
    )


def _validate_reference_contract(module, directory: Path) -> LDReferenceManifest:
    reference = module.reference
    manifest = _read_reference_manifest(directory, reference.manifest_filename)
    mismatches = []
    if manifest.format_version != reference.format_version:
        mismatches.append(
            "format_version=%s (expected %s)"
            % (manifest.format_version, reference.format_version)
        )
    if manifest.genome_build != module.genome_build:
        mismatches.append(
            "genome_build=%s (expected %s)"
            % (manifest.genome_build.value, module.genome_build.value)
        )
    if module.population not in manifest.populations:
        mismatches.append(
            "population=%s is absent from populations=%s"
            % (
                module.population.value,
                ",".join(value.value for value in manifest.populations),
            )
        )
    if manifest.orientation != reference.orientation:
        mismatches.append(
            "orientation=%s (expected %s)"
            % (manifest.orientation, reference.orientation)
        )
    if manifest.file_pattern != reference.file_pattern:
        mismatches.append(
            "file_pattern=%s (expected %s)"
            % (manifest.file_pattern, reference.file_pattern)
        )
    if manifest.reverse_file_pattern != reference.reverse_file_pattern:
        mismatches.append(
            "reverse_file_pattern=%s (expected %s)"
            % (manifest.reverse_file_pattern, reference.reverse_file_pattern)
        )
    if manifest.variant_inventory_pattern != reference.variant_inventory_pattern:
        mismatches.append(
            "variant_inventory_pattern=%s (expected %s)"
            % (
                manifest.variant_inventory_pattern,
                reference.variant_inventory_pattern,
            )
        )
    if manifest.columns != reference.columns:
        mismatches.append("columns do not match the configured seven-column contract")
    if manifest.variant_inventory_columns != reference.variant_inventory_columns:
        mismatches.append(
            "variant_inventory_columns do not match the configured contract"
        )
    if manifest.window_kb < module.window_kb:
        mismatches.append(
            "window_kb=%s is smaller than requested window_kb=%s"
            % (manifest.window_kb, module.window_kb)
        )
    required_minimum_r2 = min(module.clump_r2, module.lead_r2)
    if manifest.minimum_r2 > required_minimum_r2:
        mismatches.append(
            "minimum_r2=%s omits pairs needed at requested r2=%s"
            % (manifest.minimum_r2, required_minimum_r2)
        )
    if manifest.minimum_maf > module.minimum_reference_maf:
        mismatches.append(
            "minimum_maf=%s omits variants needed at requested "
            "minimum_reference_maf=%s"
            % (manifest.minimum_maf, module.minimum_reference_maf)
        )
    if not manifest.allele_order_preserved:
        mismatches.append("allele_order_preserved=false")
    if mismatches:
        raise LDClumpingError(
            "LD reference is scientifically incompatible: %s"
            % "; ".join(mismatches)
        )
    return manifest


def _validate_cojo_reference(
    args: argparse.Namespace,
    configuration,
    *,
    logger: PipelineLogger | None = None,
):
    """Validate the shared GCTA selection reference before VCF formatting."""
    from postgwas.modules.gcta_cojo.errors import GctaCojoError
    from postgwas.modules.gcta_cojo.service import validate_gcta_cojo_reference

    module = configuration.modules.gcta_cojo
    if module.mode != "slct":
        raise LDClumpingError(
            "--clumping-methods cojo-slct requires "
            "modules.gcta_cojo.mode: slct; the resolved mode is %s."
            % module.mode
        )
    try:
        validated = validate_gcta_cojo_reference(
            configuration,
            logger=logger,
        )
    except (GctaCojoError, OSError, ValueError) as exc:
        raise LDClumpingError(str(exc)) from exc
    reference_prefix = validated.prefix
    reference_files = validated.files
    bim_suffixes = [
        suffix for suffix in module.reference.required_extensions
        if suffix.lower() == ".bim"
    ]
    if len(bim_suffixes) != 1:
        raise LDClumpingError(
            "modules.gcta_cojo.reference.required_extensions must contain "
            "exactly one .bim suffix."
        )
    configure_reference_variant_identifiers(
        args,
        configuration.modules.formatting,
        [BimIdentifierRequirement(
            consumer="LD clumping COJO selection",
            formatter_target="gcta_gene",
            bim_file=reference_files["bim"],
            column_roles=module.reference.bim_columns,
            delimiter_pattern=module.reference.table_delimiter_pattern,
        )],
    )
    return validated


def _required_query_expressions(module, *, include_region_annotation: bool) -> list[str]:
    fields = module.vcf_fields
    selected = [
        fields.chromosome,
        fields.position,
        fields.reference_allele,
        fields.alternate_allele,
        fields.variant_id,
        fields.effect,
        fields.standard_error,
        fields.allele_frequency,
        fields.log_pvalue,
    ]
    if "region" in module.methods and include_region_annotation:
        selected.append(fields.ld_block.format(population=module.population.value))
    return selected


def _resolved_run_values(configuration) -> tuple[Path, str]:
    module = configuration.modules.ld_clumping
    output_value = module.output_directory or configuration.run.output_directory
    dataset_id = module.inputs.dataset_id or configuration.run.dataset_id
    if output_value is None:
        raise LDClumpingError(
            "Required argument not provided: --output-directory. Provide "
            "--output-directory VALUE or set modules.ld_clumping.output_directory "
            "in the run configuration."
        )
    if not dataset_id:
        raise LDClumpingError(
            "Required argument not provided: --dataset-id. Provide --dataset-id "
            "VALUE or set modules.ld_clumping.inputs.dataset_id in the run "
            "configuration."
        )
    return Path(output_value).expanduser().resolve(), str(dataset_id)


def _require_ld_clumping_arguments(
    args: argparse.Namespace,
    configuration,
    *,
    pipeline: bool,
) -> None:
    """Validate all selected-mode requirements before invoking any tool."""
    module = configuration.modules.ld_clumping
    requirements = [
        RequiredArgument(
            "--vcf", "modules.ld_clumping.inputs.vcf", module.inputs.vcf,
        ),
        RequiredArgument(
            "--dataset-id",
            "modules.ld_clumping.inputs.dataset_id",
            module.inputs.dataset_id or configuration.run.dataset_id,
        ),
        RequiredArgument(
            "--output-directory",
            "modules.ld_clumping.output_directory",
            module.output_directory or configuration.run.output_directory,
        ),
    ]
    if "standard" in module.methods:
        requirements.append(
            RequiredArgument(
                "--ld-folder",
                "modules.ld_clumping.reference.directory",
                module.reference.directory,
            )
        )
    if "cojo-slct" in module.methods:
        requirements.append(
            RequiredArgument(
                "--cojo-reference-prefix",
                "modules.gcta_cojo.reference.prefix",
                configuration.modules.gcta_cojo.reference.prefix,
            )
        )
    require_resolved_arguments(requirements)
    if (
        pipeline
        and "finemap" in getattr(args, "modules", ())
        and "standard" not in module.methods
    ):
        raise LDClumpingError(
            "Fine-mapping requires standard LD clumping; add standard to "
            "modules.ld_clumping.methods or --clumping-methods."
        )


def _validate_ld_clumping_input(
    args: argparse.Namespace,
    *,
    configuration,
    pipeline: bool,
    pipeline_entry: bool,
    logger: PipelineLogger | None,
    cached_vcf: IndexedVcfValidation | None = None,
) -> LDClumpingInputValidation:
    """Validate and count the GWAS-VCF before inspecting LD references."""
    _require_ld_clumping_arguments(args, configuration, pipeline=pipeline)
    module = configuration.modules.ld_clumping
    vcf = require_nonempty_file(
        module.inputs.vcf,
        "LD-clumping input VCF",
        error_type=LDClumpingError,
    )
    output_directory, dataset_id = _resolved_run_values(configuration)
    try:
        bcftools = resolve_executable(
            configuration.resources.executables.bcftools,
            "bcftools executable",
            error_type=LDClumpingError,
        )
    except OSError as exc:
        raise LDClumpingError(str(exc)) from exc

    query_expressions = _required_query_expressions(
        module,
        include_region_annotation=not pipeline_entry,
    )
    if "cojo-slct" in module.methods:
        query_expressions.extend(
            configuration.modules.formatting.vcf_fields.root.values()
        )
    required_fields = tuple(required_vcf_query_tags(query_expressions))
    indexed = validate_indexed_vcf(
        vcf,
        dataset_id,
        bcftools,
        genome_build_header=(
            configuration.modules.harmonisation.vcf_processing.genome_build_header
        ),
        supported_genome_builds=list(configuration.resources.genomes),
        required_fields=required_fields,
        provenance_headers=(
            configuration.modules.formatting.input_contract.provenance_headers.model_dump()
        ),
        expected_genome_build=module.genome_build.value,
        cached=cached_vcf,
        logger=logger,
        error_type=LDClumpingError,
    )
    if logger is not None:
        logger.record(
            "PASS",
            "ld_clumping_input_vcf",
            input_vcf=str(vcf),
            variants=indexed.variant_count,
            genome_build=indexed.genome_build,
            sample=indexed.sample,
            contigs=indexed.contigs,
            required_fields=required_fields,
            validation_reused=indexed.reused,
        )
    return LDClumpingInputValidation(
        configuration=configuration,
        vcf=vcf,
        output_directory=output_directory,
        dataset_id=dataset_id,
        bcftools=bcftools,
        genome_build=indexed.genome_build,
        sample=indexed.sample,
        variant_count=indexed.variant_count,
        contigs=indexed.contigs,
        required_fields=required_fields,
    )


def _validate_ld_clumping_references(
    args: argparse.Namespace,
    input_validation: LDClumpingInputValidation,
    *,
    logger: PipelineLogger | None,
) -> LDClumpingPreflight:
    """Validate references selected by the already validated GWAS-VCF run."""
    configuration = input_validation.configuration
    module = configuration.modules.ld_clumping
    tabix = None
    if "standard" in module.methods:
        try:
            tabix = resolve_executable(
                configuration.resources.executables.tabix,
                "tabix executable",
                error_type=LDClumpingError,
            )
        except OSError as exc:
            raise LDClumpingError(str(exc)) from exc

    reference_directory = None
    manifest = None
    if "standard" in module.methods:
        reference_directory = Path(module.reference.directory).expanduser().resolve()
        if not reference_directory.is_dir():
            raise LDClumpingError(
                "LD reference directory does not exist: %s" % reference_directory
            )
        try:
            nonempty = next(reference_directory.iterdir(), None) is not None
        except OSError as exc:
            raise LDClumpingError(
                "Cannot read LD reference directory %s: %s"
                % (reference_directory, exc)
            ) from exc
        if not nonempty:
            raise LDClumpingError(
                "LD reference directory is empty: %s" % reference_directory
            )
        manifest = _validate_reference_contract(module, reference_directory)

    gcta = None
    cojo_reference_prefix = None
    cojo_reference_files = None
    cojo_reference_validation = None
    cojo_identifier_observation = None
    if "cojo-slct" in module.methods:
        cojo_reference_validation = _validate_cojo_reference(
            args,
            configuration,
            logger=logger,
        )
        gcta = cojo_reference_validation.executable
        cojo_reference_prefix = cojo_reference_validation.prefix
        cojo_reference_files = cojo_reference_validation.files
        cojo_identifier_observation = dict(
            (getattr(args, "variant_id_observations", None) or {}).get(
                "gcta_gene", {}
            )
        )

    if logger is not None:
        logger.record(
            "PASS",
            "ld_clumping_preflight",
            input_vcf=str(input_validation.vcf),
            input_variants=input_validation.variant_count,
            input_sample=input_validation.sample,
            genome_build=input_validation.genome_build,
            population=module.population.value,
            methods=module.methods,
            reference_manifest=(
                None
                if reference_directory is None
                else str(reference_directory / module.reference.manifest_filename)
            ),
            reference_contract=(
                None if manifest is None else manifest.model_dump(mode="json")
            ),
            cojo_reference_prefix=(
                None
                if cojo_reference_prefix is None
                else str(cojo_reference_prefix)
            ),
            cojo_reference_files=(
                None
                if cojo_reference_files is None
                else {
                    name: str(path)
                    for name, path in cojo_reference_files.items()
                }
            ),
            cojo_reference_samples=(
                None
                if cojo_reference_validation is None
                else cojo_reference_validation.samples
            ),
            cojo_reference_identifier_observation=cojo_identifier_observation,
            gcta_version=(
                None
                if cojo_reference_validation is None
                else cojo_reference_validation.version
            ),
        )
    return LDClumpingPreflight(
        configuration=configuration,
        vcf=input_validation.vcf,
        output_directory=input_validation.output_directory,
        dataset_id=input_validation.dataset_id,
        bcftools=input_validation.bcftools,
        tabix=tabix,
        reference_directory=reference_directory,
        reference_manifest=manifest,
        gcta=gcta,
        cojo_reference_prefix=cojo_reference_prefix,
        cojo_reference_files=cojo_reference_files,
        cojo_reference_validation=cojo_reference_validation,
        input_genome_build=input_validation.genome_build,
        input_sample=input_validation.sample,
        input_variant_count=input_validation.variant_count,
        input_contigs=input_validation.contigs,
        input_required_fields=input_validation.required_fields,
        cojo_reference_identifier_observation=cojo_identifier_observation,
    )


def validate_ld_clumping_configuration(
    args: argparse.Namespace,
    *,
    configuration=None,
    pipeline: bool = False,
    pipeline_entry: bool = False,
    logger: PipelineLogger | None = None,
    cached_vcf: IndexedVcfValidation | None = None,
) -> LDClumpingPreflight:
    """Validate every selected method before creating scientific outputs."""
    configuration = configuration or resolve_ld_clumping_configuration(args)
    input_validation = _validate_ld_clumping_input(
        args,
        configuration=configuration,
        pipeline=pipeline,
        pipeline_entry=pipeline_entry,
        logger=logger,
        cached_vcf=cached_vcf,
    )
    return _validate_ld_clumping_references(
        args,
        input_validation,
        logger=logger,
    )


def preflight_ld_clumping(
    args: argparse.Namespace,
    *,
    preflight_evidence=None,
) -> PipelinePreflightEvidence:
    """Pipeline preflight using the same resolved validation as direct runs."""
    shared = require_pipeline_input_vcf(preflight_evidence)
    resources = validate_ld_clumping_configuration(
        args,
        pipeline=True,
        pipeline_entry=True,
        cached_vcf=shared["indexed"],
    )
    if resources.cojo_reference_validation is not None:
        reference = resources.cojo_reference_validation
        validate_plink_bundle_dimensions(
            reference.files,
            variants=resources.cojo_reference_identifier_observation["variants"],
            samples=reference.samples,
            error_type=LDClumpingError,
        )
    return pipeline_preflight_evidence(
        "ld_clump",
        preflight_evidence,
        resources=resources,
        deferred_checks=(
            "Validate any pipeline-added LD-block annotation on the current VCF.",
        ),
    )


def _input_validation_outcome_fields(
    validation: LDClumpingInputValidation,
) -> list[tuple]:
    """Render the validated GWAS-VCF evidence without recomputation."""
    return [
        ("analysis", "Input summary-statistics VCF"),
        ("info", "Input file", validation.vcf.name),
        ("count", "Total variants", validation.variant_count),
        (
            "genetic", "Genome build inferred from VCF",
            validation.genome_build,
        ),
        ("info", "VCF embedded dataset/sample", validation.sample),
        ("success", "VCF structural validation", "passed"),
    ]


def _reference_validation_outcome_fields(
    preflight: LDClumpingPreflight,
    module,
) -> list[tuple]:
    """Render only reference facts established by the preflight checks."""
    fields: list[tuple] = []
    for method in module.methods:
        if method == "region":
            ld_block_field = module.vcf_fields.ld_block.format(
                population=module.population.value,
            ).lstrip("%")
            fields.extend((
                ("analysis", "Annotated-region reference"),
                ("genetic", "Required VCF annotation", ld_block_field),
                ("success", "Annotation declaration", "present in VCF header"),
            ))
            continue
        if method == "standard":
            manifest = preflight.reference_manifest
            reference_directory = preflight.reference_directory
            if manifest is None or reference_directory is None:
                raise LDClumpingError(
                    "Standard LD-reference validation evidence is missing after "
                    "a successful preflight."
                )
            resources_per_chromosome = (
                3
                if manifest.orientation == "upper_triangle_dual_index"
                else 2
            )
            fields.extend((
                ("analysis", "FUMA-style indexed LD reference"),
                ("info", "Reference directory", reference_directory.name),
                (
                    "info", "Manifest file",
                    module.reference.manifest_filename,
                ),
                (
                    "info", "Reference format",
                    "version %d · %s"
                    % (
                        manifest.format_version,
                        (
                            "upper-triangle dual index"
                            if manifest.orientation == "upper_triangle_dual_index"
                            else "symmetric first endpoint"
                        ),
                    ),
                ),
                ("genetic", "Declared genome build", manifest.genome_build.value),
                (
                    "genetic", "Declared populations",
                    ", ".join(value.value for value in manifest.populations),
                ),
                (
                    "analysis", "Stored LD contract",
                    "%s kb · r² ≥ %s · MAF ≥ %s"
                    % (
                        manifest.window_kb,
                        manifest.minimum_r2,
                        manifest.minimum_maf,
                    ),
                ),
                (
                    "success", "Allele-order preservation",
                    "confirmed by manifest",
                ),
                ("success", "Manifest compatibility", "passed"),
                (
                    "info", "Chromosome resource checks",
                    "%d compressed tables + %d tabix indexes per GWS "
                    "chromosome; checked after GWS chromosomes are identified"
                    % (resources_per_chromosome, resources_per_chromosome),
                ),
            ))
            continue
        if method == "cojo-slct":
            validation = preflight.cojo_reference_validation
            observation = (
                preflight.cojo_reference_identifier_observation or {}
            )
            if validation is None or preflight.cojo_reference_prefix is None:
                raise LDClumpingError(
                    "COJO reference validation evidence is missing after a "
                    "successful preflight."
                )
            identifier_type = observation.get("variant_id_type")
            identifier_label = (
                "rsID"
                if identifier_type == "rsid"
                else "chromosome-position-allele ID"
                if identifier_type == "unique"
                else "not available"
            )
            companion_files = ", ".join(
                path.suffix.removeprefix(".").upper()
                for path in validation.files.values()
            )
            fields.extend((
                ("analysis", "PLINK LD reference for GCTA-COJO"),
                (
                    "info", "Reference prefix",
                    preflight.cojo_reference_prefix.name,
                ),
                (
                    "success", "Required companion files",
                    "%s present and non-empty" % companion_files,
                ),
                ("count", "Total BIM variants", observation.get("variants")),
                ("count", "LD-reference samples", validation.samples),
                ("genetic", "Variant identifiers", identifier_label),
                (
                    "genetic", "Declared genome build",
                    module.genome_build.value,
                ),
                (
                    "genetic", "Declared population",
                    module.population.value,
                ),
                ("success", "BIM structural validation", "passed"),
                ("success", "GCTA runtime validation", validation.version),
                (
                    "warning", "Build and population provenance",
                    "taken from configuration; genome build and population "
                    "are not independently verifiable from PLINK file contents",
                ),
            ))
            continue
        raise LDClumpingError("Unsupported LD-clumping method: %s" % method)
    return fields


def run_ld_clump_direct(
    args: argparse.Namespace,
    ctx: dict[str, Any] | None = None,
    *,
    configuration=None,
    pipeline: bool = False,
    cached_vcf: IndexedVcfValidation | None = None,
):
    """Run every selected clumping method; any selected-method failure is fatal."""
    configuration = configuration or resolve_ld_clumping_configuration(args)
    module = configuration.modules.ld_clumping
    output_directory, dataset_id = _resolved_run_values(configuration)
    output_directory.mkdir(parents=True, exist_ok=True)
    log_path = configured_output_path(
        output_directory,
        module.output_layout.canonical_log,
        dataset_id=dataset_id,
        population=module.population.value,
        error_type=LDClumpingError,
    )
    logger = PipelineLogger(
        sample_id=dataset_id,
        scope="run",
        log_dir=str(log_path.parent),
        level=configuration.logging.file_level,
        screen_level=configuration.logging.console_level,
        log_path=str(log_path),
        stage_progress=StageProgress(
            "LD clumping analysis progress",
            enabled=configuration.logging.show_progress,
            outcome_label_width=configuration.logging.terminal_label_width,
        ),
    )
    owned_paths = _owned_analysis_paths(
        output_directory, module, dataset_id,
    )
    before_outputs = {path: _path_identity(path) for path in owned_paths}
    try:
        total_steps = len(module.methods) + 3
        step_number = 1
        with logger.step(
            step_number,
            total_steps,
            "Validate the input GWAS-VCF",
            "validate_ld_clumping_input",
        ) as step:
            input_validation = _validate_ld_clumping_input(
                args,
                configuration=configuration,
                pipeline=pipeline,
                pipeline_entry=False,
                logger=logger,
                cached_vcf=cached_vcf,
            )
            step.set_rows(input_validation.variant_count)
            step.outcome(
                "Validated the GWAS-VCF structure, sample, build, and indexed "
                "record count.",
                fields=_input_validation_outcome_fields(input_validation),
                input_vcf=str(input_validation.vcf),
                variants=input_validation.variant_count,
                genome_build=input_validation.genome_build,
                sample=input_validation.sample,
                contigs=input_validation.contigs,
                required_fields=input_validation.required_fields,
            )
        step_number += 1

        with logger.step(
            step_number,
            total_steps,
            "Validate selected LD-reference inputs",
            "validate_ld_clumping_references",
        ) as step:
            preflight = _validate_ld_clumping_references(
                args,
                input_validation,
                logger=logger,
            )
            resolved_path = configured_output_path(
                preflight.output_directory,
                module.output_layout.resolved_configuration,
                dataset_id=preflight.dataset_id,
                population=module.population.value,
                error_type=LDClumpingError,
            )
            resolved_path.parent.mkdir(parents=True, exist_ok=True)
            write_resolved_configuration(
                configuration,
                resolved_path,
                modules=(
                    ("ld_clumping", "formatting", "gcta_cojo")
                    if "cojo-slct" in module.methods
                    else ("ld_clumping",)
                ),
                resource_paths=(
                    ("executables.bcftools", "executables.tabix", "executables.gcta")
                    if "cojo-slct" in module.methods
                    else ("executables.bcftools", "executables.tabix")
                ),
            )
            logger.record(
                "PARAM",
                "resolved_ld_clumping_configuration",
                methods=module.methods,
                lead_pvalue=module.lead_pvalue,
                candidate_pvalue=module.candidate_pvalue,
                clump_r2=module.clump_r2,
                lead_r2=module.lead_r2,
                window_kb=module.window_kb,
                merge_distance_bp=module.merge_distance_bp,
                cojo=module.cojo.model_dump(mode="json"),
                gcta_cojo=(
                    configuration.modules.gcta_cojo.model_dump(mode="json")
                    if "cojo-slct" in module.methods
                    else None
                ),
                remove_mhc=module.remove_mhc,
                missing_index_action=module.missing_index_action,
                minimum_reference_maf=module.minimum_reference_maf,
                missing_chromosome_action=module.missing_chromosome_action,
                html_result_table_columns=(
                    module.reporting.result_table_columns
                ),
                threads=configuration.execution.threads,
                memory_gb=configuration.execution.memory_gb,
                resolved_configuration=str(resolved_path),
            )
            step.outcome(
                "Validated the reference contract required by every selected "
                "clumping method.",
                fields=_reference_validation_outcome_fields(preflight, module),
                genome_build=module.genome_build.value,
                population=module.population.value,
                methods=module.methods,
                standard_reference_manifest=(
                    None
                    if preflight.reference_manifest is None
                    else preflight.reference_manifest.model_dump(mode="json")
                ),
                cojo_reference_prefix=(
                    None
                    if preflight.cojo_reference_prefix is None
                    else str(preflight.cojo_reference_prefix)
                ),
                cojo_reference_identifier_observation=(
                    preflight.cojo_reference_identifier_observation
                ),
            )
        step_number += 1

        outputs: dict[str, Any] = {
            "ld_clump_region": None,
            "ld_clump_standard": None,
            "ld_clump_cojo": None,
            "methods": list(module.methods),
            "genome_build": module.genome_build.value,
            "population": module.population.value,
            "canonical_log": str(log_path),
            "resolved_configuration": str(resolved_path),
        }
        if "region" in module.methods:
            with logger.step(
                step_number,
                total_steps,
                "Run annotated-region LD pruning",
                "ld_clump_by_regions",
            ) as step:
                region_result = ld_clump_by_regions(
                    sumstat_vcf=str(preflight.vcf),
                    output_directory=str(preflight.output_directory),
                    dataset_id=preflight.dataset_id,
                    population=module.population.value,
                    bcftools=preflight.bcftools,
                    threads=configuration.execution.threads,
                    include_report_tables=True,
                    configuration=module,
                    logger=logger,
                )
                outputs["ld_clump_region"] = region_result
                step.set_rows(region_result.get("significant_blocks"))
                step.outcome(
                    "Selected the strongest association in each annotated LD block.",
                    fields=[
                        (
                            "count", "Annotated variants",
                            region_result.get("annotated_variants"),
                        ),
                        ("genetic", "LD blocks", region_result.get("ld_blocks")),
                        (
                            "success", "Significant blocks",
                            region_result.get("significant_blocks"),
                        ),
                        (
                            (
                                "warning"
                                if region_result.get(
                                    "genome_wide_significant_outside_ld_regions"
                                )
                                else "info"
                            ),
                            "Significant outside LD regions",
                            region_result.get(
                                "genome_wide_significant_outside_ld_regions"
                            ),
                        ),
                    ],
                    annotated_variants=region_result.get("annotated_variants"),
                    ld_blocks=region_result.get("ld_blocks"),
                    significant_blocks=region_result.get("significant_blocks"),
                    genome_wide_significant_outside_ld_regions=(
                        region_result.get(
                            "genome_wide_significant_outside_ld_regions"
                        )
                    ),
                )
            step_number += 1
        if "standard" in module.methods:
            with logger.step(
                step_number,
                total_steps,
                "Run FUMA-style standard LD clumping",
                "ld_clump_standard",
            ) as step:
                standard_result = ld_clump_standard(
                    vcf_path=preflight.vcf,
                    output_directory=preflight.output_directory,
                    dataset_id=preflight.dataset_id,
                    threads=configuration.execution.threads,
                    memory_gb=configuration.execution.memory_gb,
                    ld_folder=preflight.reference_directory,
                    pop=module.population.value,
                    bcftools_bin=preflight.bcftools,
                    tabix_bin=preflight.tabix,
                    screen_label_width=(
                        configuration.logging.terminal_label_width
                    ),
                    emit_terminal_summary=False,
                    include_report_tables=True,
                    configuration=module,
                    reference_manifest=preflight.reference_manifest,
                    logger=logger,
                )
                outputs["ld_clump_standard"] = standard_result
                step.set_rows(standard_result.get("genomic_risk_loci"))
                step.outcome(
                    "Completed both LD clumping passes and genomic-locus definition.",
                    fields=[
                        (
                            "success",
                            "GWS variants passing LD-reference checks",
                            standard_result.get("significant_variants"),
                        ),
                        (
                            "genetic", "Independent significant SNPs",
                            standard_result.get("independent_significant_snps"),
                        ),
                        (
                            "genetic", "Lead SNPs",
                            standard_result.get("lead_snps"),
                        ),
                        (
                            "genetic", "Genomic risk loci",
                            standard_result.get("genomic_risk_loci"),
                        ),
                    ],
                    status=standard_result.get("status"),
                    reference_coverage_status=standard_result.get(
                        "reference_coverage_status"
                    ),
                    significant_variants=standard_result.get(
                        "significant_variants"
                    ),
                    independent_significant_snps=standard_result.get(
                        "independent_significant_snps"
                    ),
                    lead_snps=standard_result.get("lead_snps"),
                    genomic_risk_loci=standard_result.get("genomic_risk_loci"),
                )
            step_number += 1
        if "cojo-slct" in module.methods:
            with logger.step(
                step_number,
                total_steps,
                "Run GCTA-COJO selection and define physical loci",
                "run_cojo_slct_clumping",
            ) as step:
                cojo_module = configuration.modules.gcta_cojo
                chromosome_parallel_requested = (
                    chromosome_parallel_slct_requested(cojo_module)
                )
                step.observed(
                    "cojo_execution_setup",
                    strategy=(
                        "chromosome_parallel_after_memory_pilot"
                        if chromosome_parallel_requested
                        else "single_gcta_process"
                    ),
                    total_thread_budget=configuration.execution.threads,
                    threads_per_gcta_process=(
                        cojo_module.chromosome_execution.threads_per_worker
                        if chromosome_parallel_requested
                        else configuration.execution.threads
                    ),
                    simultaneous_gcta_processes=(
                        "resolve_after_memory_pilot"
                        if chromosome_parallel_requested
                        else 1
                    ),
                )
                from postgwas.modules.formatting.contracts import required_formats
                from postgwas.modules.formatting.service import run_formatter_direct
                from postgwas.modules.gcta_cojo.service import run_gcta_cojo_direct

                values = {
                    "dataset_id": preflight.dataset_id,
                    "population": module.population.value,
                }
                subordinate_dataset_id = module.cojo.subordinate_dataset_id.format(
                    **values,
                )
                formatter_directory = configured_output_path(
                    preflight.output_directory,
                    module.output_layout.cojo_formatter_directory,
                    error_type=LDClumpingError,
                    **values,
                )
                cojo_directory = configured_output_path(
                    preflight.output_directory,
                    module.output_layout.cojo_root_directory,
                    error_type=LDClumpingError,
                    **values,
                )
                selected_signals_path = configured_output_path(
                    preflight.output_directory,
                    module.output_layout.cojo_selected_signals,
                    error_type=LDClumpingError,
                    **values,
                )
                loci_path = configured_output_path(
                    preflight.output_directory,
                    module.output_layout.cojo_loci,
                    error_type=LDClumpingError,
                    **values,
                )
                exclusion_path = configured_output_path(
                    preflight.output_directory,
                    module.output_layout.cojo_excluded_variants,
                    error_type=LDClumpingError,
                    **values,
                )

                formatter_result = None
                if ctx is not None:
                    formatter_result = (
                        (ctx.get("formatter") or {}).get("gcta_gene")
                    )
                formatter_configuration = configuration.model_copy(deep=True)
                formatter_configuration.run.output_directory = formatter_directory
                formatter_configuration.run.dataset_id = subordinate_dataset_id
                formatter_module = formatter_configuration.modules.formatting
                formatter_module.runtime = formatter_module.runtime.model_copy(update={
                    "resolved_config_file": (
                        module.cojo.formatter_resolved_config_file
                    ),
                    "completion_manifest_file": (
                        module.cojo.formatter_completion_manifest_file
                    ),
                })
                formatter_module.formats = required_formats(
                    formatter_module, modules=("gcta_cojo",),
                )
                target_types = dict(
                    getattr(args, "variant_id_types", None)
                    or formatter_module.variant_identifiers.target_types
                )
                formatter_module.variant_identifiers = (
                    formatter_module.variant_identifiers.model_copy(update={
                        "target_types": target_types,
                    })
                )
                if formatter_result is None:
                    formatter_args = argparse.Namespace(**vars(args))
                    formatter_args.vcf = str(preflight.vcf)
                    formatter_result = run_formatter_direct(
                        formatter_args,
                        configuration=formatter_configuration,
                        emit_terminal_summary=False,
                    )["gcta_gene"]
                summary_file = formatter_result.get(
                    "summary_statistics_input_file"
                )
                if not summary_file:
                    raise LDClumpingError(
                        "The canonical formatter did not return the GCTA .ma "
                        "artifact required by cojo-slct."
                    )

                cojo_configuration = formatter_configuration.model_copy(deep=True)
                cojo_configuration.run.output_directory = cojo_directory
                cojo_module = cojo_configuration.modules.gcta_cojo
                cojo_module.output_layout = cojo_module.output_layout.model_copy(
                    update={
                        "resolved_config_file": module.cojo.gcta_resolved_config_file,
                    }
                )
                exclusion = prepare_cojo_exclusion_list(
                    reference_prefix=preflight.cojo_reference_prefix,
                    destination=exclusion_path,
                    ld_module=module,
                    cojo_module=cojo_module,
                )
                if exclusion["path"] is not None:
                    cojo_module.inputs = cojo_module.inputs.model_copy(update={
                        "exclude_snps": exclusion["path"],
                    })
                cojo_args = argparse.Namespace(**vars(args))
                cojo_args.gcta_cojo_input_file = str(summary_file)
                cojo_args.variant_id_observations = getattr(
                    args, "variant_id_observations", {},
                )
                cojo_context = RunContext({
                    "formatter": {"gcta_gene": formatter_result},
                })
                gcta_result = run_gcta_cojo_direct(
                    cojo_args,
                    cojo_context,
                    configuration=cojo_configuration,
                    emit_terminal_summary=False,
                    reference_validation=preflight.cojo_reference_validation,
                    parallel_slct_output_contract=(
                        module.cojo.parallel_output_contract
                    ),
                )
                normalized = gcta_result.artifacts["normalized_results"].path
                cojo_result = group_cojo_selected_signals(
                    normalized_result=normalized,
                    selected_signals_path=selected_signals_path,
                    loci_path=loci_path,
                    ld_module=module,
                    cojo_module=cojo_module,
                )
                cojo_result.update({
                    "summary_variants": gcta_result.metrics.get(
                        "summary_variants"
                    ),
                    "reference_overlap_variants": gcta_result.metrics.get(
                        "reference_overlap_variants"
                    ),
                    "reference_overlap_fraction": gcta_result.metrics.get(
                        "reference_overlap_fraction"
                    ),
                    "reference_missing_variants": gcta_result.metrics.get(
                        "reference_missing_variants"
                    ),
                    "reference_allele_mismatch_variants": (
                        gcta_result.metrics.get(
                            "reference_allele_mismatch_variants"
                        )
                    ),
                    "reference_samples": gcta_result.metrics.get(
                        "reference_samples"
                    ),
                    "gcta_version": gcta_result.metrics.get("gcta_version"),
                    "execution_strategy": gcta_result.metrics.get(
                        "execution_strategy"
                    ),
                    "parallel_output_contract": gcta_result.metrics.get(
                        "parallel_output_contract"
                    ),
                    "conditional_reconstruction_status": (
                        gcta_result.metrics.get(
                            "conditional_reconstruction_status"
                        )
                    ),
                    "gcta_model_window_kb": cojo_module.analysis.window_kb,
                    "gcta_significance_threshold": (
                        cojo_module.analysis.significance_threshold
                    ),
                    "gcta_collinearity_cutoff": (
                        cojo_module.analysis.collinearity_cutoff
                    ),
                    "gcta_frequency_difference_max": (
                        cojo_module.analysis.frequency_difference_max
                    ),
                    "gcta_reference_maf_min": (
                        cojo_module.analysis.reference_maf_min
                    ),
                    "gcta_chromosome": cojo_module.analysis.chromosome,
                    "gcta_command": gcta_result.metrics.get("command"),
                    "warnings": list(gcta_result.warnings),
                    "exclusions": exclusion,
                    "formatter_output": str(summary_file),
                    "cojo_output_directory": str(cojo_directory),
                    "reference_prefix": str(preflight.cojo_reference_prefix),
                })
                cojo_result["output_files"].update({
                    "formatter_ma": str(summary_file),
                })
                if exclusion["path"] is not None:
                    cojo_result["output_files"]["excluded_variants"] = (
                        exclusion["path"]
                    )
                for name, artifact in gcta_result.artifacts.items():
                    output_name = (
                        name if name.startswith("gcta_") else "gcta_%s" % name
                    )
                    cojo_result["output_files"][output_name] = str(artifact.path)
                outputs["ld_clump_cojo"] = cojo_result
                step.set_rows(cojo_result["genomic_loci"])
                step.outcome(
                    "GCTA selected conditionally independent signals; physical "
                    "locus grouping did not alter the fitted model.",
                    fields=[
                        (
                            (
                                "warning"
                                if cojo_result.get("reference_missing_variants")
                                or cojo_result.get(
                                    "reference_allele_mismatch_variants"
                                )
                                else "success"
                            ),
                            "GWAS variants matched to LD reference",
                            format_cojo_reference_match(cojo_result),
                        ),
                        (
                            "genetic", "COJO-selected signals",
                            cojo_result.get("selected_signals"),
                        ),
                        (
                            "genetic", "Physical loci",
                            cojo_result.get("genomic_loci"),
                        ),
                        (
                            "info", "COJO output scope",
                            format_cojo_output_scope(cojo_result),
                        ),
                    ],
                    selected_signals=cojo_result.get("selected_signals"),
                    genomic_loci=cojo_result.get("genomic_loci"),
                    merge_distance_bp=cojo_result.get("merge_distance_bp"),
                    summary_variants=cojo_result.get("summary_variants"),
                    reference_overlap_variants=cojo_result.get(
                        "reference_overlap_variants"
                    ),
                    reference_overlap_fraction=cojo_result.get(
                        "reference_overlap_fraction"
                    ),
                    reference_missing_variants=cojo_result.get(
                        "reference_missing_variants"
                    ),
                    reference_allele_mismatch_variants=cojo_result.get(
                        "reference_allele_mismatch_variants"
                    ),
                    reference_samples=cojo_result.get("reference_samples"),
                    execution_strategy=cojo_result.get("execution_strategy"),
                    parallel_output_contract=cojo_result.get(
                        "parallel_output_contract"
                    ),
                    conditional_reconstruction_status=cojo_result.get(
                        "conditional_reconstruction_status"
                    ),
                    gcta_command=cojo_result.get("gcta_command"),
                    warnings=cojo_result.get("warnings"),
                )
            step_number += 1
        region_result = outputs.get("ld_clump_region") or {}
        standard_result = outputs.get("ld_clump_standard") or {}
        cojo_result = outputs.get("ld_clump_cojo") or {}
        if standard_result.get("status") == "partial_reference":
            outputs["status"] = "partial_reference"
        else:
            outputs["status"] = "completed"

        report_paths = _report_paths(
            preflight.output_directory, module, preflight.dataset_id,
        )
        outputs.update(
            summary_csv=str(report_paths["summary_csv"]),
            html_report=str(report_paths["html_report"]),
        )
        with logger.step(
            step_number,
            total_steps,
            "Save CSV summary and detailed HTML report",
            "write_ld_clumping_reports",
        ) as step:
            summary = build_ld_clumping_summary(
                dataset_id=preflight.dataset_id,
                input_vcf=preflight.vcf,
                reference_directory=preflight.reference_directory,
                module=module,
                methods=module.methods,
                status=outputs["status"],
                region_result=outputs.get("ld_clump_region"),
                standard_result=outputs.get("ld_clump_standard"),
                cojo_result=outputs.get("ld_clump_cojo"),
                output_paths=_reported_output_paths(outputs),
                threads=configuration.execution.threads,
                memory_gb=configuration.execution.memory_gb,
            )
            csv_path = write_ld_clumping_summary_csv(
                summary, report_paths["summary_csv"],
            )
            html_path = write_ld_clumping_html_report(
                summary, report_paths["html_report"],
            )
            # The HTML has consumed these validated presentation rows. Keep
            # them out of returned results, canonical logs, and checkpoints.
            region_result.pop(
                "_significant_outside_ld_region_details", None,
            )
            standard_result.pop("_report_tables", None)
            standard_result.pop("_reference_exclusion_details", None)
            cojo_result.pop("_report_tables", None)
            step.set_rows(len(module.methods))
            step.outcome(
                "Published one reconciled summary across terminal, CSV, and HTML.",
                fields=[
                    ("success", "Summary CSV", csv_path.name),
                    ("success", "Detailed HTML report", html_path.name),
                ],
                summary_csv=str(csv_path),
                html_report=str(html_path),
                report_source="validated_in_memory_method_results",
            )

        terminal_summary = render_ld_clumping_summary(
            summary,
            label_width=configuration.logging.terminal_label_width,
        )
        logger.log("RESULT", terminal_summary, wrap=False, screen=False)
        if configuration.logging.show_screen:
            print(terminal_summary)

        if outputs["status"] == "partial_reference":
            logger.record(
                "WARNING",
                "ld_clumping_partial_reference",
                methods=module.methods,
                skipped_chromosomes=standard_result.get("skipped_chromosomes"),
                skipped_significant_variants=standard_result.get(
                    "skipped_significant_variants"
                ),
                outputs=outputs,
            )
        else:
            logger.record(
                "DONE",
                "ld_clumping_completed",
                methods=module.methods,
                outputs=outputs,
            )
        if ctx is not None:
            # One key, matching pipeline.runners and the fine-mapping runner.
            ctx["ld_clump"] = outputs
        return outputs
    except BaseException as exc:
        logger.error("LD clumping failed: %s: %s" % (type(exc).__name__, exc))
        _quarantine_changed_outputs(
            owned_paths,
            before_outputs,
            logger,
            validated_failure_outputs=tuple(
                getattr(exc, "validated_failure_outputs", ())
            ),
        )
        raise
    finally:
        logger.close()


__all__ = [
    "LDClumpingError",
    "LDClumpingInputValidation",
    "LDClumpingPreflight",
    "LDReferenceManifest",
    "preflight_ld_clumping",
    "resolve_ld_clumping_configuration",
    "run_ld_clump_direct",
    "validate_ld_clumping_configuration",
]
