"""Validated service boundary for LD-clumping analyses."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal
from uuid import uuid4

from pydantic import Field, field_validator
import yaml

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
from postgwas.core.pipeline_logging import PipelineLogger
from postgwas.core.required_arguments import (
    RequiredArgument,
    require_resolved_arguments,
)
from postgwas.core.vcf import (
    read_vcf_header,
    required_vcf_query_tags,
    validate_vcf_header_contract,
)
from postgwas.modules.ld_clumping.ld_prune_region import ld_clump_by_regions
from postgwas.modules.ld_clumping.ld_prune_standard import ld_clump_standard


class LDClumpingError(RuntimeError):
    """LD-clumping configuration, reference, or analysis is invalid."""


class LDReferenceManifest(StrictModel):
    """Scientific and indexing contract written by LD reference preparation."""

    format_version: int = Field(ge=1)
    genome_build: GenomeBuild
    populations: list[Population]
    orientation: Literal["symmetric_first_endpoint"]
    file_pattern: str
    columns: list[str]
    window_kb: int = Field(gt=0)
    minimum_r2: float = Field(ge=0, le=1, allow_inf_nan=False)
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
class LDClumpingPreflight:
    configuration: Any
    vcf: Path
    output_directory: Path
    dataset_id: str
    bcftools: str
    tabix: str | None
    reference_directory: Path | None
    reference_manifest: LDReferenceManifest | None


def _owned_analysis_paths(output_directory: Path, module, dataset_id: str) -> list[Path]:
    """Resolve only module-owned scientific and method-log paths."""
    values = {"dataset_id": dataset_id, "population": module.population.value}
    paths = []
    for name, pattern in module.output_layout.model_dump().items():
        if not name.startswith(("region_", "standard_")):
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
) -> None:
    """Rename outputs created or changed by a failed run so they cannot look complete."""
    for path in paths:
        if _path_identity(path) == before[path] or not path.exists():
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
    "window_kb": "window_kb",
    "merge_dist": "merge_distance_bp",
    "missing_index_action": "missing_index_action",
    "remove_mhc": "remove_mhc",
    "ld_folder": "reference.directory",
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
    return load_run_configuration_for_module(
        "ld_clumping",
        getattr(args, "run_config", None),
        module_overrides=module_overrides,
        global_overrides=global_overrides,
    )


def _read_reference_manifest(
    directory: Path,
    filename: str,
) -> LDReferenceManifest:
    path = require_nonempty_file(
        directory / filename,
        "LD reference manifest",
        error_type=LDClumpingError,
    )
    try:
        document = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise LDClumpingError(
            "Cannot read LD reference manifest %s: %s" % (path, exc)
        ) from exc
    if not isinstance(document, dict):
        raise LDClumpingError("LD reference manifest must be a YAML mapping: %s" % path)
    try:
        return LDReferenceManifest.model_validate(document)
    except Exception as exc:
        raise LDClumpingError(
            "Invalid LD reference manifest %s: %s" % (path, exc)
        ) from exc


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
    if manifest.columns != reference.columns:
        mismatches.append("columns do not match the configured seven-column contract")
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
    if not manifest.allele_order_preserved:
        mismatches.append("allele_order_preserved=false")
    if mismatches:
        raise LDClumpingError(
            "LD reference is scientifically incompatible: %s"
            % "; ".join(mismatches)
        )
    return manifest


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


def validate_ld_clumping_configuration(
    args: argparse.Namespace,
    *,
    configuration=None,
    pipeline: bool = False,
    pipeline_entry: bool = False,
    logger: PipelineLogger | None = None,
) -> LDClumpingPreflight:
    """Validate every selected method before creating scientific outputs."""
    configuration = configuration or resolve_ld_clumping_configuration(args)
    module = configuration.modules.ld_clumping
    requirements = [
        RequiredArgument(
            "--vcf", "modules.ld_clumping.inputs.vcf", module.inputs.vcf,
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
    require_resolved_arguments(requirements)
    if pipeline and "finemap" in getattr(args, "modules", ()):
        if "standard" not in module.methods:
            raise LDClumpingError(
                "Fine-mapping requires standard LD clumping; add standard to "
                "modules.ld_clumping.methods or --clumping-methods."
            )

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
        tabix = (
            resolve_executable(
                configuration.resources.executables.tabix,
                "tabix executable",
                error_type=LDClumpingError,
            )
            if "standard" in module.methods
            else None
        )
    except OSError as exc:
        raise LDClumpingError(str(exc)) from exc

    header = read_vcf_header(
        vcf, bcftools, logger=logger, error_type=LDClumpingError,
    )
    declared_build, _ = validate_vcf_header_contract(
        header=header,
        genome_build_header=(
            configuration.modules.harmonisation.vcf_processing.genome_build_header
        ),
        supported_genome_builds=list(configuration.resources.genomes),
        required_fields=required_vcf_query_tags(
            _required_query_expressions(
                module,
                include_region_annotation=not pipeline_entry,
            )
        ),
    )
    if declared_build != module.genome_build.value:
        raise LDClumpingError(
            "Input VCF genome build %s does not match configured LD-clumping "
            "genome build %s. Use build-matched VCF, LD blocks, and LD reference."
            % (declared_build, module.genome_build.value)
        )

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

    if logger is not None:
        logger.record(
            "PASS",
            "ld_clumping_preflight",
            input_vcf=str(vcf),
            genome_build=declared_build,
            population=module.population.value,
            methods=module.methods,
            reference_manifest=(
                None
                if reference_directory is None
                else str(reference_directory / module.reference.manifest_filename)
            ),
        )
    return LDClumpingPreflight(
        configuration=configuration,
        vcf=vcf,
        output_directory=output_directory,
        dataset_id=dataset_id,
        bcftools=bcftools,
        tabix=tabix,
        reference_directory=reference_directory,
        reference_manifest=manifest,
    )


def preflight_ld_clumping(args: argparse.Namespace) -> LDClumpingPreflight:
    """Pipeline preflight using the same resolved validation as direct runs."""
    return validate_ld_clumping_configuration(
        args, pipeline=True, pipeline_entry=True,
    )


def run_ld_clump_direct(
    args: argparse.Namespace,
    ctx: dict[str, Any] | None = None,
    *,
    configuration=None,
    pipeline: bool = False,
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
        screen_level="ERROR",
        log_path=str(log_path),
    )
    owned_paths = _owned_analysis_paths(
        output_directory, module, dataset_id,
    )
    before_outputs = {path: _path_identity(path) for path in owned_paths}
    try:
        preflight = validate_ld_clumping_configuration(
            args,
            configuration=configuration,
            pipeline=pipeline,
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
            modules="ld_clumping",
            resource_paths=("executables.bcftools", "executables.tabix"),
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
            remove_mhc=module.remove_mhc,
            missing_index_action=module.missing_index_action,
            threads=configuration.execution.threads,
            memory_gb=configuration.execution.memory_gb,
            resolved_configuration=str(resolved_path),
        )

        outputs: dict[str, Any] = {
            "ld_clump_region": None,
            "ld_clump_standard": None,
            "methods": list(module.methods),
            "genome_build": module.genome_build.value,
            "population": module.population.value,
            "canonical_log": str(log_path),
            "resolved_configuration": str(resolved_path),
        }
        if "region" in module.methods:
            outputs["ld_clump_region"] = ld_clump_by_regions(
                sumstat_vcf=str(preflight.vcf),
                output_directory=str(preflight.output_directory),
                dataset_id=preflight.dataset_id,
                population=module.population.value,
                bcftools=preflight.bcftools,
                threads=configuration.execution.threads,
                configuration=module,
                logger=logger,
            )
        if "standard" in module.methods:
            outputs["ld_clump_standard"] = ld_clump_standard(
                vcf_path=preflight.vcf,
                output_directory=preflight.output_directory,
                dataset_id=preflight.dataset_id,
                threads=configuration.execution.threads,
                memory_gb=configuration.execution.memory_gb,
                ld_folder=preflight.reference_directory,
                pop=module.population.value,
                bcftools_bin=preflight.bcftools,
                tabix_bin=preflight.tabix,
                configuration=module,
                logger=logger,
            )
        logger.record(
            "DONE",
            "ld_clumping_completed",
            methods=module.methods,
            outputs=outputs,
        )
        if ctx is not None:
            ctx["ld_clump_region"] = outputs["ld_clump_region"]
            ctx["ld_clump_standard"] = outputs["ld_clump_standard"]
        return outputs
    except BaseException as exc:
        logger.error("LD clumping failed: %s: %s" % (type(exc).__name__, exc))
        _quarantine_changed_outputs(owned_paths, before_outputs, logger)
        raise
    finally:
        logger.close()


__all__ = [
    "LDClumpingError",
    "LDClumpingPreflight",
    "LDReferenceManifest",
    "preflight_ld_clumping",
    "resolve_ld_clumping_configuration",
    "run_ld_clump_direct",
    "validate_ld_clumping_configuration",
]
