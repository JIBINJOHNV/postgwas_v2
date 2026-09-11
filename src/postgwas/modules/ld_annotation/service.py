"""Configuration boundary and execution service for LD-block annotation."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from postgwas.config import load_run_configuration_for_module
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.config.models import PostGWASConfig
from postgwas.core.errors import ConfigurationError
from postgwas.core.paths import (
    configured_output_path,
    require_nonempty_file,
    validate_filename_component,
)
from postgwas.core.pipeline_logging import PipelineLogger
from postgwas.core.preflight import (
    PipelinePreflightEvidence,
    pipeline_preflight_evidence,
    require_pipeline_input_vcf,
)
from postgwas.core.required_arguments import (
    RequiredArgument,
    require_resolved_arguments,
)
from postgwas.core.ui import StageProgress
from postgwas.modules.ld_annotation.annot_ldblock import (
    annotate_ldblocks,
    validate_ld_block_references,
)
from postgwas.modules.ld_annotation.reporting import render_ld_annotation_summary


MODULE_CLI_OVERRIDES = {
    "vcf": "inputs.vcf",
    "ld_region_dir": "inputs.ld_region_dir",
    "dataset_id": "inputs.dataset_id",
    "output_directory": "output_directory",
    "ld_block_populations": "populations",
}

GLOBAL_CLI_OVERRIDES = {
    "dataset_id": "run.dataset_id",
    "output_directory": "run.output_directory",
    "threads": "execution.threads",
    "memory_gb": "execution.memory_gb",
    "seed": "execution.random_seed",
    "bcftools": "resources.executables.bcftools",
}


@dataclass(frozen=True)
class LDAnnotationPreflight:
    """Resolved values proven safe before invoking annotation tools."""

    configuration: PostGWASConfig
    vcf: Path
    ld_region_dir: Path
    output_directory: Path
    dataset_id: str


def resolve_ld_annotation_configuration(
    args: argparse.Namespace,
) -> PostGWASConfig:
    """Resolve explicit CLI values over module-only or complete run YAML."""
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
        "ld_annotation",
        getattr(args, "run_config", None),
        module_overrides=module_overrides,
        global_overrides=global_overrides,
    )


def validate_ld_annotation_configuration(
    args: argparse.Namespace,
    *,
    configuration: PostGWASConfig | None = None,
) -> LDAnnotationPreflight:
    """Validate resolved run inputs before tools run or outputs are created."""
    configuration = configuration or resolve_ld_annotation_configuration(args)
    module = configuration.modules.ld_annotation
    require_resolved_arguments(
        (
            RequiredArgument(
                "--vcf",
                "modules.ld_annotation.inputs.vcf",
                module.inputs.vcf,
            ),
            RequiredArgument(
                "--ld-region-dir",
                "modules.ld_annotation.inputs.ld_region_dir",
                module.inputs.ld_region_dir,
            ),
        )
    )

    vcf = require_nonempty_file(
        module.inputs.vcf,
        "LD-annotation input VCF",
        error_type=ConfigurationError,
    )
    ld_region_value = module.inputs.ld_region_dir
    if ld_region_value is None:  # Protected by require_resolved_arguments().
        raise AssertionError("Resolved LD-block directory requirement was lost")
    ld_region_dir = ld_region_value.expanduser().resolve()
    if not ld_region_dir.is_dir():
        raise ConfigurationError(
            "LD-block directory does not exist or is not a directory: %s"
            % ld_region_dir
        )

    dataset_id = validate_filename_component(
        module.inputs.dataset_id or configuration.run.dataset_id,
        "LD-annotation dataset_id",
        error_type=ConfigurationError,
    )
    configured_output = module.output_directory or configuration.run.output_directory
    output_directory = Path(configured_output).expanduser().resolve()
    return LDAnnotationPreflight(
        configuration=configuration,
        vcf=vcf,
        ld_region_dir=ld_region_dir,
        output_directory=output_directory,
        dataset_id=dataset_id,
    )


def preflight_ld_annotation(
    args: argparse.Namespace,
    *,
    preflight_evidence=None,
) -> PipelinePreflightEvidence:
    """Pipeline preflight for resolved LD-annotation inputs."""
    entry = require_pipeline_input_vcf(preflight_evidence)
    resources = validate_ld_annotation_configuration(args)
    module = resources.configuration.modules.ld_annotation
    validate_ld_block_references(
        resources.ld_region_dir,
        genome_build=entry["indexed"].genome_build,
        populations=tuple(population.value for population in module.populations),
        bed_filename_template=module.bed_filename_template,
        vcf_contigs=entry["indexed"].contigs,
    )
    return pipeline_preflight_evidence(
        "annot_ldblock",
        preflight_evidence,
        resources=resources,
    )


def run_annot_ldblock(
    args: argparse.Namespace,
    ctx: dict[str, Any] | None = None,
    *,
    configuration: PostGWASConfig | None = None,
) -> dict[str, Any]:
    """Annotate population LD blocks using one resolved configuration."""
    preflight = validate_ld_annotation_configuration(
        args,
        configuration=configuration,
    )
    configuration = preflight.configuration
    module = configuration.modules.ld_annotation
    progress = StageProgress(
        "LD-block annotation progress",
        enabled=configuration.logging.show_progress,
        outcome_label_width=configuration.logging.terminal_label_width,
    )
    log_path = configured_output_path(
        preflight.output_directory,
        module.canonical_log_filename_template,
        dataset_id=preflight.dataset_id,
        error_type=ConfigurationError,
    )
    logger = PipelineLogger(
        sample_id=preflight.dataset_id,
        scope="run",
        log_dir=str(log_path.parent),
        level=configuration.logging.file_level,
        screen_level=configuration.logging.console_level,
        log_path=str(log_path),
        stage_progress=progress,
    )
    logger.record(
        "PARAM",
        "ld_annotation_run",
        dataset_id=preflight.dataset_id,
        input_vcf=str(preflight.vcf),
        ld_region_dir=str(preflight.ld_region_dir),
        output_directory=str(preflight.output_directory),
        populations=tuple(
            population.value for population in module.populations
        ),
        threads=configuration.execution.threads,
        bed_filename_template=module.bed_filename_template,
        info_field_template=module.info_field_template,
        html_report_filename_template=module.html_report_filename_template,
    )
    try:
        outputs = annotate_ldblocks(
            vcf_path=str(preflight.vcf),
            output_directory=str(preflight.output_directory),
            ld_dir=str(preflight.ld_region_dir),
            genome_build_header=(
                configuration.modules.harmonisation.vcf_processing.genome_build_header
            ),
            supported_genome_builds=list(configuration.resources.genomes),
            bed_filename_template=module.bed_filename_template,
            info_field_template=module.info_field_template,
            info_description_template=module.info_description_template,
            output_filename_template=module.output_filename_template,
            summary_filename_template=module.summary_filename_template,
            html_report_filename_template=(
                module.html_report_filename_template
            ),
            bcftools_bin=configuration.resources.executables.bcftools,
            provenance_headers=(
                configuration.modules.formatting.input_contract.provenance_headers.model_dump()
            ),
            populations=tuple(module.populations),
            threads=configuration.execution.threads,
            dataset_id=preflight.dataset_id,
            logger=logger,
            stage_progress=progress,
        )
    finally:
        logger.close()
    outputs["log_file"] = log_path
    print(
        render_ld_annotation_summary(
            outputs["summary"],
            outputs["summary_file"],
            label_width=configuration.logging.terminal_label_width,
            separator_column=progress.outcome_separator_column,
            annotated_vcf=outputs["annotated_vcf"],
            annotated_index=outputs["annotated_index"],
            html_report=outputs["html_report"],
            log_file=log_path,
        ),
        end="",
    )
    if ctx is not None:
        ctx["annot_ldblock"] = outputs
    return outputs


__all__ = [
    "LDAnnotationPreflight",
    "preflight_ld_annotation",
    "resolve_ld_annotation_configuration",
    "run_annot_ldblock",
    "validate_ld_annotation_configuration",
]
