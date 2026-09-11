"""Configuration boundary and execution service for GWAS-VCF filtering."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from postgwas.cli.compute import memory_limit
from postgwas.config import load_run_configuration_for_module
from postgwas.config.cli_overrides import (
    explicit_overrides,
    variant_qc_policy_override_paths,
)
from postgwas.core.paths import configured_output_path, resolve_executable
from postgwas.core.pipeline_logging import write_log_record
from postgwas.core.ui import StageProgress
from postgwas.modules.filtering.sumstat_filter import filter_gwas_vcf_bcftools


MODULE_CLI_OVERRIDES = {
    "vcf": "inputs.vcf",
    "dataset_id": "inputs.dataset_id",
    "output_directory": "output_directory",
    **variant_qc_policy_override_paths(field_overrides={
        "maximum_af_difference": "frequency_difference_max",
    }),
    "reference_af_column": "reference_population_tag",
    "mhc_chrom": "mhc_region_override.chromosome",
    "mhc_start": "mhc_region_override.start",
    "mhc_end": "mhc_region_override.end",
    "write_soft_filter_vcf": "write_soft_filter_vcf",
}

GLOBAL_CLI_OVERRIDES = {
    "dataset_id": "run.dataset_id",
    "output_directory": "run.output_directory",
    "threads": "execution.threads",
    "memory_gb": "execution.memory_gb",
    "seed": "execution.random_seed",
    "bcftools": "resources.executables.bcftools",
}


def resolve_filtering_configuration(args: argparse.Namespace):
    """Resolve CLI, module YAML, and run YAML through one validated model."""
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
        "filtering",
        getattr(args, "run_config", None),
        module_overrides=module_overrides,
        global_overrides=global_overrides,
    )


def _resolve_filtering_executable(value: str, label: str) -> str:
    try:
        return resolve_executable(value, "%s executable" % label)
    except RuntimeError as exc:
        raise RuntimeError(
            "%s. Configure resources.executables.%s in --run-config."
            % (exc, label)
        ) from exc


def run_sumstat_filter_direct(
    args: argparse.Namespace,
    ctx: dict[str, Any] | None = None,
    *,
    configuration=None,
):
    """Filter one GWAS-VCF using the completely resolved run configuration."""
    configuration = configuration or resolve_filtering_configuration(args)
    module = configuration.modules.filtering
    provenance = configuration.modules.harmonisation.vcf_processing.provenance
    required_provenance_headers = (
        configuration.modules.formatting.input_contract.provenance_headers.model_dump()
    )
    vcf = module.inputs.vcf
    dataset_id = module.inputs.dataset_id or configuration.run.dataset_id
    configured_output = module.output_directory or configuration.run.output_directory
    if vcf is None:
        raise ValueError("A harmonised GWAS-VCF is required; provide --vcf PATH.")
    if not dataset_id:
        raise ValueError("A dataset identifier is required; provide --dataset-id ID.")
    if configured_output is None:
        raise ValueError(
            "An output directory is required; provide --output-directory PATH."
        )

    output_directory = Path(configured_output).expanduser().resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    executables = configuration.resources.executables
    try:
        bcftools = _resolve_filtering_executable(executables.bcftools, "bcftools")
        tabix = _resolve_filtering_executable(executables.tabix, "tabix")
        bash = _resolve_filtering_executable(executables.bash, "bash")
    except BaseException as exc:
        log_path = configured_output_path(
            output_directory,
            module.output_layout.preflight_log,
            dataset_id=dataset_id,
        )
        write_log_record(
            log_path,
            "ERROR",
            "Filtering executable validation failed: %s: %s"
            % (type(exc).__name__, exc),
            sample_id=dataset_id,
            file_level=configuration.logging.file_level,
            screen_level=configuration.logging.console_level,
        )
        raise

    progress = StageProgress(
        "Summary-statistics filtering progress",
        enabled=configuration.logging.show_progress,
        outcome_label_width=configuration.logging.terminal_label_width,
    )
    try:
        outputs = filter_gwas_vcf_bcftools(
            vcf_path=str(vcf),
            output_folder=str(output_directory),
            output_prefix=dataset_id,
            genome_build_header=(
                configuration.modules.harmonisation.vcf_processing.genome_build_header
            ),
            supported_genome_builds=list(configuration.resources.genomes),
            required_provenance_headers=required_provenance_headers,
            provenance_headers=provenance.headers,
            provenance_missing_value=provenance.missing_value,
            pval_cutoff=module.minimum_neglog10_p,
            maf_cutoff=module.maf_min,
            allelefreq_diff_cutoff=module.frequency_difference_max,
            info_min=module.info_min,
            info_max=module.info_max,
            info_missing=module.missing_info_action,
            external_af_name=module.reference_population_tag,
            include_indels=module.include_indels,
            exclude_palindromic=module.remove_palindromic,
            palindromic_af_lower=module.palindromic_lower,
            palindromic_af_upper=module.palindromic_upper,
            remove_mhc=module.remove_mhc,
            mhc_regions={
                build.value: region.model_dump()
                for build, region in module.mhc_regions.items()
            },
            mhc_region_override=module.mhc_region_override.model_dump(),
            threads=configuration.execution.threads,
            max_mem=memory_limit(configuration.execution.memory_gb),
            lp_missing=module.missing_pvalue_action,
            af_missing=module.missing_af_action,
            empty_expression=module.empty_expression_action,
            sort_output=module.sort_output,
            write_soft_filter_vcf=module.write_soft_filter_vcf,
            filter_reason_ids=module.filter_reason_ids.model_dump(),
            vcf_fields=module.vcf_fields.model_dump(),
            output_layout=module.output_layout.model_dump(),
            bcftools_bin=bcftools,
            tabix_bin=tabix,
            bash_bin=bash,
            resolved_configuration={
                "filtering": module.model_dump(mode="json"),
                "vcf_contract": {
                    "genome_build_header": (
                        configuration.modules.harmonisation.vcf_processing.genome_build_header
                    ),
                    "supported_genome_builds": list(configuration.resources.genomes),
                    "required_provenance_headers": required_provenance_headers,
                    "provenance_headers": dict(provenance.headers),
                    "provenance_missing_value": provenance.missing_value,
                },
                "execution": {
                    "threads": configuration.execution.threads,
                    "memory_gb": configuration.execution.memory_gb,
                },
                "executables": {
                    "bcftools": bcftools,
                    "tabix": tabix,
                    "bash": bash,
                },
            },
            report_missing_counts=module.report_missing_counts,
            terminal_label_width=configuration.logging.terminal_label_width,
            stage_progress=progress,
        )
    finally:
        progress.close()
    if ctx is not None:
        ctx["sumstat_filter"] = outputs
    return outputs
