"""Configuration boundary and execution service for GWAS-VCF filtering."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from postgwas.cli.compute import memory_limit
from postgwas.config import load_run_configuration_for_module
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.core.paths import configured_output_path, resolve_executable
from postgwas.core.pipeline_logging import write_log_record
from postgwas.modules.filtering.sumstat_filter import filter_gwas_vcf_bcftools


MODULE_CLI_OVERRIDES = {
    "vcf": "inputs.vcf",
    "dataset_id": "inputs.dataset_id",
    "output_directory": "output_directory",
    "genome_build": "genome_build",
    "minimum_neglog10_p": "minimum_neglog10_p",
    "minimum_maf": "maf_min",
    "reference_af_column": "reference_population_tag",
    "maximum_af_difference": "frequency_difference_max",
    "minimum_info": "info_min",
    "maximum_info": "info_max",
    "missing_info_action": "missing_info_action",
    "include_indels": "include_indels",
    "remove_palindromic": "remove_palindromic",
    "palindromic_af_lower": "palindromic_lower",
    "palindromic_af_upper": "palindromic_upper",
    "remove_mhc": "remove_mhc",
    "mhc_chrom": "mhc.chromosome",
    "mhc_start": "mhc.start",
    "mhc_end": "mhc.end",
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
        suffix = " or provide --bcftools PATH" if label == "bcftools" else ""
        raise RuntimeError(
            "%s. Configure resources.executables.%s in --run-config%s."
            % (exc, label, suffix)
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
            module.output_layout.log_file,
            dataset_id=dataset_id,
            genome_build=module.genome_build.value,
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

    outputs = filter_gwas_vcf_bcftools(
        vcf_path=str(vcf),
        output_folder=str(output_directory),
        output_prefix=dataset_id,
        genome_build=module.genome_build.value,
        genome_build_header_tokens={
            build.value: tokens
            for build, tokens in module.genome_build_header_tokens.items()
        },
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
        mhc_chrom=module.mhc.chromosome,
        mhc_start=module.mhc.start,
        mhc_end=module.mhc.end,
        threads=configuration.execution.threads,
        max_mem=memory_limit(configuration.execution.memory_gb),
        lp_missing=module.missing_pvalue_action,
        af_missing=module.missing_af_action,
        empty_expression=module.empty_expression_action,
        sort_output=module.sort_output,
        vcf_fields=module.vcf_fields.model_dump(),
        output_layout=module.output_layout.model_dump(),
        bcftools_bin=bcftools,
        tabix_bin=tabix,
        bash_bin=bash,
        resolved_configuration={
            "filtering": module.model_dump(mode="json"),
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
    )
    if ctx is not None:
        ctx["sumstat_filter"] = outputs
    return outputs
