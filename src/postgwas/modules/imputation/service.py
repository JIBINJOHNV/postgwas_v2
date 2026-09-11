"""Validated configuration and execution boundary for PRED-LD imputation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from postgwas.config import load_run_configuration_for_module
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.core.execution.runtime import safe_thread_count
from postgwas.core.preflight import (
    PipelinePreflightEvidence,
    PreflightFileIdentity,
    capture_preflight_file_identities,
    pipeline_preflight_evidence,
    require_pipeline_input_vcf,
    require_unchanged_preflight_files,
)
from postgwas.core.required_arguments import (
    RequiredArgument,
    require_resolved_arguments,
)
from postgwas.modules.harmonisation.cli import (
    harmonisation_executable_requirements,
    run_harmonisation,
)
from postgwas.modules.harmonisation.vcf_processing import require_binaries
from postgwas.modules.imputation.engines.pred_ld.pred_ld_runner import (
    pred_ld_script_path,
    process_pred_ld_results_all_parallel,
    run_pred_ld_parallel,
    validate_pred_ld_reference,
)


@dataclass(frozen=True)
class PredLDPipelineResources:
    """External PRED-LD and harmonisation resources validated before formatting."""

    configuration: Any
    reference_directory: Path
    reference_files: tuple[Path, ...]
    pred_ld_script: Path
    harmonisation_executables: Mapping[str, str]
    resource_root: Path
    file_identities: tuple[PreflightFileIdentity, ...]


def resolve_imputation_configuration(args):
    """Resolve canonical YAML and only explicitly supplied CLI overrides."""
    module_overrides = explicit_overrides(args, {
        "imputation_engine": "engine",
        "pred_ld_input_directory": "input_directory",
        "imputation_ld_reference": "ld_reference_directory",
        "genome_build": "genome_build",
        "population": "population",
        "imputation_r2_threshold": "engines.pred_ld.minimum_r2",
        "imputation_minimum_maf": "engines.pred_ld.minimum_maf",
        "ref": "engines.pred_ld.mode",
        "corr_method": "engines.pred_ld.correlation_method",
    })
    global_overrides = explicit_overrides(args, {
        "dataset_id": "run.dataset_id",
        "output_directory": "run.output_directory",
        "resource_directory": "resources.root",
        "threads": "execution.threads",
        "memory_gb": "execution.memory_gb",
        "seed": "execution.random_seed",
        "resume": "run.resume",
        "overwrite": "run.overwrite",
    })
    return load_run_configuration_for_module(
        "imputation",
        getattr(args, "run_config", None),
        module_overrides=module_overrides,
        global_overrides=global_overrides,
    )


def _require_imputation_arguments(
    configuration,
    *,
    include_generated_input: bool,
) -> None:
    """Report every independently missing resolved imputation value together."""
    module = configuration.modules.imputation
    requirements = [
        RequiredArgument(
            "--imputation-ld-reference",
            "modules.imputation.ld_reference_directory",
            module.ld_reference_directory,
        ),
        RequiredArgument(
            "--genome-build",
            "modules.imputation.genome_build",
            module.genome_build,
        ),
        RequiredArgument(
            "--resource-directory",
            "resources.root",
            configuration.resources.root,
        ),
    ]
    if include_generated_input:
        requirements.insert(0, RequiredArgument(
            "--pred-ld-input-directory",
            "modules.imputation.input_directory",
            module.input_directory,
        ))
    require_resolved_arguments(requirements)


def _validate_imputation_resources(configuration) -> PredLDPipelineResources:
    """Validate exact PRED-LD files and downstream harmonisation executables."""
    module = configuration.modules.imputation
    _require_imputation_arguments(configuration, include_generated_input=False)
    resource_root = Path(configuration.resources.root).expanduser().resolve()
    if not resource_root.is_dir():
        raise ValueError(
            "PostGWAS resource directory does not exist or is not a directory: %s"
            % resource_root
        )
    pred_ld = module.engines.pred_ld
    safe_thread_count(
        configuration.execution.threads,
        pred_ld.memory_gb_per_worker,
        available_ram_gb=configuration.execution.memory_gb,
        enforce_memory_budget=True,
        reporter=None,
    )
    reference_directory, reference_files = validate_pred_ld_reference(
        module.ld_reference_directory,
        population=module.population.value,
        chromosomes=configuration.modules.formatting.chromosomes,
        mode=pred_ld.mode,
        subdirectory_template=pred_ld.reference_subdirectory_template,
        file_template=pred_ld.reference_file_template,
        file_kinds=pred_ld.reference_file_kinds,
    )
    script = pred_ld_script_path()
    executables = require_binaries(
        harmonisation_executable_requirements(configuration),
        plugins=(
            configuration.modules.harmonisation.vcf_processing.liftover_plugin,
        ),
    )
    resource_files = [*reference_files, script]
    resource_files.extend(
        Path(value).expanduser().resolve()
        for value in executables.values()
        if Path(value).expanduser().is_file()
    )
    return PredLDPipelineResources(
        configuration=configuration,
        reference_directory=reference_directory,
        reference_files=reference_files,
        pred_ld_script=script,
        harmonisation_executables=executables,
        resource_root=resource_root,
        file_identities=capture_preflight_file_identities(
            resource_files,
            label="PRED-LD resource",
        ),
    )


def preflight_imputation_pipeline(
    args,
    *,
    preflight_evidence=None,
) -> PipelinePreflightEvidence:
    """Validate the entry VCF and imputation resources before formatting."""
    entry_vcf = require_pipeline_input_vcf(preflight_evidence)
    configuration = resolve_imputation_configuration(args)
    module = configuration.modules.imputation
    _require_imputation_arguments(configuration, include_generated_input=False)
    observed_build = str(entry_vcf["harmonised"]["genome_build"])
    configured_build = module.genome_build.value
    if configured_build != observed_build:
        raise ValueError(
            "The harmonised GWAS-VCF declares genome build %s, but the PRED-LD "
            "reference is declared as %s. Summary-statistic positions and the "
            "LD reference must use the same build."
            % (observed_build, configured_build)
        )
    resources = _validate_imputation_resources(configuration)
    return pipeline_preflight_evidence(
        "imputation",
        preflight_evidence,
        resources=resources,
        deferred_checks=(
            "Validate the formatter-created per-chromosome PRED-LD tables.",
            "Validate imputed outputs and their exact input-dependent "
            "harmonisation reference maps before post-imputation workers start.",
        ),
    )


def _imputation_pipeline_execution_configuration(args, resources):
    """Retain preflight settings while applying generated input and stage paths."""
    module = resources.configuration.modules.imputation.model_copy(update={
        "input_directory": Path(args.pred_ld_input_directory).expanduser().resolve(),
    })
    run = resources.configuration.run.model_copy(update={
        "output_directory": Path(args.output_directory).expanduser().resolve(),
    })
    modules = resources.configuration.modules.model_copy(update={
        "imputation": module,
    })
    return resources.configuration.model_copy(
        update={"run": run, "modules": modules},
        deep=True,
    )


def _restore_argument(args, name: str, existed: bool, value) -> None:
    if existed:
        setattr(args, name, value)
    elif hasattr(args, name):
        delattr(args, name)


def run_sumstat_imputation_direct(
    args,
    ctx=None,
    *,
    configuration=None,
    pipeline_resources: PredLDPipelineResources | None = None,
):
    """Run PRED-LD, consolidate its outputs, and re-harmonise the result."""
    if pipeline_resources is not None and configuration is not None:
        raise ValueError(
            "Pass either pipeline_resources or configuration to imputation, not both."
        )
    configuration = (
        _imputation_pipeline_execution_configuration(args, pipeline_resources)
        if pipeline_resources is not None
        else configuration or resolve_imputation_configuration(args)
    )
    _require_imputation_arguments(configuration, include_generated_input=True)
    resources = pipeline_resources or _validate_imputation_resources(configuration)
    require_unchanged_preflight_files(
        resources.file_identities,
        label="PRED-LD resource",
    )

    module = configuration.modules.imputation
    pred_ld = module.engines.pred_ld
    dataset_id = configuration.run.dataset_id
    output_directory = Path(configuration.run.output_directory).expanduser().resolve()
    input_directory = Path(module.input_directory).expanduser().resolve()
    if not input_directory.is_dir():
        raise ValueError(
            "PRED-LD input directory does not exist: %s" % input_directory
        )
    run_pred_ld_parallel(
        predld_input_dir=str(input_directory),
        pred_ld_ref=str(resources.reference_directory),
        output_folder=str(output_directory),
        output_prefix=dataset_id,
        chromosomes=configuration.modules.formatting.chromosomes,
        r2threshold=pred_ld.minimum_r2,
        maf=pred_ld.minimum_maf,
        population=module.population.value,
        ref=pred_ld.mode,
        threads=configuration.execution.threads,
        memory_gb=configuration.execution.memory_gb,
        pred_ld_script=resources.pred_ld_script,
        python_executable=resources.harmonisation_executables["python"],
        memory_gb_per_worker=pred_ld.memory_gb_per_worker,
        free_memory_threshold_gb=pred_ld.free_memory_threshold_gb,
        free_memory_threshold_fraction=(
            pred_ld.free_memory_threshold_fraction
        ),
        memory_poll_seconds=pred_ld.memory_poll_seconds,
        worker_poll_seconds=pred_ld.worker_poll_seconds,
        large_chromosomes=pred_ld.large_chromosomes,
        preferred_chromosome_order=pred_ld.preferred_chromosome_order,
    )
    imputed_dataset_id = dataset_id + pred_ld.imputed_dataset_suffix
    _, _, sample_sheet_path = process_pred_ld_results_all_parallel(
        folder_path=str(output_directory),
        output_path=str(output_directory),
        output_prefix=dataset_id,
        harmonised_dataset_id=imputed_dataset_id,
        corr_method=pred_ld.correlation_method,
        threads=configuration.execution.threads,
    )

    harmonised_output = (
        output_directory.parent / module.post_harmonisation_directory
    ).resolve()
    harmonisation_run = configuration.run.model_copy(update={
        "dataset_id": imputed_dataset_id,
        "output_directory": harmonised_output,
    })
    harmonisation_configuration = configuration.model_copy(
        update={"run": harmonisation_run}, deep=True,
    )
    saved = {
        name: (hasattr(args, name), getattr(args, name, None))
        for name in ("sample_sheet", "dataset_id", "output_directory")
    }
    try:
        args.sample_sheet = sample_sheet_path
        args.dataset_id = imputed_dataset_id
        args.output_directory = str(harmonised_output)
        harmonisation_results = run_harmonisation(
            args,
            configuration=harmonisation_configuration,
            resolved_executables=resources.harmonisation_executables,
        )
    finally:
        for name, (existed, value) in saved.items():
            _restore_argument(args, name, existed, value)

    outputs = harmonisation_results.get(imputed_dataset_id)
    if not isinstance(outputs, dict):
        raise RuntimeError(
            "Post-imputation harmonisation did not return results for dataset %r"
            % imputed_dataset_id
        )
    if ctx is not None:
        ctx["imputation"] = outputs
    return outputs


__all__ = [
    "PredLDPipelineResources",
    "preflight_imputation_pipeline",
    "resolve_imputation_configuration",
    "run_sumstat_imputation_direct",
]
