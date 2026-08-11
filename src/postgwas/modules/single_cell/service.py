"""Configuration, dispatch, and publication boundary for single-cell methods."""

from __future__ import annotations

import argparse
from pathlib import Path

from postgwas.config import (
    load_configuration,
    load_run_configuration_for_module,
    write_resolved_configuration,
)
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.config.models.modules.single_cell import (
    single_cell_supporting_configurations,
)
from postgwas.core.contracts import ModuleResult
from postgwas.core.paths import (
    configured_output_path,
    validate_filename_component,
)
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.modules.single_cell.errors import SingleCellError
from postgwas.modules.single_cell.methods.base import MethodRunContext
from postgwas.modules.single_cell.methods.registry import get_single_cell_methods


def resolve_single_cell_configuration(args: argparse.Namespace):
    """Resolve the effective single-cell settings and explicit CLI overrides."""
    module_overrides = explicit_overrides(
        args,
        {
            "tools": "tools",
            "magma_gene_results_file": (
                "magma_celltype.input.gene_results_file"
            ),
            "single_cell_covariates": (
                "magma_celltype.input.covariates_file"
            ),
            "average_property": "magma_celltype.average_property",
            "cell_type_correction": "multiple_testing.methods",
            "primary_cell_type_correction": (
                "multiple_testing.primary_method"
            ),
            "cell_type_significance_threshold": (
                "multiple_testing.significance_threshold"
            ),
            "scdrs_h5ad_file": "scdrs.input.h5ad_file",
            "scdrs_gene_set_file": "scdrs.input.gene_set_file",
            "scdrs_magma_gene_results_file": (
                "scdrs.input.magma_gene_results_file"
            ),
            "scdrs_gene_id_map": "scdrs.input.gene_identifier_map_file",
            "scdrs_covariate_file": "scdrs.input.covariate_file",
            "scdrs_gene_set_source": "scdrs.magma_gene_set.source",
            "scdrs_source_gene_id_type": (
                "scdrs.magma_gene_set.source_identifier_type"
            ),
            "scdrs_target_gene_id_type": (
                "scdrs.magma_gene_set.target_identifier_type"
            ),
            "scdrs_h5ad_species": "scdrs.h5ad_species",
            "scdrs_gene_set_species": "scdrs.gene_set_species",
            "scdrs_matrix_state": "scdrs.matrix_state",
            "scdrs_group_analysis": "scdrs.downstream.group_analysis",
            "scdrs_correlation_analysis": (
                "scdrs.downstream.correlation_analysis"
            ),
            "scdrs_gene_analysis": "scdrs.downstream.gene_analysis",
            "scdrs_control_gene_sets": "scdrs.control_gene_sets",
            "ldsc_celltype_sumstats_file": "ldsc_celltype.input.sumstats_file",
            "ldsc_celltype_ldcts_file": "ldsc_celltype.input.ldcts_file",
            "ldsc_celltype_baseline_prefix": (
                "ldsc_celltype.input.baseline_ld_prefixes"
            ),
            "ldsc_celltype_weights_prefix": (
                "ldsc_celltype.input.weights_ld_prefix"
            ),
            "ldsc_celltype_merge_alleles_file": (
                "ldsc_celltype.input.merge_alleles_file"
            ),
            "ldsc_celltype_sumstats_source": "ldsc_celltype.input.source",
            "ldsc_celltype_genome_build": "ldsc_celltype.genome_build",
            "ldsc_celltype_population": "ldsc_celltype.population",
        },
    )
    global_overrides = explicit_overrides(
        args,
        {
            "dataset_id": "run.dataset_id",
            "output_directory": "run.output_directory",
            "threads": "execution.threads",
            "memory_gb": "execution.memory_gb",
            "seed": "execution.random_seed",
            "magma": "resources.executables.magma",
            "scdrs": "resources.executables.scdrs",
            "ldsc": "resources.executables.ldsc",
            "munge_sumstats": "resources.executables.munge_sumstats",
            "resume": "run.resume",
            "overwrite": "run.overwrite",
        },
    )
    return load_run_configuration_for_module(
        "single_cell",
        getattr(args, "run_config", None),
        module_overrides=module_overrides,
        global_overrides=global_overrides,
    )


def _configured_paths(output: Path, dataset: str, module) -> dict[str, Path]:
    path_names = (
        "engine_directory",
        "results_file",
        "log_file",
        "resolved_config_file",
        "completion_manifest",
        "staging_directory",
        "scdrs_engine_directory",
        "scdrs_qc_report",
        "scdrs_completion_manifest",
        "scdrs_staging_directory",
        "scdrs_generated_gene_statistics",
        "scdrs_generated_gene_set",
        "scdrs_gene_mapping_report",
        "ldsc_celltype_engine_directory",
        "ldsc_celltype_results_file",
        "ldsc_celltype_qc_report",
        "ldsc_celltype_completion_manifest",
        "ldsc_celltype_staging_directory",
        "ldsc_celltype_native_output_prefix",
        "ldsc_celltype_normalized_ldcts",
        "ldsc_celltype_munge_output_prefix",
    )
    paths = {
        name: configured_output_path(
            output,
            getattr(module.output_layout, name),
            error_type=SingleCellError,
            dataset_id=dataset,
        )
        for name in path_names
    }
    file_names = (
        "results_file",
        "log_file",
        "resolved_config_file",
        "completion_manifest",
        "scdrs_qc_report",
        "scdrs_completion_manifest",
        "scdrs_generated_gene_statistics",
        "scdrs_generated_gene_set",
        "scdrs_gene_mapping_report",
        "ldsc_celltype_results_file",
        "ldsc_celltype_qc_report",
        "ldsc_celltype_completion_manifest",
        "ldsc_celltype_normalized_ldcts",
    )
    files = {paths[name] for name in file_names}
    if len(files) != len(file_names):
        raise SingleCellError(
            "Single-cell result, QC, log, resolved configuration, and "
            "completion manifest paths must be distinct"
        )
    directory_names = (
        "staging_directory",
        "engine_directory",
        "scdrs_staging_directory",
        "scdrs_engine_directory",
        "ldsc_celltype_staging_directory",
        "ldsc_celltype_engine_directory",
    )
    directories = {paths[name] for name in directory_names}
    if len(directories) != len(directory_names) or output in directories:
        raise SingleCellError(
            "Single-cell engine and staging directories must be distinct "
            "subdirectories of the configured output directory"
        )
    prefixes = {
        paths["ldsc_celltype_native_output_prefix"],
        paths["ldsc_celltype_munge_output_prefix"],
    }
    if len(prefixes) != 2 or prefixes.intersection(files | directories):
        raise SingleCellError(
            "LDSC cell-type output prefixes must be distinct from configured "
            "files and directories"
        )
    for staging_name in (
        "staging_directory",
        "scdrs_staging_directory",
        "ldsc_celltype_staging_directory",
    ):
        staging = paths[staging_name]
        for name, path in paths.items():
            if name == staging_name:
                continue
            if path == staging or staging in path.parents:
                raise SingleCellError(
                    "modules.single_cell.output_layout.%s must not be inside "
                    "%s" % (name, staging_name)
                )
    return paths


def _fallback_log(args: argparse.Namespace, exc: BaseException) -> None:
    fallback = load_configuration()
    output = Path(
        getattr(args, "output_directory", None) or fallback.run.output_directory
    ).expanduser().resolve()
    raw_dataset = getattr(args, "dataset_id", None) or fallback.run.dataset_id
    try:
        dataset = validate_filename_component(
            raw_dataset,
            "dataset_id",
            error_type=SingleCellError,
        )
    except SingleCellError:
        dataset = fallback.run.dataset_id
    log_path = configured_output_path(
        output,
        fallback.modules.single_cell.output_layout.log_file,
        error_type=SingleCellError,
        dataset_id=dataset,
    )
    write_log_record(
        log_path,
        "ERROR",
        "Single-cell configuration failed: %s: %s"
        % (type(exc).__name__, exc),
        sample_id=dataset,
        file_level=fallback.logging.file_level,
        screen_level=fallback.logging.console_level,
    )


def preflight_single_cell_pipeline(args: argparse.Namespace) -> None:
    """Validate each selected method before upstream GWAS work starts."""
    try:
        configuration = resolve_single_cell_configuration(args)
        methods = get_single_cell_methods(
            configuration.modules.single_cell.tools
        )
        for method in methods:
            method.prepare_pipeline_args(args)
        configuration = resolve_single_cell_configuration(args)
        methods = get_single_cell_methods(
            configuration.modules.single_cell.tools
        )
        for method in methods:
            method.preflight_pipeline(args, configuration)
    except BaseException as exc:
        _fallback_log(args, exc)
        raise


def _merge_method_results(
    selected_methods: list[str],
    results: list[ModuleResult],
) -> ModuleResult:
    artifacts = {}
    metrics = {}
    warnings = []
    resumed = []
    for method, result in zip(selected_methods, results):
        overlap = set(artifacts).intersection(result.artifacts)
        if overlap:
            raise SingleCellError(
                "Selected single-cell methods produced colliding artifact keys: %s"
                % ", ".join(sorted(overlap))
            )
        artifacts.update(result.artifacts)
        metrics[method] = dict(result.metrics)
        resumed.append(bool(result.metrics.get("resumed", False)))
        warnings.extend(result.warnings)
    if len(results) == 1:
        combined_metrics = dict(results[0].metrics)
    else:
        combined_metrics = {
            "tools": metrics,
            "resumed": bool(resumed and all(resumed)),
        }
    return ModuleResult(
        "single_cell",
        artifacts=artifacts,
        metrics=combined_metrics,
        warnings=tuple(warnings),
    )


def _publish_context(ctx, result: ModuleResult) -> None:
    if ctx is None:
        return
    if hasattr(ctx, "publish"):
        ctx.publish(result)
    else:
        ctx["single_cell"] = result


def run_single_cell_direct(args: argparse.Namespace, ctx=None) -> ModuleResult:
    """Run selected methods from their exact direct or pipeline-supplied inputs."""
    try:
        configuration = resolve_single_cell_configuration(args)
        module = configuration.modules.single_cell
        output = Path(configuration.run.output_directory).expanduser().resolve()
        dataset = validate_filename_component(
            configuration.run.dataset_id,
            "dataset_id",
            error_type=SingleCellError,
        )
        paths = _configured_paths(output, dataset, module)
        methods = get_single_cell_methods(module.tools)
        preflights = {
            method.name: method.preflight_direct(args, configuration, dataset)
            for method in methods
        }
    except BaseException as exc:
        _fallback_log(args, exc)
        raise

    logger = PipelineLogger(
        dataset,
        "run",
        str(paths["log_file"].parent),
        level=configuration.logging.file_level,
        screen_level=configuration.logging.console_level,
        log_path=str(paths["log_file"]),
    )
    try:
        output.mkdir(parents=True, exist_ok=True)
        logger.record(
            "PARAM",
            "single_cell_run",
            dataset_id=dataset,
            tools=module.tools,
            overwrite=configuration.run.overwrite,
            resume=configuration.run.resume,
        )
        resource_paths = tuple(dict.fromkeys(
            resource_path
            for method in methods
            for resource_path in method.resource_paths
        ))
        resolved_modules = [
            "single_cell",
            *single_cell_supporting_configurations(module.tools),
        ]
        write_resolved_configuration(
            configuration,
            paths["resolved_config_file"],
            modules=tuple(resolved_modules),
            resource_paths=resource_paths,
        )
        context = MethodRunContext(
            args=args,
            configuration=configuration,
            output=output,
            dataset=dataset,
            paths=paths,
            logger=logger,
        )
        results = []
        for method in methods:
            logger.record("STATUS", method.name, status="STARTED")
            result = method.run(context, preflights[method.name])
            results.append(result)
            logger.record(
                "STATUS",
                method.name,
                status="COMPLETED",
                resumed=bool(result.metrics.get("resumed", False)),
            )
        combined = _merge_method_results(module.tools, results)
        logger.record(
            "STATUS",
            "single_cell_run",
            status="COMPLETED",
            tools=module.tools,
        )
        _publish_context(ctx, combined)
        return combined
    except BaseException as exc:
        logger.error(
            "Single-cell analysis failed: %s: %s"
            % (type(exc).__name__, exc)
        )
        raise
    finally:
        logger.close()


__all__ = [
    "preflight_single_cell_pipeline",
    "resolve_single_cell_configuration",
    "run_single_cell_direct",
]
