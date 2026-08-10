"""Configuration, execution, and publication boundary for single-cell analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

from postgwas.config import (
    load_configuration,
    load_run_configuration_for_module,
    write_resolved_configuration,
)
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.core.completion import (
    configuration_digest,
    validate_completion_manifest,
    write_completion_manifest,
)
from postgwas.core.contracts import Artifact, ModuleResult
from postgwas.core.paths import (
    configured_output_path,
    remove_empty_directories,
    require_nonempty_file,
    resolve_executable,
    validate_filename_component,
)
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.modules.magmacovar.service import run_magma_covar_direct
from postgwas.modules.single_cell.analysis import (
    normalize_magma_celltype_results,
    validate_magma_celltype_covariates,
    validate_magma_celltype_use_case,
)
from postgwas.modules.single_cell.errors import SingleCellError
from postgwas.modules.single_cell.scdrs import (
    ScdrsExecution,
    preflight_scdrs,
    run_scdrs,
)


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
    paths = {
        name: configured_output_path(
            output,
            getattr(module.output_layout, name),
            error_type=SingleCellError,
            dataset_id=dataset,
        )
        for name in (
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
        )
    }
    files = {
        paths["results_file"],
        paths["log_file"],
        paths["resolved_config_file"],
        paths["completion_manifest"],
        paths["scdrs_qc_report"],
        paths["scdrs_completion_manifest"],
        paths["scdrs_generated_gene_statistics"],
        paths["scdrs_generated_gene_set"],
        paths["scdrs_gene_mapping_report"],
    }
    if len(files) != 9:
        raise SingleCellError(
            "Single-cell result, QC, log, resolved configuration, and "
            "completion manifest paths must be distinct"
        )
    directories = {
        paths["staging_directory"],
        paths["engine_directory"],
        paths["scdrs_staging_directory"],
        paths["scdrs_engine_directory"],
    }
    if len(directories) != 4 or output in directories:
        raise SingleCellError(
            "Single-cell engine and staging directories must be distinct "
            "subdirectories of the configured output directory"
        )
    for staging_name in ("staging_directory", "scdrs_staging_directory"):
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
            raw_dataset, "dataset_id", error_type=SingleCellError,
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


def _magma_celltype_engine_args(
    args: argparse.Namespace,
    *,
    gene_results: Path,
    covariates: Path,
    engine_directory: Path,
) -> argparse.Namespace:
    delegated = argparse.Namespace(**vars(args))
    delegated.magma_gene_results_file = str(gene_results)
    delegated.covariates = str(covariates)
    delegated.output_directory = str(engine_directory)
    return delegated


def _magma_celltype_engine_configuration(
    configuration,
    *,
    gene_results: Path,
    covariates: Path,
    engine_directory: Path,
    use_case,
):
    """Create the exact resolved MAGMAcovar run without reparsing configuration."""
    engine_configuration = configuration.model_copy(deep=True)
    engine_configuration.run.output_directory = engine_directory
    engine = engine_configuration.modules.magmacovar
    engine.model = list(use_case.model)
    engine.direction = use_case.direction
    engine.input.gene_results_file = gene_results
    engine.input.covariates_file = covariates
    return engine_configuration


def preflight_single_cell_pipeline(args: argparse.Namespace) -> None:
    """Validate cell-type resources and tool bindings before GWAS work starts."""
    try:
        configuration = resolve_single_cell_configuration(args)
        module = configuration.modules.single_cell
        if "magma_celltype" in module.tools:
            method = module.magma_celltype
            covariates = require_nonempty_file(
                method.input.covariates_file,
                "MAGMA cell-type covariate matrix",
                error_type=SingleCellError,
            )
            validate_magma_celltype_use_case(
                module, configuration.modules.magmacovar,
            )
            validate_magma_celltype_covariates(
                covariates,
                average_property=method.average_property,
                magmacovar_config=configuration.modules.magmacovar,
            )
            resolve_executable(
                configuration.resources.executables.magma,
                "MAGMA executable",
                error_type=SingleCellError,
            )
        if "scdrs" in module.tools:
            args.scdrs_gene_set_source = "magma"
            configuration = resolve_single_cell_configuration(args)
            module = configuration.modules.single_cell
            preflight_scdrs(
                module.scdrs,
                configuration.resources.executables.scdrs,
                dataset_id=configuration.run.dataset_id,
                pipeline_pending_magma=True,
            )
    except BaseException as exc:
        _fallback_log(args, exc)
        raise


def _completion_configuration(configuration, engine_configuration) -> dict:
    return {
        "single_cell": configuration.modules.single_cell.model_dump(mode="json"),
        "effective_magmacovar": (
            engine_configuration.modules.magmacovar.model_dump(mode="json")
        ),
    }


def _magma_celltype_result(
    *,
    configuration,
    dataset: str,
    gene_property_results: Path,
    normalized_results: Path,
    metrics: dict,
    resumed: bool,
) -> ModuleResult:
    magma = configuration.modules.magma
    primary = magma.mapping.definitions[magma.mapping.primary]
    metadata = {
        "dataset_id": dataset,
        "genome_build": magma.genome_build.value,
        "population": magma.population.value,
        "gene_id_type": primary.gene_id_type,
        "method": "magma_celltype",
        "workflow": configuration.modules.single_cell.magma_celltype.workflow,
        "schema_version": (
            configuration.modules.single_cell.result_schema.normalized_schema_version
        ),
    }
    return ModuleResult(
        "single_cell",
        artifacts={
            "magma_celltype_results": Artifact(
                "single_cell_cell_type_associations",
                normalized_results,
                metadata,
            ),
            "magma_gene_property_results": Artifact(
                "magma_gene_property_results",
                gene_property_results,
                metadata,
            ),
        },
        metrics={**metrics, "resumed": resumed},
    )


def _scdrs_result(
    *,
    configuration,
    dataset: str,
    execution: ScdrsExecution,
) -> ModuleResult:
    method = configuration.modules.single_cell.scdrs
    metadata = {
        "dataset_id": dataset,
        "method": "scdrs",
        "h5ad_species": method.h5ad_species,
        "gene_set_species": method.gene_set_species,
        "matrix_source": "X",
        "matrix_state": method.matrix_state,
        "trait_ids": list(execution.traits),
    }
    artifacts = {
        "scdrs_%s" % name: Artifact(
            "scdrs_native_%s" % name,
            path,
            metadata,
        )
        for name, path in execution.outputs.items()
    }
    return ModuleResult(
        "single_cell",
        artifacts=artifacts,
        metrics={**execution.metrics, "resumed": execution.resumed},
    )


def _merge_tool_results(
    selected_tools: list[str], results: list[ModuleResult],
) -> ModuleResult:
    artifacts = {}
    metrics = {}
    warnings = []
    resumed = []
    for tool, result in zip(selected_tools, results):
        overlap = set(artifacts).intersection(result.artifacts)
        if overlap:
            raise SingleCellError(
                "Selected single-cell tools produced colliding artifact keys: %s"
                % ", ".join(sorted(overlap))
            )
        artifacts.update(result.artifacts)
        metrics[tool] = dict(result.metrics)
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


def _run_magma_celltype(
    args: argparse.Namespace,
    *,
    configuration,
    output: Path,
    dataset: str,
    paths: dict[str, Path],
    gene_results: Path,
    covariates: Path,
    use_case,
    logger,
) -> ModuleResult:
    """Delegate the published base model to the existing MAGMAcovar module."""
    module = configuration.modules.single_cell
    engine_args = _magma_celltype_engine_args(
        args,
        gene_results=gene_results,
        covariates=covariates,
        engine_directory=paths["engine_directory"],
    )
    engine_configuration = _magma_celltype_engine_configuration(
        configuration,
        gene_results=gene_results,
        covariates=covariates,
        engine_directory=paths["engine_directory"],
        use_case=use_case,
    )
    gene_property_results = Path(
        run_magma_covar_direct(
            engine_args, configuration=engine_configuration,
        )
    ).expanduser().resolve()
    digest = configuration_digest(
        _completion_configuration(configuration, engine_configuration)
    )
    completion_inputs = {
        "magma_gene_results": gene_results,
        "cell_type_covariates": covariates,
        "magma_gene_property_results": gene_property_results,
    }
    completion_outputs = {"magma_celltype_results": paths["results_file"]}
    if (
        configuration.run.resume
        and not configuration.run.overwrite
        and paths["completion_manifest"].is_file()
    ):
        manifest = validate_completion_manifest(
            paths["completion_manifest"],
            dataset_id=dataset,
            module="single_cell",
            genome_build=configuration.modules.magma.genome_build.value,
            configuration_sha256=digest,
            inputs=completion_inputs,
            outputs=completion_outputs,
            error_type=SingleCellError,
        )
        logger.record(
            "SKIP", "magma_celltype",
            reason="provenance_validated_complete_outputs",
        )
        return _magma_celltype_result(
            configuration=configuration,
            dataset=dataset,
            gene_property_results=gene_property_results,
            normalized_results=paths["results_file"],
            metrics=dict(manifest.get("metrics", {})),
            resumed=True,
        )
    if (
        paths["results_file"].exists()
        or paths["completion_manifest"].exists()
    ) and not configuration.run.overwrite:
        raise SingleCellError(
            "Existing or incomplete MAGMA cell-type output was found; use "
            "--resume for a matching complete run or --overwrite to replace it"
        )

    staging = paths["staging_directory"]
    staged_result = staging / paths["results_file"].relative_to(output)
    if staging.exists():
        owned = {staged_result.resolve()}
        unexpected = [
            path for path in staging.rglob("*")
            if path.is_file() and path.resolve() not in owned
        ]
        if unexpected:
            raise SingleCellError(
                "The configured MAGMA cell-type staging directory contains "
                "unowned files: %s"
                % ", ".join(str(path) for path in unexpected[:5])
            )
        if not configuration.run.overwrite:
            raise SingleCellError(
                "An isolated incomplete MAGMA cell-type run exists at %s; "
                "review it or use --overwrite" % staging
            )
        staged_result.unlink(missing_ok=True)
        remove_empty_directories(staged_result.parent, staging)

    metrics = normalize_magma_celltype_results(
        gene_property_results,
        covariates,
        staged_result,
        dataset_id=dataset,
        single_cell_config=module,
        magmacovar_config=configuration.modules.magmacovar,
    )
    paths["results_file"].parent.mkdir(parents=True, exist_ok=True)
    if configuration.run.overwrite:
        paths["results_file"].unlink(missing_ok=True)
        paths["completion_manifest"].unlink(missing_ok=True)
    staged_result.replace(paths["results_file"])
    remove_empty_directories(staged_result.parent, staging)
    metrics["output"] = str(paths["results_file"])
    write_completion_manifest(
        paths["completion_manifest"],
        dataset_id=dataset,
        module="single_cell",
        genome_build=configuration.modules.magma.genome_build.value,
        configuration_sha256=digest,
        inputs=completion_inputs,
        outputs=completion_outputs,
        metrics=metrics,
        error_type=SingleCellError,
    )
    logger.record(
        "OUTPUT", "magma_celltype_results",
        path=str(paths["results_file"]),
        rows=metrics["tested_cell_types"],
    )
    return _magma_celltype_result(
        configuration=configuration,
        dataset=dataset,
        gene_property_results=gene_property_results,
        normalized_results=paths["results_file"],
        metrics=metrics,
        resumed=False,
    )


def run_single_cell_direct(args: argparse.Namespace, ctx=None) -> ModuleResult:
    """Run each configured single-cell engine from its exact scientific inputs."""
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
        magma_preflight = None
        if "magma_celltype" in module.tools:
            method = module.magma_celltype
            gene_results = require_nonempty_file(
                method.input.gene_results_file,
                "MAGMA .genes.raw file",
                error_type=SingleCellError,
            )
            covariates = require_nonempty_file(
                method.input.covariates_file,
                "MAGMA cell-type covariate matrix",
                error_type=SingleCellError,
            )
            use_case = validate_magma_celltype_use_case(
                module, configuration.modules.magmacovar,
            )
            validate_magma_celltype_covariates(
                covariates,
                average_property=method.average_property,
                magmacovar_config=configuration.modules.magmacovar,
            )
            resolve_executable(
                configuration.resources.executables.magma,
                "MAGMA executable",
                error_type=SingleCellError,
            )
            method.input.gene_results_file = gene_results
            method.input.covariates_file = covariates
            magma_preflight = (gene_results, covariates, use_case)
        scdrs_preflight = None
        if "scdrs" in module.tools:
            scdrs_preflight = preflight_scdrs(
                module.scdrs,
                configuration.resources.executables.scdrs,
                dataset_id=dataset,
            )
            module.scdrs.input.h5ad_file = scdrs_preflight.h5ad_file
            module.scdrs.input.gene_set_file = scdrs_preflight.gene_set_file
            module.scdrs.input.magma_gene_results_file = (
                scdrs_preflight.magma_gene_results_file
            )
            module.scdrs.input.gene_identifier_map_file = (
                scdrs_preflight.gene_identifier_map_file
            )
            module.scdrs.input.covariate_file = scdrs_preflight.covariate_file
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
        resource_paths = tuple(
            "executables.%s" % name
            for name in ("magma", "scdrs")
            if (
                (name == "magma" and "magma_celltype" in module.tools)
                or name in module.tools
            )
        )
        write_resolved_configuration(
            configuration,
            paths["resolved_config_file"],
            modules=("single_cell",),
            resource_paths=resource_paths,
        )
        results = []
        for tool in module.tools:
            logger.record("STATUS", tool, status="STARTED")
            if tool == "magma_celltype":
                gene_results, covariates, use_case = magma_preflight
                result = _run_magma_celltype(
                    args,
                    configuration=configuration,
                    output=output,
                    dataset=dataset,
                    paths=paths,
                    gene_results=gene_results,
                    covariates=covariates,
                    use_case=use_case,
                    logger=logger,
                )
            elif tool == "scdrs":
                execution = run_scdrs(
                    preflight=scdrs_preflight,
                    configuration=configuration,
                    paths=paths,
                    dataset=dataset,
                    logger=logger,
                )
                result = _scdrs_result(
                    configuration=configuration,
                    dataset=dataset,
                    execution=execution,
                )
                logger.record(
                    "OUTPUT", "scdrs_native_outputs",
                    directory=str(paths["scdrs_engine_directory"]),
                    files=len(execution.outputs),
                    traits=len(execution.traits),
                    software_version=scdrs_preflight.software_version,
                )
            else:
                raise SingleCellError(
                    "Unsupported configured single-cell tool: %s" % tool
                )
            results.append(result)
            logger.record(
                "STATUS", tool, status="COMPLETED",
                resumed=bool(result.metrics.get("resumed", False)),
            )
        combined = _merge_tool_results(module.tools, results)
        logger.record(
            "STATUS", "single_cell_run", status="COMPLETED",
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
