"""Configuration, provenance, and publication boundary for MAGMAcovar."""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml

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
from postgwas.core.paths import (
    configured_output_path,
    remove_empty_directories,
    require_nonempty_file,
    resolve_executable,
    validate_filename_component,
)
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.modules.magmacovar.errors import MagmaCovarError
from postgwas.modules.magmacovar.main import (
    run_magma_covariates,
    validate_magma_covariate_table,
    validate_magma_covariate_output,
)


def resolve_magmacovar_configuration(args: argparse.Namespace):
    """Resolve the effective MAGMAcovar settings, including explicit CLI values."""
    module_overrides = explicit_overrides(
        args,
        {
            "magma_gene_results_file": "input.gene_results_file",
            "covariates": "input.covariates_file",
            "covariate_model": "model",
            "covariate_direction": "direction",
            "minimum_genes": "minimum_genes",
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
            "resume": "run.resume",
            "overwrite": "run.overwrite",
        },
    )
    return load_run_configuration_for_module(
        "magmacovar",
        getattr(args, "run_config", None),
        module_overrides=module_overrides,
        global_overrides=global_overrides,
    )


def _configured_paths(output: Path, dataset: str, module) -> dict[str, Path]:
    values = {
        name: configured_output_path(
            output,
            getattr(module.output_layout, name),
            error_type=MagmaCovarError,
            dataset_id=dataset,
        )
        for name in (
            "output_prefix",
            "results_file",
            "native_log_file",
            "service_log_file",
            "resolved_config_file",
            "completion_manifest",
            "staging_directory",
        )
    }
    expected_results = Path(str(values["output_prefix"]) + ".gsa.out")
    expected_native_log = Path(str(values["output_prefix"]) + ".log")
    if values["results_file"] != expected_results:
        raise MagmaCovarError(
            "modules.magmacovar.output_layout.results_file must match MAGMA's "
            "<output_prefix>.gsa.out invariant"
        )
    if values["native_log_file"] != expected_native_log:
        raise MagmaCovarError(
            "modules.magmacovar.output_layout.native_log_file must match MAGMA's "
            "<output_prefix>.log invariant"
        )
    distinct_files = {
        values[name]
        for name in (
            "results_file", "native_log_file", "service_log_file",
            "resolved_config_file", "completion_manifest",
        )
    }
    if len(distinct_files) != 5:
        raise MagmaCovarError(
            "MAGMAcovar result, native log, service log, resolved configuration, "
            "and completion manifest paths must be distinct"
        )
    staging = values["staging_directory"]
    if staging == output or staging in values["output_prefix"].parents:
        raise MagmaCovarError(
            "modules.magmacovar.output_layout.staging_directory must be an "
            "isolated directory that does not contain the final output prefix"
        )
    for name in (
        "service_log_file", "resolved_config_file", "completion_manifest",
    ):
        if values[name] == staging or staging in values[name].parents:
            raise MagmaCovarError(
                "modules.magmacovar.output_layout.%s must not be inside the "
                "staging directory" % name
            )
    return values


def _configured_artifact_paths(paths: dict[str, Path]) -> dict[str, Path]:
    prefix = paths["output_prefix"].name + "."
    return {
        paths[name].name.removeprefix(prefix): paths[name]
        for name in ("results_file", "native_log_file")
    }


def _required_artifact_paths(paths: dict[str, Path]) -> dict[str, Path]:
    artifacts = {}
    for name, path in _configured_artifact_paths(paths).items():
        if path.exists() and not path.is_file():
            raise MagmaCovarError(
                "Configured MAGMAcovar output path is not a file: %s" % path
            )
        if path.is_file():
            artifacts[name] = path
    return artifacts


def _recorded_artifact_paths(
    completion_manifest: Path, expected_outputs: dict[str, Path],
) -> dict[str, Path]:
    """Read exact owned output paths without treating arbitrary prefix files as owned."""
    try:
        document = yaml.safe_load(
            completion_manifest.read_text(encoding="utf-8")
        ) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise MagmaCovarError(
            "Cannot read MAGMAcovar completion manifest %s: %s"
            % (completion_manifest, exc)
        ) from exc
    recorded = document.get("outputs")
    if not isinstance(recorded, dict) or not recorded:
        raise MagmaCovarError(
            "MAGMAcovar completion manifest contains no recorded outputs"
        )
    if set(recorded) != set(expected_outputs):
        raise MagmaCovarError(
            "MAGMAcovar completion manifest does not contain the exact configured "
            "result and native-log outputs"
        )
    outputs = {}
    for name, fingerprint in recorded.items():
        if not isinstance(fingerprint, dict) or not fingerprint.get("path"):
            raise MagmaCovarError(
                "MAGMAcovar completion manifest contains an invalid output fingerprint"
            )
        path = Path(fingerprint["path"]).expanduser().resolve()
        if path != expected_outputs[name]:
            raise MagmaCovarError(
                "MAGMAcovar completion manifest output does not match its configured "
                "path: %s" % path
            )
        outputs[str(name)] = path
    return outputs


def _completion_digest(module, magma_executable: str) -> str:
    return configuration_digest({
        "module": module.model_dump(mode="json"),
        "magma_executable": magma_executable,
    })


def _fallback_log(args: argparse.Namespace, exc: BaseException) -> None:
    fallback = load_configuration()
    output = Path(
        getattr(args, "output_directory", None) or fallback.run.output_directory
    ).expanduser().resolve()
    raw_dataset = getattr(args, "dataset_id", None) or fallback.run.dataset_id
    try:
        dataset = validate_filename_component(
            raw_dataset, "dataset_id", error_type=MagmaCovarError,
        )
    except MagmaCovarError:
        dataset = fallback.run.dataset_id
    log_path = configured_output_path(
        output,
        fallback.modules.magmacovar.output_layout.service_log_file,
        error_type=MagmaCovarError,
        dataset_id=dataset,
    )
    write_log_record(
        log_path,
        "ERROR",
        "MAGMAcovar configuration failed: %s: %s" % (type(exc).__name__, exc),
        sample_id=dataset,
        file_level=fallback.logging.file_level,
        screen_level=fallback.logging.console_level,
    )


def preflight_magmacovar_pipeline(args: argparse.Namespace) -> None:
    """Validate external MAGMAcovar resources before upstream pipeline work."""
    try:
        configuration = resolve_magmacovar_configuration(args)
        module = configuration.modules.magmacovar
        validate_magma_covariate_table(
            module.input.covariates_file,
            minimum_genes=module.minimum_genes,
            maximum_missing_fraction=module.input.maximum_missing_fraction,
            missing_genes=module.input.missing_genes,
        )
        resolve_executable(
            configuration.resources.executables.magma,
            "MAGMA executable",
            error_type=MagmaCovarError,
        )
    except BaseException as exc:
        _fallback_log(args, exc)
        raise


def run_magma_covar_direct(
    args: argparse.Namespace, ctx=None, *, configuration=None,
) -> str:
    """Resolve configuration once and publish only validated MAGMAcovar output."""
    try:
        configuration = configuration or resolve_magmacovar_configuration(args)
        module = configuration.modules.magmacovar
        output = Path(configuration.run.output_directory).expanduser().resolve()
        dataset = validate_filename_component(
            configuration.run.dataset_id,
            "dataset_id",
            error_type=MagmaCovarError,
        )
        paths = _configured_paths(output, dataset, module)
    except BaseException as exc:
        _fallback_log(args, exc)
        raise

    logger = PipelineLogger(
        dataset,
        "run",
        str(paths["service_log_file"].parent),
        level=configuration.logging.file_level,
        screen_level=configuration.logging.console_level,
        log_path=str(paths["service_log_file"]),
    )
    try:
        output.mkdir(parents=True, exist_ok=True)
        gene_results = require_nonempty_file(
            module.input.gene_results_file,
            "MAGMA .genes.raw file",
            error_type=MagmaCovarError,
        )
        covariates = require_nonempty_file(
            module.input.covariates_file,
            "MAGMA gene-covariate file",
            error_type=MagmaCovarError,
        )
        magma = resolve_executable(
            configuration.resources.executables.magma,
            "MAGMA executable",
            error_type=MagmaCovarError,
        )
        module.input.gene_results_file = gene_results
        module.input.covariates_file = covariates
        configuration.resources.executables.magma = magma
        completion_digest = _completion_digest(module, magma)
        completion_inputs = {
            "gene_results": gene_results,
            "covariates": covariates,
            "magma_executable": Path(magma),
        }
        logger.record(
            "PARAM",
            "magmacovar_run",
            dataset_id=dataset,
            model=module.model,
            direction=module.direction,
            minimum_genes=module.minimum_genes,
            missing_values=module.input.missing_values,
            maximum_missing_fraction=module.input.maximum_missing_fraction,
            missing_genes=module.input.missing_genes,
            magma=magma,
            timeout_seconds=configuration.execution.timeout_seconds,
            overwrite=configuration.run.overwrite,
            resume=configuration.run.resume,
        )
        write_resolved_configuration(
            configuration,
            paths["resolved_config_file"],
            modules="magmacovar",
            resource_paths=("executables.magma",),
        )

        configured_artifacts = _configured_artifact_paths(paths)
        existing_artifacts = _required_artifact_paths(paths)
        if (
            configuration.run.resume
            and not configuration.run.overwrite
            and paths["completion_manifest"].is_file()
        ):
            existing_artifacts = _recorded_artifact_paths(
                paths["completion_manifest"], configured_artifacts,
            )
            validate_completion_manifest(
                paths["completion_manifest"],
                dataset_id=dataset,
                module="magmacovar",
                genome_build="not_applicable",
                configuration_sha256=completion_digest,
                inputs=completion_inputs,
                outputs=existing_artifacts,
                error_type=MagmaCovarError,
            )
            summary = validate_magma_covariate_output(
                paths["results_file"], minimum_genes=module.minimum_genes,
            )
            logger.record(
                "SKIP", "magmacovar_run",
                reason="provenance_validated_complete_outputs",
                **summary,
            )
            result = str(paths["results_file"])
            if ctx is not None:
                ctx["magma_covar"] = result
            return result

        if (
            (existing_artifacts or paths["completion_manifest"].exists())
            and not configuration.run.overwrite
        ):
            raise MagmaCovarError(
                "Existing or incomplete MAGMAcovar output was found for prefix %s; "
                "use --resume for a matching complete run or --overwrite to replace it"
                % paths["output_prefix"]
            )

        staging = paths["staging_directory"]
        relative_prefix = paths["output_prefix"].relative_to(output)
        staged_prefix = staging / relative_prefix
        staged_artifacts = {
            name: Path(str(staged_prefix) + "." + name)
            for name in configured_artifacts
        }
        if staging.exists():
            if not configuration.run.overwrite:
                raise MagmaCovarError(
                    "An isolated incomplete MAGMAcovar run exists at %s; review it "
                    "or use --overwrite" % staging
                )
            protected_inputs = {
                gene_results.resolve(), covariates.resolve(), Path(magma).resolve(),
            }
            if any(
                artifact.resolve() in protected_inputs
                for artifact in staged_artifacts.values()
            ):
                raise MagmaCovarError(
                    "The MAGMAcovar staging prefix overlaps a configured input or "
                    "executable; refusing destructive overwrite"
                )
            for artifact in staged_artifacts.values():
                if artifact.exists() and not artifact.is_file():
                    raise MagmaCovarError(
                        "Configured MAGMAcovar staged output is not a file: %s"
                        % artifact
                    )
                artifact.unlink(missing_ok=True)
            remove_empty_directories(staged_prefix.parent, staging)
            if staging.exists():
                raise MagmaCovarError(
                    "The configured MAGMAcovar staging directory contains files "
                    "not owned by this output prefix; refusing overwrite: %s"
                    % staging
                )
        staged_prefix.parent.mkdir(parents=True, exist_ok=True)
        metrics = run_magma_covariates(
            magma_bin=magma,
            gene_results_file=gene_results,
            covariates_file=covariates,
            output_prefix=staged_prefix,
            results_file=staged_artifacts["gsa.out"],
            native_log_file=staged_artifacts["log"],
            module=module,
            logger=logger,
            timeout_seconds=configuration.execution.timeout_seconds,
        )
        unexpected_staged_paths = [
            path for path in staged_prefix.parent.iterdir()
            if path not in staged_artifacts.values()
        ]
        if unexpected_staged_paths:
            raise MagmaCovarError(
                "MAGMA produced unconfigured staged artifacts; refusing to publish: %s"
                % ", ".join(str(path) for path in unexpected_staged_paths)
            )

        if configuration.run.overwrite:
            if paths["completion_manifest"].is_file():
                try:
                    existing_artifacts = _recorded_artifact_paths(
                        paths["completion_manifest"], configured_artifacts,
                    )
                except MagmaCovarError:
                    existing_artifacts = _required_artifact_paths(paths)
            for path in existing_artifacts.values():
                path.unlink(missing_ok=True)
            paths["completion_manifest"].unlink(missing_ok=True)
        published = {}
        for name, source in staged_artifacts.items():
            destination = configured_artifacts[name]
            destination.parent.mkdir(parents=True, exist_ok=True)
            source.replace(destination)
            published[name] = destination
        remove_empty_directories(staged_prefix.parent, staging, staging.parent)

        write_completion_manifest(
            paths["completion_manifest"],
            dataset_id=dataset,
            module="magmacovar",
            genome_build="not_applicable",
            configuration_sha256=completion_digest,
            inputs=completion_inputs,
            outputs=published,
            metrics={
                "gene_results_genes": metrics["inputs"]["gene_results_genes"],
                "covariate_genes": metrics["inputs"]["covariate_genes"],
                "overlapping_genes": metrics["inputs"]["overlapping_genes"],
                **metrics["output"],
            },
            error_type=MagmaCovarError,
        )
        logger.record(
            "STATUS",
            "magmacovar_run",
            status="COMPLETED",
            output=paths["results_file"],
            artifacts=len(published),
            **metrics["output"],
        )
        result = str(paths["results_file"])
        if ctx is not None:
            ctx["magma_covar"] = result
        print("\nMAGMA gene-property analysis completed: %s\n" % result)
        return result
    except BaseException as exc:
        logger.error(
            "MAGMAcovar analysis failed: %s: %s" % (type(exc).__name__, exc)
        )
        raise
    finally:
        logger.close()


__all__ = [
    "preflight_magmacovar_pipeline",
    "resolve_magmacovar_configuration",
    "run_magma_covar_direct",
]
