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
    apply_completion_restart,
    configuration_digest,
    resolve_completion_resume,
    write_completion_manifest,
)
from postgwas.core.gene_property_validation import validate_magma_covariate_table
from postgwas.core.paths import (
    configured_output_path,
    remove_empty_directories,
    require_nonempty_file,
    resolve_executable,
    validate_filename_component,
)
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.core.preflight import (
    PipelinePreflightEvidence,
    pipeline_preflight_evidence,
    require_pipeline_input_vcf,
)
from postgwas.core.ui import (
    StageProgress,
    print_screen_block,
    screen_field,
    screen_line,
)
from postgwas.modules.magmacovar.errors import MagmaCovarError
from postgwas.modules.magmacovar.contract import (
    MAGMACOVAR_EXCLUDED_PATHWAY_OPTIONS,
)
from postgwas.modules.magmacovar.main import (
    run_magma_covariates,
    validate_corrected_magma_covariate_output,
)
from postgwas.modules.magmacovar.stages import (
    MAGMACOVAR_STAGES,
    magmacovar_pipeline_stage_numbers,
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
            "covariate_missing_values": "input.missing_values",
            "covariate_max_miss": "input.maximum_missing_fraction",
            "covariate_missing_genes": "input.missing_genes",
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
            "corrected_results_file",
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
            "results_file", "corrected_results_file", "native_log_file",
            "service_log_file",
            "resolved_config_file", "completion_manifest",
        )
    }
    if len(distinct_files) != 6:
        raise MagmaCovarError(
            "MAGMAcovar raw result, corrected result, native log, service log, "
            "resolved configuration, and completion manifest paths must be distinct"
        )
    staging = values["staging_directory"]
    if staging == output or staging in values["output_prefix"].parents:
        raise MagmaCovarError(
            "modules.magmacovar.output_layout.staging_directory must be an "
            "isolated directory that does not contain the final output prefix"
        )
    for name in (
        "corrected_results_file", "service_log_file", "resolved_config_file",
        "completion_manifest",
    ):
        if values[name] == staging or staging in values[name].parents:
            raise MagmaCovarError(
                "modules.magmacovar.output_layout.%s must not be inside the "
                "staging directory" % name
            )
    return values


def _configured_artifact_paths(paths: dict[str, Path]) -> dict[str, Path]:
    return {
        "raw_results": paths["results_file"],
        "corrected_results": paths["corrected_results_file"],
        "native_log": paths["native_log_file"],
    }


def _staged_artifact_paths(
    paths: dict[str, Path], output: Path,
) -> tuple[Path, dict[str, Path]]:
    """Map each configured artifact to the same relative path in staging."""
    staging = paths["staging_directory"]
    staged_prefix = staging / paths["output_prefix"].relative_to(output)
    return staged_prefix, {
        "raw_results": staging / paths["results_file"].relative_to(output),
        "corrected_results": (
            staging / paths["corrected_results_file"].relative_to(output)
        ),
        "native_log": staging / paths["native_log_file"].relative_to(output),
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
            "raw-result, corrected-result, and native-log outputs"
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


def _result_payload(paths: dict[str, Path], summary: dict) -> dict:
    return {
        "raw_results": str(paths["results_file"]),
        "corrected_results": str(paths["corrected_results_file"]),
        "native_log": str(paths["native_log_file"]),
        **summary,
    }


def _publish_validated_outputs(
    *,
    staging: Path,
    staged_artifacts: dict[str, Path],
    configured_artifacts: dict[str, Path],
    paths: dict[str, Path],
    dataset_id: str,
    overwrite: bool,
    completion_digest: str,
    completion_inputs: dict[str, Path],
    input_metrics: dict,
    output_metrics: dict,
) -> dict[str, Path]:
    """Publish validated outputs and atomically record their completion proof."""
    unexpected_staged_paths = [
        path for path in staging.rglob("*")
        if path.is_file() and path not in staged_artifacts.values()
    ]
    if unexpected_staged_paths:
        raise MagmaCovarError(
            "MAGMA produced unconfigured staged artifacts; refusing to "
            "publish: %s"
            % ", ".join(str(path) for path in unexpected_staged_paths)
        )

    if overwrite:
        if paths["completion_manifest"].is_file():
            try:
                existing_artifacts = _recorded_artifact_paths(
                    paths["completion_manifest"], configured_artifacts,
                )
            except MagmaCovarError:
                existing_artifacts = _required_artifact_paths(paths)
        else:
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
    for parent in {path.parent for path in staged_artifacts.values()}:
        remove_empty_directories(parent, staging, staging.parent)

    completion_metrics = {
        key: value
        for key, value in input_metrics.items()
        if key not in {"property_names", "property_missingness"}
    }
    completion_metrics.update(output_metrics)
    write_completion_manifest(
        paths["completion_manifest"],
        dataset_id=dataset_id,
        module="magmacovar",
        genome_build="not_applicable",
        configuration_sha256=completion_digest,
        inputs=completion_inputs,
        outputs=published,
        metrics=completion_metrics,
        error_type=MagmaCovarError,
    )
    return published


def _format_finding_number(value: float, significant_digits: int) -> str:
    return format(float(value), ".%dg" % significant_digits)


def _flames_handoff_requested(args: argparse.Namespace) -> bool:
    requested_modules = tuple(
        getattr(
            args,
            "_pipeline_requested_modules",
            getattr(args, "modules", ()),
        )
        or ()
    )
    return "flames" in requested_modules


def _render_main_findings(
    dataset: str,
    paths: dict[str, Path],
    summary: dict,
    module,
    *,
    include_flames_handoff: bool,
) -> str:
    primary = summary["primary_correction_method"]
    primary_label = module.multiple_testing.reporting_method_labels[primary]
    method_labels = ", ".join(
        module.multiple_testing.reporting_method_labels[method]
        for method in summary["correction_methods"]
    )
    significant_by_method = summary["significant_properties_by_method"]
    top_properties = summary["top_properties"]
    highlight_method = module.reporting.highlight_method
    highlight_label = module.multiple_testing.reporting_method_labels[
        highlight_method
    ]
    p_digits = module.reporting.p_value_significant_digits
    effect_digits = module.reporting.effect_significant_digits
    label_width = 38
    lines = [
        "",
        screen_line("analysis", "MAGMAcovar main findings", indent=2),
        screen_field(
            "info", "Dataset", dataset, indent=6, label_width=label_width,
        ),
        screen_field(
            "count", "Gene properties tested", summary["tested_properties"],
            indent=6, label_width=label_width,
        ),
        screen_field(
            "genetic", "Genes per tested property",
            "%s–%s"
            % (
                f"{summary['minimum_result_genes']:,}",
                f"{summary['maximum_result_genes']:,}",
            ),
            indent=6, label_width=label_width,
        ),
        screen_field(
            "decision", "Gene-property model",
            ", ".join(module.model) if module.model else "marginal",
            indent=6, label_width=label_width,
        ),
        screen_field(
            "decision", "Property-test direction", module.direction,
            indent=6, label_width=label_width,
        ),
        screen_field(
            "decision", "Primary multiple-testing method", primary_label,
            indent=6, label_width=label_width,
        ),
        screen_field(
            "info", "Corrections reported", method_labels,
            indent=6, label_width=label_width,
        ),
        screen_field(
            "info", "Multiple-testing scope",
            "all %d COVAR properties in this run" % summary["tested_properties"],
            indent=6, label_width=label_width,
        ),
        screen_field(
            "info", "Primary significance threshold",
            summary["significance_threshold"], indent=6, label_width=label_width,
        ),
        screen_field(
            "success" if summary["primary_significant_properties"] else "info",
            "%s-significant properties" % primary_label,
            "%d / %d"
            % (
                summary["primary_significant_properties"],
                summary["tested_properties"],
            ),
            indent=6,
            label_width=label_width,
        ),
    ]
    highlight_significant = significant_by_method.get(highlight_method)
    if highlight_significant is not None and highlight_method != primary:
        lines.append(screen_field(
            "success" if highlight_significant else "info",
            "%s-significant properties" % highlight_label,
            "%d / %d" % (
                highlight_significant, summary["tested_properties"],
            ),
            indent=6,
            label_width=label_width,
        ))
    lines.append(screen_field(
        "decision", "Top-property ranking",
        "lowest raw MAGMA P; showing %d" % len(top_properties),
        indent=6,
        label_width=label_width,
    ))
    for index, finding in enumerate(top_properties, 1):
        statistics = [
            "BETA_STD=%s"
            % _format_finding_number(
                finding["standardized_beta"], effect_digits,
            ),
            "P=%s" % _format_finding_number(finding["p_value"], p_digits),
        ]
        highlighted_adjusted = finding["adjusted_p_values"].get(
            highlight_method
        )
        if highlighted_adjusted is not None:
            statistics.append(
                "%s=%s"
                % (
                    highlight_label,
                    _format_finding_number(highlighted_adjusted, p_digits),
                )
            )
        lines.append(screen_field(
            "analysis",
            "Top property %d" % index,
            "%s (%s)" % (finding["property"], "; ".join(statistics)),
            indent=6,
            label_width=label_width,
            break_long_values=True,
        ))
    lines.append(screen_field(
        "genetic", "Native MAGMA result", paths["results_file"],
        indent=6, label_width=label_width, break_long_values=True,
    ))
    if include_flames_handoff:
        lines.append(screen_field(
            "decision", "FLAMES handoff",
            "native .gsa.out with raw P values (not the corrected TSV)",
            indent=6, label_width=label_width,
        ))
    lines.extend((
        screen_field(
            "success", "Corrected results for interpretation",
            paths["corrected_results_file"], indent=6,
            label_width=label_width, break_long_values=True,
        ),
        "",
    ))
    return "\n".join(lines)


def _report_resumed_stages(
    args: argparse.Namespace,
    paths: dict[str, Path],
    manifest: dict,
    output_summary: dict,
    logger,
) -> None:
    """Display and log checkpoint-backed evidence for a resumed run."""
    progress = getattr(args, "_pipeline_stage_progress", None)
    stage_numbers = magmacovar_pipeline_stage_numbers(args)
    metrics = manifest.get("metrics") or {}
    evidence = {
        "gene_results": [
            ("analysis", "MAGMA gene-association results"),
            (
                "info",
                "Input file",
                Path(
                    metrics.get(
                        "gene_results_file", "validated checkpoint input",
                    )
                ).name,
            ),
            (
                "genetic",
                "Eligible gene identifiers",
                metrics.get("gene_results_genes", "recorded in checkpoint"),
            ),
            ("success", "Input fingerprint", "validated and unchanged"),
        ],
        "covariates": [
            ("analysis", "Gene-property covariate reference"),
            (
                "info",
                "Input file",
                Path(
                    metrics.get("covariates_file", "validated checkpoint input")
                ).name,
            ),
            (
                "count",
                "Covariate-table genes",
                metrics.get("covariate_genes", "recorded in checkpoint"),
            ),
            (
                "count",
                "Gene properties",
                metrics.get("properties", "recorded in checkpoint"),
            ),
            (
                "genetic",
                "Gene IDs shared with MAGMA results",
                metrics.get("overlapping_genes", "recorded in checkpoint"),
            ),
            ("success", "Input fingerprint", "validated and unchanged"),
        ],
        "analysis": [
            ("analysis", "MAGMA gene-property model"),
            (
                "success",
                "Native execution",
                "reused from validated checkpoint",
            ),
        ],
        "results": [
            ("analysis", "Gene-property results"),
            (
                "count",
                "Tested properties",
                output_summary["tested_properties"],
            ),
            (
                "success",
                "Native and corrected results",
                "validated and unchanged",
            ),
        ],
    }
    if progress is None or stage_numbers is None:
        for number, key in enumerate(
            ("gene_results", "covariates", "analysis", "results"),
            1,
        ):
            with logger.step(
                number,
                len(MAGMACOVAR_STAGES),
                MAGMACOVAR_STAGES[number - 1],
                "resume_validated_magmacovar_stage",
            ) as step:
                step.outcome(
                    "Reused checksum-validated MAGMAcovar evidence.",
                    fields=evidence[key],
                    resume_mode="validated_checkpoint",
                )
        return
    for key in ("gene_results", "covariates", "analysis"):
        number = int(stage_numbers[key])
        progress.start(number)
        progress.complete(number, outcome_fields=evidence[key])
    progress.start(int(stage_numbers["results"]))
    args._pipeline_stage_completion = {
        "outcome_fields": evidence["results"],
    }


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


def preflight_magmacovar_pipeline(
    args: argparse.Namespace,
    *,
    preflight_evidence=None,
) -> PipelinePreflightEvidence:
    """Validate external MAGMAcovar resources before upstream pipeline work."""
    require_pipeline_input_vcf(preflight_evidence)
    try:
        requested_modules = tuple(
            getattr(
                args,
                "_pipeline_requested_modules",
                getattr(args, "modules", ()),
            )
            or ()
        )
        excluded_arguments = [
            flag
            for destination, flag in MAGMACOVAR_EXCLUDED_PATHWAY_OPTIONS
            if getattr(args, destination, None) is not None
        ]
        if requested_modules == ("magmacovar",) and excluded_arguments:
            raise MagmaCovarError(
                "Invalid pathway argument(s) for a MAGMAcovar-only pipeline: "
                "%s. MAGMAcovar requires the upstream MAGMA "
                "gene-association result but does not run pathway analysis; "
                "use --modules magma when pathway association is required."
                % ", ".join(excluded_arguments)
            )
        configuration = resolve_magmacovar_configuration(args)
        module = configuration.modules.magmacovar
        covariates = validate_magma_covariate_table(
            module.input.covariates_file,
            minimum_genes=module.minimum_genes,
            maximum_missing_fraction=module.input.maximum_missing_fraction,
            missing_genes=module.input.missing_genes,
            error_type=MagmaCovarError,
        )
        executable = resolve_executable(
            configuration.resources.executables.magma,
            "MAGMA executable",
            error_type=MagmaCovarError,
        )
        resources = {
            "configuration": configuration,
            "covariates": covariates,
            "executable": executable,
        }
        return pipeline_preflight_evidence(
            "magmacovar",
            preflight_evidence,
            resources=resources,
            deferred_checks=(
                "Validate the pipeline-generated MAGMA gene-association result.",
            ),
        )
    except BaseException as exc:
        _fallback_log(args, exc)
        raise


def run_magma_covar_direct(
    args: argparse.Namespace, ctx=None, *, configuration=None,
) -> dict:
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
        stage_progress=StageProgress(
            "MAGMAcovar analysis progress",
            enabled=getattr(args, "_pipeline_stage_progress", None) is None,
            outcome_label_width=configuration.logging.terminal_label_width,
        ),
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
            correction_methods=module.multiple_testing.methods,
            primary_correction_method=module.multiple_testing.primary_method,
            significance_threshold=(
                module.multiple_testing.significance_threshold
            ),
            findings_highlight_method=module.reporting.highlight_method,
            top_property_count=module.reporting.top_property_count,
            p_value_significant_digits=(
                module.reporting.p_value_significant_digits
            ),
            effect_significant_digits=(
                module.reporting.effect_significant_digits
            ),
            raw_results=paths["results_file"],
            corrected_results=paths["corrected_results_file"],
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
            decision = resolve_completion_resume(
                paths["completion_manifest"],
                dataset_id=dataset,
                module="magmacovar",
                genome_build="not_applicable",
                configuration_sha256=completion_digest,
                inputs=completion_inputs,
                outputs=configured_artifacts,
                resume_policy=configuration.run.resume_policy,
                error_type=MagmaCovarError,
            )
            if decision.action == "resume":
                summary = validate_corrected_magma_covariate_output(
                    paths["results_file"],
                    paths["corrected_results_file"],
                    module=module,
                )
                logger.record(
                    "SKIP", "magmacovar_run",
                    reason="provenance_validated_complete_outputs",
                    **summary,
                )
                _report_resumed_stages(
                    args, paths, dict(decision.manifest), summary, logger,
                )
                result = _result_payload(paths, summary)
                if ctx is not None:
                    ctx["magma_covar"] = result
                print_screen_block(
                    _render_main_findings(
                        dataset,
                        paths,
                        summary,
                        module,
                        include_flames_handoff=_flames_handoff_requested(args),
                    )
                )
                return result
            apply_completion_restart(
                decision,
                output_root=output,
                manifest=paths["completion_manifest"],
                logger=logger,
                operation="magmacovar_resume",
                error_type=MagmaCovarError,
            )
            existing_artifacts = _required_artifact_paths(paths)

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
        staged_prefix, staged_artifacts = _staged_artifact_paths(paths, output)
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
            for parent in {path.parent for path in staged_artifacts.values()}:
                remove_empty_directories(parent, staging)
            if staging.exists():
                raise MagmaCovarError(
                    "The configured MAGMAcovar staging directory contains files "
                    "not owned by this output prefix; refusing overwrite: %s"
                    % staging
                )
        staged_prefix.parent.mkdir(parents=True, exist_ok=True)

        def finalize_results(input_metrics: dict, output_metrics: dict):
            return _publish_validated_outputs(
                staging=staging,
                staged_artifacts=staged_artifacts,
                configured_artifacts=configured_artifacts,
                paths=paths,
                dataset_id=dataset,
                overwrite=configuration.run.overwrite,
                completion_digest=completion_digest,
                completion_inputs=completion_inputs,
                input_metrics=input_metrics,
                output_metrics=output_metrics,
            )

        metrics = run_magma_covariates(
            magma_bin=magma,
            gene_results_file=gene_results,
            covariates_file=covariates,
            output_prefix=staged_prefix,
            results_file=staged_artifacts["raw_results"],
            corrected_results_file=staged_artifacts["corrected_results"],
            native_log_file=staged_artifacts["native_log"],
            module=module,
            logger=logger,
            timeout_seconds=configuration.execution.timeout_seconds,
            pipeline_progress=getattr(args, "_pipeline_stage_progress", None),
            pipeline_stage_numbers=magmacovar_pipeline_stage_numbers(args),
            finalize_results=finalize_results,
        )
        published = metrics["finalized"]
        pipeline_progress = getattr(args, "_pipeline_stage_progress", None)
        pipeline_stage_numbers = magmacovar_pipeline_stage_numbers(args)
        if pipeline_progress is not None and pipeline_stage_numbers is not None:
            args._pipeline_stage_completion = {
                "outcome_fields": metrics["result_fields"],
            }
        logger.record(
            "STATUS",
            "magmacovar_run",
            status="COMPLETED",
            raw_results=paths["results_file"],
            corrected_results=paths["corrected_results_file"],
            artifacts=len(published),
            **metrics["output"],
        )
        result = _result_payload(paths, metrics["output"])
        if ctx is not None:
            ctx["magma_covar"] = result
        print_screen_block(
            _render_main_findings(
                dataset,
                paths,
                metrics["output"],
                module,
                include_flames_handoff=_flames_handoff_requested(args),
            )
        )
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
