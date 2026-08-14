"""Configuration boundary and publication service for LDSC heritability."""

from __future__ import annotations

import argparse
from pathlib import Path
import tempfile
from typing import Any

from postgwas.config import (
    load_configuration,
    load_run_configuration_for_module,
    write_resolved_configuration,
)
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.core.paths import configured_output_path, resolve_executable
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.modules.ldsc.ldsc_runner import (
    LDSCError,
    ldsc_owned_output_paths,
    run_ldsc,
)


MODULE_CLI_OVERRIDES = {
    "ldsc_population": "population",
    "ldsc_genome_build": "genome_build",
    "ldsc_minimum_info": "minimum_info",
    "ldsc_minimum_maf": "minimum_maf",
    "ldsc_minimum_n": "minimum_n",
    "ldsc_chunksize": "chunksize",
    "ldsc_keep_maf": "keep_maf",
    "ldsc_intercept": "intercept",
    "ldsc_two_step": "two_step",
    "ldsc_chisq_max": "chisq_max",
    "ldsc_n_blocks": "n_blocks",
    "ldsc_use_m_5_50": "use_m_5_50",
    "ldsc_print_covariance": "print_covariance",
    "ldsc_print_delete_values": "print_delete_values",
    "samp_prev": "sample_prevalence",
    "pop_prev": "population_prevalence",
}

GLOBAL_CLI_OVERRIDES = {
    "dataset_id": "run.dataset_id",
    "output_directory": "run.output_directory",
    "threads": "execution.threads",
    "memory_gb": "execution.memory_gb",
    "seed": "execution.random_seed",
    "overwrite": "run.overwrite",
    "ldsc_executable": "resources.executables.ldsc",
    "munge_sumstats_executable": "resources.executables.munge_sumstats",
}


def resolve_ldsc_configuration(args: argparse.Namespace):
    """Resolve packaged defaults, YAML, and only explicitly supplied CLI values."""
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
        "ldsc",
        getattr(args, "run_config", None),
        module_overrides=module_overrides,
        global_overrides=global_overrides,
    )


def _formatter_sample_prevalence(ctx: dict[str, Any] | None) -> float:
    """Read, but never derive, sample prevalence from the formatter result."""
    try:
        formatter_result = ctx["formatter"]["ldsc"]  # type: ignore[index]
    except (KeyError, TypeError) as exc:
        raise LDSCError(
            "Population prevalence requires sample prevalence. Supply --samp-prev "
            "for direct execution, or run the LDSC formatter first so its returned "
            "sample_prev artifact is available."
        ) from exc
    value = formatter_result.get("sample_prev")
    if value is None:
        trait_type = formatter_result.get("trait_type", "unknown")
        raise LDSCError(
            "The formatter returned no sample prevalence (trait_type=%s); "
            "liability-scale LDSC cannot run." % trait_type
        )
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise LDSCError(
            "The formatter returned a non-numeric sample_prev value: %r" % value
        ) from exc


def _resolve_prevalence(module, ctx: dict[str, Any] | None) -> str:
    """Resolve sample prevalence only when population prevalence requests liability."""
    source = (
        "configuration_or_cli"
        if module.sample_prevalence is not None
        else "unset"
    )
    if module.population_prevalence is None:
        return source
    if module.sample_prevalence is None:
        module.sample_prevalence = _formatter_sample_prevalence(ctx)
        source = "formatter_return_value"
    return source


def _publish_staged_result(
    result: dict[str, Any], staged_prefix: Path, final_prefix: Path,
) -> dict[str, Any]:
    """Move the exact validated files out of staging and update result paths."""
    staged_paths = list(dict.fromkeys([
        Path(result["munged_sumstats"]),
        Path(result["munge_log"]),
        *(Path(value) for value in result["observed_outputs"]),
        *(Path(value) for value in result["liability_outputs"]),
    ]))
    replacements: dict[str, str] = {}
    for staged in staged_paths:
        if staged.parent != staged_prefix.parent or not staged.name.startswith(
            staged_prefix.name
        ):
            raise LDSCError(
                "Refusing to publish unexpected staged LDSC path: %s" % staged
            )
        if not staged.is_file() or staged.stat().st_size <= 0:
            raise LDSCError(
                "Refusing to publish missing or empty staged LDSC output: %s"
                % staged
            )
    for staged in staged_paths:
        suffix = staged.name[len(staged_prefix.name):]
        final = Path(str(final_prefix) + suffix)
        final.parent.mkdir(parents=True, exist_ok=True)
        staged.replace(final)
        replacements[str(staged)] = str(final)

    for key in ("munged_sumstats", "munge_log", "h2_observed", "h2_liability"):
        if result[key] is not None:
            result[key] = replacements[result[key]]
    for key in ("observed_outputs", "liability_outputs"):
        result[key] = [replacements[value] for value in result[key]]
    return result


def _record_findings(logger: PipelineLogger, result: dict[str, Any]) -> None:
    """Record the validated LDSC estimates as structured scientific results."""
    observed = result["observed_metrics"]
    logger.record(
        "RESULT", "ldsc_observed_scale",
        h2=observed["h2"],
        intercept=observed["intercept"],
        ratio=observed["ratio"],
    )
    liability = result["liability_metrics"]
    if liability is not None:
        logger.record(
            "RESULT", "ldsc_liability_scale",
            h2=liability["h2"],
            intercept=liability["intercept"],
            ratio=liability["ratio"],
            sample_prevalence=result["sample_prevalence"],
            population_prevalence=result["population_prevalence"],
            sample_prevalence_source=result["sample_prevalence_source"],
        )


def _prevalence_source_label(source: str) -> str:
    """Render the recorded prevalence provenance in plain language."""
    labels = {
        "configuration_or_cli": "CLI or YAML configuration",
        "formatter_return_value": "formatter-returned value",
        "unset": "not supplied",
    }
    return labels.get(source, source.replace("_", " "))


def _render_summary(
    result: dict[str, Any],
    dataset: str,
    log_path: Path,
    label_width: int,
) -> str:
    """Render the validated LDSC findings for the default terminal summary."""
    observed = result["observed_metrics"]
    lines = [
        "",
        screen_line("analysis", "LDSC heritability summary", indent=2),
        screen_field(
            "info", "Dataset", dataset, indent=6, label_width=label_width,
        ),
        screen_field(
            "success", "Analysis status", "COMPLETED",
            indent=6, label_width=label_width,
        ),
        "",
        screen_line("genetic", "Main findings", indent=6),
        screen_field(
            "analysis", "Observed-scale h²", observed["h2"],
            indent=10, label_width=label_width,
        ),
        screen_field(
            "analysis", "Observed intercept", observed["intercept"],
            indent=10, label_width=label_width,
        ),
        screen_field(
            "analysis", "Observed ratio", observed["ratio"],
            indent=10, label_width=label_width,
        ),
        "",
    ]
    liability = result["liability_metrics"]
    if liability is None:
        lines.append(screen_field(
            "info", "Liability-scale h²",
            "not run because population prevalence was not provided",
            indent=10, label_width=label_width,
        ))
    else:
        lines.extend([
            screen_field(
                "analysis", "Liability-scale h²", liability["h2"],
                indent=10, label_width=label_width,
            ),
            screen_field(
                "analysis", "Liability intercept", liability["intercept"],
                indent=10, label_width=label_width,
            ),
            screen_field(
                "analysis", "Liability ratio", liability["ratio"],
                indent=10, label_width=label_width,
            ),
            "",
            screen_field(
                "info", "Sample prevalence",
                "%s (%s)" % (
                    result["sample_prevalence"],
                    _prevalence_source_label(result["sample_prevalence_source"]),
                ),
                indent=10, label_width=label_width,
            ),
            screen_field(
                "info", "Population prevalence",
                result["population_prevalence"],
                indent=10, label_width=label_width,
            ),
        ])
    lines.extend([
        "",
        screen_field(
            "success", "Observed LDSC results", result["h2_observed"],
            indent=6, label_width=label_width,
        ),
    ])
    if result["h2_liability"] is not None:
        lines.append(screen_field(
            "success", "Liability LDSC results", result["h2_liability"],
            indent=6, label_width=label_width,
        ))
    lines.extend([
        screen_field(
            "info", "Full PostGWAS log", log_path,
            indent=6, label_width=label_width,
        ),
        "",
    ])
    return "\n".join(lines)


def run_ldsc_direct(
    args: argparse.Namespace,
    ctx: dict[str, Any] | None = None,
    *,
    configuration=None,
) -> dict[str, Any]:
    """Resolve, validate, run, and atomically publish one LDSC analysis."""
    if configuration is None:
        try:
            configuration = resolve_ldsc_configuration(args)
        except BaseException as exc:
            fallback = load_configuration()
            output = Path(
                getattr(args, "output_directory", None)
                or fallback.run.output_directory
            ).expanduser().resolve()
            dataset = getattr(args, "dataset_id", None) or fallback.run.dataset_id
            log_path = configured_output_path(
                output,
                fallback.modules.ldsc.output_layout.service_log_file,
                error_type=LDSCError,
                dataset_id=dataset,
            )
            write_log_record(
                log_path,
                "ERROR",
                "LDSC configuration validation failed: %s: %s"
                % (type(exc).__name__, exc),
                sample_id=dataset,
                file_level=fallback.logging.file_level,
                screen_level=fallback.logging.console_level,
            )
            raise
    module = configuration.modules.ldsc
    output = Path(configuration.run.output_directory).expanduser().resolve()
    dataset = configuration.run.dataset_id
    output.mkdir(parents=True, exist_ok=True)
    final_prefix = configured_output_path(
        output, module.output_layout.output_prefix,
        error_type=LDSCError, dataset_id=dataset,
    )
    log_path = configured_output_path(
        output, module.output_layout.service_log_file,
        error_type=LDSCError, dataset_id=dataset,
    )
    logger = PipelineLogger(
        dataset,
        "run",
        str(log_path.parent),
        level=configuration.logging.file_level,
        screen_level=configuration.logging.console_level,
        log_path=str(log_path),
    )
    try:
        try:
            prevalence_source = _resolve_prevalence(module, ctx)
            ldsc_executable = resolve_executable(
                configuration.resources.executables.ldsc,
                "LDSC executable",
                error_type=LDSCError,
            )
            munge_executable = resolve_executable(
                configuration.resources.executables.munge_sumstats,
                "LDSC munge_sumstats executable",
                error_type=LDSCError,
            )
            ldsc_input = getattr(args, "ldsc_input", None)
            merge_alleles = getattr(args, "merge_alleles", None)
            reference = getattr(args, "ref_ld_chr", None)
            weights = getattr(args, "w_ld_chr", None)
            if any(
                value is None
                for value in (ldsc_input, merge_alleles, reference, weights)
            ):
                raise LDSCError(
                    "LDSC requires ldsc_input, merge_alleles, ref_ld_chr, and w_ld_chr."
                )

            logger.record(
                "PARAM", "ldsc_run",
                population=(
                    None if module.population is None else module.population.value
                ),
                genome_build=(
                    None if module.genome_build is None else module.genome_build.value
                ),
                chromosomes=module.chromosomes,
                minimum_info=module.minimum_info,
                minimum_maf=module.minimum_maf,
                minimum_n=module.minimum_n,
                chunksize=module.chunksize,
                keep_maf=module.keep_maf,
                intercept=module.intercept,
                two_step=module.two_step,
                chisq_max=module.chisq_max,
                n_blocks=module.n_blocks,
                use_m_5_50=module.use_m_5_50,
                print_covariance=module.print_covariance,
                print_delete_values=module.print_delete_values,
                sample_prevalence=module.sample_prevalence,
                population_prevalence=module.population_prevalence,
                sample_prevalence_source=prevalence_source,
            )
            resolved_path = configured_output_path(
                output, module.output_layout.resolved_config_file,
                error_type=LDSCError, dataset_id=dataset,
            )
            write_resolved_configuration(
                configuration,
                resolved_path,
                modules="ldsc",
                resource_paths=("executables.ldsc", "executables.munge_sumstats"),
            )
            logger.record("OUTPUT", "resolved_configuration", path=str(resolved_path))

            owned_final_paths = ldsc_owned_output_paths(final_prefix, module)
            existing = [path for path in owned_final_paths if path.exists()]
            if existing and not configuration.run.overwrite:
                raise LDSCError(
                    "Existing LDSC outputs were found, so PostGWAS stopped before "
                    "changing any files. Review these outputs: %s. Then rerun with "
                    "--overwrite to replace them, or choose a different "
                    "--output-directory."
                    % ", ".join(str(path) for path in existing)
                )
            staging_root = configured_output_path(
                output, module.output_layout.staging_directory,
                error_type=LDSCError, dataset_id=dataset,
            )
            staging_root.mkdir(parents=True, exist_ok=True)
            with tempfile.TemporaryDirectory(
                prefix="run_", dir=staging_root,
            ) as directory:
                staged_prefix = Path(directory) / final_prefix.name
                result = run_ldsc(
                    sumstats_tsv=ldsc_input,
                    output_prefix=staged_prefix,
                    merge_alleles=merge_alleles,
                    reference_ld_directory=reference,
                    weights_ld_directory=weights,
                    munge_executable=munge_executable,
                    ldsc_executable=ldsc_executable,
                    configuration=module,
                    logger=logger,
                )
                if configuration.run.overwrite:
                    for path in owned_final_paths:
                        path.unlink(missing_ok=True)
                result = _publish_staged_result(result, staged_prefix, final_prefix)
            result.update({
                "resolved_configuration": str(resolved_path),
                "service_log": str(log_path),
                "sample_prevalence": module.sample_prevalence,
                "population_prevalence": module.population_prevalence,
                "sample_prevalence_source": prevalence_source,
            })
            for name in (
                "munged_sumstats", "munge_log", "h2_observed", "h2_liability",
            ):
                if result[name] is not None:
                    logger.record("OUTPUT", name, path=result[name])
            _record_findings(logger, result)
            logger.record(
                "DONE", "ldsc_run",
                observed_log=result["h2_observed"],
                liability_log=result["h2_liability"],
            )
            if ctx is not None:
                ctx["heritability"] = result
            print(_render_summary(
                result,
                dataset,
                log_path,
                configuration.logging.terminal_label_width,
            ))
            return result
        except BaseException as exc:
            logger.record(
                "FAILED", "ldsc_run",
                error_type=type(exc).__name__, error=str(exc),
            )
            raise
    finally:
        logger.close()


__all__ = ["resolve_ldsc_configuration", "run_ldsc_direct"]
