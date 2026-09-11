"""Configuration boundary and publication service for LDSC heritability."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
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
from postgwas.core.ui import print_screen_block
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.modules.ldsc.ldsc_runner import (
    LDSCError,
    LDSCReferenceValidation,
    ldsc_owned_output_paths,
    run_ldsc,
    validate_ldsc_reference_resources,
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
    "ldsc_sample_prevalence_warning_threshold": (
        "sample_prevalence_comparison.warning_absolute_difference"
    ),
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


@dataclass(frozen=True)
class LDSCPipelineResources:
    """External LDSC resources validated before formatter table creation."""

    configuration: Any
    reference: LDSCReferenceValidation
    ldsc_executable: str
    munge_executable: str
    file_identities: tuple[PreflightFileIdentity, ...]


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


def _require_ldsc_runtime_arguments(
    args: argparse.Namespace,
    *,
    include_generated_input: bool,
) -> None:
    """Report all independently missing LDSC inputs through one contract."""
    requirements = [
        RequiredArgument("--merge-alleles", None, getattr(args, "merge_alleles", None)),
        RequiredArgument("--ref-ld-chr", None, getattr(args, "ref_ld_chr", None)),
        RequiredArgument("--w-ld-chr", None, getattr(args, "w_ld_chr", None)),
    ]
    if include_generated_input:
        requirements.insert(0, RequiredArgument(
            "--ldsc-input", None, getattr(args, "ldsc_input", None),
        ))
    require_resolved_arguments(requirements)


def preflight_ldsc_pipeline(
    args: argparse.Namespace,
    *,
    preflight_evidence=None,
) -> PipelinePreflightEvidence:
    """Validate all external LDSC resources before VCF formatting starts."""
    entry_vcf = require_pipeline_input_vcf(preflight_evidence)
    configuration = resolve_ldsc_configuration(args)
    module = configuration.modules.ldsc
    _require_ldsc_runtime_arguments(args, include_generated_input=False)

    observed_build = str(entry_vcf["harmonised"]["genome_build"])
    configured_build = (
        None if module.genome_build is None else module.genome_build.value
    )
    if configured_build is not None and configured_build != observed_build:
        raise LDSCError(
            "The harmonised GWAS-VCF declares genome build %s, but the LDSC "
            "reference release is declared as %s. These resources must use "
            "the same genome build." % (observed_build, configured_build)
        )

    reference = validate_ldsc_reference_resources(
        getattr(args, "merge_alleles", None),
        getattr(args, "ref_ld_chr", None),
        getattr(args, "w_ld_chr", None),
        module,
    )
    if reference.merge_alleles is None:
        raise LDSCError("LDSC merge-alleles validation returned no file.")
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
    external_files = [
        reference.merge_alleles,
        *reference.required_files,
        ldsc_executable,
        munge_executable,
    ]
    resources = LDSCPipelineResources(
        configuration=configuration,
        reference=reference,
        ldsc_executable=ldsc_executable,
        munge_executable=munge_executable,
        file_identities=capture_preflight_file_identities(
            external_files,
            error_type=LDSCError,
            label="LDSC resource",
        ),
    )
    return pipeline_preflight_evidence(
        "heritability",
        preflight_evidence,
        resources=resources,
        deferred_checks=(
            "Validate the formatter-created LDSC summary-statistics table.",
            "Validate retained HapMap3 SNP and allele compatibility during munging.",
        ),
    )


def _ldsc_pipeline_execution_configuration(args, resources):
    """Retain preflight settings while applying the pipeline stage directory."""
    run = resources.configuration.run.model_copy(update={
        "output_directory": Path(args.output_directory).expanduser().resolve(),
    })
    return resources.configuration.model_copy(update={"run": run}, deep=True)


def _formatter_ldsc_result(
    ctx: dict[str, Any] | None,
    *,
    required: bool,
) -> dict[str, Any] | None:
    """Return the formatter's LDSC result without reading its output table."""
    try:
        formatter_result = ctx["formatter"]["ldsc"]  # type: ignore[index]
    except (KeyError, TypeError) as exc:
        if required:
            raise LDSCError(
                "Population prevalence requires sample prevalence. Supply --samp-prev "
                "for direct execution, or run the LDSC formatter first so its returned "
                "sample_prev artifact is available."
            ) from exc
        return None
    if not isinstance(formatter_result, dict):
        raise LDSCError(
            "The formatter returned an invalid LDSC result; expected a mapping "
            "containing sample_prev metadata."
        )
    return formatter_result


def _formatter_sample_prevalence(ctx: dict[str, Any] | None) -> float:
    """Read, but never derive, sample prevalence from the formatter result."""
    formatter_result = _formatter_ldsc_result(ctx, required=True)
    assert formatter_result is not None
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


def _resolve_prevalence(
    module,
    ctx: dict[str, Any] | None,
) -> tuple[str, dict[str, Any] | None]:
    """Resolve liability prevalence and compare a pipeline override to the data."""
    source = (
        "configuration_or_cli"
        if module.sample_prevalence is not None
        else "unset"
    )
    if module.sample_prevalence is None:
        if module.population_prevalence is None:
            return source, None
        module.sample_prevalence = _formatter_sample_prevalence(ctx)
        source = "gwas_vcf_case_fraction"
        return source, None

    formatter_result = _formatter_ldsc_result(ctx, required=False)
    if formatter_result is None or formatter_result.get("sample_prev") is None:
        return source, None
    gwas_vcf_case_fraction = _formatter_sample_prevalence(ctx)
    provided_sample_prevalence = float(module.sample_prevalence)
    absolute_difference_decimal = abs(
        Decimal(str(provided_sample_prevalence))
        - Decimal(str(gwas_vcf_case_fraction))
    )
    warning_threshold_decimal = Decimal(str(
        module.sample_prevalence_comparison.warning_absolute_difference
    ))
    absolute_difference = float(absolute_difference_decimal)
    return source, {
        "matches": absolute_difference_decimal == 0,
        "exceeds_warning_threshold": (
            absolute_difference_decimal > warning_threshold_decimal
        ),
        "provided_sample_prevalence": provided_sample_prevalence,
        "gwas_vcf_case_fraction": gwas_vcf_case_fraction,
        "absolute_difference": absolute_difference,
        "warning_absolute_difference": float(warning_threshold_decimal),
        "liability_scale_requested": module.population_prevalence is not None,
        "gwas_vcf_aggregation": formatter_result.get(
            "sample_prevalence_aggregation"
        ),
        "gwas_vcf_variant_count": formatter_result.get(
            "sample_prevalence_variants"
        ),
        "gwas_vcf_case_fraction_minimum": formatter_result.get(
            "sample_prevalence_minimum"
        ),
        "gwas_vcf_case_fraction_maximum": formatter_result.get(
            "sample_prevalence_maximum"
        ),
        "gwas_vcf_case_count_minimum": formatter_result.get(
            "sample_prevalence_case_count_minimum"
        ),
        "gwas_vcf_case_count_maximum": formatter_result.get(
            "sample_prevalence_case_count_maximum"
        ),
        "gwas_vcf_control_count_minimum": formatter_result.get(
            "sample_prevalence_control_count_minimum"
        ),
        "gwas_vcf_control_count_maximum": formatter_result.get(
            "sample_prevalence_control_count_maximum"
        ),
    }


def _record_prevalence_comparison(
    logger: PipelineLogger,
    comparison: dict[str, Any] | None,
) -> None:
    """Record the provided prevalence against the GWAS-VCF case fraction."""
    if comparison is None:
        return
    exceeds_threshold = bool(comparison["exceeds_warning_threshold"])
    logger.record(
        "WARNING" if exceeds_threshold else "PASS",
        "ldsc_sample_prevalence_comparison",
        result=(
            "warning_threshold_exceeded"
            if exceeds_threshold
            else "within_warning_threshold"
        ),
        matches=comparison["matches"],
        provided_sample_prevalence=comparison["provided_sample_prevalence"],
        gwas_vcf_case_fraction=comparison["gwas_vcf_case_fraction"],
        absolute_difference=comparison["absolute_difference"],
        warning_absolute_difference=comparison[
            "warning_absolute_difference"
        ],
        liability_scale_requested=comparison["liability_scale_requested"],
        gwas_vcf_aggregation=comparison["gwas_vcf_aggregation"],
        gwas_vcf_variant_count=comparison["gwas_vcf_variant_count"],
        gwas_vcf_case_fraction_minimum=comparison[
            "gwas_vcf_case_fraction_minimum"
        ],
        gwas_vcf_case_fraction_maximum=comparison[
            "gwas_vcf_case_fraction_maximum"
        ],
        gwas_vcf_case_count_minimum=comparison[
            "gwas_vcf_case_count_minimum"
        ],
        gwas_vcf_case_count_maximum=comparison[
            "gwas_vcf_case_count_maximum"
        ],
        gwas_vcf_control_count_minimum=comparison[
            "gwas_vcf_control_count_minimum"
        ],
        gwas_vcf_control_count_maximum=comparison[
            "gwas_vcf_control_count_maximum"
        ],
        action=(
            "explicit_cli_or_yaml_value_retained; the run continues; to use the "
            "GWAS-VCF case fraction stop this run and rerun without --samp-prev "
            "and with modules.ldsc.sample_prevalence unset"
        ),
    )


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
            conversion_of="ldsc_observed_scale.h2",
            regression_metrics="shared_with_observed_scale",
            sample_prevalence=result["sample_prevalence"],
            population_prevalence=result["population_prevalence"],
            sample_prevalence_source=result["sample_prevalence_source"],
        )


def _prevalence_source_label(source: str) -> str:
    """Render the recorded prevalence provenance in plain language."""
    labels = {
        "configuration_or_cli": "CLI or YAML configuration",
        "gwas_vcf_case_fraction": "GWAS-VCF case fraction",
        "unset": "not supplied",
    }
    return labels.get(source, source.replace("_", " "))


def _reported_range(minimum: Any, maximum: Any) -> str | None:
    """Render a precomputed scientific range without recalculating it."""
    if minimum is None or maximum is None:
        return None
    return "%s to %s" % (minimum, maximum)


def _prevalence_comparison_lines(
    comparison: dict[str, Any] | None,
    label_width: int,
) -> list[str]:
    """Render the pre-run prevalence comparison and selected LDSC value."""
    if comparison is None:
        return []
    exceeds_threshold = bool(comparison["exceeds_warning_threshold"])
    emphasis = "warning" if exceeds_threshold else "info"
    gwas_vcf_description = str(comparison["gwas_vcf_case_fraction"])
    aggregation = comparison["gwas_vcf_aggregation"]
    variants = comparison["gwas_vcf_variant_count"]
    provenance = [
        value
        for value in (
            None if aggregation is None else "%s aggregation" % aggregation,
            None if variants is None else "%s variants" % variants,
        )
        if value is not None
    ]
    if provenance:
        gwas_vcf_description += " (%s)" % ", ".join(provenance)
    lines = [
        "",
        screen_line(
            emphasis,
            (
                "Sample prevalence warning"
                if exceeds_threshold
                else "Sample prevalence comparison"
            ),
            indent=6,
        ),
        screen_field(
            emphasis, "Provided sample prevalence",
            comparison["provided_sample_prevalence"],
            indent=10, label_width=label_width,
        ),
        screen_field(
            "info", "GWAS-VCF case fraction", gwas_vcf_description,
            indent=10, label_width=label_width,
        ),
        screen_field(
            "info", "Absolute difference", comparison["absolute_difference"],
            indent=10, label_width=label_width,
        ),
        screen_field(
            "info", "Warning threshold",
            comparison["warning_absolute_difference"],
            indent=10, label_width=label_width,
        ),
        screen_field(
            emphasis,
            "Comparison status",
            (
                "WARNING — difference is above the configured threshold"
                if exceeds_threshold
                else "difference is at or below the configured threshold"
            ),
            indent=10,
            label_width=label_width,
        ),
    ]
    prevalence_range = _reported_range(
        comparison["gwas_vcf_case_fraction_minimum"],
        comparison["gwas_vcf_case_fraction_maximum"],
    )
    if prevalence_range is not None:
        lines.append(screen_field(
            "info", "GWAS-VCF case-fraction range", prevalence_range,
            indent=10, label_width=label_width,
        ))
    case_count_range = _reported_range(
        comparison["gwas_vcf_case_count_minimum"],
        comparison["gwas_vcf_case_count_maximum"],
    )
    if case_count_range is not None:
        lines.append(screen_field(
            "info", "GWAS-VCF case-count range", case_count_range,
            indent=10, label_width=label_width,
        ))
    control_count_range = _reported_range(
        comparison["gwas_vcf_control_count_minimum"],
        comparison["gwas_vcf_control_count_maximum"],
    )
    if control_count_range is not None:
        lines.append(screen_field(
            "info", "GWAS-VCF control-count range", control_count_range,
            indent=10, label_width=label_width,
        ))
    provided = comparison["provided_sample_prevalence"]
    gwas_vcf = comparison["gwas_vcf_case_fraction"]
    if comparison["liability_scale_requested"]:
        action = (
            "LDSC WILL USE %s (the provided sample prevalence), NOT %s (the "
            "GWAS-VCF case fraction), to calculate liability-scale h². "
            % (provided, gwas_vcf)
        )
        action += (
            "If this is intentional, allow the run to continue. Otherwise, "
            "stop now and rerun without --samp-prev "
            "and ensure the YAML sample_prevalence setting is absent or null."
        )
    else:
        action = (
            "Liability-scale h² WILL NOT RUN because population prevalence was not "
            "provided. PostGWAS retained %s (the provided sample prevalence) "
            "and did NOT replace it with %s (the GWAS-VCF case fraction). To "
            "use the GWAS-VCF value in a future liability-scale run, omit "
            "--samp-prev and ensure the YAML sample_prevalence setting is "
            "absent or null."
            % (provided, gwas_vcf)
        )
    lines.append(screen_field(
        "attention" if exceeds_threshold else "decision",
        "Important — value selected",
        action,
        indent=10,
        label_width=label_width,
    ))
    return lines


def _render_summary(
    result: dict[str, Any],
    dataset: str,
    log_path: Path,
    label_width: int,
) -> str:
    """Render the validated LDSC findings for the default terminal summary."""
    observed = result["observed_metrics"]
    prevalence_comparison = result["sample_prevalence_comparison"]
    has_prevalence_warning = bool(
        prevalence_comparison is not None
        and prevalence_comparison["exceeds_warning_threshold"]
    )
    lines = [
        "",
        screen_line("analysis", "LDSC heritability summary", indent=2),
        screen_field(
            "info", "Dataset", dataset, indent=6, label_width=label_width,
        ),
        screen_field(
            "warning" if has_prevalence_warning else "success",
            "Analysis status",
            (
                "COMPLETED WITH SCIENTIFIC WARNINGS"
                if has_prevalence_warning
                else "COMPLETED"
            ),
            indent=6, label_width=label_width,
        ),
    ]
    lines.extend([
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
    ])
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
                "analysis", "Liability-scale h² (converted)", liability["h2"],
                indent=10, label_width=label_width,
            ),
            screen_field(
                "success", "Conversion check",
                "PASS — intercept and ratio unchanged",
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
    pipeline_resources: LDSCPipelineResources | None = None,
) -> dict[str, Any]:
    """Resolve, validate, run, and atomically publish one LDSC analysis."""
    if pipeline_resources is not None and configuration is not None:
        raise LDSCError(
            "Pass either pipeline_resources or configuration to LDSC, not both."
        )
    if configuration is None:
        try:
            configuration = (
                _ldsc_pipeline_execution_configuration(args, pipeline_resources)
                if pipeline_resources is not None
                else resolve_ldsc_configuration(args)
            )
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
    if pipeline_resources is not None:
        require_unchanged_preflight_files(
            pipeline_resources.file_identities,
            error_type=LDSCError,
            label="LDSC resource",
        )
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
            prevalence_source, prevalence_comparison = _resolve_prevalence(
                module, ctx,
            )
            _record_prevalence_comparison(logger, prevalence_comparison)
            if prevalence_comparison is not None:
                print_screen_block("\n".join(_prevalence_comparison_lines(
                    prevalence_comparison,
                    configuration.logging.terminal_label_width,
                )))
            ldsc_executable = (
                pipeline_resources.ldsc_executable
                if pipeline_resources is not None
                else resolve_executable(
                    configuration.resources.executables.ldsc,
                    "LDSC executable",
                    error_type=LDSCError,
                )
            )
            munge_executable = (
                pipeline_resources.munge_executable
                if pipeline_resources is not None
                else resolve_executable(
                    configuration.resources.executables.munge_sumstats,
                    "LDSC munge_sumstats executable",
                    error_type=LDSCError,
                )
            )
            ldsc_input = getattr(args, "ldsc_input", None)
            _require_ldsc_runtime_arguments(args, include_generated_input=True)
            reference_validation = (
                pipeline_resources.reference
                if pipeline_resources is not None else None
            )
            merge_alleles = (
                reference_validation.merge_alleles
                if reference_validation is not None
                else getattr(args, "merge_alleles")
            )
            reference = (
                reference_validation.reference_directory
                if reference_validation is not None
                else getattr(args, "ref_ld_chr")
            )
            weights = (
                reference_validation.weights_directory
                if reference_validation is not None
                else getattr(args, "w_ld_chr")
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
                sample_prevalence_warning_absolute_difference=(
                    module.sample_prevalence_comparison.warning_absolute_difference
                ),
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
                    reference_validation=reference_validation,
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
                "sample_prevalence_comparison": prevalence_comparison,
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
            print_screen_block(_render_summary(
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


__all__ = [
    "LDSCPipelineResources",
    "preflight_ldsc_pipeline",
    "resolve_ldsc_configuration",
    "run_ldsc_direct",
]
