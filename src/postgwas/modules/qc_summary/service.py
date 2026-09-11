"""Configuration, provenance, and execution boundary for GWAS-VCF QC."""

from __future__ import annotations

import argparse
from collections.abc import Callable
import json
from pathlib import Path
from typing import Any, Sequence

from postgwas.config import (
    load_configuration,
    load_run_configuration_for_module,
    write_resolved_configuration,
)
from postgwas.config.cli_overrides import (
    explicit_overrides,
    variant_qc_policy_override_paths,
)
from postgwas.config.models.modules.qc_summary import QCSummaryConfig
from postgwas.core.checkpointing import software_identity
from postgwas.core.completion import (
    apply_completion_restart,
    configuration_digest,
    resolve_completion_resume,
    write_completion_manifest,
)
from postgwas.core.paths import (
    configured_output_path,
    require_nonempty_file,
    resolve_executable,
    validate_filename_component,
)
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.core.processes import run_checked_command
from postgwas.core.screen_logging import screen_recording_active
from postgwas.core.ui import StageProgress, print_screen_block, screen_line
from postgwas.core.vcf import (
    read_vcf_header,
    validate_postgwas_vcf_provenance,
    validate_vcf_header_contract,
)
from postgwas.modules.qc_summary.assessment import (
    VcfAssessmentError,
    VcfAssessmentHeaderValidation,
    qc_assessment_stage,
    qc_assessment_runtime,
    resolve_qc_field_labels,
    resolve_qc_rule_contract,
    run_vcf_qc_assessment,
    validate_qc_rule_contract,
    validate_vcf_assessment_header,
)
from postgwas.modules.qc_summary.reporting import (
    qc_final_screen_lines,
    qc_input_validation_screen_lines,
    qc_rules_applied_screen_lines,
    qc_summary_lines,
    qc_summary_screen_lines,
)


MODULE_CLI_OVERRIDES = {
    "vcf": "inputs.vcf",
    "dataset_id": "inputs.dataset_id",
    "output_directory": "output_directory",
    "reference_af_column": "reference_af_column",
    **variant_qc_policy_override_paths(prefix="rules."),
    "sample_size_reference_quantile": "rules.sample_size_reference_quantile",
    "sample_size_minimum_fraction": (
        "rules.sample_size_minimum_fraction_of_reference"
    ),
}

GLOBAL_CLI_OVERRIDES = {
    "dataset_id": "run.dataset_id",
    "output_directory": "run.output_directory",
    "threads": "execution.threads",
    "memory_gb": "execution.memory_gb",
    "seed": "execution.random_seed",
    "resume": "run.resume",
    "overwrite": "run.overwrite",
}


def resolve_qc_summary_configuration(args: argparse.Namespace):
    """Resolve direct or pipeline CLI values through the canonical QC schema."""
    module_overrides = {
        key: value
        for key, value in explicit_overrides(args, MODULE_CLI_OVERRIDES).items()
        if value is not None
    }
    global_overrides = {
        key: value
        for key, value in explicit_overrides(args, GLOBAL_CLI_OVERRIDES).items()
        if value is not None
    }
    return load_run_configuration_for_module(
        "qc_summary",
        getattr(args, "run_config", None),
        module_overrides=module_overrides,
        global_overrides=global_overrides,
    )


def _output_paths(
    output_directory: Path,
    dataset_id: str,
    genome_build: str,
    configuration: QCSummaryConfig,
) -> dict[str, Path]:
    layout = configuration.output_layout
    return {
        name: configured_output_path(
            output_directory,
            pattern,
            error_type=VcfAssessmentError,
            dataset_id=dataset_id,
            build=genome_build,
        )
        for name, pattern in {
            "summary": layout.metric_report,
            "rules": layout.rule_report,
            "json": layout.assessment_json,
            "csv": layout.summary_csv,
            "html": layout.html_report,
            "resolved_config": layout.resolved_config_file,
            "log": layout.service_log_file,
            "completion": layout.completion_manifest,
        }.items()
    }


def run_qc_assessment(
    *,
    vcf_path: str | Path,
    output_directory: str | Path,
    dataset_id: str,
    external_af_name: str,
    configuration: QCSummaryConfig,
    bcftools_bin: str,
    genome_build_header: str,
    supported_genome_builds: Sequence[str],
    threads: int,
    header_validation: VcfAssessmentHeaderValidation | None = None,
    provenance_headers: dict[str, str] | None = None,
    logger=None,
    stage_progress: StageProgress | None = None,
    stage_number_offset: int = 0,
    rule_contract: dict[str, Any] | None = None,
    metrics_completed: Callable[[dict[str, Any]], None] | None = None,
) -> dict[str, Any]:
    """Run the shared two-pass assessment used by direct and harmonisation."""
    return run_vcf_qc_assessment(
        vcf_path=vcf_path,
        output_directory=output_directory,
        dataset_id=dataset_id,
        external_af_name=external_af_name,
        configuration=configuration,
        bcftools_bin=bcftools_bin,
        genome_build_header=genome_build_header,
        supported_genome_builds=supported_genome_builds,
        threads=threads,
        header_validation=header_validation,
        provenance_headers=provenance_headers,
        logger=logger,
        stage_progress=stage_progress,
        stage_number_offset=stage_number_offset,
        rule_contract=rule_contract,
        metrics_completed=metrics_completed,
    )


def _load_completed_assessment(path: Path) -> dict[str, Any]:
    try:
        assessment = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise VcfAssessmentError(
            "Cannot read the completed QC assessment %s: %s" % (path, exc)
        ) from exc
    if not isinstance(assessment, dict) or not assessment.get("accounting_balanced"):
        raise VcfAssessmentError(
            "Completed QC assessment is missing balanced variant accounting: %s"
            % path
        )
    return assessment


def _emit_summary(
    assessment: dict[str, Any],
    configuration,
    logger,
    *,
    include_stage_cards: bool = True,
) -> None:
    detailed_report = "\n".join(qc_summary_lines(
        assessment,
        label_width=configuration.logging.terminal_label_width,
    ))
    final_report = "\n".join(qc_final_screen_lines(
        assessment, label_width=configuration.logging.terminal_label_width,
    ))
    renderer = (
        qc_summary_screen_lines if include_stage_cards else qc_final_screen_lines
    )
    screen_report = "\n".join(renderer(
        assessment, label_width=configuration.logging.terminal_label_width,
    ))
    logger.log(
        "RESULT",
        "%s\n\n%s\n%s"
        % (
            final_report,
            screen_line("analysis", "Detailed QC audit", indent=4),
            detailed_report,
        ),
        wrap=False,
        screen=False,
    )
    if configuration.logging.show_screen or screen_recording_active():
        print_screen_block("\n" + screen_report)


def _emit_stage_block(
    progress: StageProgress,
    logger: PipelineLogger,
    lines: list[str],
) -> None:
    """Write one filtering-style stage block to the screen and canonical log."""
    text = "\n".join(lines)
    logger.log("RESULT", text, wrap=False, screen=False)
    if progress.enabled:
        progress.print_block(text)


def _fallback_log(
    args: argparse.Namespace,
    exc: BaseException,
    *,
    configuration=None,
) -> None:
    fallback = configuration or load_configuration()
    module = fallback.modules.qc_summary
    output = (
        getattr(args, "output_directory", None)
        or module.output_directory
        or fallback.run.output_directory
    )
    if output is None:
        return
    candidate_dataset = (
        getattr(args, "dataset_id", None)
        or module.inputs.dataset_id
        or fallback.run.dataset_id
        or "qc_summary"
    )
    try:
        dataset = validate_filename_component(
            candidate_dataset,
            "dataset_id",
            error_type=VcfAssessmentError,
        )
    except VcfAssessmentError:
        dataset = "qc_summary"
    log_path = configured_output_path(
        output,
        module.output_layout.preflight_log_file,
        dataset_id=dataset,
    )
    write_log_record(
        log_path,
        "ERROR",
        "QC summary preflight status=FAILED exception_type=%s reason=%s"
        % (type(exc).__name__, str(exc) or "no error message was provided"),
        sample_id=dataset,
        file_level=fallback.logging.file_level,
        screen_level=fallback.logging.console_level,
    )


def run_qc_summary_direct(
    args: argparse.Namespace,
    ctx: dict[str, Any] | None = None,
    *,
    configuration=None,
) -> dict[str, Any]:
    """Assess one genotype-free GWAS-VCF and publish validated QC reports."""
    resolved_configuration = configuration
    try:
        configuration = configuration or resolve_qc_summary_configuration(args)
        resolved_configuration = configuration
        module = configuration.modules.qc_summary
        if module.inputs.vcf is None:
            raise VcfAssessmentError(
                "A harmonised GWAS-VCF is required; provide --vcf PATH."
            )
        vcf = require_nonempty_file(
            module.inputs.vcf,
            "harmonised GWAS-VCF",
            error_type=VcfAssessmentError,
        )
        dataset_id = validate_filename_component(
            module.inputs.dataset_id or configuration.run.dataset_id,
            "dataset_id",
            error_type=VcfAssessmentError,
        )
        configured_output = module.output_directory or configuration.run.output_directory
        if configured_output is None:
            raise VcfAssessmentError(
                "An output directory is required; provide --output-directory PATH."
            )
        output_directory = Path(configured_output).expanduser().resolve()
        bcftools = resolve_executable(
            configuration.resources.executables.bcftools,
            "bcftools executable",
            error_type=VcfAssessmentError,
        )
        vcf_processing = configuration.modules.harmonisation.vcf_processing
        required_provenance_headers = (
            configuration.modules.formatting.input_contract.provenance_headers.model_dump()
        )
        genome_build_header = vcf_processing.genome_build_header
        supported_genome_builds = tuple(configuration.resources.genomes)
        header = read_vcf_header(
            vcf,
            bcftools,
            error_type=VcfAssessmentError,
        )
        try:
            genome_build, _ = validate_vcf_header_contract(
                header=header,
                genome_build_header=genome_build_header,
                supported_genome_builds=supported_genome_builds,
                required_fields=(),
            )
        except ValueError as exc:
            raise VcfAssessmentError(
                "QC assessment cannot infer a supported genome build from the "
                "VCF header: %s" % exc
            ) from exc
        paths = _output_paths(output_directory, dataset_id, genome_build, module)
    except BaseException as exc:
        _fallback_log(args, exc, configuration=resolved_configuration)
        raise

    output_directory.mkdir(parents=True, exist_ok=True)
    resume_candidate = bool(
        configuration.run.resume
        and not configuration.run.overwrite
        and paths["completion"].is_file()
    )
    progress_total = 4 if resume_candidate else 6
    progress = StageProgress(
        (
            "GWAS-VCF QC validated-resume progress"
            if resume_candidate
            else "Summary-statistics GWAS-VCF QC progress"
        ),
        enabled=configuration.logging.show_progress,
        outcome_label_width=configuration.logging.terminal_label_width,
    )
    logger = PipelineLogger(
        dataset_id,
        "run",
        str(paths["log"].parent),
        level=configuration.logging.file_level,
        screen_level=configuration.logging.console_level,
        log_path=str(paths["log"]),
    )
    reports = {
        name: paths[name]
        for name in ("summary", "rules", "json", "csv", "html")
    }
    try:
        with qc_assessment_stage(
            progress,
            logger,
            1,
            progress_total,
            "Validate the input summary-statistics GWAS-VCF",
        ):
            logger.record(
                "COMMAND",
                "vcf_header_read",
                argv=[bcftools, "view", "--header-only", str(vcf)],
                genome_build=genome_build,
                genome_build_source="vcf_header",
            )
            validate_postgwas_vcf_provenance(
                header,
                required_provenance_headers,
                vcf_path=vcf,
                logger=logger,
                error_type=VcfAssessmentError,
            )
            bcftools_version_output = run_checked_command(
                [bcftools, "--version"],
                "Reading the bcftools version",
                logger=logger,
                error_type=VcfAssessmentError,
                timeout_seconds=configuration.execution.timeout_seconds,
            )
            bcftools_version = bcftools_version_output.splitlines()[0].strip()
            if not bcftools_version:
                raise VcfAssessmentError(
                    "bcftools --version returned no version information"
                )
            provenance_headers = dict(
                configuration.modules.harmonisation.vcf_processing.provenance.headers
            )
            header_validation = validate_vcf_assessment_header(
                vcf_path=vcf,
                external_af_name=module.reference_af_column,
                vcf_fields=module.vcf_fields.model_dump(),
                bcftools_bin=bcftools,
                genome_build_header=genome_build_header,
                supported_genome_builds=supported_genome_builds,
                header=header,
                provenance_headers=provenance_headers,
                logger=logger,
            )
            field_labels = resolve_qc_field_labels(
                module.vcf_fields.model_dump(), module.reference_af_column,
            )
        validation_evidence = {
            "raw_vcf": str(vcf),
            "genome_build": genome_build,
            "vcf_header_validation": header_validation.as_report(),
            "vcf_provenance": header_validation.provenance_report(),
            "field_labels": field_labels,
        }
        _emit_stage_block(progress, logger, qc_input_validation_screen_lines(
            validation_evidence,
            label_width=configuration.logging.terminal_label_width,
            include_total=False,
        ))

        with qc_assessment_stage(
            progress,
            logger,
            2,
            progress_total,
            "Resolve and validate the configured QC rules",
        ):
            rule_contract = resolve_qc_rule_contract(
                module, genome_build, module.reference_af_column,
            )
        _emit_stage_block(progress, logger, qc_rules_applied_screen_lines(
            rule_contract,
            label_width=configuration.logging.terminal_label_width,
        ))

        aggregation_runtime = qc_assessment_runtime()
        software = software_identity(bcftools)
        software["polars"] = aggregation_runtime["polars_version"]
        digest = configuration_digest({
            "qc_summary": module.model_dump(mode="json"),
            "genome_build": genome_build,
            "genome_build_source": "vcf_header",
            "genome_build_header": genome_build_header,
            "supported_genome_builds": supported_genome_builds,
            "bcftools": bcftools,
            "bcftools_version": bcftools_version,
            "aggregation": aggregation_runtime,
            "required_provenance_headers": required_provenance_headers,
            "provenance_headers": provenance_headers,
            "software": software,
        })
        logger.record(
            "PARAM",
            "qc_summary_run",
            dataset_id=dataset_id,
            genome_build=genome_build,
            genome_build_source="vcf_header",
            vcf=vcf,
            reference_af_column=module.reference_af_column,
            maximum_af_difference=module.rules.maximum_af_difference,
            bcftools=bcftools,
            bcftools_version=bcftools_version,
            resume=configuration.run.resume,
            overwrite=configuration.run.overwrite,
            threads=configuration.execution.threads,
        )
        logger.record(
            "PARAM",
            "qc_summary_runtime",
            postgwas_version=software["postgwas"],
            python_version=software["python"],
            requested_threads=configuration.execution.threads,
            thread_budget_enforcement="spawned_worker_pre_import",
            **aggregation_runtime,
        )
        logger.record(
            "PARAM", "qc_summary_rules",
            **module.rules.model_dump(mode="json"),
        )
        logger.record(
            "PARAM", "qc_summary_vcf_fields",
            **module.vcf_fields.model_dump(mode="json"),
        )
        logger.record(
            "PARAM",
            "qc_summary_resolved_policy",
            active_rule_count=rule_contract["active_rule_count"],
            active_rules=",".join(
                rule["key"] for rule in rule_contract["rules"]
            ),
            inactive_rules=",".join(
                rule["key"] for rule in rule_contract["inactive_rules"]
            ),
        )
        for number, item in enumerate(rule_contract["decision_plan"], 1):
            logger.record(
                "PARAM",
                "qc_summary_decision_plan",
                number=number,
                action=item["action"],
                rule=item["rule_key"],
                condition=item["label"],
            )
        logger.record(
            "INPUT", "qc_summary_vcf",
            path=vcf,
            size_bytes=vcf.stat().st_size,
        )
        if resume_candidate:
            with qc_assessment_stage(
                progress,
                logger,
                3,
                progress_total,
                "Validate the completed QC report set",
            ):
                decision = resolve_completion_resume(
                    paths["completion"],
                    dataset_id=dataset_id,
                    module="qc_summary",
                    genome_build=genome_build,
                    configuration_sha256=digest,
                    inputs={"vcf": vcf},
                    outputs=reports,
                    resume_policy=configuration.run.resume_policy,
                    error_type=VcfAssessmentError,
                )
            if decision.action == "resume":
                with qc_assessment_stage(
                    progress,
                    logger,
                    4,
                    progress_total,
                    "Load and verify the completed QC assessment",
                ):
                    assessment = _load_completed_assessment(paths["json"])
                    if assessment.get("vcf_header_validation") != (
                        header_validation.as_report()
                    ):
                        raise VcfAssessmentError(
                            "The completed QC assessment header evidence does "
                            "not match the validated input GWAS-VCF."
                        )
                    validate_qc_rule_contract(assessment, rule_contract)
                    assessment["decision_plan"] = rule_contract["decision_plan"]
                    presentation_by_key = {
                        rule["key"]: rule for rule in rule_contract["rules"]
                    }
                    for rule in assessment["rules"]:
                        presentation = presentation_by_key[rule["key"]]
                        rule["display_group"] = presentation["display_group"]
                        rule["display_group_kind"] = presentation[
                            "display_group_kind"
                        ]
                logger.record(
                    "SKIP",
                    "qc_summary_run",
                    reason="provenance_validated_complete_outputs",
                )
                _emit_summary(
                    assessment,
                    configuration,
                    logger,
                    include_stage_cards=not progress.enabled,
                )
                if ctx is not None:
                    ctx["qc_summary"] = assessment
                return assessment
            with qc_assessment_stage(
                progress,
                logger,
                4,
                progress_total,
                "Prepare a clean restart from the invalidated completion",
            ):
                apply_completion_restart(
                    decision,
                    output_root=output_directory,
                    manifest=paths["completion"],
                    logger=logger,
                    operation="qc_summary_resume",
                    error_type=VcfAssessmentError,
                )
            progress = StageProgress(
                "Restarted GWAS-VCF QC analysis progress",
                enabled=configuration.logging.show_progress,
                outcome_label_width=configuration.logging.terminal_label_width,
            )

        existing = [path for path in (*reports.values(), paths["completion"]) if path.exists()]
        if existing and not configuration.run.overwrite:
            raise VcfAssessmentError(
                "Existing or incomplete QC-summary outputs were found; use --resume "
                "for a matching completed run or --overwrite to replace them: %s"
                % ", ".join(str(path) for path in existing)
            )
        if configuration.run.overwrite:
            for path in (*reports.values(), paths["completion"]):
                if path.exists() and not path.is_file():
                    raise VcfAssessmentError(
                        "Configured QC output is not a file: %s" % path
                    )
                path.unlink(missing_ok=True)

        write_resolved_configuration(
            configuration,
            paths["resolved_config"],
            modules="qc_summary",
            resource_paths=("executables.bcftools",),
        )

        assessment = run_qc_assessment(
            vcf_path=vcf,
            output_directory=output_directory,
            dataset_id=dataset_id,
            external_af_name=module.reference_af_column,
            configuration=module,
            bcftools_bin=bcftools,
            genome_build_header=genome_build_header,
            supported_genome_builds=supported_genome_builds,
            threads=configuration.execution.threads,
            header_validation=header_validation,
            provenance_headers=provenance_headers,
            logger=logger,
            stage_progress=progress,
            stage_number_offset=(0 if resume_candidate else 2),
            rule_contract=rule_contract,
        )
        write_completion_manifest(
            paths["completion"],
            dataset_id=dataset_id,
            module="qc_summary",
            genome_build=genome_build,
            configuration_sha256=digest,
            inputs={"vcf": vcf},
            outputs=reports,
            metrics={
                "raw_variants": assessment["raw"]["num_records"],
                "qc_passed_variants": assessment["qc_passed"]["num_records"],
                "excluded_variants": assessment["excluded_total"],
                "variants_unique_to_one_rule": assessment[
                    "unique_rule_only_total"
                ],
                "variants_matching_multiple_rules": assessment[
                    "overlap_variants"
                ],
                "scientific_provenance_status": assessment[
                    "vcf_provenance"
                ]["status"],
                "genome_build_source": "vcf_header",
                "sample_size_reference_quantile": assessment[
                    "sample_size_reference_quantile"
                ],
                "sample_size_reference_value": assessment["raw"][
                    "effective_sample_size_reference_quantile_value"
                ],
                "sample_size_minimum_fraction_of_reference": assessment[
                    "sample_size_minimum_fraction_of_reference"
                ],
                "sample_size_minimum_threshold": assessment["raw"][
                    "effective_sample_size_minimum_threshold"
                ],
                "raw_low_neff_variants": assessment["raw"][
                    "effective_sample_size_below_minimum_threshold"
                ],
                "raw_low_neff_fraction": assessment["raw"][
                    "effective_sample_size_below_minimum_threshold_fraction"
                ],
                "qc_passed_low_neff_variants": assessment["qc_passed"][
                    "effective_sample_size_below_minimum_threshold"
                ],
                "qc_passed_low_neff_fraction": assessment["qc_passed"][
                    "effective_sample_size_below_minimum_threshold_fraction"
                ],
                "validated_vcf_fields": assessment["vcf_header_validation"][
                    "required_field_count"
                ],
                "aggregation_strategy": assessment["aggregation"][
                    "aggregation_strategy"
                ],
                "temporary_table_scans": assessment["aggregation"][
                    "temporary_table_scans"
                ],
                "streaming_collection_api": assessment["aggregation"][
                    "streaming_collection_api"
                ],
                "postgwas_version": software["postgwas"],
                "python_version": software["python"],
                "polars_version": assessment["aggregation"]["polars_version"],
                "requested_threads": assessment["aggregation"][
                    "requested_threads"
                ],
                "polars_thread_pool_size": assessment["aggregation"][
                    "polars_thread_pool_size"
                ],
                "thread_budget_enforced": assessment["aggregation"][
                    "thread_budget_enforced"
                ],
                "bcftools_version": bcftools_version,
            },
            error_type=VcfAssessmentError,
        )
        logger.record(
            "STATUS",
            "qc_summary_run",
            status="COMPLETED",
            raw=assessment["raw"]["num_records"],
            qc_passed=assessment["qc_passed"]["num_records"],
            excluded=assessment["excluded_total"],
            variants_unique_to_one_rule=assessment["unique_rule_only_total"],
            variants_matching_multiple_rules=assessment["overlap_variants"],
            scientific_provenance_status=assessment["vcf_provenance"]["status"],
            raw_low_neff=assessment["raw"][
                "effective_sample_size_below_minimum_threshold"
            ],
            qc_passed_low_neff=assessment["qc_passed"][
                "effective_sample_size_below_minimum_threshold"
            ],
            reports=len(reports),
        )
        for name, path in {
            **reports,
            "resolved_config": paths["resolved_config"],
            "completion": paths["completion"],
        }.items():
            logger.record("OUTPUT", name, path=path)
        _emit_summary(
            assessment,
            configuration,
            logger,
            include_stage_cards=not progress.enabled,
        )
        if ctx is not None:
            ctx["qc_summary"] = assessment
        return assessment
    except BaseException as exc:
        logger.record(
            "STATUS",
            "qc_summary_run",
            status=(
                "INTERRUPTED" if isinstance(exc, KeyboardInterrupt) else "FAILED"
            ),
            exception_type=type(exc).__name__,
            reason=str(exc) or "no error message was provided",
        )
        logger.error("QC summary failed: %s: %s" % (type(exc).__name__, exc))
        raise
    finally:
        progress.close()
        logger.close()


__all__ = [
    "resolve_qc_summary_configuration",
    "run_qc_assessment",
    "run_qc_summary_direct",
]
