"""Configuration, provenance, and execution boundary for GWAS-VCF QC."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from postgwas.config import (
    load_configuration,
    load_run_configuration_for_module,
    write_resolved_configuration,
)
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.config.models.modules.qc_summary import QCSummaryConfig
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
from postgwas.modules.qc_summary.assessment import (
    VcfAssessmentError,
    run_vcf_qc_assessment,
)
from postgwas.modules.qc_summary.reporting import qc_summary_lines


MODULE_CLI_OVERRIDES = {
    "vcf": "inputs.vcf",
    "dataset_id": "inputs.dataset_id",
    "output_directory": "output_directory",
    "genome_build": "target_build",
    "reference_af_column": "reference_af_column",
    "maximum_af_difference": "rules.maximum_af_difference",
}

GLOBAL_CLI_OVERRIDES = {
    "dataset_id": "run.dataset_id",
    "output_directory": "run.output_directory",
    "threads": "execution.threads",
    "memory_gb": "execution.memory_gb",
    "seed": "execution.random_seed",
    "bcftools": "resources.executables.bcftools",
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
    genome_build: str,
    external_af_name: str,
    configuration: QCSummaryConfig,
    bcftools_bin: str,
    logger=None,
) -> dict[str, Any]:
    """Run the shared single-pass assessment used by direct and harmonisation."""
    return run_vcf_qc_assessment(
        vcf_path=vcf_path,
        output_directory=output_directory,
        dataset_id=dataset_id,
        genome_build=genome_build,
        external_af_name=external_af_name,
        configuration=configuration,
        bcftools_bin=bcftools_bin,
        logger=logger,
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


def _emit_summary(assessment: dict[str, Any], configuration, logger) -> None:
    report = "\n".join(qc_summary_lines(assessment))
    logger.log("RESULT", report, wrap=False, screen=False)
    if configuration.logging.show_screen:
        print("\n" + report)


def _fallback_log(args: argparse.Namespace, exc: BaseException) -> None:
    output = getattr(args, "output_directory", None)
    dataset = getattr(args, "dataset_id", None) or "qc_summary"
    if output is None:
        return
    fallback = load_configuration()
    genome_build = (
        getattr(args, "genome_build", None)
        or fallback.modules.qc_summary.target_build.value
    )
    log_path = configured_output_path(
        output,
        fallback.modules.qc_summary.output_layout.service_log_file,
        dataset_id=dataset,
        build=genome_build,
    )
    write_log_record(
        log_path,
        "ERROR",
        "QC summary configuration or input validation failed: %s: %s"
        % (type(exc).__name__, exc),
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
    try:
        configuration = configuration or resolve_qc_summary_configuration(args)
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
        genome_build = module.target_build.value
        bcftools = resolve_executable(
            configuration.resources.executables.bcftools,
            "bcftools executable",
            error_type=VcfAssessmentError,
        )
        paths = _output_paths(output_directory, dataset_id, genome_build, module)
    except BaseException as exc:
        _fallback_log(args, exc)
        raise

    output_directory.mkdir(parents=True, exist_ok=True)
    logger = PipelineLogger(
        dataset_id,
        "run",
        str(paths["log"].parent),
        level=configuration.logging.file_level,
        screen_level=configuration.logging.console_level,
        log_path=str(paths["log"]),
    )
    reports = {name: paths[name] for name in ("summary", "rules", "json")}
    try:
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
        digest = configuration_digest({
            "qc_summary": module.model_dump(mode="json"),
            "bcftools": bcftools,
            "bcftools_version": bcftools_version,
        })
        logger.record(
            "PARAM",
            "qc_summary_run",
            dataset_id=dataset_id,
            genome_build=genome_build,
            vcf=vcf,
            reference_af_column=module.reference_af_column,
            maximum_af_difference=module.rules.maximum_af_difference,
            bcftools=bcftools,
            bcftools_version=bcftools_version,
            resume=configuration.run.resume,
            overwrite=configuration.run.overwrite,
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
            "INPUT", "qc_summary_vcf",
            path=vcf,
            size_bytes=vcf.stat().st_size,
        )
        if (
            configuration.run.resume
            and not configuration.run.overwrite
            and paths["completion"].is_file()
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
                assessment = _load_completed_assessment(paths["json"])
                logger.record(
                    "SKIP",
                    "qc_summary_run",
                    reason="provenance_validated_complete_outputs",
                )
                _emit_summary(assessment, configuration, logger)
                if ctx is not None:
                    ctx["qc_summary"] = assessment
                return assessment
            apply_completion_restart(
                decision,
                output_root=output_directory,
                manifest=paths["completion"],
                logger=logger,
                operation="qc_summary_resume",
                error_type=VcfAssessmentError,
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
            genome_build=genome_build,
            external_af_name=module.reference_af_column,
            configuration=module,
            bcftools_bin=bcftools,
            logger=logger,
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
            },
            error_type=VcfAssessmentError,
        )
        logger.record(
            "STATUS",
            "qc_summary_run",
            status="COMPLETED",
            raw=assessment["raw"]["num_records"],
            qc_passed=assessment["qc_passed"]["num_records"],
            reports=len(reports),
        )
        for name, path in {
            **reports,
            "resolved_config": paths["resolved_config"],
            "completion": paths["completion"],
        }.items():
            logger.record("OUTPUT", name, path=path)
        _emit_summary(assessment, configuration, logger)
        if ctx is not None:
            ctx["qc_summary"] = assessment
        return assessment
    except BaseException as exc:
        logger.error("QC summary failed: %s: %s" % (type(exc).__name__, exc))
        raise
    finally:
        logger.close()


__all__ = [
    "resolve_qc_summary_configuration",
    "run_qc_assessment",
    "run_qc_summary_direct",
]
