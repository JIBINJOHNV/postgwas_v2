"""Method-aware progress plans for GCTA gene and set pipelines."""

from __future__ import annotations

from typing import Any, Mapping


_PIPELINE_PLANS: dict[str, tuple[tuple[str, str], ...]] = {
    "fastbat_gene": (
        ("vcf", "Validate the input harmonised GWAS-VCF"),
        (
            "reference",
            "Validate the PLINK LD-reference files and determine the BIM "
            "variant-ID format",
        ),
        ("gene_list", "Validate the GCTA gene-coordinate file"),
        ("formatted_input", "Create and validate the BIM-compatible GCTA .ma input"),
        ("run", "Run GCTA fastBAT-gene"),
        ("raw_results", "Validate the raw GCTA results"),
        (
            "corrections",
            "Add and validate nominal, Bonferroni and Benjamini-Hochberg FDR results",
        ),
        (
            "publish",
            "Validate and publish all outputs and display the scientific summary",
        ),
    ),
    "fastbat_segment": (
        ("vcf", "Validate the input harmonised GWAS-VCF"),
        (
            "reference",
            "Validate the PLINK LD-reference files and determine the BIM "
            "variant-ID format",
        ),
        ("formatted_input", "Create and validate the BIM-compatible GCTA .ma input"),
        ("run", "Run GCTA fastBAT-segment"),
        ("raw_results", "Validate the raw GCTA results"),
        (
            "corrections",
            "Add and validate nominal, Bonferroni and Benjamini-Hochberg FDR results",
        ),
        (
            "publish",
            "Validate and publish all outputs and display the scientific summary",
        ),
    ),
    "fastbat_set_prepared": (
        ("vcf", "Validate the input harmonised GWAS-VCF"),
        (
            "reference",
            "Validate the PLINK LD-reference files and determine the BIM "
            "variant-ID format",
        ),
        ("set_list", "Validate the supplied fastBAT set list"),
        ("formatted_input", "Create and validate the BIM-compatible GCTA .ma input"),
        (
            "set_policies",
            "Apply empty-set and set-size policies and write the final fastBAT sets",
        ),
        ("run", "Run GCTA fastBAT-set"),
        ("raw_results", "Validate the raw GCTA results"),
        (
            "corrections",
            "Add and validate nominal, Bonferroni and Benjamini-Hochberg FDR results",
        ),
        (
            "publish",
            "Validate and publish all outputs and display the scientific summary",
        ),
    ),
    "fastbat_set_gmt": (
        ("vcf", "Validate the input harmonised GWAS-VCF"),
        (
            "reference",
            "Validate the PLINK LD-reference files and determine the BIM "
            "variant-ID format",
        ),
        ("gene_list", "Validate the GCTA gene-coordinate file"),
        ("gmt", "Validate the original pathway GMT file"),
        (
            "gene_overlap",
            "Compare pathway genes with the gene-coordinate reference and enforce "
            "the configured minimum coordinate-reference coverage",
        ),
        ("formatted_input", "Create and validate the BIM-compatible GCTA .ma input"),
        ("map_variants", "Map analyzable BIM variants to matched gene intervals"),
        ("candidate_memberships", "Create candidate pathway-to-variant memberships"),
        (
            "set_policies",
            "Apply duplicate-membership, empty-set and set-size policies and write "
            "the final fastBAT sets",
        ),
        ("run", "Run GCTA fastBAT-set"),
        ("raw_results", "Validate the raw GCTA results"),
        (
            "corrections",
            "Add and validate nominal, Bonferroni and Benjamini-Hochberg FDR results",
        ),
        (
            "publish",
            "Validate and publish all outputs and display the scientific summary",
        ),
    ),
    "mbat_combo": (
        ("vcf", "Validate the input harmonised GWAS-VCF"),
        (
            "reference",
            "Validate the PLINK LD-reference files and determine the BIM "
            "variant-ID format",
        ),
        ("gene_list", "Validate the GCTA gene-coordinate file"),
        ("formatted_input", "Create and validate the BIM-compatible GCTA .ma input"),
        ("run", "Run GCTA mBAT-combo"),
        ("raw_results", "Validate the raw GCTA results"),
        (
            "corrections",
            "Add and validate nominal, Bonferroni and Benjamini-Hochberg FDR results",
        ),
        (
            "publish",
            "Validate and publish all outputs and display the scientific summary",
        ),
    ),
}


def _plan_key(module) -> str:
    if module.method == "fastbat_set":
        return (
            "fastbat_set_gmt"
            if module.set_annotation.gmt_file is not None
            else "fastbat_set_prepared"
        )
    return str(module.method)


def build_gcta_gene_pipeline_progress_plan(module) -> dict[str, Any]:
    """Build the one ordered stage plan for the resolved GCTA analysis."""
    plan_key = _plan_key(module)
    entries = _PIPELINE_PLANS[plan_key]
    stage_numbers = {
        key: number for number, (key, _title) in enumerate(entries, 1)
    }
    formatter_end = stage_numbers["formatted_input"]
    return {
        "kind": "gcta_gene",
        "variant": plan_key,
        "label": "GCTA %s pipeline execution progress"
        % str(module.method).replace("_", "-").replace("fastbat", "fastBAT"),
        "stages": tuple(title for _key, title in entries),
        "stage_numbers": stage_numbers,
        "modules": {
            "formatter": (1, formatter_end),
            "gcta_gene": (formatter_end + 1, len(entries)),
        },
    }


def gcta_gene_pipeline_progress_plan(args) -> dict[str, Any]:
    """Resolve and return the method-aware plan for a GCTA-only pipeline."""
    from postgwas.modules.gcta_gene.service import _resolved_configuration

    module = _resolved_configuration(args).modules.gcta_gene
    return build_gcta_gene_pipeline_progress_plan(module)


def pipeline_stage_number(args, key: str) -> int | None:
    """Return one active GCTA pipeline-stage number, or ``None`` in direct mode."""
    plan: Mapping[str, Any] | None = getattr(
        args, "_pipeline_progress_plan", None,
    )
    if not plan or plan.get("kind") != "gcta_gene":
        return None
    numbers = plan.get("stage_numbers", {})
    return int(numbers[key]) if key in numbers else None


def start_pipeline_stage(args, key: str, logger=None) -> int | None:
    """Start one configured stage and mirror its status to the module log."""
    number = pipeline_stage_number(args, key)
    if number is None:
        return None
    plan = args._pipeline_progress_plan
    title = plan["stages"][number - 1]
    args._pipeline_stage_progress.start(number)
    if logger is not None:
        logger.record(
            "STATUS", "gcta_pipeline_stage", step=number,
            total=len(plan["stages"]), stage=title, status="RUNNING",
        )
    return number


def configure_pipeline_stage_callbacks(args, configuration, log_path) -> None:
    """Install durable final-completion and failure logging callbacks."""
    plan: Mapping[str, Any] | None = getattr(
        args, "_pipeline_progress_plan", None,
    )
    if not plan or plan.get("kind") != "gcta_gene":
        return

    from postgwas.core.pipeline_logging import PipelineLogger

    def record_status(status: str, error: BaseException | None = None) -> None:
        progress = args._pipeline_stage_progress
        number = progress.current or min(
            progress.completed + 1, len(plan["stages"]),
        )
        details = {}
        if error is not None:
            details = {
                "error_type": type(error).__name__,
                "error": str(error),
            }
        with PipelineLogger(
            configuration.run.dataset_id,
            "run",
            str(log_path.parent),
            level=configuration.logging.file_level,
            screen_level=configuration.logging.console_level,
            log_path=str(log_path),
        ) as callback_logger:
            callback_logger.record(
                "STATUS",
                "gcta_pipeline_stage",
                step=number,
                total=len(plan["stages"]),
                stage=plan["stages"][number - 1],
                status=status,
                **details,
            )

    args._pipeline_progress_failure = lambda error: record_status(
        "FAILED", error,
    )
    args._pipeline_progress_completion = lambda: record_status("COMPLETED")


def complete_pipeline_stage(
    args,
    key: str,
    *,
    outcome: str | None = None,
    outcome_fields=None,
    screen_outcome_fields=None,
    logger=None,
) -> int | None:
    """Complete one validated stage, retaining full evidence in the log.

    A caller may compact already-reported startup inventory on screen without
    discarding its canonical log fields. The default presentation is unchanged.
    """
    number = pipeline_stage_number(args, key)
    if number is None:
        return None
    plan = args._pipeline_progress_plan
    title = plan["stages"][number - 1]
    display_fields = (
        outcome_fields if screen_outcome_fields is None else screen_outcome_fields
    )
    deferred = (
        number == len(plan["stages"])
        and hasattr(args, "_pipeline_progress_completion")
    )
    if deferred:
        args._pipeline_stage_completion = {
            "outcome": outcome,
            "outcome_fields": display_fields,
        }
    else:
        args._pipeline_stage_progress.complete(
            number, outcome=outcome, outcome_fields=display_fields,
        )
    if logger is not None:
        details = {
            str(field[1]): field[2]
            for field in (outcome_fields or ())
            if len(field) == 3
        }
        logger.record(
            "RESULT", "gcta_pipeline_stage", step=number,
            total=len(plan["stages"]), stage=title,
            outcome=outcome, details=details,
        )
        logger.record(
            "STATUS", "gcta_pipeline_stage", step=number,
            total=len(plan["stages"]), stage=title,
            status=(
                "VALIDATED_PENDING_PIPELINE_CHECKPOINT"
                if deferred else "COMPLETED"
            ),
        )
    return number


__all__ = [
    "build_gcta_gene_pipeline_progress_plan",
    "configure_pipeline_stage_callbacks",
    "gcta_gene_pipeline_progress_plan",
    "pipeline_stage_number",
    "start_pipeline_stage",
    "complete_pipeline_stage",
]
