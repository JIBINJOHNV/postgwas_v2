"""Shared entry-VCF preflight callables for pipeline stages."""

from __future__ import annotations

import argparse

from postgwas.core.preflight import (
    PipelinePreflightEvidence,
    pipeline_preflight_evidence,
    require_pipeline_input_vcf,
)


def _vcf_only_preflight(
    module: str,
    preflight_evidence,
    *,
    deferred_checks: tuple[str, ...] = (),
) -> PipelinePreflightEvidence:
    """Reuse the common entry-VCF contract without another file scan."""
    return pipeline_preflight_evidence(
        module,
        preflight_evidence,
        deferred_checks=deferred_checks,
    )


def preflight_sumstat_filter(
    args: argparse.Namespace,
    *,
    preflight_evidence=None,
) -> PipelinePreflightEvidence:
    del args
    return _vcf_only_preflight("sumstat_filter", preflight_evidence)


def preflight_post_imputation_filter(
    args: argparse.Namespace,
    *,
    preflight_evidence=None,
) -> PipelinePreflightEvidence:
    del args
    return _vcf_only_preflight(
        "post_imputation_filter",
        preflight_evidence,
        deferred_checks=(
            "Validate the post-imputation VCF before applying QC filters.",
        ),
    )


def preflight_formatter(
    args: argparse.Namespace,
    *,
    preflight_evidence=None,
) -> PipelinePreflightEvidence:
    del args
    return _vcf_only_preflight(
        "formatter",
        preflight_evidence,
        deferred_checks=(
            "Validate format-specific references selected by downstream modules.",
        ),
    )


def preflight_manhattan(
    args: argparse.Namespace,
    *,
    preflight_evidence=None,
) -> PipelinePreflightEvidence:
    require_pipeline_input_vcf(preflight_evidence)
    from postgwas.modules.manhattan.service import resolve_plot_inputs

    _, indexed, rscript, runtime = resolve_plot_inputs(args)
    return pipeline_preflight_evidence(
        "manhattan",
        preflight_evidence,
        resources={"rscript": rscript, "bcftools": indexed.bcftools, "runtime": runtime},
        deferred_checks=("Validate filtered point data and rendered plot before completion.",),
    )


def preflight_qc_summary(
    args: argparse.Namespace,
    *,
    preflight_evidence=None,
) -> PipelinePreflightEvidence:
    del args
    return _vcf_only_preflight("qc_summary", preflight_evidence)


__all__ = [
    "preflight_formatter",
    "preflight_manhattan",
    "preflight_post_imputation_filter",
    "preflight_qc_summary",
    "preflight_sumstat_filter",
]
