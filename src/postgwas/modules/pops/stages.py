"""Ordered progress contract for MAGMA followed by PoPS."""

from __future__ import annotations

from typing import Any, Mapping

from postgwas.modules.magma.reporting import (
    MAGMA_GENE_ONLY_STAGE_KEYS,
    MAGMA_STAGE_TITLES,
    magma_pipeline_progress_plan,
)


POPS_STAGES = (
    "Validate the gene-score inputs and optional score companions",
    "Validate the PoPS gene-location annotation",
    "Validate the PoPS feature-matrix files",
    "Validate the PoPS feature-control files",
    "Compare input-score, annotation, and feature gene identifiers",
    "Load the gene scores used to fit PoPS",
    "Adjust the PoPS fitting scores for configured covariates",
    "Test and select predictive features",
    "Fit the PoPS prediction model",
    "Calculate genome-wide PoPS scores",
    "Validate and publish PoPS results",
)

POPS_STAGE_KEYS = (
    "target_inputs",
    "gene_annotation",
    "feature_matrix",
    "feature_controls",
    "gene_compatibility",
    "target_loading",
    "covariate_adjustment",
    "feature_selection",
    "model_fitting",
    "gene_scoring",
    "results",
)


def pops_pipeline_progress_plan(args) -> dict[str, Any] | None:
    """Build the gene-only MAGMA dependency and PoPS stage account."""
    magma_plan = magma_pipeline_progress_plan(args)
    if magma_plan is None:
        return None
    magma_stages = tuple(
        MAGMA_STAGE_TITLES[key] for key in MAGMA_GENE_ONLY_STAGE_KEYS
    )
    magma_stage_numbers = {
        key: number for number, key in enumerate(MAGMA_GENE_ONLY_STAGE_KEYS, 1)
    }
    first = len(magma_stages) + 1
    stage_numbers = {
        key: first + index for index, key in enumerate(POPS_STAGE_KEYS)
    }
    return {
        **magma_plan,
        "kind": "pops",
        "label": "PoPS pipeline execution progress",
        "stages": (*magma_stages, *POPS_STAGES),
        "magma_stage_numbers": magma_stage_numbers,
        "stage_numbers": stage_numbers,
        "modules": {
            "formatter": (
                magma_stage_numbers["vcf"],
                magma_stage_numbers["variant_inputs"],
            ),
            "magma": (
                magma_stage_numbers["variant_inputs"],
                magma_stage_numbers["publish"],
            ),
            "pops": (first, first + len(POPS_STAGES) - 1),
        },
        "deferred_completion_modules": {
            "formatter": magma_stage_numbers["variant_inputs"],
        },
    }


def pops_pipeline_stage_numbers(args) -> Mapping[str, int] | None:
    """Return the active PoPS stage mapping, or ``None`` in direct mode."""
    plan = getattr(args, "_pipeline_progress_plan", None)
    if not isinstance(plan, Mapping) or plan.get("kind") != "pops":
        return None
    return plan["stage_numbers"]


__all__ = [
    "POPS_STAGES",
    "pops_pipeline_progress_plan",
    "pops_pipeline_stage_numbers",
]
