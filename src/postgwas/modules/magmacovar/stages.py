"""Ordered progress contract for MAGMA followed by MAGMAcovar."""

from __future__ import annotations

from typing import Any, Mapping

from postgwas.modules.magma.reporting import (
    MAGMA_GENE_ONLY_STAGE_KEYS,
    MAGMA_STAGE_TITLES,
    magma_pipeline_progress_plan,
)


MAGMACOVAR_STAGES = (
    "Validate the MAGMA gene-association results",
    "Validate the gene-property covariates against the MAGMA genes",
    "Run MAGMA gene-property analysis",
    "Validate and adjust the MAGMA gene-property results",
)

MAGMACOVAR_STAGE_KEYS = (
    "gene_results",
    "covariates",
    "analysis",
    "results",
)


def magmacovar_pipeline_progress_plan(args) -> dict[str, Any] | None:
    """Build the gene-only MAGMA dependency and MAGMAcovar stage account."""
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
        key: first + index
        for index, key in enumerate(MAGMACOVAR_STAGE_KEYS)
    }
    return {
        **magma_plan,
        "kind": "magmacovar",
        "label": "MAGMAcovar pipeline execution progress",
        "stages": (*magma_stages, *MAGMACOVAR_STAGES),
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
            "magmacovar": (first, first + len(MAGMACOVAR_STAGES) - 1),
        },
        "deferred_completion_modules": {
            "formatter": magma_stage_numbers["variant_inputs"],
        },
    }


def magmacovar_pipeline_stage_numbers(args) -> Mapping[str, int] | None:
    """Return the active MAGMAcovar stage mapping, or ``None`` in direct mode."""
    plan = getattr(args, "_pipeline_progress_plan", None)
    if not isinstance(plan, Mapping) or plan.get("kind") != "magmacovar":
        return None
    return plan["stage_numbers"]


__all__ = [
    "MAGMACOVAR_STAGES",
    "magmacovar_pipeline_progress_plan",
    "magmacovar_pipeline_stage_numbers",
]
