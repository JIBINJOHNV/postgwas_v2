"""Construct one deterministic execution plan from requested module targets."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from postgwas.core.errors import PipelinePlanningError
from postgwas.pipeline.registry import ModuleRegistry, REGISTRY


@dataclass(frozen=True)
class PipelinePlan:
    requested_modules: tuple[str, ...]
    active_modules: tuple[str, ...]
    steps: tuple[str, ...]

    def __post_init__(self):
        if not self.steps:
            raise PipelinePlanningError("Pipeline plan contains no executable steps")


def _unique(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(value for value in values if value))


def _dependency_order(targets: Iterable[str], registry: ModuleRegistry) -> list[str]:
    ordered = []
    visiting = set()
    visited = set()

    def visit(name: str) -> None:
        if name in visited:
            return
        if name in visiting:
            raise PipelinePlanningError("Dependency cycle detected at module '%s'" % name)
        spec = registry.require_pipeline_enabled(name)
        visiting.add(name)
        for dependency in spec.dependencies:
            visit(dependency)
        visiting.remove(name)
        visited.add(name)
        ordered.append(name)

    for target in targets:
        visit(target)
    return ordered


def build_pipeline_plan(
    modules: Iterable[str] | None,
    *,
    apply_filter: bool = False,
    apply_imputation: bool = False,
    apply_manhattan: bool = False,
    heritability: bool = False,
    registry: ModuleRegistry = REGISTRY,
) -> PipelinePlan:
    """Validate targets and produce the pipeline step order once.

    Duplicate formatter stages are intentional: imputation needs a pre-imputation
    representation, while downstream tools must be formatted from the imputed and
    optionally post-filtered VCF.
    """

    requested = _unique(modules or ())
    if apply_filter and "sumstat_filter" not in requested:
        requested.insert(0, "sumstat_filter")
    if apply_imputation and "imputation" not in requested:
        requested.append("imputation")
    if apply_manhattan and "manhattan" not in requested:
        requested.append("manhattan")
    if heritability and "heritability" not in requested:
        requested.append("heritability")
    if not requested:
        raise PipelinePlanningError("Select at least one pipeline module")

    active = _dependency_order(requested, registry)
    active_set = set(active)
    steps = []

    if "sumstat_filter" in active_set:
        steps.append("sumstat_filter")

    imputing = "imputation" in active_set
    if imputing:
        steps.extend(("formatter", "imputation"))
        if apply_filter or "post_imputation_filter" in active_set:
            steps.append("post_imputation_filter")

    if "annot_ldblock" in active_set:
        steps.append("annot_ldblock")

    downstream = (
        "ld_clump", "magma", "magmacovar", "single_cell", "pops", "kpops", "finemap",
        "caldera", "flames",
        "heritability", "mixer", "gcta_cojo", "gcta_gene",
    )
    needs_analysis_format = any(name in active_set for name in downstream)
    formatter_requested = "formatter" in requested
    if needs_analysis_format or formatter_requested:
        if not imputing or needs_analysis_format or formatter_requested:
            steps.append("formatter")

    for name in downstream:
        if name in active_set:
            steps.append(name)

    if "manhattan" in active_set:
        steps.append("manhattan")
    if "qc_summary" in active_set:
        steps.append("qc_summary")

    for step in steps:
        registry.require_pipeline_enabled(step)

    return PipelinePlan(tuple(requested), tuple(active), tuple(steps))
