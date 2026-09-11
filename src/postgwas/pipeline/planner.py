"""Construct one deterministic execution plan from requested module targets."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

from postgwas.core.errors import PipelinePlanningError
from postgwas.pipeline.registry import ModuleRegistry, REGISTRY, resolve_reference


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


def resolve_pipeline_dependency_overrides(
    arguments: Sequence[str],
    configuration,
    *,
    registry: ModuleRegistry = REGISTRY,
) -> dict[str, tuple[str, ...]]:
    """Resolve registered method-aware dependencies for planning and export."""
    overrides = {}
    references = {
        spec.pipeline_dependency_override_factory
        for spec in (
            registry.get(name)
            for name in registry.names(include_internal=True)
        )
        if spec.pipeline_dependency_override_factory is not None
    }
    for reference in sorted(references):
        resolved = resolve_reference(reference)(arguments, configuration)
        for module_name, dependencies in resolved.items():
            value = tuple(dependencies)
            previous = overrides.get(module_name)
            if previous is not None and previous != value:
                raise PipelinePlanningError(
                    "Conflicting dynamic dependencies for module '%s': %s "
                    "and %s" % (module_name, previous, value)
                )
            overrides[module_name] = value
    return overrides


def _dependency_order(
    targets: Iterable[str],
    registry: ModuleRegistry,
    dependency_overrides: Mapping[str, Iterable[str]],
) -> list[str]:
    ordered = []
    visiting = set()
    visited = set()

    def visit(name: str) -> None:
        if name in visited:
            return
        if name in visiting:
            raise PipelinePlanningError(
                "Dependency cycle detected at module '%s'" % name
            )
        spec = registry.require_pipeline_enabled(name)
        visiting.add(name)
        dependencies = dependency_overrides.get(name, spec.dependencies)
        for dependency in dependencies:
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
    dependency_overrides: Mapping[str, Iterable[str]] | None = None,
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

    active = _dependency_order(requested, registry, dependency_overrides or {})
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
    active_positions = {name: index for index, name in enumerate(active)}
    ld_clump_before_formatter = (
        "ld_clump" in active_positions
        and "formatter" in active_positions
        and active_positions["ld_clump"] < active_positions["formatter"]
    )
    formatter_before_ld_clump = (
        "ld_clump" in active_positions
        and "formatter" in active_positions
        and active_positions["formatter"] < active_positions["ld_clump"]
    )
    if ld_clump_before_formatter:
        steps.append("ld_clump")

    analysis_format_consumers = tuple(
        name for name in downstream if name != "ld_clump"
    )
    needs_analysis_format = any(
        name in active_set for name in analysis_format_consumers
    ) or formatter_before_ld_clump
    formatter_requested = "formatter" in requested
    if needs_analysis_format or formatter_requested:
        if not imputing or needs_analysis_format or formatter_requested:
            steps.append("formatter")

    for name in downstream:
        if name in active_set and not (
            name == "ld_clump" and ld_clump_before_formatter
        ):
            steps.append(name)

    if "manhattan" in active_set:
        steps.append("manhattan")
    if "qc_summary" in active_set:
        steps.append("qc_summary")

    for step in steps:
        registry.require_pipeline_enabled(step)

    return PipelinePlan(tuple(requested), tuple(active), tuple(steps))
