"""Registered adapter for the official LDSC cell-type-specific workflow."""

from __future__ import annotations

import argparse

from postgwas.core.contracts import Artifact, ModuleResult
from postgwas.modules.single_cell.methods.base import MethodRunContext
from postgwas.modules.single_cell.methods.ldsc_celltype.runner import (
    LdscCelltypeExecution,
    LdscCelltypePreflight,
    preflight_ldsc_celltype,
    run_ldsc_celltype,
)


def _result(
    *,
    configuration,
    dataset: str,
    execution: LdscCelltypeExecution,
) -> ModuleResult:
    method = configuration.modules.single_cell.ldsc_celltype
    metadata = {
        "dataset_id": dataset,
        "method": "ldsc_celltype",
        "workflow": method.workflow,
        "genome_build": method.genome_build.value,
        "population": method.population.value,
        "schema_version": (
            configuration.modules.single_cell.result_schema.normalized_schema_version
        ),
        "coefficient_interpretation": "additional_per_snp_heritability",
        "p_value_alternative": "coefficient_greater_than_zero",
    }
    kinds = {
        "normalized_results": "single_cell_cell_type_associations",
        "native_results": "ldsc_cell_type_results",
        "native_log": "ldsc_native_log",
        "normalized_ldcts": "ldsc_cell_type_manifest",
        "qc_report": "ldsc_cell_type_qc",
        "munged_sumstats": "ldsc_munged_sumstats",
        "munge_log": "ldsc_munge_log",
    }
    return ModuleResult(
        "single_cell",
        artifacts={
            "ldsc_celltype_%s" % name: Artifact(kinds[name], path, metadata)
            for name, path in execution.outputs.items()
        },
        metrics={**execution.metrics, "resumed": execution.resumed},
    )


class LdscCelltypeMethod:
    """Preflight, run, and publish LDSC ``--h2-cts`` analyses."""

    name = "ldsc_celltype"
    resource_paths = ("executables.ldsc", "executables.munge_sumstats")

    def prepare_pipeline_args(self, args: argparse.Namespace) -> None:
        args.ldsc_celltype_sumstats_source = "formatter"

    def preflight_pipeline(self, args: argparse.Namespace, configuration) -> LdscCelltypePreflight:
        return preflight_ldsc_celltype(
            configuration.modules.single_cell.ldsc_celltype,
            configuration.resources.executables.ldsc,
            configuration.resources.executables.munge_sumstats,
            pipeline_pending_sumstats=True,
        )

    def preflight_direct(
        self,
        args: argparse.Namespace,
        configuration,
        dataset: str,
    ) -> LdscCelltypePreflight:
        method = configuration.modules.single_cell.ldsc_celltype
        preflight = preflight_ldsc_celltype(
            method,
            configuration.resources.executables.ldsc,
            configuration.resources.executables.munge_sumstats,
        )
        method.input.sumstats_file = preflight.sumstats_file
        method.input.ldcts_file = preflight.ldcts_file
        method.input.merge_alleles_file = preflight.merge_alleles_file
        method.input.baseline_ld_prefixes = list(preflight.baseline_ld_prefixes)
        method.input.weights_ld_prefix = preflight.weights_ld_prefix
        return preflight

    def run(
        self,
        context: MethodRunContext,
        preflight: LdscCelltypePreflight,
    ) -> ModuleResult:
        execution = run_ldsc_celltype(
            preflight=preflight,
            configuration=context.configuration,
            paths=context.paths,
            dataset=context.dataset,
            logger=context.logger,
        )
        return _result(
            configuration=context.configuration,
            dataset=context.dataset,
            execution=execution,
        )


METHOD = LdscCelltypeMethod()


__all__ = ["LdscCelltypeMethod", "METHOD"]
