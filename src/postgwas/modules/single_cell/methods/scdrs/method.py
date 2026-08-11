"""Registered scDRS method adapter for shared single-cell orchestration."""

from __future__ import annotations

import argparse

from postgwas.core.contracts import Artifact, ModuleResult
from postgwas.modules.single_cell.methods.base import MethodRunContext
from postgwas.modules.single_cell.methods.scdrs.runner import (
    ScdrsExecution,
    ScdrsPreflight,
    preflight_scdrs,
    run_scdrs,
)


def _result(
    *,
    configuration,
    dataset: str,
    execution: ScdrsExecution,
) -> ModuleResult:
    method = configuration.modules.single_cell.scdrs
    metadata = {
        "dataset_id": dataset,
        "method": "scdrs",
        "h5ad_species": method.h5ad_species,
        "gene_set_species": method.gene_set_species,
        "matrix_source": "X",
        "matrix_state": method.matrix_state,
        "trait_ids": list(execution.traits),
    }
    artifacts = {
        "scdrs_%s" % name: Artifact(
            "scdrs_native_%s" % name,
            path,
            metadata,
        )
        for name, path in execution.outputs.items()
    }
    return ModuleResult(
        "single_cell",
        artifacts=artifacts,
        metrics={**execution.metrics, "resumed": execution.resumed},
    )


class ScdrsMethod:
    """Native scDRS preflight, execution, and publication adapter."""

    name = "scdrs"
    resource_paths = ("executables.scdrs",)

    def prepare_pipeline_args(self, args: argparse.Namespace) -> None:
        args.scdrs_gene_set_source = "magma"

    def preflight_pipeline(self, args: argparse.Namespace, configuration) -> None:
        preflight_scdrs(
            configuration.modules.single_cell.scdrs,
            configuration.resources.executables.scdrs,
            dataset_id=configuration.run.dataset_id,
            pipeline_pending_magma=True,
        )

    def preflight_direct(
        self,
        args: argparse.Namespace,
        configuration,
        dataset: str,
    ) -> ScdrsPreflight:
        method = configuration.modules.single_cell.scdrs
        preflight = preflight_scdrs(
            method,
            configuration.resources.executables.scdrs,
            dataset_id=dataset,
        )
        method.input.h5ad_file = preflight.h5ad_file
        method.input.gene_set_file = preflight.gene_set_file
        method.input.magma_gene_results_file = preflight.magma_gene_results_file
        method.input.gene_identifier_map_file = (
            preflight.gene_identifier_map_file
        )
        method.input.covariate_file = preflight.covariate_file
        return preflight

    def run(
        self,
        context: MethodRunContext,
        preflight: ScdrsPreflight,
    ) -> ModuleResult:
        execution = run_scdrs(
            preflight=preflight,
            configuration=context.configuration,
            paths=context.paths,
            dataset=context.dataset,
            logger=context.logger,
        )
        context.logger.record(
            "OUTPUT",
            "scdrs_native_outputs",
            directory=str(context.paths["scdrs_engine_directory"]),
            files=len(execution.outputs),
            traits=len(execution.traits),
            software_version=preflight.software_version,
        )
        return _result(
            configuration=context.configuration,
            dataset=context.dataset,
            execution=execution,
        )


METHOD = ScdrsMethod()


__all__ = ["METHOD", "ScdrsMethod"]
