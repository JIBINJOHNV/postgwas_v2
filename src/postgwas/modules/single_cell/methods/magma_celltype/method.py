"""MAGMA cell-typing adapter built on the existing MAGMAcovar service."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from postgwas.core.completion import (
    apply_completion_restart,
    configuration_digest,
    resolve_completion_resume,
    write_completion_manifest,
)
from postgwas.core.contracts import Artifact, ModuleResult
from postgwas.core.paths import (
    remove_empty_directories,
    require_nonempty_file,
    resolve_executable,
)
from postgwas.modules.magmacovar.service import run_magma_covar_direct
from postgwas.modules.single_cell.errors import SingleCellError
from postgwas.modules.single_cell.methods.base import MethodRunContext
from postgwas.modules.single_cell.methods.magma_celltype.analysis import (
    normalize_magma_celltype_results,
    validate_magma_celltype_covariates,
    validate_magma_celltype_use_case,
)


@dataclass(frozen=True)
class MagmaCelltypePreflight:
    """Exact inputs and model definition for one MAGMA cell-typing run."""

    gene_results: Path
    covariates: Path
    use_case: object


def _engine_args(
    args: argparse.Namespace,
    *,
    gene_results: Path,
    covariates: Path,
    engine_directory: Path,
) -> argparse.Namespace:
    delegated = argparse.Namespace(**vars(args))
    delegated.magma_gene_results_file = str(gene_results)
    delegated.covariates = str(covariates)
    delegated.output_directory = str(engine_directory)
    return delegated


def _engine_configuration(
    configuration,
    *,
    gene_results: Path,
    covariates: Path,
    engine_directory: Path,
    use_case,
):
    """Create the exact resolved MAGMAcovar run without reparsing configuration."""
    engine_configuration = configuration.model_copy(deep=True)
    engine_configuration.run.output_directory = engine_directory
    engine = engine_configuration.modules.magmacovar
    engine.model = list(use_case.model)
    engine.direction = use_case.direction
    engine.input.gene_results_file = gene_results
    engine.input.covariates_file = covariates
    return engine_configuration


def _completion_configuration(configuration, engine_configuration) -> dict:
    single_cell = configuration.modules.single_cell
    layout = single_cell.output_layout
    return {
        "magma_celltype": single_cell.magma_celltype.model_dump(mode="json"),
        "multiple_testing": single_cell.multiple_testing.model_dump(mode="json"),
        "result_schema": single_cell.result_schema.model_dump(mode="json"),
        "output_layout": {
            name: getattr(layout, name)
            for name in (
                "engine_directory",
                "results_file",
                "completion_manifest",
                "staging_directory",
            )
        },
        "magma_provenance": configuration.modules.magma.model_dump(mode="json"),
        "effective_magmacovar": (
            engine_configuration.modules.magmacovar.model_dump(mode="json")
        ),
    }


def _result(
    *,
    configuration,
    dataset: str,
    gene_property_results: Path,
    normalized_results: Path,
    metrics: dict,
    resumed: bool,
) -> ModuleResult:
    magma = configuration.modules.magma
    primary = magma.mapping.definitions[magma.mapping.primary]
    metadata = {
        "dataset_id": dataset,
        "genome_build": magma.genome_build.value,
        "population": magma.population.value,
        "gene_id_type": primary.gene_id_type,
        "method": "magma_celltype",
        "workflow": configuration.modules.single_cell.magma_celltype.workflow,
        "schema_version": (
            configuration.modules.single_cell.result_schema.normalized_schema_version
        ),
    }
    return ModuleResult(
        "single_cell",
        artifacts={
            "magma_celltype_results": Artifact(
                "single_cell_cell_type_associations",
                normalized_results,
                metadata,
            ),
            "magma_gene_property_results": Artifact(
                "magma_gene_property_results",
                gene_property_results,
                metadata,
            ),
        },
        metrics={**metrics, "resumed": resumed},
    )


class MagmaCelltypeMethod:
    """Registered MAGMA cell-typing method implementation."""

    name = "magma_celltype"
    resource_paths = ("executables.magma",)

    def prepare_pipeline_args(self, args: argparse.Namespace) -> None:
        return None

    def preflight_pipeline(self, args: argparse.Namespace, configuration) -> dict:
        module = configuration.modules.single_cell
        method = module.magma_celltype
        covariates = require_nonempty_file(
            method.input.covariates_file,
            "MAGMA cell-type covariate matrix",
            error_type=SingleCellError,
        )
        validate_magma_celltype_use_case(
            module,
            configuration.modules.magmacovar,
        )
        validated = validate_magma_celltype_covariates(
            covariates,
            average_property=method.average_property,
            magmacovar_config=configuration.modules.magmacovar,
        )
        resolve_executable(
            configuration.resources.executables.magma,
            "MAGMA executable",
            error_type=SingleCellError,
        )
        return validated

    def preflight_direct(
        self,
        args: argparse.Namespace,
        configuration,
        dataset: str,
    ) -> MagmaCelltypePreflight:
        module = configuration.modules.single_cell
        method = module.magma_celltype
        gene_results = require_nonempty_file(
            method.input.gene_results_file,
            "MAGMA .genes.raw file",
            error_type=SingleCellError,
        )
        covariates = require_nonempty_file(
            method.input.covariates_file,
            "MAGMA cell-type covariate matrix",
            error_type=SingleCellError,
        )
        use_case = validate_magma_celltype_use_case(
            module,
            configuration.modules.magmacovar,
        )
        validate_magma_celltype_covariates(
            covariates,
            average_property=method.average_property,
            magmacovar_config=configuration.modules.magmacovar,
        )
        resolve_executable(
            configuration.resources.executables.magma,
            "MAGMA executable",
            error_type=SingleCellError,
        )
        method.input.gene_results_file = gene_results
        method.input.covariates_file = covariates
        return MagmaCelltypePreflight(
            gene_results=gene_results,
            covariates=covariates,
            use_case=use_case,
        )

    def run(
        self,
        context: MethodRunContext,
        preflight: MagmaCelltypePreflight,
    ) -> ModuleResult:
        """Delegate the published base model to the existing MAGMAcovar module."""
        configuration = context.configuration
        module = configuration.modules.single_cell
        paths = context.paths
        engine_args = _engine_args(
            context.args,
            gene_results=preflight.gene_results,
            covariates=preflight.covariates,
            engine_directory=paths["engine_directory"],
        )
        engine_configuration = _engine_configuration(
            configuration,
            gene_results=preflight.gene_results,
            covariates=preflight.covariates,
            engine_directory=paths["engine_directory"],
            use_case=preflight.use_case,
        )
        magmacovar_result = run_magma_covar_direct(
            engine_args,
            configuration=engine_configuration,
        )
        gene_property_results = Path(
            magmacovar_result["raw_results"]
        ).expanduser().resolve()
        digest = configuration_digest(
            _completion_configuration(configuration, engine_configuration)
        )
        completion_inputs = {
            "magma_gene_results": preflight.gene_results,
            "cell_type_covariates": preflight.covariates,
            "magma_gene_property_results": gene_property_results,
        }
        completion_outputs = {
            "magma_celltype_results": paths["results_file"],
        }
        if (
            configuration.run.resume
            and not configuration.run.overwrite
            and paths["completion_manifest"].is_file()
        ):
            decision = resolve_completion_resume(
                paths["completion_manifest"],
                dataset_id=context.dataset,
                module="single_cell",
                genome_build=configuration.modules.magma.genome_build.value,
                configuration_sha256=digest,
                inputs=completion_inputs,
                outputs=completion_outputs,
                resume_policy=configuration.run.resume_policy,
                error_type=SingleCellError,
            )
            if decision.action == "resume":
                context.logger.record(
                    "SKIP",
                    self.name,
                    reason="provenance_validated_complete_outputs",
                )
                return _result(
                    configuration=configuration,
                    dataset=context.dataset,
                    gene_property_results=gene_property_results,
                    normalized_results=paths["results_file"],
                    metrics=dict(decision.manifest.get("metrics", {})),
                    resumed=True,
                )
            apply_completion_restart(
                decision,
                output_root=context.output,
                manifest=paths["completion_manifest"],
                logger=context.logger,
                operation=self.name,
                error_type=SingleCellError,
            )
        if (
            paths["results_file"].exists()
            or paths["completion_manifest"].exists()
        ) and not configuration.run.overwrite:
            raise SingleCellError(
                "Existing or incomplete MAGMA cell-type output was found; use "
                "--resume for a matching complete run or --overwrite to replace it"
            )

        staging = paths["staging_directory"]
        staged_result = staging / paths["results_file"].relative_to(context.output)
        if staging.exists():
            owned = {staged_result.resolve()}
            unexpected = [
                path
                for path in staging.rglob("*")
                if path.is_file() and path.resolve() not in owned
            ]
            if unexpected:
                raise SingleCellError(
                    "The configured MAGMA cell-type staging directory contains "
                    "unowned files: %s"
                    % ", ".join(str(path) for path in unexpected[:5])
                )
            if not configuration.run.overwrite:
                raise SingleCellError(
                    "An isolated incomplete MAGMA cell-type run exists at %s; "
                    "review it or use --overwrite" % staging
                )
            staged_result.unlink(missing_ok=True)
            remove_empty_directories(staged_result.parent, staging)

        metrics = normalize_magma_celltype_results(
            gene_property_results,
            preflight.covariates,
            staged_result,
            dataset_id=context.dataset,
            single_cell_config=module,
            magmacovar_config=configuration.modules.magmacovar,
        )
        paths["results_file"].parent.mkdir(parents=True, exist_ok=True)
        if configuration.run.overwrite:
            paths["results_file"].unlink(missing_ok=True)
            paths["completion_manifest"].unlink(missing_ok=True)
        staged_result.replace(paths["results_file"])
        remove_empty_directories(staged_result.parent, staging)
        metrics["output"] = str(paths["results_file"])
        write_completion_manifest(
            paths["completion_manifest"],
            dataset_id=context.dataset,
            module="single_cell",
            genome_build=configuration.modules.magma.genome_build.value,
            configuration_sha256=digest,
            inputs=completion_inputs,
            outputs=completion_outputs,
            metrics=metrics,
            error_type=SingleCellError,
        )
        context.logger.record(
            "OUTPUT",
            "magma_celltype_results",
            path=str(paths["results_file"]),
            rows=metrics["tested_cell_types"],
        )
        return _result(
            configuration=configuration,
            dataset=context.dataset,
            gene_property_results=gene_property_results,
            normalized_results=paths["results_file"],
            metrics=metrics,
            resumed=False,
        )


METHOD = MagmaCelltypeMethod()


__all__ = ["METHOD", "MagmaCelltypeMethod", "MagmaCelltypePreflight"]
