"""Validated execution of native scDRS score and downstream analyses."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

from postgwas.core.completion import (
    apply_completion_restart,
    configuration_digest,
    resolve_completion_resume,
    write_completion_manifest,
)
from postgwas.core.io.reports import write_delimited_report, write_yaml_report
from postgwas.core.paths import (
    configured_output_path,
    remove_owned_directory,
    require_nonempty_file,
    resolve_executable,
    validate_filename_component,
)
from postgwas.core.processes import run_checked_command
from postgwas.core.single_cell_validation import (
    H5adValidationPolicy,
    normalize_gene_identifier,
    read_cell_covariates,
    read_gene_identifier_crosswalk,
    read_scdrs_gene_sets,
    validate_h5ad_expression,
)
from postgwas.modules.single_cell.errors import SingleCellError


# scDRS 1.0.3 exposes no CLI seed option and calls score_cell without a seed;
# the supported upstream implementation therefore uses score_cell's fixed default.
_SCDRS_INTERNAL_RANDOM_SEED = 0


@dataclass(frozen=True)
class ScdrsPreflight:
    """Resolved inputs and scientific summaries needed for one scDRS run."""

    h5ad_file: Path
    gene_set_file: Path | None
    magma_gene_results_file: Path | None
    gene_identifier_map_file: Path | None
    covariate_file: Path | None
    executable: str
    software_version: str
    traits: tuple[str, ...]
    h5ad_summary: Mapping[str, Any]
    gene_set_summary: Mapping[str, Any]
    covariate_summary: Mapping[str, Any] | None
    gene_universe: tuple[str, ...]
    cell_names: tuple[str, ...]


@dataclass(frozen=True)
class ScdrsExecution:
    """Native scDRS outputs and execution metrics."""

    outputs: Mapping[str, Path]
    traits: tuple[str, ...]
    metrics: Mapping[str, Any]
    resumed: bool


def validate_scdrs_h5ad(path: str | Path, method) -> dict[str, Any]:
    """Apply the configured scDRS contract through the shared H5AD validator."""
    policy = H5adValidationPolicy(
        matrix_state=method.matrix_state,
        filter_data=method.filter_data,
        minimum_genes_per_cell=method.minimum_genes_per_cell,
        minimum_cells_per_gene=method.minimum_cells_per_gene,
        matrix_chunk_rows=method.validation.matrix_chunk_rows,
        raw_count_integer_tolerance=method.validation.raw_count_integer_tolerance,
        group_analysis=tuple(method.downstream.group_analysis),
        correlation_analysis=tuple(method.downstream.correlation_analysis),
        adjust_proportion_column=method.adjust_proportion_column,
        allow_missing_annotation_values=method.validation.allow_missing_annotation_values,
    )
    return validate_h5ad_expression(path, policy=policy, error_type=SingleCellError)


def validate_scdrs_gene_sets(
    path: str | Path, method, *, gene_universe: tuple[str, ...],
) -> dict[str, Any]:
    """Check effective gene-set overlap after the shared native file parser."""
    parsed = read_scdrs_gene_sets(
        path, **method.gene_set_format.model_dump(), error_type=SingleCellError,
    )
    universe = set(gene_universe)
    records = []
    for row in parsed:
        trait, genes, weighted_genes = row.trait, row.genes, row.weighted_genes
        effective = set(genes).intersection(universe)
        effective_count = len(effective)
        effective_fraction = effective_count / len(universe)
        validation = method.validation
        if effective_count < validation.minimum_effective_genes:
            raise SingleCellError(
                "scDRS trait %s has %d genes in the filtered H5AD universe; "
                "the configured minimum is %d"
                % (
                    trait, effective_count,
                    validation.minimum_effective_genes,
                )
            )
        if effective_fraction >= validation.maximum_effective_gene_fraction:
            raise SingleCellError(
                "scDRS trait %s covers %.6g of the filtered H5AD gene "
                "universe; it must be below the configured maximum %.6g"
                % (
                    trait, effective_fraction,
                    validation.maximum_effective_gene_fraction,
                )
            )
        records.append({
            "trait": trait,
            "input_genes": len(genes),
            "weighted_genes": weighted_genes,
            "effective_genes": effective_count,
            "effective_gene_fraction": effective_fraction,
        })
    return {"traits": records, "trait_count": len(records)}


def validate_scdrs_covariates(
    path: str | Path, method, *, cell_names: tuple[str, ...],
) -> dict[str, Any]:
    """Match the shared parsed covariate evidence to this exact H5AD cell set."""
    configured = method.covariate_format
    evidence = read_cell_covariates(
        path,
        delimiter=configured.delimiter,
        constant_column=configured.constant_column,
        constant_tolerance=configured.constant_tolerance,
        error_type=SingleCellError,
    )
    cells = evidence.cells
    h5ad_cells = set(cell_names)
    missing = h5ad_cells - cells
    unexpected = cells - h5ad_cells
    if configured.require_exact_cell_ids and (missing or unexpected):
        raise SingleCellError(
            "scDRS covariate cell identifiers must exactly match adata.obs_names "
            "(missing=%d, unexpected=%d)" % (len(missing), len(unexpected))
        )
    overlap = len(cells.intersection(h5ad_cells))
    if overlap <= 0:
        raise SingleCellError(
            "scDRS covariates do not overlap adata.obs_names"
        )
    overlap_fraction = overlap / len(h5ad_cells)
    if overlap_fraction <= configured.minimum_cell_overlap_fraction:
        raise SingleCellError(
            "scDRS covariates overlap %.6g of H5AD cells; upstream scDRS "
            "requires overlap above the configured fraction %.6g"
            % (
                overlap_fraction,
                configured.minimum_cell_overlap_fraction,
            )
        )
    return {
        "rows": evidence.rows,
        "covariates": list(evidence.columns),
        "overlapping_cells": overlap,
        "overlap_fraction": overlap_fraction,
        "missing_h5ad_cells": len(missing),
        "unexpected_covariate_cells": len(unexpected),
        "constant_column": configured.constant_column,
        "maximum_constant_deviation": evidence.maximum_constant_deviation,
    }


def preflight_scdrs(
    method,
    executable_value: str | Path,
    *,
    dataset_id: str,
    pipeline_pending_magma: bool = False,
) -> ScdrsPreflight:
    """Resolve exact atlas inputs and either validate or defer the GWAS gene set."""
    h5ad_file = require_nonempty_file(
        method.input.h5ad_file, "scDRS H5AD file", error_type=SingleCellError,
    )
    executable = resolve_executable(
        executable_value, "scDRS executable", error_type=SingleCellError,
    )
    version_output = run_checked_command(
        [executable, *method.version_probe_arguments],
        "scDRS version probe",
        error_type=SingleCellError,
    )
    version_lines = [
        line.strip() for line in version_output.splitlines() if line.strip()
    ]
    if not version_lines:
        raise SingleCellError("scDRS version probe returned no version")
    software_version = version_lines[-1]
    if software_version not in method.supported_versions:
        raise SingleCellError(
            "Unsupported scDRS version %r; supported versions are %s"
            % (software_version, ", ".join(method.supported_versions))
        )
    h5ad_summary = validate_scdrs_h5ad(h5ad_file, method)
    gene_set_file = None
    magma_gene_results_file = None
    gene_identifier_map_file = None
    construction = method.magma_gene_set
    if construction.source == "file":
        gene_set_file = require_nonempty_file(
            method.input.gene_set_file,
            "scDRS .gs gene-set file",
            error_type=SingleCellError,
        )
        gene_set_summary = validate_scdrs_gene_sets(
            gene_set_file,
            method,
            gene_universe=h5ad_summary["gene_universe"],
        )
        traits = tuple(
            record["trait"] for record in gene_set_summary["traits"]
        )
    else:
        trait = validate_filename_component(
            construction.trait_pattern.format(dataset_id=dataset_id),
            "generated scDRS trait",
            error_type=SingleCellError,
        )
        if "@" in trait:
            raise SingleCellError(
                "Generated scDRS trait contains the reserved '@' character"
            )
        traits = (trait,)
        if not pipeline_pending_magma:
            magma_gene_results_file = require_nonempty_file(
                method.input.magma_gene_results_file,
                "MAGMA gene result used to construct scDRS gene sets",
                error_type=SingleCellError,
            )
        if construction.mapping_mode == "crosswalk":
            gene_identifier_map_file = require_nonempty_file(
                method.input.gene_identifier_map_file,
                "scDRS gene-identifier crosswalk",
                error_type=SingleCellError,
            )
            _read_identifier_crosswalk(gene_identifier_map_file, construction)
        gene_set_summary = {
            "source": "magma",
            "status": (
                "awaiting_pipeline_magma"
                if pipeline_pending_magma else "ready_for_construction"
            ),
            "trait_count": 1,
            "traits": [{"trait": trait}],
        }
    covariate_file = None
    covariate_summary = None
    if method.input.covariate_file is not None:
        covariate_file = require_nonempty_file(
            method.input.covariate_file,
            "scDRS covariate file",
            error_type=SingleCellError,
        )
        covariate_summary = validate_scdrs_covariates(
            covariate_file,
            method,
            cell_names=h5ad_summary["cell_names"],
        )
    return ScdrsPreflight(
        h5ad_file=h5ad_file,
        gene_set_file=gene_set_file,
        magma_gene_results_file=magma_gene_results_file,
        gene_identifier_map_file=gene_identifier_map_file,
        covariate_file=covariate_file,
        executable=executable,
        software_version=software_version,
        traits=traits,
        h5ad_summary={
            key: value for key, value in h5ad_summary.items()
            if key not in {"gene_universe", "cell_names"}
        },
        gene_set_summary=gene_set_summary,
        covariate_summary=covariate_summary,
        gene_universe=h5ad_summary["gene_universe"],
        cell_names=h5ad_summary["cell_names"],
    )


def _read_magma_gene_statistics(path: Path, construction) -> pd.DataFrame:
    try:
        frame = pd.read_csv(
            path,
            sep=construction.table_delimiter_pattern,
            comment="#",
            engine="python",
            dtype={construction.gene_id_column: str},
        )
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise SingleCellError(
            "Cannot read MAGMA gene statistics %s: %s" % (path, exc)
        ) from exc
    required = {construction.gene_id_column, construction.z_score_column}
    missing = sorted(required - set(frame.columns))
    if missing or frame.empty:
        raise SingleCellError(
            "MAGMA gene statistics are empty or missing columns: %s"
            % ", ".join(missing or sorted(required))
        )
    genes = frame[construction.gene_id_column].map(normalize_gene_identifier)
    if genes.isna().any():
        raise SingleCellError(
            "MAGMA gene statistics contain missing or empty gene identifiers"
        )
    if genes.duplicated().any():
        raise SingleCellError(
            "MAGMA gene statistics contain duplicate normalized gene identifiers"
        )
    z_scores = pd.to_numeric(frame[construction.z_score_column], errors="coerce")
    if z_scores.isna().any() or not np.isfinite(z_scores.to_numpy(dtype=float)).all():
        raise SingleCellError(
            "MAGMA gene statistics contain missing or non-finite Z scores"
        )
    return pd.DataFrame({"source_gene": genes, "z_score": z_scores.astype(float)})


def _read_identifier_crosswalk(path: Path, construction) -> tuple[dict[str, str], dict]:
    return read_gene_identifier_crosswalk(
        path,
        delimiter=construction.mapping_delimiter,
        source_column=construction.mapping_source_column,
        target_column=construction.mapping_target_column,
        error_type=SingleCellError,
    )


def prepare_scdrs_gene_set(
    preflight: ScdrsPreflight,
    method,
    *,
    statistics_file: Path,
    gene_set_file: Path,
    mapping_report: Path,
    logger,
) -> ScdrsPreflight:
    """Convert validated MAGMA Z statistics to one auditable native scDRS .gs."""
    construction = method.magma_gene_set
    if construction.source != "magma":
        return preflight
    if preflight.magma_gene_results_file is None:
        raise SingleCellError(
            "MAGMA gene results are required to construct the pipeline scDRS gene set"
        )
    statistics = _read_magma_gene_statistics(
        preflight.magma_gene_results_file, construction,
    )
    mapping_metrics = {
        "mapping_mode": construction.mapping_mode,
        "source_identifier_type": construction.source_identifier_type,
        "target_identifier_type": construction.target_identifier_type,
        "input_magma_genes": len(statistics),
    }
    if construction.mapping_mode == "crosswalk":
        if preflight.gene_identifier_map_file is None:
            raise SingleCellError(
                "A gene-identifier crosswalk is required for scDRS mapping_mode "
                "crosswalk"
            )
        mapping, crosswalk_metrics = _read_identifier_crosswalk(
            preflight.gene_identifier_map_file, construction,
        )
        mapping_metrics.update(crosswalk_metrics)
        statistics["target_gene"] = statistics["source_gene"].map(mapping)
        unmapped = int(statistics["target_gene"].isna().sum())
        mapping_metrics["unmapped_magma_genes"] = unmapped
        if unmapped and construction.unmapped_policy == "error":
            raise SingleCellError(
                "%d MAGMA genes are absent from the configured identifier "
                "crosswalk and unmapped_policy is error" % unmapped
            )
        statistics = statistics.loc[
            statistics["target_gene"].notna()
        ].copy()
    else:
        statistics["target_gene"] = statistics["source_gene"]
        mapping_metrics.update({
            "crosswalk_rows": None,
            "unmapped_magma_genes": 0,
        })
    if statistics["target_gene"].duplicated().any():
        raise SingleCellError(
            "MAGMA-to-H5AD identifier conversion produced duplicate target genes"
        )
    if len(statistics) < construction.minimum_genes:
        raise SingleCellError(
            "Only %d MAGMA genes remain after identifier mapping; scDRS "
            "munge-gs requires at least the configured minimum %d"
            % (len(statistics), construction.minimum_genes)
        )
    trait = preflight.traits[0]
    if trait == construction.statistics_gene_column:
        raise SingleCellError(
            "Generated scDRS trait must differ from statistics_gene_column"
        )
    records = [
        {
            construction.statistics_gene_column: row.target_gene,
            trait: float(row.z_score),
        }
        for row in statistics.itertuples(index=False)
    ]
    write_delimited_report(
        records,
        statistics_file,
        fieldnames=[construction.statistics_gene_column, trait],
        delimiter=method.gene_set_format.delimiter,
        null_value=construction.statistics_null_value,
    )
    run_checked_command(
        [
            preflight.executable,
            "munge-gs",
            "--out-file", str(gene_set_file),
            "--zscore-file", str(statistics_file),
            "--weight", construction.weight,
            "--n-min", str(construction.minimum_genes),
            "--n-max", str(construction.maximum_genes),
        ],
        "Construct scDRS gene set from MAGMA Z scores",
        logger=logger,
        error_type=SingleCellError,
        expected_outputs=[gene_set_file],
    )
    gene_set_summary = validate_scdrs_gene_sets(
        gene_set_file,
        method,
        gene_universe=preflight.gene_universe,
    )
    mapping_metrics.update({
        "mapped_magma_genes": len(statistics),
        "mapped_genes_in_h5ad": int(
            statistics["target_gene"].isin(set(preflight.gene_universe)).sum()
        ),
        "unmapped_policy": construction.unmapped_policy,
        "ambiguous_policy": construction.ambiguous_policy,
        "selection_weight": construction.weight,
        "selection_minimum_genes": construction.minimum_genes,
        "selection_maximum_genes": construction.maximum_genes,
        "generated_gene_set": gene_set_summary,
    })
    write_yaml_report(
        {
            "schema_version": 1,
            "status": "COMPLETED",
            "source_gene_results": str(preflight.magma_gene_results_file),
            "crosswalk": (
                None if preflight.gene_identifier_map_file is None
                else str(preflight.gene_identifier_map_file)
            ),
            **mapping_metrics,
        },
        mapping_report,
    )
    logger.record("RESULT", "scdrs_gene_mapping", **mapping_metrics)
    return replace(
        preflight,
        gene_set_file=gene_set_file,
        gene_set_summary=gene_set_summary,
    )


def expected_scdrs_outputs(
    directory: Path,
    traits: tuple[str, ...],
    method,
    output_layout,
) -> dict[str, Path]:
    """Resolve every native output implied by the configured scDRS analyses."""
    outputs: dict[str, Path] = {}
    for trait_index, trait in enumerate(traits, start=1):
        values = {"trait": trait}
        outputs["score_%d" % trait_index] = configured_output_path(
            directory,
            output_layout.scdrs_score_file_pattern,
            error_type=SingleCellError,
            **values,
        )
        outputs["full_score_%d" % trait_index] = configured_output_path(
            directory,
            output_layout.scdrs_full_score_file_pattern,
            error_type=SingleCellError,
            **values,
        )
        for annotation_index, annotation in enumerate(
            method.downstream.group_analysis, start=1,
        ):
            validate_filename_component(
                annotation,
                "scDRS group annotation",
                error_type=SingleCellError,
            )
            outputs[
                "group_%d_%d" % (trait_index, annotation_index)
            ] = configured_output_path(
                directory,
                output_layout.scdrs_group_file_pattern,
                error_type=SingleCellError,
                trait=trait,
                annotation=annotation,
            )
        if method.downstream.correlation_analysis:
            outputs["correlation_%d" % trait_index] = configured_output_path(
                directory,
                output_layout.scdrs_correlation_file_pattern,
                error_type=SingleCellError,
                **values,
            )
        if method.downstream.gene_analysis:
            outputs["gene_%d" % trait_index] = configured_output_path(
                directory,
                output_layout.scdrs_gene_file_pattern,
                error_type=SingleCellError,
                **values,
            )
    if len(outputs) != len(set(outputs.values())):
        raise SingleCellError("Configured scDRS output patterns collide")
    return outputs


def _bool_argument(value: bool) -> str:
    return "True" if value else "False"


def _compute_command(preflight: ScdrsPreflight, method, output: Path) -> list[str]:
    command = [
        preflight.executable,
        "compute-score",
        "--h5ad-file", str(preflight.h5ad_file),
        "--h5ad-species", method.h5ad_species,
        "--gs-file", str(preflight.gene_set_file),
        "--gs-species", method.gene_set_species,
        "--out-folder", str(output),
        "--weight-opt", method.weight_option,
        "--flag-filter-data", _bool_argument(method.filter_data),
        "--flag-raw-count", _bool_argument(method.matrix_state == "raw_counts"),
        "--n-ctrl", str(method.control_gene_sets),
        "--min-genes", str(method.minimum_genes_per_cell),
        "--min-cells", str(method.minimum_cells_per_gene),
        "--flag-return-ctrl-raw-score",
        _bool_argument(method.return_control_raw_score),
        "--flag-return-ctrl-norm-score",
        _bool_argument(method.return_control_normalized_score),
    ]
    if preflight.covariate_file is not None:
        command.extend(["--cov-file", str(preflight.covariate_file)])
    if method.adjust_proportion_column is not None:
        command.extend(["--adj-prop", method.adjust_proportion_column])
    return command


def _downstream_command(
    preflight: ScdrsPreflight,
    method,
    output: Path,
    output_layout,
) -> list[str]:
    full_scores = configured_output_path(
        output,
        output_layout.scdrs_full_score_file_pattern,
        error_type=SingleCellError,
        trait="@",
    )
    command = [
        preflight.executable,
        "perform-downstream",
        "--h5ad-file", str(preflight.h5ad_file),
        "--score-file", str(full_scores),
        "--out-folder", str(output),
        "--flag-filter-data", _bool_argument(method.filter_data),
        "--flag-raw-count", _bool_argument(method.matrix_state == "raw_counts"),
        "--min-genes", str(method.minimum_genes_per_cell),
        "--min-cells", str(method.minimum_cells_per_gene),
        "--knn-n-neighbors", str(method.downstream.knn_neighbors),
        "--knn-n-pcs", str(method.downstream.knn_principal_components),
    ]
    if method.downstream.group_analysis:
        command.extend([
            "--group-analysis", ",".join(method.downstream.group_analysis),
        ])
    if method.downstream.correlation_analysis:
        command.extend([
            "--corr-analysis",
            ",".join(method.downstream.correlation_analysis),
        ])
    if method.downstream.gene_analysis:
        command.append("--gene-analysis")
    return command


def run_scdrs(
    *,
    preflight: ScdrsPreflight,
    configuration,
    paths: Mapping[str, Path],
    dataset: str,
    logger,
) -> ScdrsExecution:
    """Run scDRS atomically, or resume only checksum-matched native outputs."""
    method = configuration.modules.single_cell.scdrs
    layout = configuration.modules.single_cell.output_layout
    engine = paths["scdrs_engine_directory"]
    staging = paths["scdrs_staging_directory"]
    qc_report = paths["scdrs_qc_report"]
    manifest_path = paths["scdrs_completion_manifest"]
    if engine not in qc_report.parents:
        raise SingleCellError(
            "modules.single_cell.output_layout.scdrs_qc_report must be inside "
            "scdrs_engine_directory"
        )
    generated_paths = {
        paths["scdrs_generated_gene_statistics"],
        paths["scdrs_generated_gene_set"],
        paths["scdrs_gene_mapping_report"],
    }
    if any(engine not in path.parents for path in generated_paths):
        raise SingleCellError(
            "Configured generated scDRS gene-set artifacts must be inside "
            "scdrs_engine_directory"
        )
    outputs = expected_scdrs_outputs(
        engine, preflight.traits, method, layout,
    )
    outputs["qc_report"] = qc_report
    construction = method.magma_gene_set
    if construction.source == "magma":
        outputs.update({
            "generated_gene_statistics": paths[
                "scdrs_generated_gene_statistics"
            ],
            "generated_gene_set": paths["scdrs_generated_gene_set"],
            "gene_mapping_report": paths["scdrs_gene_mapping_report"],
        })
        if len(outputs) != len(set(outputs.values())):
            raise SingleCellError("Configured scDRS output paths collide")
        completion_inputs = {
            "h5ad": preflight.h5ad_file,
            "magma_gene_results": preflight.magma_gene_results_file,
        }
        if preflight.gene_identifier_map_file is not None:
            completion_inputs["gene_identifier_map"] = (
                preflight.gene_identifier_map_file
            )
    else:
        completion_inputs = {
            "h5ad": preflight.h5ad_file,
            "gene_sets": preflight.gene_set_file,
        }
    if preflight.covariate_file is not None:
        completion_inputs["covariates"] = preflight.covariate_file
    digest = configuration_digest({
        "scdrs": method.model_dump(mode="json"),
        "magma_gene_set_source": (
            configuration.modules.magma.model_dump(mode="json")
            if construction.source == "magma" else None
        ),
        "output_layout": {
            key: value for key, value in layout.model_dump().items()
            if key.startswith("scdrs_")
        },
        "executable": preflight.executable,
    })
    if (
        configuration.run.resume
        and not configuration.run.overwrite
        and manifest_path.is_file()
    ):
        decision = resolve_completion_resume(
            manifest_path,
            dataset_id=dataset,
            module="single_cell.scdrs",
            genome_build="not_applicable_gene_level_input",
            configuration_sha256=digest,
            inputs=completion_inputs,
            outputs=outputs,
            resume_policy=configuration.run.resume_policy,
            error_type=SingleCellError,
        )
        if decision.action == "resume":
            return ScdrsExecution(
                outputs=outputs,
                traits=preflight.traits,
                metrics=dict(decision.manifest.get("metrics", {})),
                resumed=True,
            )
        apply_completion_restart(
            decision,
            output_root=configuration.run.output_directory,
            manifest=manifest_path,
            logger=logger,
            operation="scdrs_resume",
            error_type=SingleCellError,
        )
        if engine.is_dir() and not any(
            path.is_file() or path.is_symlink()
            for path in engine.rglob("*")
        ):
            remove_owned_directory(
                engine,
                configuration.run.output_directory,
                "empty stale scDRS engine directory",
                error_type=SingleCellError,
            )
    if (
        engine.exists()
        or manifest_path.exists()
    ) and not configuration.run.overwrite:
        raise SingleCellError(
            "Existing or incomplete scDRS output was found; use --resume for "
            "a matching complete run or --overwrite to replace it"
        )
    if staging.exists():
        if not configuration.run.overwrite:
            raise SingleCellError(
                "An isolated incomplete scDRS run exists at %s; review it or "
                "use --overwrite" % staging
            )
        remove_owned_directory(
            staging,
            configuration.run.output_directory,
            "incomplete scDRS staging directory",
            error_type=SingleCellError,
        )
    if configuration.run.overwrite:
        remove_owned_directory(
            engine,
            configuration.run.output_directory,
            "scDRS engine directory",
            error_type=SingleCellError,
        )
        manifest_path.unlink(missing_ok=True)

    staging.mkdir(parents=True, exist_ok=False)
    staged_outputs = expected_scdrs_outputs(
        staging, preflight.traits, method, layout,
    )
    execution_preflight = preflight
    if construction.source == "magma":
        staged_statistics = staging / paths[
            "scdrs_generated_gene_statistics"
        ].relative_to(engine)
        staged_gene_set = staging / paths[
            "scdrs_generated_gene_set"
        ].relative_to(engine)
        staged_mapping_report = staging / paths[
            "scdrs_gene_mapping_report"
        ].relative_to(engine)
        execution_preflight = prepare_scdrs_gene_set(
            preflight,
            method,
            statistics_file=staged_statistics,
            gene_set_file=staged_gene_set,
            mapping_report=staged_mapping_report,
            logger=logger,
        )
        staged_outputs.update({
            "generated_gene_statistics": staged_statistics,
            "generated_gene_set": staged_gene_set,
            "gene_mapping_report": staged_mapping_report,
        })
    run_checked_command(
        _compute_command(execution_preflight, method, staging),
        "scDRS cell-score computation",
        logger=logger,
        error_type=SingleCellError,
        expected_outputs=[
            path for name, path in staged_outputs.items()
            if name.startswith(("score_", "full_score_"))
        ],
    )
    downstream_requested = bool(
        method.downstream.group_analysis
        or method.downstream.correlation_analysis
        or method.downstream.gene_analysis
    )
    if downstream_requested:
        run_checked_command(
            _downstream_command(execution_preflight, method, staging, layout),
            "scDRS downstream analyses",
            logger=logger,
            error_type=SingleCellError,
            expected_outputs=[
                path for name, path in staged_outputs.items()
                if not name.startswith(("score_", "full_score_"))
            ],
        )

    staged_qc = staging / qc_report.relative_to(engine)
    metrics = {
        "traits": len(execution_preflight.traits),
        "trait_ids": list(execution_preflight.traits),
        "cells": preflight.h5ad_summary["cells"],
        "genes": preflight.h5ad_summary["genes"],
        "cells_after_configured_filter": (
            preflight.h5ad_summary["cells_after_configured_filter"]
        ),
        "genes_after_configured_filter": (
            preflight.h5ad_summary["genes_after_configured_filter"]
        ),
        "control_gene_sets": method.control_gene_sets,
        "software_version": preflight.software_version,
        "matrix_source": "X",
        "matrix_state": method.matrix_state,
        "native_output_files": len(
            [
                name for name in staged_outputs
                if name not in {
                    "generated_gene_statistics", "generated_gene_set",
                    "gene_mapping_report",
                }
            ]
        ),
        "gene_set_source": construction.source,
        "upstream_internal_random_seed": _SCDRS_INTERNAL_RANDOM_SEED,
        "postgwas_execution_seed_applied": False,
    }
    magma_provenance = None
    if construction.source == "magma":
        magma = configuration.modules.magma
        primary_mapping = magma.mapping.primary
        definition = magma.mapping.definitions[primary_mapping]
        magma_provenance = {
            "genome_build": magma.genome_build.value,
            "population": magma.population.value,
            "gene_window_upstream_kb": magma.gene_window_upstream_kb,
            "gene_window_downstream_kb": magma.gene_window_downstream_kb,
            "gene_model": magma.gene_model,
            "primary_mapping": primary_mapping,
            "mapping_method": definition.method,
            "gene_id_type": definition.gene_id_type,
            "annotation_source": definition.source_name,
            "annotation_version": definition.source_version,
            "result_statistic_type": definition.result_statistic_type,
        }
    write_yaml_report(
        {
            "schema_version": 1,
            "status": "COMPLETED",
            "dataset_id": dataset,
            "method": "scdrs",
            "software_version": preflight.software_version,
            "inputs": {
                "h5ad": str(preflight.h5ad_file),
                "gene_sets": (
                    str(paths["scdrs_generated_gene_set"])
                    if construction.source == "magma"
                    else str(execution_preflight.gene_set_file)
                ),
                "magma_gene_results": (
                    None if preflight.magma_gene_results_file is None
                    else str(preflight.magma_gene_results_file)
                ),
                "gene_identifier_map": (
                    None if preflight.gene_identifier_map_file is None
                    else str(preflight.gene_identifier_map_file)
                ),
                "covariates": (
                    None if preflight.covariate_file is None
                    else str(preflight.covariate_file)
                ),
            },
            "h5ad_validation": dict(preflight.h5ad_summary),
            "gene_set_validation": dict(execution_preflight.gene_set_summary),
            "covariate_validation": (
                None if preflight.covariate_summary is None
                else dict(preflight.covariate_summary)
            ),
            "magma_provenance": magma_provenance,
            "scientific_settings": {
                "h5ad_species": method.h5ad_species,
                "gene_set_species": method.gene_set_species,
                "gene_set_source": construction.source,
                "matrix_source": "X",
                "matrix_state": method.matrix_state,
                "filter_data": method.filter_data,
                "minimum_genes_per_cell": method.minimum_genes_per_cell,
                "minimum_cells_per_gene": method.minimum_cells_per_gene,
                "control_gene_sets": method.control_gene_sets,
                "upstream_internal_random_seed": (
                    _SCDRS_INTERNAL_RANDOM_SEED
                ),
                "postgwas_execution_seed_applied": False,
                "weight_option": method.weight_option,
                "group_analysis": list(method.downstream.group_analysis),
                "correlation_analysis": list(
                    method.downstream.correlation_analysis
                ),
                "gene_analysis": method.downstream.gene_analysis,
            },
            "metrics": metrics,
        },
        staged_qc,
    )
    engine.parent.mkdir(parents=True, exist_ok=True)
    staging.replace(engine)
    write_completion_manifest(
        manifest_path,
        dataset_id=dataset,
        module="single_cell.scdrs",
        genome_build="not_applicable_gene_level_input",
        configuration_sha256=digest,
        inputs=completion_inputs,
        outputs=outputs,
        metrics=metrics,
        error_type=SingleCellError,
    )
    return ScdrsExecution(
        outputs=outputs,
        traits=execution_preflight.traits,
        metrics=metrics,
        resumed=False,
    )


__all__ = [
    "ScdrsExecution",
    "ScdrsPreflight",
    "expected_scdrs_outputs",
    "preflight_scdrs",
    "run_scdrs",
    "validate_scdrs_covariates",
    "validate_scdrs_gene_sets",
    "validate_scdrs_h5ad",
]
