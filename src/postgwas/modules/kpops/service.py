"""Validated service boundary around the upstream K-POPS v1.0.0 script."""

from __future__ import annotations

import argparse
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

from postgwas.config import (
    load_configuration,
    load_run_configuration_for_module,
    write_resolved_configuration,
)
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.core.completion import (
    configuration_digest,
    validate_completion_manifest,
    write_completion_manifest,
)
from postgwas.core.paths import (
    configured_output_path,
    require_nonempty_file,
    remove_empty_directories,
    resolve_executable,
    validate_filename_component,
)
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.core.processes import run_checked_command
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.modules.kpops.errors import KPopsError


def _resolved_configuration(args: argparse.Namespace):
    module_overrides = explicit_overrides(args, {
        "kpops_genome_build": "genome_build",
        "kpops_script": "script_path",
        "kpops_gene_annotation_file": "gene_annotation_file",
        "kernel_matrix_prefix": "kernel_matrix_prefix",
        "magma_association_prefix": "magma_association_prefix",
        "training_chromosomes": "training_chromosomes",
        "kpops_device": "device",
        "top_contributor_gene_count": "top_contributor_gene_count",
        "anchor_genes": "anchor_genes",
        "anchor_gene_type": "anchor_gene_type",
        "use_magma_covariates": "use_magma_covariates",
        "save_attribution_files": "save_attribution_files",
        "kpops_verbose": "verbose",
    })
    global_overrides = explicit_overrides(args, {
        "dataset_id": "run.dataset_id",
        "output_directory": "run.output_directory",
        "threads": "execution.threads",
        "memory_gb": "execution.memory_gb",
        "seed": "execution.random_seed",
        "resume": "run.resume",
        "overwrite": "run.overwrite",
    })
    return load_run_configuration_for_module(
        "kpops", getattr(args, "run_config", None),
        module_overrides=module_overrides, global_overrides=global_overrides,
    )


def _prefix_file(prefix: str | None, suffix: str, label: str) -> Path:
    if prefix is None:
        raise KPopsError("%s prefix is required" % label)
    return require_nonempty_file(str(Path(prefix).expanduser()) + suffix, label, error_type=KPopsError)


def _read_table(path: Path, delimiter: str, label: str) -> pd.DataFrame:
    try:
        table = pd.read_csv(path, sep=delimiter)
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise KPopsError("Cannot read %s %s: %s" % (label, path, exc)) from exc
    if table.empty:
        raise KPopsError("%s contains no data rows: %s" % (label, path))
    return table


def _require_columns(table: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise KPopsError("%s is missing required columns: %s" % (label, ", ".join(missing)))


def validate_kpops_configuration(args: argparse.Namespace, *, pipeline: bool = False):
    configuration = _resolved_configuration(args)
    module = configuration.modules.kpops
    schema = module.input_schema
    if module.genome_build is None:
        raise KPopsError(
            "Declare modules.kpops.genome_build; K-POPS does not encode or infer a genome build"
        )
    if pipeline and module.genome_build != configuration.modules.magma.genome_build:
        raise KPopsError(
            "K-POPS genome build %s does not match pipeline MAGMA genome build %s"
            % (module.genome_build.value, configuration.modules.magma.genome_build.value)
        )
    script = Path(resolve_executable(module.script_path, "K-POPS", error_type=KPopsError))
    annotation_path = require_nonempty_file(
        module.gene_annotation_file, "K-POPS gene annotation", error_type=KPopsError,
    )
    annotation = _read_table(annotation_path, schema.table_delimiter_pattern, "K-POPS gene annotation")
    annotation_columns = [
        schema.gene_id_column, schema.gene_name_column,
        schema.chromosome_column, schema.tss_column,
    ]
    _require_columns(annotation, annotation_columns, "K-POPS gene annotation")
    gene_ids = annotation[schema.gene_id_column].astype(str)
    gene_names = annotation[schema.gene_name_column].astype(str)
    chromosomes = annotation[schema.chromosome_column].astype(str)
    tss = pd.to_numeric(annotation[schema.tss_column], errors="coerce")
    if gene_ids.duplicated().any():
        raise KPopsError("K-POPS annotation gene IDs must be unique")
    if (
        gene_ids.str.strip().eq("").any() or gene_names.str.strip().eq("").any()
        or chromosomes.str.strip().eq("").any() or tss.isna().any()
        or not np.isfinite(tss.to_numpy()).all() or (tss < 0).any()
    ):
        raise KPopsError("K-POPS annotation contains invalid gene, chromosome, or TSS values")
    explicit_chromosomes = set(module.covariate_projection_chromosomes or ())
    if not {"loco", "all"}.intersection(module.training_chromosomes):
        explicit_chromosomes.update(module.training_chromosomes)
    missing_chromosomes = sorted(explicit_chromosomes - set(chromosomes))
    if missing_chromosomes:
        raise KPopsError(
            "Configured K-POPS chromosomes are absent from the annotation: %s"
            % ", ".join(missing_chromosomes)
        )
    anchor_values = set(gene_ids if module.anchor_gene_type == "ENSGID" else gene_names)
    unknown_anchors = sorted(set(module.anchor_genes) - anchor_values)
    if unknown_anchors:
        raise KPopsError(
            "Configured K-POPS anchor genes are absent from the annotation: %s"
            % ", ".join(unknown_anchors[:10])
        )
    if module.anchor_gene_type == "NAME" and module.anchor_genes:
        name_counts = gene_names.value_counts()
        ambiguous_anchors = sorted(
            name for name in module.anchor_genes if name_counts.get(name, 0) > 1
        )
        if ambiguous_anchors:
            raise KPopsError(
                "Configured K-POPS gene-name anchors map to multiple Ensembl IDs: %s"
                % ", ".join(ambiguous_anchors[:10])
            )
    kernel_genes_path = _prefix_file(
        module.kernel_matrix_prefix, schema.kernel_genes_suffix, "K-POPS kernel genes",
    )
    kernel_path = _prefix_file(
        module.kernel_matrix_prefix, schema.kernel_matrix_suffix, "K-POPS kernel matrix",
    )
    kernel_genes = np.atleast_1d(np.loadtxt(kernel_genes_path, dtype=str)).reshape(-1)
    if len(kernel_genes) < module.minimum_gene_count or len(kernel_genes) != len(set(kernel_genes)):
        raise KPopsError("K-POPS kernel genes must be unique and meet minimum_gene_count")
    expected_bytes = len(kernel_genes) ** 2 * schema.kernel_float_bytes
    if kernel_path.stat().st_size != expected_bytes:
        raise KPopsError(
            "K-POPS kernel has %d bytes; expected %d for %d float32 genes"
            % (kernel_path.stat().st_size, expected_bytes, len(kernel_genes))
        )
    kernel = np.memmap(
        kernel_path, dtype=np.float32, mode="r",
        shape=(len(kernel_genes), len(kernel_genes)), order="C",
    )
    for start in range(0, len(kernel_genes), module.kernel_validation_chunk_rows):
        stop = min(len(kernel_genes), start + module.kernel_validation_chunk_rows)
        block = np.asarray(kernel[start:stop, :])
        if not np.isfinite(block).all():
            raise KPopsError("K-POPS kernel contains non-finite values")
        if not np.allclose(
            block, np.asarray(kernel[:, start:stop]).T,
            rtol=module.kernel_symmetry_relative_tolerance,
            atol=module.kernel_symmetry_absolute_tolerance,
        ):
            raise KPopsError("K-POPS kernel matrix must be symmetric")
    annotation_set = set(gene_ids)
    missing_annotation = set(kernel_genes) - annotation_set
    if missing_annotation:
        raise KPopsError("%d kernel genes are absent from the K-POPS annotation" % len(missing_annotation))
    outcome = None
    if not pipeline:
        outcome = _validate_magma(module, set(kernel_genes), annotation_set)
        outcome_genes = outcome["genes"]
        chromosome_by_gene = dict(zip(gene_ids, chromosomes))
        if module.training_chromosomes == ["loco"]:
            counts = pd.Series([chromosome_by_gene[gene] for gene in outcome_genes]).value_counts()
            minimum_training = len(outcome_genes) - int(counts.max())
        elif module.training_chromosomes == ["all"]:
            minimum_training = len(outcome_genes)
        else:
            selected = set(module.training_chromosomes)
            minimum_training = sum(chromosome_by_gene[gene] in selected for gene in outcome_genes)
        if minimum_training < max(2, module.top_contributor_gene_count):
            raise KPopsError(
                "K-POPS has only %d genes in its smallest training fit, fewer than "
                "the required %d"
                % (minimum_training, max(2, module.top_contributor_gene_count))
            )
    return configuration, {
        "script": script,
        "annotation": annotation_path,
        "annotation_gene_names": dict(zip(gene_ids, gene_names)),
        "annotation_gene_chromosomes": dict(zip(gene_ids, chromosomes)),
        "kernel": kernel_path,
        "kernel_genes": kernel_genes_path,
        "kernel_gene_count": len(kernel_genes),
        "kernel_gene_ids": set(kernel_genes),
        "outcome": outcome,
    }


def _validate_magma(module, kernel_genes: set[str], annotation_genes: set[str]) -> dict:
    schema = module.input_schema
    genes_out = _prefix_file(
        module.magma_association_prefix, schema.magma_genes_out_suffix, "MAGMA gene results",
    )
    genes_raw = _prefix_file(
        module.magma_association_prefix, schema.magma_genes_raw_suffix, "MAGMA raw gene results",
    )
    table = _read_table(genes_out, schema.table_delimiter_pattern, "MAGMA gene results")
    _require_columns(table, [schema.magma_gene_id_column, schema.magma_score_column], "MAGMA gene results")
    genes = table[schema.magma_gene_id_column].astype(str)
    scores = pd.to_numeric(table[schema.magma_score_column], errors="coerce")
    if genes.duplicated().any() or scores.isna().any() or not np.isfinite(scores.to_numpy()).all():
        raise KPopsError("MAGMA gene IDs must be unique and Z statistics must be finite")
    absent_kernel = set(genes) - kernel_genes
    absent_annotation = set(genes) - annotation_genes
    if absent_kernel or absent_annotation:
        raise KPopsError(
            "Every MAGMA gene must occur in the kernel and annotation "
            "(missing kernel=%d, missing annotation=%d)"
            % (len(absent_kernel), len(absent_annotation))
        )
    if len(genes) < module.minimum_gene_count:
        raise KPopsError("MAGMA gene count is below configured minimum_gene_count")
    if module.use_magma_covariates:
        try:
            raw_lines = genes_raw.read_text(encoding="utf-8").splitlines()[2:]
        except OSError as exc:
            raise KPopsError("Cannot read MAGMA raw gene results: %s" % exc) from exc
        rows = [line.split() for line in raw_lines if line.strip()]
        if len(rows) != len(genes) or any(len(row) < 8 for row in rows):
            raise KPopsError(
                "MAGMA raw gene results must contain one upstream covariance row "
                "per MAGMA gene after its two header lines"
            )
        if [row[0] for row in rows] != genes.tolist():
            raise KPopsError("MAGMA .genes.raw and .genes.out gene order must match")
        try:
            covariate_values = np.asarray(
                [[float(row[4]), float(row[5]), float(row[7])] for row in rows],
                dtype=float,
            )
        except ValueError as exc:
            raise KPopsError("MAGMA raw NSNPS, NPARAM, and MAC values must be numeric") from exc
        if not np.isfinite(covariate_values).all() or (covariate_values <= 0).any():
            raise KPopsError("MAGMA raw NSNPS, NPARAM, and MAC values must be finite and positive")
        chromosome_blocks = [row[1] for row in rows]
        completed = set()
        previous = None
        for chromosome in chromosome_blocks:
            if chromosome != previous:
                if chromosome in completed:
                    raise KPopsError("MAGMA raw chromosomes must occur in contiguous blocks")
                if previous is not None:
                    completed.add(previous)
                previous = chromosome
    return {
        "genes_out": genes_out, "genes_raw": genes_raw,
        "gene_count": len(genes), "genes": set(genes),
    }


def preflight_kpops_pipeline(args: argparse.Namespace) -> None:
    validate_kpops_configuration(args, pipeline=True)


def _arguments(configuration, output_prefix: Path, script: Path) -> list[str]:
    module = configuration.modules.kpops
    values = [
        configuration.resources.executables.python, script,
        "--gene_annot_path", module.gene_annotation_file,
        "--kernel_mat_prefix", module.kernel_matrix_prefix,
        "--magma_prefix", module.magma_association_prefix,
        "--out_prefix", str(output_prefix),
        "--random_seed", str(configuration.execution.random_seed),
        "--device", module.device,
        "--top_n_contributor_gene", str(module.top_contributor_gene_count),
        "--anchor_genes_type", module.anchor_gene_type,
        "--training_chromosomes", *module.training_chromosomes,
    ]
    values.append("--use_magma_covariates" if module.use_magma_covariates else "--ignore_magma_covariates")
    values.append(
        "--project_out_covariates_remove_hla"
        if module.remove_hla_during_covariate_projection
        else "--project_out_covariates_keep_hla"
    )
    values.append("--training_remove_hla" if module.remove_hla_during_training else "--training_keep_hla")
    values.append("--testing_remove_hla" if module.remove_hla_during_testing else "--testing_keep_hla")
    values.append("--save_attr_files" if module.save_attribution_files else "--no_save_attr_files")
    values.append("--verbose" if module.verbose else "--no_verbose")
    if module.covariate_projection_chromosomes is not None:
        values.extend(["--project_out_covariates_chromosomes", *module.covariate_projection_chromosomes])
    if module.anchor_genes:
        values.extend(["--anchor_genes", ",".join(module.anchor_genes)])
    return [str(value) for value in values]


def _output_paths(root: Path, dataset: str, module) -> tuple[Path, dict[str, Path]]:
    prefix = configured_output_path(root, module.output_layout.output_prefix, error_type=KPopsError, dataset_id=dataset)
    paths = {
        "predictions": Path(str(prefix) + module.output_layout.predictions_suffix),
        "coefficients": Path(str(prefix) + module.output_layout.coefficients_suffix),
    }
    if module.save_attribution_files:
        paths.update({
            "attribution": Path(str(prefix) + module.output_layout.attribution_suffix),
            "attribution_rows": Path(str(prefix) + module.output_layout.attribution_rows_suffix),
            "attribution_columns": Path(str(prefix) + module.output_layout.attribution_columns_suffix),
        })
    return prefix, paths


def _completion_inputs(resources: dict) -> dict[str, Path]:
    outcome = resources["outcome"]
    return {
        "kpops_script": resources["script"],
        "gene_annotation": resources["annotation"],
        "kernel_matrix": resources["kernel"],
        "kernel_genes": resources["kernel_genes"],
        "magma_genes_out": outcome["genes_out"],
        "magma_genes_raw": outcome["genes_raw"],
    }


def _completion_configuration(configuration) -> str:
    return configuration_digest({
        "module": configuration.modules.kpops.model_dump(
            mode="json", exclude={"reporting"},
        ),
        "random_seed": configuration.execution.random_seed,
        "python": configuration.resources.executables.python,
    })


def _validate_predictions(path: Path, module) -> pd.DataFrame:
    table = _read_table(
        path,
        module.input_schema.predictions_table_delimiter,
        "K-POPS predictions",
    )
    source_id = module.input_schema.predictions_gene_id_column
    published_id = module.input_schema.published_gene_id_column
    if source_id in table.columns and published_id not in table.columns:
        table = table.rename(columns={source_id: published_id})
        table.to_csv(path, sep=module.input_schema.published_table_delimiter, index=False)
    elif (
        published_id not in table.columns
        and module.input_schema.predictions_gene_id_from_index
        and not isinstance(table.index, pd.RangeIndex)
    ):
        table.insert(0, published_id, table.index.astype(str))
        table = table.reset_index(drop=True)
        table.to_csv(path, sep=module.input_schema.published_table_delimiter, index=False)
    required = [published_id, module.input_schema.predictions_score_column]
    _require_columns(table, required, "K-POPS predictions")
    if table[required[0]].astype(str).duplicated().any():
        raise KPopsError("K-POPS predictions contain duplicate gene IDs")
    scores = pd.to_numeric(table[required[1]], errors="coerce")
    if scores.notna().sum() < module.minimum_gene_count or not np.isfinite(scores.dropna()).all():
        raise KPopsError("K-POPS predictions contain too few finite scores")
    return table


def _training_design(module) -> str:
    if module.training_chromosomes == ["loco"]:
        return "leave-one-chromosome-out (LOCO)"
    if module.training_chromosomes == ["all"]:
        return "all target genes"
    return "configured chromosome(s): %s" % ", ".join(module.training_chromosomes)


def _summarise_outputs(path: Path, module, resources: dict) -> dict:
    """Validate K-POPS results and derive scientific reporting metrics."""
    table = _validate_predictions(path, module)
    schema = module.input_schema
    gene_column = schema.published_gene_id_column
    score_column = schema.predictions_score_column
    gene_ids = table[gene_column].astype(str)
    scores = pd.to_numeric(table[score_column], errors="coerce")
    kernel_genes = resources["kernel_gene_ids"]
    unknown_genes = set(gene_ids) - kernel_genes
    if unknown_genes:
        raise KPopsError(
            "%d K-POPS prediction genes are absent from the kernel"
            % len(unknown_genes)
        )

    scored = table.loc[scores.notna()].assign(
        score_value=scores[scores.notna()].to_numpy(),
        gene_value=gene_ids[scores.notna()].to_numpy(),
    )
    ranked = scored.nlargest(module.reporting.top_gene_count, "score_value")
    gene_names = resources["annotation_gene_names"]
    top_genes = [
        {
            "rank": rank,
            "gene_id": row.gene_value,
            "gene_name": gene_names.get(row.gene_value, row.gene_value),
            "score": float(row.score_value),
        }
        for rank, row in enumerate(ranked.itertuples(index=False), 1)
    ]

    scored_genes = set(scored["gene_value"])
    gene_chromosomes = resources["annotation_gene_chromosomes"]
    compatible_chromosomes = {
        gene_chromosomes[gene] for gene in kernel_genes
    }
    scored_chromosomes = {
        gene_chromosomes[gene] for gene in scored_genes
    }
    missing_chromosomes = sorted(compatible_chromosomes - scored_chromosomes)
    compatible_count = len(kernel_genes)
    warnings = []
    if len(scored_genes) < compatible_count:
        warnings.append(
            "K-POPS scores cover %d of %d compatible kernel genes (%.1f%%); "
            "genes without finite scores are excluded from the ranking."
            % (
                len(scored_genes), compatible_count,
                100 * len(scored_genes) / compatible_count,
            )
        )
    if missing_chromosomes:
        warnings.append(
            "No finite K-POPS scores were available on chromosome(s): %s."
            % ", ".join(missing_chromosomes)
        )
    return {
        "genes_in_output": len(table),
        "genes_scored": len(scored_genes),
        "compatible_genes": compatible_count,
        "target_genes": resources["outcome"]["gene_count"],
        "training_design": _training_design(module),
        "top_genes": top_genes,
        "warnings": warnings,
    }


def _render_summary(
    summary: dict, dataset: str, module, paths: dict[str, Path],
    log_path: Path, label_width: int,
) -> str:
    status = (
        "COMPLETED WITH SCIENTIFIC WARNINGS"
        if summary["warnings"] else "COMPLETED"
    )
    lines = [
        "",
        screen_line("analysis", "K-POPS gene-prioritisation summary", indent=2),
        screen_field("info", "Dataset", dataset, indent=6, label_width=label_width),
        screen_field(
            "warning" if summary["warnings"] else "success",
            "Analysis status", status, indent=6, label_width=label_width,
        ),
        "",
        screen_line("genetic", "Scientific findings", indent=6),
        screen_field(
            "count", "Genes represented in output",
            f'{summary["genes_in_output"]:,}', indent=10,
            label_width=label_width,
        ),
        screen_field(
            "count", "Genes with finite K-POPS scores",
            "%s of %s" % (
                f'{summary["genes_scored"]:,}',
                f'{summary["compatible_genes"]:,}',
            ),
            indent=10, label_width=label_width,
        ),
        screen_field(
            "count", "MAGMA target genes", f'{summary["target_genes"]:,}',
            indent=10, label_width=label_width,
        ),
        screen_field(
            "analysis", "Training design", summary["training_design"],
            indent=10, label_width=label_width,
        ),
        "",
        screen_line("decision", "Top prioritized genes", indent=6),
    ]
    precision = module.reporting.score_decimal_places
    for gene in summary["top_genes"]:
        lines.append(screen_field(
            "genetic", "%d. %s" % (gene["rank"], gene["gene_name"]),
            "%s · K-POPS score %.*f" % (
                gene["gene_id"], precision, gene["score"],
            ),
            indent=10, label_width=label_width,
        ))
    lines.extend([
        "",
        screen_line("decision", "How to interpret", indent=6),
        screen_field(
            "info", "K-POPS significance cutoff",
            "None. K-POPS scores are relative rankings, not p-values.",
            indent=10, label_width=label_width,
        ),
        screen_field(
            "info", "Relevant-gene count",
            "Not statistically defined; report a prespecified top-ranked set and "
            "support it with independent genetic or functional evidence.",
            indent=10, label_width=label_width,
        ),
        screen_field(
            "info", "Causal interpretation",
            "A high K-POPS score prioritizes a gene; it does not establish causality.",
            indent=10, label_width=label_width,
        ),
    ])
    if summary["warnings"]:
        lines.extend(["", screen_line("warning", "Scientific warnings", indent=6)])
        lines.extend(
            screen_field(
                "warning", "Warning %d" % number, warning,
                indent=10, label_width=label_width,
            )
            for number, warning in enumerate(summary["warnings"], 1)
        )
    lines.extend([
        "",
        screen_field(
            "success", "Complete K-POPS results", paths["predictions"],
            indent=6, label_width=label_width,
        ),
        screen_field("info", "Full log", log_path, indent=6, label_width=label_width),
        "",
    ])
    return "\n".join(lines)


def _record_summary(logger: PipelineLogger, summary: dict) -> None:
    logger.record("OBSERVED", "kpops_summary", **{
        key: value for key, value in summary.items() if key != "top_genes"
    })
    for warning in summary["warnings"]:
        logger.warning(warning)


def run_kpops_direct(args: argparse.Namespace, ctx=None):
    try:
        configuration, resources = validate_kpops_configuration(args)
    except BaseException as exc:
        fallback = load_configuration()
        output = Path(
            getattr(args, "output_directory", None) or fallback.run.output_directory
        ).expanduser().resolve()
        raw_dataset = getattr(args, "dataset_id", None) or fallback.run.dataset_id
        try:
            dataset = validate_filename_component(raw_dataset, "dataset_id", error_type=KPopsError)
        except KPopsError:
            dataset = fallback.run.dataset_id
        log_path = configured_output_path(
            output, fallback.modules.kpops.output_layout.service_log_file,
            error_type=KPopsError, dataset_id=dataset,
        )
        write_log_record(
            log_path, "ERROR",
            "K-POPS configuration or input validation failed: %s: %s"
            % (type(exc).__name__, exc),
            sample_id=dataset, file_level=fallback.logging.file_level,
            screen_level=fallback.logging.console_level,
        )
        raise
    module = configuration.modules.kpops
    output = Path(configuration.run.output_directory).expanduser().resolve()
    dataset = validate_filename_component(configuration.run.dataset_id, "dataset_id", error_type=KPopsError)
    output.mkdir(parents=True, exist_ok=True)
    prefix, final_paths = _output_paths(output, dataset, module)
    completion = configured_output_path(output, module.output_layout.completion_manifest, error_type=KPopsError, dataset_id=dataset)
    log_path = configured_output_path(output, module.output_layout.service_log_file, error_type=KPopsError, dataset_id=dataset)
    logger = PipelineLogger(
        dataset, "run", str(log_path.parent), level=configuration.logging.file_level,
        screen_level=configuration.logging.console_level, log_path=str(log_path),
    )
    try:
        resolved_path = configured_output_path(output, module.output_layout.resolved_config_file, error_type=KPopsError, dataset_id=dataset)
        write_resolved_configuration(configuration, resolved_path, modules="kpops")
        completion_inputs = _completion_inputs(resources)
        completion_digest = _completion_configuration(configuration)
        if configuration.run.resume and not configuration.run.overwrite and completion.is_file() and all(path.is_file() for path in final_paths.values()):
            validate_completion_manifest(
                completion, dataset_id=dataset, module="kpops",
                genome_build=module.genome_build.value,
                configuration_sha256=completion_digest,
                inputs=completion_inputs, outputs=final_paths,
                error_type=KPopsError,
            )
            summary = _summarise_outputs(
                final_paths["predictions"], module, resources,
            )
            result = {
                "status": "success",
                "kpops_file": str(final_paths["predictions"]),
                "published_files": [str(path) for path in final_paths.values()],
                "completion_manifest": str(completion),
                "summary": summary,
            }
            if ctx is not None:
                ctx["kpops"] = result
            logger.record("SKIP", "kpops_run", reason="validated_complete_outputs")
            _record_summary(logger, summary)
            print(_render_summary(
                summary, dataset, module, final_paths, log_path,
                configuration.logging.terminal_label_width,
            ))
            return result
        existing = [path for path in final_paths.values() if path.exists()]
        if existing and not configuration.run.overwrite:
            raise KPopsError("Existing K-POPS outputs require --overwrite: %s" % ", ".join(map(str, existing)))
        staging_root = configured_output_path(output, module.output_layout.staging_directory, error_type=KPopsError, dataset_id=dataset)
        staging_root.mkdir(parents=True, exist_ok=True)
        try:
            with tempfile.TemporaryDirectory(prefix="run_", dir=staging_root) as directory:
                staged_prefix = Path(directory) / prefix.name
                staged_paths = {
                    name: Path(str(staged_prefix) + str(path).removeprefix(str(prefix)))
                    for name, path in final_paths.items()
                }
                python = resolve_executable(configuration.resources.executables.python, "Python", error_type=KPopsError)
                command = _arguments(configuration, staged_prefix, resources["script"])
                command[0] = python
                run_checked_command(
                    command, "K-POPS", logger=logger, error_type=KPopsError,
                    timeout_seconds=configuration.execution.timeout_seconds,
                    expected_outputs=list(staged_paths.values()),
                )
                summary = _summarise_outputs(
                    staged_paths["predictions"], module, resources,
                )
                if configuration.run.overwrite:
                    completion.unlink(missing_ok=True)
                    for path in final_paths.values():
                        path.unlink(missing_ok=True)
                for name, source in staged_paths.items():
                    final_paths[name].parent.mkdir(parents=True, exist_ok=True)
                    source.replace(final_paths[name])
        finally:
            remove_empty_directories(staging_root, staging_root.parent)
        result = {
            "status": "success", "kpops_file": str(final_paths["predictions"]),
            "published_files": [str(path) for path in final_paths.values()],
            "summary": summary,
        }
        write_completion_manifest(
            completion, dataset_id=dataset, module="kpops",
            genome_build=module.genome_build.value,
            configuration_sha256=completion_digest,
            inputs=completion_inputs, outputs=final_paths,
            metrics={
                "prediction_count": summary["genes_in_output"],
                "finite_score_count": summary["genes_scored"],
                "target_gene_count": summary["target_genes"],
            },
            error_type=KPopsError,
        )
        result["completion_manifest"] = str(completion)
        _record_summary(logger, summary)
        logger.record(
            "STATUS", "kpops_run", status="COMPLETED",
            predictions=summary["genes_in_output"],
            finite_scores=summary["genes_scored"],
        )
        if ctx is not None:
            ctx["kpops"] = result
        print(_render_summary(
            summary, dataset, module, final_paths, log_path,
            configuration.logging.terminal_label_width,
        ))
        return result
    except BaseException as exc:
        logger.error("K-POPS analysis failed: %s: %s" % (type(exc).__name__, exc))
        raise
    finally:
        logger.close()


__all__ = ["preflight_kpops_pipeline", "run_kpops_direct", "validate_kpops_configuration"]
