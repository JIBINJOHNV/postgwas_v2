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
    apply_completion_restart,
    configuration_digest,
    resolve_completion_resume,
    write_completion_manifest,
)
from postgwas.core.gene_annotation_validation import validate_gene_tss_annotation
from postgwas.core.io.tables import read_pandas_table, require_table_columns
from postgwas.core.io.reports import write_yaml_report
from postgwas.core.matrix_validation import read_unique_names, validate_binary_kernel
from postgwas.core.paths import (
    configured_output_path,
    require_nonempty_file,
    remove_empty_directories,
    resolve_executable,
    validate_filename_component,
)
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.core.preflight import (
    PipelinePreflightEvidence,
    capture_preflight_file_identities,
    pipeline_preflight_evidence,
    require_unchanged_preflight_files,
    require_pipeline_input_vcf,
)
from postgwas.core.processes import run_checked_command
from postgwas.core.required_arguments import (
    RequiredArgument,
    require_resolved_arguments,
)
from postgwas.core.ui import MeasuredProgress, StageProgress
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.modules.kpops.errors import KPopsError


_KPOPS_PROGRESS_STAGES = (
    "Validate MAGMA gene-association input files",
    "Validate the K-POPS gene-annotation resource",
    "Validate the K-POPS kernel and software resources",
    "Resolve the compatible MAGMA gene universe",
    "Fit K-POPS models and calculate gene scores",
    "Validate and publish K-POPS results",
)


def _resolved_configuration(args: argparse.Namespace, *, pipeline: bool = False):
    # Pipeline PoPS options must never overwrite K-POPS's independently
    # configured scientific settings, including for programmatic callers.
    option_prefix = "kpops_" if pipeline else ""
    module_overrides = explicit_overrides(args, {
        "kpops_genome_build": "genome_build",
        "kpops_script": "script_path",
        "kpops_gene_annotation_file": "gene_annotation_file",
        "kernel_matrix_prefix": "kernel_matrix_prefix",
        "magma_association_prefix": "magma_association_prefix",
        option_prefix + "gene_universe_policy": "gene_universe_policy",
        option_prefix + "training_chromosomes": "training_chromosomes",
        "kpops_device": "device",
        "top_contributor_gene_count": "top_contributor_gene_count",
        "anchor_genes": "anchor_genes",
        "anchor_gene_type": "anchor_gene_type",
        option_prefix + "use_magma_covariates": "use_magma_covariates",
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


def validate_kpops_configuration(
    args: argparse.Namespace,
    *,
    pipeline: bool = False,
    configuration=None,
    progress: StageProgress | None = None,
):
    configuration = configuration or _resolved_configuration(args, pipeline=pipeline)
    module = configuration.modules.kpops
    schema = module.input_schema
    requirements = [
        RequiredArgument(
            "--kpops-genome-build",
            "modules.kpops.genome_build",
            module.genome_build,
        ),
        RequiredArgument(
            "--kpops-gene-annotation-file",
            "modules.kpops.gene_annotation_file",
            module.gene_annotation_file,
        ),
        RequiredArgument(
            "--kernel-matrix-prefix",
            "modules.kpops.kernel_matrix_prefix",
            module.kernel_matrix_prefix,
        ),
    ]
    if not pipeline:
        requirements.append(RequiredArgument(
            "--magma-association-prefix",
            "modules.kpops.magma_association_prefix",
            module.magma_association_prefix,
        ))
    require_resolved_arguments(requirements)
    if pipeline and module.genome_build != configuration.modules.magma.genome_build:
        raise KPopsError(
            "K-POPS genome build %s does not match pipeline MAGMA genome build %s"
            % (module.genome_build.value, configuration.modules.magma.genome_build.value)
        )
    outcome = None
    if not pipeline:
        if progress is not None:
            progress.active_stage_number = 1
            progress.start_step(1, len(_KPOPS_PROGRESS_STAGES), _KPOPS_PROGRESS_STAGES[0])
        outcome = _validate_magma_inputs(module)
        if progress is not None:
            progress.complete_step(
                1, len(_KPOPS_PROGRESS_STAGES), _KPOPS_PROGRESS_STAGES[0],
                outcome_fields=[
                    ("info", "MAGMA gene-results file", outcome["genes_out"]),
                    ("info", "MAGMA raw-results file", outcome["genes_raw"]),
                    ("count", "MAGMA genes validated", outcome["original_gene_count"]),
                    ("success", "Gene order and statistics", "valid"),
                ],
            )

    if progress is not None:
        progress.active_stage_number = 2
        progress.start_step(2, len(_KPOPS_PROGRESS_STAGES), _KPOPS_PROGRESS_STAGES[1])
    annotation_path = require_nonempty_file(
        module.gene_annotation_file, "K-POPS gene annotation", error_type=KPopsError,
    )
    annotation = validate_gene_tss_annotation(
        annotation_path, delimiter=schema.table_delimiter_pattern,
        id_column=schema.gene_id_column, name_column=schema.gene_name_column,
        chromosome_column=schema.chromosome_column, tss_column=schema.tss_column,
        require_names=True, require_nonempty_identifiers=True,
        label="K-POPS gene annotation", error_type=KPopsError,
    )
    gene_ids = pd.Series(list(annotation["gene_names"]), dtype=str)
    gene_names = pd.Series(list(annotation["gene_names"].values()), dtype=str)
    chromosomes = pd.Series(list(annotation["gene_chromosomes"].values()), dtype=str)
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
    if module.anchor_gene_type == "ENSGID":
        anchor_gene_ids = set(module.anchor_genes)
    else:
        gene_id_by_name = dict(zip(gene_names, gene_ids))
        anchor_gene_ids = {
            gene_id_by_name[name] for name in module.anchor_genes
        }
    if progress is not None:
        progress.complete_step(
            2, len(_KPOPS_PROGRESS_STAGES), _KPOPS_PROGRESS_STAGES[1],
            outcome_fields=[
                ("info", "Gene-annotation file", annotation_path),
                ("genetic", "Declared genome build", module.genome_build.value),
                ("count", "Annotation genes", len(gene_ids)),
                ("success", "Annotation structure", "valid"),
            ],
        )

    if progress is not None:
        progress.active_stage_number = 3
        progress.start_step(3, len(_KPOPS_PROGRESS_STAGES), _KPOPS_PROGRESS_STAGES[2])
    script = Path(resolve_executable(module.script_path, "K-POPS", error_type=KPopsError))
    python = Path(resolve_executable(
        configuration.resources.executables.python,
        "Python",
        error_type=KPopsError,
    ))
    kernel_genes_path = _prefix_file(
        module.kernel_matrix_prefix, schema.kernel_genes_suffix, "K-POPS kernel genes",
    )
    kernel_path = _prefix_file(
        module.kernel_matrix_prefix, schema.kernel_matrix_suffix, "K-POPS kernel matrix",
    )
    kernel_genes = read_unique_names(kernel_genes_path, "K-POPS kernel genes", error_type=KPopsError)
    if len(kernel_genes) < module.minimum_gene_count:
        raise KPopsError("K-POPS kernel genes must be unique and meet minimum_gene_count")
    # The supported upstream binary kernel format is C-order float32.
    validate_binary_kernel(
        kernel_path, dimension=len(kernel_genes), dtype=np.float32,
        bytes_per_value=schema.kernel_float_bytes,
        chunk_rows=module.kernel_validation_chunk_rows,
        relative_tolerance=module.kernel_symmetry_relative_tolerance,
        absolute_tolerance=module.kernel_symmetry_absolute_tolerance,
        label="K-POPS kernel", error_type=KPopsError,
    )
    annotation_set = set(gene_ids)
    annotation_chromosomes = dict(zip(gene_ids, chromosomes))
    missing_annotation = set(kernel_genes) - annotation_set
    if missing_annotation:
        raise KPopsError("%d kernel genes are absent from the K-POPS annotation" % len(missing_annotation))
    if progress is not None:
        progress.complete_step(
            3, len(_KPOPS_PROGRESS_STAGES), _KPOPS_PROGRESS_STAGES[2],
            outcome_fields=[
                ("info", "K-POPS script", script),
                ("info", "Python executable", python),
                ("info", "Kernel gene-order file", kernel_genes_path),
                ("info", "Kernel matrix file", kernel_path),
                ("count", "Kernel genes", len(kernel_genes)),
                ("success", "Kernel dimensions and symmetry", "valid"),
            ],
        )
    if not pipeline:
        if progress is not None:
            progress.active_stage_number = 4
            progress.start_step(4, len(_KPOPS_PROGRESS_STAGES), _KPOPS_PROGRESS_STAGES[3])
        outcome = _validate_magma_compatibility(
            module, outcome, set(kernel_genes), annotation_chromosomes,
            anchor_gene_ids,
        )
        if progress is not None:
            compatibility = outcome["compatibility"]
            progress.complete_step(
                4, len(_KPOPS_PROGRESS_STAGES), _KPOPS_PROGRESS_STAGES[3],
                outcome_fields=[
                    ("decision", "Gene-universe policy", module.gene_universe_policy),
                    ("success", "Compatible target genes", compatibility["retained_target_genes"]),
                    (
                        "warning" if compatibility["excluded_target_genes"] else "success",
                        "Excluded target genes", compatibility["excluded_target_genes"],
                    ),
                ],
            )
    resources = {
        "script": script,
        "python": python,
        "annotation": annotation_path,
        "annotation_gene_count": len(gene_ids),
        "annotation_gene_names": dict(zip(gene_ids, gene_names)),
        "annotation_gene_chromosomes": dict(zip(gene_ids, chromosomes)),
        "kernel": kernel_path,
        "kernel_genes": kernel_genes_path,
        "kernel_gene_count": len(kernel_genes),
        "kernel_gene_ids": set(kernel_genes),
        "anchor_gene_ids": anchor_gene_ids,
        "outcome": outcome,
    }
    resources["preflight_file_identities"] = capture_preflight_file_identities(
        (
            resources["script"],
            resources["python"],
            resources["annotation"],
            resources["kernel"],
            resources["kernel_genes"],
        ),
        error_type=KPopsError,
        label="K-POPS resource",
    )
    return configuration, resources


def _validate_magma_inputs(module) -> dict:
    schema = module.input_schema
    genes_out = _prefix_file(
        module.magma_association_prefix, schema.magma_genes_out_suffix, "MAGMA gene results",
    )
    genes_raw = _prefix_file(
        module.magma_association_prefix, schema.magma_genes_raw_suffix, "MAGMA raw gene results",
    )
    table = read_pandas_table(genes_out, schema.table_delimiter_pattern, "MAGMA gene results", error_type=KPopsError)
    require_table_columns(table, [schema.magma_gene_id_column, schema.magma_score_column], "MAGMA gene results", error_type=KPopsError)
    genes = table[schema.magma_gene_id_column].astype(str)
    scores = pd.to_numeric(table[schema.magma_score_column], errors="coerce")
    if genes.duplicated().any() or scores.isna().any() or not np.isfinite(scores.to_numpy()).all():
        raise KPopsError("MAGMA gene IDs must be unique and Z statistics must be finite")

    try:
        raw_lines = genes_raw.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise KPopsError("Cannot read MAGMA raw gene results: %s" % exc) from exc
    header_count = schema.magma_raw_header_lines
    headers = raw_lines[:header_count]
    raw_data_lines = [
        line for line in raw_lines[header_count:] if line.strip()
    ]
    rows = [line.split() for line in raw_data_lines]
    raw_indexes = (
        schema.magma_raw_gene_id_index,
        schema.magma_raw_chromosome_index,
        schema.magma_raw_nsnp_index,
        schema.magma_raw_nparam_index,
        schema.magma_raw_mac_index,
    )
    minimum_fields = max(raw_indexes) + 1
    if (
        len(headers) != header_count
        or len(rows) != len(genes)
        or any(len(row) < minimum_fields for row in rows)
    ):
        raise KPopsError(
            "MAGMA .genes.raw must contain the configured header and one "
            "sufficiently populated row per MAGMA gene"
        )
    raw_gene_ids = [row[schema.magma_raw_gene_id_index] for row in rows]
    gene_ids = genes.tolist()
    if raw_gene_ids != gene_ids:
        raise KPopsError("MAGMA .genes.raw and .genes.out gene order must match")

    raw_chromosomes = [row[schema.magma_raw_chromosome_index] for row in rows]
    if module.use_magma_covariates:
        try:
            covariate_values = np.asarray(
                [[
                    float(row[schema.magma_raw_nsnp_index]),
                    float(row[schema.magma_raw_nparam_index]),
                    float(row[schema.magma_raw_mac_index]),
                ] for row in rows],
                dtype=float,
            )
        except ValueError as exc:
            raise KPopsError("MAGMA raw NSNPS, NPARAM, and MAC values must be numeric") from exc
        if not np.isfinite(covariate_values).all() or (covariate_values <= 0).any():
            raise KPopsError("MAGMA raw NSNPS, NPARAM, and MAC values must be finite and positive")
        completed = set()
        previous = None
        for chromosome in raw_chromosomes:
            if chromosome != previous:
                if chromosome in completed:
                    raise KPopsError("MAGMA raw chromosomes must occur in contiguous blocks")
                if previous is not None:
                    completed.add(previous)
                previous = chromosome
    return {
        "genes_out": genes_out, "genes_raw": genes_raw,
        "original_gene_count": len(gene_ids), "original_gene_ids": gene_ids,
        "table": table, "raw_headers": headers,
        "raw_data_lines": raw_data_lines, "raw_chromosomes": raw_chromosomes,
    }


def _validate_magma_compatibility(
    module,
    outcome: dict,
    kernel_genes: set[str],
    annotation_chromosomes: dict[str, str],
    anchor_gene_ids: set[str],
) -> dict:
    gene_ids = outcome["original_gene_ids"]
    raw_chromosomes = outcome["raw_chromosomes"]
    annotation_genes = set(annotation_chromosomes)
    absent_kernel = sorted(set(gene_ids) - kernel_genes)
    absent_annotation = sorted(set(gene_ids) - annotation_genes)
    shared = kernel_genes & annotation_genes
    retained_gene_ids = [gene for gene in gene_ids if gene in shared]
    retained = set(retained_gene_ids)
    excluded_gene_ids = [gene for gene in gene_ids if gene not in retained]
    absent_kernel_set = set(absent_kernel)
    absent_annotation_set = set(absent_annotation)
    absent_from_both = absent_kernel_set & absent_annotation_set
    compatibility = {
        "policy": module.gene_universe_policy,
        "original_target_genes": len(gene_ids),
        "retained_target_genes": len(retained_gene_ids),
        "excluded_target_genes": len(excluded_gene_ids),
        "retained_percent": 100 * len(retained_gene_ids) / len(gene_ids),
        "excluded_percent": 100 * len(excluded_gene_ids) / len(gene_ids),
        "absent_from_kernel": len(absent_kernel),
        "absent_from_annotation": len(absent_annotation),
        "absent_from_both": len(absent_from_both),
        "absent_only_from_kernel": len(absent_kernel_set - absent_annotation_set),
        "absent_only_from_annotation": len(
            absent_annotation_set - absent_kernel_set
        ),
        "retained_gene_ids": retained_gene_ids,
        "excluded_gene_ids": excluded_gene_ids,
        "missing_kernel_gene_ids": absent_kernel,
        "missing_annotation_gene_ids": absent_annotation,
    }
    chromosome_mismatches = [
        (gene, chromosome, annotation_chromosomes[gene])
        for gene, chromosome in zip(gene_ids, raw_chromosomes)
        if gene in retained and chromosome != annotation_chromosomes[gene]
    ]
    if chromosome_mismatches:
        examples = ", ".join(
            "%s (MAGMA=%s, K-POPS=%s)" % values
            for values in chromosome_mismatches[:module.reporting.top_gene_count]
        )
        raise KPopsError(
            "MAGMA raw metadata and the K-POPS annotation disagree on "
            "chromosome for %d retained genes. Examples: %s. Gene-universe "
            "intersection cannot repair a chromosome disagreement."
            % (len(chromosome_mismatches), examples)
        )
    if excluded_gene_ids and module.gene_universe_policy == "strict":
        raise KPopsError(
            "K-POPS input gene universes are incompatible: MAGMA genes=%d, "
            "shared genes=%d, missing from kernel=%d, missing from annotation=%d. "
            "Example excluded genes: %s. Supply matching resources or set "
            "modules.kpops.gene_universe_policy to intersect in the run "
            "configuration to retain and audit the shared genes."
            % (
                len(gene_ids), len(retained_gene_ids), len(absent_kernel),
                len(absent_annotation),
                ", ".join(excluded_gene_ids[:module.reporting.top_gene_count]),
            )
        )
    if len(retained_gene_ids) < module.minimum_gene_count:
        raise KPopsError(
            "K-POPS retains only %d shared MAGMA genes, below the configured "
            "minimum_gene_count of %d"
            % (len(retained_gene_ids), module.minimum_gene_count)
        )
    missing_target_anchors = sorted(anchor_gene_ids - retained)
    if missing_target_anchors:
        raise KPopsError(
            "Configured K-POPS anchor genes are absent from the retained "
            "MAGMA target universe: %s"
            % ", ".join(missing_target_anchors[:10])
        )
    if module.training_chromosomes == ["loco"]:
        counts = pd.Series(
            [annotation_chromosomes[gene] for gene in retained]
        ).value_counts()
        minimum_training = len(retained) - int(counts.max())
    elif module.training_chromosomes == ["all"]:
        minimum_training = len(retained)
    else:
        selected = set(module.training_chromosomes)
        minimum_training = sum(
            annotation_chromosomes[gene] in selected for gene in retained
        )
    if minimum_training < max(2, module.top_contributor_gene_count):
        raise KPopsError(
            "K-POPS has only %d genes in its smallest training fit, fewer than "
            "the required %d"
            % (minimum_training, max(2, module.top_contributor_gene_count))
        )
    outcome.update({
        "gene_count": len(retained_gene_ids), "genes": retained,
        "compatibility": compatibility,
    })
    return outcome


def _prepare_intersected_magma(
    module,
    outcome: dict,
    paths: dict[str, Path],
    published_paths: dict[str, Path],
    derived_prefix: Path,
) -> tuple[str, dict]:
    """Create temporary aligned K-POPS inputs and a durable exclusion audit."""
    compatibility = outcome["compatibility"]
    retained = set(compatibility["retained_gene_ids"])
    missing_kernel = set(compatibility["missing_kernel_gene_ids"])
    missing_annotation = set(compatibility["missing_annotation_gene_ids"])
    audit_rows = []
    chromosome_statistics = {}
    for row_number, (gene, chromosome) in enumerate(
        zip(outcome["original_gene_ids"], outcome["raw_chromosomes"]), 1,
    ):
        in_kernel = gene not in missing_kernel
        in_annotation = gene not in missing_annotation
        is_retained = gene in retained
        if is_retained:
            decision = "retained"
        elif not in_kernel and not in_annotation:
            decision = "missing_kernel_and_annotation"
        elif not in_kernel:
            decision = "missing_kernel"
        else:
            decision = "missing_annotation"
        audit_rows.append({
            "magma_row": row_number,
            "gene_id": gene,
            "chromosome": chromosome,
            "present_in_kpops_annotation": in_annotation,
            "present_in_kernel": in_kernel,
            "retained_for_kpops": is_retained,
            "decision": decision,
        })
        counts = chromosome_statistics.setdefault(
            chromosome, {"original": 0, "retained": 0, "excluded": 0},
        )
        counts["original"] += 1
        counts["retained" if is_retained else "excluded"] += 1

    pd.DataFrame(audit_rows).to_csv(
        paths["gene_compatibility_table"],
        sep=module.input_schema.published_table_delimiter,
        index=False,
    )
    excluded = compatibility["excluded_target_genes"]
    report = {
        key: value
        for key, value in compatibility.items()
        if not key.endswith("_gene_ids")
    }
    report.update({
        "decision": (
            "MAGMA targets were restricted to genes present in both the K-POPS "
            "annotation and kernel."
            if excluded
            else "All MAGMA targets were already compatible; no restriction was required."
        ),
        "analysis_effect": (
            "K-POPS fitting uses only retained MAGMA target genes; rankings may "
            "differ from a run built from one fully matched gene annotation release."
        ),
        "original_files_unchanged": True,
        "derived_magma_inputs": (
            "Temporary aligned .genes.out and .genes.raw inputs are removed after "
            "the K-POPS run and must not be reused by another analysis."
        ),
        "excluded_gene_examples": compatibility["excluded_gene_ids"][
            :module.reporting.top_gene_count
        ],
        "exclusion_reason_counts": {
            reason: sum(row["decision"] == reason for row in audit_rows)
            for reason in (
                "retained", "missing_kernel", "missing_annotation",
                "missing_kernel_and_annotation",
            )
        },
        "chromosomes": chromosome_statistics,
        "files": {
            "original_genes_out": str(outcome["genes_out"]),
            "original_genes_raw": str(outcome["genes_raw"]),
            "gene_audit_table": str(published_paths["gene_compatibility_table"]),
            "compatibility_report": str(
                published_paths["gene_compatibility_report"]
            ),
        },
    })
    write_yaml_report(report, paths["gene_compatibility_report"])
    if not excluded:
        return str(module.magma_association_prefix), report

    schema = module.input_schema
    genes_out = Path(str(derived_prefix) + schema.magma_genes_out_suffix)
    genes_raw = Path(str(derived_prefix) + schema.magma_genes_raw_suffix)
    output_ids = outcome["table"][schema.magma_gene_id_column].astype(str)
    outcome["table"].loc[output_ids.isin(retained)].to_csv(
        genes_out, sep=schema.published_table_delimiter, index=False,
    )
    selected_raw_lines = [
        line
        for gene, line in zip(
            outcome["original_gene_ids"], outcome["raw_data_lines"],
        )
        if gene in retained
    ]
    genes_raw.write_text(
        "\n".join(outcome["raw_headers"] + selected_raw_lines) + "\n",
        encoding="utf-8",
    )
    derived = read_pandas_table(
        genes_out, schema.table_delimiter_pattern,
        "derived compatible MAGMA gene results", error_type=KPopsError)
    derived_ids = derived[schema.magma_gene_id_column].astype(str).tolist()
    if derived_ids != compatibility["retained_gene_ids"]:
        raise KPopsError(
            "Derived compatible MAGMA gene results failed order validation"
        )
    derived_raw_ids = [
        line.split()[schema.magma_raw_gene_id_index]
        for line in selected_raw_lines
    ]
    if derived_raw_ids != compatibility["retained_gene_ids"]:
        raise KPopsError(
            "Derived compatible MAGMA raw results failed order validation"
        )
    return str(derived_prefix), report


def preflight_kpops_pipeline(
    args: argparse.Namespace,
    *,
    preflight_evidence=None,
) -> PipelinePreflightEvidence:
    require_pipeline_input_vcf(preflight_evidence)
    resources = validate_kpops_configuration(args, pipeline=True)
    return pipeline_preflight_evidence(
        "kpops",
        preflight_evidence,
        resources=resources,
        deferred_checks=(
            "Validate the pipeline-generated MAGMA gene-association result.",
        ),
    )


def _reuse_kpops_pipeline_resources(
    configuration,
    pipeline_resources,
    progress: StageProgress,
):
    """Reuse unchanged external resources and validate the generated MAGMA input."""
    if not (
        isinstance(pipeline_resources, tuple)
        and len(pipeline_resources) == 2
        and isinstance(pipeline_resources[1], dict)
    ):
        raise KPopsError(
            "K-POPS pipeline execution requires completed resource-preflight evidence"
        )
    preflight_configuration, validated = pipeline_resources
    current_module = configuration.modules.kpops
    preflight_module = preflight_configuration.modules.kpops
    excluded = {"magma_association_prefix"}
    if current_module.model_dump(mode="json", exclude=excluded) != preflight_module.model_dump(
        mode="json", exclude=excluded,
    ):
        raise KPopsError(
            "K-POPS configuration changed after pipeline preflight; restart so "
            "resources are validated against the executed configuration"
        )
    if (
        configuration.resources.executables.python
        != preflight_configuration.resources.executables.python
    ):
        raise KPopsError(
            "K-POPS Python configuration changed after pipeline preflight; restart "
            "so the executable is validated before analysis"
        )
    identities = validated.get("preflight_file_identities")
    if not identities:
        raise KPopsError("K-POPS resource-preflight evidence is incomplete")
    resources = dict(validated)
    progress.active_stage_number = 1
    progress.start_step(1, len(_KPOPS_PROGRESS_STAGES), _KPOPS_PROGRESS_STAGES[0])
    outcome = _validate_magma_inputs(current_module)
    progress.complete_step(
        1,
        len(_KPOPS_PROGRESS_STAGES),
        _KPOPS_PROGRESS_STAGES[0],
        outcome_fields=[
            ("info", "MAGMA gene-results file", outcome["genes_out"]),
            ("info", "MAGMA raw-results file", outcome["genes_raw"]),
            ("count", "MAGMA genes validated", outcome["original_gene_count"]),
            ("success", "Gene order and statistics", "valid"),
        ],
    )
    annotation_path = Path(resources["annotation"]).resolve()
    validation_stages = (
        (2, tuple(
            identity for identity in identities
            if identity.path == annotation_path
        ), [
            ("info", "Gene-annotation file", resources["annotation"]),
            ("count", "Annotation genes", resources["annotation_gene_count"]),
            ("decision", "Validation evidence", "reused; file identity unchanged"),
        ]),
        (3, tuple(
            identity for identity in identities
            if identity.path != annotation_path
        ), [
            ("info", "Kernel matrix file", resources["kernel"]),
            ("count", "Kernel genes", resources["kernel_gene_count"]),
            ("decision", "Validation evidence", "reused; file identities unchanged"),
        ]),
    )
    if any(not stage_identities for _, stage_identities, _ in validation_stages):
        raise KPopsError("K-POPS resource-preflight evidence is incomplete")
    for stage_number, stage_identities, fields in validation_stages:
        progress.active_stage_number = stage_number
        progress.start_step(
            stage_number,
            len(_KPOPS_PROGRESS_STAGES),
            _KPOPS_PROGRESS_STAGES[stage_number - 1],
        )
        require_unchanged_preflight_files(
            stage_identities,
            error_type=KPopsError,
            label=(
                "K-POPS gene-annotation resource"
                if stage_number == 2
                else "K-POPS kernel or software resource"
            ),
        )
        progress.complete_step(
            stage_number,
            len(_KPOPS_PROGRESS_STAGES),
            _KPOPS_PROGRESS_STAGES[stage_number - 1],
            outcome_fields=fields,
        )

    progress.active_stage_number = 4
    progress.start_step(4, len(_KPOPS_PROGRESS_STAGES), _KPOPS_PROGRESS_STAGES[3])
    outcome = _validate_magma_compatibility(
        current_module,
        outcome,
        resources["kernel_gene_ids"],
        resources["annotation_gene_chromosomes"],
        resources["anchor_gene_ids"],
    )
    compatibility = outcome["compatibility"]
    progress.complete_step(
        4,
        len(_KPOPS_PROGRESS_STAGES),
        _KPOPS_PROGRESS_STAGES[3],
        outcome_fields=[
            ("decision", "Gene-universe policy", current_module.gene_universe_policy),
            ("success", "Compatible target genes", compatibility["retained_target_genes"]),
            (
                "warning" if compatibility["excluded_target_genes"] else "success",
                "Excluded target genes",
                compatibility["excluded_target_genes"],
            ),
        ],
    )
    resources["outcome"] = outcome
    return configuration, resources


def _arguments(
    configuration,
    output_prefix: Path,
    script: Path,
    *,
    magma_prefix: str | None = None,
    emit_model_progress: bool = False,
) -> list[str]:
    module = configuration.modules.kpops
    values = [
        configuration.resources.executables.python, script,
        "--gene_annot_path", module.gene_annotation_file,
        "--kernel_mat_prefix", module.kernel_matrix_prefix,
        "--magma_prefix", magma_prefix or module.magma_association_prefix,
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
    values.append("--verbose" if module.verbose or emit_model_progress else "--no_verbose")
    if module.covariate_projection_chromosomes is not None:
        values.extend(["--project_out_covariates_chromosomes", *module.covariate_projection_chromosomes])
    if module.anchor_genes:
        values.extend(["--anchor_genes", ",".join(module.anchor_genes)])
    return [str(value) for value in values]


def _model_fit_count(module, resources: dict) -> int:
    """Return the exact number of model fits performed by upstream K-POPS."""
    if module.training_chromosomes == ["loco"]:
        return len(set(resources["annotation_gene_chromosomes"].values()))
    return 1


def _model_progress_observer(
    native_log: Path,
    measured: MeasuredProgress,
    total: int,
    logger: PipelineLogger,
):
    """Translate observed upstream model-fit messages into measured progress."""
    state = {"offset": 0, "pending": "", "completed": 0}
    completion_message = "Computing PoPS scores."

    def observe() -> None:
        try:
            with native_log.open("r", encoding="utf-8", errors="replace") as handle:
                handle.seek(state["offset"])
                chunk = handle.read()
                state["offset"] = handle.tell()
        except FileNotFoundError:
            return
        if not chunk:
            return
        text = state["pending"] + chunk
        lines = text.splitlines(keepends=True)
        state["pending"] = ""
        if lines and not lines[-1].endswith(("\n", "\r")):
            state["pending"] = lines.pop()
        observed = sum(completion_message in line for line in lines)
        newly_completed = min(observed, total - state["completed"])
        for _ in range(newly_completed):
            state["completed"] += 1
            measured.update(
                state["completed"],
                total=total,
                title="Fit K-POPS chromosome models",
            )
            logger.record(
                "OBSERVED",
                "kpops_model_fit_progress",
                completed_models=state["completed"],
                total_models=total,
                upstream_event=completion_message,
            )

    return observe


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
    if module.gene_universe_policy == "intersect":
        paths.update({
            "gene_compatibility_table": Path(
                str(prefix) + module.output_layout.gene_compatibility_table_suffix
            ),
            "gene_compatibility_report": Path(
                str(prefix) + module.output_layout.gene_compatibility_report_suffix
            ),
        })
    return prefix, paths


def _completion_inputs(resources: dict) -> dict[str, Path]:
    outcome = resources["outcome"]
    return {
        "kpops_script": resources["script"],
        "python_executable": resources["python"],
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
    table = read_pandas_table(
        path,
        module.input_schema.predictions_table_delimiter,
        "K-POPS predictions", error_type=KPopsError)
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
    require_table_columns(table, required, "K-POPS predictions", error_type=KPopsError)
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


def _summarise_outputs(
    path: Path,
    module,
    resources: dict,
    *,
    paths: dict[str, Path] | None = None,
) -> dict:
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
    compatibility = {
        key: value
        for key, value in resources["outcome"]["compatibility"].items()
        if not key.endswith("_gene_ids")
    }
    if compatibility["excluded_target_genes"]:
        warnings.append(
            "Gene-universe intersection excluded %d of %d MAGMA target genes "
            "and retained %d (%.1f%%). K-POPS fitting used only the retained "
            "genes; review the compatibility audit before interpreting rankings."
            % (
                compatibility["excluded_target_genes"],
                compatibility["original_target_genes"],
                compatibility["retained_target_genes"],
                compatibility["retained_percent"],
            )
        )
    if module.gene_universe_policy == "intersect" and paths is not None:
        compatibility["audit_table"] = str(paths["gene_compatibility_table"])
        compatibility["report"] = str(paths["gene_compatibility_report"])
    return {
        "genes_in_output": len(table),
        "genes_scored": len(scored_genes),
        "compatible_genes": compatible_count,
        "target_genes": resources["outcome"]["gene_count"],
        "original_target_genes": resources["outcome"]["original_gene_count"],
        "excluded_target_genes": compatibility["excluded_target_genes"],
        "gene_compatibility": compatibility,
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
            "count", "Original MAGMA target genes",
            f'{summary["original_target_genes"]:,}',
            indent=10, label_width=label_width,
        ),
        screen_field(
            "success", "MAGMA target genes retained",
            f'{summary["target_genes"]:,}',
            indent=10, label_width=label_width,
        ),
        screen_field(
            "warning" if summary["excluded_target_genes"] else "success",
            "MAGMA target genes excluded",
            f'{summary["excluded_target_genes"]:,}',
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
    compatibility = summary["gene_compatibility"]
    if "audit_table" in compatibility:
        lines.extend([
            "",
            screen_line("decision", "Gene-universe audit", indent=6),
            screen_field(
                "info", "Gene-level decisions", compatibility["audit_table"],
                indent=10, label_width=label_width,
            ),
            screen_field(
                "info", "Compatibility report", compatibility["report"],
                indent=10, label_width=label_width,
            ),
        ])
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


def run_kpops_direct(
    args: argparse.Namespace,
    ctx=None,
    *,
    pipeline_resources=None,
):
    progress = StageProgress(
        "K-POPS analysis progress",
        enabled=True,
    )
    total_progress_stages = len(_KPOPS_PROGRESS_STAGES)
    active_progress_step = 0
    try:
        configuration = _resolved_configuration(args, pipeline=ctx is not None)
        if pipeline_resources is None:
            configuration, resources = validate_kpops_configuration(
                args,
                configuration=configuration,
                progress=progress,
            )
        else:
            configuration, resources = _reuse_kpops_pipeline_resources(
                configuration,
                pipeline_resources,
                progress,
            )
        compatibility = resources["outcome"]["compatibility"]
    except BaseException as exc:
        active_progress_step = getattr(progress, "active_stage_number", 1)
        progress.fail_step(
            active_progress_step,
            total_progress_stages,
            _KPOPS_PROGRESS_STAGES[active_progress_step - 1],
        )
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
        progress.close()
        raise
    try:
        module = configuration.modules.kpops
        output = Path(configuration.run.output_directory).expanduser().resolve()
        dataset = validate_filename_component(
            configuration.run.dataset_id,
            "dataset_id",
            error_type=KPopsError,
        )
        output.mkdir(parents=True, exist_ok=True)
        prefix, final_paths = _output_paths(output, dataset, module)
        completion = configured_output_path(
            output,
            module.output_layout.completion_manifest,
            error_type=KPopsError,
            dataset_id=dataset,
        )
        log_path = configured_output_path(
            output,
            module.output_layout.service_log_file,
            error_type=KPopsError,
            dataset_id=dataset,
        )
        logger = PipelineLogger(
            dataset,
            "run",
            str(log_path.parent),
            level=configuration.logging.file_level,
            screen_level=configuration.logging.console_level,
            log_path=str(log_path),
        )
        active_progress_step = 5
    except BaseException:
        progress.close()
        raise
    try:
        resolved_path = configured_output_path(output, module.output_layout.resolved_config_file, error_type=KPopsError, dataset_id=dataset)
        write_resolved_configuration(configuration, resolved_path, modules="kpops")
        logger.record(
            "OBSERVED", "kpops_magma_input_validation",
            genes_out=resources["outcome"]["genes_out"],
            genes_raw=resources["outcome"]["genes_raw"],
            genes=resources["outcome"]["original_gene_count"], status="PASSED",
        )
        logger.record(
            "OBSERVED", "kpops_resource_validation",
            script=resources["script"], python=resources["python"],
            annotation=resources["annotation"],
            annotation_genes=resources["annotation_gene_count"],
            kernel=resources["kernel"], kernel_genes=resources["kernel_genes"],
            kernel_gene_count=resources["kernel_gene_count"], status="PASSED",
        )
        logger.record(
            "PARAM", "kpops_run",
            genome_build=module.genome_build.value,
            gene_universe_policy=module.gene_universe_policy,
            training_chromosomes=module.training_chromosomes,
            use_magma_covariates=module.use_magma_covariates,
        )
        logger.record(
            "OBSERVED", "kpops_gene_compatibility",
            **{
                key: value
                for key, value in compatibility.items()
                if not key.endswith("_gene_ids")
            },
        )
        completion_inputs = _completion_inputs(resources)
        completion_digest = _completion_configuration(configuration)
        if (
            configuration.run.resume
            and not configuration.run.overwrite
            and completion.is_file()
        ):
            decision = resolve_completion_resume(
                completion, dataset_id=dataset, module="kpops",
                genome_build=module.genome_build.value,
                configuration_sha256=completion_digest,
                inputs=completion_inputs, outputs=final_paths,
                resume_policy=configuration.run.resume_policy,
                error_type=KPopsError,
            )
            if decision.action == "resume":
                active_progress_step = 5
                progress.start_step(
                    active_progress_step,
                    total_progress_stages,
                    _KPOPS_PROGRESS_STAGES[active_progress_step - 1],
                )
                progress.complete_step(
                    active_progress_step,
                    total_progress_stages,
                    _KPOPS_PROGRESS_STAGES[active_progress_step - 1],
                    outcome_fields=[
                        (
                            "decision", "Model execution",
                            "reused checksum-validated outputs",
                        ),
                    ],
                )
                active_progress_step = 6
                progress.start_step(
                    active_progress_step,
                    total_progress_stages,
                    _KPOPS_PROGRESS_STAGES[active_progress_step - 1],
                )
                summary = _summarise_outputs(
                    final_paths["predictions"], module, resources,
                    paths=final_paths,
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
                progress.complete_step(
                    active_progress_step,
                    total_progress_stages,
                    _KPOPS_PROGRESS_STAGES[active_progress_step - 1],
                    outcome_fields=[
                        (
                            "count", "Genes with finite K-POPS scores",
                            summary["genes_scored"],
                        ),
                        (
                            "success", "Published result files",
                            len(final_paths),
                        ),
                    ],
                )
                active_progress_step = 0
                print(_render_summary(
                    summary, dataset, module, final_paths, log_path,
                    configuration.logging.terminal_label_width,
                ))
                return result
            apply_completion_restart(
                decision,
                output_root=output,
                manifest=completion,
                logger=logger,
                operation="kpops_resume",
                error_type=KPopsError,
            )
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
                effective_magma_prefix = module.magma_association_prefix
                if module.gene_universe_policy == "intersect":
                    effective_magma_prefix, compatibility_report = (
                        _prepare_intersected_magma(
                            module,
                            resources["outcome"],
                            staged_paths,
                            final_paths,
                            Path(directory)
                            / module.output_layout.compatible_magma_prefix,
                        )
                    )
                    logger.record(
                        "ACTION", "kpops_gene_universe_intersection",
                        **{
                            key: value
                            for key, value in compatibility_report.items()
                            if key not in {"chromosomes", "files"}
                        },
                    )
                    for chromosome, counts in compatibility_report[
                        "chromosomes"
                    ].items():
                        logger.record(
                            "OBSERVED", "kpops_gene_compatibility_chromosome",
                            chromosome=chromosome, **counts,
                        )
                    logger.record(
                        "OUTPUT", "kpops_gene_compatibility_files",
                        **compatibility_report["files"],
                    )
                    if compatibility_report["excluded_target_genes"]:
                        logger.warning(
                            "K-POPS gene-universe intersection retained %d of %d "
                            "MAGMA target genes (%.1f%%) and excluded %d; original "
                            "MAGMA files were not modified."
                            % (
                                compatibility_report["retained_target_genes"],
                                compatibility_report["original_target_genes"],
                                compatibility_report["retained_percent"],
                                compatibility_report["excluded_target_genes"],
                            )
                        )
                active_progress_step = 5
                progress.start_step(
                    active_progress_step,
                    total_progress_stages,
                    _KPOPS_PROGRESS_STAGES[active_progress_step - 1],
                )
                command = _arguments(
                    configuration,
                    staged_prefix,
                    resources["script"],
                    magma_prefix=effective_magma_prefix,
                    emit_model_progress=True,
                )
                command[0] = str(resources["python"])
                model_fit_count = _model_fit_count(module, resources)
                model_progress = MeasuredProgress(
                    "K-POPS model-fitting progress",
                    enabled=configuration.logging.show_progress,
                )
                model_progress.start(
                    "Fit K-POPS chromosome models",
                    total=model_fit_count,
                )
                with tempfile.NamedTemporaryFile(
                    prefix="kpops_native_",
                    suffix=".log",
                    dir=directory,
                    delete=False,
                ) as native_handle:
                    native_log = Path(native_handle.name)
                observe_model_progress = _model_progress_observer(
                    native_log,
                    model_progress,
                    model_fit_count,
                    logger,
                )
                try:
                    run_checked_command(
                        command,
                        "K-POPS",
                        logger=logger,
                        error_type=KPopsError,
                        stdout_path=native_log,
                        stderr_to_stdout=True,
                        timeout_seconds=configuration.execution.timeout_seconds,
                        expected_outputs=list(staged_paths.values()),
                        progress_callback=observe_model_progress,
                        progress_refresh_seconds=(
                            configuration.logging.progress_refresh_seconds
                        ),
                    )
                    observe_model_progress()
                except BaseException:
                    model_progress.fail(
                        title="Fit K-POPS chromosome models",
                    )
                    raise
                else:
                    model_progress.complete(
                        model_fit_count,
                        total=model_fit_count,
                        title="Fit K-POPS chromosome models",
                    )
                finally:
                    if native_log.is_file():
                        native_output = native_log.read_text(
                            encoding="utf-8",
                            errors="replace",
                        ).strip()
                        if native_output:
                            logger.info("K-POPS output: %s" % native_output)
                progress.complete_step(
                    active_progress_step,
                    total_progress_stages,
                    _KPOPS_PROGRESS_STAGES[active_progress_step - 1],
                    outcome_fields=[
                        (
                            "decision", "Training design",
                            " ".join(module.training_chromosomes),
                        ),
                        ("info", "Compute device", module.device),
                        (
                            "success", "Scored target genes submitted",
                            compatibility["retained_target_genes"],
                        ),
                    ],
                )
                active_progress_step = 6
                progress.start_step(
                    active_progress_step,
                    total_progress_stages,
                    _KPOPS_PROGRESS_STAGES[active_progress_step - 1],
                )
                summary = _summarise_outputs(
                    staged_paths["predictions"], module, resources,
                    paths=final_paths,
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
                "original_target_gene_count": summary["original_target_genes"],
                "excluded_target_gene_count": summary["excluded_target_genes"],
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
        progress.complete_step(
            active_progress_step,
            total_progress_stages,
            _KPOPS_PROGRESS_STAGES[active_progress_step - 1],
            outcome_fields=[
                (
                    "success", "Output validation",
                    "all required files are non-empty",
                ),
                (
                    "count", "Genes with finite K-POPS scores",
                    summary["genes_scored"],
                ),
                (
                    "warning" if summary["warnings"] else "success",
                    "Scientific warnings", len(summary["warnings"]),
                ),
                ("success", "Published result files", len(final_paths)),
            ],
        )
        active_progress_step = 0
        if ctx is not None:
            ctx["kpops"] = result
        print(_render_summary(
            summary, dataset, module, final_paths, log_path,
            configuration.logging.terminal_label_width,
        ))
        return result
    except BaseException as exc:
        if active_progress_step:
            progress.fail_step(
                active_progress_step,
                total_progress_stages,
                _KPOPS_PROGRESS_STAGES[active_progress_step - 1],
            )
        logger.error("K-POPS analysis failed: %s: %s" % (type(exc).__name__, exc))
        raise
    finally:
        progress.close()
        logger.close()


__all__ = ["preflight_kpops_pipeline", "run_kpops_direct", "validate_kpops_configuration"]
