"""Validated service boundary for the upstream CALDERA R implementation."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import tempfile

import numpy as np
import pandas as pd

from postgwas.config import (
    load_configuration,
    load_run_configuration_for_module,
    write_resolved_configuration,
)
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.config.models.common import GenomeBuild
from postgwas.core.completion import (
    apply_completion_restart,
    configuration_digest,
    resolve_completion_resume,
    write_completion_manifest,
)
from postgwas.core.io.tables import read_pandas_table, require_table_columns
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
from postgwas.core.r_runtime import resolve_r_runtime
from postgwas.core.reference_resources import contained_resource_path
from postgwas.core.required_arguments import (
    RequiredArgument,
    require_resolved_arguments,
)
from postgwas.core.ui import StageProgress
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.modules.caldera.errors import CalderaError


_CALDERA_PROGRESS_STAGES = (
    "Validate PoPS gene-prioritisation input",
    "Validate credible-set input",
    "Validate CALDERA reference and software resources",
    "Run CALDERA gene prioritisation",
    "Validate and publish CALDERA results",
)

# Required by the pinned upstream caldera() implementation, not optional model
# settings: https://github.com/kheilbron/caldera/blob/81a8a0308741ae986660f711bbf6a7abbd3bca19/z_caldera.R.
_CALDERA_REQUIRED_R_PACKAGES = ("data.table", "dplyr")


def _resolved_configuration(args: argparse.Namespace):
    module_overrides = explicit_overrides(args, {
        "caldera_genome_build": "genome_build",
        "caldera_repository": "repository_path",
        "caldera_adapter_script": "adapter_script_path",
        "pops_file": "pops_file",
        "credible_set_file": "credible_set_file",
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
        "caldera", getattr(args, "run_config", None),
        module_overrides=module_overrides, global_overrides=global_overrides,
    )


def _repository_file(repository: Path, relative: str, label: str) -> Path:
    candidate = contained_resource_path(
        repository, relative, label=label,
        root_label="the configured CALDERA repository", error_type=CalderaError,
    )
    return require_nonempty_file(candidate, label, error_type=CalderaError)


def _installed_repository(configuration, module) -> Path:
    python = Path(resolve_executable(
        configuration.resources.executables.python,
        "Python",
        error_type=CalderaError,
    )).resolve()
    environment = python.parent.parent
    repository = (environment / module.installed_repository_relative_path).resolve()
    if environment != repository and environment not in repository.parents:
        raise CalderaError(
            "Configured CALDERA installation path leaves the active Python environment"
        )
    return repository


def _validate_pops(path: Path, module) -> dict:
    schema = module.input_schema
    table = read_pandas_table(
        path, schema.table_delimiter, "PoPS predictions", engine="c", error_type=CalderaError,
    )
    required = [schema.pops_gene_id_column, schema.pops_score_column]
    require_table_columns(table, required, "PoPS predictions", error_type=CalderaError)
    identifiers = table[required[0]].astype(str)
    scores = pd.to_numeric(table[required[1]], errors="coerce")
    if identifiers.str.strip().eq("").any() or identifiers.duplicated().any():
        raise CalderaError("PoPS gene identifiers must be non-empty and unique")
    if scores.isna().any() or not np.isfinite(scores.to_numpy()).all():
        raise CalderaError("PoPS scores must all be finite")
    return {"genes": len(table)}


def _validate_credible_sets(path: Path, module) -> dict:
    table = read_pandas_table(
        path, None, "CALDERA credible sets", engine="python", error_type=CalderaError,
    )
    return _validate_credible_set_table(table, module)


def _validate_credible_set_table(table: pd.DataFrame, module) -> dict:
    schema = module.input_schema
    required = [
        schema.credible_set_locus_column,
        schema.credible_set_chromosome_column,
        schema.credible_set_position_column,
        schema.credible_set_probability_column,
    ]
    require_table_columns(table, required, "CALDERA credible sets", error_type=CalderaError)
    coding_columns = {
        schema.credible_set_coding_gene_column,
        schema.credible_set_rsid_column,
    }
    if module.genome_build == GenomeBuild.GRCH38 and not coding_columns.intersection(table.columns):
        raise CalderaError(
            "GRCh38 CALDERA credible sets require the configured c_gene or rsid "
            "column because upstream coding-variant coordinates are GRCh37"
        )
    if table[required].isna().any().any():
        raise CalderaError(
            "Credible-set locus, chromosome, position, and PIP values cannot be missing; "
            "PostGWAS will not allow CALDERA to discard them silently"
        )
    loci = table[required[0]].astype(str).str.strip()
    chromosomes = table[required[1]].astype(str).str.strip().str.replace(r"(?i)^chr", "", regex=True)
    positions = pd.to_numeric(table[required[2]], errors="coerce")
    probabilities = pd.to_numeric(table[required[3]], errors="coerce")
    if loci.eq("").any() or chromosomes.eq("").any():
        raise CalderaError("Credible-set locus and chromosome labels must not be empty")
    chromosome_numbers = pd.to_numeric(chromosomes, errors="coerce")
    if (
        chromosome_numbers.isna().any()
        or not np.isfinite(chromosome_numbers).all()
        or (chromosome_numbers < 1).any()
        or (chromosome_numbers % 1 != 0).any()
    ):
        raise CalderaError(
            "Credible-set chromosomes must be positive integer labels accepted by CALDERA"
        )
    if positions.isna().any() or not np.isfinite(positions).all() or (positions < 1).any() or (positions % 1 != 0).any():
        raise CalderaError("Credible-set positions must be positive integers")
    if probabilities.isna().any() or not np.isfinite(probabilities).all() or ((probabilities < 0) | (probabilities > 1)).any():
        raise CalderaError("Credible-set PIPs must be finite values in [0, 1]")
    normalized = table.copy()
    normalized[required[0]] = loci
    normalized[required[1]] = chromosomes
    normalized[required[2]] = positions.astype(np.int64)
    normalized[required[3]] = probabilities.astype(float)
    chromosome_counts = normalized.groupby(required[0], sort=False)[required[1]].nunique()
    if (chromosome_counts != 1).any():
        raise CalderaError("Each credible-set locus must contain exactly one chromosome")
    pip_mass = normalized.groupby(required[0], sort=False)[required[3]].sum()
    # Upstream CALDERA rejects a locus only when its total PIP is below the
    # boundary. At an exact boundary, its empty `which(total > boundary)` result
    # leaves all variants in the credible set, so exact coverage is valid.
    insufficient = pip_mass[pip_mass < module.minimum_credible_set_coverage]
    if not insufficient.empty:
        raise CalderaError(
            "%d credible sets are below the configured %.3g cumulative PIP "
            "required by upstream CALDERA"
            % (len(insufficient), module.minimum_credible_set_coverage)
        )
    retained = 0
    for _, group in normalized.groupby(required[0], sort=False):
        cumulative = group.sort_values(required[3], ascending=False, kind="stable")[required[3]].cumsum()
        above = np.flatnonzero(
            cumulative.to_numpy() > module.minimum_credible_set_coverage
        )
        retained += int(above[0]) + 1 if len(above) else len(group)
    return {
        "rows": len(normalized), "loci": normalized[required[0]].nunique(),
        "upstream_retained_rows": retained,
        "upstream_trimmed_rows": len(normalized) - retained,
    }


def validate_caldera_configuration(
    args: argparse.Namespace,
    *,
    pipeline: bool = False,
    configuration=None,
    pipeline_pops: str | Path | None = None,
    pipeline_credible_sets_directory: str | Path | None = None,
    validated_pipeline_resources=None,
    progress: StageProgress | None = None,
):
    configuration = configuration or _resolved_configuration(args)
    module = configuration.modules.caldera
    if not pipeline:
        require_resolved_arguments([
            RequiredArgument("--pops-file", "modules.caldera.pops_file", module.pops_file),
            RequiredArgument(
                "--credible-set-file",
                "modules.caldera.credible_set_file",
                module.credible_set_file,
            ),
        ])
    if pipeline and module.genome_build != configuration.modules.fine_mapping.genome_build:
        raise CalderaError(
            "CALDERA genome build %s does not match pipeline fine-mapping build %s"
            % (module.genome_build.value, configuration.modules.fine_mapping.genome_build.value)
        )
    if pipeline and module.genome_build == GenomeBuild.GRCH38:
        raise CalderaError(
            "CALDERA pipeline conversion cannot preserve c_gene or rsid from the "
            "current fine-mapping interchange; use direct mode with an annotated "
            "GRCh38 credible-set file"
        )
    files = {}
    pops_value = pipeline_pops if pipeline else module.pops_file
    credible_value = None if pipeline else module.credible_set_file
    validate_runtime_inputs = not pipeline or (
        pipeline_pops is not None and pipeline_credible_sets_directory is not None
    )
    if validate_runtime_inputs:
        if progress is not None:
            progress.active_stage_number = 1
            progress.start_step(1, len(_CALDERA_PROGRESS_STAGES), _CALDERA_PROGRESS_STAGES[0])
        pops = require_nonempty_file(
            pops_value, "PoPS predictions", error_type=CalderaError,
        )
        pops_metrics = _validate_pops(pops, module)
        files.update({"pops": pops, "pops_metrics": pops_metrics})
        if progress is not None:
            progress.complete_step(
                1, len(_CALDERA_PROGRESS_STAGES), _CALDERA_PROGRESS_STAGES[0],
                outcome_fields=[
                    ("info", "PoPS predictions", pops),
                    ("count", "Unique genes", pops_metrics["genes"]),
                    ("success", "Gene identifiers and scores", "valid"),
                ],
            )

        if progress is not None:
            progress.active_stage_number = 2
            progress.start_step(2, len(_CALDERA_PROGRESS_STAGES), _CALDERA_PROGRESS_STAGES[1])
        if pipeline:
            credible_table, credible_metrics = _collect_finemap_credible_sets(
                pipeline_credible_sets_directory, module,
            )
            files["pipeline_credible_set_table"] = credible_table
        else:
            credible_sets = require_nonempty_file(
                credible_value, "CALDERA credible sets", error_type=CalderaError,
            )
            credible_metrics = _validate_credible_sets(credible_sets, module)
            files["credible_sets"] = credible_sets
        files["credible_set_metrics"] = credible_metrics
        if progress is not None:
            fields = [
                ("count", "Credible-set loci", credible_metrics["loci"]),
                ("count", "Credible-set variants", credible_metrics["rows"]),
                ("count", "Variants retained upstream", credible_metrics["upstream_retained_rows"]),
                ("warning" if credible_metrics["upstream_trimmed_rows"] else "success", "Variants trimmed upstream", credible_metrics["upstream_trimmed_rows"]),
            ]
            if pipeline:
                fields.insert(0, ("info", "Fine-mapping index", credible_metrics["source_index"]))
                fields.insert(3, ("count", "Credible-set files", credible_metrics["credible_set_files"]))
            else:
                fields.insert(0, ("info", "Credible-set file", credible_sets))
            progress.complete_step(
                2, len(_CALDERA_PROGRESS_STAGES), _CALDERA_PROGRESS_STAGES[1],
                outcome_fields=fields,
            )

    preflight_configuration = None
    validated_external = None
    if validated_pipeline_resources is not None:
        if not (
            isinstance(validated_pipeline_resources, tuple)
            and len(validated_pipeline_resources) == 2
            and isinstance(validated_pipeline_resources[1], dict)
        ):
            raise CalderaError(
                "CALDERA pipeline execution requires completed resource-preflight evidence"
            )
        if progress is not None:
            progress.active_stage_number = 3
            progress.start_step(
                3,
                len(_CALDERA_PROGRESS_STAGES),
                _CALDERA_PROGRESS_STAGES[2],
            )
        preflight_configuration, validated_external = validated_pipeline_resources
        excluded = {"pops_file", "credible_set_file"}
        if module.model_dump(mode="json", exclude=excluded) != (
            preflight_configuration.modules.caldera.model_dump(
                mode="json", exclude=excluded,
            )
        ):
            raise CalderaError(
                "CALDERA configuration changed after pipeline preflight; restart so "
                "resources are validated against the executed configuration"
            )
        if (
            configuration.resources.executables.rscript
            != preflight_configuration.resources.executables.rscript
        ):
            raise CalderaError(
                "CALDERA Rscript configuration changed after pipeline preflight; "
                "restart so the executable is validated before analysis"
            )
        identities = validated_external.get("preflight_file_identities")
        if not identities:
            raise CalderaError("CALDERA resource-preflight evidence is incomplete")
        require_unchanged_preflight_files(
            identities,
            error_type=CalderaError,
            label="CALDERA resource",
        )

    if progress is not None and validated_external is None:
        progress.active_stage_number = 3
        progress.start_step(3, len(_CALDERA_PROGRESS_STAGES), _CALDERA_PROGRESS_STAGES[2])
    if validated_external is None:
        repository = (
            Path(module.repository_path).expanduser().resolve()
            if module.repository_path is not None
            else _installed_repository(configuration, module)
        )
        if not repository.is_dir():
            raise CalderaError("CALDERA repository does not exist: %s" % repository)
        assembly = module.assembly_by_genome_build[module.genome_build]
        files.update({
            "repository": repository,
            "upstream_script": _repository_file(repository, module.upstream_script_relative_path, "CALDERA upstream script"),
            "coding_variants": _repository_file(repository, module.coding_variants_relative_path, "CALDERA coding variants"),
            "gene_locations": _repository_file(
                repository,
                module.gene_locations_relative_pattern.format(assembly=assembly),
                "CALDERA gene locations",
            ),
            "model": _repository_file(repository, module.model_relative_path, "CALDERA trained model"),
            "adapter": require_nonempty_file(
                module.adapter_script_path
                if module.adapter_script_path is not None
                else Path(__file__).with_name("run_caldera.R"),
                "PostGWAS CALDERA adapter",
                error_type=CalderaError,
            ),
        })
        files["rscript"] = Path(resolve_executable(
            configuration.resources.executables.rscript,
            "Rscript",
            error_type=CalderaError,
        ))
        try:
            files["r_runtime"] = resolve_r_runtime(
                files["rscript"], configuration.execution.timeout_seconds,
                _CALDERA_REQUIRED_R_PACKAGES, label="CALDERA",
            )
        except (RuntimeError, OSError) as exc:
            raise CalderaError(str(exc)) from exc
        files["preflight_file_identities"] = capture_preflight_file_identities(
            (
                files["upstream_script"],
                files["coding_variants"],
                files["gene_locations"],
                files["model"],
                files["adapter"],
                files["rscript"],
            ),
            error_type=CalderaError,
            label="CALDERA resource",
        )
    else:
        for key in (
            "repository",
            "upstream_script",
            "coding_variants",
            "gene_locations",
            "model",
            "adapter",
            "rscript",
            "r_runtime",
            "preflight_file_identities",
        ):
            if key not in validated_external:
                raise CalderaError(
                    "CALDERA resource-preflight evidence is missing %s" % key
                )
            files[key] = validated_external[key]
        repository = files["repository"]
    if progress is not None:
        progress.complete_step(
            3, len(_CALDERA_PROGRESS_STAGES), _CALDERA_PROGRESS_STAGES[2],
            outcome_fields=[
                ("info", "CALDERA repository", repository),
                ("info", "Upstream analysis script", files["upstream_script"]),
                ("info", "Coding-variant resource", files["coding_variants"]),
                ("info", "Gene-location resource", files["gene_locations"]),
                ("info", "Trained model", files["model"]),
                ("info", "PostGWAS adapter", files["adapter"]),
                ("info", "Rscript executable", files["rscript"]),
                ("info", "R runtime", files["r_runtime"]["version"]),
                ("genetic", "Declared genome build", module.genome_build.value),
                (
                    "decision" if validated_external is not None else "success",
                    "All reference files",
                    "validation reused; file identities unchanged"
                    if validated_external is not None
                    else "present and non-empty",
                ),
            ],
        )
    return configuration, files


def preflight_caldera_pipeline(
    args: argparse.Namespace,
    *,
    preflight_evidence=None,
) -> PipelinePreflightEvidence:
    require_pipeline_input_vcf(preflight_evidence)
    resources = validate_caldera_configuration(args, pipeline=True)
    return pipeline_preflight_evidence(
        "caldera",
        preflight_evidence,
        resources=resources,
        deferred_checks=(
            "Validate pipeline-generated PoPS and fine-mapping inputs.",
        ),
    )


def _collect_finemap_credible_sets(
    directory: str | Path,
    module,
) -> tuple[pd.DataFrame, dict]:
    """Read and validate the pipeline interchange without creating outputs."""
    root = Path(directory).expanduser().resolve()
    if not root.is_dir():
        raise CalderaError("Fine-mapping credible-set directory does not exist: %s" % root)
    mapping = module.pipeline_input
    index_path, index, source_paths = _finemap_source_paths(root, module)
    pattern = re.compile(mapping.variant_identifier_pattern)
    records = []
    seen_loci = set()
    for (row_number, row), source in zip(index.iterrows(), source_paths):
        credible = read_pandas_table(
            source, None, "fine-mapping credible-set file", engine="python", error_type=CalderaError,
        )
        require_table_columns(
            credible, [mapping.credible_variant_column, mapping.credible_probability_column],
            "fine-mapping credible-set file", error_type=CalderaError,
        )
        base_locus = str(row[mapping.index_locus_column]).strip()
        locus = base_locus + mapping.locus_separator + source.name
        if not base_locus or locus in seen_loci:
            raise CalderaError("Fine-mapping credible-set index contains an empty or duplicate locus/file pair")
        seen_loci.add(locus)
        for _, variant in credible.iterrows():
            match = pattern.fullmatch(str(variant[mapping.credible_variant_column]).strip())
            if match is None:
                raise CalderaError("Cannot parse fine-mapping variant identifier: %s" % variant[mapping.credible_variant_column])
            records.append({
                module.input_schema.credible_set_locus_column: locus,
                module.input_schema.credible_set_chromosome_column: match.group("chromosome"),
                module.input_schema.credible_set_position_column: int(match.group("position")),
                module.input_schema.credible_set_probability_column: variant[mapping.credible_probability_column],
            })
    if not records:
        raise CalderaError("Fine-mapping credible-set index contains no variants")
    table = pd.DataFrame(records)
    metrics = _validate_credible_set_table(table, module)
    metrics.update({
        "source_index": str(index_path), "credible_set_files": len(index),
        "source_files": [str(path) for path in source_paths],
    })
    return table, metrics


def convert_finemap_credible_sets(directory: str | Path, destination: Path, module) -> dict:
    """Convert the validated fine-mapping interchange into CALDERA's CS schema."""
    table, metrics = _collect_finemap_credible_sets(directory, module)
    destination.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(destination, sep=module.input_schema.table_delimiter, index=False)
    return metrics


def _finemap_source_paths(root: Path, module) -> tuple[Path, pd.DataFrame, list[Path]]:
    mapping = module.pipeline_input
    index_path = next((root / name for name in mapping.index_file_names if (root / name).is_file()), None)
    if index_path is None:
        raise CalderaError(
            "Fine-mapping output lacks a configured credible-set index: %s"
            % ", ".join(mapping.index_file_names)
        )
    index = read_pandas_table(
        index_path, None, "fine-mapping credible-set index", engine="python", error_type=CalderaError,
    )
    require_table_columns(
        index, [mapping.index_filename_column, mapping.index_locus_column],
        "fine-mapping credible-set index", error_type=CalderaError,
    )
    source_paths = []
    for row_number, row in index.iterrows():
        filename = str(row[mapping.index_filename_column]).strip()
        source = Path(filename).expanduser()
        if not source.is_absolute():
            source = root / source
        source = source.resolve()
        if root != source and root not in source.parents:
            raise CalderaError("Credible-set index row %d leaves its output directory" % (row_number + 1))
        source = require_nonempty_file(source, "fine-mapping credible-set file", error_type=CalderaError)
        source_paths.append(source)
    return index_path.resolve(), index, source_paths


def _completion_configuration(configuration) -> str:
    return configuration_digest({
        "module": configuration.modules.caldera.model_dump(mode="json"),
        "rscript": configuration.resources.executables.rscript,
    })


def _base_completion_inputs(resources: dict, pops: Path) -> dict[str, Path]:
    return {
        "pops_predictions": pops,
        "caldera_upstream_script": resources["upstream_script"],
        "caldera_coding_variants": resources["coding_variants"],
        "caldera_gene_locations": resources["gene_locations"],
        "caldera_model": resources["model"],
        "postgwas_adapter": resources["adapter"],
    }


def _validate_results(path: Path, module) -> dict:
    table = read_pandas_table(
        path, module.input_schema.table_delimiter, "CALDERA results",
        engine="c", error_type=CalderaError,
    )
    schema = module.result_schema
    require_table_columns(table, list(schema.model_dump().values()), "CALDERA results", error_type=CalderaError)
    probability = pd.to_numeric(table[schema.normalized_probability_column], errors="coerce")
    raw = pd.to_numeric(table[schema.raw_probability_column], errors="coerce")
    if probability.isna().any() or raw.isna().any() or not np.isfinite(probability).all() or not np.isfinite(raw).all():
        raise CalderaError("CALDERA probabilities must be finite")
    if ((probability < 0) | (probability > 1) | (raw < 0) | (raw > 1)).any():
        raise CalderaError("CALDERA probabilities must be in [0, 1]")
    loci = table[schema.locus_column].astype(str).str.strip()
    locus_positions = table[schema.locus_position_column].astype(str).str.strip()
    gene_names = table[schema.gene_name_column].astype(str).str.strip()
    genes = table[schema.gene_id_column].astype(str).str.strip()
    locus_gene_counts = pd.to_numeric(table[schema.locus_gene_count_column], errors="coerce")
    distance = pd.to_numeric(table[schema.distance_column], errors="coerce")
    pops = pd.to_numeric(table[schema.pops_score_column], errors="coerce")
    coding = pd.to_numeric(table[schema.coding_probability_column], errors="coerce")
    if (
        loci.eq("").any() or locus_positions.eq("").any()
        or gene_names.eq("").any() or genes.eq("").any()
        or table.assign(_locus=loci, _gene=genes)
        .duplicated(["_locus", "_gene"]).any()
    ):
        raise CalderaError("CALDERA locus/gene identifiers must be non-empty and unique per locus")
    numeric = np.column_stack([locus_gene_counts, distance, pops, coding])
    if not np.isfinite(numeric).all():
        raise CalderaError("CALDERA gene counts, distances, PoPS scores, and coding PIPs must be finite")
    if (
        (locus_gene_counts < 1).any() or (locus_gene_counts % 1 != 0).any()
        or (distance < 0).any() or (coding < 0).any() or (coding > 1).any()
    ):
        raise CalderaError("CALDERA gene counts, distances, or coding PIPs are invalid")
    imputed = table[schema.imputed_pops_column].astype(str).str.upper()
    if not imputed.isin({"TRUE", "FALSE"}).all():
        raise CalderaError("CALDERA impute_pops values must be boolean")
    observed_counts = table.assign(_locus=loci).groupby("_locus", sort=False).size()
    declared_counts = table.assign(_locus=loci, _count=locus_gene_counts).groupby("_locus", sort=False)["_count"].nunique()
    first_counts = table.assign(_locus=loci, _count=locus_gene_counts).groupby("_locus", sort=False)["_count"].first()
    if (declared_counts != 1).any() or not np.array_equal(observed_counts.to_numpy(), first_counts.to_numpy()):
        raise CalderaError("CALDERA n_genes must equal the published row count in every locus")
    sums = table.assign(_locus=loci, _probability=probability).groupby("_locus", sort=False)["_probability"].sum()
    if not np.allclose(sums.to_numpy(), 1.0, rtol=1e-8, atol=1e-10):
        raise CalderaError("CALDERA normalized probabilities must sum to one per locus")
    return {"genes": len(table), "loci": loci.nunique()}


def _render_caldera_summary(
    dataset: str,
    result_metrics: dict,
    credible_metrics: dict,
    results_path: Path,
    log_path: Path,
    label_width: int,
) -> str:
    lines = [
        "",
        screen_line("analysis", "CALDERA gene-prioritisation summary", indent=2),
        screen_field("info", "Dataset", dataset, indent=6, label_width=label_width),
        "",
        screen_line("genetic", "Scientific findings", indent=6),
        screen_field("count", "Credible-set loci analysed", credible_metrics["loci"], indent=10, label_width=label_width),
        screen_field("count", "Credible-set variants supplied", credible_metrics["rows"], indent=10, label_width=label_width),
        screen_field("count", "Variants retained by CALDERA", credible_metrics["upstream_retained_rows"], indent=10, label_width=label_width),
        screen_field("count", "Prioritised locus-gene rows", result_metrics["genes"], indent=10, label_width=label_width),
        screen_field("count", "Loci reported", result_metrics["loci"], indent=10, label_width=label_width),
        "",
        screen_field("success", "Complete CALDERA results", results_path, indent=6, label_width=label_width),
        screen_field("info", "Full log", log_path, indent=6, label_width=label_width),
        "",
    ]
    return "\n".join(lines)


def run_caldera_direct(
    args: argparse.Namespace,
    ctx=None,
    *,
    pipeline_resources=None,
):
    pipeline_directory = getattr(args, "finemap_credible_sets_directory", None)
    progress = StageProgress("CALDERA analysis progress", enabled=True)
    active_progress_step = 0
    try:
        configuration = _resolved_configuration(args)
        configuration, resources = validate_caldera_configuration(
            args,
            pipeline=pipeline_directory is not None,
            configuration=configuration,
            pipeline_pops=getattr(args, "pops_file", None),
            pipeline_credible_sets_directory=pipeline_directory,
            validated_pipeline_resources=pipeline_resources,
            progress=progress,
        )
    except BaseException as exc:
        active_progress_step = getattr(progress, "active_stage_number", 1)
        progress.fail_step(
            active_progress_step,
            len(_CALDERA_PROGRESS_STAGES),
            _CALDERA_PROGRESS_STAGES[active_progress_step - 1],
        )
        fallback = load_configuration()
        output = Path(
            getattr(args, "output_directory", None) or fallback.run.output_directory
        ).expanduser().resolve()
        raw_dataset = getattr(args, "dataset_id", None) or fallback.run.dataset_id
        try:
            dataset = validate_filename_component(raw_dataset, "dataset_id", error_type=CalderaError)
        except CalderaError:
            dataset = fallback.run.dataset_id
        log_path = configured_output_path(
            output, fallback.modules.caldera.output_layout.service_log_file,
            error_type=CalderaError, dataset_id=dataset,
        )
        write_log_record(
            log_path, "ERROR",
            "CALDERA configuration or input validation failed: %s: %s"
            % (type(exc).__name__, exc),
            sample_id=dataset, file_level=fallback.logging.file_level,
            screen_level=fallback.logging.console_level,
        )
        progress.close()
        raise
    module = configuration.modules.caldera
    output = Path(configuration.run.output_directory).expanduser().resolve()
    dataset = validate_filename_component(configuration.run.dataset_id, "dataset_id", error_type=CalderaError)
    output.mkdir(parents=True, exist_ok=True)
    results_path = configured_output_path(output, module.output_layout.results_file, error_type=CalderaError, dataset_id=dataset)
    converted_path = configured_output_path(output, module.output_layout.pipeline_credible_sets_file, error_type=CalderaError, dataset_id=dataset)
    completion = configured_output_path(output, module.output_layout.completion_manifest, error_type=CalderaError, dataset_id=dataset)
    log_path = configured_output_path(output, module.output_layout.service_log_file, error_type=CalderaError, dataset_id=dataset)
    logger = PipelineLogger(
        dataset, "run", str(log_path.parent), level=configuration.logging.file_level,
        screen_level=configuration.logging.console_level, log_path=str(log_path),
    )
    active_progress_step = 4
    try:
        logger.record(
            "OBSERVED", "caldera_pops_input_validation",
            path=resources["pops"], genes=resources["pops_metrics"]["genes"],
            status="PASSED",
        )
        logger.record(
            "OBSERVED", "caldera_credible_set_input_validation",
            rows=resources["credible_set_metrics"]["rows"],
            loci=resources["credible_set_metrics"]["loci"],
            retained=resources["credible_set_metrics"]["upstream_retained_rows"],
            trimmed=resources["credible_set_metrics"]["upstream_trimmed_rows"],
            status="PASSED",
        )
        logger.record(
            "OBSERVED", "caldera_resource_validation",
            repository=resources["repository"],
            upstream_script=resources["upstream_script"],
            coding_variants=resources["coding_variants"],
            gene_locations=resources["gene_locations"],
            model=resources["model"], adapter=resources["adapter"],
            rscript=resources["rscript"], status="PASSED",
        )
        logger.record(
            "OBSERVED", "caldera_r_runtime",
            version=resources["r_runtime"]["version"],
            packages=resources["r_runtime"]["package_versions"],
            library_paths=resources["r_runtime"]["library_paths"],
        )
        completion_digest = _completion_configuration(configuration)
        if pipeline_directory is not None:
            pops = resources["pops"]
            pops_metrics = resources["pops_metrics"]
            credible_metrics = resources["credible_set_metrics"]
            source_index = Path(credible_metrics["source_index"])
            source_files = [Path(path) for path in credible_metrics["source_files"]]
            completion_inputs = _base_completion_inputs(resources, pops)
            completion_inputs["fine_mapping_index"] = source_index
            completion_inputs.update({
                "fine_mapping_credible_set_%d" % number: path
                for number, path in enumerate(source_files, 1)
            })
        else:
            pops = resources["pops"]
            pops_metrics = resources["pops_metrics"]
            completion_inputs = _base_completion_inputs(resources, pops)
            completion_inputs["credible_sets"] = resources["credible_sets"]
        resolved_path = configured_output_path(output, module.output_layout.resolved_config_file, error_type=CalderaError, dataset_id=dataset)
        write_resolved_configuration(configuration, resolved_path, modules="caldera")
        expected = {"results": results_path}
        if pipeline_directory is not None:
            expected["pipeline_credible_sets"] = converted_path
        if (
            configuration.run.resume
            and not configuration.run.overwrite
            and completion.is_file()
        ):
            decision = resolve_completion_resume(
                completion, dataset_id=dataset, module="caldera",
                genome_build=module.genome_build.value,
                configuration_sha256=completion_digest,
                inputs=completion_inputs, outputs=expected,
                resume_policy=configuration.run.resume_policy,
                error_type=CalderaError,
            )
            if decision.action == "resume":
                active_progress_step = 4
                progress.start_step(
                    active_progress_step,
                    len(_CALDERA_PROGRESS_STAGES),
                    _CALDERA_PROGRESS_STAGES[active_progress_step - 1],
                )
                progress.complete_step(
                    active_progress_step,
                    len(_CALDERA_PROGRESS_STAGES),
                    _CALDERA_PROGRESS_STAGES[active_progress_step - 1],
                    outcome_fields=[("decision", "Model execution", "reused checksum-validated outputs")],
                )
                active_progress_step = 5
                progress.start_step(
                    active_progress_step,
                    len(_CALDERA_PROGRESS_STAGES),
                    _CALDERA_PROGRESS_STAGES[active_progress_step - 1],
                )
                if pipeline_directory is not None:
                    _validate_credible_sets(converted_path, module)
                metrics = _validate_results(results_path, module)
                result = {"status": "success", "caldera_file": str(results_path), "completion_manifest": str(completion), "published_files": [str(path) for path in expected.values()]}
                if ctx is not None:
                    ctx["caldera"] = result
                logger.record("SKIP", "caldera_run", reason="validated_complete_outputs")
                progress.complete_step(
                    active_progress_step,
                    len(_CALDERA_PROGRESS_STAGES),
                    _CALDERA_PROGRESS_STAGES[active_progress_step - 1],
                    outcome_fields=[
                        ("count", "Prioritised gene rows", metrics["genes"]),
                        ("count", "Loci reported", metrics["loci"]),
                        ("success", "Published result files", len(result["published_files"])),
                    ],
                )
                active_progress_step = 0
                print(_render_caldera_summary(
                    dataset,
                    metrics,
                    credible_metrics,
                    results_path,
                    log_path,
                    configuration.logging.terminal_label_width,
                ))
                return result
            apply_completion_restart(
                decision,
                output_root=output,
                manifest=completion,
                logger=logger,
                operation="caldera_resume",
                error_type=CalderaError,
            )
        existing = [path for path in expected.values() if path.exists()]
        if existing and not configuration.run.overwrite:
            raise CalderaError("Existing CALDERA outputs require --overwrite: %s" % ", ".join(map(str, existing)))
        staging_root = configured_output_path(output, module.output_layout.staging_directory, error_type=CalderaError, dataset_id=dataset)
        staging_root.mkdir(parents=True, exist_ok=True)
        try:
            with tempfile.TemporaryDirectory(prefix="run_", dir=staging_root) as directory:
                work = Path(directory)
                staged_results = work / results_path.name
                if pipeline_directory is not None:
                    staged_credible_sets = work / converted_path.name
                    resources["pipeline_credible_set_table"].to_csv(
                        staged_credible_sets,
                        sep=module.input_schema.table_delimiter,
                        index=False,
                    )
                else:
                    staged_credible_sets = resources["credible_sets"]
                    credible_metrics = resources["credible_set_metrics"]
                active_progress_step = 4
                progress.start_step(
                    active_progress_step,
                    len(_CALDERA_PROGRESS_STAGES),
                    _CALDERA_PROGRESS_STAGES[active_progress_step - 1],
                )
                rscript = resources["rscript"]
                assembly = module.assembly_by_genome_build[module.genome_build]
                run_checked_command(
                    [
                        rscript, "--vanilla", resources["adapter"], resources["repository"],
                        resources["upstream_script"], pops, staged_credible_sets,
                        str(assembly), staged_results, module.input_schema.table_delimiter,
                    ],
                    "CALDERA", logger=logger, error_type=CalderaError,
                    timeout_seconds=configuration.execution.timeout_seconds,
                    env=resources["r_runtime"]["environment"],
                    expected_outputs=[staged_results],
                )
                progress.complete_step(
                    active_progress_step,
                    len(_CALDERA_PROGRESS_STAGES),
                    _CALDERA_PROGRESS_STAGES[active_progress_step - 1],
                    outcome_fields=[
                        ("genetic", "Genome build", module.genome_build.value),
                        ("count", "PoPS genes submitted", pops_metrics["genes"]),
                        ("count", "Credible-set loci submitted", credible_metrics["loci"]),
                    ],
                )
                active_progress_step = 5
                progress.start_step(
                    active_progress_step,
                    len(_CALDERA_PROGRESS_STAGES),
                    _CALDERA_PROGRESS_STAGES[active_progress_step - 1],
                )
                result_metrics = _validate_results(staged_results, module)
                if configuration.run.overwrite:
                    completion.unlink(missing_ok=True)
                    results_path.unlink(missing_ok=True)
                    if pipeline_directory is not None:
                        converted_path.unlink(missing_ok=True)
                results_path.parent.mkdir(parents=True, exist_ok=True)
                staged_results.replace(results_path)
                published = [results_path]
                if pipeline_directory is not None:
                    converted_path.parent.mkdir(parents=True, exist_ok=True)
                    staged_credible_sets.replace(converted_path)
                    published.append(converted_path)
        finally:
            remove_empty_directories(staging_root, staging_root.parent)
        result = {"status": "success", "caldera_file": str(results_path), "published_files": [str(path) for path in published]}
        write_completion_manifest(
            completion, dataset_id=dataset, module="caldera",
            genome_build=module.genome_build.value,
            configuration_sha256=completion_digest,
            inputs=completion_inputs, outputs=expected,
            metrics={
                "pops": pops_metrics, "credible_sets": credible_metrics,
                "results": result_metrics,
            },
            error_type=CalderaError,
        )
        result["completion_manifest"] = str(completion)
        logger.record(
            "STATUS", "caldera_run", status="COMPLETED",
            genes=result_metrics["genes"], loci=result_metrics["loci"],
            credible_set_rows=credible_metrics["rows"],
            upstream_trimmed_rows=credible_metrics["upstream_trimmed_rows"],
        )
        progress.complete_step(
            active_progress_step,
            len(_CALDERA_PROGRESS_STAGES),
            _CALDERA_PROGRESS_STAGES[active_progress_step - 1],
            outcome_fields=[
                ("success", "Output validation", "probabilities and locus counts valid"),
                ("count", "Prioritised gene rows", result_metrics["genes"]),
                ("count", "Loci reported", result_metrics["loci"]),
                ("success", "Published result files", len(published)),
            ],
        )
        active_progress_step = 0
        if ctx is not None:
            ctx["caldera"] = result
        print(_render_caldera_summary(
            dataset,
            result_metrics,
            credible_metrics,
            results_path,
            log_path,
            configuration.logging.terminal_label_width,
        ))
        return result
    except BaseException as exc:
        if active_progress_step:
            progress.fail_step(
                active_progress_step,
                len(_CALDERA_PROGRESS_STAGES),
                _CALDERA_PROGRESS_STAGES[active_progress_step - 1],
            )
        logger.error("CALDERA analysis failed: %s: %s" % (type(exc).__name__, exc))
        raise
    finally:
        progress.close()
        logger.close()


__all__ = [
    "convert_finemap_credible_sets", "preflight_caldera_pipeline",
    "run_caldera_direct", "validate_caldera_configuration",
]
