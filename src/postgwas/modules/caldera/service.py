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
from postgwas.modules.caldera.errors import CalderaError


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


def _read_table(path: Path, label: str, delimiter=None) -> pd.DataFrame:
    try:
        table = pd.read_csv(
            path,
            sep=delimiter,
            engine="python" if delimiter is None else "c",
        )
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise CalderaError("Cannot read %s %s: %s" % (label, path, exc)) from exc
    if table.empty:
        raise CalderaError("%s contains no data rows: %s" % (label, path))
    return table


def _require_columns(table: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise CalderaError("%s is missing required columns: %s" % (label, ", ".join(missing)))


def _repository_file(repository: Path, relative: str, label: str) -> Path:
    candidate = (repository / relative).resolve()
    if repository != candidate and repository not in candidate.parents:
        raise CalderaError("%s leaves the configured CALDERA repository" % label)
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
    table = _read_table(path, "PoPS predictions", schema.table_delimiter)
    required = [schema.pops_gene_id_column, schema.pops_score_column]
    _require_columns(table, required, "PoPS predictions")
    identifiers = table[required[0]].astype(str)
    scores = pd.to_numeric(table[required[1]], errors="coerce")
    if identifiers.str.strip().eq("").any() or identifiers.duplicated().any():
        raise CalderaError("PoPS gene identifiers must be non-empty and unique")
    if scores.isna().any() or not np.isfinite(scores.to_numpy()).all():
        raise CalderaError("PoPS scores must all be finite")
    return {"genes": len(table)}


def _validate_credible_sets(path: Path, module) -> dict:
    schema = module.input_schema
    table = _read_table(path, "CALDERA credible sets")
    required = [
        schema.credible_set_locus_column,
        schema.credible_set_chromosome_column,
        schema.credible_set_position_column,
        schema.credible_set_probability_column,
    ]
    _require_columns(table, required, "CALDERA credible sets")
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


def validate_caldera_configuration(args: argparse.Namespace, *, pipeline: bool = False):
    configuration = _resolved_configuration(args)
    module = configuration.modules.caldera
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
    repository = (
        Path(module.repository_path).expanduser().resolve()
        if module.repository_path is not None
        else _installed_repository(configuration, module)
    )
    if not repository.is_dir():
        raise CalderaError("CALDERA repository does not exist: %s" % repository)
    assembly = module.assembly_by_genome_build[module.genome_build]
    files = {
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
    }
    if not pipeline:
        pops = require_nonempty_file(module.pops_file, "PoPS predictions", error_type=CalderaError)
        credible_sets = require_nonempty_file(module.credible_set_file, "CALDERA credible sets", error_type=CalderaError)
        files.update({"pops": pops, "credible_sets": credible_sets})
        files["pops_metrics"] = _validate_pops(pops, module)
        files["credible_set_metrics"] = _validate_credible_sets(credible_sets, module)
    return configuration, files


def preflight_caldera_pipeline(args: argparse.Namespace) -> None:
    validate_caldera_configuration(args, pipeline=True)


def convert_finemap_credible_sets(directory: str | Path, destination: Path, module) -> dict:
    """Convert the validated fine-mapping interchange into CALDERA's CS schema."""
    root = Path(directory).expanduser().resolve()
    if not root.is_dir():
        raise CalderaError("Fine-mapping credible-set directory does not exist: %s" % root)
    mapping = module.pipeline_input
    index_path, index, source_paths = _finemap_source_paths(root, module)
    pattern = re.compile(mapping.variant_identifier_pattern)
    records = []
    seen_loci = set()
    for (row_number, row), source in zip(index.iterrows(), source_paths):
        credible = _read_table(source, "fine-mapping credible-set file")
        _require_columns(credible, [mapping.credible_variant_column, mapping.credible_probability_column], "fine-mapping credible-set file")
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
    destination.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(destination, sep=module.input_schema.table_delimiter, index=False)
    metrics = _validate_credible_sets(destination, module)
    metrics.update({
        "source_index": str(index_path), "credible_set_files": len(index),
        "source_files": [str(path) for path in source_paths],
    })
    return metrics


def _finemap_source_paths(root: Path, module) -> tuple[Path, pd.DataFrame, list[Path]]:
    mapping = module.pipeline_input
    index_path = next((root / name for name in mapping.index_file_names if (root / name).is_file()), None)
    if index_path is None:
        raise CalderaError(
            "Fine-mapping output lacks a configured credible-set index: %s"
            % ", ".join(mapping.index_file_names)
        )
    index = _read_table(index_path, "fine-mapping credible-set index")
    _require_columns(index, [mapping.index_filename_column, mapping.index_locus_column], "fine-mapping credible-set index")
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
    table = _read_table(path, "CALDERA results", module.input_schema.table_delimiter)
    schema = module.result_schema
    _require_columns(table, list(schema.model_dump().values()), "CALDERA results")
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


def run_caldera_direct(args: argparse.Namespace, ctx=None):
    pipeline_directory = getattr(args, "finemap_credible_sets_directory", None)
    try:
        configuration, resources = validate_caldera_configuration(
            args, pipeline=pipeline_directory is not None,
        )
    except BaseException as exc:
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
    try:
        completion_digest = _completion_configuration(configuration)
        if pipeline_directory is not None:
            pops = require_nonempty_file(getattr(args, "pops_file", None), "PoPS predictions", error_type=CalderaError)
            pops_metrics = _validate_pops(pops, module)
            source_index, _, source_files = _finemap_source_paths(
                Path(pipeline_directory).expanduser().resolve(), module,
            )
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
        if configuration.run.resume and not configuration.run.overwrite and completion.is_file() and all(path.is_file() for path in expected.values()):
            validate_completion_manifest(
                completion, dataset_id=dataset, module="caldera",
                genome_build=module.genome_build.value,
                configuration_sha256=completion_digest,
                inputs=completion_inputs, outputs=expected,
                error_type=CalderaError,
            )
            if pipeline_directory is not None:
                _validate_credible_sets(converted_path, module)
            metrics = _validate_results(results_path, module)
            result = {"status": "success", "caldera_file": str(results_path), "completion_manifest": str(completion), "published_files": [str(path) for path in expected.values()]}
            if ctx is not None:
                ctx["caldera"] = result
            logger.record("SKIP", "caldera_run", reason="validated_complete_outputs")
            return result
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
                    credible_metrics = convert_finemap_credible_sets(pipeline_directory, staged_credible_sets, module)
                else:
                    staged_credible_sets = resources["credible_sets"]
                    credible_metrics = resources["credible_set_metrics"]
                rscript = resolve_executable(configuration.resources.executables.rscript, "Rscript", error_type=CalderaError)
                assembly = module.assembly_by_genome_build[module.genome_build]
                run_checked_command(
                    [
                        rscript, resources["adapter"], resources["repository"],
                        resources["upstream_script"], pops, staged_credible_sets,
                        str(assembly), staged_results, module.input_schema.table_delimiter,
                    ],
                    "CALDERA", logger=logger, error_type=CalderaError,
                    timeout_seconds=configuration.execution.timeout_seconds,
                    expected_outputs=[staged_results],
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
        if ctx is not None:
            ctx["caldera"] = result
        return result
    except BaseException as exc:
        logger.error("CALDERA analysis failed: %s: %s" % (type(exc).__name__, exc))
        raise
    finally:
        logger.close()


__all__ = [
    "convert_finemap_credible_sets", "preflight_caldera_pipeline",
    "run_caldera_direct", "validate_caldera_configuration",
]
