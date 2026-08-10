"""Strict PostGWAS service boundary for the unchanged upstream FLAMES program."""

from __future__ import annotations

import argparse
import re
import tempfile
from pathlib import Path

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
    remove_empty_directories,
    require_nonempty_file,
    resolve_executable,
    validate_filename_component,
)
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.core.processes import run_checked_command
from postgwas.modules.flames.errors import FlamesError


def _resolved_configuration(args: argparse.Namespace):
    module_overrides = explicit_overrides(args, {
        "flames_genome_build": "genome_build",
        "credible_sets_directory": "credible_sets_directory",
        "magma_gene_results_file": "magma_gene_results_file",
        "magma_covariate_results_file": "magma_covariate_results_file",
        "pops_scores_file": "pops_scores_file",
        "annotation_resource_directory": "annotation_resource_directory",
        "flames_model_directory": "model_directory",
        "flames_vep_mode": "vep_mode",
        "vep_command": "vep_command",
        "vep_cache": "vep_cache",
        "flames_cadd_mode": "cadd_mode",
        "cadd_file": "cadd_file",
    })
    global_overrides = explicit_overrides(args, {
        "dataset_id": "run.dataset_id",
        "output_directory": "run.output_directory",
        "threads": "execution.threads",
        "memory_gb": "execution.memory_gb",
        "seed": "execution.random_seed",
        "resume": "run.resume",
        "overwrite": "run.overwrite",
        "tabix": "resources.executables.tabix",
    })
    return load_run_configuration_for_module(
        "flames", getattr(args, "run_config", None),
        module_overrides=module_overrides,
        global_overrides=global_overrides,
    )


def _read_table(
    path: Path,
    label: str,
    separator: str,
    *,
    comment=None,
    allow_empty: bool = False,
) -> pd.DataFrame:
    try:
        table = pd.read_csv(
            path,
            sep=separator,
            comment=comment,
            dtype=str,
            engine="python" if len(separator) != 1 else "c",
        )
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise FlamesError("Cannot read %s %s: %s" % (label, path, exc)) from exc
    if table.empty and not allow_empty:
        raise FlamesError("%s contains no data rows: %s" % (label, path))
    return table


def _require_columns(table: pd.DataFrame, columns, label: str) -> None:
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise FlamesError(
            "%s is missing required columns: %s" % (label, ", ".join(missing))
        )


def _module_resource(root: Path, relative: str, label: str) -> Path:
    candidate = (root / relative).resolve()
    if root != candidate and root not in candidate.parents:
        raise FlamesError("%s leaves the FLAMES module directory" % label)
    return candidate


def _validate_runtime(configuration, module, logger=None) -> str:
    python = resolve_executable(
        configuration.resources.executables.python,
        "Python",
        error_type=FlamesError,
    )
    imports = ",".join(module.upstream.runtime_imports)
    run_checked_command(
        [python, "-c", "import %s" % imports],
        "FLAMES Python dependency preflight",
        logger=logger,
        error_type=FlamesError,
        timeout_seconds=configuration.execution.timeout_seconds,
    )
    return python


def _validate_upstream_resources(configuration, logger=None) -> dict:
    module = configuration.modules.flames
    package_root = Path(__file__).parent.resolve()
    script = require_nonempty_file(
        _module_resource(
            package_root, module.upstream.script_relative_path,
            "FLAMES upstream script",
        ),
        "FLAMES upstream script",
        error_type=FlamesError,
    )
    if module.model_directory is None:
        model_directory = _module_resource(
            package_root,
            module.upstream.model_directory_relative_path,
            "FLAMES model directory",
        )
    else:
        model_directory = Path(module.model_directory).expanduser().resolve()
    if not model_directory.is_dir():
        raise FlamesError("FLAMES model directory does not exist: %s" % model_directory)
    model = require_nonempty_file(
        model_directory / module.upstream.model_file,
        "FLAMES model",
        error_type=FlamesError,
    )
    features = require_nonempty_file(
        model_directory / module.upstream.feature_file,
        "FLAMES model feature manifest",
        error_type=FlamesError,
    )
    if module.annotation_resource_directory is None:
        raise FlamesError(
            "FLAMES annotation_resource_directory is required; set it in YAML or "
            "provide --flames-annotation-directory"
        )
    annotation_directory = Path(
        module.annotation_resource_directory
    ).expanduser().resolve()
    if not annotation_directory.is_dir():
        raise FlamesError(
            "FLAMES annotation resource directory does not exist: %s"
            % annotation_directory
        )
    missing_directories = []
    empty_directories = []
    for relative in module.upstream.required_annotation_directories:
        required_directory = annotation_directory / relative
        if not required_directory.is_dir():
            missing_directories.append(relative)
        elif next(
            (entry for entry in required_directory.rglob("*") if entry.is_file()),
            None,
        ) is None:
            empty_directories.append(relative)
    if missing_directories:
        raise FlamesError(
            "FLAMES annotation bundle %s is missing configured directories: %s"
            % (
                module.upstream.annotation_bundle_id,
                ", ".join(missing_directories),
            )
        )
    if empty_directories:
        raise FlamesError(
            "FLAMES annotation bundle %s has configured resource directories "
            "without any files: %s"
            % (
                module.upstream.annotation_bundle_id,
                ", ".join(empty_directories),
            )
        )
    annotation_files = {}
    for pattern in module.upstream.required_annotation_file_patterns:
        relative = pattern.format(genome_build=module.genome_build.value.upper())
        annotation_files[relative] = require_nonempty_file(
            annotation_directory / relative,
            "FLAMES annotation resource",
            error_type=FlamesError,
        )
    vep_command = None
    vep_cache = None
    if module.vep_mode == "local":
        vep_command = resolve_executable(
            module.vep_command, "VEP", error_type=FlamesError,
        )
        vep_cache = Path(module.vep_cache).expanduser().resolve()
        if not vep_cache.is_dir():
            raise FlamesError("VEP cache directory does not exist: %s" % vep_cache)
    cadd_file = None
    tabix = None
    cadd_index = None
    if module.cadd_mode == "local":
        cadd_file = require_nonempty_file(
            module.cadd_file, "CADD score file", error_type=FlamesError,
        )
        cadd_index = require_nonempty_file(
            Path(str(cadd_file) + module.upstream.cadd_index_suffix),
            "CADD tabix index",
            error_type=FlamesError,
        )
        tabix = resolve_executable(
            configuration.resources.executables.tabix,
            "tabix",
            error_type=FlamesError,
        )
    python = _validate_runtime(configuration, module, logger=logger)
    return {
        "python": python,
        "script": script,
        "model_directory": model_directory,
        "model": model,
        "features": features,
        "annotation_directory": annotation_directory,
        "annotation_files": annotation_files,
        "vep_command": vep_command,
        "vep_cache": vep_cache,
        "cadd_file": cadd_file,
        "cadd_index": cadd_index,
        "tabix": tabix,
    }


def validate_fine_mapping_index(fine_mapping_directory, module):
    """Validate the indexed PostGWAS fine-mapping interchange without inference."""
    schema = module.input_schema
    root = Path(fine_mapping_directory).expanduser().resolve()
    if not root.is_dir():
        raise FlamesError("Fine-mapping FLAMES input is not a directory: %s" % root)
    index_path = require_nonempty_file(
        root / schema.index_file,
        "Fine-mapping FLAMES index",
        error_type=FlamesError,
    )
    index = _read_table(
        index_path, "fine-mapping FLAMES index", schema.table_delimiter,
    )
    required = [
        schema.index_filename_column,
        schema.index_locus_column,
        schema.index_annotation_column,
    ]
    _require_columns(index, required, "Fine-mapping FLAMES index")
    if index[required].isna().any().any() or any(
        index[column].astype(str).str.strip().eq("").any() for column in required
    ):
        raise FlamesError("Fine-mapping FLAMES index contains empty required values")
    if index[schema.index_filename_column].duplicated().any():
        raise FlamesError(
            "Fine-mapping FLAMES index contains duplicate credible-set files"
        )
    if index[schema.index_locus_column].duplicated().any():
        raise FlamesError(
            "Fine-mapping FLAMES index must contain one credible set per locus; "
            "duplicate locus identifiers would make per-locus score normalization "
            "ambiguous"
        )
    variant_pattern = re.compile(schema.variant_identifier_pattern, re.IGNORECASE)
    credible_paths = []
    credible_metrics = []
    for row_number, filename in enumerate(index[schema.index_filename_column], 1):
        source = Path(filename).expanduser()
        if not source.is_absolute():
            source = root / source
        source = require_nonempty_file(
            source.resolve(), "Indexed credible-set file", error_type=FlamesError,
        )
        credible = _read_table(
            source,
            "credible-set file",
            schema.whitespace_delimiter_pattern,
            comment="#",
        )
        required_credible = [
            schema.credible_set_index_column,
            schema.credible_set_variant_column,
            schema.credible_set_probability_column,
        ]
        _require_columns(credible, required_credible, "Credible-set file")
        if credible[required_credible].isna().any().any():
            raise FlamesError(
                "Credible-set row identifiers, variants, and probabilities cannot "
                "be missing: %s" % source
            )
        row_identifiers = (
            credible[schema.credible_set_index_column].astype(str).str.strip()
        )
        if row_identifiers.eq("").any() or row_identifiers.duplicated().any():
            raise FlamesError(
                "Credible-set row identifiers must be non-empty and unique: %s"
                % source
            )
        variants = credible[schema.credible_set_variant_column].astype(str).str.strip()
        if variants.eq("").any() or variants.duplicated().any():
            raise FlamesError(
                "Credible-set variant identifiers must be non-empty and unique: %s"
                % source
            )
        matches = [variant_pattern.fullmatch(value) for value in variants]
        if any(match is None for match in matches):
            raise FlamesError(
                "Credible-set variants do not match the configured FLAMES "
                "identifier contract: %s"
                % source
            )
        chromosomes = np.asarray([int(match.group(1)) for match in matches])
        positions = np.asarray([int(match.group(2)) for match in matches])
        if (
            (chromosomes < schema.chromosome_minimum).any()
            or (chromosomes > schema.chromosome_maximum).any()
            or (positions < 1).any()
            or (positions > schema.position_maximum).any()
        ):
            raise FlamesError(
                "Credible-set chromosome or position is outside the configured "
                "FLAMES range: %s" % source
            )
        if len(set(chromosomes.tolist())) != 1:
            raise FlamesError(
                "Each FLAMES credible set must contain exactly one chromosome: %s"
                % source
            )
        probabilities = pd.to_numeric(
            credible[schema.credible_set_probability_column], errors="coerce",
        )
        if (
            probabilities.isna().any()
            or not np.isfinite(probabilities.to_numpy()).all()
            or ((probabilities < 0) | (probabilities > 1)).any()
        ):
            raise FlamesError(
                "Credible-set probabilities must be finite values in [0, 1]: %s"
                % source
            )
        mass = float(probabilities.sum())
        if mass + module.probability_tolerance < module.minimum_credible_set_coverage:
            raise FlamesError(
                "Credible-set cumulative PIP %.12g is below the configured %.12g: %s"
                % (mass, module.minimum_credible_set_coverage, source)
            )
        if mass > 1.0 + module.probability_tolerance:
            raise FlamesError(
                "Credible-set cumulative PIP %.12g exceeds one; PostGWAS will not "
                "permit upstream FLAMES to rescale it silently: %s" % (mass, source)
            )
        credible_paths.append(source)
        credible_metrics.append({
            "row_number": row_number,
            "variants": len(credible),
            "pip_mass": mass,
            "chromosome": int(chromosomes[0]),
        })
    return {
        "root": root,
        "index_path": index_path,
        "index": index,
        "credible_paths": credible_paths,
        "credible_metrics": credible_metrics,
    }


def _validated_gene_ids(values: pd.Series, pattern: re.Pattern, label: str) -> set[str]:
    identifiers = values.astype(str).str.strip()
    if identifiers.eq("").any() or identifiers.duplicated().any():
        raise FlamesError("%s gene identifiers must be non-empty and unique" % label)
    invalid = identifiers[~identifiers.str.fullmatch(pattern)]
    if not invalid.empty:
        raise FlamesError(
            "%s must use Ensembl gene identifiers compatible with upstream FLAMES; "
            "first invalid value: %s" % (label, invalid.iloc[0])
        )
    return set(identifiers)


def _validate_scientific_inputs(module, resources: dict) -> dict:
    schema = module.input_schema
    required_values = {
        "credible_sets_directory": module.credible_sets_directory,
        "MAGMA gene results": module.magma_gene_results_file,
        "MAGMA gene-property results": module.magma_covariate_results_file,
        "PoPS scores": module.pops_scores_file,
    }
    missing = [name for name, value in required_values.items() if value is None]
    if missing:
        raise FlamesError("Missing required FLAMES inputs: %s" % ", ".join(missing))
    interchange = validate_fine_mapping_index(module.credible_sets_directory, module)
    paths = {
        "magma": require_nonempty_file(
            module.magma_gene_results_file,
            "MAGMA gene results",
            error_type=FlamesError,
        ),
        "magma_covariate": require_nonempty_file(
            module.magma_covariate_results_file,
            "MAGMA gene-property results",
            error_type=FlamesError,
        ),
        "pops": require_nonempty_file(
            module.pops_scores_file, "PoPS scores", error_type=FlamesError,
        ),
    }
    magma = _read_table(
        paths["magma"], "MAGMA gene results", schema.whitespace_delimiter_pattern,
    )
    _require_columns(
        magma,
        [schema.magma_gene_column, schema.magma_z_column],
        "MAGMA gene results",
    )
    gene_pattern = re.compile(schema.ensembl_gene_pattern)
    magma_genes = _validated_gene_ids(
        magma[schema.magma_gene_column], gene_pattern, "MAGMA"
    )
    magma_z = pd.to_numeric(magma[schema.magma_z_column], errors="coerce")
    if magma_z.isna().any() or not np.isfinite(magma_z.to_numpy()).all():
        raise FlamesError("MAGMA Z statistics must all be finite")
    covariate = _read_table(
        paths["magma_covariate"], "MAGMA gene-property results",
        schema.whitespace_delimiter_pattern, comment="#",
    )
    covariate_columns = [schema.magma_covariate_column, schema.magma_covariate_p_column]
    _require_columns(covariate, covariate_columns, "MAGMA gene-property results")
    variables = covariate[schema.magma_covariate_column].astype(str).str.strip()
    p_values = pd.to_numeric(
        covariate[schema.magma_covariate_p_column], errors="coerce"
    )
    if variables.eq("").any() or variables.duplicated().any():
        raise FlamesError("MAGMA gene-property variables must be non-empty and unique")
    if (
        p_values.isna().any()
        or not np.isfinite(p_values.to_numpy()).all()
        or ((p_values < 0) | (p_values > 1)).any()
    ):
        raise FlamesError(
            "MAGMA gene-property P values must be finite values in [0, 1]"
        )
    pops = _read_table(paths["pops"], "PoPS scores", schema.table_delimiter)
    _require_columns(
        pops,
        [schema.pops_gene_column, schema.pops_score_column],
        "PoPS scores",
    )
    pops_genes = _validated_gene_ids(
        pops[schema.pops_gene_column], gene_pattern, "PoPS"
    )
    pops_scores = pd.to_numeric(pops[schema.pops_score_column], errors="coerce")
    if pops_scores.isna().any() or not np.isfinite(pops_scores.to_numpy()).all():
        raise FlamesError("PoPS scores must all be finite")
    shared = magma_genes & pops_genes
    if not shared:
        raise FlamesError(
            "MAGMA and PoPS contain no matching Ensembl gene identifiers; verify "
            "that both used the same gene annotation"
        )
    resources.update(paths)
    resources.update({
        "interchange": interchange,
        "magma_genes": magma_genes,
        "pops_genes": pops_genes,
        "input_metrics": {
            "credible_sets": len(interchange["index"]),
            "credible_set_variants": sum(
                item["variants"] for item in interchange["credible_metrics"]
            ),
            "magma_genes": len(magma_genes),
            "pops_genes": len(pops_genes),
            "shared_magma_pops_genes": len(shared),
            "magma_covariates": len(covariate),
        },
    })
    return resources


def _feature_names(path: Path) -> list[str]:
    names = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not names or len(names) != len(set(names)):
        raise FlamesError(
            "FLAMES feature manifest must contain unique non-empty names: %s" % path
        )
    return names


def _validate_annotations(paths, resources: dict, module) -> dict:
    schema = module.result_schema
    features = _feature_names(resources["features"])
    required = [schema.gene_column, schema.symbol_column, *features]
    genes = 0
    zero_features = set(features)
    for path in paths:
        table = _read_table(
            path, "FLAMES annotated locus", module.input_schema.table_delimiter
        )
        _require_columns(table, required, "FLAMES annotated locus")
        gene_ids = table[schema.gene_column].astype(str).str.strip()
        symbols = table[schema.symbol_column].astype(str).str.strip()
        if gene_ids.eq("").any() or symbols.eq("").any() or gene_ids.duplicated().any():
            raise FlamesError(
                "FLAMES annotated loci require unique non-empty genes and symbols: %s"
                % path
            )
        missing_magma = sorted(set(gene_ids) - resources["magma_genes"])
        missing_pops = sorted(set(gene_ids) - resources["pops_genes"])
        if missing_magma or missing_pops:
            raise FlamesError(
                "Annotated FLAMES genes are absent from compatible upstream inputs "
                "(%d missing MAGMA, %d missing PoPS): %s"
                % (len(missing_magma), len(missing_pops), path)
            )
        numeric = table[features].apply(pd.to_numeric, errors="coerce")
        if numeric.isna().any().any() or not np.isfinite(numeric.to_numpy()).all():
            raise FlamesError("FLAMES model features must all be finite: %s" % path)
        zero_features &= {
            column
            for column in features
            if np.allclose(numeric[column].to_numpy(), 0.0)
        }
        genes += len(table)
    return {
        "annotation_files": len(paths),
        "annotated_genes": genes,
        "features": len(features),
        "features_zero_in_every_locus": sorted(zero_features),
    }


def _validate_scores(raw_path: Path, prediction_path: Path, module) -> dict:
    schema = module.result_schema
    raw = _read_table(
        raw_path, "FLAMES raw scores", module.input_schema.table_delimiter
    )
    required_raw = list(schema.model_dump().values())
    _require_columns(raw, required_raw, "FLAMES raw scores")
    identity = [schema.locus_column, schema.gene_column]
    text_columns = [
        schema.filename_column,
        schema.locus_column,
        schema.symbol_column,
        schema.gene_column,
    ]
    if raw[text_columns].isna().any().any() or any(
        raw[column].astype(str).str.strip().eq("").any()
        for column in text_columns
    ):
        raise FlamesError("FLAMES raw score identifiers must be non-empty")
    if raw.duplicated(identity).any():
        raise FlamesError("FLAMES raw scores require unique locus/gene pairs")
    numeric_columns = [
        schema.xgboost_score_column, schema.pops_score_column,
        schema.raw_score_column, schema.scaled_score_column,
        schema.precision_column, schema.highest_column, schema.causal_column,
    ]
    numeric = raw[numeric_columns].apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any() or not np.isfinite(numeric.to_numpy()).all():
        raise FlamesError("FLAMES score columns must all be finite")
    if (numeric[schema.raw_score_column] < 0).any():
        raise FlamesError("FLAMES raw scores must be non-negative")
    unit_columns = [
        schema.xgboost_score_column, schema.scaled_score_column,
        schema.precision_column, schema.highest_column, schema.causal_column,
    ]
    if ((numeric[unit_columns] < 0) | (numeric[unit_columns] > 1)).any().any():
        raise FlamesError(
            "FLAMES probability, indicator, and scaled-score columns must be in "
            "[0, 1]"
        )
    for column in (schema.highest_column, schema.causal_column):
        if not numeric[column].isin([0, 1]).all():
            raise FlamesError("FLAMES %s must be binary" % column)
    sums = raw.assign(_score=numeric[schema.scaled_score_column]).groupby(
        schema.locus_column, sort=False,
    )["_score"].sum()
    if not np.allclose(sums.to_numpy(), 1.0, rtol=1e-8, atol=1e-10):
        raise FlamesError("FLAMES scaled scores must sum to one within every locus")
    raw_sums = raw.assign(_raw=numeric[schema.raw_score_column]).groupby(
        schema.locus_column, sort=False,
    )["_raw"].transform("sum")
    if (raw_sums <= 0).any():
        raise FlamesError(
            "Each FLAMES locus must have a positive total raw score for normalization"
        )
    expected_scaled = numeric[schema.raw_score_column] / raw_sums
    if not np.allclose(
        numeric[schema.scaled_score_column].to_numpy(),
        expected_scaled.to_numpy(),
        rtol=1e-8,
        atol=1e-10,
    ):
        raise FlamesError(
            "FLAMES scaled scores must equal raw score divided by the locus total"
        )
    predictions = _read_table(
        prediction_path,
        "FLAMES prioritized genes",
        module.input_schema.table_delimiter,
        allow_empty=True,
    )
    prediction_required = [
        schema.locus_column, schema.filename_column, schema.symbol_column,
        schema.gene_column, schema.scaled_score_column, schema.raw_score_column,
        schema.precision_column,
    ]
    _require_columns(predictions, prediction_required, "FLAMES prioritized genes")
    if predictions[prediction_required].isna().any().any() or any(
        predictions[column].astype(str).str.strip().eq("").any()
        for column in (
            schema.locus_column,
            schema.filename_column,
            schema.symbol_column,
            schema.gene_column,
        )
    ):
        raise FlamesError("FLAMES prioritized-gene fields must be non-empty")
    if predictions.duplicated(identity).any():
        raise FlamesError("FLAMES prioritized genes contain duplicate locus/gene pairs")
    raw_keys = set(map(tuple, raw[identity].astype(str).to_numpy()))
    prediction_keys = set(map(tuple, predictions[identity].astype(str).to_numpy()))
    if not prediction_keys.issubset(raw_keys):
        raise FlamesError("FLAMES prioritized genes are not a subset of raw scores")
    causal_keys = set(
        map(
            tuple,
            raw.loc[numeric[schema.causal_column] == 1, identity]
            .astype(str)
            .to_numpy(),
        )
    )
    if prediction_keys != causal_keys:
        raise FlamesError(
            "FLAMES prioritized genes do not match causal rows in raw scores"
        )
    if not predictions.empty:
        comparison_columns = [
            schema.filename_column,
            schema.symbol_column,
            schema.scaled_score_column,
            schema.raw_score_column,
            schema.precision_column,
        ]
        comparison = predictions.merge(
            raw[[*identity, *comparison_columns]],
            on=identity,
            how="left",
            validate="one_to_one",
            suffixes=("_prediction", "_raw"),
        )
        for column in (schema.filename_column, schema.symbol_column):
            if not comparison["%s_prediction" % column].equals(
                comparison["%s_raw" % column]
            ):
                raise FlamesError(
                    "FLAMES prioritized-gene %s values do not match raw scores"
                    % column
                )
        for column in (
            schema.scaled_score_column,
            schema.raw_score_column,
            schema.precision_column,
        ):
            prediction_values = pd.to_numeric(
                comparison["%s_prediction" % column], errors="coerce"
            )
            raw_values = pd.to_numeric(
                comparison["%s_raw" % column], errors="coerce"
            )
            if (
                prediction_values.isna().any()
                or not np.isfinite(prediction_values.to_numpy()).all()
                or not np.allclose(
                    prediction_values.to_numpy(),
                    raw_values.to_numpy(),
                    rtol=1e-8,
                    atol=1e-10,
                )
            ):
                raise FlamesError(
                    "FLAMES prioritized-gene %s values do not match raw scores"
                    % column
                )
    return {
        "scored_genes": len(raw),
        "loci": raw[schema.locus_column].nunique(),
        "prioritized_genes": len(predictions),
    }


def _rewrite_annotation_paths(path: Path, mapping: dict[str, str], module) -> int:
    table = _read_table(
        path,
        "FLAMES score output",
        module.input_schema.table_delimiter,
        allow_empty=True,
    )
    column = module.result_schema.filename_column
    _require_columns(table, [column], "FLAMES score output")
    original = table[column].astype(str)
    rewritten = original.map(mapping)
    if rewritten.isna().any():
        unknown = original[rewritten.isna()].iloc[0]
        raise FlamesError(
            "FLAMES score output references an unknown annotation file: %s" % unknown
        )
    table[column] = rewritten
    table.to_csv(path, sep=module.input_schema.table_delimiter, index=False)
    return len(table)


def _build_commands(
    configuration,
    resources,
    stage_index: Path,
    work: Path,
    score_name: str,
):
    module = configuration.modules.flames
    schema = module.input_schema
    annotation = [
        resources["python"], resources["script"], "annotate",
        "--annotation_dir", resources["annotation_directory"],
        "--pops", resources["pops"],
        "--magma_z", resources["magma"],
        "--magma_tissue", resources["magma_covariate"],
        "--indexfile", stage_index,
        "--build", module.genome_build.value,
        "--prob_col", schema.credible_set_probability_column,
        "--SNP_col", schema.credible_set_variant_column,
        "--filter", str(module.locus_window_bp),
    ]
    if module.vep_mode == "local":
        annotation.extend([
            "--cmd_vep", resources["vep_command"],
            "--vep_cache", resources["vep_cache"],
        ])
    if module.cadd_mode == "local":
        annotation.extend([
            "--tabix", resources["tabix"],
            "--CADD_file", resources["cadd_file"],
        ])
    scoring = [
        resources["python"], resources["script"], "FLAMES",
        "--indexfile", stage_index,
        "--outdir", work,
        "--filename", score_name,
        "--distance", str(module.locus_window_bp),
        "--modelpath", resources["model_directory"],
    ]
    return [str(value) for value in annotation], [str(value) for value in scoring]


def _completion_inputs(resources: dict) -> dict[str, Path]:
    interchange = resources["interchange"]
    inputs = {
        "fine_mapping_index": interchange["index_path"],
        "magma_gene_results": resources["magma"],
        "magma_gene_property_results": resources["magma_covariate"],
        "pops_scores": resources["pops"],
        "upstream_script": resources["script"],
        "model": resources["model"],
        "model_features": resources["features"],
    }
    inputs.update({
        "credible_set_%d" % number: path
        for number, path in enumerate(interchange["credible_paths"], 1)
    })
    inputs.update({
        "annotation_resource_%d" % number: path
        for number, path in enumerate(resources["annotation_files"].values(), 1)
    })
    if resources["cadd_file"] is not None:
        inputs["cadd_scores"] = resources["cadd_file"]
        inputs["cadd_index"] = resources["cadd_index"]
    return inputs


def _completion_configuration(configuration) -> str:
    return configuration_digest({
        "module": configuration.modules.flames.model_dump(mode="json"),
        "python": configuration.resources.executables.python,
        "tabix": configuration.resources.executables.tabix,
    })


def preflight_flames_pipeline(args: argparse.Namespace) -> None:
    configuration = _resolved_configuration(args)
    module = configuration.modules.flames
    configured_builds = {
        "fine_mapping": configuration.modules.fine_mapping.genome_build,
        "magma": configuration.modules.magma.genome_build,
        "pops": configuration.modules.pops.genome_build,
    }
    mismatches = {
        name: genome_build
        for name, genome_build in configured_builds.items()
        if genome_build is not None and genome_build != module.genome_build
    }
    if mismatches:
        details = ", ".join(
            "%s=%s" % (name, genome_build.value)
            for name, genome_build in mismatches.items()
        )
        raise FlamesError(
            "FLAMES, fine-mapping, MAGMA, and configured PoPS genome builds "
            "must match exactly: FLAMES=%s; %s"
            % (module.genome_build.value, details)
        )
    _validate_upstream_resources(configuration)


def run_flames_direct(args: argparse.Namespace, ctx=None):
    """Validate, execute, and publish one provenance-tracked FLAMES analysis."""
    try:
        configuration = _resolved_configuration(args)
    except BaseException as exc:
        fallback = load_configuration()
        output = Path(
            getattr(args, "output_directory", None) or fallback.run.output_directory
        ).expanduser().resolve()
        raw_dataset = getattr(args, "dataset_id", None) or fallback.run.dataset_id
        try:
            dataset = validate_filename_component(
                raw_dataset, "dataset_id", error_type=FlamesError,
            )
        except FlamesError:
            dataset = fallback.run.dataset_id
        log_path = configured_output_path(
            output,
            fallback.modules.flames.output_layout.service_log_file,
            error_type=FlamesError,
            dataset_id=dataset,
        )
        write_log_record(
            log_path,
            "ERROR",
            "FLAMES configuration validation failed: %s: %s"
            % (type(exc).__name__, exc),
            sample_id=dataset,
            file_level=fallback.logging.file_level,
            screen_level=fallback.logging.console_level,
        )
        raise
    module = configuration.modules.flames
    output = Path(configuration.run.output_directory).expanduser().resolve()
    dataset = validate_filename_component(
        configuration.run.dataset_id, "dataset_id", error_type=FlamesError,
    )
    output.mkdir(parents=True, exist_ok=True)
    log_path = configured_output_path(
        output, module.output_layout.service_log_file,
        error_type=FlamesError, dataset_id=dataset,
    )
    logger = PipelineLogger(
        dataset, "run", str(log_path.parent),
        level=configuration.logging.file_level,
        screen_level=configuration.logging.console_level,
        log_path=str(log_path),
    )
    staging_root = configured_output_path(
        output, module.output_layout.staging_directory,
        error_type=FlamesError, dataset_id=dataset,
    )
    try:
        resources = _validate_upstream_resources(configuration, logger=logger)
        resources = _validate_scientific_inputs(module, resources)
        logger.record(
            "PARAM", "flames_run",
            genome_build=module.genome_build.value,
            locus_window_bp=module.locus_window_bp,
            credible_set_coverage=module.minimum_credible_set_coverage,
            vep_mode=module.vep_mode,
            cadd_mode=module.cadd_mode,
            annotation_bundle=module.upstream.annotation_bundle_id,
        )
        logger.record("OBSERVED", "flames_inputs", **resources["input_metrics"])
        resolved_path = configured_output_path(
            output, module.output_layout.resolved_config_file,
            error_type=FlamesError, dataset_id=dataset,
        )
        write_resolved_configuration(configuration, resolved_path, modules="flames")
        logger.record("OUTPUT", "resolved_configuration", path=str(resolved_path))
        interchange = resources["interchange"]
        score_base = configured_output_path(
            output, module.output_layout.score_basename,
            error_type=FlamesError, dataset_id=dataset,
        )
        raw_path = Path(str(score_base) + module.output_layout.raw_score_suffix)
        prediction_path = Path(
            str(score_base) + module.output_layout.prediction_suffix
        )
        final_index = configured_output_path(
            output, module.output_layout.index_file,
            error_type=FlamesError, dataset_id=dataset,
        )
        final_annotations = [
            configured_output_path(
                output, module.output_layout.annotation_file,
                error_type=FlamesError, dataset_id=dataset, row_number=number,
            )
            for number in range(1, len(interchange["index"]) + 1)
        ]
        completion = configured_output_path(
            output, module.output_layout.completion_manifest,
            error_type=FlamesError, dataset_id=dataset,
        )
        expected = {
            "raw_scores": raw_path,
            "prioritized_genes": prediction_path,
            "flames_index": final_index,
        }
        expected.update({
            "annotation_%d" % number: path
            for number, path in enumerate(final_annotations, 1)
        })
        completion_inputs = _completion_inputs(resources)
        digest = _completion_configuration(configuration)
        if (
            configuration.run.resume
            and not configuration.run.overwrite
            and completion.is_file()
            and all(path.is_file() for path in expected.values())
        ):
            validate_completion_manifest(
                completion,
                dataset_id=dataset,
                module="flames",
                genome_build=module.genome_build.value,
                configuration_sha256=digest,
                inputs=completion_inputs,
                outputs=expected,
                error_type=FlamesError,
            )
            annotation_metrics = _validate_annotations(
                final_annotations, resources, module,
            )
            score_metrics = _validate_scores(raw_path, prediction_path, module)
            result = {
                "status": "success",
                "flames_raw_file": str(raw_path),
                "flames_predictions_file": str(prediction_path),
                "flames_index": str(final_index),
                "annotations": [str(path) for path in final_annotations],
                "completion_manifest": str(completion),
                "published_files": [str(path) for path in expected.values()],
                "metrics": {"annotations": annotation_metrics, "scores": score_metrics},
            }
            logger.record("SKIP", "flames_run", reason="validated_complete_outputs")
            if ctx is not None:
                ctx["flames"] = result
            return result
        existing = [path for path in expected.values() if path.exists()]
        if existing and not configuration.run.overwrite:
            raise FlamesError(
                "Existing FLAMES outputs require --overwrite or a valid resumable "
                "completion manifest: %s" % ", ".join(map(str, existing))
            )
        staging_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="run_", dir=staging_root) as directory:
            work = Path(directory)
            staged_annotations = [
                work / final.name for final in final_annotations
            ]
            staged_index = work / module.input_schema.index_file
            index = interchange["index"].copy()
            index[module.input_schema.index_filename_column] = [
                str(path) for path in interchange["credible_paths"]
            ]
            index[module.input_schema.index_annotation_column] = [
                str(path) for path in staged_annotations
            ]
            index.to_csv(
                staged_index, sep=module.input_schema.table_delimiter, index=False,
            )
            score_name = score_base.name
            annotation_command, scoring_command = _build_commands(
                configuration, resources, staged_index, work, score_name,
            )
            dry_run = bool(getattr(args, "dry_run", False))
            run_checked_command(
                annotation_command,
                "FLAMES annotation",
                logger=logger,
                error_type=FlamesError,
                timeout_seconds=configuration.execution.timeout_seconds,
                expected_outputs=staged_annotations,
                dry_run=dry_run,
            )
            if dry_run:
                run_checked_command(
                    scoring_command,
                    "FLAMES scoring",
                    logger=logger,
                    error_type=FlamesError,
                    timeout_seconds=configuration.execution.timeout_seconds,
                    dry_run=True,
                )
                logger.record("STATUS", "flames_run", status="DRY_RUN_VALIDATED")
                return {
                    "status": "dry_run",
                    "resolved_configuration": str(resolved_path),
                    "log": str(log_path),
                }
            annotation_metrics = _validate_annotations(
                staged_annotations, resources, module,
            )
            staged_raw = work / (
                score_name + module.output_layout.raw_score_suffix
            )
            staged_prediction = work / (
                score_name + module.output_layout.prediction_suffix
            )
            run_checked_command(
                scoring_command,
                "FLAMES scoring",
                logger=logger,
                error_type=FlamesError,
                timeout_seconds=configuration.execution.timeout_seconds,
                expected_outputs=[staged_raw, staged_prediction],
            )
            score_metrics = _validate_scores(
                staged_raw, staged_prediction, module,
            )
            path_mapping = {
                str(staged): str(final)
                for staged, final in zip(staged_annotations, final_annotations)
            }
            rewritten_rows = _rewrite_annotation_paths(
                staged_raw, path_mapping, module,
            ) + _rewrite_annotation_paths(staged_prediction, path_mapping, module)
            published_index = index.copy()
            published_index[module.input_schema.index_annotation_column] = [
                str(path) for path in final_annotations
            ]
            staged_published_index = work / final_index.name
            published_index.to_csv(
                staged_published_index,
                sep=module.input_schema.table_delimiter,
                index=False,
            )
            if configuration.run.overwrite:
                completion.unlink(missing_ok=True)
                for path in expected.values():
                    path.unlink(missing_ok=True)
            for staged, final in zip(staged_annotations, final_annotations):
                final.parent.mkdir(parents=True, exist_ok=True)
                staged.replace(final)
            raw_path.parent.mkdir(parents=True, exist_ok=True)
            staged_raw.replace(raw_path)
            prediction_path.parent.mkdir(parents=True, exist_ok=True)
            staged_prediction.replace(prediction_path)
            final_index.parent.mkdir(parents=True, exist_ok=True)
            staged_published_index.replace(final_index)
        result = {
            "status": "success",
            "flames_raw_file": str(raw_path),
            "flames_predictions_file": str(prediction_path),
            "flames_index": str(final_index),
            "annotations": [str(path) for path in final_annotations],
            "published_files": [str(path) for path in expected.values()],
            "metrics": {"annotations": annotation_metrics, "scores": score_metrics},
        }
        write_completion_manifest(
            completion,
            dataset_id=dataset,
            module="flames",
            genome_build=module.genome_build.value,
            configuration_sha256=digest,
            inputs=completion_inputs,
            outputs=expected,
            metrics={
                "inputs": resources["input_metrics"],
                "credible_sets": interchange["credible_metrics"],
                "annotations": annotation_metrics,
                "scores": score_metrics,
                "rewritten_provenance_rows": rewritten_rows,
            },
            error_type=FlamesError,
        )
        result["completion_manifest"] = str(completion)
        logger.record(
            "STATUS", "flames_run", status="COMPLETED",
            loci=score_metrics["loci"],
            scored_genes=score_metrics["scored_genes"],
            prioritized_genes=score_metrics["prioritized_genes"],
            annotation_files=annotation_metrics["annotation_files"],
            zero_features=annotation_metrics["features_zero_in_every_locus"],
        )
        if ctx is not None:
            ctx["flames"] = result
        return result
    except BaseException as exc:
        logger.error("FLAMES analysis failed: %s: %s" % (type(exc).__name__, exc))
        raise
    finally:
        remove_empty_directories(staging_root, staging_root.parent)
        logger.close()


__all__ = [
    "preflight_flames_pipeline",
    "run_flames_direct",
    "validate_fine_mapping_index",
]
