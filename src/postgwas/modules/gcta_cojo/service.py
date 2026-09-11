"""Validated direct and pipeline execution service for GCTA-COJO."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
import re
import shutil
import sqlite3
import tempfile
from collections.abc import Mapping
from typing import Any, get_args

import yaml

from postgwas.config import (
    load_configuration,
    load_run_configuration_for_module,
    write_resolved_configuration,
)
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.config.models.modules.gcta_cojo import (
    GctaCojoMode,
    GctaCojoParallelOutputContract,
)
from postgwas.core.completion import (
    apply_completion_restart,
    configuration_digest,
    resolve_completion_resume,
    write_completion_manifest,
)
from postgwas.core.checkpointing import software_identity
from postgwas.core.contracts import Artifact, ModuleResult
from postgwas.core.paths import (
    configured_output_path,
    require_nonempty_file,
    resolve_executable,
    validate_filename_component,
)
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.core.plink import (
    count_plink_samples,
    validate_plink_bundle_dimensions,
    validate_plink_files,
)
from postgwas.core.preflight import (
    PipelinePreflightEvidence,
    PreflightFileIdentity,
    capture_preflight_file_identities,
    pipeline_preflight_evidence,
    require_pipeline_input_vcf,
    require_unchanged_preflight_files,
)
from postgwas.core.required_arguments import RequiredArgument, require_resolved_arguments
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.modules.formatting.reference_identifiers import (
    BimIdentifierRequirement,
    configure_reference_variant_identifiers,
)
from postgwas.modules.gcta_cojo.adapters import (
    build_cojo_command,
    require_supported_gcta,
    run_cojo_command,
)
from postgwas.modules.gcta_cojo.errors import GctaCojoError
from postgwas.modules.gcta_cojo.parallel import (
    chromosome_parallel_slct_requested,
    execute_parallel_slct,
    numeric_chromosome_workloads,
    should_parallelize_slct,
)
from postgwas.modules.gcta_cojo.results import (
    normalize_cojo_results,
    write_empty_cojo_results,
)


_GCTA_FORMATTER_TARGET = "gcta_gene"
_GCTA_FORMATTER_RESULT_KEY = "summary_statistics_input_file"


@dataclass(frozen=True)
class GctaCojoReferenceValidation:
    prefix: Path
    files: dict[str, Path]
    samples: int
    executable: str
    version: str
    software: dict[str, str]


@dataclass(frozen=True)
class GctaCojoPipelineResources:
    """External COJO resources validated before formatter input creation."""

    configuration: Any
    reference: GctaCojoReferenceValidation
    identifier_observation: Mapping[str, Any]
    snp_list_metrics: Mapping[str, Mapping[str, Any]]
    file_identities: tuple[PreflightFileIdentity, ...]


_GCTA_COJO_ANALYSIS_LABELS = {
    "slct": "GCTA-COJO stepwise signal selection",
    "top_snps": "GCTA-COJO top-SNP selection",
    "joint": "GCTA-COJO specified-SNP joint analysis",
    "cond": "GCTA-COJO conditional analysis",
}

# Both files are upper-cased while entering the private validation database.
# GCTA can reverse A1/A2 coding, so compatibility is equality of the complete
# unordered allele pair rather than equality of allele order.
_ALLELE_COMPATIBILITY_SQL = (
    "((s.allele1=b.allele1 AND s.allele2=b.allele2) OR "
    "(s.allele1=b.allele2 AND s.allele2=b.allele1))"
)


def gcta_cojo_analysis_label(module) -> str:
    """Render the resolved scientific analysis without changing its mode key."""
    label = _GCTA_COJO_ANALYSIS_LABELS[module.mode]
    if module.mode == "top_snps":
        return "%s (%d SNPs requested)" % (
            label, module.analysis.top_snp_count,
        )
    return label


def gcta_cojo_pipeline_title(args) -> str:
    """Return the progress title from the same mode resolved for execution."""
    module = resolve_gcta_cojo_configuration(args).modules.gcta_cojo
    return "Run %s." % gcta_cojo_analysis_label(module)


def resolve_gcta_cojo_configuration(args):
    """Resolve CLI/YAML GCTA-COJO settings once and validate mode controls."""
    module_overrides = explicit_overrides(args, {
        "cojo_mode": "mode",
        "gcta_cojo_input_file": "input_file",
        "genome_build": "genome_build",
        "cojo_reference_prefix": "reference.prefix",
        "cojo_reference_population": "reference.population",
        "condition_snps": "inputs.condition_snps",
        "joint_snps": "inputs.joint_snps",
        "cojo_extract": "inputs.extract_snps",
        "cojo_exclude": "inputs.exclude_snps",
        "cojo_p": "analysis.significance_threshold",
        "cojo_top_snps": "analysis.top_snp_count",
        "cojo_window_kb": "analysis.window_kb",
        "cojo_collinear": "analysis.collinearity_cutoff",
        "cojo_diff_freq": "analysis.frequency_difference_max",
        "cojo_maf": "analysis.reference_maf_min",
        "cojo_gc": "analysis.genomic_control",
        "cojo_gc_lambda": "analysis.genomic_control_lambda",
        "cojo_chromosome": "analysis.chromosome",
        "cojo_minimum_reference_overlap": (
            "input_validation.minimum_reference_overlap_fraction"
        ),
        "cojo_top_results": "reporting.top_result_count",
        "cojo_finding_threshold": "reporting.finding_threshold",
        "cojo_p_value_digits": "reporting.p_value_significant_digits",
    })
    if "inputs.condition_snps" in module_overrides and "mode" not in module_overrides:
        module_overrides["mode"] = "cond"
    if "inputs.joint_snps" in module_overrides and "mode" not in module_overrides:
        module_overrides["mode"] = "joint"
    if "analysis.top_snp_count" in module_overrides and "mode" not in module_overrides:
        module_overrides["mode"] = "top_snps"
    if "analysis.genomic_control_lambda" in module_overrides:
        module_overrides["analysis.genomic_control"] = True
    global_overrides = explicit_overrides(args, {
        "dataset_id": "run.dataset_id",
        "output_directory": "run.output_directory",
        "threads": "execution.threads",
        "memory_gb": "execution.memory_gb",
        "seed": "execution.random_seed",
        "gcta": "resources.executables.gcta",
        "resume": "run.resume",
        "overwrite": "run.overwrite",
    })
    configuration = load_run_configuration_for_module(
        "gcta_cojo",
        getattr(args, "run_config", None),
        module_overrides=module_overrides,
        global_overrides=global_overrides,
    )
    mode = configuration.modules.gcta_cojo.mode
    if "analysis.significance_threshold" in module_overrides and mode != "slct":
        raise GctaCojoError(
            "--cojo-p is valid only with --cojo-mode slct; resolved mode is %s."
            % mode
        )
    if "analysis.top_snp_count" in module_overrides and mode != "top_snps":
        raise GctaCojoError(
            "--cojo-top-snps is valid only with --cojo-mode top_snps; resolved "
            "mode is %s." % mode
        )
    return configuration


def _require_gcta_cojo_arguments(
    configuration,
    input_file,
    *,
    include_generated_input: bool = True,
) -> None:
    """Report every missing direct or pipeline requirement in one error."""
    module = configuration.modules.gcta_cojo
    requirements = [
        RequiredArgument(
            "--cojo-reference-prefix",
            "modules.gcta_cojo.reference.prefix",
            module.reference.prefix,
        ),
        RequiredArgument(
            "--genome-build", "modules.gcta_cojo.genome_build",
            module.genome_build,
        ),
        RequiredArgument(
            "--cojo-reference-population",
            "modules.gcta_cojo.reference.population",
            module.reference.population,
        ),
    ]
    if include_generated_input:
        requirements.insert(0, RequiredArgument(
            "--cojo-file", "modules.gcta_cojo.input_file", input_file,
        ))
    if module.mode == "cond":
        requirements.append(RequiredArgument(
            "--condition-snps", "modules.gcta_cojo.inputs.condition_snps",
            module.inputs.condition_snps,
        ))
    if module.mode == "joint":
        requirements.append(RequiredArgument(
            "--joint-snps", "modules.gcta_cojo.inputs.joint_snps",
            module.inputs.joint_snps,
        ))
    require_resolved_arguments(requirements)


def _validate_build_and_population(module, configuration) -> None:
    if module.genome_build is None:
        raise GctaCojoError(
            "Declare --genome-build for the GWAS summary statistics and PLINK "
            "LD reference."
        )
    if module.genome_build not in configuration.resources.genomes:
        raise GctaCojoError(
            "Unknown configured genome build: %s" % module.genome_build
        )
    if module.reference.population is None:
        raise GctaCojoError(
            "Declare --cojo-reference-population; LD-reference ancestry cannot "
            "be inferred safely."
        )
    if module.reference.population not in configuration.resources.populations:
        raise GctaCojoError(
            "Unknown configured reference population: %s"
            % module.reference.population
        )


def validate_gcta_cojo_reference(
    configuration,
    *,
    logger=None,
) -> GctaCojoReferenceValidation:
    """Validate the complete reusable COJO reference and software contract."""
    module = configuration.modules.gcta_cojo
    _validate_build_and_population(module, configuration)
    if module.reference.prefix is None:
        raise GctaCojoError(
            "Required argument not provided: --cojo-reference-prefix. Provide "
            "--cojo-reference-prefix VALUE or set "
            "modules.gcta_cojo.reference.prefix in the run configuration."
        )
    prefix = Path(module.reference.prefix).expanduser().resolve()
    reference_files = validate_plink_files(
        prefix, module.reference.required_extensions, error_type=GctaCojoError,
    )
    reference_samples = count_plink_samples(
        reference_files["fam"], allow_blank_rows=True, error_type=GctaCojoError,
    )
    executable = resolve_executable(
        configuration.resources.executables.gcta,
        "GCTA executable",
        error_type=GctaCojoError,
    )
    version = require_supported_gcta(
        executable, module, logger, configuration.execution.timeout_seconds,
    )
    return GctaCojoReferenceValidation(
        prefix=prefix,
        files=reference_files,
        samples=reference_samples,
        executable=executable,
        version=version,
        software=software_identity(executable),
    )


def _validate_configured_snp_lists(module) -> dict[str, dict[str, Any]]:
    """Validate external SNP-list syntax before a pipeline creates its .ma file."""
    connection = sqlite3.connect(":memory:")
    metrics: dict[str, dict[str, Any]] = {}
    try:
        for name, value in (
            ("condition", module.inputs.condition_snps),
            ("joint", module.inputs.joint_snps),
            ("extract", module.inputs.extract_snps),
            ("exclude", module.inputs.exclude_snps),
        ):
            if value is None:
                continue
            path = require_nonempty_file(
                value, "%s SNP list" % name, error_type=GctaCojoError,
            )
            variants = _load_snp_list(connection, name, path, module)
            metrics[name] = {"path": str(path), "variants": variants}
    finally:
        connection.close()
    return metrics


def preflight_gcta_cojo_pipeline(
    args,
    *,
    preflight_evidence=None,
) -> PipelinePreflightEvidence:
    """Validate entry data and all external COJO resources before formatting."""
    entry_vcf = require_pipeline_input_vcf(preflight_evidence)
    configuration = resolve_gcta_cojo_configuration(args)
    module = configuration.modules.gcta_cojo
    _require_gcta_cojo_arguments(
        configuration,
        None,
        include_generated_input=False,
    )
    observed_build = str(entry_vcf["harmonised"]["genome_build"])
    if str(module.genome_build) != observed_build:
        raise GctaCojoError(
            "The harmonised GWAS-VCF declares genome build %s, but GCTA-COJO "
            "resolves to %s. The GWAS and PLINK LD reference must use the "
            "same build." % (observed_build, module.genome_build)
        )
    reference = validate_gcta_cojo_reference(configuration)
    configure_reference_variant_identifiers(
        args,
        configuration.modules.formatting,
        [BimIdentifierRequirement(
            consumer=gcta_cojo_analysis_label(module),
            formatter_target=_GCTA_FORMATTER_TARGET,
            bim_file=reference.files["bim"],
            column_roles=module.reference.bim_columns,
            delimiter_pattern=module.reference.table_delimiter_pattern,
        )],
    )
    identifier_observation = dict(
        args.variant_id_observations[_GCTA_FORMATTER_TARGET]
    )
    validate_plink_bundle_dimensions(
        reference.files, variants=identifier_observation["variants"],
        samples=reference.samples, error_type=GctaCojoError,
    )
    snp_list_metrics = _validate_configured_snp_lists(module)
    resource_paths = [*reference.files.values(), reference.executable]
    resource_paths.extend(
        metric["path"] for metric in snp_list_metrics.values()
    )
    resources = GctaCojoPipelineResources(
        configuration=configuration,
        reference=reference,
        identifier_observation=identifier_observation,
        snp_list_metrics=snp_list_metrics,
        file_identities=capture_preflight_file_identities(
            resource_paths,
            error_type=GctaCojoError,
            label="GCTA-COJO resource",
        ),
    )
    return pipeline_preflight_evidence(
        "gcta_cojo",
        preflight_evidence,
        resources=resources,
        deferred_checks=(
            "Validate the formatter-created GCTA .ma table.",
            "Validate exact SNP-ID and allele-pair compatibility with the BIM.",
            "Validate configured SNP-list membership in the .ma and BIM.",
        ),
    )


def _gcta_cojo_pipeline_execution_configuration(args, resources):
    """Retain validated settings while applying the orchestrator stage path."""
    configuration = resources.configuration
    run = configuration.run.model_copy(update={
        "output_directory": Path(args.output_directory).expanduser().resolve(),
    })
    return configuration.model_copy(update={"run": run}, deep=True)


def _create_validation_database(
    output: Path,
    configuration,
) -> tuple[sqlite3.Connection, Path]:
    temporary_root = configuration.execution.temporary_directory or output
    Path(temporary_root).mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        prefix=configuration.modules.gcta_cojo.input_validation.temporary_database_prefix,
        suffix=configuration.modules.gcta_cojo.input_validation.temporary_database_suffix,
        dir=temporary_root,
        delete=False,
    )
    database_path = Path(handle.name)
    handle.close()
    connection = sqlite3.connect(database_path)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute(
        "CREATE TABLE bim (variant_id TEXT PRIMARY KEY, chromosome TEXT, "
        "allele1 TEXT, allele2 TEXT)"
    )
    connection.execute(
        "CREATE TABLE summary (variant_id TEXT PRIMARY KEY, allele1 TEXT, allele2 TEXT)"
    )
    return connection, database_path


def _insert_bim(connection, bim_file: Path, module) -> int:
    roles = {name: index for index, name in enumerate(module.reference.bim_columns)}
    expected = len(roles)
    batch, count = [], 0
    pattern = re.compile(module.reference.table_delimiter_pattern)
    with bim_file.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = pattern.split(raw.strip())
            if len(fields) != expected:
                raise GctaCojoError(
                    "PLINK BIM line %d has %d fields; expected %d."
                    % (line_number, len(fields), expected)
                )
            variant = fields[roles["variant_id"]].strip()
            chromosome = fields[roles["chromosome"]].strip()
            alleles = (
                fields[roles["allele1"]].strip().upper(),
                fields[roles["allele2"]].strip().upper(),
            )
            if (
                not variant
                or not chromosome
                or not all(alleles)
                or alleles[0] == alleles[1]
            ):
                raise GctaCojoError(
                    "PLINK BIM contains an invalid variant or allele at line %d."
                    % line_number
                )
            batch.append((variant, chromosome, *alleles))
            count += 1
            if len(batch) >= module.input_validation.sqlite_batch_size:
                try:
                    connection.executemany(
                        "INSERT INTO bim VALUES (?, ?, ?, ?)", batch,
                    )
                except sqlite3.IntegrityError as exc:
                    raise GctaCojoError(
                        "PLINK BIM contains duplicate variant identifiers near line %d."
                        % line_number
                    ) from exc
                batch.clear()
    if batch:
        try:
            connection.executemany("INSERT INTO bim VALUES (?, ?, ?, ?)", batch)
        except sqlite3.IntegrityError as exc:
            raise GctaCojoError(
                "PLINK BIM contains duplicate variant identifiers."
            ) from exc
    connection.commit()
    if count == 0:
        raise GctaCojoError("PLINK BIM contains no variants: %s" % bim_file)
    return count


def _insert_summary(connection, summary_file: Path, configuration) -> int:
    module = configuration.modules.gcta_cojo
    schema = configuration.modules.formatting.exports[
        _GCTA_FORMATTER_TARGET
    ].outputs["summary_statistics"]
    expected_header = list(schema.columns.values())
    official_header = list(module.input_validation.official_summary_header)
    roles = {
        name: index
        for index, name in enumerate(module.input_validation.summary_column_roles)
    }
    if len(expected_header) != len(roles):
        raise GctaCojoError(
            "Resolved formatter GCTA .ma schema has %d columns; expected %d "
            "configured scientific roles."
            % (len(expected_header), len(roles))
        )
    pattern = re.compile(module.input_validation.table_delimiter_pattern)
    batch, count = [], 0
    with summary_file.open("r", encoding="utf-8", errors="replace") as handle:
        header = handle.readline()
        observed_header = pattern.split(header.strip()) if header else []
        if observed_header not in (expected_header, official_header):
            raise GctaCojoError(
                "GCTA .ma header does not match the resolved schema. "
                "Expected %s or %s; observed %s."
                % (
                    " ".join(expected_header),
                    " ".join(official_header),
                    " ".join(observed_header),
                )
            )
        for line_number, raw in enumerate(handle, 2):
            if not raw.strip():
                continue
            fields = pattern.split(raw.strip())
            if len(fields) != len(expected_header):
                raise GctaCojoError(
                    "GCTA .ma line %d has %d fields; expected %d."
                    % (line_number, len(fields), len(expected_header))
                )
            variant = fields[roles["variant_id"]].strip()
            alleles = (
                fields[roles["effect_allele"]].upper(),
                fields[roles["other_allele"]].upper(),
            )
            if not variant or not all(alleles) or alleles[0] == alleles[1]:
                raise GctaCojoError(
                    "GCTA .ma contains an invalid variant or allele at line %d."
                    % line_number
                )
            numeric = {}
            for role in (
                "effect_allele_frequency", "effect", "standard_error",
                "p_value", "sample_size",
            ):
                try:
                    numeric[role] = float(fields[roles[role]])
                except ValueError as exc:
                    raise GctaCojoError(
                        "GCTA .ma %s is non-numeric at line %d."
                        % (role, line_number)
                    ) from exc
                if not math.isfinite(numeric[role]):
                    raise GctaCojoError(
                        "GCTA .ma %s is non-finite at line %d."
                        % (role, line_number)
                    )
            frequency = numeric["effect_allele_frequency"]
            if frequency <= 0 or frequency >= 1:
                raise GctaCojoError(
                    "GCTA .ma effect-allele frequency must be within (0, 1) at "
                    "line %d." % line_number
                )
            if numeric["standard_error"] <= 0:
                raise GctaCojoError(
                    "GCTA .ma standard error must be positive at line %d."
                    % line_number
                )
            if numeric["p_value"] < 0 or numeric["p_value"] > 1:
                raise GctaCojoError(
                    "GCTA .ma P-value must be within [0, 1] at line %d."
                    % line_number
                )
            if (
                numeric["sample_size"]
                < module.input_validation.minimum_summary_sample_size
            ):
                raise GctaCojoError(
                    "GCTA .ma sample size must be at least %s at line %d."
                    % (
                        module.input_validation.minimum_summary_sample_size,
                        line_number,
                    )
                )
            batch.append((variant, *alleles))
            count += 1
            if len(batch) >= module.input_validation.sqlite_batch_size:
                try:
                    connection.executemany(
                        "INSERT INTO summary VALUES (?, ?, ?)", batch,
                    )
                except sqlite3.IntegrityError as exc:
                    raise GctaCojoError(
                        "GCTA .ma contains duplicate variant IDs near line %d."
                        % line_number
                    ) from exc
                batch.clear()
    if batch:
        try:
            connection.executemany("INSERT INTO summary VALUES (?, ?, ?)", batch)
        except sqlite3.IntegrityError as exc:
            raise GctaCojoError(
                "GCTA .ma contains duplicate variant identifiers."
            ) from exc
    connection.commit()
    if count == 0:
        raise GctaCojoError("GCTA .ma contains no variants: %s" % summary_file)
    return count


def _load_snp_list(connection, name: str, path: Path, module) -> int:
    connection.execute("CREATE TABLE %s (variant_id TEXT PRIMARY KEY)" % name)
    pattern = re.compile(module.input_validation.table_delimiter_pattern)
    batch, count = [], 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = pattern.split(raw.strip())
            if len(fields) != 1 or not fields[0]:
                raise GctaCojoError(
                    "%s line %d must contain exactly one SNP identifier."
                    % (path, line_number)
                )
            batch.append((fields[0],))
            count += 1
            if len(batch) >= module.input_validation.sqlite_batch_size:
                try:
                    connection.executemany(
                        "INSERT INTO %s VALUES (?)" % name, batch,
                    )
                except sqlite3.IntegrityError as exc:
                    raise GctaCojoError(
                        "SNP list contains duplicate identifiers near line %d: %s"
                        % (line_number, path)
                    ) from exc
                batch.clear()
    if batch:
        try:
            connection.executemany("INSERT INTO %s VALUES (?)" % name, batch)
        except sqlite3.IntegrityError as exc:
            raise GctaCojoError(
                "SNP list contains duplicate identifiers: %s" % path
            ) from exc
    connection.commit()
    if count == 0:
        raise GctaCojoError("SNP list contains no identifiers: %s" % path)
    return count


def _allele_mismatch_examples(connection, limit: int) -> list[dict[str, object]]:
    rows = connection.execute(
        "SELECT s.variant_id, s.allele1, s.allele2, b.allele1, b.allele2 "
        "FROM summary s INNER JOIN bim b USING (variant_id) WHERE NOT %s "
        "ORDER BY s.rowid LIMIT ?" % _ALLELE_COMPATIBILITY_SQL,
        (limit,),
    ).fetchall()
    return [
        {
            "variant_id": variant,
            "summary_alleles": [summary_a1, summary_a2],
            "bim_alleles": [bim_a1, bim_a2],
        }
        for variant, summary_a1, summary_a2, bim_a1, bim_a2 in rows
    ]


def _format_allele_mismatch_examples(examples) -> str:
    return "; ".join(
        "%s (.ma %s/%s; BIM %s/%s)"
        % (
            example["variant_id"],
            *example["summary_alleles"],
            *example["bim_alleles"],
        )
        for example in examples
    )


def _validate_snp_list_membership(
    connection,
    table: str,
    path: Path,
    *,
    strict: bool,
):
    missing_summary = connection.execute(
        "SELECT COUNT(*) FROM %s l LEFT JOIN summary s USING (variant_id) "
        "WHERE s.variant_id IS NULL" % table
    ).fetchone()[0]
    missing_reference = connection.execute(
        "SELECT COUNT(*) FROM %s l LEFT JOIN bim b USING (variant_id) "
        "WHERE b.variant_id IS NULL" % table
    ).fetchone()[0]
    incompatible_alleles = connection.execute(
        "SELECT COUNT(*) FROM %s l INNER JOIN summary s USING (variant_id) "
        "INNER JOIN bim b USING (variant_id) WHERE NOT %s"
        % (table, _ALLELE_COMPATIBILITY_SQL)
    ).fetchone()[0]
    if strict and (missing_summary or missing_reference or incompatible_alleles):
        raise GctaCojoError(
            "SNP list %s is incompatible with the GCTA .ma and LD reference: "
            "%d IDs are absent from the .ma file, %d are absent from the BIM, "
            "and %d have incompatible .ma/BIM allele pairs. Every conditioning, "
            "joint, or extract SNP must be usable by GCTA."
            % (
                path,
                missing_summary,
                missing_reference,
                incompatible_alleles,
            )
        )
    return {
        "variants": connection.execute(
            "SELECT COUNT(*) FROM %s" % table
        ).fetchone()[0],
        "missing_summary": missing_summary,
        "missing_reference": missing_reference,
        "incompatible_alleles": incompatible_alleles,
    }


def _write_effective_exclusion(connection, output: Path, module) -> tuple[Path, int]:
    """Write the deterministic union of user and allele-mismatch exclusions."""
    validation = module.input_validation
    handle = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        prefix=validation.temporary_exclusion_prefix,
        suffix=validation.temporary_exclusion_suffix,
        dir=output,
        delete=False,
    )
    path = Path(handle.name)
    mismatch_query = (
        "SELECT s.variant_id FROM summary s INNER JOIN bim b USING (variant_id) "
        "WHERE NOT %s" % _ALLELE_COMPATIBILITY_SQL
    )
    if module.inputs.exclude_snps is not None:
        query = (
            "%s UNION SELECT variant_id FROM exclude ORDER BY variant_id"
            % mismatch_query
        )
    else:
        query = mismatch_query + " ORDER BY s.variant_id"
    count = 0
    try:
        with handle:
            for (variant,) in connection.execute(query):
                handle.write(variant + "\n")
                count += 1
    except BaseException:
        path.unlink(missing_ok=True)
        raise
    if count == 0:
        path.unlink(missing_ok=True)
        raise GctaCojoError(
            "Allele-mismatch exclusion generation produced no SNP identifiers."
        )
    return path, count


def _eligible_chromosome_counts(connection, module) -> dict[str, int]:
    """Count GCTA-usable selection candidates without discarding raw labels."""
    joins = [
        "summary s INNER JOIN bim b USING (variant_id)",
    ]
    conditions = [_ALLELE_COMPATIBILITY_SQL]
    if module.inputs.extract_snps is not None:
        joins.append("INNER JOIN extract x USING (variant_id)")
    if module.inputs.exclude_snps is not None:
        joins.append("LEFT JOIN exclude e USING (variant_id)")
        conditions.append("e.variant_id IS NULL")
    rows = connection.execute(
        "SELECT b.chromosome, COUNT(*) FROM %s WHERE %s "
        "GROUP BY b.chromosome ORDER BY b.chromosome"
        % (" ".join(joins), " AND ".join(conditions))
    ).fetchall()
    return {str(chromosome): int(count) for chromosome, count in rows}


def _validate_inputs(
    summary_file: Path,
    reference_files: dict[str, Path],
    output: Path,
    configuration,
) -> tuple[dict, list[str], Path | None]:
    module = configuration.modules.gcta_cojo
    connection, database_path = _create_validation_database(output, configuration)
    warnings = []
    try:
        bim_variants = _insert_bim(connection, reference_files["bim"], module)
        summary_variants = _insert_summary(connection, summary_file, configuration)
        reference_id_overlap = connection.execute(
            "SELECT COUNT(*) FROM summary INNER JOIN bim USING (variant_id)"
        ).fetchone()[0]
        incompatible_alleles = connection.execute(
            "SELECT COUNT(*) FROM summary s INNER JOIN bim b USING (variant_id) "
            "WHERE NOT %s" % _ALLELE_COMPATIBILITY_SQL
        ).fetchone()[0]
        reference_overlap = reference_id_overlap - incompatible_alleles
        missing_reference = summary_variants - reference_id_overlap
        overlap_fraction = reference_overlap / summary_variants
        id_overlap_fraction = reference_id_overlap / summary_variants
        mismatch_examples = _allele_mismatch_examples(
            connection,
            module.input_validation.maximum_allele_mismatch_examples,
        )
        mismatch_example_text = _format_allele_mismatch_examples(
            mismatch_examples,
        )
        if overlap_fraction < module.input_validation.minimum_reference_overlap_fraction:
            raise GctaCojoError(
                "Only %.2f%% of GCTA .ma variants are usable with the PLINK BIM; "
                "the configured minimum is %.2f%%. Usable means that both the SNP "
                "ID and complete allele pair match. Breakdown — usable: %d; IDs "
                "absent from the BIM: %d; shared variant IDs with incompatible "
                "alleles: %d.%s "
                "Check --genome-build, the LD-reference prefix, variant identifiers, "
                "and reference preparation."
                % (
                    100 * overlap_fraction,
                    100 * module.input_validation.minimum_reference_overlap_fraction,
                    reference_overlap,
                    missing_reference,
                    incompatible_alleles,
                    (
                        " Examples: %s." % mismatch_example_text
                        if mismatch_example_text else ""
                    ),
                )
            )
        if missing_reference:
            warnings.append(
                "GWAS variants absent from the LD reference: %s of %s "
                "formatted variants are missing from the BIM and cannot enter "
                "the COJO genotype/LD model."
                % (
                    f"{missing_reference:,}",
                    f"{summary_variants:,}",
                )
            )
        if incompatible_alleles:
            warnings.append(
                "GWAS/reference allele mismatch: %s of %s variants with shared "
                "IDs have incompatible allele pairs. PostGWAS excludes them "
                "from the GCTA genotype/LD model while retaining the complete "
                "formatted GWAS input for phenotypic-variance estimation.%s"
                % (
                    f"{incompatible_alleles:,}",
                    f"{reference_id_overlap:,}",
                    (
                        " Examples: %s." % mismatch_example_text
                        if mismatch_example_text else ""
                    ),
                )
            )

        list_metrics = {}
        list_fields = (
            ("condition", module.inputs.condition_snps, True),
            ("joint", module.inputs.joint_snps, True),
            ("extract", module.inputs.extract_snps, True),
            ("exclude", module.inputs.exclude_snps, False),
        )
        for name, value, strict in list_fields:
            if value is None:
                continue
            path = require_nonempty_file(
                value, "%s SNP list" % name, error_type=GctaCojoError,
            )
            _load_snp_list(connection, name, path, module)
            observed = _validate_snp_list_membership(
                connection, name, path, strict=strict,
            )
            list_metrics[name] = {"path": str(path), **observed}
            if not strict and (
                observed["missing_summary"] or observed["missing_reference"]
            ):
                warnings.append(
                    "Configured exclusion IDs with no effect: %s are absent "
                    "from the formatted GWAS input and %s are absent from the "
                    "BIM; IDs absent from the relevant input cannot remove an "
                    "additional COJO candidate."
                    % (
                        f"{observed['missing_summary']:,}",
                        f"{observed['missing_reference']:,}",
                    )
                )
        if module.inputs.condition_snps is not None:
            if module.inputs.exclude_snps is not None:
                excluded_conditioning = connection.execute(
                    "SELECT COUNT(*) FROM condition INNER JOIN exclude "
                    "USING (variant_id)"
                ).fetchone()[0]
                if excluded_conditioning:
                    raise GctaCojoError(
                        "%d conditioning SNPs also occur in --cojo-exclude. "
                        "Conditioning SNPs must not be excluded."
                        % excluded_conditioning
                    )
            if module.inputs.extract_snps is not None:
                missing = connection.execute(
                    "SELECT COUNT(*) FROM condition c LEFT JOIN extract e "
                    "USING (variant_id) WHERE e.variant_id IS NULL"
                ).fetchone()[0]
                if missing:
                    raise GctaCojoError(
                        "%d conditioning SNPs are absent from --cojo-extract. "
                        "Include every conditioning SNP so GCTA retains the "
                        "complete conditioning model." % missing
                    )
        if (
            module.inputs.joint_snps is not None
            and module.inputs.exclude_snps is not None
        ):
            excluded_joint = connection.execute(
                "SELECT COUNT(*) FROM joint INNER JOIN exclude USING (variant_id)"
            ).fetchone()[0]
            if excluded_joint:
                raise GctaCojoError(
                    "%d joint-model SNPs also occur in --cojo-exclude. Joint "
                    "SNPs must not be excluded." % excluded_joint
                )
        effective_exclusion = None
        effective_exclusion_variants = (
            list_metrics.get("exclude", {}).get("variants", 0)
        )
        if incompatible_alleles:
            effective_exclusion, effective_exclusion_variants = (
                _write_effective_exclusion(connection, output, module)
            )
        eligible_chromosome_counts = _eligible_chromosome_counts(
            connection, module,
        )
        return {
            "summary_variants": summary_variants,
            "reference_variants": bim_variants,
            "reference_id_overlap_variants": reference_id_overlap,
            "reference_id_overlap_fraction": id_overlap_fraction,
            "reference_overlap_variants": reference_overlap,
            "reference_overlap_fraction": overlap_fraction,
            "reference_missing_variants": missing_reference,
            "reference_allele_mismatch_variants": incompatible_alleles,
            "reference_allele_mismatch_examples": mismatch_examples,
            "effective_exclude_variants": effective_exclusion_variants,
            "eligible_chromosome_variant_counts": eligible_chromosome_counts,
            "snp_lists": list_metrics,
        }, warnings, effective_exclusion
    finally:
        connection.close()
        database_path.unlink(missing_ok=True)


def _staged_path(final_prefix: Path, final_path: Path, staged_prefix: Path) -> Path:
    suffix = final_path.name[len(final_prefix.name):]
    return staged_prefix.with_name(staged_prefix.name + suffix)


def _publish_staged_outputs(
    staged_prefix: Path,
    final_prefix: Path,
    overwrite: bool,
) -> list[Path]:
    sources = sorted(staged_prefix.parent.glob(staged_prefix.name + "*"))
    if not sources:
        raise GctaCojoError("GCTA-COJO produced no output files.")
    destinations = [
        final_prefix.with_name(final_prefix.name + source.name[len(staged_prefix.name):])
        for source in sources
    ]
    existing = [str(path) for path in destinations if path.exists()]
    if existing and not overwrite:
        raise GctaCojoError(
            "Output already exists: %s. Use --resume or --overwrite."
            % ", ".join(existing)
        )
    final_prefix.parent.mkdir(parents=True, exist_ok=True)
    for source, destination in zip(sources, destinations):
        if destination.exists():
            destination.unlink()
        source.replace(destination)
    return destinations


def _input_manifest_paths(summary_file, reference_files, module):
    inputs = {"summary_statistics": summary_file, **{
        "reference_%s" % name: path for name, path in reference_files.items()
    }}
    for name, value in (
        ("condition_snps", module.inputs.condition_snps),
        ("joint_snps", module.inputs.joint_snps),
        ("extract_snps", module.inputs.extract_snps),
        ("exclude_snps", module.inputs.exclude_snps),
    ):
        if value is not None:
            inputs[name] = Path(value).expanduser().resolve()
    return inputs


def _read_gcta_warnings(log_path: Path, module) -> list[str]:
    if not log_path.is_file():
        return []
    pattern = re.compile(module.log_parsing.warning_pattern)
    warnings = []
    with log_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            message = raw.strip()
            if message and pattern.search(message) and message not in warnings:
                warnings.append(message)
                if len(warnings) >= module.log_parsing.maximum_warning_messages:
                    break
    return warnings


def _gcta_reported_no_signals(log_path: Path, module) -> bool:
    if not log_path.is_file():
        return False
    pattern = re.compile(module.log_parsing.no_signals_pattern)
    with log_path.open("r", encoding="utf-8", errors="replace") as handle:
        return any(pattern.search(line) for line in handle)


def _read_validated_snp_list(path: str | Path, module) -> set[str]:
    pattern = re.compile(module.input_validation.table_delimiter_pattern)
    identifiers = set()
    with Path(path).expanduser().resolve().open(
        "r", encoding="utf-8", errors="replace",
    ) as handle:
        for raw in handle:
            if not raw.strip():
                continue
            identifiers.add(pattern.split(raw.strip())[0])
    return identifiers


def _read_output_snp_ids(path: Path, module) -> set[str]:
    path = require_nonempty_file(
        path, "GCTA modeled-SNP output", error_type=GctaCojoError,
    )
    pattern = re.compile(module.input_validation.table_delimiter_pattern)
    identifier_column = module.results.schemas[module.mode].identifier_column
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        header = handle.readline()
        columns = pattern.split(header.strip()) if header else []
        if identifier_column not in columns:
            raise GctaCojoError(
                "GCTA modeled-SNP output %s lacks identifier column %s."
                % (path, identifier_column)
            )
        identifier_index = columns.index(identifier_column)
        identifiers = set()
        rows = 0
        for line_number, raw in enumerate(handle, 2):
            if not raw.strip():
                continue
            fields = pattern.split(raw.strip())
            if len(fields) != len(columns):
                raise GctaCojoError(
                    "GCTA modeled-SNP output %s line %d has %d fields; expected %d."
                    % (path, line_number, len(fields), len(columns))
                )
            identifier = fields[identifier_index]
            if identifier in identifiers:
                raise GctaCojoError(
                    "GCTA modeled-SNP output contains duplicate SNP ID %s: %s"
                    % (identifier, path)
                )
            identifiers.add(identifier)
            rows += 1
    if rows == 0:
        raise GctaCojoError("GCTA modeled-SNP output contains no SNPs: %s" % path)
    return identifiers


def _validate_modeled_snp_survival(
    requested_path: str | Path,
    observed_path: Path,
    module,
    label: str,
) -> None:
    requested = _read_validated_snp_list(requested_path, module)
    observed = _read_output_snp_ids(observed_path, module)
    missing = requested - observed
    unexpected = observed - requested
    if missing or unexpected:
        raise GctaCojoError(
            "GCTA %s SNP set does not exactly match the requested list: %d "
            "requested SNPs are missing and %d unexpected SNPs were returned. "
            "Check --cojo-maf, --cojo-diff-freq, --cojo-chromosome, extract/"
            "exclude lists, alleles, and reference genotypes."
            % (label, len(missing), len(unexpected))
        )


def _completion_records_no_signals(path: Path) -> bool:
    try:
        document = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise GctaCojoError(
            "Cannot read GCTA-COJO completion manifest %s: %s" % (path, exc)
        ) from exc
    metrics = document.get("metrics")
    return isinstance(metrics, dict) and metrics.get("no_signals") is True


def _expected_output_contract(
    module,
    paths: dict[str, Path],
    *,
    no_signals: bool,
    allele_mismatches: int,
    include_conditional_results: bool = True,
):
    outputs = {
        "normalized_result": paths["normalized_result"],
        "gcta_log": paths["gcta_log"],
    }
    if no_signals:
        if allele_mismatches:
            outputs["allele_mismatch_exclusion"] = paths[
                "allele_mismatch_exclusion"
            ]
        return outputs
    outputs["primary_result"] = paths["primary_result"]
    if module.mode != "cond":
        outputs["ld_matrix"] = paths["ld_matrix"]
    if (
        module.mode in {"slct", "top_snps"}
        and include_conditional_results
    ):
        outputs["conditional_results"] = paths["conditional_results"]
    if module.mode == "cond":
        outputs["condition_snps"] = paths["condition_snps"]
    if allele_mismatches:
        outputs["allele_mismatch_exclusion"] = paths[
            "allele_mismatch_exclusion"
        ]
    return outputs


def _publish_completed_execution(
    *,
    staging: Path,
    staged_prefix: Path,
    final_prefix: Path,
    staged_normalized: Path,
    normalized: Path,
    output_paths: dict[str, Path],
    module,
    configuration,
    input_metrics: dict,
    no_signals: bool,
    result_metrics: dict,
    include_conditional_results: bool,
    completion: Path,
    dataset: str,
    digest: str,
    manifest_inputs: dict,
    step,
) -> list[Path]:
    """Atomically publish one validated single or chromosome-combined run."""
    published = _publish_staged_outputs(
        staged_prefix, final_prefix, configuration.run.overwrite,
    )
    normalized.parent.mkdir(parents=True, exist_ok=True)
    if normalized.exists():
        if not configuration.run.overwrite:
            raise GctaCojoError("Output already exists: %s" % normalized)
        normalized.unlink()
    staged_normalized.replace(normalized)
    result_metrics["normalized_result"] = str(normalized)
    step.set_rows(result_metrics["result_rows"], removed=0)
    step.output("cojo_outputs", files=[str(path) for path in published])
    step.output("normalized_result", path=str(normalized))
    shutil.rmtree(staging)
    expected_outputs = _expected_output_contract(
        module,
        output_paths,
        no_signals=no_signals,
        allele_mismatches=input_metrics["reference_allele_mismatch_variants"],
        include_conditional_results=include_conditional_results,
    )
    missing_expected = [
        str(path) for path in expected_outputs.values()
        if not path.is_file() or path.stat().st_size <= 0
    ]
    if missing_expected:
        raise GctaCojoError(
            "GCTA-COJO completed but expected outputs are missing or empty: %s"
            % ", ".join(missing_expected)
        )
    write_completion_manifest(
        completion,
        dataset_id=dataset,
        module="gcta_cojo",
        genome_build=module.genome_build,
        configuration_sha256=digest,
        inputs=manifest_inputs,
        outputs=expected_outputs,
        metrics=result_metrics,
        error_type=GctaCojoError,
    )
    step.output("completion_manifest", path=str(completion))
    return published


def _format_p(value, digits):
    return ("%.*g" % (digits, float(value))) if value is not None else "NA"


def _render_success_summary(
    dataset,
    module,
    metrics,
    artifacts,
    warnings,
    log_path,
    label_width,
) -> str:
    lines = ["", screen_line("analysis", "GCTA-COJO completed", indent=2)]
    fields = [
        ("info", "Dataset", dataset),
        ("analysis", "COJO analysis", gcta_cojo_analysis_label(module)),
        ("genetic", "LD population", module.reference.population),
        ("count", "Reference samples", "%s" % f"{metrics['reference_samples']:,}"),
        (
            "genetic", "BIM identifier type",
            metrics["variant_id_type"],
        ),
        (
            "count", "GCTA-usable GWAS/reference variants",
            "%s/%s variants (%.2f%%)"
            % (
                f"{metrics['reference_overlap_variants']:,}",
                f"{metrics['summary_variants']:,}",
                100 * metrics["reference_overlap_fraction"],
            ),
        ),
    ]
    if metrics.get("reference_missing_variants"):
        fields.append((
            "warning",
            "GWAS variants absent from BIM",
            f"{metrics['reference_missing_variants']:,}",
        ))
    if metrics.get("reference_allele_mismatch_variants"):
        fields.append((
            "warning",
            "Shared IDs with incompatible alleles",
            f"{metrics['reference_allele_mismatch_variants']:,}",
        ))
    if metrics.get("execution_strategy") == "chromosome_parallel":
        fields.append((
            "analysis",
            "Parallel chromosome analyses",
            "up to %d simultaneously · %d total"
            % (
                metrics["parallel_workers"],
                len(metrics["chromosomes"]),
            ),
        ))
        if metrics.get("retried_chromosomes"):
            fields.append((
                "warning",
                "Automatically retried chromosomes",
                ", ".join(
                    str(value) for value in metrics["retried_chromosomes"]
                ),
            ))
    if module.mode in {"slct", "top_snps", "joint"}:
        fields.append((
            "success", "Independent/joint signals",
            "%s reported; %s at p <= %s"
            % (
                f"{metrics['result_rows']:,}",
                f"{metrics['significant_rows']:,}",
                _format_p(
                    metrics["finding_threshold"],
                    module.reporting.p_value_significant_digits,
                ),
            ),
        ))
    else:
        fields.extend([
            (
                "success", "Conditional findings",
                "%s/%s estimable variants significant; %s newly significant"
                % (
                    f"{metrics['significant_rows']:,}",
                    f"{metrics['estimated_rows']:,}",
                    f"{metrics['newly_significant_rows']:,}",
                ),
            ),
            (
                "warning" if metrics["not_estimable_rows"] else "info",
                "Not estimable",
                "%s variants (GCTA multivariate-collinearity rule)"
                % f"{metrics['not_estimable_rows']:,}",
            ),
        ])
    for kind, label, value in fields:
        lines.append(screen_field(
            kind, label, value, indent=6, label_width=label_width,
        ))
    if metrics.get("top_findings"):
        p_column = metrics["cojo_p_value_column"]
        for index, finding in enumerate(metrics["top_findings"], 1):
            lines.append(screen_field(
                "analysis",
                "Top finding %d" % index,
                "%s (chr%s:%s; %s=%s)"
                % (
                    finding["variant_id"], finding["chromosome"], finding["position"], p_column,
                    _format_p(
                        finding["cojo_p_value"],
                        module.reporting.p_value_significant_digits,
                    ),
                ),
                indent=6,
                label_width=label_width,
            ))
    if warnings:
        lines.append(screen_field(
            "warning", "Warnings", "%d; details follow and are saved in the log" % len(warnings),
            indent=6, label_width=label_width,
        ))
        for warning in warnings:
            lines.append(screen_field(
                "warning", "Warning", warning,
                indent=6, label_width=label_width,
            ))
    else:
        lines.append(screen_field(
            "info", "Warnings", "none",
            indent=6, label_width=label_width,
        ))
    lines.append(screen_field(
        "info", "Normalized result", artifacts["normalized_results"].path,
        indent=6, label_width=label_width,
    ))
    lines.append(screen_field(
        "info", "Full PostGWAS log", log_path,
        indent=6, label_width=label_width,
    ))
    lines.append("")
    return "\n".join(lines)


def _render_failure_summary(dataset, mode, reason, log_path, label_width) -> str:
    return "\n".join([
        "",
        screen_line("error", "GCTA-COJO failed", indent=2),
        screen_field("info", "Dataset", dataset, indent=6, label_width=label_width),
        screen_field("analysis", "COJO mode", mode, indent=6, label_width=label_width),
        screen_field("error", "Reason", reason, indent=6, label_width=label_width),
        screen_field(
            "info", "Full PostGWAS log", log_path,
            indent=6, label_width=label_width,
        ),
        "",
    ])


def _return_dry_run_result(
    *,
    ctx,
    logger,
    emit_terminal_summary: bool,
    dataset: str,
    module,
    configuration,
    reference_samples: int,
    observation: dict,
    software: dict,
    version: str,
    input_metrics: dict,
    execution_metrics: dict,
    warnings: list[str],
    final_prefix: Path,
    log_path: Path,
) -> ModuleResult:
    metrics = {
        "mode": module.mode,
        "dry_run": True,
        "gcta_version": version,
        "postgwas_version": software["postgwas"],
        "python_version": software["python"],
        "reference_samples": reference_samples,
        "variant_id_type": observation["variant_id_type"],
        **input_metrics,
        **execution_metrics,
    }
    result = ModuleResult(
        "gcta_cojo",
        metrics=metrics,
        warnings=tuple(warnings),
    )
    if ctx is not None:
        ctx.publish(result)
    logger.record(
        "STATUS", "gcta_cojo", status="VALIDATED", mode=module.mode,
    )
    if emit_terminal_summary:
        print("\n".join([
            "",
            screen_line("analysis", "GCTA-COJO command validated", indent=2),
            screen_field(
                "info", "Dataset", dataset, indent=6,
                label_width=configuration.logging.terminal_label_width,
            ),
            screen_field(
                "info", "Validated output prefix", final_prefix, indent=6,
                label_width=configuration.logging.terminal_label_width,
            ),
            screen_field(
                "info", "Full PostGWAS log", log_path, indent=6,
                label_width=configuration.logging.terminal_label_width,
            ),
            "",
        ]))
    return result


def run_gcta_cojo_direct(
    args,
    ctx=None,
    *,
    configuration=None,
    emit_terminal_summary: bool = True,
    reference_validation: GctaCojoReferenceValidation | None = None,
    pipeline_resources: GctaCojoPipelineResources | None = None,
    parallel_slct_output_contract: GctaCojoParallelOutputContract = "complete",
) -> ModuleResult:
    """Validate an existing .ma artifact and run COJO in direct or pipeline mode."""
    if pipeline_resources is not None and configuration is not None:
        raise GctaCojoError(
            "Pass either pipeline_resources or configuration to GCTA-COJO, not both."
        )
    if pipeline_resources is not None and reference_validation is not None:
        raise GctaCojoError(
            "Pass either pipeline_resources or reference_validation to "
            "GCTA-COJO, not both."
        )
    try:
        configuration = (
            _gcta_cojo_pipeline_execution_configuration(args, pipeline_resources)
            if pipeline_resources is not None
            else configuration or resolve_gcta_cojo_configuration(args)
        )
    except BaseException as exc:
        fallback = load_configuration()
        output = Path(
            getattr(args, "output_directory", None) or fallback.run.output_directory
        ).expanduser().resolve()
        dataset = str(getattr(args, "dataset_id", None) or fallback.run.dataset_id)
        try:
            dataset = validate_filename_component(
                dataset, "run.dataset_id", error_type=GctaCojoError,
            )
        except GctaCojoError:
            dataset = fallback.run.dataset_id
        mode = str(
            getattr(args, "cojo_mode", None) or fallback.modules.gcta_cojo.mode
        )
        if mode not in get_args(GctaCojoMode):
            mode = fallback.modules.gcta_cojo.mode
        log_path = configured_output_path(
            output,
            fallback.modules.gcta_cojo.output_layout.log_file,
            error_type=GctaCojoError,
            dataset_id=dataset,
            mode=mode,
        )
        write_log_record(
            log_path,
            "ERROR",
            "GCTA-COJO configuration failed: %s: %s" % (type(exc).__name__, exc),
            sample_id=dataset,
            file_level=fallback.logging.file_level,
            screen_level=fallback.logging.console_level,
        )
        if emit_terminal_summary:
            print(_render_failure_summary(
                dataset, mode, str(exc), log_path,
                fallback.logging.terminal_label_width,
            ))
        raise

    module = configuration.modules.gcta_cojo
    if pipeline_resources is not None:
        require_unchanged_preflight_files(
            pipeline_resources.file_identities,
            error_type=GctaCojoError,
            label="GCTA-COJO resource",
        )
    output = Path(configuration.run.output_directory).expanduser().resolve()
    dataset = configuration.run.dataset_id
    output.mkdir(parents=True, exist_ok=True)
    log_path = configured_output_path(
        output,
        module.output_layout.log_file,
        error_type=GctaCojoError,
        dataset_id=dataset,
        mode=module.mode,
    )
    logger = PipelineLogger(
        dataset,
        "run",
        str(log_path.parent),
        level=configuration.logging.file_level,
        screen_level=configuration.logging.console_level,
        log_path=str(log_path),
    )
    temporary_exclusion = None
    try:
        if parallel_slct_output_contract not in get_args(
            GctaCojoParallelOutputContract
        ):
            raise GctaCojoError(
                "Invalid chromosome-parallel COJO output contract: %s. "
                "Expected one of: %s."
                % (
                    parallel_slct_output_contract,
                    ", ".join(get_args(GctaCojoParallelOutputContract)),
                )
            )
        input_file = getattr(args, "gcta_cojo_input_file", None) or module.input_file
        formatter_input = None
        if ctx is not None:
            formatter = ctx.get("formatter", {})
            formatter_input = formatter.get(_GCTA_FORMATTER_TARGET, {}).get(
                _GCTA_FORMATTER_RESULT_KEY
            )
        if input_file is None and ctx is not None:
            input_file = formatter_input
        _require_gcta_cojo_arguments(configuration, input_file)
        formatter_generated_input = (
            formatter_input is not None
            and Path(input_file).expanduser().resolve()
            == Path(formatter_input).expanduser().resolve()
        )
        reference = (
            pipeline_resources.reference
            if pipeline_resources is not None
            else reference_validation or validate_gcta_cojo_reference(
                configuration, logger=logger,
            )
        )
        reference_prefix = reference.prefix
        reference_files = reference.files
        reference_samples = reference.samples
        warnings = []
        reference_sample_warning_threshold = (
            module.input_validation.reference_sample_size_warning_threshold
        )
        if reference_samples <= reference_sample_warning_threshold:
            warnings.append(
                "LD reference sample size is below GCTA's COJO recommendation: "
                "%s samples available; more than %s recommended. Selected and "
                "conditional signals may be less stable."
                % (
                    f"{reference_samples:,}",
                    f"{reference_sample_warning_threshold:,}",
                )
            )

        executable = reference.executable
        version = reference.version
        software = reference.software

        summary_file = require_nonempty_file(
            input_file, "GCTA .ma input", error_type=GctaCojoError,
        )

        input_metrics, input_warnings, temporary_exclusion = _validate_inputs(
            summary_file, reference_files, output, configuration,
        )
        warnings.extend(input_warnings)
        chromosome_counts = input_metrics["eligible_chromosome_variant_counts"]
        parallel_slct = should_parallelize_slct(module, chromosome_counts)
        effective_output_contract = (
            parallel_slct_output_contract if parallel_slct else "complete"
        )
        include_conditional_results = (
            module.mode not in {"slct", "top_snps"}
            or effective_output_contract == "complete"
        )
        if (
            chromosome_parallel_slct_requested(module)
            and len(chromosome_counts) > 1
            and numeric_chromosome_workloads(chromosome_counts) is None
        ):
            warnings.append(
                "Chromosome-parallel COJO was disabled for this run because the "
                "GCTA-usable BIM chromosomes are not unique positive numeric "
                "codes. PostGWAS will run one genome-wide GCTA command so no "
                "chromosome is silently omitted."
            )
        observation = (
            getattr(args, "variant_id_observations", None) or {}
        ).get(_GCTA_FORMATTER_TARGET)
        if observation is None:
            formatting = configuration.modules.formatting
            configure_reference_variant_identifiers(
                args,
                formatting,
                [BimIdentifierRequirement(
                    consumer=gcta_cojo_analysis_label(module),
                    formatter_target=_GCTA_FORMATTER_TARGET,
                    bim_file=reference_files["bim"],
                    column_roles=module.reference.bim_columns,
                    delimiter_pattern=module.reference.table_delimiter_pattern,
                )],
            )
            observation = args.variant_id_observations[_GCTA_FORMATTER_TARGET]

        resolved_path = configured_output_path(
            output,
            module.output_layout.resolved_config_file,
            error_type=GctaCojoError,
            dataset_id=dataset,
            mode=module.mode,
        )
        write_resolved_configuration(
            configuration,
            resolved_path,
            # Direct and pipeline paths share the formatter's canonical GCTA
            # .ma schema; direct mode also accepts GCTA's official header.
            modules=("formatting", "gcta_cojo"),
            resource_paths=(
                "executables.gcta",
                "genomes.%s" % module.genome_build,
                "populations.%s" % module.reference.population,
            ),
        )
        logger.record(
            "PARAM", "gcta_cojo",
            mode=module.mode,
            genome_build=module.genome_build,
            reference_population=module.reference.population,
            reference_prefix=str(reference_prefix),
            analysis=module.analysis.model_dump(mode="json"),
            inputs=module.inputs.model_dump(mode="json"),
            input_validation=module.input_validation.model_dump(mode="json"),
            reporting=module.reporting.model_dump(mode="json"),
            chromosome_execution=module.chromosome_execution.model_dump(
                mode="json",
            ),
            execution_strategy=(
                "chromosome_parallel" if parallel_slct else "single_command"
            ),
            requested_parallel_output_contract=parallel_slct_output_contract,
            effective_output_contract=effective_output_contract,
            threads=configuration.execution.threads,
            gcta_version=version,
            software=software,
        )
        logger.record(
            "OBSERVED", "reference_variant_identifiers", **observation,
        )
        logger.record(
            "OBSERVED", "cojo_input_validation",
            reference_samples=reference_samples,
            **input_metrics,
        )
        if pipeline_resources is not None:
            logger.record(
                "PASS",
                "gcta_cojo_pipeline_resource_preflight",
                reference_files={
                    name: str(path) for name, path in reference_files.items()
                },
                gcta_executable=executable,
                gcta_version=version,
                reference_identifier_observation=dict(
                    pipeline_resources.identifier_observation
                ),
                configured_snp_lists=dict(
                    pipeline_resources.snp_list_metrics
                ),
                deferred_checks_completed=list(
                    (
                        "formatter_created_gcta_ma",
                        "summary_bim_id_and_allele_compatibility",
                        "configured_snp_list_membership",
                    )
                ),
            )
        for warning in warnings:
            logger.warning(warning)

        layout = module.output_layout
        final_prefix = configured_output_path(
            output, layout.output_prefix, error_type=GctaCojoError,
            dataset_id=dataset, mode=module.mode,
        )
        final_result = configured_output_path(
            output, layout.primary_results[module.mode], error_type=GctaCojoError,
            dataset_id=dataset, mode=module.mode,
        )
        normalized = configured_output_path(
            output, layout.normalized_result, error_type=GctaCojoError,
            dataset_id=dataset, mode=module.mode,
        )
        completion = configured_output_path(
            output, layout.completion_manifest, error_type=GctaCojoError,
            dataset_id=dataset, mode=module.mode,
        )
        gcta_log = configured_output_path(
            output, layout.gcta_log_result, error_type=GctaCojoError,
            dataset_id=dataset, mode=module.mode,
        )
        ld_matrix = configured_output_path(
            output, layout.ld_matrix_result, error_type=GctaCojoError,
            dataset_id=dataset, mode=module.mode,
        )
        conditional_results = configured_output_path(
            output, layout.conditional_results, error_type=GctaCojoError,
            dataset_id=dataset, mode=module.mode,
        )
        condition_snps_result = configured_output_path(
            output, layout.condition_snps_result, error_type=GctaCojoError,
            dataset_id=dataset, mode=module.mode,
        )
        allele_mismatch_exclusion = configured_output_path(
            output,
            layout.allele_mismatch_exclusion,
            error_type=GctaCojoError,
            dataset_id=dataset,
            mode=module.mode,
        )
        output_paths = {
            "primary_result": final_result,
            "normalized_result": normalized,
            "gcta_log": gcta_log,
            "ld_matrix": ld_matrix,
            "conditional_results": conditional_results,
            "condition_snps": condition_snps_result,
            "allele_mismatch_exclusion": allele_mismatch_exclusion,
        }
        manifest_inputs = _input_manifest_paths(
            summary_file, reference_files, module,
        )
        digest = configuration_digest({
            "module": module.model_dump(mode="json"),
            "gcta_version": version,
            "software": software,
            "threads": configuration.execution.threads,
            "effective_output_contract": effective_output_contract,
        })
        resumed = False
        result_metrics = None
        if (
            configuration.run.resume
            and not configuration.run.overwrite
            and completion.is_file()
        ):
            expected_outputs = _expected_output_contract(
                module,
                output_paths,
                no_signals=_completion_records_no_signals(completion),
                allele_mismatches=input_metrics[
                    "reference_allele_mismatch_variants"
                ],
                include_conditional_results=include_conditional_results,
            )
            decision = resolve_completion_resume(
                completion,
                dataset_id=dataset,
                module="gcta_cojo",
                genome_build=module.genome_build,
                configuration_sha256=digest,
                inputs=manifest_inputs,
                outputs=expected_outputs,
                resume_policy=configuration.run.resume_policy,
                error_type=GctaCojoError,
            )
            if decision.action == "resume":
                result_metrics = dict(decision.manifest.get("metrics") or {})
                resumed = True
                logger.record(
                    "SKIP", "gcta_cojo_execution", reason="validated_resume",
                    manifest=str(completion),
                )
            else:
                apply_completion_restart(
                    decision,
                    output_root=output,
                    manifest=completion,
                    logger=logger,
                    operation="gcta_cojo_resume",
                    error_type=GctaCojoError,
                )
        if not resumed:
            if (
                not include_conditional_results
                and conditional_results.exists()
            ):
                if not configuration.run.overwrite:
                    raise GctaCojoError(
                        "A conditional-results file exists but is not part of "
                        "the resolved selection-only output contract: %s. Use "
                        "--overwrite to replace this stale output safely."
                        % conditional_results
                    )
                conditional_results.unlink()
            if (
                not input_metrics["reference_allele_mismatch_variants"]
                and allele_mismatch_exclusion.exists()
            ):
                if not configuration.run.overwrite:
                    raise GctaCojoError(
                        "Output already exists: %s. Use --resume or --overwrite."
                        % allele_mismatch_exclusion
                    )
                allele_mismatch_exclusion.unlink()
            expected_outputs = _expected_output_contract(
                module,
                output_paths,
                no_signals=False,
                allele_mismatches=input_metrics[
                    "reference_allele_mismatch_variants"
                ],
                include_conditional_results=include_conditional_results,
            )
            existing = [
                str(path) for path in expected_outputs.values() if path.exists()
            ]
            if existing and not configuration.run.overwrite:
                raise GctaCojoError(
                    "Output already exists: %s. Use --resume or --overwrite."
                    % ", ".join(existing)
                )
            staging = configured_output_path(
                output, layout.staging_directory, error_type=GctaCojoError,
                dataset_id=dataset, mode=module.mode,
            )
            if staging.exists():
                if not configuration.run.overwrite:
                    raise GctaCojoError(
                        "Incomplete staging output exists: %s. Review it, then use "
                        "--overwrite to replace it." % staging
                    )
                shutil.rmtree(staging)
            staging.mkdir(parents=True)
            staged_prefix = staging / final_prefix.name
            staged_result = _staged_path(final_prefix, final_result, staged_prefix)
            staged_gcta_log = _staged_path(
                final_prefix, gcta_log, staged_prefix,
            )
            staged_ld_matrix = _staged_path(
                final_prefix, ld_matrix, staged_prefix,
            )
            staged_conditional_results = _staged_path(
                final_prefix, conditional_results, staged_prefix,
            )
            staged_condition_snps = _staged_path(
                final_prefix, condition_snps_result, staged_prefix,
            )
            effective_exclude = None
            if temporary_exclusion is not None:
                staged_allele_mismatch_exclusion = _staged_path(
                    final_prefix, allele_mismatch_exclusion, staged_prefix,
                )
                shutil.copyfile(
                    temporary_exclusion, staged_allele_mismatch_exclusion,
                )
                effective_exclude = staged_allele_mismatch_exclusion
            staged_normalized = staging / normalized.relative_to(output)
            dry_run = bool(getattr(args, "dry_run", False))
            if parallel_slct:
                chromosome_count = len(
                    input_metrics["eligible_chromosome_variant_counts"]
                )
                with logger.step(
                    1,
                    1,
                    "Run GCTA-COJO slct across %d chromosomes"
                    % chromosome_count,
                    "run_cojo_chromosomes",
                    rows_in=input_metrics["summary_variants"],
                ) as step:
                    parallel_result = execute_parallel_slct(
                        configuration=configuration,
                        module=module,
                        logger=logger,
                        executable=executable,
                        summary_file=summary_file,
                        reference_prefix=reference_prefix,
                        staging=staging,
                        combined_prefix=staged_prefix,
                        combined_primary=staged_result,
                        combined_ld=staged_ld_matrix,
                        combined_conditional=staged_conditional_results,
                        combined_condition_snps=staged_condition_snps,
                        combined_log=staged_gcta_log,
                        normalized_destination=staged_normalized,
                        effective_exclude=effective_exclude,
                        chromosome_counts=input_metrics[
                            "eligible_chromosome_variant_counts"
                        ],
                        output_contract=effective_output_contract,
                        dry_run=dry_run,
                    )
                    if parallel_result.dry_run:
                        step.output(
                            "validated_commands",
                            output_prefix=str(final_prefix),
                            chromosomes=parallel_result.metrics["chromosomes"],
                            workers=parallel_result.metrics["parallel_workers"],
                        )
                        shutil.rmtree(staging)
                        return _return_dry_run_result(
                            ctx=ctx,
                            logger=logger,
                            emit_terminal_summary=emit_terminal_summary,
                            dataset=dataset,
                            module=module,
                            configuration=configuration,
                            reference_samples=reference_samples,
                            observation=observation,
                            software=software,
                            version=version,
                            input_metrics=input_metrics,
                            execution_metrics=parallel_result.metrics,
                            warnings=warnings,
                            final_prefix=final_prefix,
                            log_path=log_path,
                        )
                    result_metrics = parallel_result.metrics
                    no_signals = parallel_result.no_signals
                    _publish_completed_execution(
                        staging=staging,
                        staged_prefix=staged_prefix,
                        final_prefix=final_prefix,
                        staged_normalized=staged_normalized,
                        normalized=normalized,
                        output_paths=output_paths,
                        module=module,
                        configuration=configuration,
                        input_metrics=input_metrics,
                        no_signals=no_signals,
                        result_metrics=result_metrics,
                        include_conditional_results=include_conditional_results,
                        completion=completion,
                        dataset=dataset,
                        digest=digest,
                        manifest_inputs=manifest_inputs,
                        step=step,
                    )
            else:
                command = build_cojo_command(
                    executable,
                    summary_file,
                    reference_prefix,
                    staged_prefix,
                    configuration.execution.threads,
                    module,
                    exclude_snps=effective_exclude,
                )
                with logger.step(
                    1, 1, "Run GCTA-COJO %s" % module.mode, "run_cojo_command",
                    rows_in=input_metrics["summary_variants"],
                ) as step:
                    run_cojo_command(
                        command,
                        configuration,
                        logger,
                        dry_run=dry_run,
                    )
                    if dry_run:
                        step.output(
                            "validated_command", output_prefix=str(final_prefix),
                        )
                        shutil.rmtree(staging)
                        return _return_dry_run_result(
                            ctx=ctx,
                            logger=logger,
                            emit_terminal_summary=emit_terminal_summary,
                            dataset=dataset,
                            module=module,
                            configuration=configuration,
                            reference_samples=reference_samples,
                            observation=observation,
                            software=software,
                            version=version,
                            input_metrics=input_metrics,
                            execution_metrics={
                                "execution_strategy": "single_command",
                                "parallel_output_contract": "not_applicable",
                                "conditional_reconstruction_status": (
                                    "native_single_command"
                                    if module.mode in {"slct", "top_snps"}
                                    else "not_applicable"
                                ),
                                "command": command,
                            },
                            warnings=warnings,
                            final_prefix=final_prefix,
                            log_path=log_path,
                        )
                    require_nonempty_file(
                        staged_gcta_log,
                        "GCTA log",
                        error_type=GctaCojoError,
                    )
                    no_signals = (
                        module.mode in {"slct", "top_snps"}
                        and not staged_result.is_file()
                        and _gcta_reported_no_signals(staged_gcta_log, module)
                    )
                    if no_signals:
                        result_metrics, _ = write_empty_cojo_results(
                            staged_normalized, module,
                        )
                    else:
                        require_nonempty_file(
                            staged_result,
                            "GCTA primary result",
                            error_type=GctaCojoError,
                        )
                        if module.mode != "cond":
                            require_nonempty_file(
                                staged_ld_matrix,
                                "GCTA LD matrix",
                                error_type=GctaCojoError,
                            )
                        if module.mode in {"slct", "top_snps"}:
                            require_nonempty_file(
                                staged_conditional_results,
                                "GCTA conditional result",
                                error_type=GctaCojoError,
                            )
                        if module.mode == "cond":
                            require_nonempty_file(
                                staged_condition_snps,
                                "GCTA conditioning-SNP result",
                                error_type=GctaCojoError,
                            )
                            _validate_modeled_snp_survival(
                                module.inputs.condition_snps,
                                staged_condition_snps,
                                module,
                                "conditioning",
                            )
                        if module.mode == "joint":
                            _validate_modeled_snp_survival(
                                module.inputs.joint_snps,
                                staged_result,
                                module,
                                "joint-model",
                            )
                        result_metrics, _ = normalize_cojo_results(
                            staged_result, staged_normalized, module,
                        )
                    result_metrics.update({
                        "execution_strategy": "single_command",
                        "parallel_output_contract": "not_applicable",
                        "conditional_reconstruction_status": (
                            "native_single_command"
                            if module.mode in {"slct", "top_snps"}
                            else "not_applicable"
                        ),
                        "command": command,
                    })
                    _publish_completed_execution(
                        staging=staging,
                        staged_prefix=staged_prefix,
                        final_prefix=final_prefix,
                        staged_normalized=staged_normalized,
                        normalized=normalized,
                        output_paths=output_paths,
                        module=module,
                        configuration=configuration,
                        input_metrics=input_metrics,
                        no_signals=no_signals,
                        result_metrics=result_metrics,
                        include_conditional_results=include_conditional_results,
                        completion=completion,
                        dataset=dataset,
                        digest=digest,
                        manifest_inputs=manifest_inputs,
                        step=step,
                    )

        if result_metrics and result_metrics.get("no_signals"):
            no_signals_warning = (
                "GCTA-COJO completed successfully but selected no independent "
                "signals under the resolved model."
            )
            warnings.append(no_signals_warning)
            logger.warning(no_signals_warning)
        gcta_warnings = _read_gcta_warnings(gcta_log, module)
        for warning in gcta_warnings:
            if warning not in warnings:
                warnings.append(warning)
                logger.warning("GCTA: %s" % warning)
        metrics = {
            "mode": module.mode,
            "dry_run": False,
            "resumed": resumed,
            "gcta_version": version,
            "postgwas_version": software["postgwas"],
            "python_version": software["python"],
            "reference_samples": reference_samples,
            "variant_id_type": observation["variant_id_type"],
            **input_metrics,
            **(result_metrics or {}),
        }
        metadata = {
            "dataset_id": dataset,
            "mode": module.mode,
            "genome_build": module.genome_build,
            "reference_population": module.reference.population,
            "variant_id_type": observation["variant_id_type"],
            "effect_allele": "A1",
            "effect_allele_origin": (
                "VCF_ALT" if formatter_generated_input else "input_ma"
            ),
            "gcta_version": version,
            "postgwas_version": software["postgwas"],
            "python_version": software["python"],
        }
        artifacts = {
            "normalized_results": Artifact(
                "gcta_cojo_results", normalized, metadata,
            ),
            "gcta_log": Artifact("gcta_log", gcta_log, metadata),
            "completion_manifest": Artifact(
                "gcta_cojo_completion_manifest", completion, metadata,
            ),
        }
        if final_result.is_file():
            artifacts["raw_results"] = Artifact(
                "gcta_cojo_raw", final_result, metadata,
            )
        if module.mode != "cond" and ld_matrix.is_file():
            artifacts["ld_matrix"] = Artifact(
                "gcta_cojo_ld_matrix", ld_matrix, metadata,
            )
        if (
            module.mode in {"slct", "top_snps"}
            and include_conditional_results
            and conditional_results.is_file()
        ):
            artifacts["conditional_results"] = Artifact(
                "gcta_cojo_conditional_results", conditional_results, metadata,
            )
        if module.mode == "cond" and condition_snps_result.is_file():
            artifacts["condition_snps"] = Artifact(
                "gcta_cojo_condition_snps", condition_snps_result, metadata,
            )
        if (
            input_metrics["reference_allele_mismatch_variants"]
            and allele_mismatch_exclusion.is_file()
        ):
            artifacts["allele_mismatch_exclusion"] = Artifact(
                "gcta_cojo_allele_mismatch_exclusion",
                allele_mismatch_exclusion,
                metadata,
            )
        result = ModuleResult(
            "gcta_cojo",
            artifacts=artifacts,
            metrics=metrics,
            warnings=tuple(warnings),
        )
        if ctx is not None:
            ctx.publish(result)
        logger.record(
            "STATUS", "gcta_cojo", status="COMPLETED", mode=module.mode,
            resumed=resumed, warnings=len(warnings),
            significant_findings=metrics["significant_rows"],
        )
        if emit_terminal_summary:
            print(_render_success_summary(
                dataset,
                module,
                metrics,
                artifacts,
                warnings,
                log_path,
                configuration.logging.terminal_label_width,
            ))
        return result
    except BaseException as exc:
        if not logger.summary()["failed"]:
            logger.error(
                "GCTA-COJO failed: %s: %s" % (type(exc).__name__, exc)
            )
        if emit_terminal_summary:
            print(_render_failure_summary(
                dataset,
                module.mode,
                str(exc),
                log_path,
                configuration.logging.terminal_label_width,
            ))
        raise
    finally:
        if temporary_exclusion is not None:
            temporary_exclusion.unlink(missing_ok=True)
        logger.close()


__all__ = [
    "GctaCojoPipelineResources",
    "GctaCojoReferenceValidation",
    "gcta_cojo_analysis_label",
    "gcta_cojo_pipeline_title",
    "preflight_gcta_cojo_pipeline",
    "resolve_gcta_cojo_configuration",
    "run_gcta_cojo_direct",
    "validate_gcta_cojo_reference",
]
