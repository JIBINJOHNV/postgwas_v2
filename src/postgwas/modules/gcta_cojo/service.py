"""Validated direct and pipeline execution service for GCTA-COJO."""

from __future__ import annotations

from pathlib import Path
import re
import shutil
import sqlite3
import tempfile
from typing import get_args

from postgwas.config import (
    load_configuration,
    load_run_configuration_for_module,
    write_resolved_configuration,
)
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.config.models.modules.gcta_cojo import GctaCojoMode
from postgwas.core.completion import (
    configuration_digest,
    validate_completion_manifest,
    write_completion_manifest,
)
from postgwas.core.contracts import Artifact, ModuleResult
from postgwas.core.paths import (
    configured_output_path,
    require_nonempty_file,
    resolve_executable,
    validate_filename_component,
)
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
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
from postgwas.modules.gcta_cojo.results import normalize_cojo_results


_GCTA_FORMATTER_TARGET = "gcta_gene"
_GCTA_FORMATTER_RESULT_KEY = "summary_statistics_input_file"


def _resolved_configuration(args):
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
    return load_run_configuration_for_module(
        "gcta_cojo",
        getattr(args, "run_config", None),
        module_overrides=module_overrides,
        global_overrides=global_overrides,
    )


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


def _validate_reference_files(reference_prefix: Path, module) -> dict[str, Path]:
    paths = {
        suffix.lstrip("."): Path(str(reference_prefix) + suffix)
        for suffix in module.reference.required_extensions
    }
    missing = [
        str(path) for path in paths.values()
        if not path.is_file() or path.stat().st_size <= 0
    ]
    if missing:
        raise GctaCojoError(
            "PLINK LD reference files are missing or empty: %s"
            % ", ".join(missing)
        )
    return paths


def _count_reference_samples(fam_file: Path) -> int:
    count = 0
    with fam_file.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            # Six fields are a PLINK FAM format invariant.
            if len(raw.split()) != 6:
                raise GctaCojoError(
                    "PLINK FAM line %d has %d fields; expected 6: %s"
                    % (line_number, len(raw.split()), fam_file)
                )
            count += 1
    if count == 0:
        raise GctaCojoError("PLINK FAM contains no reference samples: %s" % fam_file)
    return count


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
        "CREATE TABLE bim (variant_id TEXT PRIMARY KEY, allele1 TEXT, allele2 TEXT)"
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
            alleles = (
                fields[roles["allele1"]].strip().upper(),
                fields[roles["allele2"]].strip().upper(),
            )
            if not variant or not all(alleles) or alleles[0] == alleles[1]:
                raise GctaCojoError(
                    "PLINK BIM contains an invalid variant or allele at line %d."
                    % line_number
                )
            batch.append((variant, *alleles))
            count += 1
            if len(batch) >= module.input_validation.sqlite_batch_size:
                try:
                    connection.executemany("INSERT INTO bim VALUES (?, ?, ?)", batch)
                except sqlite3.IntegrityError as exc:
                    raise GctaCojoError(
                        "PLINK BIM contains duplicate variant identifiers near line %d."
                        % line_number
                    ) from exc
                batch.clear()
    if batch:
        try:
            connection.executemany("INSERT INTO bim VALUES (?, ?, ?)", batch)
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
    pattern = re.compile(module.input_validation.table_delimiter_pattern)
    batch, count = [], 0
    with summary_file.open("r", encoding="utf-8", errors="replace") as handle:
        header = handle.readline()
        observed_header = pattern.split(header.strip()) if header else []
        if observed_header != expected_header:
            raise GctaCojoError(
                "GCTA .ma header does not match the resolved schema. "
                "Expected %s; observed %s."
                % (" ".join(expected_header), " ".join(observed_header))
            )
        identifier_index = expected_header.index(
            schema.columns[configuration.modules.formatting.canonical_columns.resolved_variant_id]
        )
        a1_index = expected_header.index(
            schema.columns[configuration.modules.formatting.canonical_columns.alternate_allele]
        )
        a2_index = expected_header.index(
            schema.columns[configuration.modules.formatting.canonical_columns.reference_allele]
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
            variant = fields[identifier_index].strip()
            alleles = fields[a1_index].upper(), fields[a2_index].upper()
            if not variant or not all(alleles) or alleles[0] == alleles[1]:
                raise GctaCojoError(
                    "GCTA .ma contains an invalid variant or allele at line %d."
                    % line_number
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


def _validate_snp_list_membership(connection, table: str, path: Path, *, strict: bool):
    missing_summary = connection.execute(
        "SELECT COUNT(*) FROM %s l LEFT JOIN summary s USING (variant_id) "
        "WHERE s.variant_id IS NULL" % table
    ).fetchone()[0]
    missing_reference = connection.execute(
        "SELECT COUNT(*) FROM %s l LEFT JOIN bim b USING (variant_id) "
        "WHERE b.variant_id IS NULL" % table
    ).fetchone()[0]
    if strict and (missing_summary or missing_reference):
        raise GctaCojoError(
            "SNP list %s is incompatible with the GCTA .ma and LD reference: "
            "%d IDs are absent from the .ma file and %d are absent from the BIM."
            % (path, missing_summary, missing_reference)
        )
    return {
        "variants": connection.execute(
            "SELECT COUNT(*) FROM %s" % table
        ).fetchone()[0],
        "missing_summary": missing_summary,
        "missing_reference": missing_reference,
    }


def _validate_inputs(
    summary_file: Path,
    reference_files: dict[str, Path],
    output: Path,
    configuration,
) -> tuple[dict, list[str]]:
    module = configuration.modules.gcta_cojo
    connection, database_path = _create_validation_database(output, configuration)
    warnings = []
    try:
        bim_variants = _insert_bim(connection, reference_files["bim"], module)
        summary_variants = _insert_summary(connection, summary_file, configuration)
        overlap = connection.execute(
            "SELECT COUNT(*) FROM summary INNER JOIN bim USING (variant_id)"
        ).fetchone()[0]
        incompatible_alleles = connection.execute(
            "SELECT COUNT(*) FROM summary s INNER JOIN bim b USING (variant_id) "
            "WHERE NOT ((UPPER(s.allele1)=UPPER(b.allele1) AND "
            "UPPER(s.allele2)=UPPER(b.allele2)) OR "
            "(UPPER(s.allele1)=UPPER(b.allele2) AND "
            "UPPER(s.allele2)=UPPER(b.allele1)))"
        ).fetchone()[0]
        if incompatible_alleles:
            raise GctaCojoError(
                "%d variants share an ID between the GCTA .ma and BIM but "
                "have incompatible alleles." % incompatible_alleles
            )
        overlap_fraction = overlap / summary_variants
        if overlap_fraction < module.input_validation.minimum_reference_overlap_fraction:
            raise GctaCojoError(
                "Only %.2f%% of GCTA .ma variants overlap the PLINK BIM; "
                "the configured minimum is %.2f%%. Check genome build, ancestry, "
                "and variant identifiers."
                % (
                    100 * overlap_fraction,
                    100 * module.input_validation.minimum_reference_overlap_fraction,
                )
            )
        if overlap < summary_variants:
            warnings.append(
                "%d of %d GCTA .ma variants are absent from the LD-reference "
                "BIM and will not enter COJO."
                % (summary_variants - overlap, summary_variants)
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
                    "Exclude list contains %d IDs absent from the .ma file and %d "
                    "absent from the BIM; those exclusions have no effect."
                    % (
                        observed["missing_summary"],
                        observed["missing_reference"],
                    )
                )
        return {
            "summary_variants": summary_variants,
            "reference_variants": bim_variants,
            "reference_overlap_variants": overlap,
            "reference_overlap_fraction": overlap_fraction,
            "snp_lists": list_metrics,
        }, warnings
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
        ("analysis", "COJO mode", module.mode),
        ("genetic", "LD population", module.reference.population),
        ("count", "Reference samples", "%s" % f"{metrics['reference_samples']:,}"),
        (
            "genetic", "BIM identifier type",
            metrics["variant_id_type"],
        ),
        (
            "count", "GWAS/reference overlap",
            "%s/%s variants (%.2f%%)"
            % (
                f"{metrics['reference_overlap_variants']:,}",
                f"{metrics['summary_variants']:,}",
                100 * metrics["reference_overlap_fraction"],
            ),
        ),
    ]
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
        p_column = metrics["adjusted_p_value_column"]
        for index, finding in enumerate(metrics["top_findings"], 1):
            lines.append(screen_field(
                "analysis",
                "Top finding %d" % index,
                "%s (chr%s:%s; %s=%s)"
                % (
                    finding["variant_id"], finding["chromosome"], finding["position"], p_column,
                    _format_p(
                        finding["adjusted_p_value"],
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


def run_gcta_cojo_direct(args, ctx=None) -> ModuleResult:
    """Validate an existing .ma artifact and run COJO in direct or pipeline mode."""
    try:
        configuration = _resolved_configuration(args)
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
        print(_render_failure_summary(
            dataset, mode, str(exc), log_path,
            fallback.logging.terminal_label_width,
        ))
        raise

    module = configuration.modules.gcta_cojo
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
    try:
        _validate_build_and_population(module, configuration)
        if module.reference.prefix is None:
            raise GctaCojoError(
                "Provide --cojo-reference-prefix or modules.gcta_cojo.reference.prefix."
            )
        reference_prefix = Path(module.reference.prefix).expanduser().resolve()
        reference_files = _validate_reference_files(reference_prefix, module)
        reference_samples = _count_reference_samples(reference_files["fam"])
        warnings = []
        if reference_samples <= module.input_validation.reference_sample_size_warning_threshold:
            warnings.append(
                "The LD reference has %d samples; GCTA recommends more than %d for "
                "COJO. Treat selected and conditional signals cautiously."
                % (
                    reference_samples,
                    module.input_validation.reference_sample_size_warning_threshold,
                )
            )

        executable = resolve_executable(
            configuration.resources.executables.gcta,
            "GCTA executable",
            error_type=GctaCojoError,
        )
        version = require_supported_gcta(
            executable, module, logger, configuration.execution.timeout_seconds,
        )

        input_file = getattr(args, "gcta_cojo_input_file", None) or module.input_file
        if input_file is None and ctx is not None:
            formatter = ctx.get("formatter", {})
            input_file = formatter.get(_GCTA_FORMATTER_TARGET, {}).get(
                _GCTA_FORMATTER_RESULT_KEY
            )
        if input_file is None:
            raise GctaCojoError(
                "Direct mode requires an existing GCTA .ma file; provide "
                "--cojo-file PATH. To start from a GWAS-VCF, use "
                "postgwas pipeline --modules gcta_cojo --vcf PATH."
            )
        summary_file = require_nonempty_file(
            input_file, "GCTA .ma input", error_type=GctaCojoError,
        )

        input_metrics, input_warnings = _validate_inputs(
            summary_file, reference_files, output, configuration,
        )
        warnings.extend(input_warnings)
        observation = (
            getattr(args, "variant_id_observations", None) or {}
        ).get(_GCTA_FORMATTER_TARGET)
        if observation is None:
            formatting = configuration.modules.formatting
            configure_reference_variant_identifiers(
                args,
                formatting,
                [BimIdentifierRequirement(
                    consumer="GCTA COJO",
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
            # The formatter is not executed in direct mode, but its shared GCTA
            # .ma schema is the canonical input contract validated above.
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
            threads=configuration.execution.threads,
        )
        logger.record(
            "OBSERVED", "reference_variant_identifiers", **observation,
        )
        logger.record(
            "OBSERVED", "cojo_input_validation",
            reference_samples=reference_samples,
            **input_metrics,
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
        expected_outputs = {
            "primary_result": final_result,
            "normalized_result": normalized,
            "gcta_log": gcta_log,
        }
        if module.mode != "cond":
            expected_outputs["ld_matrix"] = ld_matrix
        manifest_inputs = _input_manifest_paths(
            summary_file, reference_files, module,
        )
        digest = configuration_digest({
            "module": module.model_dump(mode="json"),
            "gcta_version": version,
            "threads": configuration.execution.threads,
        })
        resumed = False
        result_metrics = None
        if (
            configuration.run.resume
            and not configuration.run.overwrite
            and completion.is_file()
        ):
            manifest = validate_completion_manifest(
                completion,
                dataset_id=dataset,
                module="gcta_cojo",
                genome_build=module.genome_build,
                configuration_sha256=digest,
                inputs=manifest_inputs,
                outputs=expected_outputs,
                error_type=GctaCojoError,
            )
            result_metrics = dict(manifest.get("metrics") or {})
            resumed = True
            logger.record(
                "SKIP", "gcta_cojo_execution", reason="validated_resume",
                manifest=str(completion),
            )
        else:
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
            command = build_cojo_command(
                executable,
                summary_file,
                reference_prefix,
                staged_prefix,
                configuration.execution.threads,
                module,
            )
            with logger.step(
                1, 1, "Run GCTA-COJO %s" % module.mode, "run_cojo_command",
                rows_in=input_metrics["summary_variants"],
            ) as step:
                run_cojo_command(
                    command,
                    staged_result,
                    configuration,
                    logger,
                    dry_run=bool(getattr(args, "dry_run", False)),
                )
                if getattr(args, "dry_run", False):
                    step.output("validated_command", output_prefix=str(final_prefix))
                    shutil.rmtree(staging)
                    result = ModuleResult(
                        "gcta_cojo",
                        metrics={
                            "mode": module.mode,
                            "dry_run": True,
                            "gcta_version": version,
                            "reference_samples": reference_samples,
                            "variant_id_type": observation["variant_id_type"],
                            **input_metrics,
                        },
                        warnings=tuple(warnings),
                    )
                    if ctx is not None:
                        ctx.publish(result)
                    logger.record(
                        "STATUS", "gcta_cojo", status="VALIDATED", mode=module.mode,
                    )
                    print("\n".join([
                        "",
                        screen_line("analysis", "GCTA-COJO command validated", indent=2),
                        screen_field(
                            "info", "Dataset", dataset, indent=6,
                            label_width=configuration.logging.terminal_label_width,
                        ),
                        screen_field(
                            "info", "Full PostGWAS log", log_path, indent=6,
                            label_width=configuration.logging.terminal_label_width,
                        ),
                        "",
                    ]))
                    return result
                staged_normalized = staging / normalized.relative_to(output)
                result_metrics, _ = normalize_cojo_results(
                    staged_result, staged_normalized, module,
                )
                published = _publish_staged_outputs(
                    staged_prefix, final_prefix, configuration.run.overwrite,
                )
                normalized.parent.mkdir(parents=True, exist_ok=True)
                if normalized.exists():
                    if not configuration.run.overwrite:
                        raise GctaCojoError("Output already exists: %s" % normalized)
                    normalized.unlink()
                staged_normalized.replace(normalized)
                step.set_rows(result_metrics["result_rows"], removed=0)
                step.output("cojo_outputs", files=[str(path) for path in published])
                step.output("normalized_result", path=str(normalized))
            shutil.rmtree(staging)
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
            "effect_allele": "A1=ALT",
            "gcta_version": version,
        }
        artifacts = {
            "raw_results": Artifact("gcta_cojo_raw", final_result, metadata),
            "normalized_results": Artifact(
                "gcta_cojo_results", normalized, metadata,
            ),
            "gcta_log": Artifact("gcta_log", gcta_log, metadata),
            "completion_manifest": Artifact(
                "gcta_cojo_completion_manifest", completion, metadata,
            ),
        }
        if module.mode != "cond":
            artifacts["ld_matrix"] = Artifact(
                "gcta_cojo_ld_matrix", ld_matrix, metadata,
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
        print(_render_failure_summary(
            dataset,
            module.mode,
            str(exc),
            log_path,
            configuration.logging.terminal_label_width,
        ))
        raise
    finally:
        logger.close()


__all__ = ["run_gcta_cojo_direct"]
