"""Validated, atomic execution of LDSC cell-type-specific analyses."""

from __future__ import annotations

import bz2
from dataclasses import dataclass
import gzip
from pathlib import Path
import re
from typing import Any, Mapping

from postgwas.core.completion import (
    apply_completion_restart,
    configuration_digest,
    resolve_completion_resume,
    write_completion_manifest,
)
from postgwas.core.io.reports import write_yaml_report
from postgwas.core.paths import (
    remove_empty_directories,
    remove_owned_directory,
    require_nonempty_file,
    resolve_executable,
)
from postgwas.core.processes import run_checked_command
from postgwas.modules.ldsc.ldsc_runner import (
    build_h2_cts_command,
    build_munge_sumstats_command,
)
from postgwas.modules.single_cell.errors import SingleCellError
from postgwas.modules.single_cell.methods.ldsc_celltype.analysis import (
    normalize_ldsc_celltype_results,
)


@dataclass(frozen=True)
class LdscCelltypeManifestEntry:
    """One tested annotation and its required control annotation prefixes."""

    label: str
    prefixes: tuple[str, ...]


@dataclass(frozen=True)
class LdscCelltypePreflight:
    """Resolved inputs and complete reference inventory for one LDSC run."""

    source: str
    sumstats_file: Path | None
    ldcts_file: Path
    merge_alleles_file: Path | None
    baseline_ld_prefixes: tuple[str, ...]
    weights_ld_prefix: str
    entries: tuple[LdscCelltypeManifestEntry, ...]
    executable: str
    munge_executable: str | None
    software_version: str
    reference_inventory: tuple[Mapping[str, Any], ...]


@dataclass(frozen=True)
class LdscCelltypeExecution:
    """Published LDSC outputs and execution metrics."""

    outputs: Mapping[str, Path]
    metrics: Mapping[str, Any]
    resumed: bool


def _resolve_prefix(value: str, *, relative_to: Path) -> str:
    candidate = Path(value).expanduser()
    if not candidate.is_absolute():
        candidate = relative_to / candidate
    return str(candidate.resolve())


def _read_ldcts(path: Path, method) -> tuple[LdscCelltypeManifestEntry, ...]:
    """Parse the exact two-field syntax consumed by LDSC ``--ref-ld-chr-cts``."""
    configured = method.reference_format
    entries: list[LdscCelltypeManifestEntry] = []
    labels: set[str] = set()
    tested_prefixes: set[str] = set()
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise SingleCellError(
            "Cannot read LDSC .ldcts file %s: %s" % (path, exc)
        ) from exc
    for line_number, line in enumerate(lines, start=1):
        fields = line.split()
        if len(fields) != 2:
            raise SingleCellError(
                "LDSC .ldcts line %d must contain exactly a cell-type name and "
                "a comma-delimited prefix list" % line_number
            )
        label, raw_prefixes = fields
        if label in labels:
            raise SingleCellError(
                "Duplicate LDSC .ldcts cell-type name: %s" % label
            )
        labels.add(label)
        values = raw_prefixes.split(configured.prefix_separator)
        if any(not value.strip() for value in values):
            raise SingleCellError(
                "LDSC .ldcts line %d contains an empty LD-score prefix" % line_number
            )
        if len(values) < configured.minimum_prefixes_per_cell_type:
            raise SingleCellError(
                "LDSC .ldcts line %d requires at least %d prefixes: the tested "
                "annotation followed by its all-genes control"
                % (line_number, configured.minimum_prefixes_per_cell_type)
            )
        prefixes = tuple(
            _resolve_prefix(value.strip(), relative_to=path.parent)
            for value in values
        )
        if len(prefixes) != len(set(prefixes)):
            raise SingleCellError(
                "LDSC .ldcts line %d contains duplicate LD-score prefixes"
                % line_number
            )
        if prefixes[0] in tested_prefixes:
            raise SingleCellError(
                "LDSC .ldcts tested annotation prefix is repeated: %s"
                % prefixes[0]
            )
        tested_prefixes.add(prefixes[0])
        entries.append(LdscCelltypeManifestEntry(label, prefixes))
    if not entries:
        raise SingleCellError("LDSC .ldcts file must contain at least one cell type")
    return tuple(entries)


def _chromosome_prefix(prefix: str, chromosome: int, placeholder: str) -> str:
    occurrences = prefix.count(placeholder)
    if occurrences > 1:
        raise SingleCellError(
            "LDSC reference prefix contains the chromosome placeholder more "
            "than once: %s" % prefix
        )
    if occurrences:
        return prefix.replace(placeholder, str(chromosome))
    return "%s%d" % (prefix, chromosome)


def _reference_file_record(path: Path, **metadata: Any) -> dict[str, Any]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise SingleCellError(
            "LDSC reference file does not exist or is empty: %s" % path
        )
    stat = path.stat()
    return {
        **metadata,
        "path": str(path),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def _validate_prefix(
    prefix: str,
    *,
    role: str,
    method,
    require_m_file: bool,
) -> list[dict[str, Any]]:
    configured = method.reference_format
    records = []
    for chromosome in range(1, configured.chromosomes + 1):
        chromosome_prefix = _chromosome_prefix(
            prefix, chromosome, configured.chromosome_placeholder,
        )
        ldscore = Path(chromosome_prefix + configured.ldscore_suffix)
        records.append(_reference_file_record(
            ldscore,
            role=role,
            prefix=prefix,
            chromosome=chromosome,
            resource="ld_score",
        ))
        if require_m_file:
            m_file = Path(chromosome_prefix + configured.m_suffix)
            records.append(_reference_file_record(
                m_file,
                role=role,
                prefix=prefix,
                chromosome=chromosome,
                resource="regression_snp_count",
            ))
    return records


def _validate_references(
    baseline_ld_prefixes: tuple[str, ...],
    weights_ld_prefix: str,
    entries: tuple[LdscCelltypeManifestEntry, ...],
    method,
) -> tuple[Mapping[str, Any], ...]:
    records: list[dict[str, Any]] = []
    for prefix in baseline_ld_prefixes:
        records.extend(_validate_prefix(
            prefix, role="baseline", method=method, require_m_file=True,
        ))
    records.extend(_validate_prefix(
        weights_ld_prefix,
        role="regression_weights",
        method=method,
        require_m_file=False,
    ))
    validated_cts_prefixes: set[str] = set()
    for entry in entries:
        for index, prefix in enumerate(entry.prefixes):
            if prefix in validated_cts_prefixes:
                continue
            validated_cts_prefixes.add(prefix)
            records.extend(_validate_prefix(
                prefix,
                role=("cell_type_annotation" if index == 0 else "all_genes_control"),
                method=method,
                require_m_file=True,
            ))
    return tuple(records)


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    if path.suffix == ".bz2":
        return bz2.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def validate_munged_sumstats(path: str | Path, method) -> dict[str, Any]:
    """Validate columns required by the LDSC regression before execution."""
    sumstats = require_nonempty_file(
        path, "munged LDSC summary statistics", error_type=SingleCellError,
    )
    try:
        with _open_text(sumstats) as handle:
            header_line = handle.readline()
            first_data_line = next((line for line in handle if line.strip()), "")
    except (OSError, UnicodeError) as exc:
        raise SingleCellError(
            "Cannot read munged LDSC summary statistics %s: %s" % (sumstats, exc)
        ) from exc
    columns = header_line.split()
    required = method.result_format.sumstats_required_columns
    missing = [column for column in required if column not in columns]
    if missing:
        raise SingleCellError(
            "Munged LDSC summary statistics are missing required columns: %s"
            % ", ".join(missing)
        )
    if len(columns) != len(set(columns)):
        raise SingleCellError("Munged LDSC summary-statistic columns must be unique")
    if not first_data_line:
        raise SingleCellError("Munged LDSC summary statistics contain no variants")
    return {"columns": columns, "required_columns": list(required)}


def _probe_version(executable: str, method) -> str:
    output = run_checked_command(
        [executable, *method.version_probe_arguments],
        "LDSC version probe",
        error_type=SingleCellError,
    )
    matches = re.findall(method.version_pattern, output)
    if len(matches) != 1:
        raise SingleCellError(
            "LDSC version probe did not return exactly one parseable version"
        )
    return matches[0]


def preflight_ldsc_celltype(
    method,
    ldsc_executable_value: str | Path,
    munge_executable_value: str | Path,
    *,
    pipeline_pending_sumstats: bool = False,
) -> LdscCelltypePreflight:
    """Resolve and validate all inputs needed by the official LDSC workflow."""
    ldcts_file = require_nonempty_file(
        method.input.ldcts_file, "LDSC .ldcts manifest", error_type=SingleCellError,
    )
    if not method.input.baseline_ld_prefixes:
        raise SingleCellError("At least one baseline LD-score prefix is required")
    if method.input.weights_ld_prefix is None:
        raise SingleCellError("A regression-weights LD-score prefix is required")
    cwd = Path.cwd()
    baseline_ld_prefixes = tuple(
        _resolve_prefix(prefix, relative_to=cwd)
        for prefix in method.input.baseline_ld_prefixes
    )
    weights_ld_prefix = _resolve_prefix(
        method.input.weights_ld_prefix, relative_to=cwd,
    )
    entries = _read_ldcts(ldcts_file, method)
    reference_inventory = _validate_references(
        baseline_ld_prefixes, weights_ld_prefix, entries, method,
    )
    executable = resolve_executable(
        ldsc_executable_value, "LDSC executable", error_type=SingleCellError,
    )
    software_version = _probe_version(executable, method)
    munge_executable = None
    merge_alleles_file = None
    if method.input.source == "formatter":
        munge_executable = resolve_executable(
            munge_executable_value,
            "LDSC munge_sumstats executable",
            error_type=SingleCellError,
        )
        merge_alleles_file = require_nonempty_file(
            method.input.merge_alleles_file,
            "LDSC HapMap3 merge-alleles file",
            error_type=SingleCellError,
        )
    sumstats_file = None
    if not pipeline_pending_sumstats:
        sumstats_file = require_nonempty_file(
            method.input.sumstats_file,
            (
                "formatter LDSC summary statistics"
                if method.input.source == "formatter"
                else "munged LDSC summary statistics"
            ),
            error_type=SingleCellError,
        )
        if method.input.source == "munged":
            validate_munged_sumstats(sumstats_file, method)
    return LdscCelltypePreflight(
        source=method.input.source,
        sumstats_file=sumstats_file,
        ldcts_file=ldcts_file,
        merge_alleles_file=merge_alleles_file,
        baseline_ld_prefixes=baseline_ld_prefixes,
        weights_ld_prefix=weights_ld_prefix,
        entries=entries,
        executable=executable,
        munge_executable=munge_executable,
        software_version=software_version,
        reference_inventory=reference_inventory,
    )


def _write_normalized_ldcts(
    path: Path,
    entries: tuple[LdscCelltypeManifestEntry, ...],
    separator: str,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        for entry in entries:
            handle.write("%s\t%s\n" % (entry.label, separator.join(entry.prefixes)))
    return path


def _expected_outputs(
    paths: Mapping[str, Path], method, source: str,
) -> dict[str, Path]:
    native_prefix = paths["ldsc_celltype_native_output_prefix"]
    outputs = {
        "normalized_results": paths["ldsc_celltype_results_file"],
        "native_results": Path(
            str(native_prefix) + method.result_format.native_results_suffix
        ),
        "native_log": Path(
            str(native_prefix) + method.result_format.native_log_suffix
        ),
        "normalized_ldcts": paths["ldsc_celltype_normalized_ldcts"],
        "qc_report": paths["ldsc_celltype_qc_report"],
    }
    if source == "formatter":
        munge_prefix = paths["ldsc_celltype_munge_output_prefix"]
        outputs.update({
            "munged_sumstats": Path(
                str(munge_prefix) + method.result_format.munged_sumstats_suffix
            ),
            "munge_log": Path(
                str(munge_prefix) + method.result_format.munge_log_suffix
            ),
        })
    if len(outputs) != len(set(outputs.values())):
        raise SingleCellError("Configured LDSC cell-type output paths collide")
    return outputs


def _staged_path(path: Path, *, output: Path, staging: Path) -> Path:
    try:
        relative = path.relative_to(output)
    except ValueError as exc:
        raise SingleCellError(
            "Configured LDSC cell-type output must be inside the output root: %s"
            % path
        ) from exc
    return staging / relative


def run_ldsc_celltype(
    *,
    preflight: LdscCelltypePreflight,
    configuration,
    paths: Mapping[str, Path],
    dataset: str,
    logger,
) -> LdscCelltypeExecution:
    """Run LDSC ``--h2-cts`` atomically or validate a completed resume."""
    single_cell = configuration.modules.single_cell
    method = single_cell.ldsc_celltype
    output = Path(configuration.run.output_directory).expanduser().resolve()
    engine = paths["ldsc_celltype_engine_directory"]
    staging = paths["ldsc_celltype_staging_directory"]
    manifest_path = paths["ldsc_celltype_completion_manifest"]
    outputs = _expected_outputs(paths, method, preflight.source)
    for name, path in outputs.items():
        if name != "normalized_results" and engine not in path.parents:
            raise SingleCellError(
                "Configured LDSC cell-type %s must be inside its engine directory"
                % name
            )
    if preflight.sumstats_file is None:
        raise SingleCellError("LDSC cell-type execution requires summary statistics")
    completion_inputs = {
        "summary_statistics": preflight.sumstats_file,
        "ldcts_manifest": preflight.ldcts_file,
    }
    if preflight.merge_alleles_file is not None:
        completion_inputs["merge_alleles"] = preflight.merge_alleles_file
    digest = configuration_digest({
        "ldsc_celltype": method.model_dump(mode="json"),
        "ldsc_munging": (
            configuration.modules.ldsc.model_dump(mode="json")
            if preflight.source == "formatter" else None
        ),
        "multiple_testing": single_cell.multiple_testing.model_dump(mode="json"),
        "result_schema": single_cell.result_schema.model_dump(mode="json"),
        "output_layout": {
            key: value
            for key, value in single_cell.output_layout.model_dump().items()
            if key.startswith("ldsc_celltype_")
        },
        "ldsc_executable": preflight.executable,
        "munge_executable": preflight.munge_executable,
        "software_version": preflight.software_version,
        "reference_inventory": list(preflight.reference_inventory),
    })
    if (
        configuration.run.resume
        and not configuration.run.overwrite
        and manifest_path.is_file()
    ):
        decision = resolve_completion_resume(
            manifest_path,
            dataset_id=dataset,
            module="single_cell.ldsc_celltype",
            genome_build=method.genome_build.value,
            configuration_sha256=digest,
            inputs=completion_inputs,
            outputs=outputs,
            resume_policy=configuration.run.resume_policy,
            error_type=SingleCellError,
        )
        if decision.action == "resume":
            return LdscCelltypeExecution(
                outputs=outputs,
                metrics=dict(decision.manifest.get("metrics", {})),
                resumed=True,
            )
        apply_completion_restart(
            decision,
            output_root=output,
            manifest=manifest_path,
            logger=logger,
            operation="ldsc_celltype_resume",
            error_type=SingleCellError,
        )
        if engine.is_dir() and not any(
            path.is_file() or path.is_symlink()
            for path in engine.rglob("*")
        ):
            remove_owned_directory(
                engine,
                output,
                "empty stale LDSC cell-type engine directory",
                error_type=SingleCellError,
            )
    if (
        engine.exists()
        or outputs["normalized_results"].exists()
        or manifest_path.exists()
    ) and not configuration.run.overwrite:
        raise SingleCellError(
            "Existing or incomplete LDSC cell-type output was found; use "
            "--resume for a matching completed run or --overwrite to replace it"
        )
    if staging.exists():
        if not configuration.run.overwrite:
            raise SingleCellError(
                "An isolated incomplete LDSC cell-type run exists at %s; "
                "review it or use --overwrite" % staging
            )
        remove_owned_directory(
            staging,
            output,
            "incomplete LDSC cell-type staging directory",
            error_type=SingleCellError,
        )
    if configuration.run.overwrite:
        remove_owned_directory(
            engine,
            output,
            "LDSC cell-type engine directory",
            error_type=SingleCellError,
        )
        outputs["normalized_results"].unlink(missing_ok=True)
        manifest_path.unlink(missing_ok=True)

    staging.mkdir(parents=True, exist_ok=False)
    staged_outputs = {
        name: _staged_path(path, output=output, staging=staging)
        for name, path in outputs.items()
    }
    staged_prefix = _staged_path(
        paths["ldsc_celltype_native_output_prefix"],
        output=output,
        staging=staging,
    )
    staged_prefix.parent.mkdir(parents=True, exist_ok=True)
    _write_normalized_ldcts(
        staged_outputs["normalized_ldcts"],
        preflight.entries,
        method.reference_format.prefix_separator,
    )
    analysis_sumstats = preflight.sumstats_file
    if preflight.source == "formatter":
        if (
            preflight.munge_executable is None
            or preflight.merge_alleles_file is None
        ):
            raise SingleCellError(
                "Formatter-source LDSC execution requires munging resources"
            )
        staged_munge_prefix = _staged_path(
            paths["ldsc_celltype_munge_output_prefix"],
            output=output,
            staging=staging,
        )
        staged_munge_prefix.parent.mkdir(parents=True, exist_ok=True)
        ldsc_config = configuration.modules.ldsc
        run_checked_command(
            build_munge_sumstats_command(
                preflight.munge_executable,
                sumstats=preflight.sumstats_file,
                output_prefix=staged_munge_prefix,
                merge_alleles=preflight.merge_alleles_file,
                minimum_info=ldsc_config.minimum_info,
                minimum_maf=ldsc_config.minimum_maf,
                minimum_n=ldsc_config.minimum_n,
                chunksize=ldsc_config.chunksize,
                keep_maf=ldsc_config.keep_maf,
            ),
            "Munge formatter output for LDSC cell-type analysis",
            logger=logger,
            error_type=SingleCellError,
            expected_outputs=[
                staged_outputs["munged_sumstats"], staged_outputs["munge_log"],
            ],
        )
        analysis_sumstats = staged_outputs["munged_sumstats"]
        validate_munged_sumstats(analysis_sumstats, method)
    run_checked_command(
        build_h2_cts_command(
            preflight.executable,
            sumstats=analysis_sumstats,
            baseline_ld_prefixes=preflight.baseline_ld_prefixes,
            weights_ld_prefix=preflight.weights_ld_prefix,
            ldcts_file=staged_outputs["normalized_ldcts"],
            output_prefix=staged_prefix,
            prefix_separator=method.reference_format.prefix_separator,
        ),
        "LDSC cell-type-specific heritability regression",
        logger=logger,
        error_type=SingleCellError,
        expected_outputs=[
            staged_outputs["native_results"], staged_outputs["native_log"],
        ],
    )
    metrics = normalize_ldsc_celltype_results(
        staged_outputs["native_results"],
        staged_outputs["normalized_results"],
        expected_cell_types=tuple(entry.label for entry in preflight.entries),
        dataset_id=dataset,
        single_cell_config=single_cell,
    )
    metrics.update({
        "software_version": preflight.software_version,
        "genome_build": method.genome_build.value,
        "population": method.population.value,
        "sumstats_source": preflight.source,
        "baseline_prefixes": len(preflight.baseline_ld_prefixes),
        "reference_files_validated": len(preflight.reference_inventory),
        "chromosomes_validated": method.reference_format.chromosomes,
    })
    write_yaml_report(
        {
            "schema_version": 1,
            "status": "COMPLETED",
            "dataset_id": dataset,
            "method": "ldsc_celltype",
            "workflow": method.workflow,
            "software_version": preflight.software_version,
            "inputs": {
                "summary_statistics": str(preflight.sumstats_file),
                "summary_statistics_source": preflight.source,
                "ldcts_manifest": str(preflight.ldcts_file),
                "merge_alleles": (
                    None if preflight.merge_alleles_file is None
                    else str(preflight.merge_alleles_file)
                ),
            },
            "scientific_settings": {
                "genome_build": method.genome_build.value,
                "population": method.population.value,
                "tested_coefficient": "first_ldcts_prefix",
                "control_annotations": "remaining_ldcts_prefixes",
                "coefficient_interpretation": "additional_per_snp_heritability",
                "p_value_alternative": "coefficient_greater_than_zero",
                "multiple_testing_methods": list(
                    single_cell.multiple_testing.methods
                ),
            },
            "cell_types": [entry.label for entry in preflight.entries],
            "reference_inventory": list(preflight.reference_inventory),
            "metrics": metrics,
        },
        staged_outputs["qc_report"],
    )
    engine.parent.mkdir(parents=True, exist_ok=True)
    outputs["normalized_results"].parent.mkdir(parents=True, exist_ok=True)
    staged_engine = _staged_path(engine, output=output, staging=staging)
    staged_engine.replace(engine)
    staged_outputs["normalized_results"].replace(outputs["normalized_results"])
    remove_empty_directories(
        staged_outputs["normalized_results"].parent,
        staged_engine.parent,
        staging,
    )
    metrics["output"] = str(outputs["normalized_results"])
    write_completion_manifest(
        manifest_path,
        dataset_id=dataset,
        module="single_cell.ldsc_celltype",
        genome_build=method.genome_build.value,
        configuration_sha256=digest,
        inputs=completion_inputs,
        outputs=outputs,
        metrics=metrics,
        error_type=SingleCellError,
    )
    logger.record(
        "OUTPUT",
        "ldsc_celltype_results",
        path=str(outputs["normalized_results"]),
        cell_types=metrics["tested_cell_types"],
        reference_files=metrics["reference_files_validated"],
    )
    return LdscCelltypeExecution(outputs=outputs, metrics=metrics, resumed=False)


__all__ = [
    "LdscCelltypeExecution",
    "LdscCelltypeManifestEntry",
    "LdscCelltypePreflight",
    "preflight_ldsc_celltype",
    "run_ldsc_celltype",
    "validate_munged_sumstats",
]
