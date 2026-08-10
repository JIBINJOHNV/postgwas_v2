"""Provenance-validated formatter completion and resume support."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import yaml

from postgwas.config import (
    resolved_configuration_values,
    select_configuration_values,
)
from postgwas.core.io.delimiters import open_text
from postgwas.core.io.reports import write_yaml_report
from postgwas.core.paths import configured_output_path
from postgwas.core.resource_preparation import sha256

from .contracts import GCTA_SAMPLE_SIZE_MODE, NAMED_OUTPUT_RESULT_KEYS
from .table import FormattingError


_MANIFEST_SCHEMA_VERSION = 1


def formatter_resolved_paths(configuration, selected: list[str]) -> list[str]:
    """Return YAML-configured metadata paths for the selected targets."""
    metadata = configuration.modules.formatting.resolved_config
    missing = sorted(set(selected) - set(metadata.format_fields))
    if missing:
        raise FormattingError(
            "No resolved-metadata fields are configured for: %s"
            % ", ".join(missing)
        )
    return list(dict.fromkeys([
        *metadata.common_fields,
        *(path for target in selected for path in metadata.format_fields[target]),
    ]))


def _file_fingerprint(path: Path) -> dict[str, Any]:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size <= 0:
        raise FormattingError("Resume file is missing or empty: %s" % resolved)
    return {
        "path": str(resolved),
        "size": resolved.stat().st_size,
        "sha256": sha256(resolved),
    }


def _configuration_values_digest(
    formatting: Mapping[str, Any], bcftools: object, selected: list[str],
) -> str:
    payload = {
        "formats": selected,
        "formatting": dict(formatting),
        "bcftools": str(bcftools),
    }
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _configuration_digest(configuration, selected: list[str]) -> str:
    formatting = resolved_configuration_values(
        configuration,
        modules=("formatting",),
        module_paths={
            "formatting": formatter_resolved_paths(configuration, selected),
        },
    )["modules"]["formatting"]
    return _configuration_values_digest(
        formatting, configuration.resources.executables.bcftools, selected,
    )


def _result_files(value: Any, *, key: str | None = None):
    if key == "log_file":
        return
    if isinstance(value, Mapping):
        for child_key, child in value.items():
            yield from _result_files(child, key=str(child_key))
    elif isinstance(value, (list, tuple)):
        for child in value:
            yield from _result_files(child)
    elif isinstance(value, (str, Path)):
        candidate = Path(value).expanduser()
        if candidate.is_file():
            yield candidate.resolve()
        elif candidate.is_dir():
            yield from (
                path.resolve() for path in sorted(candidate.rglob("*"))
                if path.is_file()
            )


def _output_fingerprints(results: Mapping[str, Any]) -> list[dict[str, Any]]:
    paths = sorted(set(_result_files(results)), key=str)
    if not paths:
        raise FormattingError("Formatter completion has no output files to validate.")
    return [_file_fingerprint(path) for path in paths]


def _manifest_path(output_directory: Path, dataset_id: str, module) -> Path:
    return configured_output_path(
        output_directory,
        module.runtime.completion_manifest_file,
        error_type=FormattingError,
        dataset_id=dataset_id,
    )


def write_formatter_completion_manifest(
    *, output_directory: Path, dataset_id: str, vcf: Path,
    selected: list[str], configuration, results: Mapping[str, Any],
) -> Path:
    """Atomically record the exact input, configuration, results, and outputs."""
    module = configuration.modules.formatting
    manifest = {
        "schema_version": _MANIFEST_SCHEMA_VERSION,
        "status": "COMPLETED",
        "dataset_id": dataset_id,
        "input_vcf": _file_fingerprint(vcf),
        "configuration_sha256": _configuration_digest(configuration, selected),
        "configuration_scope": "selected_formatter_targets",
        "formats": selected,
        "results": dict(results),
        "outputs": _output_fingerprints(results),
    }
    return write_yaml_report(
        manifest, _manifest_path(output_directory, dataset_id, module),
    )


def _validate_fingerprint(record: Mapping[str, Any], expected_path: Path) -> None:
    current = _file_fingerprint(expected_path)
    for field in ("path", "size", "sha256"):
        if current[field] != record.get(field):
            raise FormattingError(
                "Cannot resume because %s changed (%s mismatch): %s"
                % (expected_path, field, current["path"])
            )


def _load_manifest(path: Path) -> Mapping[str, Any]:
    try:
        value = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise FormattingError("Cannot read formatter completion manifest: %s" % exc) from exc
    if not isinstance(value, Mapping):
        raise FormattingError("Formatter completion manifest must be a mapping: %s" % path)
    return value


def _validate_manifest(
    manifest: Mapping[str, Any], *, output_directory: Path,
    dataset_id: str, vcf: Path,
    selected: list[str], configuration,
) -> dict[str, Any]:
    if manifest.get("schema_version") != _MANIFEST_SCHEMA_VERSION:
        raise FormattingError("Unsupported formatter completion manifest schema.")
    if manifest.get("status") != "COMPLETED":
        raise FormattingError("Formatter completion manifest is not complete.")
    if manifest.get("dataset_id") != dataset_id:
        raise FormattingError("Formatter completion manifest dataset does not match.")
    if manifest.get("formats") != selected:
        raise FormattingError("Formatter completion manifest formats do not match.")
    recorded_digest = manifest.get("configuration_sha256")
    current_digest = _configuration_digest(configuration, selected)
    migrated_digest = recorded_digest != current_digest
    if migrated_digest:
        _matching_recorded_configuration(
            output_directory,
            dataset_id,
            selected,
            configuration,
            expected_digest=recorded_digest,
        )
    input_record = manifest.get("input_vcf")
    if not isinstance(input_record, Mapping):
        raise FormattingError("Formatter completion manifest has no input fingerprint.")
    _validate_fingerprint(input_record, vcf)
    outputs = manifest.get("outputs")
    if not isinstance(outputs, list) or not outputs:
        raise FormattingError("Formatter completion manifest has no output fingerprints.")
    for record in outputs:
        if not isinstance(record, Mapping) or not record.get("path"):
            raise FormattingError("Formatter completion manifest has an invalid output.")
        _validate_fingerprint(record, Path(str(record["path"])))
    results = manifest.get("results")
    if not isinstance(results, dict):
        raise FormattingError("Formatter completion manifest has no reusable results.")
    resumed = dict(results)
    for result in resumed.values():
        if isinstance(result, dict):
            result["resumed"] = True
    if migrated_digest:
        migrated = dict(manifest)
        migrated["configuration_sha256"] = current_digest
        migrated["configuration_scope"] = "selected_formatter_targets"
        write_yaml_report(
            migrated,
            _manifest_path(
                output_directory,
                dataset_id,
                configuration.modules.formatting,
            ),
        )
    return resumed


def _completed_log_block(log_path: Path, vcf: Path) -> bool:
    if not log_path.is_file():
        return False
    lines = log_path.read_text(encoding="utf-8").splitlines()
    completed = [
        index for index, line in enumerate(lines)
        if "formatter_run status=COMPLETED" in line
    ]
    if not completed:
        return False
    end = completed[-1]
    starts = [
        index for index, line in enumerate(lines[: end + 1])
        if "INPUT    formatter_run" in line
    ]
    if not starts:
        return False
    block = "\n".join(lines[starts[-1] : end + 1])
    return "vcf=%s" % vcf.expanduser().resolve() in block


def _table_rows_and_header(
    path: Path, expected: list[str], delimiter: str,
) -> tuple[int, list[str]]:
    with open_text(path) as handle:
        header = handle.readline().rstrip("\r\n").split(delimiter)
        rows = 0
        for line_number, line in enumerate(handle, 2):
            if not line.strip():
                continue
            if len(line.rstrip("\r\n").split(delimiter)) != len(header):
                raise FormattingError(
                    "Existing formatter output has an incomplete record at line "
                    "%s: %s" % (line_number, path)
                )
            rows += 1
    if header != expected or rows < 1:
        raise FormattingError(
            "Existing formatter output failed schema validation: %s" % path
        )
    return rows, header


def _matching_recorded_configuration(
    output_directory: Path, dataset_id: str, selected: list[str], configuration,
    *, expected_digest: object = None,
) -> bool:
    module = configuration.modules.formatting
    resolved_path = configured_output_path(
        output_directory, module.runtime.resolved_config_file,
        error_type=FormattingError, dataset_id=dataset_id,
    )
    if not resolved_path.is_file():
        return False
    try:
        resolved = yaml.safe_load(resolved_path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise FormattingError(
            "Cannot validate the prior formatter configuration: %s" % exc
        ) from exc
    recorded_module = resolved.get("modules", {}).get("formatting")
    if not isinstance(recorded_module, Mapping):
        raise FormattingError(
            "The prior resolved configuration has no formatter section."
        )
    paths = formatter_resolved_paths(configuration, selected)
    recorded_values = select_configuration_values(recorded_module, paths)
    current_module = resolved_configuration_values(
        configuration,
        modules=("formatting",),
        module_paths={"formatting": paths},
    )["modules"]["formatting"]
    recorded_values.get("runtime", {}).pop("completion_manifest_file", None)
    current_module.get("runtime", {}).pop("completion_manifest_file", None)
    if recorded_values != current_module:
        raise FormattingError(
            "Existing formatter outputs use a different resolved configuration."
        )
    if expected_digest is not None:
        try:
            recorded_bcftools = resolved["resources"]["executables"]["bcftools"]
        except (KeyError, TypeError) as exc:
            raise FormattingError(
                "The prior resolved configuration has no bcftools resource."
            ) from exc
        digest_values = dict(recorded_module)
        digest_without_projection = dict(digest_values)
        digest_without_projection.pop("resolved_config", None)
        compatible_digests = {
            _configuration_values_digest(
                digest_values, recorded_bcftools, selected,
            ),
            _configuration_values_digest(
                digest_without_projection, recorded_bcftools, selected,
            ),
        }
        recorded_digest = str(expected_digest or "")
        if recorded_digest not in compatible_digests:
            raise FormattingError(
                "Formatter completion metadata does not match its recorded "
                "resolved configuration."
            )
    return True


def _adopt_existing_named_outputs(
    *, output_directory: Path, dataset_id: str, vcf: Path,
    selected: list[str], configuration, log_path: Path,
) -> dict[str, Any] | None:
    if (
        len(selected) != 1
        or selected[0] not in NAMED_OUTPUT_RESULT_KEYS
        or not _completed_log_block(log_path, vcf)
        or not _matching_recorded_configuration(
            output_directory, dataset_id, selected, configuration,
        )
    ):
        return None
    target = selected[0]
    result_keys = NAMED_OUTPUT_RESULT_KEYS[target]
    module = configuration.modules.formatting
    schema = module.exports[target]
    outputs = {
        name: configured_output_path(
            output_directory, output_schema.output_file,
            error_type=FormattingError, dataset_id=dataset_id,
        )
        for name, output_schema in schema.outputs.items()
    }
    if not all(path.is_file() for path in outputs.values()):
        return None
    if vcf.stat().st_mtime_ns > min(
        path.stat().st_mtime_ns for path in outputs.values()
    ):
        raise FormattingError("Input VCF is newer than existing formatter outputs.")
    validated = {
        name: _table_rows_and_header(
            path,
            list(schema.outputs[name].columns.values()),
            module.runtime.table_delimiter,
        )
        for name, path in outputs.items()
    }
    row_counts = [value[0] for value in validated.values()]
    if schema.validation.required_columns and len(set(row_counts)) != 1:
        raise FormattingError(
            "Existing %s formatter outputs contain different variant counts."
            % target
        )
    result = {
        target: {
            **{
                result_keys[name]: str(path)
                for name, path in outputs.items()
            },
            "rows_in": max(row_counts),
            "rows_out": min(row_counts),
            "rows_excluded": None,
            "p_values_bounded": None,
            "columns": {name: value[1] for name, value in validated.items()},
            "log_file": str(log_path),
            "resumed": True,
            "adopted_existing_outputs": True,
        }
    }
    if target == "gcta_gene":
        result[target]["sample_size_mode"] = GCTA_SAMPLE_SIZE_MODE
    write_formatter_completion_manifest(
        output_directory=output_directory,
        dataset_id=dataset_id,
        vcf=vcf,
        selected=selected,
        configuration=configuration,
        results=result,
    )
    return result


def resume_formatter_outputs(
    *, output_directory: Path, dataset_id: str, vcf: Path,
    selected: list[str], configuration, log_path: Path,
) -> dict[str, Any] | None:
    """Return validated prior formatter results, or ``None`` when none exist."""
    module = configuration.modules.formatting
    path = _manifest_path(output_directory, dataset_id, module)
    if path.is_file():
        return _validate_manifest(
            _load_manifest(path), output_directory=output_directory,
            dataset_id=dataset_id, vcf=vcf,
            selected=selected, configuration=configuration,
        )
    return _adopt_existing_named_outputs(
        output_directory=output_directory, dataset_id=dataset_id, vcf=vcf,
        selected=selected, configuration=configuration, log_path=log_path,
    )


__all__ = [
    "formatter_resolved_paths",
    "resume_formatter_outputs",
    "write_formatter_completion_manifest",
]
