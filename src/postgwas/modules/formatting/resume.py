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

from .contracts import (
    CUSTOM_OUTPUT_TARGET,
    GCTA_SAMPLE_SIZE_MODE,
    NAMED_OUTPUT_RESULT_KEYS,
    PARTITIONED_OUTPUT_RESULT_KEYS,
    SINGLE_OUTPUT_RESULT_KEYS,
    formatter_result_targets,
)
from .table import FormattingError


_MANIFEST_SCHEMA_VERSION = 1
_NON_PATH_RESULT_KEYS = frozenset({"columns", "field_roles"})


def formatter_resolved_paths(configuration, selected: list[str]) -> list[str]:
    """Return YAML-configured metadata paths for the selected targets."""
    metadata = configuration.modules.formatting.resolved_config
    missing = sorted(set(selected) - set(metadata.format_fields))
    if missing:
        raise FormattingError(
            "No resolved-metadata fields are configured for: %s"
            % ", ".join(missing)
        )
    paths = [
        *metadata.common_fields,
        *(path for target in selected for path in metadata.format_fields[target]),
    ]
    if configuration.modules.formatting.custom_output.active:
        paths.extend(metadata.custom_fields)
    return list(dict.fromkeys(paths))


def formatter_output_paths(
    output_directory: Path,
    dataset_id: str,
    selected: list[str],
    module,
) -> dict[str, Path]:
    """Resolve the canonical configured path for every selected artifact."""
    destinations: dict[str, Path] = {}
    for target in selected:
        schema = module.exports[target]
        if schema.output_file is not None:
            destinations[target] = configured_output_path(
                output_directory,
                schema.output_file,
                error_type=FormattingError,
                dataset_id=dataset_id,
            )
        for name, output_schema in schema.outputs.items():
            destinations["%s.%s" % (target, name)] = configured_output_path(
                output_directory,
                output_schema.output_file,
                error_type=FormattingError,
                dataset_id=dataset_id,
            )
        if schema.partition_file is not None:
            if schema.output_directory is None:
                raise FormattingError(
                    "Partitioned formatter target %s requires output_directory."
                    % target
                )
            partition_root = configured_output_path(
                output_directory,
                schema.output_directory,
                error_type=FormattingError,
                dataset_id=dataset_id,
            )
            for partition in module.chromosomes:
                destinations[
                    "%s[partition=%s]" % (target, partition)
                ] = configured_output_path(
                    partition_root,
                    schema.partition_file,
                    error_type=FormattingError,
                    dataset_id=dataset_id,
                    partition=partition,
                )
    if module.custom_output.active:
        destinations[CUSTOM_OUTPUT_TARGET] = configured_output_path(
            output_directory,
            module.custom_output.output_file,
            error_type=FormattingError,
            dataset_id=dataset_id,
        )
    return destinations


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


def _output_fingerprints(
    results: Mapping[str, Any], selected: list[str], module,
) -> list[dict[str, Any]]:
    artifacts, _ = _recorded_artifact_paths(results, selected, module)
    paths = sorted(set(artifacts.values()), key=str)
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
        "selected_targets": formatter_result_targets(module, selected),
        "results": dict(results),
        "outputs": _output_fingerprints(results, selected, module),
    }
    return write_yaml_report(
        manifest, _manifest_path(output_directory, dataset_id, module),
    )


def _validate_fingerprint(
    record: Mapping[str, Any],
    expected_path: Path,
    *,
    validate_path: bool = True,
) -> dict[str, Any]:
    current = _file_fingerprint(expected_path)
    fields = ("path", "size", "sha256") if validate_path else ("size", "sha256")
    for field in fields:
        if current[field] != record.get(field):
            raise FormattingError(
                "Cannot resume because %s changed (%s mismatch): %s"
                % (expected_path, field, current["path"])
            )
    return current


def _required_result_path(
    result: Mapping[str, Any], key: str, label: str,
) -> Path:
    value = result.get(key)
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise FormattingError(
            "Formatter completion manifest has no reusable path for %s." % label
        )
    path = Path(value).expanduser()
    if not path.is_absolute():
        raise FormattingError(
            "Formatter completion manifest artifact paths must be absolute: %s"
            % label
        )
    return path.resolve()


def _recorded_artifact_paths(
    results: Mapping[str, Any], selected: list[str], module,
) -> tuple[dict[str, Path], dict[str, Path]]:
    targets = formatter_result_targets(module, selected)
    if set(results) != set(targets):
        raise FormattingError(
            "Formatter completion manifest results do not match the selected outputs."
        )
    artifacts: dict[str, Path] = {}
    auxiliary: dict[str, Path] = {}
    for target in targets:
        result = results.get(target)
        if not isinstance(result, Mapping):
            raise FormattingError(
                "Formatter completion manifest has no reusable results for %s."
                % target
            )
        if target == CUSTOM_OUTPUT_TARGET:
            result_key = SINGLE_OUTPUT_RESULT_KEYS[CUSTOM_OUTPUT_TARGET]
            artifacts[target] = _required_result_path(result, result_key, target)
            auxiliary["%s.log_file" % target] = _required_result_path(
                result, "log_file", "%s.log_file" % target,
            )
            continue
        schema = module.exports[target]
        if schema.output_file is not None:
            result_key = SINGLE_OUTPUT_RESULT_KEYS.get(target)
            if result_key is None:
                raise FormattingError(
                    "Formatter resume has no result-path contract for %s." % target
                )
            artifacts[target] = _required_result_path(result, result_key, target)
        for name in schema.outputs:
            label = "%s.%s" % (target, name)
            try:
                result_key = NAMED_OUTPUT_RESULT_KEYS[target][name]
            except KeyError as exc:
                raise FormattingError(
                    "Formatter resume has no result-path contract for %s." % label
                ) from exc
            artifacts[label] = _required_result_path(result, result_key, label)
            detailed_outputs = result.get("outputs")
            if detailed_outputs is not None:
                if not isinstance(detailed_outputs, Mapping):
                    raise FormattingError(
                        "Formatter completion manifest has invalid output details "
                        "for %s." % target
                    )
                detail = detailed_outputs.get(name)
                if detail is not None:
                    if not isinstance(detail, Mapping):
                        raise FormattingError(
                            "Formatter completion manifest has invalid output "
                            "details for %s." % label
                        )
                    detailed_path = _required_result_path(
                        detail, "path", "%s.details" % label,
                    )
                    if detailed_path != artifacts[label]:
                        raise FormattingError(
                            "Cannot resume because manifest artifact %s has "
                            "conflicting recorded paths: %s and %s."
                            % (label, artifacts[label], detailed_path)
                        )
        if schema.partition_file is not None:
            try:
                result_keys = PARTITIONED_OUTPUT_RESULT_KEYS[target]
            except KeyError as exc:
                raise FormattingError(
                    "Formatter resume has no partition-path contract for %s."
                    % target
                ) from exc
            files = result.get(result_keys["files"])
            if not isinstance(files, list) or len(files) != len(module.chromosomes):
                raise FormattingError(
                    "Formatter completion manifest has invalid partition paths for %s."
                    % target
                )
            for partition, value in zip(module.chromosomes, files, strict=True):
                label = "%s[partition=%s]" % (target, partition)
                artifacts[label] = _required_result_path(
                    {"path": value}, "path", label,
                )
            auxiliary["%s.output_directory" % target] = _required_result_path(
                result,
                result_keys["directory"],
                "%s.output_directory" % target,
            )
        auxiliary["%s.log_file" % target] = _required_result_path(
            result, "log_file", "%s.log_file" % target,
        )
    return artifacts, auxiliary


def _recorded_root(
    recorded: Path,
    expected: Path,
    output_directory: Path,
    label: str,
) -> Path:
    relative = expected.relative_to(output_directory)
    relative_parts = relative.parts
    if (
        len(recorded.parts) <= len(relative_parts)
        or recorded.parts[-len(relative_parts):] != relative_parts
    ):
        raise FormattingError(
            "Cannot resume because manifest artifact %s does not match the "
            "currently configured path: recorded %s; expected %s."
            % (label, recorded, expected)
        )
    return Path(*recorded.parts[:-len(relative_parts)])


def _absolute_result_paths(value: Any, *, key: str | None = None):
    if key in _NON_PATH_RESULT_KEYS:
        return
    if isinstance(value, Mapping):
        for child_key, child in value.items():
            yield from _absolute_result_paths(child, key=str(child_key))
    elif isinstance(value, (list, tuple)):
        for child in value:
            yield from _absolute_result_paths(child)
    elif isinstance(value, (str, Path)):
        candidate = Path(value).expanduser()
        if candidate.is_absolute():
            yield candidate.resolve()


def _replace_result_paths(
    value: Any,
    replacements: Mapping[Path, Path],
    *,
    key: str | None = None,
) -> Any:
    if key in _NON_PATH_RESULT_KEYS:
        return value
    if isinstance(value, Mapping):
        return {
            child_key: _replace_result_paths(
                child, replacements, key=str(child_key),
            )
            for child_key, child in value.items()
        }
    if isinstance(value, list):
        return [_replace_result_paths(child, replacements) for child in value]
    if isinstance(value, tuple):
        return tuple(_replace_result_paths(child, replacements) for child in value)
    if isinstance(value, (str, Path)):
        candidate = Path(value).expanduser()
        if candidate.is_absolute():
            replacement = replacements.get(candidate.resolve())
            if replacement is not None:
                return str(replacement)
    return value


def _rebase_manifest_results(
    results: Mapping[str, Any],
    outputs: list[Any],
    *,
    output_directory: Path,
    dataset_id: str,
    selected: list[str],
    module,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    expected = formatter_output_paths(
        output_directory, dataset_id, selected, module,
    )
    recorded, auxiliary = _recorded_artifact_paths(results, selected, module)
    if set(recorded) != set(expected):
        raise FormattingError(
            "Formatter completion manifest artifacts do not match the configured "
            "formatter outputs."
        )

    expected_auxiliary = {
        "%s.log_file" % target: configured_output_path(
            output_directory,
            module.runtime.log_file,
            error_type=FormattingError,
            dataset_id=dataset_id,
        )
        for target in formatter_result_targets(module, selected)
    }
    for target in selected:
        schema = module.exports[target]
        if schema.partition_file is not None:
            expected_auxiliary[
                "%s.output_directory" % target
            ] = configured_output_path(
                output_directory,
                schema.output_directory,
                error_type=FormattingError,
                dataset_id=dataset_id,
            )

    replacements: dict[Path, Path] = {}
    roots: set[Path] = set()
    for label, current in {**expected, **expected_auxiliary}.items():
        prior = (recorded | auxiliary)[label]
        roots.add(_recorded_root(prior, current, output_directory, label))
        if prior in replacements and replacements[prior] != current:
            raise FormattingError(
                "Formatter completion manifest reuses one path for multiple artifacts: "
                "%s" % prior
            )
        replacements[prior] = current
    if len(roots) != 1:
        raise FormattingError(
            "Cannot resume because formatter manifest artifacts come from multiple "
            "run directories: %s" % ", ".join(map(str, sorted(roots, key=str)))
        )

    permitted_paths = set(replacements)
    unexpected = sorted(
        set(_absolute_result_paths(results)) - permitted_paths,
        key=str,
    )
    if unexpected:
        raise FormattingError(
            "Cannot resume because formatter results contain an unconfigured artifact "
            "path: %s" % unexpected[0]
        )

    records_by_path: dict[Path, Mapping[str, Any]] = {}
    for record in outputs:
        if not isinstance(record, Mapping) or not record.get("path"):
            raise FormattingError("Formatter completion manifest has an invalid output.")
        path = Path(str(record["path"])).expanduser()
        if not path.is_absolute():
            raise FormattingError(
                "Formatter completion manifest output fingerprints must use absolute "
                "paths."
            )
        resolved = path.resolve()
        if resolved in records_by_path:
            raise FormattingError(
                "Formatter completion manifest repeats an output fingerprint: %s"
                % resolved
            )
        records_by_path[resolved] = record
    if set(records_by_path) != set(recorded.values()):
        raise FormattingError(
            "Formatter completion manifest output fingerprints do not match its "
            "result artifacts."
        )
    validated_outputs = [
        _validate_fingerprint(
            records_by_path[prior], expected[label], validate_path=False,
        )
        for label, prior in recorded.items()
    ]

    rebased = _replace_result_paths(results, replacements)
    current_root = output_directory.resolve()
    outside = sorted({
        path
        for path in _absolute_result_paths(rebased)
        if path != current_root and current_root not in path.parents
    }, key=str)
    if outside:
        raise FormattingError(
            "Cannot resume because a returned formatter artifact is outside the "
            "current output directory: %s" % outside[0]
        )
    return dict(rebased), sorted(validated_outputs, key=lambda item: item["path"])


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
    recorded_targets = manifest.get("selected_targets")
    current_targets = formatter_result_targets(
        configuration.modules.formatting, selected,
    )
    if recorded_targets is not None and recorded_targets != current_targets:
        raise FormattingError(
            "Formatter completion manifest selected outputs do not match."
        )
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
    results = manifest.get("results")
    if not isinstance(results, dict):
        raise FormattingError("Formatter completion manifest has no reusable results.")
    rebased, rebased_outputs = _rebase_manifest_results(
        results,
        outputs,
        output_directory=output_directory,
        dataset_id=dataset_id,
        selected=selected,
        module=configuration.modules.formatting,
    )
    resumed = {
        target: dict(rebased[target])
        for target in formatter_result_targets(
            configuration.modules.formatting, selected,
        )
    }
    for result in resumed.values():
        result["resumed"] = True
    migrated = dict(manifest)
    migrated["configuration_sha256"] = current_digest
    migrated["configuration_scope"] = "selected_formatter_targets"
    migrated["selected_targets"] = current_targets
    migrated["results"] = rebased
    migrated["outputs"] = rebased_outputs
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
        configuration.modules.formatting.custom_output.active
        or len(selected) != 1
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
    "formatter_output_paths",
    "formatter_resolved_paths",
    "resume_formatter_outputs",
    "write_formatter_completion_manifest",
]
