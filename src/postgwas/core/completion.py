"""Generic checksummed completion manifests for atomic scientific modules."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Type

import yaml

from postgwas.core.io.reports import write_yaml_report
from postgwas.core.resource_preparation import sha256


MANIFEST_SCHEMA_VERSION = 1


def file_fingerprint(path: str | Path, *, error_type: Type[Exception] = RuntimeError) -> dict[str, Any]:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size <= 0:
        raise error_type("Completion file is missing or empty: %s" % resolved)
    return {"path": str(resolved), "size": resolved.stat().st_size, "sha256": sha256(resolved)}


def configuration_digest(value: Any) -> str:
    if hasattr(value, "model_dump"):
        value = value.model_dump(mode="json")
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def write_completion_manifest(
    path: str | Path,
    *,
    dataset_id: str,
    module: str,
    genome_build: str,
    configuration_sha256: str,
    inputs: Mapping[str, str | Path],
    outputs: Mapping[str, str | Path],
    metrics: Mapping[str, Any],
    error_type: Type[Exception] = RuntimeError,
) -> Path:
    document = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "status": "COMPLETED",
        "dataset_id": dataset_id,
        "module": module,
        "genome_build": genome_build,
        "configuration_sha256": configuration_sha256,
        "inputs": {name: file_fingerprint(value, error_type=error_type) for name, value in inputs.items()},
        "outputs": {name: file_fingerprint(value, error_type=error_type) for name, value in outputs.items()},
        "metrics": dict(metrics),
    }
    return write_yaml_report(document, path)


def validate_completion_manifest(
    path: str | Path,
    *,
    dataset_id: str,
    module: str,
    genome_build: str,
    configuration_sha256: str,
    inputs: Mapping[str, str | Path],
    outputs: Mapping[str, str | Path],
    error_type: Type[Exception] = RuntimeError,
) -> Mapping[str, Any]:
    manifest_path = Path(path)
    try:
        document = yaml.safe_load(manifest_path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise error_type("Cannot read completion manifest %s: %s" % (manifest_path, exc)) from exc
    expected = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "status": "COMPLETED",
        "dataset_id": dataset_id,
        "module": module,
        "genome_build": genome_build,
        "configuration_sha256": configuration_sha256,
    }
    for field, value in expected.items():
        if document.get(field) != value:
            raise error_type("Completion manifest %s does not match %s" % (field, value))
    for group_name, current_paths in (("inputs", inputs), ("outputs", outputs)):
        recorded = document.get(group_name)
        if not isinstance(recorded, Mapping) or set(recorded) != set(current_paths):
            raise error_type("Completion manifest %s set does not match" % group_name)
        for name, value in current_paths.items():
            if not isinstance(recorded[name], Mapping):
                raise error_type(
                    "Completion manifest %s %s fingerprint is invalid"
                    % (group_name[:-1], name)
                )
            current = file_fingerprint(value, error_type=error_type)
            for field in ("path", "size", "sha256"):
                if recorded[name].get(field) != current[field]:
                    raise error_type(
                        "Cannot resume because %s %s changed (%s mismatch)"
                        % (group_name[:-1], name, field)
                    )
    return document


__all__ = [
    "configuration_digest", "file_fingerprint", "validate_completion_manifest",
    "write_completion_manifest",
]
