"""Generic checksummed completion manifests for atomic scientific modules."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Literal, Mapping, Type

import yaml

from postgwas.core.io.reports import write_yaml_report
from postgwas.core.resource_preparation import sha256


MANIFEST_SCHEMA_VERSION = 1


@dataclass(frozen=True)
class CompletionResumeDecision:
    """Validated decision for one completed module's declared outputs."""

    action: Literal["resume", "restart"]
    manifest: Mapping[str, Any]
    reason: str | None = None
    warning: str | None = None
    missing_outputs: tuple[str, ...] = ()


@dataclass(frozen=True)
class CompletionOutputDecision:
    """Availability decision made before output fingerprint validation."""

    action: Literal["validate", "restart"]
    missing_outputs: tuple[str, ...] = ()


def assess_completion_outputs(
    outputs: Mapping[str, str | Path],
    *,
    error_type: Type[Exception] = RuntimeError,
) -> CompletionOutputDecision:
    """Allow restart only when none of the declared output paths remain."""
    if not outputs:
        raise error_type("Completion output contract contains no declared files")
    resolved = {
        name: Path(value).expanduser().resolve()
        for name, value in outputs.items()
    }
    present = {
        name: path.exists() or path.is_symlink()
        for name, path in resolved.items()
    }
    missing = tuple(name for name, exists in present.items() if not exists)
    if len(missing) == len(resolved):
        return CompletionOutputDecision(
            action="restart",
            missing_outputs=missing,
        )
    if missing:
        existing = tuple(name for name, exists in present.items() if exists)
        raise error_type(
            "Cannot resume because the declared output set is partial; existing: "
            "%s; missing: %s. Use --overwrite after reviewing the surviving "
            "outputs."
            % (", ".join(existing), ", ".join(missing))
        )
    return CompletionOutputDecision(action="validate")


def _policy_action(policy, field: str) -> str:
    """Read one schema-validated policy action without coupling to config models."""
    return str(getattr(policy, field, "error"))


def completion_restart_decision(
    document: Mapping[str, Any],
    *,
    policy,
    policy_field: str,
    reason: str,
    warning: str,
    error_type: Type[Exception],
    missing_outputs: tuple[str, ...] = (),
) -> CompletionResumeDecision:
    """Apply one validated global policy action to a checkpoint mismatch."""
    if _policy_action(policy, policy_field) != "warn_and_restart":
        raise error_type(warning)
    return CompletionResumeDecision(
        action="restart",
        manifest=document,
        reason=reason,
        warning=warning,
        missing_outputs=missing_outputs,
    )


def record_completion_restart(
    logger,
    operation: str,
    *,
    missing_outputs=(),
    manifest: str | Path,
    decision: CompletionResumeDecision | None = None,
) -> None:
    """Record the canonical reason for regenerating a stale completed step."""
    reason = decision.reason if decision is not None else None
    warning = decision.warning if decision is not None else None
    if warning:
        logger.warning(warning)
    logger.record(
        "ACTION",
        operation,
        action="restart",
        reason=reason or "all_declared_outputs_missing",
        missing_outputs=(
            decision.missing_outputs if decision is not None
            else missing_outputs
        ),
        manifest=str(manifest),
    )


def apply_completion_restart(
    decision: CompletionResumeDecision,
    *,
    output_root: str | Path,
    manifest: str | Path,
    logger,
    operation: str,
    error_type: Type[Exception] = RuntimeError,
) -> tuple[Path, ...]:
    """Remove only intact, manifest-owned outputs for a validated restart."""
    if decision.action != "restart":
        raise error_type("Completion restart can be applied only to restart decisions")
    root = Path(output_root).expanduser().resolve()
    manifest_path = Path(manifest).expanduser()
    resolved_manifest = manifest_path.resolve()
    if (
        resolved_manifest == root
        or root not in resolved_manifest.parents
        or manifest_path.is_symlink()
    ):
        raise error_type(
            "Refusing to replace a completion manifest outside its configured "
            "output root: %s" % resolved_manifest
        )
    recorded_outputs = decision.manifest.get("outputs")
    if isinstance(recorded_outputs, Mapping):
        output_records = recorded_outputs.items()
    elif isinstance(recorded_outputs, list):
        output_records = (
            ("output_%d" % number, record)
            for number, record in enumerate(recorded_outputs, 1)
        )
    else:
        output_records = ()
    output_records = tuple(output_records)
    if not output_records:
        raise error_type("Completion manifest contains no replaceable output records")

    removable: list[Path] = []
    for name, record in output_records:
        if not isinstance(record, Mapping) or not isinstance(record.get("path"), str):
            raise error_type(
                "Completion manifest output %s fingerprint is invalid" % name
            )
        raw_path = Path(record["path"]).expanduser()
        path = raw_path.resolve()
        if path == root or root not in path.parents or raw_path.is_symlink():
            raise error_type(
                "Refusing to replace completion output %s outside its configured "
                "output root or through a symlink: %s" % (name, path)
            )
        if not raw_path.exists():
            continue
        if not raw_path.is_file():
            raise error_type(
                "Refusing automatic restart because output %s is not a regular "
                "file: %s" % (name, path)
            )
        current = {
            "path": str(path),
            "size": raw_path.stat().st_size,
            "sha256": sha256(raw_path),
        }
        for field in ("size", "sha256"):
            if record.get(field) != current[field]:
                raise error_type(
                    "Refusing automatic restart because output %s was modified "
                    "after its checkpoint (%s mismatch): %s"
                    % (name, field, path)
                )
        removable.append(raw_path)

    record_completion_restart(
        logger,
        operation,
        manifest=resolved_manifest,
        decision=decision,
    )
    for path in removable:
        path.unlink()
    if manifest_path.exists():
        manifest_path.unlink()
    return tuple(path.resolve() for path in removable)


def file_fingerprint(
    path: str | Path,
    *,
    error_type: Type[Exception] = RuntimeError,
) -> dict[str, Any]:
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
        "inputs": {
            name: file_fingerprint(value, error_type=error_type)
            for name, value in inputs.items()
        },
        "outputs": {
            name: file_fingerprint(value, error_type=error_type)
            for name, value in outputs.items()
        },
        "metrics": dict(metrics),
    }
    return write_yaml_report(document, path)


def resolve_completion_resume(
    path: str | Path,
    *,
    dataset_id: str,
    module: str,
    genome_build: str,
    configuration_sha256: str,
    inputs: Mapping[str, str | Path],
    outputs: Mapping[str, str | Path],
    resume_policy=None,
    error_type: Type[Exception] = RuntimeError,
) -> CompletionResumeDecision:
    """Validate provenance and decide whether to resume or safely restart.

    Policy-controlled restarts are returned for changed parameters, changed
    inputs, changed output contracts, and missing checkpoint outputs. A checksum
    mismatch in an existing output is always an error because that file is no
    longer provably identical to the PostGWAS-owned checkpoint artifact.
    """
    manifest_path = Path(path)
    try:
        document = yaml.safe_load(manifest_path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise error_type(
            "Cannot read completion manifest %s: %s" % (manifest_path, exc)
        ) from exc
    identity = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "status": "COMPLETED",
        "dataset_id": dataset_id,
        "module": module,
        "genome_build": genome_build,
    }
    for field, value in identity.items():
        if document.get(field) != value:
            raise error_type(
                "Completion manifest %s does not match %s" % (field, value)
            )
    if document.get("configuration_sha256") != configuration_sha256:
        return completion_restart_decision(
            document,
            policy=resume_policy,
            policy_field="changed_parameters",
            reason="changed_parameters",
            warning=(
                "Resolved parameters changed since this checkpoint; the previous "
                "PostGWAS-owned results will be replaced and the run restarted."
            ),
            error_type=error_type,
        )
    recorded_inputs = document.get("inputs")
    if (
        not isinstance(recorded_inputs, Mapping)
        or set(recorded_inputs) != set(inputs)
    ):
        return completion_restart_decision(
            document,
            policy=resume_policy,
            policy_field="changed_inputs",
            reason="changed_inputs",
            warning=(
                "The configured input set changed since this checkpoint; the "
                "previous PostGWAS-owned results will be replaced and the run "
                "restarted."
            ),
            error_type=error_type,
        )
    for name, value in inputs.items():
        if not isinstance(recorded_inputs[name], Mapping):
            raise error_type(
                "Completion manifest input %s fingerprint is invalid" % name
            )
        current = file_fingerprint(value, error_type=error_type)
        for field in ("path", "size", "sha256"):
            if recorded_inputs[name].get(field) != current[field]:
                return completion_restart_decision(
                    document,
                    policy=resume_policy,
                    policy_field="changed_inputs",
                    reason="changed_inputs",
                    warning=(
                        "Input %s changed since this checkpoint (%s mismatch); "
                        "the previous PostGWAS-owned results will be replaced and "
                        "the run restarted." % (name, field)
                    ),
                    error_type=error_type,
                )

    recorded_outputs = document.get("outputs")
    if (
        not isinstance(recorded_outputs, Mapping)
        or set(recorded_outputs) != set(outputs)
    ):
        return completion_restart_decision(
            document,
            policy=resume_policy,
            policy_field="unvalidated_outputs",
            reason="changed_output_contract",
            warning=(
                "The expected output contract changed since this checkpoint; "
                "the previous PostGWAS-owned results will be replaced and the "
                "run restarted."
            ),
            error_type=error_type,
        )
    resolved_outputs = {
        name: Path(value).expanduser().resolve()
        for name, value in outputs.items()
    }
    for name, current_path in resolved_outputs.items():
        recorded = recorded_outputs[name]
        if not isinstance(recorded, Mapping):
            raise error_type(
                "Completion manifest output %s fingerprint is invalid" % name
            )
        if recorded.get("path") != str(current_path):
            return completion_restart_decision(
                document,
                policy=resume_policy,
                policy_field="unvalidated_outputs",
                reason="changed_output_contract",
                warning=(
                    "Output %s changed since this checkpoint (path mismatch); "
                    "the previous PostGWAS-owned results will be replaced and "
                    "the run restarted." % name
                ),
                error_type=error_type,
            )

    present = {
        name: path.exists() or path.is_symlink()
        for name, path in resolved_outputs.items()
    }
    missing = tuple(name for name, exists in present.items() if not exists)
    if missing:
        return completion_restart_decision(
            document,
            policy=resume_policy,
            policy_field="unvalidated_outputs",
            reason="incomplete_outputs",
            warning=(
                "The checkpoint output set is incomplete (missing: %s); no later "
                "stage will trust it. PostGWAS will restart from the earliest "
                "safe checkpoint and replace only its owned outputs."
                % ", ".join(missing)
            ),
            missing_outputs=missing,
            error_type=error_type,
        )

    for name, value in resolved_outputs.items():
        current = file_fingerprint(value, error_type=error_type)
        for field in ("size", "sha256"):
            if recorded_outputs[name].get(field) != current[field]:
                raise error_type(
                    "Refusing automatic restart because output %s changed since "
                    "its checkpoint (%s mismatch): %s. Review the file and use "
                    "--overwrite only when replacement is intended."
                    % (name, field, value)
                )
    return CompletionResumeDecision(action="resume", manifest=document)


def validate_completion_manifest(
    path: str | Path,
    *,
    dataset_id: str,
    module: str,
    genome_build: str,
    configuration_sha256: str,
    inputs: Mapping[str, str | Path],
    outputs: Mapping[str, str | Path],
    resume_policy=None,
    error_type: Type[Exception] = RuntimeError,
) -> Mapping[str, Any]:
    """Return a fully reusable manifest; reject a restart-only decision."""
    decision = resolve_completion_resume(
        path,
        dataset_id=dataset_id,
        module=module,
        genome_build=genome_build,
        configuration_sha256=configuration_sha256,
        inputs=inputs,
        outputs=outputs,
        resume_policy=resume_policy,
        error_type=error_type,
    )
    if decision.action == "restart":
        raise error_type(
            "The completion checkpoint requires a policy-controlled restart "
            "and cannot be returned as a reusable manifest."
        )
    return decision.manifest


__all__ = [
    "CompletionOutputDecision", "CompletionResumeDecision",
    "apply_completion_restart", "assess_completion_outputs",
    "completion_restart_decision", "configuration_digest", "file_fingerprint",
    "record_completion_restart",
    "resolve_completion_resume",
    "validate_completion_manifest", "write_completion_manifest",
]
