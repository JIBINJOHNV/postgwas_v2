"""Validated execution checkpoints shared by direct commands and pipelines."""

from __future__ import annotations

from argparse import Namespace
from dataclasses import dataclass
from datetime import datetime, timezone
from enum import Enum
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path
import platform
import shutil
import sys
from typing import Any, Literal, Mapping, Sequence, Type

from rich.console import Console
import yaml

from postgwas.core.completion import (
    apply_completion_restart,
    completion_restart_decision,
    configuration_digest,
    file_fingerprint,
)
from postgwas.core.contracts import Artifact, ModuleResult
from postgwas.core.io.reports import write_yaml_report
from postgwas.core.resource_preparation import sha256


EXECUTION_CHECKPOINT_SCHEMA_VERSION = 1
_TYPE_KEY = "__postgwas_type__"
_RUN_CONTROL_DESTINATIONS = frozenset({"resume", "overwrite"})
_RUN_CONTROL_OPTIONS = frozenset({"--resume", "--no-resume", "--overwrite"})


class NoOwnedCheckpointOutputs(RuntimeError):
    """A successful legacy command exposed no safely attributable output files."""


@dataclass(frozen=True)
class CheckpointDecision:
    """One pre-execution checkpoint decision and its validated document."""

    action: Literal["run", "resume"]
    document: Mapping[str, Any] | None = None


class CheckpointAuditLogger:
    """Append-only global checkpoint audit log with terminal warnings."""

    def __init__(self, path: str | Path, *, console: Console | None = None):
        self.path = Path(path).expanduser().resolve()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.console = console or Console()

    def _write(self, marker: str, subject: str, values: Mapping[str, Any]) -> None:
        stamp = datetime.now(timezone.utc).isoformat()
        detail = ""
        if values:
            detail = " " + json.dumps(
                dict(values), sort_keys=True, ensure_ascii=False, default=str,
            )
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write("[%s] %s %s%s\n" % (stamp, marker, subject, detail))
            handle.flush()

    def warning(self, message: str) -> None:
        self._write("WARNING", "checkpoint_restart", {"message": message})
        self.console.print("\n[bold yellow]WARNING:[/bold yellow] %s" % message)

    warn = warning

    def record(self, marker: str, subject: str, **values: Any) -> None:
        self._write(str(marker).upper(), str(subject), values)


def encode_checkpoint_value(value: Any) -> Any:
    """Encode supported pipeline state as safe YAML data without pickle."""
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, Path):
        return {_TYPE_KEY: "path", "value": str(value)}
    if isinstance(value, Enum):
        return {
            _TYPE_KEY: "enum_value",
            "value": encode_checkpoint_value(value.value),
        }
    if isinstance(value, Artifact):
        return {
            _TYPE_KEY: "artifact",
            "kind": value.kind,
            "path": str(value.path),
            "metadata": encode_checkpoint_value(dict(value.metadata)),
        }
    if isinstance(value, ModuleResult):
        return {
            _TYPE_KEY: "module_result",
            "module": value.module,
            "artifacts": encode_checkpoint_value(dict(value.artifacts)),
            "metrics": encode_checkpoint_value(dict(value.metrics)),
            "warnings": encode_checkpoint_value(value.warnings),
        }
    if isinstance(value, tuple):
        return {
            _TYPE_KEY: "tuple",
            "items": [encode_checkpoint_value(item) for item in value],
        }
    if isinstance(value, list):
        return [encode_checkpoint_value(item) for item in value]
    if isinstance(value, Mapping):
        if not all(isinstance(key, str) for key in value):
            raise TypeError("Checkpoint mappings require string keys")
        return {
            key: encode_checkpoint_value(item)
            for key, item in value.items()
        }
    item = getattr(value, "item", None)
    if callable(item):
        return encode_checkpoint_value(item())
    raise TypeError(
        "Checkpoint state contains unsupported value %s" % type(value).__name__
    )


def decode_checkpoint_value(value: Any) -> Any:
    """Decode state written by :func:`encode_checkpoint_value`."""
    if isinstance(value, list):
        return [decode_checkpoint_value(item) for item in value]
    if not isinstance(value, Mapping):
        return value
    kind = value.get(_TYPE_KEY)
    if kind == "path":
        return Path(value["value"])
    if kind == "enum_value":
        return decode_checkpoint_value(value["value"])
    if kind == "tuple":
        return tuple(decode_checkpoint_value(item) for item in value["items"])
    if kind == "artifact":
        return Artifact(
            kind=str(value["kind"]),
            path=Path(value["path"]),
            metadata=decode_checkpoint_value(value["metadata"]),
        )
    if kind == "module_result":
        return ModuleResult(
            module=str(value["module"]),
            artifacts=decode_checkpoint_value(value["artifacts"]),
            metrics=decode_checkpoint_value(value["metrics"]),
            warnings=tuple(decode_checkpoint_value(value["warnings"])),
        )
    if kind is not None:
        raise ValueError("Unknown checkpoint state type: %s" % kind)
    return {
        str(key): decode_checkpoint_value(item)
        for key, item in value.items()
    }


def checkpoint_namespace(namespace: Namespace) -> dict[str, Any]:
    """Return the public, safely serialisable state of one CLI namespace."""
    return {
        key: value for key, value in vars(namespace).items()
        if not key.startswith("_")
    }


def checkpoint_content_arguments(arguments: Sequence[str]) -> list[str]:
    """Remove execution controls that cannot change scientific run content."""
    return [
        str(argument) for argument in arguments
        if argument not in _RUN_CONTROL_OPTIONS
    ]


def checkpoint_content_configuration(configuration: Any) -> dict[str, Any]:
    """Return resolved configuration without resume/overwrite run controls."""
    document = configuration.model_dump(mode="json")
    run = document.get("run")
    if isinstance(run, dict):
        for destination in _RUN_CONTROL_DESTINATIONS:
            run.pop(destination, None)
    return document


def checkpoint_content_namespace(namespace: Namespace) -> dict[str, Any]:
    """Return public CLI state without resume/overwrite run controls."""
    return {
        key: value
        for key, value in checkpoint_namespace(namespace).items()
        if key not in _RUN_CONTROL_DESTINATIONS
    }


def preserve_checkpoint_run_controls(namespace: Namespace) -> dict[str, Any]:
    """Capture controls that must survive restoration of prior pipeline state."""
    return {
        key: getattr(namespace, key)
        for key in _RUN_CONTROL_DESTINATIONS
        if hasattr(namespace, key)
    }


def restore_checkpoint_namespace(namespace: Namespace, values: Mapping[str, Any]) -> None:
    """Restore public CLI state after a pipeline stage is safely skipped."""
    for key in tuple(vars(namespace)):
        if not key.startswith("_"):
            delattr(namespace, key)
    for key, value in values.items():
        setattr(namespace, key, value)


def _iter_scalars(value: Any, *, key: str = ""):
    if hasattr(value, "model_dump"):
        value = value.model_dump(mode="python")
    if isinstance(value, Namespace):
        value = checkpoint_namespace(value)
    if isinstance(value, Mapping):
        for child_key, child in value.items():
            yield from _iter_scalars(child, key=str(child_key))
    elif isinstance(value, (list, tuple, set)):
        for child in value:
            yield from _iter_scalars(child, key=key)
    elif isinstance(value, (str, Path)):
        yield key, value


def _inside(path: Path, root: Path) -> bool:
    return path == root or root in path.parents


def discover_input_files(
    *values: Any,
    excluded_paths: Sequence[str | Path] = (),
    excluded_roots: Sequence[str | Path] = (),
) -> dict[str, Path]:
    """Discover explicit files, directory contents, and reference-prefix files."""
    ignored_paths = {Path(value).expanduser().resolve() for value in excluded_paths}
    ignored_roots = tuple(
        Path(value).expanduser().resolve() for value in excluded_roots
    )
    discovered: set[Path] = set()
    for key, raw_value in _iter_scalars(values):
        if key in {"output_directory", "output_layout"}:
            continue
        text = str(raw_value).strip()
        if not text or text in {".", ".."} or "{" in text or "}" in text:
            continue
        candidate = Path(text).expanduser()
        try:
            resolved = candidate.resolve()
        except OSError:
            continue
        if resolved in ignored_paths or any(
            _inside(resolved, root) for root in ignored_roots
        ):
            continue
        if candidate.is_symlink():
            continue
        if resolved.is_file():
            discovered.add(resolved)
            continue
        if resolved.is_dir():
            for child in resolved.rglob("*"):
                if child.is_file() and not child.is_symlink():
                    discovered.add(child.resolve())
            continue
        if not (candidate.is_absolute() or candidate.parent != Path(".")):
            continue
        parent = resolved.parent
        if parent.is_dir():
            for child in parent.glob(resolved.name + "*"):
                if child.is_file() and not child.is_symlink():
                    discovered.add(child.resolve())
    return {str(path): path for path in sorted(discovered)}


def software_identity(*values: Any) -> dict[str, Any]:
    """Fingerprint PostGWAS, Python, and executable values used by a stage."""
    try:
        package_version = version("postgwas")
    except PackageNotFoundError:
        package_version = "source-tree"
    executables: dict[str, dict[str, Any]] = {}
    for _, raw_value in _iter_scalars(values):
        text = str(raw_value).strip()
        if not text or any(character.isspace() for character in text):
            continue
        resolved = shutil.which(text)
        if not resolved:
            continue
        path = Path(resolved).resolve()
        if path.is_file() and not path.is_symlink():
            executables[str(path)] = {
                "size": path.stat().st_size,
                "sha256": sha256(path),
            }
    return {
        "postgwas": package_version,
        "python": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "executable": str(Path(sys.executable).resolve()),
        "executables": executables,
    }


def _snapshot_files(root: Path) -> dict[Path, tuple[int, int]]:
    if not root.exists():
        return {}
    if root.is_symlink() or not root.is_dir():
        raise RuntimeError("Checkpoint artifact root must be a real directory: %s" % root)
    return {
        path.resolve(): (path.stat().st_size, path.stat().st_mtime_ns)
        for path in root.rglob("*")
        if path.is_file() and not path.is_symlink()
    }


def _declared_files(value: Any, root: Path) -> set[Path]:
    files: set[Path] = set()
    for _, raw_value in _iter_scalars(value):
        candidate = Path(str(raw_value)).expanduser()
        if candidate.is_symlink():
            continue
        resolved = candidate.resolve()
        if not _inside(resolved, root):
            continue
        if resolved.is_file():
            files.add(resolved)
        elif resolved.is_dir():
            files.update(
                child.resolve() for child in resolved.rglob("*")
                if child.is_file() and not child.is_symlink()
            )
    return files


def _owned_fingerprint(path: Path) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise RuntimeError("Checkpoint output must be a regular file: %s" % path)
    return {
        "path": str(path.resolve()),
        "size": path.stat().st_size,
        "sha256": sha256(path),
    }


class ExecutionCheckpoint:
    """Validate, resume, or safely restart one direct command or pipeline stage."""

    def __init__(
        self,
        *,
        manifest_path: str | Path,
        output_root: str | Path,
        artifact_root: str | Path,
        identity: Mapping[str, Any],
        configuration: Any,
        inputs: Mapping[str, str | Path],
        policy: Any,
        resume: bool,
        overwrite: bool,
        logger: CheckpointAuditLogger,
        software: Mapping[str, Any],
        upstream_checkpoints: Mapping[str, str | Path] | None = None,
        error_type: Type[Exception] = RuntimeError,
    ):
        self.manifest_path = Path(manifest_path).expanduser().resolve()
        self.output_root = Path(output_root).expanduser().resolve()
        self.artifact_root = Path(artifact_root).expanduser().resolve()
        self.identity = dict(identity)
        self.configuration_sha256 = configuration_digest(configuration)
        self.input_paths = {
            str(name): Path(value).expanduser().resolve()
            for name, value in inputs.items()
        }
        self.policy = policy
        self.resume = bool(resume)
        self.overwrite = bool(overwrite)
        self.logger = logger
        logger_path = getattr(logger, "path", getattr(logger, "log_path", None))
        if logger_path is None:
            raise TypeError("Checkpoint logger must expose path or log_path")
        self.logger_path = Path(logger_path).expanduser().resolve()
        self.software = dict(software)
        self.upstream_paths = {
            str(name): Path(value).expanduser().resolve()
            for name, value in (upstream_checkpoints or {}).items()
        }
        self.error_type = error_type
        self._manifest_existed_before = self.manifest_path.is_file()
        self._restart_applied = False
        self._before = _snapshot_files(self.artifact_root)

    @property
    def should_record_partial(self) -> bool:
        """Return whether a failed attempt may replace the prior checkpoint."""
        return (
            self.overwrite
            or not self._manifest_existed_before
            or self._restart_applied
        )

    def _input_fingerprints(self) -> dict[str, dict[str, Any]]:
        return {
            name: file_fingerprint(path, error_type=self.error_type)
            for name, path in self.input_paths.items()
        }

    def _upstream_fingerprints(self) -> dict[str, dict[str, Any]]:
        return {
            name: file_fingerprint(path, error_type=self.error_type)
            for name, path in self.upstream_paths.items()
        }

    def _load(self) -> Mapping[str, Any]:
        try:
            document = yaml.safe_load(
                self.manifest_path.read_text(encoding="utf-8")
            ) or {}
        except (OSError, yaml.YAMLError) as exc:
            raise self.error_type(
                "Cannot read execution checkpoint %s: %s"
                % (self.manifest_path, exc)
            ) from exc
        if not isinstance(document, Mapping):
            raise self.error_type(
                "Execution checkpoint root must be a mapping: %s"
                % self.manifest_path
            )
        return document

    def _restart(
        self,
        document: Mapping[str, Any],
        *,
        policy_field: str,
        reason: str,
        warning: str,
    ) -> CheckpointDecision:
        decision = completion_restart_decision(
            document,
            policy=self.policy,
            policy_field=policy_field,
            reason=reason,
            warning=warning,
            error_type=self.error_type,
        )
        outputs = document.get("outputs")
        if isinstance(outputs, Mapping) and outputs:
            apply_completion_restart(
                decision,
                output_root=self.output_root,
                manifest=self.manifest_path,
                logger=self.logger,
                operation=str(self.identity.get("stage", "execution")),
                error_type=self.error_type,
            )
        else:
            self.logger.warning(warning)
            self.logger.record(
                "ACTION",
                str(self.identity.get("stage", "execution")),
                action="restart",
                reason=reason,
                manifest=str(self.manifest_path),
            )
            self.manifest_path.unlink(missing_ok=True)
        self._restart_applied = True
        self._before = _snapshot_files(self.artifact_root)
        return CheckpointDecision(action="run")

    def prepare(self) -> CheckpointDecision:
        """Return a validated reuse decision, applying safe restart when needed."""
        if self.overwrite:
            self.logger.record(
                "ACTION", str(self.identity.get("stage", "execution")),
                action="forced_restart", reason="overwrite_enabled",
            )
            return CheckpointDecision(action="run")
        if not self.resume:
            self.logger.record(
                "ACTION", str(self.identity.get("stage", "execution")),
                action="run", reason="resume_disabled",
            )
            return CheckpointDecision(action="run")
        if not self.manifest_path.is_file():
            return CheckpointDecision(action="run")

        document = self._load()
        expected_identity = {
            "schema_version": EXECUTION_CHECKPOINT_SCHEMA_VERSION,
            **self.identity,
        }
        for field, expected in expected_identity.items():
            if document.get(field) != expected:
                raise self.error_type(
                    "Execution checkpoint %s does not match %s"
                    % (field, expected)
                )
        if document.get("configuration_sha256") != self.configuration_sha256:
            return self._restart(
                document,
                policy_field="changed_parameters",
                reason="changed_parameters",
                warning=(
                    "Resolved parameters changed for %s; PostGWAS will replace "
                    "only checksum-validated owned outputs and restart from this "
                    "safe boundary." % self.identity.get("stage", "the execution")
                ),
            )
        if document.get("software") != self.software:
            return self._restart(
                document,
                policy_field="changed_parameters",
                reason="changed_software",
                warning=(
                    "The software identity changed for %s; PostGWAS will restart "
                    "from this safe boundary."
                    % self.identity.get("stage", "the execution")
                ),
            )

        current_inputs = self._input_fingerprints()
        recorded_inputs = document.get("inputs")
        if isinstance(recorded_inputs, Mapping):
            for name, recorded in recorded_inputs.items():
                recorded_path = (
                    Path(recorded.get("path", "")).expanduser()
                    if isinstance(recorded, Mapping)
                    else Path("")
                )
                if not recorded_path.is_file():
                    raise self.error_type(
                        "Tracked checkpoint input is missing; existing validated "
                        "outputs were preserved: %s" % name
                    )
        if (
            not isinstance(recorded_inputs, Mapping)
            or set(recorded_inputs) != set(current_inputs)
        ):
            return self._restart(
                document,
                policy_field="changed_inputs",
                reason="changed_inputs",
                warning=(
                    "The tracked input set changed for %s; PostGWAS will restart "
                    "from this safe boundary."
                    % self.identity.get("stage", "the execution")
                ),
            )
        for name, fingerprint in current_inputs.items():
            if recorded_inputs.get(name) != fingerprint:
                return self._restart(
                    document,
                    policy_field="changed_inputs",
                    reason="changed_inputs",
                    warning=(
                        "Tracked input %s changed for %s; PostGWAS will restart "
                        "from this safe boundary."
                        % (name, self.identity.get("stage", "the execution"))
                    ),
                )

        current_upstream = self._upstream_fingerprints()
        if document.get("upstream_checkpoints") != current_upstream:
            return self._restart(
                document,
                policy_field="changed_inputs",
                reason="changed_upstream_checkpoint",
                warning=(
                    "An upstream checkpoint changed for %s; PostGWAS will "
                    "invalidate this downstream checkpoint and restart safely."
                    % self.identity.get("stage", "the execution")
                ),
            )

        if document.get("status") != "COMPLETED":
            return self._restart(
                document,
                policy_field="unvalidated_outputs",
                reason="partial_results",
                warning=(
                    "A real partial result exists for %s, but it is not a completed "
                    "scientific checkpoint; PostGWAS will discard only its verified "
                    "owned files and restart this stage."
                    % self.identity.get("stage", "the execution")
                ),
            )
        if document.get("scientific_validation") != "passed_by_stage":
            return self._restart(
                document,
                policy_field="unvalidated_outputs",
                reason="unvalidated_outputs",
                warning=(
                    "Scientific output validation is not recorded for %s; "
                    "PostGWAS will restart from this safe boundary."
                    % self.identity.get("stage", "the execution")
                ),
            )

        outputs = document.get("outputs")
        if not isinstance(outputs, Mapping) or not outputs:
            raise self.error_type("Execution checkpoint contains no owned outputs")
        missing = []
        for name, recorded in outputs.items():
            if not isinstance(recorded, Mapping) or not isinstance(recorded.get("path"), str):
                raise self.error_type(
                    "Execution checkpoint output %s fingerprint is invalid" % name
                )
            raw_path = Path(recorded["path"]).expanduser()
            path = raw_path.resolve()
            if raw_path.is_symlink() or not _inside(path, self.output_root):
                raise self.error_type(
                    "Execution checkpoint output is unsafe: %s" % path
                )
            if not path.is_file():
                missing.append(str(name))
                continue
            current = _owned_fingerprint(path)
            if dict(recorded) != current:
                raise self.error_type(
                    "Refusing automatic restart because checkpoint output %s was "
                    "modified after validation: %s" % (name, path)
                )
        if missing:
            return self._restart(
                document,
                policy_field="unvalidated_outputs",
                reason="incomplete_outputs",
                warning=(
                    "Checkpoint outputs are incomplete for %s (missing: %s); "
                    "PostGWAS will restart from this safe boundary."
                    % (self.identity.get("stage", "the execution"), ", ".join(missing))
                ),
            )
        self.logger.record(
            "RESUME", str(self.identity.get("stage", "execution")),
            manifest=str(self.manifest_path),
        )
        return CheckpointDecision(action="resume", document=document)

    def _output_fingerprints(
        self,
        *,
        declared_values: Any,
        include_modified: bool,
    ) -> dict[str, dict[str, Any]]:
        after = _snapshot_files(self.artifact_root)
        owned = {
            path for path, signature in after.items()
            if path not in self._before
            or (include_modified and self._before.get(path) != signature)
        }
        owned.update(_declared_files(declared_values, self.artifact_root))
        owned.discard(self.manifest_path)
        owned.discard(self.logger_path)
        return {
            str(path.relative_to(self.output_root)): _owned_fingerprint(path)
            for path in sorted(owned)
            if _inside(path, self.output_root)
        }

    def write(
        self,
        *,
        status: Literal["COMPLETED", "PARTIAL"],
        declared_values: Any = None,
        state: Any = None,
        metrics: Mapping[str, Any] | None = None,
    ) -> Path:
        """Atomically record one completed or isolated partial execution."""
        outputs = self._output_fingerprints(
            declared_values=declared_values,
            include_modified=self.overwrite,
        )
        if status == "COMPLETED" and not outputs:
            raise NoOwnedCheckpointOutputs(
                "Cannot create a reusable checkpoint because no owned output "
                "files were identified below %s" % self.artifact_root
            )
        document = {
            "schema_version": EXECUTION_CHECKPOINT_SCHEMA_VERSION,
            **self.identity,
            "status": status,
            "scientific_validation": (
                "passed_by_stage" if status == "COMPLETED" else "incomplete"
            ),
            "configuration_sha256": self.configuration_sha256,
            "inputs": self._input_fingerprints(),
            "software": self.software,
            "upstream_checkpoints": self._upstream_fingerprints(),
            "outputs": outputs,
            "state": encode_checkpoint_value(state),
            "metrics": encode_checkpoint_value(dict(metrics or {})),
        }
        return write_yaml_report(document, self.manifest_path)


__all__ = [
    "CheckpointAuditLogger",
    "CheckpointDecision",
    "ExecutionCheckpoint",
    "NoOwnedCheckpointOutputs",
    "checkpoint_content_arguments",
    "checkpoint_content_configuration",
    "checkpoint_content_namespace",
    "checkpoint_namespace",
    "decode_checkpoint_value",
    "discover_input_files",
    "encode_checkpoint_value",
    "preserve_checkpoint_run_controls",
    "restore_checkpoint_namespace",
    "software_identity",
]
