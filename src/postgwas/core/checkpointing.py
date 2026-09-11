"""Validated execution checkpoints shared by direct commands and pipelines."""

from __future__ import annotations

from argparse import Namespace
from dataclasses import dataclass
from datetime import datetime, timezone
from enum import Enum
import errno
from importlib.metadata import PackageNotFoundError, version
from importlib.machinery import BYTECODE_SUFFIXES
import json
from pathlib import Path
import platform
import shutil
import stat
import sys
from typing import Any, Iterator, Literal, Mapping, Sequence, Type

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
from postgwas.core.ui import print_screen_message


EXECUTION_CHECKPOINT_SCHEMA_VERSION = 1
_TYPE_KEY = "__postgwas_type__"
_INPUT_FINGERPRINT_CACHE_KEY = "input_fingerprint_cache"
_RUN_CONTROL_DESTINATIONS = frozenset({"resume", "overwrite"})
_RUN_CONTROL_OPTIONS = frozenset({"--resume", "--overwrite"})


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
        print_screen_message("warning", "WARNING: %s" % message, console=self.console)

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


def _iter_scalars(
    value: Any, *, key: str = "", excluded_key_suffixes: tuple[str, ...] = (),
):
    if any(key == name or key.endswith("_" + name) for name in excluded_key_suffixes):
        return
    if hasattr(value, "model_dump"):
        value = value.model_dump(mode="python")
    if isinstance(value, Namespace):
        value = checkpoint_namespace(value)
    if isinstance(value, Mapping):
        for child_key, child in value.items():
            yield from _iter_scalars(
                child, key=str(child_key), excluded_key_suffixes=excluded_key_suffixes,
            )
    elif isinstance(value, (list, tuple, set)):
        for child in value:
            yield from _iter_scalars(child, key=key, excluded_key_suffixes=excluded_key_suffixes)
    elif isinstance(value, (str, Path)):
        yield key, value


def _inside(path: Path, root: Path) -> bool:
    return path == root or root in path.parents


def _optional_input_signature(path: str | Path) -> dict[str, Any] | None:
    """Return the filesystem identity used to validate a cached content hash."""
    resolved = Path(path).expanduser().resolve()
    try:
        metadata = resolved.stat()
    except OSError:
        return None
    if not stat.S_ISREG(metadata.st_mode) or metadata.st_size <= 0:
        return None
    return {
        "path": str(resolved),
        "device": metadata.st_dev,
        "inode": metadata.st_ino,
        "size": metadata.st_size,
        "mtime_ns": metadata.st_mtime_ns,
        "ctime_ns": metadata.st_ctime_ns,
    }


def _required_input_signature(
    path: str | Path,
    *,
    error_type: Type[Exception],
) -> dict[str, Any]:
    signature = _optional_input_signature(path)
    if signature is None:
        raise error_type(
            "Completion file is missing or empty: %s"
            % Path(path).expanduser().resolve()
        )
    return signature


def _valid_sha256(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _cache_matches_signature(
    entry: Any,
    signature: Mapping[str, Any],
) -> bool:
    return (
        isinstance(entry, Mapping)
        and _valid_sha256(entry.get("sha256"))
        and all(entry.get(field) == value for field, value in signature.items())
    )


def _scientific_fingerprint(entry: Mapping[str, Any]) -> dict[str, Any]:
    """Keep checkpoint comparison based on canonical path, size, and SHA-256."""
    return {
        "path": entry["path"],
        "size": entry["size"],
        "sha256": entry["sha256"],
    }


def discover_input_files(
    *values: Any,
    excluded_paths: Sequence[str | Path] = (),
    excluded_roots: Sequence[str | Path] = (),
) -> dict[str, Path]:
    """Discover explicit files, directory contents, and reference-prefix files.

    An excluded directory is never entered while expanding a parent directory.
    A file below it remains discoverable when that file is supplied explicitly;
    this preserves legitimate upstream inputs located below an output root while
    preventing accidental categorical-value matches from tracking live outputs.
    """
    ignored_paths = {Path(value).expanduser().resolve() for value in excluded_paths}
    ignored_roots = tuple(
        Path(value).expanduser().resolve() for value in excluded_roots
    )
    ignored_expansion_roots = tuple(
        path for path in ignored_paths if path.is_dir()
    )

    def excluded(path: Path) -> bool:
        return path in ignored_paths or any(
            _inside(path, root) for root in ignored_roots
        )

    def excluded_from_expansion(path: Path, source: Path) -> bool:
        return excluded(path) or any(
            _inside(path, root) and not _inside(source, root)
            for root in ignored_expansion_roots
        )

    def expanded_files(source: Path) -> Iterator[Path]:
        pending = [source]
        while pending:
            directory = pending.pop()
            for child in directory.iterdir():
                if child.is_symlink():
                    continue
                resolved_child = child.resolve()
                if excluded_from_expansion(resolved_child, source):
                    continue
                if child.is_dir():
                    pending.append(child)
                elif child.is_file():
                    yield resolved_child

    discovered: set[Path] = set()
    # Output-layout mappings contain relative names, not input resources. Prune
    # the entire subtree before its parent key is lost, including namespaced
    # resolved copies in a pipeline Namespace. Explicit upstream artifact paths
    # outside these settings remain inputs and retain checksum validation.
    for _, raw_value in _iter_scalars(
        values, excluded_key_suffixes=("output_directory", "output_layout"),
    ):
        text = str(raw_value).strip()
        if not text or text in {".", ".."} or "{" in text or "}" in text:
            continue
        candidate = Path(text).expanduser()
        try:
            resolved = candidate.resolve()
            if excluded(resolved) or candidate.is_symlink():
                continue
            if resolved.is_file():
                discovered.add(resolved)
                continue
            if resolved.is_dir():
                discovered.update(expanded_files(resolved))
                continue
            if not (candidate.is_absolute() or candidate.parent != Path(".")):
                continue
            parent = resolved.parent
            if parent.is_dir():
                for child in parent.glob(resolved.name + "*"):
                    resolved_child = child.resolve()
                    if (
                        not excluded_from_expansion(resolved_child, resolved)
                        and child.is_file()
                        and not child.is_symlink()
                    ):
                        discovered.add(resolved_child)
        except OSError as exc:
            # Configuration and runtime state also contain ordinary strings.
            # A value that the operating system cannot even probe as a path
            # (for example the long LS_COLORS shell setting) is not an input.
            if exc.errno != errno.ENAMETOOLONG:
                raise
            continue
    return {str(path): path for path in sorted(discovered)}


def _package_content_identity(root: Path) -> dict[str, Any]:
    """Identify installed source and resources, including editable installs.

    Python's regenerable bytecode cache is excluded by protocol; all other
    package files (including R scripts and bundled resources) determine identity.
    Relative names make this independent of installation-directory spelling.
    Do not memoize: developers can change source without changing its version.
    """
    fingerprints = {
        path.relative_to(root).as_posix(): sha256(path)
        for path in sorted(root.rglob("*"))
        if path.is_file()
        and "__pycache__" not in path.relative_to(root).parts
        and path.suffix not in BYTECODE_SUFFIXES
    }
    if not fingerprints:
        raise RuntimeError("Cannot fingerprint PostGWAS package content: %s" % root)
    return {"files": len(fingerprints), "sha256": configuration_digest(fingerprints)}


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
        try:
            resolved = shutil.which(text)
            if not resolved:
                continue
            path = Path(resolved).resolve()
            if path.is_file() and not path.is_symlink():
                executables[str(path)] = {
                    "size": path.stat().st_size,
                    "sha256": sha256(path),
                }
        except OSError as exc:
            if exc.errno != errno.ENAMETOOLONG:
                raise
            continue
    return {
        "postgwas": package_version,
        "postgwas_content": _package_content_identity(Path(__file__).resolve().parents[1]),
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
        try:
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
        except OSError as exc:
            if exc.errno != errno.ENAMETOOLONG:
                raise
            continue
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
        self._input_start_signatures = {
            name: _optional_input_signature(path)
            for name, path in self.input_paths.items()
        }
        self._input_fingerprint_cache: dict[str, dict[str, Any]] = {}
        self._before = _snapshot_files(self.artifact_root)

    @property
    def should_record_partial(self) -> bool:
        """Return whether a failed attempt may replace the prior checkpoint."""
        return (
            self.overwrite
            or not self._manifest_existed_before
            or self._restart_applied
        )

    def _hydrate_input_fingerprint_cache(
        self,
        document: Mapping[str, Any],
    ) -> None:
        """Reuse trusted prior hashes only while the complete file identity matches."""
        stored_cache = document.get(_INPUT_FINGERPRINT_CACHE_KEY)
        recorded_inputs = document.get("inputs")
        try:
            manifest_mtime_ns = self.manifest_path.stat().st_mtime_ns
        except OSError:
            manifest_mtime_ns = -1
        reused = 0
        migrated = 0
        for name in self.input_paths:
            signature = self._input_start_signatures.get(name)
            if signature is None:
                continue
            entry = (
                stored_cache.get(name)
                if isinstance(stored_cache, Mapping)
                else None
            )
            if _cache_matches_signature(entry, signature):
                self._input_fingerprint_cache[name] = dict(entry)
                reused += 1
                continue

            # Older checkpoints contain a SHA-256 calculated immediately before
            # the atomic manifest write. Migrate it only when both filesystem
            # change clocks predate that manifest and path and size still match.
            recorded = (
                recorded_inputs.get(name)
                if isinstance(recorded_inputs, Mapping)
                else None
            )
            if (
                isinstance(recorded, Mapping)
                and recorded.get("path") == signature["path"]
                and recorded.get("size") == signature["size"]
                and _valid_sha256(recorded.get("sha256"))
                and signature["mtime_ns"] <= manifest_mtime_ns
                and signature["ctime_ns"] <= manifest_mtime_ns
            ):
                self._input_fingerprint_cache[name] = {
                    **signature,
                    "sha256": recorded["sha256"],
                }
                migrated += 1
        if reused or migrated:
            self.logger.record(
                "CACHE",
                str(self.identity.get("stage", "execution")),
                action="hydrate_input_fingerprints",
                reused=reused,
                migrated=migrated,
            )

    def _input_fingerprints(self) -> dict[str, dict[str, Any]]:
        fingerprints: dict[str, dict[str, Any]] = {}
        reused = 0
        calculated = 0
        for name, path in self.input_paths.items():
            signature = _required_input_signature(
                path, error_type=self.error_type,
            )
            if signature != self._input_start_signatures.get(name):
                raise self.error_type(
                    "Tracked input changed while the execution was active: %s"
                    % signature["path"]
                )
            cached = self._input_fingerprint_cache.get(name)
            if _cache_matches_signature(cached, signature):
                fingerprints[name] = _scientific_fingerprint(cached)
                reused += 1
                continue

            before = signature
            fingerprint = file_fingerprint(path, error_type=self.error_type)
            after = _required_input_signature(path, error_type=self.error_type)
            if after != before:
                raise self.error_type(
                    "Tracked input changed while its SHA-256 was calculated: %s"
                    % after["path"]
                )
            cached = {**after, "sha256": fingerprint["sha256"]}
            self._input_fingerprint_cache[name] = cached
            fingerprints[name] = _scientific_fingerprint(cached)
            calculated += 1
        if reused or calculated:
            self.logger.record(
                "CHECKSUM",
                str(self.identity.get("stage", "execution")),
                action="fingerprint_inputs",
                reused=reused,
                calculated=calculated,
            )
        return fingerprints

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
        document = None
        if self.manifest_path.is_file():
            try:
                document = self._load()
            except self.error_type:
                if self.resume and not self.overwrite:
                    raise
            if document is not None:
                self._hydrate_input_fingerprint_cache(document)
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

        if document is None:
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
            _INPUT_FINGERPRINT_CACHE_KEY: {
                name: self._input_fingerprint_cache[name]
                for name in sorted(self._input_fingerprint_cache)
            },
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
