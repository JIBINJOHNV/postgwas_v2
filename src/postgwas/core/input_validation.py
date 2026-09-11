"""Current-invocation evidence and reuse for shared input validators.

Validators own their scientific contracts; this module does not infer a format,
genome build, value range, or completeness policy. A successful availability
check is never promoted to successful content or compatibility validation.
Session evidence is transient and must not be restored from checkpoints.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, replace
import json
import math
from pathlib import Path
from typing import Any, Iterator, TypeVar


_Result = TypeVar("_Result")
# These are report protocol states, not configurable scientific policies.
_STATUSES = frozenset({"passed", "failed", "warning", "deferred", "blocked"})
_SESSION: ContextVar[InputValidationSession | None] = ContextVar(
    "postgwas_input_validation_session", default=None,
)
_RECORDER: ContextVar[InputValidationSession | None] = ContextVar(
    "postgwas_input_validation_recorder", default=None,
)
_CONSUMER: ContextVar[str | None] = ContextVar(
    "postgwas_input_validation_consumer", default=None,
)


def _json_value(value: Any) -> Any:
    """Copy compact evidence without stringifying dataframes or arbitrary objects."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("Input-validation evidence must contain finite numbers")
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        if not all(isinstance(key, str) for key in value):
            raise TypeError("Input-validation evidence keys must be strings")
        return {key: _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    raise TypeError(
        "Input-validation evidence requires JSON scalars or containers, not %s"
        % type(value).__name__
    )


@dataclass(frozen=True)
class PreflightFileIdentity:
    """Cheap file identity for one invocation, not a durable content fingerprint."""

    path: Path
    size_bytes: int
    modified_ns: int
    device: int | None = None
    inode: int | None = None
    changed_ns: int | None = None


def capture_preflight_file_identities(
    paths,
    *,
    error_type=RuntimeError,
    label: str = "Validated resource",
) -> tuple[PreflightFileIdentity, ...]:
    """Capture current identities using the existing preflight compatibility API."""
    identities = []
    observed = set()
    for value in paths:
        path = Path(value).expanduser().resolve()
        if path in observed:
            continue
        observed.add(path)
        try:
            metadata = path.stat()
        except OSError as exc:
            message = "%s is unavailable: %s" % (label, path)
            record_file_validation(
                path, label, checks=("regular file", "non-empty file"),
                status="failed", message=message,
            )
            raise error_type(message) from exc
        if not path.is_file() or metadata.st_size <= 0:
            message = "%s is missing or empty: %s" % (label, path)
            record_file_validation(
                path, label, checks=("regular file", "non-empty file"),
                status="failed", message=message,
            )
            raise error_type(message)
        record_file_validation(
            path, label, checks=("regular file", "non-empty file"),
        )
        identities.append(PreflightFileIdentity(
            path=path,
            size_bytes=int(metadata.st_size),
            modified_ns=int(metadata.st_mtime_ns),
            device=int(metadata.st_dev),
            inode=int(metadata.st_ino),
            changed_ns=int(metadata.st_ctime_ns),
        ))
    return tuple(identities)


def require_unchanged_preflight_files(
    identities,
    *,
    error_type=RuntimeError,
    label: str = "Validated resource",
) -> tuple[Path, ...]:
    """Reject resources replaced between early preflight and their consumer."""
    expected = tuple(identities)
    current = capture_preflight_file_identities(
        (identity.path for identity in expected),
        error_type=error_type,
        label=label,
    )
    if current != expected:
        message = (
            "%s changed after pipeline preflight; restart so validation and "
            "execution use the same files." % label
        )
        for before, after in zip(expected, current):
            if before != after:
                record_file_validation(
                    before.path, label, status="failed",
                    checks=("unchanged since validation",), message=message,
                )
        raise error_type(message)
    return tuple(identity.path for identity in current)


@dataclass(frozen=True)
class FileValidationRecord:
    """One specifically named check outcome and the consumers requiring it."""

    path: str | None
    role: str
    checks: tuple[str, ...]
    status: str
    metrics: Mapping[str, Any]
    message: str
    consumers: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "role": self.role,
            "checks": list(self.checks),
            "status": self.status,
            "metrics": _json_value(self.metrics),
            "message": self.message,
            "consumers": list(self.consumers),
        }


@dataclass(frozen=True)
class _CachedValidation:
    result: Any
    records: tuple[int, ...]
    dependencies: tuple[PreflightFileIdentity, ...]


class InputValidationSession:
    """Collect checks and optionally reuse exact contracts during one invocation.

    Nested module scopes are supported. Context-local state avoids attributing
    checks to another invocation; workers must explicitly return their evidence.
    Cached operations must be read-only and their contract must include the
    validator identity and every option that influences the result. Their file
    dependencies include those observed by nested cached validators.
    ``pipeline=False`` permits exact-contract reuse in direct commands without
    enabling pipeline-specific runtime branches. ``reuse_checks=False`` is an
    evidence-only session. Nothing is reused between invocations.
    """

    def __init__(self, *, reuse_checks: bool = True, pipeline: bool = True) -> None:
        self._reuse_checks = reuse_checks
        self._pipeline = pipeline
        self._records: list[FileValidationRecord] = []
        self._record_indexes: dict[str, int] = {}
        self._cache: dict[tuple[Any, ...], _CachedValidation] = {}
        self._identities: dict[Path, PreflightFileIdentity] = {}
        self._observations: list[set[int]] = []
        self._dependency_observations: list[dict[Path, PreflightFileIdentity]] = []
        self._token = None
        self._recorder_token = None
        self._consumer_token = None
        self._finished = False

    def __enter__(self) -> InputValidationSession:
        if self._recorder_token is not None:
            raise RuntimeError("An input-validation session is already active")
        if self._finished:
            raise RuntimeError(
                "Create a fresh input-validation session for each invocation"
            )
        self._recorder_token = _RECORDER.set(self)
        self._token = _SESSION.set(self if self._reuse_checks and self._pipeline else None)
        self._consumer_token = _CONSUMER.set(None)
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        _CONSUMER.reset(self._consumer_token)
        if self._token is not None:
            _SESSION.reset(self._token)
        _RECORDER.reset(self._recorder_token)
        self._recorder_token = None
        self._consumer_token = None
        self._token = None
        self._finished = True

    @contextmanager
    def scope(self, module: str) -> Iterator[None]:
        """Associate newly performed and reused checks with one module."""
        if _RECORDER.get() is not self:
            raise RuntimeError("Input-validation scope requires its active session")
        name = str(module).strip()
        if not name:
            raise ValueError("Input-validation scope requires a module name")
        token = _CONSUMER.set(name)
        try:
            yield
        finally:
            _CONSUMER.reset(token)

    @property
    def records(self) -> tuple[FileValidationRecord, ...]:
        """Return detached report records so consumers cannot mutate evidence."""
        return tuple(
            replace(record, metrics=_json_value(record.metrics))
            for record in self._records
        )

    def to_dict(self) -> dict[str, Any]:
        return {"files": [record.to_dict() for record in self._records]}

    def _observe(self, index: int) -> None:
        record = self._records[index]
        consumer = _CONSUMER.get()
        if consumer is not None and consumer not in record.consumers:
            self._records[index] = replace(
                record, consumers=(*record.consumers, consumer),
            )
        for observations in self._observations:
            observations.add(index)

    def _record(self, record: FileValidationRecord) -> None:
        key = json.dumps(record.to_dict(), sort_keys=True, allow_nan=False)
        index = self._record_indexes.get(key)
        if index is None:
            index = len(self._records)
            self._record_indexes[key] = index
            self._records.append(record)
        self._observe(index)

    def _validate_once(
        self,
        paths: Iterable[str | Path],
        contract: Any,
        operation: Callable[[], _Result],
        *,
        error_type=RuntimeError,
    ) -> _Result:
        # The contract is explicit: closure identity cannot prove equal policies.
        encoded_contract = json.dumps(
            _json_value(contract), sort_keys=True, allow_nan=False,
        )
        identities = capture_preflight_file_identities(paths, error_type=error_type)
        if not identities:
            # Without an input identity there is no evidence for safe file reuse.
            return operation()
        for identity in identities:
            previous = self._identities.get(identity.path)
            if previous is not None and previous != identity and self._pipeline:
                message = (
                    "Validated resource changed after pipeline preflight; "
                    "restart so validation and execution use the same file: %s"
                    % identity.path
                )
                record_file_validation(
                    identity.path, "Input resource", status="failed",
                    checks=("unchanged since validation",), message=message,
                )
                raise error_type(message)
            self._identities[identity.path] = identity
        dependencies = {identity.path: identity for identity in identities}
        for observed in self._dependency_observations:
            observed.update(dependencies)
        key = (identities, encoded_contract)
        cached = self._cache.get(key)
        if cached is not None and not self._pipeline:
            # A direct module may publish a newly generated file at a path it
            # checked earlier. Revalidate the new version, including nested
            # dependencies; never replay a check of the previous contents.
            try:
                current = capture_preflight_file_identities(
                    (identity.path for identity in cached.dependencies), error_type=error_type,
                )
            except Exception:
                # An absent/unreadable dependency must reach the direct
                # validator's own error policy, not fail inside cache lookup.
                current = None
            if current != cached.dependencies:
                cached = None
        if cached is not None:
            # A parent check can depend on nested validators of other files.
            # Replaying their evidence requires validating those identities too.
            require_unchanged_preflight_files(cached.dependencies, error_type=error_type)
            for observed in self._dependency_observations:
                observed.update((identity.path, identity)
                                for identity in cached.dependencies)
            for index in cached.records:
                self._observe(index)
            return cached.result

        observations: set[int] = set()
        self._observations.append(observations)
        self._dependency_observations.append(dependencies)
        try:
            result = operation()
            require_unchanged_preflight_files(dependencies.values(), error_type=error_type)
        except Exception as exc:
            record_file_validation(
                None, "Input resource validation", status="failed",
                metrics={"paths": [str(path) for path in dependencies]},
                message=str(exc),
            )
            raise
        finally:
            self._dependency_observations.pop()
            self._observations.pop()
        self._cache[key] = _CachedValidation(
            result, tuple(sorted(observations)), tuple(dependencies.values()),
        )
        return result


def current_validation_session() -> InputValidationSession | None:
    """Return current transient evidence, never a restored checkpoint."""
    return _SESSION.get()


def current_validation_recorder() -> InputValidationSession | None:
    """Return the evidence owner without implying pipeline execution or reuse."""
    return _RECORDER.get()


@contextmanager
def validation_scope(module: str) -> Iterator[None]:
    """Attribute recorded checks without changing a consumer's validation policy."""
    session = _RECORDER.get()
    if session is None:
        yield
    else:
        with session.scope(module):
            yield


def record_file_validation(
    path: str | Path | None,
    role: str,
    *,
    checks: Iterable[str] = (),
    metrics: Mapping[str, Any] | None = None,
    status: str = "passed",
    message: str = "",
) -> None:
    """Record the checks actually performed; do nothing outside a session."""
    session = _RECORDER.get()
    if session is None:
        return
    if status not in _STATUSES:
        raise ValueError("Unknown input-validation status: %s" % status)
    session._record(FileValidationRecord(
        path=None if path is None else str(Path(path).expanduser().resolve()),
        role=str(role),
        checks=tuple(str(check) for check in checks),
        status=status,
        metrics=_json_value({} if metrics is None else metrics),
        message=str(message),
        consumers=(),
    ))


def validate_once(
    paths: Iterable[str | Path],
    contract: Any,
    operation: Callable[[], _Result],
    *,
    error_type=RuntimeError,
) -> _Result:
    """Reuse a read-only check only for unchanged files and the exact contract."""
    session = _RECORDER.get()
    paths = tuple(paths)
    cacheable = True
    if session is not None and not session._pipeline:
        # Direct callers own their missing/empty-file policy and error type.
        # When identity capture cannot support reuse, execute that original
        # validator rather than imposing a stronger pipeline contract.
        try:
            cacheable = all(Path(path).is_file() and Path(path).stat().st_size > 0 for path in paths)
        except OSError:
            cacheable = False
    if session is None or not session._reuse_checks or not cacheable:
        # Evidence-only callers retain their original execution semantics.
        if _RECORDER.get() is None:
            return operation()
        dependencies = tuple(paths)
        try:
            return operation()
        except Exception as exc:
            record_file_validation(
                None, "Input resource validation", status="failed",
                metrics={"paths": [str(Path(path).expanduser().resolve()) for path in dependencies]},
                message=str(exc),
            )
            raise
    return session._validate_once(paths, contract, operation, error_type=error_type)


__all__ = [
    "FileValidationRecord",
    "InputValidationSession",
    "PreflightFileIdentity",
    "capture_preflight_file_identities",
    "current_validation_session",
    "current_validation_recorder",
    "record_file_validation",
    "require_unchanged_preflight_files",
    "validate_once",
    "validation_scope",
]
