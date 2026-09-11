"""Small, dependency-free contracts shared by all PostGWAS execution modes."""

from __future__ import annotations

from collections.abc import Iterator, MutableMapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping


@dataclass(frozen=True)
class Artifact:
    """One file or directory produced by a module.

    Metadata carries scientific compatibility information such as genome build,
    ancestry, sample identifier, schema version, and effect-allele convention.
    """

    kind: str
    path: Path
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        object.__setattr__(self, "path", Path(self.path))


@dataclass(frozen=True)
class ModuleResult:
    """Uniform result returned by a scientific module."""

    module: str
    artifacts: Mapping[str, Artifact] = field(default_factory=dict)
    metrics: Mapping[str, Any] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()


class RunContext(MutableMapping[str, Any]):
    """Run-scoped result and transient validation store shared by modules.

    Validation evidence is deliberately excluded from ``snapshot()`` because
    pipeline preflight is rerun for every invocation, including resume. This
    permits module-native evidence objects and prevents stale validation from
    being restored from a checkpoint.
    """

    def __init__(
        self,
        initial: Mapping[str, Any] | None = None,
        *,
        validations: Mapping[str, Any] | None = None,
    ):
        self._values = dict(initial or {})
        self._validations = dict(validations or {})

    def __getitem__(self, key: str) -> Any:
        return self._values[key]

    def __setitem__(self, key: str, value: Any) -> None:
        self._values[key] = value

    def __delitem__(self, key: str) -> None:
        del self._values[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._values)

    def __len__(self) -> int:
        return len(self._values)

    def publish(self, result: ModuleResult) -> None:
        self._values[result.module] = result

    def publish_validation(self, module: str, evidence: Any) -> None:
        """Retain current-invocation preflight evidence for one module."""
        name = str(module).strip()
        if not name:
            raise ValueError("Validation evidence requires a module name")
        if evidence is None:
            raise ValueError("Validation evidence must not be None")
        self._validations[name] = evidence

    def validation(self, module: str, default: Any = None) -> Any:
        """Return transient preflight evidence without checkpointing it."""
        return self._validations.get(str(module), default)

    def validation_modules(self) -> tuple[str, ...]:
        """Return modules with reusable evidence in stable insertion order."""
        return tuple(self._validations)

    def snapshot(self) -> dict[str, Any]:
        return dict(self._values)
