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
    """Run-scoped result store shared by pipeline modules."""

    def __init__(self, initial: Mapping[str, Any] | None = None):
        self._values = dict(initial or {})

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

    def snapshot(self) -> dict[str, Any]:
        return dict(self._values)
