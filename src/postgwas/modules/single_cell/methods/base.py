"""Common contracts for independently implemented single-cell methods."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Protocol

from postgwas.core.contracts import ModuleResult


@dataclass(frozen=True)
class MethodRunContext:
    """Resolved shared state supplied to one selected single-cell method."""

    args: argparse.Namespace
    configuration: Any
    output: Path
    dataset: str
    paths: Mapping[str, Path]
    logger: Any


class SingleCellMethod(Protocol):
    """Interface implemented by every registered single-cell method."""

    name: str
    resource_paths: tuple[str, ...]

    def prepare_pipeline_args(self, args: argparse.Namespace) -> None:
        """Apply method-owned pipeline bindings before configuration resolution."""

    def preflight_pipeline(self, args: argparse.Namespace, configuration) -> object:
        """Return evidence for resources available before upstream stages run."""

    def preflight_direct(
        self,
        args: argparse.Namespace,
        configuration,
        dataset: str,
    ) -> object:
        """Resolve and validate the exact inputs used by this method."""

    def run(self, context: MethodRunContext, preflight: object) -> ModuleResult:
        """Execute one preflighted method and publish validated artifacts."""


__all__ = ["MethodRunContext", "SingleCellMethod"]
