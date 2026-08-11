"""Canonical registry of implemented single-cell method adapters."""

from __future__ import annotations

from types import MappingProxyType
from typing import Iterable, Mapping

from postgwas.config.models.modules.single_cell import SINGLE_CELL_TOOLS
from postgwas.modules.single_cell.errors import SingleCellError
from postgwas.modules.single_cell.methods.base import SingleCellMethod
from postgwas.modules.single_cell.methods.ldsc_celltype.method import (
    METHOD as LDSC_CELLTYPE_METHOD,
)
from postgwas.modules.single_cell.methods.magma_celltype.method import (
    METHOD as MAGMA_CELLTYPE_METHOD,
)
from postgwas.modules.single_cell.methods.scdrs.method import METHOD as SCDRS_METHOD


_METHODS = (MAGMA_CELLTYPE_METHOD, SCDRS_METHOD, LDSC_CELLTYPE_METHOD)
METHOD_REGISTRY: Mapping[str, SingleCellMethod] = MappingProxyType({
    method.name: method for method in _METHODS
})

if tuple(METHOD_REGISTRY) != tuple(SINGLE_CELL_TOOLS):
    raise RuntimeError(
        "The single-cell method registry must exactly match SingleCellTool"
    )


def get_single_cell_methods(names: Iterable[str]) -> tuple[SingleCellMethod, ...]:
    """Resolve configured method names in user-selected execution order."""
    methods = []
    for name in names:
        try:
            methods.append(METHOD_REGISTRY[name])
        except KeyError as exc:
            raise SingleCellError(
                "Unsupported configured single-cell tool: %s" % name
            ) from exc
    return tuple(methods)


__all__ = ["METHOD_REGISTRY", "get_single_cell_methods"]
