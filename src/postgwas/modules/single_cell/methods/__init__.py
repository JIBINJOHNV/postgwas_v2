"""Independently registered GWAS/single-cell integration methods."""

from postgwas.modules.single_cell.methods.registry import (
    METHOD_REGISTRY,
    get_single_cell_methods,
)

__all__ = ["METHOD_REGISTRY", "get_single_cell_methods"]
