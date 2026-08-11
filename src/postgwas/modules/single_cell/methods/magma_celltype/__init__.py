"""MAGMA cell-type enrichment method."""

from postgwas.modules.single_cell.methods.magma_celltype.analysis import (
    normalize_magma_celltype_results,
    validate_magma_celltype_covariates,
    validate_magma_celltype_use_case,
)
from postgwas.modules.single_cell.methods.magma_celltype.method import METHOD

__all__ = [
    "METHOD",
    "normalize_magma_celltype_results",
    "validate_magma_celltype_covariates",
    "validate_magma_celltype_use_case",
]
