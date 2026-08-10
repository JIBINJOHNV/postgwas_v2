"""Domain-specific failures for single-cell integration."""


class SingleCellError(RuntimeError):
    """Single-cell configuration, input, analysis, or output is invalid."""


__all__ = ["SingleCellError"]
