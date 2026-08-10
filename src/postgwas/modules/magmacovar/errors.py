"""Domain-specific failures for MAGMA gene-property analysis."""


class MagmaCovarError(RuntimeError):
    """MAGMAcovar configuration, input, execution, or output is invalid."""


__all__ = ["MagmaCovarError"]
