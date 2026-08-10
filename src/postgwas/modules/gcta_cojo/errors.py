"""Module-specific failures for GCTA-COJO."""


class GctaCojoError(RuntimeError):
    """A GCTA-COJO input, scientific-validation, or execution failure."""


__all__ = ["GctaCojoError"]
