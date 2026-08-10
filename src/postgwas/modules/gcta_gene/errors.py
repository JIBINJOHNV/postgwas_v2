"""Errors raised by the GCTA gene-association module."""


class GctaGeneError(RuntimeError):
    """A fastBAT or mBAT-combo run cannot proceed scientifically or technically."""


__all__ = ["GctaGeneError"]
