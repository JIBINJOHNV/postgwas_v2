"""Typed errors shared without importing the scientific analysis stack."""


class ConcordanceValidationError(RuntimeError):
    """The accuracy audit could not be completed reliably."""
