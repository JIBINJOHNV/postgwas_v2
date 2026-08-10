"""User-facing errors raised by the PostGWAS FLAMES integration."""


class FlamesError(RuntimeError):
    """FLAMES configuration, scientific validation, or execution failed."""


__all__ = ["FlamesError"]
