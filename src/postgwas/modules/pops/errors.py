"""User-facing errors raised by the PostGWAS PoPS integration."""


class PopsError(RuntimeError):
    """PoPS configuration, input validation, or execution failed."""


__all__ = ["PopsError"]
