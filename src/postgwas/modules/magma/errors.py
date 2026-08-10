"""Public MAGMA module exception."""


class MagmaError(RuntimeError):
    """MAGMA validation or execution could not complete safely."""


__all__ = ["MagmaError"]
