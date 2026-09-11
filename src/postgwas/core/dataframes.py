"""Small, reusable Polars operations with explicit null semantics."""

from __future__ import annotations

import inspect
from typing import Any, Callable, Type

import polars as pl


FLOAT_DTYPES = (pl.Float32, pl.Float64)


def streaming_collect_options(
    collect_method: Callable[..., Any],
    *,
    error_type: Type[Exception] = RuntimeError,
) -> tuple[dict[str, Any], str]:
    """Select an explicitly supported Polars streaming collection API."""
    parameters = inspect.signature(collect_method).parameters
    if "engine" in parameters:
        return {"engine": "streaming"}, "engine=streaming"
    if "streaming" in parameters:
        return {"streaming": True}, "streaming=True"
    raise error_type(
        "Installed Polars %s does not expose a supported streaming LazyFrame "
        "collection option. Install a supported PostGWAS analysis environment."
        % getattr(pl, "__version__", "unknown")
    )


def collect_streaming(
    frame: pl.LazyFrame,
    *,
    error_type: Type[Exception] = RuntimeError,
) -> pl.DataFrame:
    """Collect a lazy frame using a declared Polars streaming option."""
    options, _ = streaming_collect_options(
        frame.collect,
        error_type=error_type,
    )
    return frame.collect(**options)


def numeric_column(frame: pl.DataFrame, column: str) -> pl.Expr:
    """Return a Float64 expression; unreadable values become null."""
    if frame.schema.get(column) in FLOAT_DTYPES:
        return pl.col(column)
    return pl.col(column).cast(pl.Float64, strict=False)


def count_matching_rows(
    frame: pl.DataFrame,
    condition: pl.Expr,
    *,
    null_is_match: bool = False,
) -> int:
    """Count a boolean condition using caller-selected null polarity."""
    if frame.height == 0:
        return 0
    value = frame.select(
        condition.fill_null(null_is_match).cast(pl.UInt64).sum()
    ).item()
    return int(value or 0)


def count_non_null(frame: pl.DataFrame, column: str) -> int:
    """Count populated cells in one column."""
    return int(frame.select(pl.col(column).is_not_null().sum()).item() or 0)


def chromosome_expression(
    column: str,
    dtype: pl.DataType,
    *,
    strip_chr_prefix: bool = True,
    strip_leading_zero: bool = False,
) -> pl.Expr:
    """Normalize a chromosome column without turning integral floats into ``7.0``."""
    expression = pl.col(column)
    if dtype in FLOAT_DTYPES:
        expression = expression.cast(pl.Int64, strict=False)
    expression = expression.cast(pl.String).str.strip_chars()
    if strip_chr_prefix:
        expression = expression.str.replace(r"(?i)^chr", "")
    if strip_leading_zero:
        expression = expression.str.replace(r"^0", "")
    return expression.alias(column)


def position_expression(column: str, dtype: pl.DataType) -> pl.Expr:
    """Normalize a position column to Int64; unreadable cells become null."""
    expression = pl.col(column)
    if dtype == pl.String:
        expression = expression.str.strip_chars()
    return expression.cast(pl.Int64, strict=False).alias(column)


def validate_cast_retention(
    before: int,
    after: int,
    label: str,
    *,
    error_type: Type[Exception] = ValueError,
    warn: Callable[[str], None] | None = None,
) -> None:
    """Reject a destructive cast and report partial value loss consistently."""
    if before > 0 and after == 0:
        raise error_type(
            "Every %s value became empty during numeric conversion (%d values "
            "before, 0 after). Check the source column format." % (label, before)
        )
    if after < before and warn is not None:
        warn(
            "%d of %d %s values could not be converted and became empty."
            % (before - after, before, label)
        )


__all__ = [
    "FLOAT_DTYPES",
    "chromosome_expression",
    "collect_streaming",
    "count_matching_rows",
    "count_non_null",
    "numeric_column",
    "position_expression",
    "streaming_collect_options",
    "validate_cast_retention",
]
