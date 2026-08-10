"""Small scientific transformations shared across PostGWAS modules."""

from __future__ import annotations

import math

import numpy as np
import polars as pl


def negative_log10_to_raw_p(
    column: str,
    minimum_p_value: float,
    *,
    output_name: str = "P",
) -> pl.Expr:
    """Return a numerically bounded ``10**(-LP)`` Polars expression."""
    minimum = float(minimum_p_value)
    if not 0 < minimum <= 1:
        raise ValueError("minimum_p_value must be greater than 0 and at most 1")
    maximum_lp = -math.log10(minimum)
    value = pl.col(column).cast(pl.Float64, strict=False)
    return (
        pl.when(value > maximum_lp)
        .then(pl.lit(minimum))
        .otherwise(10.0 ** (-value))
        .alias(output_name)
    )


def adjust_p_values(values, method: str) -> np.ndarray:
    """Apply a standard family-wise or false-discovery-rate correction."""
    p_values = np.asarray(values, dtype=float)
    if p_values.ndim != 1 or np.any(~np.isfinite(p_values)) or np.any(
        (p_values < 0) | (p_values > 1)
    ):
        raise ValueError("p-values must be a one-dimensional finite array in [0, 1]")
    count = len(p_values)
    if count == 0:
        return p_values.copy()
    if method == "bonferroni":
        return np.minimum(p_values * count, 1.0)
    if method == "sidak":
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.minimum(-np.expm1(count * np.log1p(-p_values)), 1.0)

    order = np.argsort(p_values, kind="mergesort")
    ranked = p_values[order]
    if method == "holm":
        adjusted = np.maximum.accumulate((count - np.arange(count)) * ranked)
    elif method == "fdr_bh":
        adjusted = np.minimum.accumulate(
            (ranked * count / np.arange(1, count + 1))[::-1]
        )[::-1]
    else:
        raise ValueError("unsupported p-value correction method: %s" % method)
    result = np.empty(count, dtype=float)
    result[order] = np.minimum(adjusted, 1.0)
    return result


__all__ = ["adjust_p_values", "negative_log10_to_raw_p"]
