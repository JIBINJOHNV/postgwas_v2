"""Small scientific transformations shared across PostGWAS modules."""

from __future__ import annotations

import math

import numpy as np
import polars as pl


def negative_log10_to_raw_p(
    column: str,
    minimum_p_value: float | None = None,
    *,
    output_name: str = "P",
) -> pl.Expr:
    """Return ``10**(-LP)``, optionally bounded for a target file contract.

    Harmonisation omits ``minimum_p_value`` so IEEE-754 underflow remains
    visible and can be distinguished from a study-reported zero using the
    preserved source token. Formatter callers may supply their canonical
    target-format minimum; that explicit output policy retains the established
    bounded-export behaviour without imposing the bound on internal analysis.
    """
    value = pl.col(column).cast(pl.Float64, strict=False)
    if minimum_p_value is None:
        return (10.0 ** (-value)).alias(output_name)

    minimum = float(minimum_p_value)
    if not 0 < minimum <= 1:
        raise ValueError("minimum_p_value must be greater than 0 and at most 1")
    maximum_lp = -math.log10(minimum)
    return (
        pl.when(value > maximum_lp)
        .then(pl.lit(minimum))
        .otherwise(10.0 ** (-value))
        .alias(output_name)
    )


def normal_z_magnitude_from_ln_p(values, tail: int) -> np.ndarray:
    """Return the normal-test ``|Z|`` associated with natural-log p-values.

    For a ``tail``-sided test, ``p / tail = Phi(-|Z|)``. ``ndtri_exp`` accepts
    ``log(Phi)`` directly, so values far below the smallest positive Float64
    probability retain their correct normal quantile.
    """
    from scipy.special import ndtri_exp

    tail_count = int(tail)
    if tail_count <= 0:
        raise ValueError("tail must be a positive integer")
    log_probabilities = np.asarray(values, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return -ndtri_exp(log_probabilities - math.log(tail_count))


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


__all__ = [
    "adjust_p_values",
    "negative_log10_to_raw_p",
    "normal_z_magnitude_from_ln_p",
]
