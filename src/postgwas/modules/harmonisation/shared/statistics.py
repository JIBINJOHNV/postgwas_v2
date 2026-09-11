"""Numerically stable statistical identities shared by harmonisation steps."""

import math

import polars as pl


def valid_frequency_mask(value):
    """Return a null-safe mask for finite probabilities in the closed unit interval.

    Frequencies of exactly zero or one are valid probabilities at this general
    validation layer.  Scientific steps that cannot use degenerate frequencies,
    such as Z-based effect reconstruction, apply their stricter policy later.
    """
    numeric = value.cast(pl.Float64, strict=False)
    return (
        numeric.is_not_null()
        & numeric.is_finite()
        & numeric.is_between(0.0, 1.0)
    ).fill_null(False)


def frequency_maf_screen(frame, column, cutoff):
    """Describe whether one raw frequency column is MAF-like.

    The screen is deliberately representation-only: it uses finite numeric
    probabilities and asks what fraction are at or below one half.  It does not
    decide allele semantics.  Callers must confirm a MAF-like result against an
    independent, allele-aligned reference before labelling the column EAF or
    MAF.
    """
    if column not in frame.columns:
        return {
            "column_present": False,
            "rows": frame.height,
            "usable": 0,
            "unusable": frame.height,
            "at_or_below_0_5": 0,
            "low_fraction_of_usable": 0.0,
            "decision_cutoff": float(cutoff),
            "maf_like": False,
        }
    value = pl.col(column).cast(pl.Float64, strict=False)
    usable = valid_frequency_mask(value)
    evidence = frame.select([
        pl.len().alias("rows"),
        usable.sum().alias("usable"),
        (usable & (value <= 0.5)).sum().alias("at_or_below_0_5"),
    ]).to_dicts()[0]
    usable_count = int(evidence["usable"] or 0)
    low_count = int(evidence["at_or_below_0_5"] or 0)
    fraction = low_count / usable_count if usable_count else 0.0
    return {
        "column_present": True,
        "rows": int(evidence["rows"] or 0),
        "usable": usable_count,
        "unusable": int(evidence["rows"] or 0) - usable_count,
        "at_or_below_0_5": low_count,
        "low_fraction_of_usable": fraction,
        "decision_cutoff": float(cutoff),
        "maf_like": bool(usable_count and fraction > float(cutoff)),
    }


def two_sided_negative_log10_p_from_z(values):
    """Return ``-log10(2 * P(N(0,1) > |Z|))`` without tail underflow."""
    import numpy as np
    from scipy.special import log_ndtr

    z_values = np.asarray(values, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return -(
            math.log(2.0) + log_ndtr(-np.abs(z_values))
        ) / math.log(10.0)


__all__ = [
    "frequency_maf_screen",
    "two_sided_negative_log10_p_from_z",
    "valid_frequency_mask",
]
