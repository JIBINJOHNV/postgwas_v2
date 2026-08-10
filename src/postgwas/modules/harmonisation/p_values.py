"""P-value scale: raw p, -log10 p or -ln p (per-chromosome step 08).

Two entry points matter:

``detect_p_value_type(df, pval_col, policies)``
    Pure detection.  Reads the frame, changes nothing, returns
    ``(decision, evidence)`` where ``decision`` is ``'raw'``, ``'neglog10'`` or
    ``'negln'``.  The dataset-level "resolve study properties" step calls this
    ONCE on the full input so every chromosome applies the same answer.

``harmonise_p_values(chromosome, df, sample_column_dict, ...)``
    Applies a decision - the study-level one when handed in through
    ``decision=``, otherwise one detected from this chromosome alone - and
    leaves a plain p-value column behind, exactly as before.

The natural-log case is new.  The old classifier had two outcomes, and both of
its last two branches returned the same one, so a ``-ln p`` column was labelled
``-log10 p`` and had ``10 ** (-x)`` applied to it: every p-value wrong by about
four orders of magnitude, silently.  The medians are the discriminator - about
0.5 for raw p, 0.30 for -log10 p, 0.69 for -ln p - and a median that sits
between the two log scales is reported as ambiguous and raises rather than
being guessed at.

Two internal boolean columns, ``__pval_clipped_high`` and ``__pval_clipped_low``,
mark values forced back into the valid range. ``__pval_reported_zero`` separately
preserves the fact that a raw input value was exactly zero after the value itself
is raised to the configured floating-point floor. A reported p-value of exactly 1
is valid, remains unchanged, and is never flagged. The working columns survive
until effect validation and are dropped before export.
"""

import math
from typing import Any, Dict, Optional, Tuple, Union

import polars as pl

from postgwas.core.values import format_count as _count, format_number as _fmt
from postgwas.core.statistics import negative_log10_to_raw_p

from .shared.runtime import NullStepContext, resolve_policies, step_context


__all__ = [
    "AmbiguousPValueTypeError",
    "CLIPPED_COLUMNS",
    "CLIPPED_HIGH_COLUMN",
    "CLIPPED_LOW_COLUMN",
    "INTERNAL_PVALUE_COLUMNS",
    "REPORTED_ZERO_COLUMN",
    "POLICY_KEYS",
    "PValueTypeError",
    "STEP_LABEL",
    "convert_negative_log10_to_p_value",
    "convert_p_value_to_negative_log10",
    "harmonise_p_values",
    "detect_p_value_type",
]


# The label written into the reject file's `reject_step` column, and the
# position of this step in the per-chromosome pipeline.
STEP_LABEL = "07 pvalue_detection"
STEP_NUMBER = 7
STEP_TOTAL = 16
STEP_TITLE = "P-value scale"

POLICY_KEYS = [
    "pvalue.type",
    "pvalue.mlogp_detect_threshold",
    "pvalue.mlogp_detect_proportion",
    "pvalue.tolerance_above_one",
    "pvalue.out_of_range",
    "pvalue.clip_low",
    "pvalue.clip_high",
    "pvalue.mlogp_min",
    "pvalue.mlogp_max",
    "pvalue.verify_per_chromosome",
]

# Internal flags consumed by the effect-statistics validation step and dropped
# before export.
CLIPPED_HIGH_COLUMN = "__pval_clipped_high"
CLIPPED_LOW_COLUMN = "__pval_clipped_low"
CLIPPED_COLUMNS = (CLIPPED_HIGH_COLUMN, CLIPPED_LOW_COLUMN)
REPORTED_ZERO_COLUMN = "__pval_reported_zero"
INTERNAL_PVALUE_COLUMNS = CLIPPED_COLUMNS + (REPORTED_ZERO_COLUMN,)

# Constants of the arithmetic, not settings.  Under the null hypothesis half the
# p-values are below 0.5, so the median of -log10 p is log10(2) and the median of
# -ln p is ln(2).  Real data is enriched and sits a little above these, which is
# why the decision is made on the ratio and why the middle ground is refused.
LN10 = math.log(10.0)
MEDIAN_NEGLOG10 = math.log10(2.0)        # 0.30103
MEDIAN_NEGLN = math.log(2.0)             # 0.69315
MEDIAN_MIDPOINT = math.sqrt(MEDIAN_NEGLOG10 * MEDIAN_NEGLN)   # 0.45673
MEDIAN_AMBIGUITY_FACTOR = 1.15           # a median within x1.15 of the midpoint
                                         # is too close to call
# How far, as a ratio, the observed median may sit from an expected one and
# still count as that scale.  Half the separation between the two expectations,
# so the two acceptance windows meet at the midpoint and cover 0.198 to 1.051.
# A median outside both - a file of top hits only, say - fits neither scale and
# is reported as undecidable rather than guessed at.
MEDIAN_ACCEPT_LOG_DISTANCE = 0.5 * abs(math.log(MEDIAN_NEGLN / MEDIAN_NEGLOG10))
# A p-value can never exceed 1, so anything past this is evidence of a log
# scale on its own.  This is the `max_val > 2` rule the old detector used.
MAX_RAW_PVALUE_EVIDENCE = 2.0

_SCALE_TEXT = {
    "raw": "plain p-values",
    "neglog10": "-log10 p",
    "negln": "-ln p (natural log)",
}


class PValueTypeError(ValueError):
    """The p-value column cannot be used as configured."""


class AmbiguousPValueTypeError(PValueTypeError):
    """The column is on a log scale but which log cannot be determined."""


class _Unset(object):
    def __repr__(self):
        return "<from policies>"


_UNSET = _Unset()


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def _normalise_decision(value):
    """Return a canonical p-value scale decision."""
    if value is None:
        return None
    if isinstance(value, dict):
        for key in ("pvalue_type", "decision"):
            if value.get(key):
                value = value[key]
                break
        else:
            return None
    text = str(value).strip().lower()
    if text in ("", "auto", "na", "none"):
        return None
    if text == "raw":
        return "raw"
    if text == "neglog10":
        return "neglog10"
    if text == "negln":
        return "negln"
    raise PValueTypeError(
        "unknown p-value scale %r - expected 'raw', 'neglog10' or 'negln'" % (value,)
    )


def _drop_or_reject(df, mask, reason, detail, plain, check, ctx, rejects, warn=True):
    """Remove the rows the mask selects, recorded either way.

    ``mask`` must already be null-free (nulls are handled before this is
    called).  With a reject collector the rows go into the reject file and the
    collector does its own before/after logging; without one they are filtered
    out and the QC action is logged here so the count is never silent.
    """
    before = df.height
    if rejects is not None:
        return rejects.reject(df, mask, STEP_LABEL, reason, detail=detail)
    df = df.filter(~mask)
    ctx.qc(check, plain, before, df.height, reason=reason, warn=warn, step=STEP_LABEL)
    return df


# ---------------------------------------------------------------------------
# 1.  Pure detection  (study level, plan 5.1)
# ---------------------------------------------------------------------------
def detect_p_value_type(
    df: pl.DataFrame,
    pval_col: str,
    policies: Optional[Any] = None,
    statistics: Optional[Dict[str, Any]] = None,
) -> Tuple[str, Dict[str, Any]]:
    """Decide the scale of ``pval_col``.

    Returns ``(decision, evidence)`` with ``decision`` in
    ``{'raw', 'neglog10', 'negln'}``.  Raises ``AmbiguousPValueTypeError`` when
    the column is clearly on a log scale but the median lies between the two
    log scales, rather than picking one and being wrong by a factor of
    ``ln(10)`` in the exponent. The optional ``statistics`` input lets the
    dataset-level resolver reuse the canonical aggregates from its combined scan.
    """
    policies = resolve_policies(policies)

    if not pval_col or pval_col not in df.columns:
        raise PValueTypeError(
            "P-value column %r is not in the data frame; the columns present "
            "are: %s" % (pval_col, ", ".join(df.columns))
        )

    threshold = policies.pvalue.mlogp_detect_threshold
    proportion = policies.pvalue.mlogp_detect_proportion

    row = statistics
    if row is None:
        row = df.select(
            _pvalue_type_statistic_expressions(pval_col, threshold)
        ).to_dicts()[0]

    n_rows = int(row["n_rows"] or 0)
    n_usable = int(row["n_usable"] or 0)
    if n_usable == 0:
        raise PValueTypeError(
            "No numeric P-values found for detection in column %r (%d rows)."
            % (pval_col, n_rows)
        )

    above_fraction = int(row["n_above_threshold"] or 0) / n_usable
    within_fraction = int(row["n_within"] or 0) / n_usable
    median = row["median"]
    maximum = row["max"]

    evidence = {
        "pval_col": pval_col,
        "n_rows": n_rows,
        "n_non_null": int(row["n_non_null"] or 0),
        "n_usable": n_usable,
        "n_negative": int(row["n_negative"] or 0),
        "above_threshold_fraction": above_fraction,
        "within_range_fraction": within_fraction,
        "detect_threshold": threshold,
        "detect_proportion": proportion,
        "min": row["min"],
        "max": maximum,
        "mean": row["mean"],
        "median": median,
        "expected_median_neglog10": MEDIAN_NEGLOG10,
        "expected_median_negln": MEDIAN_NEGLN,
    }

    configured = policies.pvalue.type
    if configured != "auto":
        evidence["decision_source"] = "config (pvalue.type)"
        return configured, evidence

    evidence["decision_source"] = "automatic detection"

    # Raw against log is decided by the SHARE of values above 1, not by the
    # largest one.  A p-value cannot exceed 1, while a -log10 p column has
    # roughly 10% of its values above 1 and a -ln p column roughly 37%, so the
    # share separates the hypotheses cleanly.  The old detector also declared a
    # log scale whenever the maximum exceeded 2, which turned a raw column with
    # a single corrupt cell into a log column and applied 10**(-p) to all of
    # it; the maximum is kept as evidence and the offending values are dealt
    # with by the range gate instead.
    on_log_scale = above_fraction >= proportion
    evidence["on_log_scale"] = on_log_scale
    evidence["max_above_raw_range"] = bool(
        maximum is not None and maximum > MAX_RAW_PVALUE_EVIDENCE
    )

    if not on_log_scale:
        return "raw", evidence

    # A log scale: which log?
    if median is None or median <= 0.0:
        raise AmbiguousPValueTypeError(
            "Column %r is on a log scale (%.4f%% of values exceed %s, maximum %s) "
            "but its median is %s, so -log10 p and -ln p cannot be told apart. "
            "Set pvalue.type in the config."
            % (
                pval_col,
                above_fraction * 100.0,
                _fmt(threshold),
                _fmt(maximum),
                _fmt(median),
            )
        )

    band_low = MEDIAN_MIDPOINT / MEDIAN_AMBIGUITY_FACTOR
    band_high = MEDIAN_MIDPOINT * MEDIAN_AMBIGUITY_FACTOR
    accept_low = MEDIAN_NEGLOG10 / math.exp(MEDIAN_ACCEPT_LOG_DISTANCE)
    accept_high = MEDIAN_NEGLN * math.exp(MEDIAN_ACCEPT_LOG_DISTANCE)
    evidence["ambiguous_median_band"] = [band_low, band_high]
    evidence["plausible_median_range"] = [accept_low, accept_high]

    if median < accept_low or median > accept_high:
        # Neither scale predicts this median.  The usual cause is a file that
        # holds only the top hits, where the median says nothing about the
        # scale.  Refuse rather than pick.
        raise AmbiguousPValueTypeError(
            "Column %r is on a log scale, but its median of %s fits neither "
            "-log10 p (expected about %s) nor -ln p (expected about %s); "
            "anything outside %s to %s is not decidable from the median, which "
            "usually means the file holds only the strongest associations. "
            "Guessing would make every p-value wrong by a factor of ln(10) in "
            "the exponent, so nothing was guessed. Set pvalue.type to "
            "'neglog10' or 'negln' in the config."
            % (
                pval_col,
                _fmt(median, "%.4f"),
                _fmt(MEDIAN_NEGLOG10, "%.4f"),
                _fmt(MEDIAN_NEGLN, "%.4f"),
                _fmt(accept_low, "%.4f"),
                _fmt(accept_high, "%.4f"),
            )
        )

    if band_low <= median <= band_high:
        raise AmbiguousPValueTypeError(
            "Column %r is on a log scale, but its median of %s sits between the "
            "%s expected of -log10 p and the %s expected of -ln p (the "
            "undecidable band is %s to %s). Guessing would make every p-value "
            "wrong by a factor of ln(10) in the exponent, so nothing was "
            "guessed. Set pvalue.type to 'neglog10' or 'negln' in the config."
            % (
                pval_col,
                _fmt(median, "%.4f"),
                _fmt(MEDIAN_NEGLOG10, "%.4f"),
                _fmt(MEDIAN_NEGLN, "%.4f"),
                _fmt(band_low, "%.4f"),
                _fmt(band_high, "%.4f"),
            )
        )

    if median < MEDIAN_MIDPOINT:
        return "neglog10", evidence
    return "negln", evidence


def _pvalue_type_statistic_expressions(
    pval_col: str,
    threshold: float,
    prefix: str = "",
):
    """Return the canonical aggregate expressions used by p-value detection."""
    values = pl.col(pval_col).cast(pl.Float64, strict=False)
    usable = values.filter(values.is_finite())
    return [
        pl.len().alias(prefix + "n_rows"),
        values.is_not_null().sum().alias(prefix + "n_non_null"),
        values.is_finite().sum().alias(prefix + "n_usable"),
        (usable > threshold).sum().alias(prefix + "n_above_threshold"),
        ((usable >= 0.0) & (usable <= threshold)).sum().alias(prefix + "n_within"),
        (usable < 0.0).sum().alias(prefix + "n_negative"),
        usable.min().alias(prefix + "min"),
        usable.max().alias(prefix + "max"),
        usable.mean().alias(prefix + "mean"),
        usable.median().alias(prefix + "median"),
    ]


def _detection_sentence(evidence):
    return (
        "%s of %s usable values exceed %s (median %s, maximum %s)"
        % (
            "%.3f%%" % (evidence.get("above_threshold_fraction", 0.0) * 100.0),
            _count(evidence.get("n_usable")),
            _fmt(evidence.get("detect_threshold")),
            _fmt(evidence.get("median"), "%.4f"),
            _fmt(evidence.get("max")),
        )
    )


# ---------------------------------------------------------------------------
# 3.  Raw p  ->  -log10 p
# ---------------------------------------------------------------------------
def convert_p_value_to_negative_log10(
    df: pl.DataFrame,
    sample_column_dict,
    output_col: str = "LP",
    min_p: Optional[float] = None,
    max_p: Optional[float] = None,
    logger: Optional[Any] = None,
    policies: Optional[Any] = None,
    ctx: Optional[Any] = None,
    rejects: Optional[Any] = None,
):
    """Turn plain p-values into ``-log10 p`` in ``output_col``.

    ``min_p`` / ``max_p`` default to ``pvalue.clip_low`` and
    ``pvalue.clip_high``. The canonical upper bound is 1 because p=1 is valid.
    """
    policies = resolve_policies(policies)
    ctx = ctx if ctx is not None else NullStepContext(df.height)
    if min_p is None:
        min_p = policies.pvalue.clip_low
    if max_p is None:
        max_p = policies.pvalue.clip_high

    pval_col = sample_column_dict["pval_col"]
    if pval_col not in df.columns:
        raise PValueTypeError("Column '%s' not found in DataFrame." % (pval_col,))

    values = pl.col(pval_col).cast(pl.Float64, strict=False)
    counts = df.select(
        [
            values.is_null().sum().alias("n_null"),
            ((values <= 0.0) | (values > 1.0)).sum().alias("n_out"),
        ]
    ).to_dicts()[0]
    n_null = int(counts["n_null"] or 0)
    n_out = int(counts["n_out"] or 0)

    if n_null:
        df = _drop_or_reject(
            df,
            pl.col(pval_col).is_null(),
            "pval_null",
            None,
            "%s variants have no p-value, so -log10 p cannot be computed for "
            "them." % (_count(n_null),),
            "p-value present",
            ctx,
            rejects,
        )
    if n_out:
        df = _drop_or_reject(
            df,
            ((pl.col(pval_col) <= 0.0) | (pl.col(pval_col) > 1.0)).fill_null(False),
            "pval_out_of_range",
            "p is above 0 and at most 1",
            "%s variants have a p-value of 0 or below, or above 1, which has "
            "no usable logarithm." % (_count(n_out),),
            "p-value range",
            ctx,
            rejects,
        )

    df = df.with_columns(pl.col(pval_col).clip(min_p, max_p).alias("__p_clipped"))
    df = df.with_columns(
        (-pl.col("__p_clipped").log10()).clip(0.0, None).alias(output_col)
    ).drop("__p_clipped")

    stats = df.select(
        pl.col(output_col).min().alias("min_LP"),
        pl.col(output_col).max().alias("max_LP"),
        pl.col(output_col).mean().alias("mean_LP"),
        pl.col(output_col).median().alias("median_LP"),
    ).to_dicts()[0]

    ctx.info(
        "Converted plain p-values to -log10 p in '%s': minimum %s, maximum %s, "
        "median %s."
        % (
            output_col,
            _fmt(stats["min_LP"]),
            _fmt(stats["max_LP"]),
            _fmt(stats["median_LP"]),
        )
    )

    sample_column_dict["pval_col"] = output_col
    return df, stats, sample_column_dict


# ---------------------------------------------------------------------------
# 4.  -log10 p  ->  raw p
# ---------------------------------------------------------------------------
def convert_negative_log10_to_p_value(
    df: pl.DataFrame,
    sample_column_dict,
    output_col: str = "PVAL",
    min_lp: Optional[float] = None,
    max_lp: Union[float, None, "_Unset"] = _UNSET,
    logger: Optional[Any] = None,
    policies: Optional[Any] = None,
    ctx: Optional[Any] = None,
):
    """Turn ``-log10 p`` back into plain p-values in ``output_col``.

    ``min_lp`` defaults to ``pvalue.mlogp_min`` (0.0) and ``max_lp`` to
    ``pvalue.mlogp_max`` (300.0, or ``null`` to turn the cap off).  Both caps
    are now reported: capping used to be silent, which made -log10 p of 320,
    500 and 1200 indistinguishable in the output and left no trace of it.

    Adds ``__pval_clipped_high`` / ``__pval_clipped_low`` so the
    effect-statistics validation step can see which p-values were forced.
    """
    policies = resolve_policies(policies)
    ctx = ctx if ctx is not None else NullStepContext(df.height)
    if min_lp is None:
        min_lp = policies.pvalue.mlogp_min
    if isinstance(max_lp, _Unset):
        max_lp = policies.pvalue.mlogp_max

    pval_col = sample_column_dict["pval_col"]
    if pval_col not in df.columns:
        raise PValueTypeError("Column '%s' not found in DataFrame." % (pval_col,))

    values = pl.col(pval_col).cast(pl.Float64, strict=False)
    low_mask = (values < min_lp) if min_lp is not None else pl.lit(False)
    high_mask = (values > max_lp) if max_lp is not None else pl.lit(False)

    counts = df.select(
        [
            low_mask.sum().alias("n_below"),
            high_mask.sum().alias("n_above"),
            values.max().alias("max_lp_seen"),
        ]
    ).to_dicts()[0]
    n_below = int(counts["n_below"] or 0)
    n_above = int(counts["n_above"] or 0)
    n_clipped = n_below + n_above

    df = df.with_columns(pl.col(pval_col).clip(min_lp, max_lp).alias("__lp_clipped"))

    if n_above:
        # A -log10 p above the cap is a genuine top hit being flattened.
        ctx.qc(
            "-log10 p cap",
            "%s variants have a -log10 p above the cap of %s (the largest seen "
            "was %s). Their p-values were all set to %s, so they can no longer "
            "be told apart. Policy 'pvalue.mlogp_max' controls this; null turns "
            "the cap off."
            % (
                _count(n_above),
                _fmt(max_lp),
                _fmt(counts["max_lp_seen"]),
                _fmt(10.0 ** (-max_lp) if max_lp is not None else None, "%.3e"),
            ),
            df.height,
            df.height,
            changed=n_above,
            warn=True,
            step=STEP_LABEL,
        )
    if n_below:
        ctx.qc(
            "-log10 p floor",
            "%s variants have a -log10 p below %s, which means a p-value above "
            "1. They were raised to %s, i.e. p = %s."
            % (
                _count(n_below),
                _fmt(min_lp),
                _fmt(min_lp),
                _fmt(10.0 ** (-min_lp) if min_lp is not None else None),
            ),
            df.height,
            df.height,
            changed=n_below,
            warn=True,
            step=STEP_LABEL,
        )

    df = df.with_columns(
        negative_log10_to_raw_p(
            "__lp_clipped", policies.pvalue.clip_low, output_name="__raw_p",
        )
    )

    # 10 ** -LP underflows to exactly 0 for very large LP; that is a p-value
    # clipped at the bottom just as surely as the cap above.
    underflow = pl.col("__raw_p") == 0
    df = df.with_columns(
        [
            pl.when(underflow)
            .then(pl.lit(policies.pvalue.clip_low))
            .otherwise(pl.col("__raw_p"))
            .alias(output_col),
            (low_mask.fill_null(False)).alias(CLIPPED_HIGH_COLUMN),
            (high_mask.fill_null(False) | underflow).alias(CLIPPED_LOW_COLUMN),
        ]
    ).drop(["__lp_clipped", "__raw_p"])

    mlogp_stats = df.select(
        pl.col(output_col).min().alias("min_PVAL"),
        pl.col(output_col).max().alias("max_PVAL"),
        pl.col(output_col).mean().alias("mean_PVAL"),
        pl.col(output_col).median().alias("median_PVAL"),
    ).to_dicts()[0]

    mlogp_stats["mlogp_values_total_clipped"] = n_clipped
    mlogp_stats["mlogp_values_clipped_at_max"] = n_above
    mlogp_stats["mlogp_values_clipped_at_min"] = n_below
    mlogp_stats["mlogp_max_observed"] = counts["max_lp_seen"]

    ctx.info(
        "Converted -log10 p to plain p-values in '%s': minimum %s, maximum %s, "
        "median %s."
        % (
            output_col,
            _fmt(mlogp_stats["min_PVAL"], "%.3e"),
            _fmt(mlogp_stats["max_PVAL"], "%.3e"),
            _fmt(mlogp_stats["median_PVAL"], "%.3e"),
        )
    )

    sample_column_dict["pval_col"] = output_col
    return df, mlogp_stats, sample_column_dict


# ---------------------------------------------------------------------------
# 5.  The step itself
# ---------------------------------------------------------------------------
def _harmonise_raw_pvalues(
    df, pval_col, output_col, policies, ctx, rejects, qc_info, chromosome,
    entry_counts=None,
):
    """Bring plain p-values into the allowed range, recording every removal.

    Three gates, in this order:

    1. a missing p-value is unusable                    -> ``pval_null``
    2. below 0 or above ``pvalue.tolerance_above_one``  -> not a p-value at all
    3. outside ``[clip_low, clip_high]``                -> ``pvalue.out_of_range``

    Gates 1 and 2 remove the variant, which is what this module always did -
    silently.  Gate 3 is the one ``pvalue.out_of_range`` governs, and its
    default, ``clip``, is exactly the clip that was hardcoded here before.
    """
    tolerance = policies.pvalue.tolerance_above_one
    clip_low = policies.pvalue.clip_low
    clip_high = policies.pvalue.clip_high
    action = policies.pvalue.out_of_range

    if entry_counts is None:
        values = pl.col(pval_col)
        counts = df.select(
            [
                values.is_null().sum().alias("n_null"),
                (values < 0.0).sum().alias("n_lt0"),
                (values > tolerance).sum().alias("n_gt_tolerance"),
                ((values >= 0.0) & (values < clip_low)).sum().alias("n_below_clip"),
                ((values > clip_high) & (values <= tolerance)).sum().alias("n_above_clip"),
            ]
        ).to_dicts()[0]
    else:
        counts = entry_counts
    n_null = int(counts["n_null"] or 0)
    n_lt0 = int(counts["n_lt0"] or 0)
    n_gt_tolerance = int(counts["n_gt_tolerance"] or 0)
    n_below_clip = int(counts["n_below_clip"] or 0)
    n_above_clip = int(counts["n_above_clip"] or 0)
    n_unsalvageable = n_lt0 + n_gt_tolerance
    n_out_of_range = n_below_clip + n_above_clip

    total_before = df.height

    if action == "fail" and (n_unsalvageable or n_out_of_range):
        raise PValueTypeError(
            "Chromosome %s: %s p-values are outside 0 to 1 (%s below 0, %s above "
            "%s) and %s more are outside the kept range %s to %s. Policy "
            "'pvalue.out_of_range' is 'fail'."
            % (
                chromosome,
                _count(n_unsalvageable),
                _count(n_lt0),
                _count(n_gt_tolerance),
                _fmt(tolerance),
                _count(n_out_of_range),
                _fmt(clip_low),
                _fmt(clip_high),
            )
        )

    # 1. missing
    if n_null:
        df = _drop_or_reject(
            df,
            pl.col(pval_col).is_null(),
            "pval_null",
            None,
            "%s variants have no p-value at all." % (_count(n_null),),
            "p-value present",
            ctx,
            rejects,
        )

    # 2. impossible values
    unsalvageable = (
        (pl.col(pval_col) < 0.0) | (pl.col(pval_col) > tolerance)
    ).fill_null(False)
    if n_unsalvageable:
        plain = (
            "%s variants have a p-value below 0 or above %s (%s below 0, %s "
            "above the tolerance), which no rounding error explains."
            % (
                _count(n_unsalvageable),
                _fmt(tolerance),
                _count(n_lt0),
                _count(n_gt_tolerance),
            )
        )
        if action == "null":
            before = df.height
            df = df.with_columns(
                pl.when(unsalvageable)
                .then(None)
                .otherwise(pl.col(pval_col))
                .alias(pval_col)
            )
            ctx.qc(
                "p-value range",
                "%s Policy 'pvalue.out_of_range' is 'null', so their p-value was "
                "blanked out and the variants kept." % (plain,),
                before,
                df.height,
                changed=n_unsalvageable,
                warn=True,
                step=STEP_LABEL,
            )
        else:
            df = _drop_or_reject(
                df,
                unsalvageable,
                "pval_out_of_range",
                "outside 0 to %s" % (tolerance,),
                plain,
                "p-value range",
                ctx,
                rejects,
            )

    # 3. inside the tolerance but outside the kept range
    too_low = ((pl.col(pval_col) >= 0.0) & (pl.col(pval_col) < clip_low)).fill_null(
        False
    )
    too_high = (pl.col(pval_col) > clip_high).fill_null(False)
    if action == "reject" and n_out_of_range:
        df = _drop_or_reject(
            df,
            too_low | too_high,
            "pval_out_of_range",
            "outside %s to %s" % (clip_low, clip_high),
            "%s variants have a p-value outside the kept range %s to %s (%s too "
            "small, %s above 1)."
            % (
                _count(n_out_of_range),
                _fmt(clip_low),
                _fmt(clip_high),
                _count(n_below_clip),
                _count(n_above_clip),
            ),
            "p-value range",
            ctx,
            rejects,
        )
        df = df.with_columns(
            [
                pl.lit(False).alias(CLIPPED_HIGH_COLUMN),
                pl.lit(False).alias(CLIPPED_LOW_COLUMN),
            ]
        )
    elif action == "null" and n_out_of_range:
        before = df.height
        df = df.with_columns(
            [
                pl.when(too_low | too_high)
                .then(None)
                .otherwise(pl.col(pval_col))
                .alias(pval_col),
                pl.lit(False).alias(CLIPPED_HIGH_COLUMN),
                pl.lit(False).alias(CLIPPED_LOW_COLUMN),
            ]
        )
        ctx.qc(
            "p-value range",
            "%s variants have a p-value outside %s to %s. Policy "
            "'pvalue.out_of_range' is 'null', so their p-value was blanked out "
            "and the variants kept."
            % (_count(n_out_of_range), _fmt(clip_low), _fmt(clip_high)),
            before,
            df.height,
            changed=n_out_of_range,
            warn=True,
            step=STEP_LABEL,
        )
    else:
        # 'clip' - the default, and what this module always did.  The flags are
        # what step 11 uses to find values corrected from above 1.
        before = df.height
        df = df.with_columns(
            [
                too_high.alias(CLIPPED_HIGH_COLUMN),
                too_low.alias(CLIPPED_LOW_COLUMN),
            ]
        )
        df = df.with_columns(pl.col(pval_col).clip(clip_low, clip_high).alias(pval_col))
        if n_out_of_range:
            ctx.qc(
                "p-value range",
                "%s variants have a p-value outside %s to %s (%s too small, %s "
                "above 1). Policy 'pvalue.out_of_range' is 'clip', so they were "
                "forced back into the mathematically valid range."
                % (
                    _count(n_out_of_range),
                    _fmt(clip_low),
                    _fmt(clip_high),
                    _count(n_below_clip),
                    _count(n_above_clip),
                ),
                before,
                df.height,
                changed=n_out_of_range,
                warn=True,
                step=STEP_LABEL,
            )

    df = df.with_columns(pl.col(pval_col).alias(output_col))

    qc_info.update(
        {
            "raw_total_variants_with_pvalues": total_before,
            "after_filter_variants_with_pvalues": df.height,
            "variants_with_null_pvalues_removed": n_null,
            "variants_with_lt0_pvalues_removed": n_lt0,
            "variants_with_gt1_05_pvalues_removed": n_gt_tolerance,
            "variants_with_zero_pvalues_replaced": n_below_clip,
            "variants_with_pvalues_clipped_low": n_below_clip,
            "variants_with_pvalues_clipped_high": n_above_clip,
            "pvalue_out_of_range_action": action,
            "pvalue_clip_low": clip_low,
            "pvalue_clip_high": clip_high,
            "pvalue_tolerance_above_one": tolerance,
        }
    )
    return df, qc_info


def harmonise_p_values(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict,
    output_col: str = "LP",
    proportion_threshold: Optional[float] = None,
    max_allowed_lp: Union[float, None, "_Unset"] = _UNSET,
    logger: Optional[Any] = None,
    policies: Optional[Any] = None,
    rejects: Optional[Any] = None,
    decision: Optional[Any] = None,
    step_number: Optional[int] = None,
    step_total: Optional[int] = None,
):
    """Detect the p-value scale (or apply the study-level decision) and leave a
    plain p-value column behind.

    ``proportion_threshold`` and ``max_allowed_lp`` now default to
    ``pvalue.mlogp_detect_proportion`` (0.001) and ``pvalue.mlogp_max`` (300.0),
    the same numbers that used to be written into the signature.  Passing
    ``max_allowed_lp=None`` turns the cap off.

    Returns ``(df, qc_info, sample_column_dict)`` exactly as before.
    """
    policies = resolve_policies(policies)
    if proportion_threshold is not None:
        policies = policies.with_overrides(
            {"pvalue.mlogp_detect_proportion": proportion_threshold}
        )
    number = STEP_NUMBER if step_number is None else step_number
    total_steps = STEP_TOTAL if step_total is None else step_total

    pval_col = sample_column_dict.get("pval_col")
    if not pval_col or pval_col not in df.columns:
        raise PValueTypeError(
            "P-value column ('pval_col') is missing from the config or not "
            "present in the data frame for chromosome %s." % (chromosome,)
        )

    with step_context(
        logger,
        None,
        number=number,
        total=total_steps,
        title=STEP_TITLE,
        operation="p_values.harmonise_p_values",
        rows_in=df.height,
        policy_keys=POLICY_KEYS,
    ) as ctx:
        df = df.with_columns(pl.col(pval_col).cast(pl.Float64, strict=False))

        initial_stats = df.select(
            [
                pl.len().alias("initial_total_variants_with_pvalues"),
                pl.col(pval_col).is_null().sum().alias("initial_variants_with_null_pvalues"),
                (pl.col(pval_col) == 0).sum().alias("initial_variants_with_zero_pvalues"),
                (pl.col(pval_col) < 0).sum().alias("initial_variants_with_lt0_pvalues"),
                (pl.col(pval_col) > 1).sum().alias("initial_variants_with_gt1_pvalues"),
                (pl.col(pval_col) > policies.pvalue.tolerance_above_one)
                .sum()
                .alias("initial_variants_with_gt1_05_pvalues"),
                (
                    (pl.col(pval_col) >= 0.0)
                    & (pl.col(pval_col) < policies.pvalue.clip_low)
                ).sum().alias("__n_below_clip"),
                (
                    (pl.col(pval_col) > policies.pvalue.clip_high)
                    & (
                        pl.col(pval_col)
                        <= policies.pvalue.tolerance_above_one
                    )
                ).sum().alias("__n_above_clip"),
                pl.col(pval_col).min().alias("initial_min_pvalues"),
                pl.col(pval_col).max().alias("initial_max_pvalues"),
                pl.col(pval_col).mean().alias("initial_mean_pvalues"),
                pl.col(pval_col).median().alias("initial_median_pvalues"),
            ]
        ).to_dicts()[0]
        raw_entry_counts = {
            "n_null": initial_stats["initial_variants_with_null_pvalues"],
            "n_lt0": initial_stats["initial_variants_with_lt0_pvalues"],
            "n_gt_tolerance": initial_stats["initial_variants_with_gt1_05_pvalues"],
            "n_below_clip": initial_stats.pop("__n_below_clip"),
            "n_above_clip": initial_stats.pop("__n_above_clip"),
        }
        qc_info = initial_stats.copy()
        qc_info["pvalue_column"] = pval_col

        ctx.info(
            "P-values on entry: %s missing, %s exactly zero, %s below zero, %s "
            "above one; minimum %s, median %s, maximum %s."
            % (
                _count(initial_stats["initial_variants_with_null_pvalues"]),
                _count(initial_stats["initial_variants_with_zero_pvalues"]),
                _count(initial_stats["initial_variants_with_lt0_pvalues"]),
                _count(initial_stats["initial_variants_with_gt1_pvalues"]),
                _fmt(initial_stats["initial_min_pvalues"], "%.3e"),
                _fmt(initial_stats["initial_median_pvalues"]),
                _fmt(initial_stats["initial_max_pvalues"]),
            )
        )

        # ------------------------------------------------------------------
        # Which decision applies?
        # ------------------------------------------------------------------
        study_decision = _normalise_decision(decision)
        evidence = None
        if study_decision is not None:
            scale = study_decision
            decision_source = "study-level decision (made once, before the split)"
            if policies.pvalue.verify_per_chromosome:
                # Verification may only warn.  A chromosome whose own values are
                # undecidable must not fail a run whose scale is already known.
                try:
                    local, evidence = detect_p_value_type(df, pval_col, policies)
                except PValueTypeError as exc:
                    local = None
                    ctx.warn(
                        "Re-checking the p-value scale on this chromosome alone "
                        "was inconclusive, which changes nothing: the "
                        "study-level decision '%s' was applied. %s"
                        % (scale, exc)
                    )
                else:
                    if local != scale:
                        ctx.warn(
                            "Re-checking on this chromosome alone suggests '%s', "
                            "but the study-level decision '%s' is authoritative "
                            "and was applied. %s"
                            % (local, scale, _detection_sentence(evidence))
                        )
                    else:
                        ctx.info(
                            "Re-checked on this chromosome: agrees with the "
                            "study-level decision '%s'." % (scale,)
                        )
                qc_info["pvalue_scale_chromosome_check"] = local
        else:
            scale, evidence = detect_p_value_type(df, pval_col, policies)
            decision_source = evidence.get("decision_source", "automatic detection")

        qc_info["detected_scale"] = scale
        qc_info["pvalue_decision_source"] = decision_source
        sample_column_dict["pvalue_type"] = scale
        if evidence is not None:
            qc_info.update(
                {
                    "n_values": evidence.get("n_usable"),
                    "over_1.2_fraction": evidence.get("above_threshold_fraction"),
                    "within_[0,1.05]_fraction": evidence.get("within_range_fraction"),
                    "median": evidence.get("median"),
                    "max": evidence.get("max"),
                }
            )

        ctx.decide(
            "P-value column '%s'" % (pval_col,),
            _detection_sentence(evidence) if evidence is not None else decision_source,
            "%s (%s)" % (_SCALE_TEXT[scale], decision_source),
        )
        if scale == "raw" and evidence is not None and evidence.get("max_above_raw_range"):
            ctx.warn(
                "The column is plain p-values, but its largest value is %s, "
                "which no p-value can be. Only %s of values exceed %s, so those "
                "are corrupt cells rather than a log scale; they are handled by "
                "the range check below."
                % (
                    _fmt(evidence.get("max")),
                    "%.4f%%" % (evidence.get("above_threshold_fraction", 0.0) * 100.0),
                    _fmt(evidence.get("detect_threshold")),
                )
            )

        # ------------------------------------------------------------------
        # Apply
        # ------------------------------------------------------------------
        if scale == "raw":
            # Preserve exact-zero provenance before the raw value is raised to
            # pvalue.clip_low. This is a raw-input fact: an underflow while
            # reconstructing a supplied -log10 p retains more information and
            # must not be confused with a study that literally reported zero.
            df = df.with_columns(
                (pl.col(pval_col) == 0).fill_null(False).alias(REPORTED_ZERO_COLUMN)
            )
            df, qc_info = _harmonise_raw_pvalues(
                df, pval_col, output_col, policies, ctx, rejects, qc_info,
                chromosome, entry_counts=raw_entry_counts,
            )
            sample_column_dict["pval_col"] = output_col
        else:
            df = df.with_columns(pl.lit(False).alias(REPORTED_ZERO_COLUMN))
            total_before = df.height
            scratch_col = None
            if scale == "negln":
                # Put it on the -log10 scale first, in a scratch column, so one
                # code path and one set of policies handle both logs and the
                # study's own column is left as it arrived.
                scratch_col = "__lp_from_negln"
                df = df.with_columns((pl.col(pval_col) / LN10).alias(scratch_col))
                sample_column_dict["pval_col"] = scratch_col
                ctx.info(
                    "The column is on the natural-log scale, so it was divided "
                    "by ln(10) = %s to become -log10 p before conversion. "
                    "Reading it directly as -log10 p would make every p-value "
                    "wrong by about four orders of magnitude." % (_fmt(LN10),)
                )
            df, mlogp_stats, sample_column_dict = convert_negative_log10_to_p_value(
                df=df,
                sample_column_dict=sample_column_dict,
                output_col="PVAL",
                min_lp=policies.pvalue.mlogp_min,
                max_lp=max_allowed_lp,
                policies=policies,
                ctx=ctx,
            )
            if scratch_col is not None and scratch_col in df.columns:
                df = df.drop(scratch_col)
            qc_info.update(
                {
                    "mlogp_total_variants_with_pvalues": total_before,
                    "after_mlogp_to_pval_variants_with_pvalues": df.height,
                }
            )
            qc_info.update(mlogp_stats)

        pval_col = sample_column_dict["pval_col"]

        # Both branches must leave the flags behind for the effect-statistics
        # step, even when nothing was clipped.
        for column in INTERNAL_PVALUE_COLUMNS:
            if column not in df.columns:
                df = df.with_columns(pl.lit(False).alias(column))

        final_stats = df.select(
            [
                pl.len().alias("final_total_variants_with_pvalues"),
                pl.col(pval_col).is_null().sum().alias("final_variants_with_null_pvalues"),
                (pl.col(pval_col) == 0).sum().alias("final_variants_with_zero_pvalues"),
                (pl.col(pval_col) < 0).sum().alias("final_variants_with_lt0_pvalues"),
                (pl.col(pval_col) > 1).sum().alias("final_variants_with_gt1_pvalues"),
                (pl.col(pval_col) > policies.pvalue.tolerance_above_one)
                .sum()
                .alias("final_variants_with_gt1_05_pvalues"),
                pl.col(pval_col).min().alias("final_min_pvalues"),
                pl.col(pval_col).max().alias("final_max_pvalues"),
                pl.col(pval_col).mean().alias("final_mean_pvalues"),
                pl.col(pval_col).median().alias("final_median_pvalues"),
                pl.col(CLIPPED_HIGH_COLUMN).sum().alias("pvalues_flagged_clipped_high"),
                pl.col(CLIPPED_LOW_COLUMN).sum().alias("pvalues_flagged_clipped_low"),
            ]
        ).to_dicts()[0]
        qc_info.update(final_stats)

        ctx.set_rows(df.height)

    return df, qc_info, sample_column_dict
