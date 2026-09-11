"""P-value scale: raw p or -log10 p (per-chromosome step 07).

Two entry points matter:

``detect_p_value_type(df, pval_col, policies)``
    Pure detection.  Reads the frame, changes nothing, returns
    ``(decision, evidence)`` where ``decision`` is ``'raw'`` or ``'neglog10'``.
    The dataset-level "resolve study properties" step calls this ONCE on the
    full input so every chromosome applies the same answer.

``harmonise_p_values(chromosome, df, sample_column_dict, ...)``
    Applies a decision - the study-level one when handed in through
    ``decision=``, otherwise one detected from this chromosome alone - and
    leaves one plain p-value column named by ``pvalue.output_column``. The
    GWAS-to-VCF boundary consumes that raw value and creates FORMAT/LP later.

Automatic detection recognises only the two standard GWAS representations.
Negative values are incompatible with both. A small configured fraction is
warned about and rejected row by row; a larger fraction stops the study before
chromosome processing because it is evidence of a wrong scale or mapping. A
-log10 decision requires agreement between the configured count, fraction,
and study-wide median policies; therefore one sentinel cannot select a scale
for the file. Signed logarithms and natural-log representations must be
converted before PostGWAS reads them.

The canonical public value remains a plain raw p-value. Private source-text,
natural-log-p and exact-raw-text columns preserve values that cannot be
represented as a positive Float64. These are numerical provenance, not a
second public p-value field. A positive source token such as ``1e-400`` is
therefore never confused with a study-reported zero.
"""

from decimal import Decimal, InvalidOperation, localcontext
import math
import sys
from typing import Any, Dict, Optional, Tuple

import polars as pl

from postgwas.core.values import format_count as _count, format_number as _fmt
from postgwas.core.statistics import negative_log10_to_raw_p

from .shared.runtime import (
    NullStepContext,
    reject_rows,
    resolve_policies,
    step_context,
)


__all__ = [
    "AmbiguousPValueTypeError",
    "CLIPPED_COLUMNS",
    "CLIPPED_HIGH_COLUMN",
    "CLIPPED_LOW_COLUMN",
    "EXACT_RAW_P_COLUMN",
    "ExcessiveNegativePValueError",
    "FLOAT_UNDERFLOW_COLUMN",
    "INTERNAL_PVALUE_COLUMNS",
    "LOG_P_COLUMN",
    "P_VALUE_EXACT_RAW_KEY",
    "P_VALUE_LOG_KEY",
    "P_VALUE_SOURCE_TEXT_KEY",
    "REPORTED_ZERO_COLUMN",
    "SOURCE_TEXT_COLUMN",
    "POLICY_KEYS",
    "PValueTypeError",
    "STEP_LABEL",
    "assess_negative_pvalues",
    "convert_negative_log10_to_p_value",
    "harmonise_p_values",
    "detect_p_value_type",
    "preserve_pvalue_source_text",
]


# The label written into the reject file's `reject_step` column, and the
# position of this step in the per-chromosome pipeline.
STEP_LABEL = "07 pvalue_detection"
STEP_NUMBER = 7
STEP_TOTAL = 16
STEP_TITLE = "P-value scale"

POLICY_KEYS = [
    "pvalue.type",
    "pvalue.output_column",
    "pvalue.mlogp_detect_threshold",
    "pvalue.mlogp_detect_proportion",
    "pvalue.mlogp_detect_min_count",
    "pvalue.mlogp_expected_median",
    "pvalue.mlogp_median_tolerance",
    "pvalue.max_negative_fraction",
    "pvalue.tolerance_above_one",
    "pvalue.out_of_range",
    "pvalue.clip_low",
    "pvalue.clip_high",
    "pvalue.mlogp_min",
    "pvalue.verify_per_chromosome",
]

# Private calculation and provenance columns. Boolean correction flags are
# consumed by effect validation; mapped source/log/exact columns remain until
# the adapter boundary and are excluded from its explicit scientific mapping.
CLIPPED_HIGH_COLUMN = "__pval_clipped_high"
CLIPPED_LOW_COLUMN = "__pval_clipped_low"
CLIPPED_COLUMNS = (CLIPPED_HIGH_COLUMN, CLIPPED_LOW_COLUMN)
REPORTED_ZERO_COLUMN = "__pval_reported_zero"
SOURCE_TEXT_COLUMN = "__pval_source_text"
LOG_P_COLUMN = "__pval_ln"
EXACT_RAW_P_COLUMN = "__pval_exact_raw"
FLOAT_UNDERFLOW_COLUMN = "__pval_float_underflow"

# Temporary source-status columns. They are materialised before either scale
# branch changes the selected input column, then consumed and removed by the
# shared canonical range/missingness gate. Keeping the categories separate is
# what prevents a valid Float64-underflow probability from being mistaken for
# a missing or non-finite source value after -log10 conversion.
SOURCE_MISSING_STATUS_COLUMN = "__pval_source_missing"
SOURCE_NONFINITE_STATUS_COLUMN = "__pval_source_nonfinite"
SOURCE_NEGATIVE_STATUS_COLUMN = "__pval_source_negative"
SOURCE_STATUS_COLUMNS = (
    SOURCE_MISSING_STATUS_COLUMN,
    SOURCE_NONFINITE_STATUS_COLUMN,
    SOURCE_NEGATIVE_STATUS_COLUMN,
)

P_VALUE_SOURCE_TEXT_KEY = "pvalue_source_text_col"
P_VALUE_LOG_KEY = "pvalue_log_col"
P_VALUE_EXACT_RAW_KEY = "pvalue_exact_raw_col"

INTERNAL_PVALUE_COLUMNS = CLIPPED_COLUMNS + (
    REPORTED_ZERO_COLUMN,
    SOURCE_TEXT_COLUMN,
    LOG_P_COLUMN,
    EXACT_RAW_P_COLUMN,
    FLOAT_UNDERFLOW_COLUMN,
) + SOURCE_STATUS_COLUMNS

# Decimal128 precision is far beyond the Float32 precision of GWAS-VCF LP and
# is used only for exceptional source tokens that Float64 cannot represent.
_DECIMAL_PRECISION = 34
_DECIMAL_EMIN = -999999999
_DECIMAL_EMAX = 999999999

_SCALE_TEXT = {
    "raw": "plain p-values",
    "neglog10": "-log10 p",
}


class PValueTypeError(ValueError):
    """The p-value column cannot be used as configured."""


class AmbiguousPValueTypeError(PValueTypeError):
    """Automatic evidence does not safely identify raw p or -log10 p."""


class ExcessiveNegativePValueError(PValueTypeError):
    """Too much of the study p-value column is negative to continue safely."""


def assess_negative_pvalues(
    statistics: Dict[str, Any],
    pval_col: str,
    policies: Optional[Any] = None,
) -> Dict[str, Any]:
    """Validate the study-wide fraction of negative usable p-value cells.

    The denominator is the finite numeric subset. Missing, unparseable, NaN and
    infinite cells cannot dilute evidence of a signed-log or wrongly mapped
    p-value column.
    """
    policies = resolve_policies(policies)
    n_usable = int(statistics.get("n_usable") or 0)
    n_negative = int(statistics.get("n_negative") or 0)
    negative_fraction = n_negative / n_usable if n_usable else 0.0
    maximum_fraction = float(policies.pvalue.max_negative_fraction)
    result = {
        "n_negative": n_negative,
        "negative_fraction": negative_fraction,
        "max_negative_fraction": maximum_fraction,
    }
    if n_negative and negative_fraction > maximum_fraction:
        value_word = "value" if n_negative == 1 else "values"
        raise ExcessiveNegativePValueError(
            "P-value column %r contains %s negative usable %s out of %s "
            "(%s), exceeding pvalue.max_negative_fraction=%s (%s). PostGWAS "
            "supports raw p-values in [0,1] or non-negative -log10(p). This "
            "negative fraction usually indicates signed log10(p), ln(p), an "
            "incorrect p-value-column mapping, or mixed/corrupt data. Convert "
            "the complete column to raw p or non-negative -log10(p), correct "
            "p_value_column/p_value_type if necessary, and rerun. No chromosome "
            "processing was started."
            % (
                pval_col,
                _count(n_negative),
                value_word,
                _count(n_usable),
                "%.4f%%" % (negative_fraction * 100.0),
                _fmt(maximum_fraction),
                "%.2f%%" % (maximum_fraction * 100.0),
            )
        )
    return result


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def preserve_pvalue_source_text(df: pl.DataFrame, sample_column_dict):
    """Keep the configured source token before any numeric normalization.

    The stable private mapping is carried through chromosome partition I/O and
    final validation. Direct library callers receive the same protection when
    they call :func:`harmonise_p_values` without the dataset reader.
    """
    pval_col = sample_column_dict.get("pval_col")
    if not pval_col or pval_col not in df.columns:
        return df, sample_column_dict

    configured = sample_column_dict.get(P_VALUE_SOURCE_TEXT_KEY)
    if configured and configured in df.columns:
        return df, sample_column_dict
    if SOURCE_TEXT_COLUMN in df.columns and SOURCE_TEXT_COLUMN != pval_col:
        raise PValueTypeError(
            "Input uses reserved internal p-value provenance column %r. Rename "
            "that input column before harmonisation." % SOURCE_TEXT_COLUMN
        )

    df = df.with_columns(
        pl.col(pval_col)
        .cast(pl.String, strict=False)
        .str.strip_chars()
        .replace("", None)
        .alias(SOURCE_TEXT_COLUMN)
    )
    sample_column_dict[P_VALUE_SOURCE_TEXT_KEY] = SOURCE_TEXT_COLUMN
    return df, sample_column_dict


def _attach_pvalue_source_status(df, pval_col, source_col):
    """Materialise source missing/non-finite/negative status before conversion.

    The ordinary path is fully vectorised. Decimal inspection is restricted to
    negative-signed numeric zeros so an extreme token such as ``-1e-400`` is
    not rounded to ``-0.0`` and then silently accepted as ``-log10(p) = 0``.
    A literal ``-0`` remains mathematically zero and is not marked negative.
    """
    collisions = [name for name in SOURCE_STATUS_COLUMNS if name in df.columns]
    if collisions:
        raise PValueTypeError(
            "Input uses reserved internal p-value source-status column(s): %s. "
            "Rename those input columns before harmonisation."
            % ", ".join(collisions)
        )

    row_column = "__pval_source_status_row"
    if row_column in df.columns:
        raise PValueTypeError(
            "Input uses reserved internal p-value provenance column %r. Rename "
            "that input column before harmonisation." % row_column
        )

    indexed = df.with_row_index(row_column)
    numeric = pl.col(pval_col).cast(pl.Float64, strict=False)
    signed_zero_candidates = indexed.filter(
        numeric.is_finite()
        & (numeric == 0.0)
        & pl.col(source_col).cast(pl.String, strict=False)
        .str.starts_with("-")
        .fill_null(False)
    ).select(row_column, source_col)
    negative_underflow_indices = []
    for index, token in signed_zero_candidates.iter_rows():
        value = _decimal_value(token)
        if value is not None and value < 0:
            negative_underflow_indices.append(index)

    source_negative = (
        (numeric.is_finite() & (numeric < 0.0))
        | pl.col(row_column).is_in(negative_underflow_indices)
    ).fill_null(False)
    return indexed.with_columns(
        numeric.is_null().alias(SOURCE_MISSING_STATUS_COLUMN),
        (
            numeric.is_not_null() & ~numeric.is_finite()
        ).fill_null(False).alias(SOURCE_NONFINITE_STATUS_COLUMN),
        source_negative.alias(SOURCE_NEGATIVE_STATUS_COLUMN),
    ).drop(row_column)


def _decimal_value(token):
    if token is None:
        return None
    try:
        value = Decimal(str(token).strip())
    except (InvalidOperation, ValueError):
        return None
    return value if value.is_finite() else None


def _decimal_ln(value: Decimal) -> float:
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        ctx.Emin = _DECIMAL_EMIN
        ctx.Emax = _DECIMAL_EMAX
        return float(value.ln())


def _raw_text_from_neglog10(value: Decimal) -> str:
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        ctx.Emin = _DECIMAL_EMIN
        ctx.Emax = _DECIMAL_EMAX
        return format(ctx.power(Decimal(10), -value), "E")


def _attach_raw_source_provenance(df, pval_col, source_col):
    """Attach exact log-p and underflow provenance for a raw-p source."""
    row_column = "__pval_provenance_row"
    if row_column in df.columns:
        raise PValueTypeError(
            "Input uses reserved internal p-value provenance column %r. Rename "
            "that input column before harmonisation." % row_column
        )
    indexed = df.with_row_index(row_column)
    numeric = pl.col(pval_col).cast(pl.Float64, strict=False)
    # Only literal zero, Float64 underflow, and subnormal values need Decimal.
    # Ordinary GWAS rows stay entirely vectorised and do not acquire Python
    # objects, which keeps this step practical for tens of millions of rows.
    candidates = indexed.filter(
        numeric.is_not_null()
        & numeric.is_finite()
        & (numeric.abs() < sys.float_info.min)
    ).select(row_column, source_col, pval_col)

    reported_zero_indices = []
    underflow_indices = []
    negative_underflow_indices = []
    log_override_indices = []
    log_override_values = []
    for index, token, numeric_value in candidates.iter_rows():
        value = _decimal_value(token)
        if value is None:
            continue
        if value.is_zero():
            reported_zero_indices.append(index)
        elif value < 0:
            negative_underflow_indices.append(index)
        elif value <= 1:
            log_override_indices.append(index)
            log_override_values.append(_decimal_ln(value))
            if float(numeric_value) == 0.0:
                underflow_indices.append(index)

    if log_override_indices:
        overrides = pl.DataFrame({
            row_column: pl.Series(
                row_column,
                log_override_indices,
                dtype=indexed.schema[row_column],
            ),
            "__pval_log_override": log_override_values,
        })
        indexed = indexed.join(
            overrides, on=row_column, how="left", coalesce=True,
        )
        exact_log = pl.col("__pval_log_override")
    else:
        exact_log = pl.lit(None, dtype=pl.Float64)

    underflow = pl.col(row_column).is_in(underflow_indices)
    reported_zero = pl.col(row_column).is_in(reported_zero_indices)
    negative_underflow = pl.col(row_column).is_in(negative_underflow_indices)
    numeric_log = (
        pl.when(
            numeric.is_finite() & (numeric > 0.0) & (numeric <= 1.0)
        )
        .then(numeric.log())
        .otherwise(None)
    )
    result = indexed.with_columns(
        pl.when(underflow)
        .then(None)
        .when(negative_underflow)
        .then(pl.lit(-sys.float_info.min))
        .otherwise(numeric)
        .alias(pval_col),
        pl.coalesce(exact_log, numeric_log).alias(LOG_P_COLUMN),
        reported_zero.alias(REPORTED_ZERO_COLUMN),
        underflow.alias(FLOAT_UNDERFLOW_COLUMN),
        pl.when(underflow)
        .then(pl.col(source_col))
        .otherwise(None)
        .cast(pl.String)
        .alias(EXACT_RAW_P_COLUMN),
    )
    drop_columns = [row_column]
    if "__pval_log_override" in result.columns:
        drop_columns.append("__pval_log_override")
    return result.drop(drop_columns)


def _synchronise_log_from_raw_p(df, pval_col):
    """Make log provenance follow every representable validated raw value."""
    numeric = pl.col(pval_col).cast(pl.Float64, strict=False)
    existing = pl.col(LOG_P_COLUMN).cast(pl.Float64, strict=False)
    return df.with_columns(
        pl.when(existing.is_finite() & (existing <= 0.0))
        .then(existing)
        .when(numeric.is_finite() & (numeric > 0.0) & (numeric <= 1.0))
        .then(numeric.log())
        .otherwise(None)
        .alias(LOG_P_COLUMN)
    )


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
    raise PValueTypeError(
        "unknown p-value scale %r - expected 'raw' or 'neglog10'" % (value,)
    )


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
    ``{'raw', 'neglog10'}``. Raises ``AmbiguousPValueTypeError`` when the
    configured count, fraction, and median evidence do not agree. The optional
    ``statistics`` input lets the dataset-level resolver reuse the canonical
    aggregates from its combined scan.
    """
    policies = resolve_policies(policies)

    if not pval_col or pval_col not in df.columns:
        raise PValueTypeError(
            "P-value column %r is not in the data frame; the columns present "
            "are: %s" % (pval_col, ", ".join(df.columns))
        )

    threshold = policies.pvalue.mlogp_detect_threshold
    proportion = policies.pvalue.mlogp_detect_proportion
    minimum_count = policies.pvalue.mlogp_detect_min_count
    expected_median = policies.pvalue.mlogp_expected_median
    median_tolerance = policies.pvalue.mlogp_median_tolerance

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

    n_above_threshold = int(row["n_above_threshold"] or 0)
    n_negative = int(row["n_negative"] or 0)
    above_fraction = n_above_threshold / n_usable
    median = row["median"]
    maximum = row["max"]
    minimum = row["min"]
    median_low = expected_median - median_tolerance
    median_high = expected_median + median_tolerance

    evidence = {
        "pval_col": pval_col,
        "n_rows": n_rows,
        "n_non_null": int(row["n_non_null"] or 0),
        "n_usable": n_usable,
        "n_negative": n_negative,
        "n_above_threshold": n_above_threshold,
        "above_threshold_fraction": above_fraction,
        "detect_threshold": threshold,
        "detect_proportion": proportion,
        "detect_min_count": minimum_count,
        "min": minimum,
        "max": maximum,
        "mean": row["mean"],
        "median": median,
        "expected_median_neglog10": expected_median,
        "median_tolerance": median_tolerance,
        "median_acceptance_range": [median_low, median_high],
    }
    evidence.update(assess_negative_pvalues(row, pval_col, policies))

    configured = policies.pvalue.type
    if configured != "auto":
        evidence["decision_source"] = "config (pvalue.type)"
        return configured, evidence

    evidence["decision_source"] = "automatic detection"

    fraction_supports_log = above_fraction >= proportion
    count_supports_log = n_above_threshold >= minimum_count
    evidence["fraction_supports_neglog10"] = fraction_supports_log
    evidence["count_supports_neglog10"] = count_supports_log
    evidence["on_log_scale"] = fraction_supports_log and count_supports_log
    evidence["has_values_above_detection_threshold"] = n_above_threshold > 0

    if fraction_supports_log and not count_supports_log:
        raise AmbiguousPValueTypeError(
            "Column %r has %s usable values above the configured -log10 evidence "
            "threshold %s (%s), but pvalue.mlogp_detect_min_count requires at least "
            "%s. A single sentinel or corrupt value cannot select the scale for the "
            "study. Validate the outlying cells or set pvalue.type only after confirming "
            "the source representation. No chromosome processing was started."
            % (
                pval_col,
                _count(n_above_threshold),
                _fmt(threshold),
                "%.4f%%" % (above_fraction * 100.0),
                _count(minimum_count),
            )
        )

    if not evidence["on_log_scale"]:
        return "raw", evidence

    # Range evidence alone is insufficient: require the configured study-wide
    # median so -ln p and selected/top-hit-only files are never guessed.
    if median is None or median < median_low or median > median_high:
        raise AmbiguousPValueTypeError(
            "Column %r has -log10-like range evidence (%s values; %.4f%% exceed %s), "
            "but its study-wide median %s does not agree with the configured -log10 "
            "median %s within tolerance %s (accepted range %s to %s). This can indicate "
            "-ln(p), a selected/top-hit-only file, or mixed/invalid data. Convert the "
            "column to raw p or -log10 p, or explicitly set pvalue.type='neglog10' only "
            "when the source metadata proves that representation. No chromosome processing "
            "was started."
            % (
                pval_col,
                _count(n_above_threshold),
                above_fraction * 100.0,
                _fmt(threshold),
                _fmt(median),
                _fmt(expected_median),
                _fmt(median_tolerance),
                _fmt(median_low),
                _fmt(median_high),
            )
        )
    return "neglog10", evidence


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
        (usable < 0.0).sum().alias(prefix + "n_negative"),
        usable.min().alias(prefix + "min"),
        usable.max().alias(prefix + "max"),
        usable.mean().alias(prefix + "mean"),
        usable.median().alias(prefix + "median"),
    ]


def _detection_sentence(evidence):
    return (
        "%s of %s usable values (%s) exceed %s (median %s, maximum %s)"
        % (
            _count(evidence.get("n_above_threshold")),
            _count(evidence.get("n_usable")),
            "%.3f%%" % (evidence.get("above_threshold_fraction", 0.0) * 100.0),
            _fmt(evidence.get("detect_threshold")),
            _fmt(evidence.get("median"), "%.4f"),
            _fmt(evidence.get("max")),
        )
    )


# ---------------------------------------------------------------------------
# 3.  -log10 p  ->  raw p
# ---------------------------------------------------------------------------
def convert_negative_log10_to_p_value(
    df: pl.DataFrame,
    sample_column_dict,
    output_col: str,
    min_lp: Optional[float] = None,
    logger: Optional[Any] = None,
    policies: Optional[Any] = None,
    ctx: Optional[Any] = None,
):
    """Turn ``-log10 p`` back into plain p-values in ``output_col``.

    Negative LP values are refused because they imply p > 1; the orchestration
    path masks source-invalid rows for conversion, retains their source-status
    provenance, and applies the shared range policy afterward. Non-negative
    values below ``pvalue.mlogp_min`` follow that explicit floor policy.
    Positive LP values are never capped. When their raw probability is smaller
    than Float64 can represent, the public raw column is null and the exact raw
    token plus natural-log probability are retained in private provenance
    columns for calculations and export.
    """
    policies = resolve_policies(policies)
    ctx = ctx if ctx is not None else NullStepContext(df.height)
    if min_lp is None:
        min_lp = policies.pvalue.mlogp_min

    pval_col = sample_column_dict["pval_col"]
    if pval_col not in df.columns:
        raise PValueTypeError("Column '%s' not found in DataFrame." % (pval_col,))
    source_col = sample_column_dict.get(P_VALUE_SOURCE_TEXT_KEY)
    if not source_col or source_col not in df.columns:
        df, sample_column_dict = preserve_pvalue_source_text(
            df, sample_column_dict
        )
        source_col = sample_column_dict[P_VALUE_SOURCE_TEXT_KEY]

    values = pl.col(pval_col).cast(pl.Float64, strict=False)
    negative_mask = (values.is_finite() & (values < 0.0)).fill_null(False)
    low_mask = (
        values.is_finite() & (values >= 0.0) & (values < min_lp)
        if min_lp is not None
        else pl.lit(False)
    )
    counts = df.select(
        [
            negative_mask.sum().alias("n_negative"),
            low_mask.sum().alias("n_below"),
            values.max().alias("max_lp_seen"),
        ]
    ).to_dicts()[0]
    n_negative = int(counts["n_negative"] or 0)
    if n_negative:
        raise PValueTypeError(
            "Cannot convert %s negative value(s) in -log10 p column %r. "
            "Negative values are incompatible with -log10(p) and must be "
            "handled by the study-level pvalue.max_negative_fraction gate."
            % (_count(n_negative), pval_col)
        )
    n_below = int(counts["n_below"] or 0)

    adjusted_lp = (
        pl.when(low_mask)
        .then(pl.lit(min_lp))
        .otherwise(values)
        .alias("__lp_adjusted")
    )
    df = df.with_columns(adjusted_lp)
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
        negative_log10_to_raw_p("__lp_adjusted", output_name="__raw_p")
    )

    underflow_expr = (
        pl.col("__lp_adjusted").is_finite()
        & (pl.col("__lp_adjusted") >= 0.0)
        & (pl.col("__raw_p") == 0.0)
    ).fill_null(False)
    row_column = "__pval_neglog_row"
    if row_column in df.columns:
        raise PValueTypeError(
            "Input uses reserved internal p-value provenance column %r. Rename "
            "that input column before harmonisation." % row_column
        )
    df = df.with_row_index(row_column)
    underflow_rows = df.filter(underflow_expr).select(
        row_column, source_col, "__lp_adjusted",
    )
    exact_indices = []
    exact_values = []
    for index, source_value, adjusted_value in underflow_rows.iter_rows():
        original = _decimal_value(source_value)
        if original is None or (
            min_lp is not None and original < Decimal(str(min_lp))
        ):
            original = _decimal_value(adjusted_value)
        if original is not None:
            exact_indices.append(index)
            exact_values.append(_raw_text_from_neglog10(original))
    n_underflow = underflow_rows.height
    if exact_indices:
        exact_overrides = pl.DataFrame({
            row_column: pl.Series(
                row_column, exact_indices, dtype=df.schema[row_column],
            ),
            "__pval_exact_override": pl.Series(
                "__pval_exact_override", exact_values, dtype=pl.String,
            ),
        })
        df = df.join(
            exact_overrides, on=row_column, how="left", coalesce=True,
        )
        exact_raw = pl.col("__pval_exact_override")
    else:
        exact_raw = pl.lit(None, dtype=pl.String)

    df = df.with_columns(
        pl.when(~pl.col("__lp_adjusted").is_finite() | underflow_expr)
        .then(None)
        .otherwise(pl.col("__raw_p"))
        .alias(output_col),
        low_mask.fill_null(False).alias(CLIPPED_HIGH_COLUMN),
        pl.lit(False).alias(CLIPPED_LOW_COLUMN),
        underflow_expr.alias(FLOAT_UNDERFLOW_COLUMN),
        pl.when(pl.col("__lp_adjusted").is_finite())
        .then(-pl.col("__lp_adjusted") * math.log(10.0))
        .otherwise(None)
        .alias(LOG_P_COLUMN),
        exact_raw.alias(EXACT_RAW_P_COLUMN),
        pl.lit(False).alias(REPORTED_ZERO_COLUMN),
    )
    df = df.drop([
        column
        for column in (
            "__lp_adjusted", "__raw_p", row_column, "__pval_exact_override",
        )
        if column in df.columns
    ])

    mlogp_stats = df.select(
        pl.col(output_col).min().alias("min_PVAL"),
        pl.col(output_col).max().alias("max_PVAL"),
        pl.col(output_col).mean().alias("mean_PVAL"),
        pl.col(output_col).median().alias("median_PVAL"),
    ).to_dicts()[0]

    mlogp_stats["mlogp_values_total_clipped"] = n_below
    mlogp_stats["mlogp_values_clipped_at_min"] = n_below
    mlogp_stats["mlogp_max_observed"] = counts["max_lp_seen"]
    mlogp_stats["pvalue_float_underflow"] = n_underflow
    ctx.info(
        "Converted -log10 p to plain p-values in '%s': minimum %s, maximum %s, "
        "median %s. %s values were smaller than positive Float64 and retain "
        "their exact raw text and log probability without clipping."
        % (
            output_col,
            _fmt(mlogp_stats["min_PVAL"], "%.3e"),
            _fmt(mlogp_stats["max_PVAL"], "%.3e"),
            _fmt(mlogp_stats["median_PVAL"], "%.3e"),
            _count(n_underflow),
        )
    )

    sample_column_dict["pval_col"] = output_col
    sample_column_dict[P_VALUE_LOG_KEY] = LOG_P_COLUMN
    sample_column_dict[P_VALUE_EXACT_RAW_KEY] = EXACT_RAW_P_COLUMN
    return df, mlogp_stats, sample_column_dict


# ---------------------------------------------------------------------------
# 5.  The step itself
# ---------------------------------------------------------------------------
def _harmonise_canonical_pvalues(
    df, pval_col, input_scale, policies, ctx, rejects, qc_info, chromosome,
):
    """Validate one canonical raw-p column after either input-scale branch.

    Three gates, in this order:

    1. source/conversion missing without exact provenance -> ``pval_null``
    2. source non-finite, scale-invalid, or impossible raw value -> range policy
    3. literal raw zero or rounding excess above ``clip_high`` -> range policy

    Source-status flags are materialised before conversion, so missing and
    non-finite -log10 cells cannot collapse into the same anonymous null.
    Positive probabilities below Float64 retain exact raw/log provenance and
    are explicitly exempt from the missing gate.
    """
    tolerance = policies.pvalue.tolerance_above_one
    clip_low = policies.pvalue.clip_low
    clip_high = policies.pvalue.clip_high
    action = policies.pvalue.out_of_range

    missing_clip_columns = [
        name for name in CLIPPED_COLUMNS if name not in df.columns
    ]
    if missing_clip_columns:
        df = df.with_columns(
            [pl.lit(False).alias(name) for name in missing_clip_columns]
        )

    exact_underflow = (
        pl.col(FLOAT_UNDERFLOW_COLUMN).fill_null(False)
        & pl.col(LOG_P_COLUMN).is_finite()
        & pl.col(EXACT_RAW_P_COLUMN).is_not_null()
    )
    values = pl.col(pval_col).cast(pl.Float64, strict=False)
    source_missing = pl.col(SOURCE_MISSING_STATUS_COLUMN).fill_null(False)
    source_nonfinite = pl.col(SOURCE_NONFINITE_STATUS_COLUMN).fill_null(False)
    source_negative = pl.col(SOURCE_NEGATIVE_STATUS_COLUMN).fill_null(False)
    scale_invalid = (
        source_negative if input_scale == "neglog10" else pl.lit(False)
    )
    source_invalid = source_nonfinite | scale_invalid
    missing = (
        source_missing
        | (values.is_null() & ~exact_underflow & ~source_invalid)
    ).fill_null(False)
    below_zero = (values.is_finite() & (values < 0.0)).fill_null(False)
    above_tolerance = (
        values.is_finite() & (values > tolerance)
    ).fill_null(False)
    unsalvageable = (source_invalid | below_zero | above_tolerance).fill_null(False)
    reported_zero = pl.col(REPORTED_ZERO_COLUMN).fill_null(False)
    too_high = (
        values.is_finite()
        & (values > clip_high)
        & (values <= tolerance)
    ).fill_null(False)
    counts = df.select(
        [
            missing.sum().alias("n_null"),
            source_missing.sum().alias("n_source_missing"),
            source_nonfinite.sum().alias("n_source_nonfinite"),
            source_negative.sum().alias("n_source_negative"),
            below_zero.sum().alias("n_lt0"),
            above_tolerance.sum().alias("n_gt_tolerance"),
            unsalvageable.sum().alias("n_unsalvageable"),
            reported_zero.sum().alias("n_below_clip"),
            too_high.sum().alias("n_above_clip"),
        ]
    ).to_dicts()[0]
    n_null = int(counts["n_null"] or 0)
    n_source_missing = int(counts["n_source_missing"] or 0)
    n_source_nonfinite = int(counts["n_source_nonfinite"] or 0)
    n_source_negative = int(counts["n_source_negative"] or 0)
    n_lt0 = int(counts["n_lt0"] or 0)
    n_gt_tolerance = int(counts["n_gt_tolerance"] or 0)
    n_unsalvageable = int(counts["n_unsalvageable"] or 0)
    n_below_clip = int(counts["n_below_clip"] or 0)
    n_above_clip = int(counts["n_above_clip"] or 0)
    n_out_of_range = n_below_clip + n_above_clip

    total_before = df.height

    if action == "fail" and (n_unsalvageable or n_out_of_range):
        raise PValueTypeError(
            "Chromosome %s: %s %s-source or canonical p-values are unusable "
            "(%s non-finite source values, %s negative source values, %s "
            "canonical values below 0, %s above %s), and %s more are literal "
            "zero or above %s. Policy "
            "'pvalue.out_of_range' is 'fail'."
            % (
                chromosome,
                _count(n_unsalvageable),
                input_scale,
                _count(n_source_nonfinite),
                _count(n_source_negative),
                _count(n_lt0),
                _count(n_gt_tolerance),
                _fmt(tolerance),
                _count(n_out_of_range),
                _fmt(clip_high),
            )
        )

    # 1. missing
    if n_null:
        df, _ = reject_rows(
            df,
            missing,
            step_label=STEP_LABEL,
            reason="pval_null",
            context=ctx,
            collector=rejects,
            detail=None,
            description=(
                "%s variants have a missing or unparseable source p-value, or "
                "conversion produced no canonical probability and no exact "
                "log-probability provenance." % (_count(n_null),)
            ),
            check_name="p-value present",
            warn_on_remove=True,
        )

    # 2. impossible values. Under 'clip' they are rejected rather than forced
    # into range: no finite correction can recover a non-finite source value or
    # a negative -log10(p) without inventing p=1.
    if n_unsalvageable:
        plain = (
            "%s variants have an unusable %s source or canonical p-value (%s "
            "non-finite source values, %s negative source values, %s canonical "
            "values below 0, %s above %s), which no rounding correction can "
            "repair."
            % (
                _count(n_unsalvageable),
                input_scale,
                _count(n_source_nonfinite),
                _count(n_source_negative),
                _count(n_lt0),
                _count(n_gt_tolerance),
                _fmt(tolerance),
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
            df, _ = reject_rows(
                df,
                unsalvageable,
                step_label=STEP_LABEL,
                reason="pval_out_of_range",
                context=ctx,
                collector=rejects,
                detail="non-finite, scale-invalid, or outside 0 to %s" % tolerance,
                description=plain,
                check_name="p-value range",
                warn_on_remove=True,
            )

    # 3. literal zero or a rounding excess above one. Positive p-values are
    # valid at every representable magnitude and are not floored.
    too_low = pl.col(REPORTED_ZERO_COLUMN).fill_null(False)
    if action == "reject" and n_out_of_range:
        df, _ = reject_rows(
            df,
            too_low | too_high,
            step_label=STEP_LABEL,
            reason="pval_out_of_range",
            context=ctx,
            collector=rejects,
            detail="literal zero or above %s" % clip_high,
            description=(
                "%s variants have a literal-zero p-value or a value above %s (%s "
                "zero, %s above 1)."
                % (
                    _count(n_out_of_range),
                    _fmt(clip_high),
                    _count(n_below_clip),
                    _count(n_above_clip),
                )
            ),
            check_name="p-value range",
            warn_on_remove=True,
        )
    elif action == "null" and n_out_of_range:
        before = df.height
        df = df.with_columns(
            pl.when(too_low | too_high)
            .then(None)
            .otherwise(pl.col(pval_col))
            .alias(pval_col)
        )
        ctx.qc(
            "p-value range",
            "%s variants have a literal-zero p-value or a value above %s. Policy "
            "'pvalue.out_of_range' is 'null', so their p-value was blanked out "
            "and the variants kept."
            % (_count(n_out_of_range), _fmt(clip_high)),
            before,
            df.height,
            changed=n_out_of_range,
            warn=True,
            step=STEP_LABEL,
        )
    else:
        # 'clip' supplies an explicit approximation only for literal zero and
        # corrects small rounding excesses above one.
        before = df.height
        df = df.with_columns(
            [
                (
                    pl.col(CLIPPED_HIGH_COLUMN).fill_null(False) | too_high
                ).alias(CLIPPED_HIGH_COLUMN),
                (
                    pl.col(CLIPPED_LOW_COLUMN).fill_null(False) | too_low
                ).alias(CLIPPED_LOW_COLUMN),
            ]
        )
        df = df.with_columns(
            pl.when(too_low)
            .then(pl.lit(clip_low))
            .when(too_high)
            .then(pl.lit(clip_high))
            .otherwise(pl.col(pval_col))
            .alias(pval_col)
        )
        if n_out_of_range:
            ctx.qc(
                "p-value range",
                "%s variants have a literal-zero p-value or a value above %s "
                "(%s zero, %s above 1). Policy 'pvalue.out_of_range' is 'clip', "
                "so zero was approximated by %s and rounding excesses were set "
                "to %s."
                % (
                    _count(n_out_of_range),
                    _fmt(clip_high),
                    _count(n_below_clip),
                    _count(n_above_clip),
                    _fmt(clip_low),
                    _fmt(clip_high),
                ),
                before,
                df.height,
                changed=n_out_of_range,
                warn=True,
                step=STEP_LABEL,
            )

    df = _synchronise_log_from_raw_p(df, pval_col)

    qc_info.update(
        {
            "pvalue_validation_input_scale": input_scale,
            "pvalue_validation_initial_variants": total_before,
            "after_filter_variants_with_pvalues": df.height,
            "variants_with_null_pvalues_removed": n_null,
            "pvalue_source_missing_or_unparseable": n_source_missing,
            "pvalue_source_non_finite": n_source_nonfinite,
            "pvalue_source_negative": n_source_negative,
            "variants_with_lt0_pvalues_removed": (
                n_source_negative if action in ("clip", "reject") else 0
            ),
            "variants_with_non_finite_pvalues_removed": (
                n_source_nonfinite if action in ("clip", "reject") else 0
            ),
            "variants_with_gt1_05_pvalues_removed": (
                n_gt_tolerance if action in ("clip", "reject") else 0
            ),
            "variants_with_zero_pvalues_replaced": n_below_clip,
            "variants_with_pvalues_clipped_low": n_below_clip,
            "variants_with_pvalues_clipped_high": n_above_clip,
            "pvalue_out_of_range_action": action,
            "pvalue_clip_low": clip_low,
            "pvalue_clip_high": clip_high,
            "pvalue_tolerance_above_one": tolerance,
            "pvalue_float_underflow": int(
                df.select(pl.col(FLOAT_UNDERFLOW_COLUMN).sum()).item() or 0
            ),
        }
    )
    if input_scale == "raw":
        qc_info["raw_total_variants_with_pvalues"] = total_before
    else:
        qc_info.update(
            {
                "mlogp_missing_or_unparseable": n_source_missing,
                "mlogp_non_finite": n_source_nonfinite,
                "mlogp_negative": n_source_negative,
            }
        )
    return df, qc_info


def harmonise_p_values(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict,
    proportion_threshold: Optional[float] = None,
    logger: Optional[Any] = None,
    policies: Optional[Any] = None,
    rejects: Optional[Any] = None,
    decision: Optional[Any] = None,
    step_number: Optional[int] = None,
    step_total: Optional[int] = None,
):
    """Detect the p-value scale (or apply the study-level decision) and leave a
    plain p-value column in the configured ``pvalue.output_column``.

    ``proportion_threshold`` overrides the canonical detection fraction only
    for this call. Positive -log10 values are never capped.

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
    output_col = policies.pvalue.output_column
    if output_col in df.columns and output_col != pval_col:
        raise PValueTypeError(
            "Configured canonical raw p-value output column %r "
            "(pvalue.output_column) already exists, but the selected input "
            "p-value column is %r. Refusing to overwrite unrelated data; rename "
            "the existing column or configure a distinct pvalue.output_column."
            % (output_col, pval_col)
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
        df, sample_column_dict = preserve_pvalue_source_text(
            df, sample_column_dict
        )
        source_col = sample_column_dict[P_VALUE_SOURCE_TEXT_KEY]
        df = df.with_columns(pl.col(pval_col).cast(pl.Float64, strict=False))
        df = _attach_pvalue_source_status(df, pval_col, source_col)

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
                pl.col(pval_col).min().alias("initial_min_pvalues"),
                pl.col(pval_col).max().alias("initial_max_pvalues"),
                pl.col(pval_col).mean().alias("initial_mean_pvalues"),
                pl.col(pval_col).median().alias("initial_median_pvalues"),
            ]
        ).to_dicts()[0]
        qc_info = initial_stats.copy()
        qc_info["pvalue_column"] = pval_col
        qc_info["pvalue_input_column"] = pval_col
        qc_info["pvalue_output_column"] = output_col

        ctx.info(
            "Numeric p-value-column values on entry, before scale and source-token "
            "interpretation: %s missing, %s numeric zero, %s below zero, %s above "
            "one; minimum %s, median %s, maximum %s. A numeric zero can still be a "
            "positive raw source token below Float64."
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
                    "pvalue_above_detection_threshold_count": evidence.get(
                        "n_above_threshold"
                    ),
                    "pvalue_above_detection_threshold_fraction": evidence.get(
                        "above_threshold_fraction"
                    ),
                    "pvalue_detection_minimum_count": evidence.get("detect_min_count"),
                    "pvalue_expected_neglog10_median": evidence.get(
                        "expected_median_neglog10"
                    ),
                    "pvalue_neglog10_median_tolerance": evidence.get(
                        "median_tolerance"
                    ),
                    "median": evidence.get("median"),
                    "max": evidence.get("max"),
                }
            )

        ctx.decide(
            "P-value column '%s'" % (pval_col,),
            _detection_sentence(evidence) if evidence is not None else decision_source,
            "%s (%s)" % (_SCALE_TEXT[scale], decision_source),
        )

        if (
            scale == "raw"
            and evidence is not None
            and evidence.get("has_values_above_detection_threshold")
        ):
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
            df = _attach_raw_source_provenance(df, pval_col, source_col)
            source_counts = df.select(
                pl.col(REPORTED_ZERO_COLUMN).sum().alias("reported_zero"),
                pl.col(FLOAT_UNDERFLOW_COLUMN).sum().alias("underflow"),
            ).row(0, named=True)
            qc_info["initial_variants_with_zero_pvalues"] = int(
                source_counts["reported_zero"] or 0
            )
            qc_info["initial_positive_pvalues_below_float64"] = int(
                source_counts["underflow"] or 0
            )
            if output_col != pval_col:
                df = df.with_columns(pl.col(pval_col).alias(output_col))
            sample_column_dict["pval_col"] = output_col
        else:
            total_before = df.height
            # The low-level converter deliberately refuses negative LP values.
            # Mask source-invalid cells for conversion while retaining their
            # materialised status so the shared post-conversion gate can apply
            # pvalue.out_of_range without ever turning a negative LP into p=1.
            converter_invalid = (
                pl.col(SOURCE_NONFINITE_STATUS_COLUMN).fill_null(False)
                | pl.col(SOURCE_NEGATIVE_STATUS_COLUMN).fill_null(False)
            )
            df = df.with_columns(
                pl.when(converter_invalid)
                .then(None)
                .otherwise(pl.col(pval_col))
                .alias(pval_col)
            )
            df, mlogp_stats, sample_column_dict = convert_negative_log10_to_p_value(
                df=df,
                sample_column_dict=sample_column_dict,
                output_col=output_col,
                min_lp=policies.pvalue.mlogp_min,
                policies=policies,
                ctx=ctx,
            )
            qc_info.update(
                {
                    "mlogp_total_variants_with_pvalues": total_before,
                    "after_mlogp_to_pval_variants_with_pvalues": df.height,
                }
            )
            qc_info.update(mlogp_stats)

        pval_col = sample_column_dict["pval_col"]
        df, qc_info = _harmonise_canonical_pvalues(
            df,
            pval_col,
            scale,
            policies,
            ctx,
            rejects,
            qc_info,
            chromosome,
        )

        sample_column_dict[P_VALUE_LOG_KEY] = LOG_P_COLUMN
        sample_column_dict[P_VALUE_EXACT_RAW_KEY] = EXACT_RAW_P_COLUMN

        # Both branches must leave the flags behind for the effect-statistics
        # step, even when nothing was clipped.
        for column in (
            CLIPPED_HIGH_COLUMN,
            CLIPPED_LOW_COLUMN,
            REPORTED_ZERO_COLUMN,
            FLOAT_UNDERFLOW_COLUMN,
        ):
            if column not in df.columns:
                df = df.with_columns(pl.lit(False).alias(column))
        for column in (SOURCE_TEXT_COLUMN, EXACT_RAW_P_COLUMN):
            if column not in df.columns:
                df = df.with_columns(pl.lit(None, dtype=pl.String).alias(column))
        if LOG_P_COLUMN not in df.columns:
            df = df.with_columns(
                pl.lit(None, dtype=pl.Float64).alias(LOG_P_COLUMN)
            )

        final_stats = df.select(
            [
                pl.len().alias("final_total_variants_with_pvalues"),
                (
                    pl.col(pval_col).is_null()
                    & ~(
                        pl.col(FLOAT_UNDERFLOW_COLUMN).fill_null(False)
                        & pl.col(LOG_P_COLUMN).is_finite()
                        & pl.col(EXACT_RAW_P_COLUMN).is_not_null()
                    )
                )
                .sum()
                .alias("final_variants_with_null_pvalues"),
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
                pl.col(FLOAT_UNDERFLOW_COLUMN)
                .sum()
                .alias("pvalues_preserved_below_float64"),
            ]
        ).to_dicts()[0]
        qc_info.update(final_stats)

        df = df.drop(
            [name for name in SOURCE_STATUS_COLUMNS if name in df.columns]
        )

        ctx.set_rows(df.height)

    return df, qc_info, sample_column_dict
