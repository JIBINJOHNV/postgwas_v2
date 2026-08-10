"""Step 08 (per chromosome) — standard error from BETA and the p-value.

Public entry point: ``derive_standard_error_from_effect_and_p_value``.

    SE = |BETA / z|,   z = norm.isf(p / tail)

**The clamp can fabricate standard errors that look healthy.** A value above 1
may be clamped just below 1, producing a very small Z score and an enormous SE.
A raw p-value reported as exactly zero is a different form of information loss:
it is a censored or underflow value whose true magnitude is unknown. The default
``pvalue.zero_missing_se = fail`` therefore stops only when SE is still missing
after the Z-based recovery step. An explicit ``approximate`` override uses
``pvalue.clip_low`` and marks the affected rows; ``reject`` records them in the
rejected-variant output. A supplied SE or an SE recovered from a usable Z is
preserved and does not trigger this policy.

GWAS Catalog requires accompanying precision metadata when a submitted p-value
is zero: https://www.ebi.ac.uk/gwas/docs/summary-statistics-format

The column is always present on the returned frame — including when the study
already supplied a standard error and nothing was fabricated — so the
validation step can read it unconditionally.

Logging and rejected variants
-----------------------------
Pass ``logger`` (a ``PipelineLogger``), ``policies`` (a ``Policies``) and
``rejects`` (a ``RejectCollector``); all three are optional.  ``POLICY_KEYS``
lists every policy this step reads.
"""

import numpy as np
import polars as pl
from postgwas.core.values import optional_text

from functools import partial
from scipy.stats import norm
from typing import Any, Dict, Optional, Tuple

from postgwas.core.values import format_number

from .effect_from_z import SE_UNAVAILABLE_FROM_Z_COLUMN
from .p_values import (
    CLIPPED_HIGH_COLUMN,
    CLIPPED_LOW_COLUMN,
    REPORTED_ZERO_COLUMN,
)
from .shared.runtime import log_info, log_warning, reject_rows, resolve_policies


__all__ = [
    "derive_standard_error_from_effect_and_p_value",
    "SE_FROM_CLIPPED_PVAL_COL",
    "SE_FROM_ZERO_P_APPROXIMATION_COL",
    "POLICY_KEYS",
    "STEP_LABEL",
]


#: Internal boolean column marking a standard error derived from a clamped
#: p-value.  Read by the effect-statistics validation step.
SE_FROM_CLIPPED_PVAL_COL = "__se_from_clipped_pval"

#: Internal boolean provenance for the explicit approximation policy. Unlike
#: the generic clipped-p-value marker, this tag is not rejected later because
#: the user deliberately selected the approximation. It remains in QC/logs and
#: is removed with all other internal columns before export.
SE_FROM_ZERO_P_APPROXIMATION_COL = "__se_from_zero_p_approximation"

#: Every policy this step reads.  Pass to ``logger.step(policy_keys=...)``.
POLICY_KEYS = (
    "pvalue.se_tail",
    "pvalue.out_of_range",
    "pvalue.clip_low",
    "pvalue.zero_missing_se",
    "validation.se_from_clipped_pval",
)

#: The free-text step label recorded in the reject file.
STEP_LABEL = "08 se_from_beta_pval"

#: Numerical guard used only for invalid values above 1. A valid p-value equal
#: to 1 is not clamped; when SE must be inferred it naturally yields z=0 and an
#: undefined SE, which is handled as missing rather than fabricated.
_PVAL_CLIP_HIGH = 0.999999999

# Shared mechanics stay outside the scientific implementation.
_info = log_info
_warn = log_warning
_fmt = partial(format_number, pattern="%.6f", missing="n/a")


def _flag(df: pl.DataFrame, name: str, default: bool = False) -> pl.Expr:
    """Return a null-safe boolean expression for an optional working column."""
    if name not in df.columns:
        return pl.lit(default)
    value = pl.col(name)
    if df.schema[name] != pl.Boolean:
        value = value.cast(pl.Boolean, strict=False)
    return value.fill_null(default)


def _se_statistics(df: pl.DataFrame, column: str) -> Dict[str, Any]:
    return df.select(
        [
            pl.col(column).min().alias("min"),
            pl.col(column).max().alias("max"),
            pl.col(column).mean().alias("mean"),
            pl.col(column).median().alias("median"),
            pl.len().alias("total"),
        ]
    ).to_dicts()[0]


# public entry point
# ---------------------------------------------------------------------------
def derive_standard_error_from_effect_and_p_value(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    logger: Optional[Any] = None,
    policies: Optional[Any] = None,
    rejects: Optional[Any] = None,
) -> Tuple[pl.DataFrame, Dict, dict]:
    """Derive SE from BETA and the p-value when the study supplies no SE."""
    pol = resolve_policies(policies)

    qc_info = {"initial_variants": df.height}  # type: Dict[str, Any]

    if logger is not None:
        logger.settings(list(POLICY_KEYS))

    beta_col = optional_text(sample_column_dict.get("beta_col"))
    se_col = optional_text(sample_column_dict.get("se_col"))
    pval_col = optional_text(sample_column_dict.get("pval_col"))

    # -------------------------------------------------------
    # VALIDATION
    # -------------------------------------------------------
    if (
        beta_col is None
        or pval_col is None
        or beta_col not in df.columns
        or pval_col not in df.columns
    ):
        raise ValueError(
            "Chromosome {}: the standard error must be derived from the effect size "
            "and the p-value, but the effect size column ('{}') and the p-value column "
            "('{}') are not both present.".format(chromosome, beta_col, pval_col)
        )

    clip_low = float(pol.get("pvalue.clip_low"))
    zero_missing_se_action = str(pol.get("pvalue.zero_missing_se"))
    has_se = se_col is not None and se_col in df.columns

    # A raw zero is recorded by the p-value step before it raises the value to
    # pvalue.clip_low. Standalone callers that have not run that step yet still
    # receive the same behaviour by testing the raw p-value directly.
    reported_zero = (
        _flag(df, REPORTED_ZERO_COLUMN)
        if REPORTED_ZERO_COLUMN in df.columns
        else (pl.col(pval_col).cast(pl.Float64, strict=False) == 0).fill_null(False)
    )
    needs_se = _flag(df, SE_UNAVAILABLE_FROM_Z_COLUMN) if has_se else pl.lit(True)
    zero_missing_se = reported_zero & needs_se
    n_zero_missing_se = int(
        df.select(zero_missing_se.sum().alias("count")).item() or 0
    )
    n_zero_rejected = 0
    n_zero_approximated = 0

    if n_zero_missing_se and zero_missing_se_action == "fail":
        raise ValueError(
            "Chromosome {}: {:,} highly significant variant(s) have a reported "
            "p-value of 0, but neither a supplied standard error nor a usable Z "
            "score is available. A zero p-value is a censored or underflow value, "
            "not an exact probability, so PostGWAS cannot reconstruct a "
            "scientifically valid SE. No variants were silently discarded. Provide "
            "SE, provide Z, supply negative-log10 p-values or precision metadata, or "
            "explicitly choose approximation/rejection in the harmonisation "
            "configuration. To override this run, use '--zero-p-se-action "
            "approximate' or '--zero-p-se-action reject'.".format(
                chromosome, n_zero_missing_se
            )
        )

    if n_zero_missing_se and zero_missing_se_action == "reject":
        df, n_zero_rejected = reject_rows(
            df,
            zero_missing_se,
            step_label=STEP_LABEL,
            reason="zero_p_missing_se",
            check_name="zero p-value without SE or usable Z",
            description=(
                "The reported p-value is zero and neither a supplied standard error "
                "nor a usable Z score is available, so an exact standard error "
                "cannot be reconstructed."
            ),
            collector=rejects,
            context=logger,
            detail="policy pvalue.zero_missing_se=reject",
        )
        reported_zero = (
            _flag(df, REPORTED_ZERO_COLUMN)
            if REPORTED_ZERO_COLUMN in df.columns
            else (pl.col(pval_col).cast(pl.Float64, strict=False) == 0).fill_null(False)
        )
        needs_se = _flag(df, SE_UNAVAILABLE_FROM_Z_COLUMN) if has_se else pl.lit(True)
        zero_missing_se = reported_zero & needs_se

    # -------------------------------------------------------
    # SKIP WHEN EVERY ROW ALREADY HAS AN SE
    # -------------------------------------------------------
    if has_se:
        if n_zero_missing_se and zero_missing_se_action == "approximate":
            df = df.with_columns(
                [
                    pl.col(beta_col).cast(pl.Float64, strict=False),
                    pl.col(pval_col).cast(pl.Float64, strict=False),
                    pl.col(se_col).cast(pl.Float64, strict=False),
                ]
            )
            p_safe = (
                pl.when(zero_missing_se)
                .then(pl.lit(clip_low))
                .when(pl.col(pval_col) <= 0)
                .then(pl.lit(clip_low))
                .otherwise(pl.col(pval_col))
            )
            p_values = df.select(p_safe.alias("p")).get_column("p").to_numpy()
            with np.errstate(divide="ignore", invalid="ignore"):
                candidate = np.abs(
                    df[beta_col].to_numpy()
                    / norm.isf(p_values / int(pol.get("pvalue.se_tail")))
                )
            df = df.with_columns(
                pl.Series("__zero_p_approximate_se", candidate, dtype=pl.Float64)
            ).with_columns(
                pl.when(zero_missing_se)
                .then(pl.col("__zero_p_approximate_se"))
                .otherwise(pl.col(se_col))
                .alias(se_col)
            ).drop("__zero_p_approximate_se")
            n_zero_approximated = n_zero_missing_se

        df = df.with_columns(
            [
                pl.lit(False).alias(SE_FROM_CLIPPED_PVAL_COL),
                (
                    zero_missing_se
                    if zero_missing_se_action == "approximate"
                    else pl.lit(False)
                ).alias(SE_FROM_ZERO_P_APPROXIMATION_COL),
            ]
        )
        source = (
            "The study supplied no standard error, but one was recovered as BETA / Z "
            "where Z was usable"
            if SE_UNAVAILABLE_FROM_Z_COLUMN in df.columns
            else "The study supplies a standard error in column '{}'".format(se_col)
        )
        _info(
            logger,
            "{}, so no other standard error is derived from the p-value.".format(source),
        )
        if n_zero_approximated:
            _warn(
                logger,
                "{:,} variant(s) reported p=0 and had no usable Z-derived SE. The "
                "explicit approximation policy used pvalue.clip_low={} and marked "
                "them in '{}'. Their true significance and exact SE remain unknown."
                .format(
                    n_zero_approximated,
                    clip_low,
                    SE_FROM_ZERO_P_APPROXIMATION_COL,
                ),
            )
            if logger is not None:
                logger.qc(
                    "zero p-value SE approximation",
                    "The raw p-value was exactly zero and no supplied or usable "
                    "Z-derived SE was available. The explicit approximation policy "
                    "used pvalue.clip_low; the exact SE remains unknown.",
                    df.height,
                    df.height,
                    changed=n_zero_approximated,
                    warn=True,
                    step=STEP_LABEL,
                )

        stats = _se_statistics(df, se_col)
        n_undefined = int(
            df.select(
                (pl.col(se_col).is_null() | ~pl.col(se_col).is_finite()).sum()
            ).item()
            or 0
        )
        qc_info.update(
            {
                "status": "SE already present"
                if not n_zero_approximated
                else "missing Z-derived SE approximated from zero p-value floor",
                "se_min": stats["min"],
                "se_max": stats["max"],
                "se_mean": stats["mean"],
                "se_median": stats["median"],
                "total_variants": stats["total"],
                "se_from_clipped_pval": 0,
                "se_from_zero_p_approximation": n_zero_approximated,
                "zero_p_missing_se": n_zero_missing_se,
                "zero_p_missing_se_action": zero_missing_se_action,
                "zero_p_missing_se_rejected": n_zero_rejected,
                "se_undefined": n_undefined,
                "variants_removed_total": qc_info["initial_variants"] - df.height,
                "variants_remaining": df.height,
            }
        )
        return df, qc_info, sample_column_dict

    # -------------------------------------------------------
    # CAST
    # -------------------------------------------------------
    df = df.with_columns(
        [
            pl.col(beta_col).cast(pl.Float64, strict=False),
            pl.col(pval_col).cast(pl.Float64, strict=False),
        ]
    )

    se_tail = int(pol.get("pvalue.se_tail"))

    out_of_range_action = pol.get("pvalue.out_of_range")
    n_initial = qc_info["initial_variants"]
    # -------------------------------------------------------
    # QC on the p-value before anything is clamped
    # -------------------------------------------------------
    counts = df.select(
        [
            (pl.col(pval_col) == 0).fill_null(False).sum().alias("zero"),
            (pl.col(pval_col) < 0).fill_null(False).sum().alias("negative"),
            (pl.col(pval_col) > 1).fill_null(False).sum().alias("above_one"),
            pl.col(pval_col).is_null().sum().alias("missing"),
        ]
    ).to_dicts()[0]

    zero_pval = int(counts["zero"] or 0)
    negative_pval = int(counts["negative"] or 0)
    high_pval = int(counts["above_one"] or 0)
    null_pval = int(counts["missing"] or 0)
    # An explicit zero approximation owns this one otherwise-invalid value;
    # the generic range action continues to govern negative and high values.
    invalid_zero_pval = (
        0 if zero_missing_se_action == "approximate" else zero_pval
    )
    invalid_pval = invalid_zero_pval + negative_pval + high_pval

    qc_info.update(
        {
            "zero_pval": zero_pval,
            "negative_pval": negative_pval,
            "high_pval": high_pval,
            "null_pval": null_pval,
            "invalid_pval": invalid_pval,
            "se_tail": se_tail,
        }
    )

    if zero_pval:
        _warn(
            logger,
            "{:,} variants report a p-value of exactly 0. The true significance is "
            "unknown and below the floating-point floor, so the standard error derived "
            "for them corresponds to a Z score of about 37 rather than their real "
            "value.".format(zero_pval),
        )
    if negative_pval:
        _warn(
            logger,
            "{:,} variants report a negative p-value, which is impossible.".format(
                negative_pval
            ),
        )
    if null_pval:
        _warn(
            logger,
            "{:,} variants have no p-value, so no standard error can be derived for "
            "them and their SE is left empty.".format(null_pval),
        )

    # -------------------------------------------------------
    # OUT-OF-RANGE P-VALUES
    # -------------------------------------------------------
    zero_approximation = (
        reported_zero
        if zero_missing_se_action == "approximate"
        else pl.lit(False)
    )
    cond_low = (pl.col(pval_col) <= 0) & ~zero_approximation
    cond_high = pl.col(pval_col) > 1
    n_pval_removed = 0

    if out_of_range_action == "fail":
        if invalid_pval:
            raise ValueError(
                "Chromosome {}: {:,} variants have a p-value outside (0, 1] in column "
                "'{}' and policy pvalue.out_of_range is 'fail'.".format(
                    chromosome, invalid_pval, pval_col
                )
            )
        df = df.with_columns(
            [
                pl.lit(False).alias(SE_FROM_CLIPPED_PVAL_COL),
                pl.col(pval_col).alias("_pval_safe"),
            ]
        )
    elif out_of_range_action == "reject":
        df, n_pval_removed = reject_rows(
            df,
            cond_low | cond_high,
            step_label=STEP_LABEL,
            reason="pval_out_of_range",
            check_name="p-value inside (0, 1]",
            description=(
                "The p-value is not inside (0, 1], so no standard error can "
                "honestly be derived from it."
            ),
            collector=rejects,
            context=logger,
            detail="column '{}'".format(pval_col),
        )
        df = df.with_columns(
            [
                pl.lit(False).alias(SE_FROM_CLIPPED_PVAL_COL),
                pl.col(pval_col).alias("_pval_safe"),
            ]
        )
    elif out_of_range_action == "null":
        df = df.with_columns(
            [
                pl.lit(False).alias(SE_FROM_CLIPPED_PVAL_COL),
                pl.when(cond_low | cond_high)
                .then(None)
                .otherwise(pl.col(pval_col))
                .alias("_pval_safe"),
            ]
        )
        if invalid_pval and logger is not None:
            logger.qc(
                "p-value inside (0, 1]",
                "The p-value is not inside (0, 1]. Policy pvalue.out_of_range is "
                "'null', so it was emptied and no standard error was derived for "
                "those variants.",
                df.height,
                df.height,
                changed=invalid_pval,
                warn=True,
                step=STEP_LABEL,
            )
    else:
        # Clip invalid values to the configured valid range.
        df = df.with_columns(
            [
                cond_high.fill_null(False).alias(SE_FROM_CLIPPED_PVAL_COL),
                pl.when(cond_low)
                .then(pl.lit(clip_low))
                .when(cond_high)
                .then(pl.lit(_PVAL_CLIP_HIGH))
                .otherwise(pl.col(pval_col))
                .alias("_pval_safe"),
            ]
        )
        if invalid_pval and logger is not None:
            logger.qc(
                "p-value inside (0, 1]",
                "The p-value is not inside (0, 1]. Policy pvalue.out_of_range is "
                "'clip', so it was forced back into range before the standard error "
                "was derived. Standard errors built on a p-value above 1 are "
                "fabricated - they come out around 8,000,000 - so those variants are "
                "tagged and removed later under validation.se_from_clipped_pval.",
                df.height,
                df.height,
                changed=invalid_pval,
                warn=True,
                step=STEP_LABEL,
            )

    if zero_missing_se_action == "approximate":
        # pvalue.out_of_range=null may already have blanked the numeric value,
        # but REPORTED_ZERO_COLUMN retains its origin. The explicit
        # approximation must therefore restore only those rows to clip_low.
        df = df.with_columns(
            pl.when(zero_approximation)
            .then(pl.lit(clip_low))
            .otherwise(pl.col("_pval_safe"))
            .alias("_pval_safe")
        )
    upstream_clipped = _flag(df, CLIPPED_HIGH_COLUMN) | _flag(
        df, CLIPPED_LOW_COLUMN
    )
    df = df.with_columns(
        [
            (
                _flag(df, SE_FROM_CLIPPED_PVAL_COL)
                | (upstream_clipped & ~zero_approximation)
            ).alias(SE_FROM_CLIPPED_PVAL_COL),
            zero_approximation.alias(SE_FROM_ZERO_P_APPROXIMATION_COL),
        ]
    )
    n_clipped_pval = int(
        df.select(pl.col(SE_FROM_CLIPPED_PVAL_COL).sum()).item() or 0
    )
    n_zero_approximated = int(
        df.select(pl.col(SE_FROM_ZERO_P_APPROXIMATION_COL).sum()).item() or 0
    )
    if n_clipped_pval:
        _warn(
            logger,
            "{:,} variants require an SE derived from a p-value that was clipped. "
            "They are tagged in '{}' and removed by the effect-statistics validation step "
            "unless validation.se_from_clipped_pval is set to 'keep'.".format(
                n_clipped_pval, SE_FROM_CLIPPED_PVAL_COL
            ),
        )
    if n_zero_approximated:
        _warn(
            logger,
            "{:,} variant(s) reported p=0 with no supplied SE or usable Z. The "
            "explicit approximation policy used pvalue.clip_low={} and marked them "
            "in '{}'. Their true significance and exact SE remain unknown.".format(
                n_zero_approximated,
                clip_low,
                SE_FROM_ZERO_P_APPROXIMATION_COL,
            ),
        )
        if logger is not None:
            logger.qc(
                "zero p-value SE approximation",
                "The raw p-value was exactly zero and no supplied or usable "
                "Z-derived SE was available. The explicit approximation policy "
                "used pvalue.clip_low; the exact SE remains unknown.",
                df.height,
                df.height,
                changed=n_zero_approximated,
                warn=True,
                step=STEP_LABEL,
            )

    _info(
        logger,
        "Standard errors are calculated as |BETA / z| with z from the {}-sided normal "
        "tail (norm.isf), using effect size column '{}' and p-value column '{}'.".format(
            se_tail, beta_col, pval_col
        ),
    )

    # -------------------------------------------------------
    # Compute Z, then SE
    # -------------------------------------------------------
    # isf is the right tail, better conditioned than ppf on the left tail; the
    # magnitude is identical, and only the magnitude survives the abs().
    p_safe = df["_pval_safe"].to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        z = norm.isf(p_safe / se_tail)
        se = np.abs(df[beta_col].to_numpy() / z)

    n_undefined = int((~np.isfinite(se)).sum())

    df = df.with_columns(pl.Series("SE", se, dtype=pl.Float64))
    # A z of exactly 0 (p == tail/2) makes SE infinite or NaN.  Store those as
    # a missing value rather than NaN so the standard error gate removes them
    # with the reason 'se_null' instead of exporting a NaN.
    df = df.with_columns(
        pl.when(pl.col("SE").is_finite())
        .then(pl.col("SE"))
        .otherwise(None)
        .alias("SE")
    ).drop("_pval_safe")

    sample_column_dict["se_col"] = "SE"

    if n_undefined:
        _warn(
            logger,
            "{:,} variants have no usable standard error because the p-value gives a "
            "Z score of exactly zero or is missing. Their SE is empty, so the standard "
            "error check removes them rather than exporting a NaN.".format(n_undefined),
        )

    if logger is not None:
        logger.qc(
            "standard error derived from the p-value",
            "SE was derived as |BETA / z| for every variant.",
            df.height,
            df.height,
            changed=df.height - n_undefined,
            step=STEP_LABEL,
        )

    # -------------------------------------------------------
    # QC SUMMARY
    # -------------------------------------------------------
    stats = _se_statistics(df, "SE")

    qc_info.update(
        {
            "status": "SE calculated",
            "se_min": stats["min"],
            "se_max": stats["max"],
            "se_mean": stats["mean"],
            "se_median": stats["median"],
            "total_variants": stats["total"],
            "se_from_clipped_pval": n_clipped_pval,
            "se_from_zero_p_approximation": n_zero_approximated,
            "zero_p_missing_se": n_zero_missing_se,
            "zero_p_missing_se_action": zero_missing_se_action,
            "zero_p_missing_se_rejected": n_zero_rejected,
            "se_undefined": n_undefined,
            "variants_removed_total": n_initial - df.height,
            "variants_remaining": df.height,
        }
    )

    _info(
        logger,
        "Standard error summary: min={}, max={}, mean={}, median={} (n={:,}).".format(
            _fmt(stats["min"]),
            _fmt(stats["max"]),
            _fmt(stats["mean"]),
            _fmt(stats["median"]),
            stats["total"],
        ),
    )

    return df, qc_info, sample_column_dict
