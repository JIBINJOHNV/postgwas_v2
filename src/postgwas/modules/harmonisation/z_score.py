"""Step 10 (per chromosome) — Z score from BETA and SE.

Public entry point: ``derive_z_score_from_effect_and_standard_error``.

    Z = BETA / SE

A study that already supplies a Z score keeps it; the BETA/SE quality checks
still run when both of those columns exist, and are skipped when they do not.
That is the point of the step: a file with a Z score but no effect size and no
standard error is exactly the case the pipeline exists to serve, and it must not
abort.

Logging and rejected variants
-----------------------------
Pass ``logger`` (a ``PipelineLogger``), ``policies`` (a ``Policies``) and
``rejects`` (a ``RejectCollector``). The pipeline supplies the collector so
every removed variant is written once to the unified reject file.

``POLICY_KEYS`` lists every policy this step reads.
"""

import polars as pl

from functools import partial
from typing import Any, Dict, Optional, Tuple

from postgwas.core.dataframes import count_matching_rows
from postgwas.core.values import format_number, optional_text

from .shared.runtime import log_info, log_warning, reject_rows, resolve_policies


__all__ = ["derive_z_score_from_effect_and_standard_error", "POLICY_KEYS", "STEP_LABEL"]


#: Every policy this step reads.  Pass to ``logger.step(policy_keys=...)``.
POLICY_KEYS = (
    "validation.beta_invalid",
    "validation.beta_zero",
    "validation.se_invalid",
    "validation.se_division_floor",
    "validation.z_invalid",
)

#: The free-text step label recorded in the reject file.
STEP_LABEL = "09 z_from_beta_se"

# Shared mechanics stay outside the scientific implementation.
_info = log_info
_warn = log_warning
_fmt = partial(format_number, pattern="%.6f", missing="n/a")
_count = partial(count_matching_rows, null_is_match=True)


# public entry point
# ---------------------------------------------------------------------------
def derive_z_score_from_effect_and_standard_error(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    logger: Optional[Any] = None,
    policies: Optional[Any] = None,
    rejects: Optional[Any] = None,
    validation_state: Optional[Dict[str, Any]] = None,
) -> Tuple[pl.DataFrame, Dict, dict]:
    """Compute Z = BETA / SE, or keep the Z score the study already supplied.

    Zero-effect handling is controlled exclusively by
    ``validation.beta_zero``. The default is ``keep`` because a beta of exactly
    zero with a valid standard error is a legitimate null estimate giving Z=0.

    ``validation_state`` is an internal, row-safe hand-off to the immediately
    following effect-validation step.  It records which basic checks were
    already handled and retains the exact resulting DataFrame object, avoiding
    a second scan without trusting masks from an earlier row layout.
    """
    pol = resolve_policies(policies)

    if validation_state is not None:
        validation_state.clear()

    qc_info = {"initial_variants": df.height}  # type: Dict[str, Any]

    if logger is not None:
        logger.settings(list(POLICY_KEYS))

    beta_col = optional_text(sample_column_dict.get("beta_col"))
    se_col = optional_text(sample_column_dict.get("se_col"))

    # ==========================================================
    # EXISTING Z COLUMN
    # ==========================================================
    existing_z_col = optional_text(sample_column_dict.get("imp_z_col"))
    has_existing_z = existing_z_col is not None and existing_z_col in df.columns

    if has_existing_z:
        df = df.with_columns(
            pl.col(existing_z_col).cast(pl.Float64, strict=False).alias(existing_z_col)
        )
        sample_column_dict["imp_z_col"] = existing_z_col
        _info(
            logger,
            "The study supplies a Z score in column '{}', so it is used as given and "
            "not recalculated.".format(existing_z_col),
        )

    has_beta_se = (
        beta_col is not None
        and se_col is not None
        and beta_col in df.columns
        and se_col in df.columns
    )
    qc_info.update({
        "z_source": (
            "study_column:%s" % existing_z_col
            if has_existing_z else "calculated_from_beta_se"
        ),
        "beta_column": beta_col,
        "se_column": se_col,
    })

    # BETA and SE are only *required* when there is no Z score to fall back on.
    if not has_beta_se and not has_existing_z:
        raise ValueError(
            "Chromosome {}: a Z score cannot be produced. The effect size column "
            "('{}') and the standard error column ('{}') are not both present, and the "
            "file supplies no Z score column either.".format(chromosome, beta_col, se_col)
        )

    division_floor = float(pol.get("validation.se_division_floor"))
    beta_zero_action = pol.get("validation.beta_zero")
    beta_invalid_action = pol.get("validation.beta_invalid")
    se_invalid_action = pol.get("validation.se_invalid")
    z_invalid_action = pol.get("validation.z_invalid")

    n_initial = qc_info["initial_variants"]
    n_beta_null = 0
    n_beta_zero = 0
    n_se_null = 0
    n_se_non_finite = 0
    n_se_invalid = 0

    if not has_beta_se:
        _info(
            logger,
            "The effect size and standard error columns are not both present, so the "
            "BETA and SE quality checks were skipped. The Z score column '{}' supplied "
            "by the study is used as it stands.".format(existing_z_col),
        )
    else:
        df = df.with_columns(
            [
                pl.col(beta_col).cast(pl.Float64, strict=False),
                pl.col(se_col).cast(pl.Float64, strict=False),
            ]
        )
        if not has_existing_z:
            _info(
                logger,
                "The study supplies no Z score, so it is calculated as BETA / SE from "
                "column '{}' divided by column '{}'.".format(beta_col, se_col),
            )

        # --- BETA missing or infinite ---------------------------------------
        cond_beta_null = pl.col(beta_col).is_null() | ~pl.col(beta_col).is_finite()
        if beta_invalid_action == "fail":
            n_beta_null = _count(df, cond_beta_null)
            if n_beta_null:
                raise ValueError(
                    "Chromosome {}: {:,} variants have a missing or non-finite effect "
                    "size in column '{}' and policy validation.beta_invalid is "
                    "'fail'.".format(chromosome, n_beta_null, beta_col)
                )
        elif beta_invalid_action == "keep":
            n_beta_null = _count(df, cond_beta_null)
            if n_beta_null:
                _warn(
                    logger,
                    "{:,} variants have a missing or non-finite effect size, but policy "
                    "validation.beta_invalid is 'keep', so they were left in "
                    "place.".format(n_beta_null),
                )
        else:
            df, n_beta_null = reject_rows(
                df,
                cond_beta_null,
                step_label=STEP_LABEL,
                reason="beta_invalid",
                check_name="effect size usable",
                description=(
                    "The effect size is missing or not a finite number, so no Z "
                    "score can be calculated from it."
                ),
                collector=rejects,
                context=logger,
                detail="column '{}'".format(beta_col),
            )

        # --- SE missing, then infinite, then too small for stable division ---
        # Order matters: `inf > 0` is True, so finiteness must be tested first.
        cond_se_null = pl.col(se_col).is_null()
        cond_se_non_finite = pl.col(se_col).is_not_null() & ~pl.col(se_col).is_finite()
        cond_se_non_positive = (
            pl.col(se_col).is_not_null()
            & pl.col(se_col).is_finite()
            & (pl.col(se_col) <= division_floor)
        )

        if se_invalid_action == "fail":
            n_se_null = _count(df, cond_se_null)
            n_se_non_finite = _count(df, cond_se_non_finite)
            n_se_invalid = _count(df, cond_se_non_positive)
            if n_se_null + n_se_non_finite + n_se_invalid:
                raise ValueError(
                    "Chromosome {}: {:,} variants have a standard error that is "
                    "missing, infinite or at or below the configured division floor "
                    "({}) in column '{}' and policy "
                    "validation.se_invalid is 'fail'.".format(
                        chromosome,
                        n_se_null + n_se_non_finite + n_se_invalid,
                        division_floor,
                        se_col,
                    )
                )
        elif se_invalid_action == "keep":
            n_se_null = _count(df, cond_se_null)
            n_se_non_finite = _count(df, cond_se_non_finite)
            n_se_invalid = _count(df, cond_se_non_positive)
            if n_se_null + n_se_non_finite + n_se_invalid:
                _warn(
                    logger,
                    "{:,} variants have a standard error that is missing, infinite or "
                    "at or below the configured division floor, but policy "
                    "validation.se_invalid is 'keep', so they "
                    "were left in place and their Z score will be empty.".format(
                        n_se_null + n_se_non_finite + n_se_invalid
                    ),
                )
        else:
            df, n_se_null = reject_rows(
                df,
                cond_se_null,
                step_label=STEP_LABEL,
                reason="se_null",
                check_name="standard error present",
                description=(
                    "The standard error is missing, so no Z score can be calculated."
                ),
                collector=rejects,
                context=logger,
                detail="column '{}'".format(se_col),
            )
            df, n_se_non_finite = reject_rows(
                df,
                cond_se_non_finite,
                step_label=STEP_LABEL,
                reason="se_non_finite",
                check_name="standard error finite",
                description=(
                    "The standard error is infinite or not a number, which no later "
                    "positivity test would catch because infinity compares as "
                    "greater than zero."
                ),
                collector=rejects,
                context=logger,
                detail="column '{}'".format(se_col),
            )

        # --- BETA exactly zero ----------------------------------------------
        cond_beta_zero = pl.col(beta_col) == 0
        if beta_zero_action == "fail":
            n_beta_zero = _count(df, cond_beta_zero)
            if n_beta_zero:
                raise ValueError(
                    "Chromosome {}: {:,} variants have an effect size of exactly 0 in "
                    "column '{}' and policy validation.beta_zero is 'fail'.".format(
                        chromosome, n_beta_zero, beta_col
                    )
                )
        elif beta_zero_action == "reject":
            df, n_beta_zero = reject_rows(
                df,
                cond_beta_zero,
                step_label=STEP_LABEL,
                reason="beta_zero",
                check_name="effect size exactly zero",
                description=(
                    "The effect size is exactly 0. Policy validation.beta_zero is "
                    "'reject', so these variants were removed and recorded."
                ),
                collector=rejects,
                context=logger,
                detail="column '{}'".format(beta_col),
            )
        else:
            n_beta_zero = _count(df, cond_beta_zero)
            if logger is not None and n_beta_zero:
                logger.qc(
                    "effect size exactly zero",
                    "An effect size of exactly 0 with a valid standard error is a "
                    "legitimate null estimate giving a Z score of 0. Policy "
                    "validation.beta_zero is 'keep', so these variants were kept.",
                    df.height,
                    df.height,
                    reason="beta_zero",
                    matched=n_beta_zero,
                    outcome="kept",
                    step=STEP_LABEL,
                )

        # --- SE at or below the configured division floor -------------------
        if se_invalid_action == "reject":
            df, n_se_invalid = reject_rows(
                df,
                cond_se_non_positive,
                step_label=STEP_LABEL,
                reason="se_non_positive",
                check_name="standard error large enough for stable division",
                description=(
                    "The standard error is at or below validation.se_division_floor "
                    "(%s), so BETA / SE would be numerically unstable."
                    % division_floor
                ),
                collector=rejects,
                context=logger,
                detail="column '{}'".format(se_col),
            )

    # ==========================================================
    # COMPUTE Z
    # ==========================================================
    if not has_existing_z:
        df = df.with_columns(
            pl.when(pl.col(se_col) > division_floor)
            .then(pl.col(beta_col) / pl.col(se_col))
            .otherwise(None)
            .alias("imp_z_col")
        )
        sample_column_dict["imp_z_col"] = "imp_z_col"

    z_col = sample_column_dict["imp_z_col"]

    cond_z_invalid = pl.col(z_col).is_null() | ~pl.col(z_col).is_finite()
    n_z_invalid = 0
    if z_invalid_action == "fail":
        n_z_invalid = _count(df, cond_z_invalid)
        if n_z_invalid:
            raise ValueError(
                "Chromosome {}: {:,} variants have a missing or non-finite Z score in "
                "column '{}' and policy validation.z_invalid is 'fail'.".format(
                    chromosome, n_z_invalid, z_col
                )
            )
    elif z_invalid_action == "keep":
        n_z_invalid = _count(df, cond_z_invalid)
        if n_z_invalid:
            _warn(
                logger,
                "{:,} variants have a missing or non-finite Z score, but policy "
                "validation.z_invalid is 'keep', so they were left in place.".format(
                    n_z_invalid
                ),
            )
    else:
        df, n_z_invalid = reject_rows(
            df,
            cond_z_invalid,
            step_label=STEP_LABEL,
            reason="z_invalid",
            check_name="Z score usable",
            description="The Z score is missing or not a finite number.",
            collector=rejects,
            context=logger,
            detail="column '{}'".format(z_col),
        )

    # Every removal above came out of `df`, so the difference is the truth and
    # cannot drift from the per-reason counts.  Under a 'keep' or 'fail' policy
    # the per-reason numbers below are variants *detected*, not removed.
    n_final = df.height
    n_removed = n_initial - n_final
    pct_removed = (n_removed / n_initial * 100) if n_initial > 0 else 0.0

    _info(
        logger,
        "BETA/SE quality checks:\n"
        "Variants in: {:,}\n"
        "Removed: {:,} ({:.2f}%)\n"
        "  unusable effect size: {:,}\n"
        "  missing standard error: {:,}\n"
        "  non-finite standard error: {:,}\n"
        "  effect size exactly zero: {:,} ({})\n"
        "  standard error at/below division floor: {:,}\n"
        "  unusable Z score: {:,}\n"
        "Variants out: {:,}".format(
            n_initial,
            n_removed,
            pct_removed,
            n_beta_null,
            n_se_null,
            n_se_non_finite,
            n_beta_zero,
            "removed" if beta_zero_action == "reject" else "kept",
            n_se_invalid,
            n_z_invalid,
            n_final,
        ),
    )

    qc_info.update(
        {
            "variants_removed_due_to_null_beta": n_beta_null,
            "variants_removed_due_to_null_se": n_se_null + n_se_non_finite,
            "variants_removed_due_to_missing_se": n_se_null,
            "variants_removed_due_to_non_finite_se": n_se_non_finite,
            "variants_removed_due_to_zero_beta": (
                n_beta_zero if beta_zero_action == "reject" else 0
            ),
            "variants_with_zero_beta": n_beta_zero,
            "variants_removed_due_to_invalid_se": n_se_invalid,
            "variants_removed_due_to_invalid_z": n_z_invalid,
            "variants_removed_invalid_beta_se": n_removed,
            "variants_after_filter_invalid_beta_se": n_final,
            "beta_zero_removed_flag": beta_zero_action == "reject",
            "beta_zero_action": beta_zero_action,
            "se_division_floor": division_floor,
        }
    )

    if n_final == 0:
        raise RuntimeError(
            "Chromosome {}: no valid variants are left after the BETA, SE and Z score "
            "checks. Started with {:,} variants and removed all of them.".format(
                chromosome, n_initial
            )
        )

    # ==========================================================
    # Z SUMMARY
    # ==========================================================
    z_summary = df.select(
        [
            pl.col(z_col).min().alias("z_min"),
            pl.col(z_col).max().alias("z_max"),
            pl.col(z_col).mean().alias("z_mean"),
            pl.col(z_col).std().alias("z_std"),
            pl.len().alias("total"),
        ]
    ).to_dicts()[0]

    _info(
        logger,
        "Z score summary: min={}, max={}, mean={}, std={} (n={:,}).".format(
            _fmt(z_summary["z_min"]),
            _fmt(z_summary["z_max"]),
            _fmt(z_summary["z_mean"]),
            _fmt(z_summary["z_std"]),
            z_summary["total"],
        ),
    )

    qc_info.update(
        {
            "z_min": z_summary["z_min"],
            "z_max": z_summary["z_max"],
            "z_mean": z_summary["z_mean"],
            "z_std": z_summary["z_std"],
            "total_variants": z_summary["total"],
        }
    )

    if validation_state is not None:
        completed = ["z_invalid"]
        if has_beta_se:
            completed.extend(("se_invalid", "beta_invalid", "beta_zero"))
        validation_state.update(
            {
                "frame": df,
                "row_count": df.height,
                "columns": {
                    "se": se_col if has_beta_se else None,
                    "beta": beta_col if has_beta_se else None,
                    "zscore": z_col,
                },
                "policies": {
                    "validation.se_invalid": se_invalid_action,
                    "validation.beta_invalid": beta_invalid_action,
                    "validation.beta_zero": beta_zero_action,
                    "validation.z_invalid": z_invalid_action,
                },
                "se_division_floor": division_floor,
                "completed": tuple(completed),
            }
        )

    return df, qc_info, sample_column_dict
