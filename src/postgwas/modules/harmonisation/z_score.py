"""Step 09 (per chromosome) — produce Z from BETA and SE.

Public entry point: ``derive_z_score_from_effect_and_standard_error``.

    Z = BETA / SE

A study that already supplies a Z score keeps it. Otherwise this step computes
``Z = BETA / SE`` only where BETA and SE are finite and SE is above the
configured numerical division floor. The immediately following effect-
statistics step owns the six basic BETA, SE and Z checks and applies their
policies exactly once.

Logging and rejected variants
-----------------------------
Pass ``logger`` (a ``PipelineLogger``) and ``policies`` (a ``Policies``). This
calculation step never removes a variant; invalid-value rejection and its
provenance belong to step 10.

``POLICY_KEYS`` lists every policy this step reads.
"""

from typing import Any, Dict, Optional, Tuple

import polars as pl

from postgwas.core.values import optional_text

from .shared.runtime import log_info, resolve_policies


__all__ = ["derive_z_score_from_effect_and_standard_error", "POLICY_KEYS", "STEP_LABEL"]


#: Every policy this step reads.  Pass to ``logger.step(policy_keys=...)``.
POLICY_KEYS = (
    "validation.se_division_floor",
)

#: The free-text step label recorded in the reject file.
STEP_LABEL = "09 z_from_beta_se"

# Shared mechanics stay outside the scientific implementation.
_info = log_info


# public entry point
# ---------------------------------------------------------------------------
def derive_z_score_from_effect_and_standard_error(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    logger: Optional[Any] = None,
    policies: Optional[Any] = None,
) -> Tuple[pl.DataFrame, Dict, dict]:
    """Compute Z = BETA / SE, or retain a study-supplied Z score.

    This step performs no policy-driven filtering. A finite zero BETA produces
    Z=0; unusable inputs produce a null calculated Z and are handled once by
    :func:`effect_validation.validate_effect_statistics` in step 10.
    """
    pol = resolve_policies(policies)
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

    # BETA and SE are required only when there is no Z score to retain.
    if not has_beta_se and not has_existing_z:
        raise ValueError(
            "Chromosome {}: a Z score cannot be produced. The effect size column "
            "('{}') and the standard error column ('{}') are not both present, and the "
            "file supplies no Z score column either.".format(chromosome, beta_col, se_col)
        )

    if has_beta_se:
        df = df.with_columns(
            pl.col(beta_col).cast(pl.Float64, strict=False),
            pl.col(se_col).cast(pl.Float64, strict=False),
        )
    else:
        _info(
            logger,
            "The effect size and standard error columns are not both present. The "
            "study-supplied Z score in column '{}' is retained for the single "
            "validation gate in step 10.".format(existing_z_col),
        )

    division_floor = float(pol.get("validation.se_division_floor"))
    if not has_existing_z:
        beta = pl.col(beta_col)
        se = pl.col(se_col)
        df = df.with_columns(
            pl.when(
                beta.is_not_null()
                & beta.is_finite()
                & se.is_not_null()
                & se.is_finite()
                & (se > division_floor)
            )
            .then(beta / se)
            .otherwise(None)
            .alias("imp_z_col")
        )
        sample_column_dict["imp_z_col"] = "imp_z_col"
        _info(
            logger,
            "The study supplies no Z score, so Z was calculated as BETA / SE from "
            "columns '{}' and '{}'. Rows with unusable inputs were left with a null "
            "calculated Z for the single validation gate in step 10.".format(
                beta_col, se_col
            ),
        )

    qc_info.update({
        "z_column": sample_column_dict["imp_z_col"],
        "se_division_floor": division_floor,
        "total_variants": df.height,
    })

    return df, qc_info, sample_column_dict
