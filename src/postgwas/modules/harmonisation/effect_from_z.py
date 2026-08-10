"""Step 06 (per chromosome) — derive BETA and SE from an imputed Z score.

Public entry point: ``derive_effect_and_standard_error_from_z``.

The standardized-effect approximation from Zhu et al. (2016), PMID 27019110,
is::

    BETA = z / sqrt(2 p (1 - p) (Neff + z^2))
    SE   = 1 / sqrt(2 p (1 - p) (Neff + z^2))

where ``p`` is the effect allele frequency.  The denominator is *exactly* zero
when ``p`` is 0 or 1, and polars then yields ``inf`` — not null, not an error —
so a degenerate frequency silently becomes an infinite standard error that no
later positivity check can catch.  Those variants are therefore removed here,
at source, under the reason ``eaf_degenerate``.

This does not reconstruct a study's original regression coefficient.  It
estimates an additive, standardized effect under unit phenotype variance and
genotype variance ``2p(1-p)``.  It assumes that ``p`` represents the analysed
samples and that the supplied sample size describes the same association test.
For a binary trait, ``Neff`` is the balanced-design effective sample size from
the preceding sample-size step; the result is an approximation, not a direct
recovery of the logistic-regression log odds ratio or a liability-scale effect.

Reference: Zhu Z et al. Nature Genetics 2016; PMID 27019110;
https://pubmed.ncbi.nlm.nih.gov/27019110/

Logging and rejected variants
-----------------------------
Pass ``logger`` (a ``PipelineLogger``), ``policies`` (a ``Policies``) and
``rejects`` (a ``RejectCollector``). The pipeline supplies the collector so
every removed variant is written once to the unified reject file.

Give the logger the *same* ``Policies`` object you pass here — the settings
block is read from ``logger.policies``.  ``POLICY_KEYS`` lists every policy this
step reads, so a caller opening the step itself can pass
``policy_keys=POLICY_KEYS`` to ``logger.step(...)``.
"""

import polars as pl

from functools import partial
from typing import Any, Dict, Optional, Tuple

from postgwas.core.dataframes import count_matching_rows
from postgwas.core.values import format_number, optional_text

from .shared.runtime import (
    log_info,
    log_warning,
    reject_rows,
    resolve_policies,
)


__all__ = [
    "derive_effect_and_standard_error_from_z",
    "SE_UNAVAILABLE_FROM_Z_COLUMN",
    "POLICY_KEYS",
    "STEP_LABEL",
]


#: Internal row-level provenance used by the following standard-error step.
#: True means the study supplied no SE and BETA / Z could not recover one for
#: this row because Z or BETA was unusable. It is dropped before export with
#: the other ``__`` working columns.
SE_UNAVAILABLE_FROM_Z_COLUMN = "__se_unavailable_from_z"


#: Every policy this step reads.  Pass to ``logger.step(policy_keys=...)``.
POLICY_KEYS = (
    "eaf.degenerate",
    "validation.z_invalid",
    "sample_size.trait_type",
    "sample_size.min_value",
)

#: The free-text step label recorded in the reject file.
STEP_LABEL = "06 beta_se_from_z"


# Shared mechanics stay outside the scientific implementation.
_info = log_info
_warn = log_warning
_fmt = partial(format_number, pattern="%10.6f", missing="       n/a")
_count = partial(count_matching_rows, null_is_match=True)


# public entry point
# ---------------------------------------------------------------------------
def derive_effect_and_standard_error_from_z(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    logger: Optional[Any] = None,
    policies: Optional[Any] = None,
    rejects: Optional[Any] = None,
) -> Tuple[pl.DataFrame, Dict, dict]:
    """Fill in whichever of BETA and SE the study did not supply.

    Four combinations of BETA/SE presence are handled explicitly:

    ==============  =========  ==================================================
    BETA present    SE present  action
    ==============  =========  ==================================================
    yes             yes        nothing to do
    yes             no         SE = BETA / Z when a Z score exists, otherwise
                               left to ``derive_standard_error_from_effect_and_p_value``
    no              yes        BETA = Z * SE  (needs a Z score)
    no              no         both derived from Z, EAF and Neff
    ==============  =========  ==================================================
    """
    sample_column_dict = sample_column_dict.copy()  # SAFE
    pol = resolve_policies(policies)

    qc_info = {"initial_variants": df.height}  # type: Dict[str, Any]

    if logger is not None:
        logger.settings(list(POLICY_KEYS))

    # -------------------------------------------------------
    # Column mappings
    # -------------------------------------------------------
    imp_z_col = optional_text(sample_column_dict.get("imp_z_col"))
    beta_col = optional_text(sample_column_dict.get("beta_or_col"))
    se_col = optional_text(sample_column_dict.get("se_col"))
    eaf_col = optional_text(sample_column_dict.get("eaf_col"))

    has_z = imp_z_col is not None and imp_z_col in df.columns
    has_beta = beta_col is not None and beta_col in df.columns
    has_se = se_col is not None and se_col in df.columns

    for name, column in (
        ("beta_or_col", beta_col),
        ("se_col", se_col),
        ("imp_z_col", imp_z_col),
    ):
        if column is not None and column not in df.columns:
            _warn(
                logger,
                "The config sets {} = '{}' but that column is not in the data for "
                "chromosome {}, so it is treated as absent.".format(name, column, chromosome),
            )

    qc_info.update({
        "has_beta": has_beta,
        "has_se": has_se,
        "has_zscore": has_z,
        "effect_column": beta_col if has_beta else None,
        "se_column": se_col if has_se else None,
        "z_column": imp_z_col if has_z else None,
    })

    # -------------------------------------------------------
    # Combination 1 of 4 — BETA and SE both supplied.
    # -------------------------------------------------------
    if has_beta and has_se:
        if has_z:
            _info(
                logger,
                "The study supplies both an effect size ('{}') and a standard error "
                "('{}'), so nothing is derived from the Z score column '{}'.".format(
                    beta_col, se_col, imp_z_col
                ),
            )
        else:
            _info(
                logger,
                "The study supplies both an effect size ('{}') and a standard error "
                "('{}'), so nothing is derived here. There is no Z score column; one "
                "is calculated from BETA and SE in a later step.".format(beta_col, se_col),
            )
        qc_info.update(
            {
                "status": "Beta and SE columns are present",
                "beta_computed": False,
                "se_computed": False,
                "variants_removed_total": 0,
                "variants_remaining": df.height,
            }
        )
        return df, qc_info, sample_column_dict

    # -------------------------------------------------------
    # Nothing to work from: no effect size and no Z score.
    # -------------------------------------------------------
    if not has_beta and not has_z:
        qc_info["status"] = "skipped_no_zscore_or_effectsize_column"
        raise KeyError(
            "Chromosome {}: the summary statistics have neither an effect size "
            "column nor a Z score column, so BETA and SE cannot be produced. "
            "Provide beta_or_col (with se_col) or imp_z_col.".format(chromosome)
        )

    # -------------------------------------------------------
    # Combination 2 of 4 — BETA supplied, SE missing.
    # -------------------------------------------------------
    if has_beta and not has_se:
        if not has_z:
            _info(
                logger,
                "The study supplies an effect size ('{}') but no standard error, and "
                "there is no Z score column to recover one from. SE is calculated in a "
                "later step from BETA and the p-value.".format(beta_col),
            )
            qc_info.update(
                {
                    "status": "Beta column is present, but SE column is missing",
                    "beta_computed": False,
                    "se_computed": False,
                    "variants_removed_total": 0,
                    "variants_remaining": df.height,
                }
            )
            return df, qc_info, sample_column_dict

        df = df.with_columns(
            [
                pl.col(beta_col).cast(pl.Float64, strict=False),
                pl.col(imp_z_col).cast(pl.Float64, strict=False),
            ]
        )
        # SE = BETA / Z is exact: Z is defined as BETA / SE.  A zero or
        # non-finite Z leaves SE null rather than infinite, so the standard
        # error gate removes it with a reason instead of exporting inf.
        usable_for_se = (
                (pl.col(imp_z_col) != 0)
                & pl.col(imp_z_col).is_finite()
                & pl.col(beta_col).is_finite()
        )
        df = df.with_columns(
            [
            pl.when(usable_for_se)
            .then(pl.col(beta_col) / pl.col(imp_z_col))
            .otherwise(None)
            .alias("SE"),
            (~usable_for_se).fill_null(True).alias(SE_UNAVAILABLE_FROM_Z_COLUMN),
            ]
        )
        n_undefined = int(df.select(pl.col("SE").is_null().sum()).item() or 0)

        sample_column_dict["se_col"] = "SE"
        _info(
            logger,
            "The study supplies an effect size ('{}') but no standard error, so SE was "
            "recovered exactly as BETA / Z from the Z score column '{}'.".format(
                beta_col, imp_z_col
            ),
        )
        if logger is not None:
            logger.qc(
                "standard error recovered from Z",
                "SE was derived as BETA / Z for every variant with a usable Z score.",
                df.height,
                df.height,
                changed=df.height - n_undefined,
                step=STEP_LABEL,
                warn=n_undefined > 0,
            )
        if n_undefined > 0:
            _warn(
                logger,
                "{:,} variants have a zero or non-finite Z score, so no standard error "
                "could be derived for them. The following standard-error step applies "
                "pvalue.zero_missing_se to any row whose raw p-value was zero; other "
                "empty SE values reach the configured standard-error validation gate."
                .format(n_undefined),
            )
        qc_info.update(
            {
                "status": "SE derived from BETA and Z",
                "beta_computed": False,
                "se_computed": True,
                "se_undefined_from_z": n_undefined,
                "variants_removed_total": 0,
                "variants_remaining": df.height,
            }
        )
        return df, qc_info, sample_column_dict

    # -------------------------------------------------------
    # Combination 3 of 4 — SE supplied, BETA missing.
    # (has_z is guaranteed here: the no-BETA-no-Z case raised above.)
    # -------------------------------------------------------
    if has_se and not has_beta:
        df = df.with_columns(
            [
                pl.col(se_col).cast(pl.Float64, strict=False),
                pl.col(imp_z_col).cast(pl.Float64, strict=False),
            ]
        )
        # BETA = Z * SE is exact and needs no allele frequency.
        df = df.with_columns(
            pl.when(pl.col(imp_z_col).is_finite() & pl.col(se_col).is_finite())
            .then(pl.col(imp_z_col) * pl.col(se_col))
            .otherwise(None)
            .alias("BETA")
        )
        n_undefined = int(df.select(pl.col("BETA").is_null().sum()).item() or 0)

        sample_column_dict["beta_or_col"] = "BETA"
        sample_column_dict["beta_col"] = "BETA"
        _info(
            logger,
            "The study supplies a standard error ('{}') but no effect size, so BETA was "
            "recovered exactly as Z * SE from the Z score column '{}'. No allele "
            "frequency is needed for this recovery.".format(se_col, imp_z_col),
        )
        if logger is not None:
            logger.qc(
                "effect size recovered from Z",
                "BETA was derived as Z * SE for every variant with a usable Z score "
                "and standard error.",
                df.height,
                df.height,
                changed=df.height - n_undefined,
                step=STEP_LABEL,
                warn=n_undefined > 0,
            )
        if n_undefined > 0:
            _warn(
                logger,
                "{:,} variants have a non-finite Z score or standard error, so no effect "
                "size could be derived for them; their BETA is empty and the effect size "
                "check removes them.".format(n_undefined),
            )
        qc_info.update(
            {
                "status": "BETA derived from Z and SE",
                "beta_computed": True,
                "se_computed": False,
                "beta_undefined_from_z": n_undefined,
                "variants_removed_total": 0,
                "variants_remaining": df.height,
            }
        )
        return df, qc_info, sample_column_dict

    # -------------------------------------------------------
    # Combination 4 of 4 — neither BETA nor SE; derive both from Z.
    # -------------------------------------------------------
    _info(
        logger,
        "Neither an effect size nor a standard error was supplied, so both are derived "
        "from the Z score column '{}' using the effect allele frequency '{}' and the "
        "effective sample size.".format(imp_z_col, eaf_col),
    )

    if eaf_col is None or eaf_col not in df.columns:
        raise KeyError(
            "Chromosome {}: BETA and SE must be derived from the Z score column '{}', "
            "which needs an effect allele frequency, but the configured column '{}' is "
            "not present.".format(chromosome, imp_z_col, eaf_col)
        )

    if "Neff" not in df.columns:
        raise KeyError(
            "Chromosome {}: BETA and SE must be derived from the Z score column '{}', "
            "which needs the effective sample size, but the 'Neff' column is not "
            "present. It is produced by the sample size step.".format(chromosome, imp_z_col)
        )

    # -------------------------------------------------------
    # Cast numeric
    # -------------------------------------------------------
    df = df.with_columns(
        [
            pl.col(imp_z_col).cast(pl.Float64, strict=False),
            pl.col(eaf_col).cast(pl.Float64, strict=False),
            pl.col("Neff").cast(pl.Float64, strict=False),
        ]
    )

    z_invalid_action = pol.get("validation.z_invalid")
    degenerate_action = pol.get("eaf.degenerate")
    trait_type = str(pol.get("sample_size.trait_type"))
    min_sample_size = pol.get("sample_size.min_value")

    n_initial = qc_info["initial_variants"]

    # --- 1. Z score ---------------------------------------------------------
    cond_z_invalid = pl.col(imp_z_col).is_null() | ~pl.col(imp_z_col).is_finite()
    n_z_null = 0
    if z_invalid_action == "fail":
        n_z_null = _count(df, cond_z_invalid)
        if n_z_null:
            raise ValueError(
                "Chromosome {}: {:,} variants have a missing or non-finite Z score in "
                "column '{}' and policy validation.z_invalid is 'fail'.".format(
                    chromosome, n_z_null, imp_z_col
                )
            )
    elif z_invalid_action == "keep":
        n_z_null = _count(df, cond_z_invalid)
        if n_z_null:
            _warn(
                logger,
                "{:,} variants have a missing or non-finite Z score, but policy "
                "validation.z_invalid is 'keep', so they were left in place and their "
                "BETA and SE will be empty.".format(n_z_null),
            )
    else:
        df, n_z_null = reject_rows(
            df,
            cond_z_invalid,
            step_label=STEP_LABEL, reason="z_invalid",
            check_name="Z score usable",
            description="The Z score is missing or not a finite number, so no effect size can be derived from it.",
            collector=rejects, context=logger,
            detail="column '{}'".format(imp_z_col),
        )

    # --- 2. Effect allele frequency: missing --------------------------------
    cond_eaf_null = pl.col(eaf_col).is_null() | ~pl.col(eaf_col).is_finite()
    df, n_eaf_null = reject_rows(
        df,
        cond_eaf_null,
        step_label=STEP_LABEL, reason="eaf_null",
        check_name="effect allele frequency present",
        description="The effect allele frequency is missing or not a finite number, so the denominator of the effect size formula cannot be evaluated.",
        collector=rejects, context=logger,
        detail="column '{}'".format(eaf_col),
    )

    # --- 3. Effect allele frequency: degenerate -----------------------------
    # At exactly 0 or 1 the denominator 2*p*(1-p)*(Neff + z^2) is exactly zero
    # and polars yields inf, so the predicate must be <= 0 / >= 1, not < 0 / > 1.
    cond_eaf_lt_zero = pl.col(eaf_col) <= 0
    cond_eaf_gt_one = pl.col(eaf_col) >= 1
    cond_eaf_degenerate = cond_eaf_lt_zero | cond_eaf_gt_one

    edge = df.select(
        [
            cond_eaf_lt_zero.fill_null(False).sum().alias("le_zero"),
            cond_eaf_gt_one.fill_null(False).sum().alias("ge_one"),
        ]
    ).to_dicts()[0]
    n_eaf_lt_zero = int(edge["le_zero"] or 0)
    n_eaf_gt_one = int(edge["ge_one"] or 0)
    n_eaf_invalid = 0

    if degenerate_action == "fail":
        if n_eaf_lt_zero + n_eaf_gt_one:
            raise ValueError(
                "Chromosome {}: {:,} variants have an effect allele frequency of exactly "
                "0 or 1 (or outside that range) in column '{}' and policy eaf.degenerate "
                "is 'fail'.".format(chromosome, n_eaf_lt_zero + n_eaf_gt_one, eaf_col)
            )
    elif degenerate_action == "keep":
        if n_eaf_lt_zero + n_eaf_gt_one:
            _warn(
                logger,
                "{:,} variants have an effect allele frequency of exactly 0 or 1, but "
                "policy eaf.degenerate is 'keep', so they were left in place. Their "
                "standard error will be infinite.".format(n_eaf_lt_zero + n_eaf_gt_one),
            )
    elif degenerate_action == "null":
        if n_eaf_lt_zero + n_eaf_gt_one:
            before = df.height
            df = df.with_columns(
                pl.when(cond_eaf_degenerate.fill_null(False))
                .then(None)
                .otherwise(pl.col(eaf_col))
                .alias(eaf_col)
            )
            if logger is not None:
                logger.qc(
                    "degenerate effect allele frequency",
                    "Frequencies of exactly 0 or 1 make the standard error formula "
                    "divide by zero. Policy eaf.degenerate is 'null', so the frequency "
                    "was emptied and BETA and SE are left empty for those variants.",
                    before,
                    df.height,
                    changed=n_eaf_lt_zero + n_eaf_gt_one,
                    warn=True,
                    step=STEP_LABEL,
                )
    else:
        df, n_eaf_invalid = reject_rows(
            df,
            cond_eaf_degenerate,
            step_label=STEP_LABEL, reason="eaf_degenerate",
            check_name="degenerate effect allele frequency",
            description="The effect allele frequency is exactly 0 or 1, which makes the standard error formula divide by zero and yield an infinite value.",
            collector=rejects, context=logger,
            detail="{:,} at or below 0, {:,} at or above 1".format(
                n_eaf_lt_zero, n_eaf_gt_one
            ),
        )

    # --- 4. Effective sample size -------------------------------------------
    cond_neff_null = pl.col("Neff").is_null() | ~pl.col("Neff").is_finite()
    df, n_neff_null = reject_rows(
        df,
        cond_neff_null,
        step_label=STEP_LABEL, reason="neff_invalid",
        check_name="effective sample size present",
        description="The effective sample size is missing or not a finite number, so the denominator of the effect size formula cannot be evaluated.",
        collector=rejects, context=logger,
        detail="Neff missing or non-finite",
    )

    # sample_size.min_value sets the minimum usable effective sample size.
    cond_neff_invalid = (pl.col("Neff") <= 0) | (pl.col("Neff") < min_sample_size)
    df, n_neff_invalid = reject_rows(
        df,
        cond_neff_invalid,
        step_label=STEP_LABEL, reason="neff_invalid",
        check_name="effective sample size usable",
        description="The effective sample size is not a positive number, so the effect size formula would divide by zero or return a complex value.",
        collector=rejects, context=logger,
        detail="Neff at or below {}".format(min_sample_size)
        if min_sample_size > 0
        else "Neff at or below 0",
    )

    # Every removal above came out of `df`, so the difference is the truth and
    # cannot drift from the per-reason counts.  Under a 'keep' or 'fail' policy
    # the per-reason numbers are variants *detected*, not removed.
    n_final = df.height
    n_removed = n_initial - n_final
    pct_removed = (n_removed / n_initial * 100) if n_initial > 0 else 0.0

    qc_info.update(
        {
            "removed_z_null": n_z_null,
            "removed_eaf_null": n_eaf_null,
            "removed_eaf_invalid": n_eaf_invalid,
            "removed_eaf_lt_zero": n_eaf_lt_zero,
            "removed_eaf_gt_one": n_eaf_gt_one,
            "removed_neff_null": n_neff_null,
            "removed_neff_invalid": n_neff_invalid,
            "variants_removed_total": n_removed,
            "variants_remaining": n_final,
        }
    )

    if n_final == 0:
        raise ValueError(
            "Chromosome {}: zero usable variants remain, so BETA and SE cannot be "
            "derived from the Z score. Started with {:,} variants and removed all of "
            "them ({:,} unusable Z score, {:,} missing frequency, {:,} degenerate "
            "frequency, {:,} unusable effective sample size). Check the Z score, "
            "frequency and sample size columns for this chromosome.".format(
                chromosome,
                n_initial,
                n_z_null,
                n_eaf_null,
                n_eaf_invalid,
                n_neff_null + n_neff_invalid,
            )
        )

    # -------------------------------------------------------
    # COMPUTE BETA AND SE
    # -------------------------------------------------------
    df = df.with_columns(
        (
            2
            * pl.col(eaf_col)
            * (1 - pl.col(eaf_col))
            * (pl.col("Neff") + pl.col(imp_z_col) ** 2)
        ).alias("_den")
    )

    df = df.with_columns(
        [
            (pl.col(imp_z_col) / pl.col("_den").sqrt()).alias("BETA"),
            (1 / pl.col("_den").sqrt()).alias("SE"),
        ]
    ).drop("_den")

    sample_column_dict["beta_or_col"] = "BETA"
    sample_column_dict["beta_col"] = "BETA"
    sample_column_dict["se_col"] = "SE"
    beta_col = "BETA"
    se_col = "SE"
    qc_info["beta_computed"] = True
    qc_info["se_computed"] = True
    qc_info["effect_estimate_scale"] = "standardized_effect_estimate"
    qc_info["effect_estimate_citation"] = "Zhu et al. 2016; PMID 27019110"
    qc_info["effect_estimate_trait_type"] = trait_type
    qc_info["effect_estimate_uses_neff"] = True
    qc_info["effect_estimate_eaf_assumption"] = "same_analysed_samples_as_z"
    qc_info["effect_estimate_binary_interpretation"] = (
        "approximate_standardized_not_log_odds_or_liability_scale"
        if trait_type == "binary"
        else None
    )

    _info(
        logger,
        "A standardized BETA estimate and its SE were calculated as "
        "z / sqrt(2p(1-p)(Neff + z^2)) and "
        "1 / sqrt(2p(1-p)(Neff + z^2)) (Zhu et al. 2016; PMID 27019110).",
    )
    _warn(
        logger,
        "The standardized BETA and SE approximation assumes that EAF was measured in "
        "the same analysed samples as the Z score. If EAF comes from an external "
        "reference panel, ancestry or sample-frequency differences can make the "
        "derived estimates inaccurate.",
    )
    if trait_type == "binary":
        _warn(
            logger,
            "For this binary trait, the approximation uses the case-control effective "
            "sample size. The derived standardized BETA is not a direct reconstruction "
            "of the study's logistic-regression log odds ratio and is not a "
            "liability-scale effect.",
        )

    _info(
        logger,
        "Z to BETA/SE conversion:\n"
        "Variants in: {:,}\n"
        "Removed: {:,} ({:.2f}%)\n"
        "  unusable Z score: {:,}\n"
        "  missing frequency: {:,}\n"
        "  degenerate frequency (exactly 0 or 1): {:,} "
        "(at or below 0: {:,}; at or above 1: {:,})\n"
        "  missing effective sample size: {:,}\n"
        "  non-positive effective sample size: {:,}\n"
        "Variants out: {:,}".format(
            n_initial,
            n_removed,
            pct_removed,
            n_z_null,
            n_eaf_null,
            n_eaf_invalid,
            n_eaf_lt_zero,
            n_eaf_gt_one,
            n_neff_null,
            n_neff_invalid,
            n_final,
        ),
    )

    # -------------------------------------------------------
    # FINAL SUMMARY
    # -------------------------------------------------------
    summary = df.select(
        [
            pl.col(beta_col).min().alias("beta_min"),
            pl.col(beta_col).max().alias("beta_max"),
            pl.col(beta_col).mean().alias("beta_mean"),
            pl.col(beta_col).std().alias("beta_std"),
            pl.col(se_col).min().alias("se_min"),
            pl.col(se_col).max().alias("se_max"),
            pl.col(se_col).mean().alias("se_mean"),
            pl.col(se_col).std().alias("se_std"),
            pl.len().alias("total"),
        ]
    ).to_dicts()[0]

    qc_info.update(summary)

    _info(
        logger,
        "BETA and SE summary:\n"
        "BETA -> Min: {} | Max: {} | Mean: {} | Std: {}\n"
        "SE   -> Min: {} | Max: {} | Mean: {} | Std: {} | Total: {:,}".format(
            _fmt(summary["beta_min"]),
            _fmt(summary["beta_max"]),
            _fmt(summary["beta_mean"]),
            _fmt(summary["beta_std"]),
            _fmt(summary["se_min"]),
            _fmt(summary["se_max"]),
            _fmt(summary["se_mean"]),
            _fmt(summary["se_std"]),
            summary["total"],
        ),
    )

    return df, qc_info, sample_column_dict
