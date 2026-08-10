"""Effect-size and standard-error scales (dataset step 07 / chromosome step 03).

Three entry points:

``detect_effect_type(df, effect_col, policies)``
    Pure detection.  Reads the frame, changes nothing, and returns
    ``(decision, evidence)``.  The dataset-level "resolve study properties" step
    calls this ONCE on the full input so every chromosome applies the same
    answer instead of re-deriving it 23 times and throwing it away.

``detect_or_standard_error_scale(...)``
    Pure study-level validation for odds-ratio studies. It compares the Z score
    or two-sided p-value implied when SE is interpreted as log-OR SE with the
    corresponding value implied when SE is interpreted as raw OR-scale SE.

``harmonise_effect_estimates(chromosome, df, sample_column_dict, ...)``
    Applies a decision - the study-level one when it is handed in through
    ``decision=``, otherwise one detected from this chromosome alone - and
    converts odds ratios and any raw OR-scale SE to log-odds units before allele
    orientation. This ordering is required because swapping effect alleles
    changes the sign of log(OR), but does not change its standard error.

GWAS-SSF defines ``standard_error`` as the standard error of the reported effect
and GWAS-VCF requires SE to accompany its effect estimate. For the uncommon case
where uncertainty is supplied on the raw OR scale, the first-order delta method
gives ``SE[log(OR)] = SE[OR] / OR``. See:
https://www.ebi.ac.uk/gwas/docs/summary-statistics-format and
https://github.com/MRCIEU/gwas-vcf-specification.

Everything numeric is a policy loaded from the canonical harmonisation YAML
(see ``policies.py``, group ``effect``).
"""

import math
from typing import Any, Dict, Optional, Tuple

import polars as pl

from postgwas.core.values import format_number as _fmt, optional_text

from .shared.runtime import (
    emit_high_visibility_warning,
    resolve_policies,
    step_context,
)
from .shared.statistics import two_sided_negative_log10_p_from_z


__all__ = [
    "EffectTypeError",
    "POLICY_KEYS",
    "STEP_LABEL",
    "detect_effect_type",
    "detect_or_standard_error_scale",
    "harmonise_effect_estimates",
]


# The label written into the reject file's `reject_step` column, and the
# position of this step in the per-chromosome pipeline.
STEP_LABEL = "03 effect_type"
STEP_NUMBER = 3
STEP_TOTAL = 16
STEP_TITLE = "Effect type (beta or odds ratio)"

# Everything this step reads out of the policy object.  The logger prints these
# with their plain-English explanations at the top of the step.
POLICY_KEYS = [
    "effect.type",
    "effect.or_detection_max_non_positive_fraction",
    "effect.or_detection_require_median_range",
    "effect.or_non_positive",
    "effect.se_scale",
    "effect.se_scale_min_variants",
    "effect.se_scale_min_agreement_fraction",
    "effect.se_scale_min_agreement_margin",
    "effect.cast_strict",
    "effect.verify_per_chromosome",
]


class EffectTypeError(ValueError):
    """The effect column cannot be classified or applied."""


def _normalise_se_scale(value):
    """Return the canonical configured standard-error scale."""
    if value is None:
        return None
    text = str(value).strip().lower()
    if text in ("", "auto", "na", "none"):
        return None
    if text in ("log_odds", "as_given"):
        return text
    raise EffectTypeError(
        "unknown standard-error scale %r - expected 'auto', 'log_odds' or "
        "'as_given'" % (value,)
    )


def _decision_se_scale(decision, configured):
    """Use a study-level SE-scale decision, falling back to explicit policy."""
    value = decision.get("se_scale") if isinstance(decision, dict) else None
    resolved = _normalise_se_scale(value)
    if resolved is not None:
        return resolved
    return _normalise_se_scale(configured)


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def _normalise_decision(value):
    """Return a canonical study-level effect decision."""
    if value is None:
        return None
    if isinstance(value, dict):
        for key in ("effect_type", "decision"):
            if value.get(key):
                value = value[key]
                break
        else:
            return None
    text = str(value).strip().lower()
    if text in ("", "auto", "na", "none"):
        return None
    if text == "beta":
        return "beta"
    if text == "odds_ratio":
        return "odds_ratio"
    raise EffectTypeError(
        "unknown effect type %r - expected 'beta' or 'odds_ratio'" % (value,)
    )


# ---------------------------------------------------------------------------
# 1.  Pure detection  (study level, plan 5.1)
# ---------------------------------------------------------------------------
def detect_effect_type(
    df: pl.DataFrame,
    effect_col: str,
    policies: Optional[Any] = None,
    statistics: Optional[Dict[str, Any]] = None,
) -> Tuple[str, Dict[str, Any]]:
    """Decide whether ``effect_col`` holds beta values or odds ratios.

    Pure: reads the frame, writes nothing, returns ``(decision, evidence)``
    where ``decision`` is ``'beta'`` or ``'odds_ratio'`` and ``evidence`` is the
    dict of statistics behind it (for the log, the QC summary and the manifest).

    Automatic classification is deliberately conservative. An OR candidate
    must have a median in the configured near-one interval and no more than the
    configured combined fraction of zeros and negative values. Distributions
    with contradictory evidence fail instead of risking an invalid logarithm.
    Fractions use only finite numeric values as their denominator. The optional
    ``statistics`` input lets the dataset-level resolver reuse the same canonical
    aggregates without scanning the full study again.
    """
    policies = resolve_policies(policies)

    if not effect_col or effect_col not in df.columns:
        raise EffectTypeError(
            "Effect size column %r is not in the data frame; the columns present "
            "are: %s" % (effect_col, ", ".join(df.columns))
        )

    row = statistics
    if row is None:
        row = df.select(_effect_type_statistic_expressions(effect_col)).to_dicts()[0]

    n_rows = int(row["n_rows"] or 0)
    n_non_null = int(row["n_non_null"] or 0)
    n_usable = int(row["n_usable"] or 0)
    n_negative = int(row["n_negative"] or 0)
    n_zero = int(row["n_zero"] or 0)
    n_non_positive = n_negative + n_zero

    negative_fraction = (n_negative / n_usable) if n_usable else 0.0
    non_positive_fraction = (n_non_positive / n_usable) if n_usable else 0.0
    negative_fraction_all_rows = (n_negative / n_rows) if n_rows else 0.0

    evidence = {
        "effect_col": effect_col,
        "n_rows": n_rows,
        "n_non_null": n_non_null,
        "n_usable": n_usable,
        "n_null_or_non_finite": n_rows - n_usable,
        "negative_value_count": n_negative,
        "zero_value_count": n_zero,
        "non_positive_value_count": n_non_positive,
        "negative_fraction": round(negative_fraction, 6),
        "non_positive_fraction": round(non_positive_fraction, 6),
        "negative_fraction_all_rows": round(negative_fraction_all_rows, 6),
        "min": row["min"],
        "max": row["max"],
        "mean": row["mean"],
        "median": row["median"],
        "std": row["std"],
        "max_non_positive_fraction": (
            policies.effect.or_detection_max_non_positive_fraction
        ),
        "median_range": policies.effect.or_detection_require_median_range,
        "median_range_checked": False,
        "median_range_confirmed": None,
    }

    configured = policies.effect.type
    if configured != "auto":
        evidence["decision_source"] = "config (effect.type)"
        return configured, evidence

    evidence["decision_source"] = "automatic detection"

    if n_usable == 0:
        raise EffectTypeError(
            "Effect size column %r has no usable values (%d rows, %d non-null, "
            "none of them a finite number), so it cannot be classified as beta "
            "or odds ratio.  Set effect.type in the config, or fix the column."
            % (effect_col, n_rows, n_non_null)
        )

    threshold = float(policies.effect.or_detection_max_non_positive_fraction)
    median_range = policies.effect.or_detection_require_median_range
    low, high = float(median_range[0]), float(median_range[1])
    median = row["median"]
    median_confirmed = median is not None and low <= median <= high
    evidence["median_range_checked"] = True
    evidence["median_range_confirmed"] = bool(median_confirmed)

    if non_positive_fraction <= threshold and median_confirmed:
        return "odds_ratio", evidence

    fraction_text = "%.3f%%" % (non_positive_fraction * 100.0)
    threshold_text = "%.3f%%" % (threshold * 100.0)
    range_text = "%s to %s" % (_fmt(low), _fmt(high))
    if non_positive_fraction > threshold and not median_confirmed:
        return "beta", evidence
    if non_positive_fraction > threshold:
        raise EffectTypeError(
            "Effect column %r has conflicting evidence: its median %s is inside "
            "the configured odds-ratio range %s, but %s of usable values are zero "
            "or negative, above the allowed OR contamination limit %s. Specify "
            "effect_type as beta or odds_ratio in the sample sheet."
            % (effect_col, _fmt(median), range_text, fraction_text, threshold_text)
        )
    raise EffectTypeError(
        "Effect column %r is ambiguous: only %s of usable values are zero or "
        "negative, but its median %s is outside the configured odds-ratio range "
        "%s. This may be an all-positive beta column or an unusual selected OR "
        "dataset. Specify effect_type as beta or odds_ratio in the sample sheet."
        % (effect_col, fraction_text, _fmt(median), range_text)
    )


def _effect_type_statistic_expressions(
    effect_col: str,
    prefix: str = "",
):
    """Return the canonical aggregate expressions used by effect detection."""
    values = pl.col(effect_col).cast(pl.Float64, strict=False)
    usable = values.filter(values.is_finite())
    return [
        pl.len().alias(prefix + "n_rows"),
        values.is_not_null().sum().alias(prefix + "n_non_null"),
        values.is_finite().sum().alias(prefix + "n_usable"),
        (usable < 0.0).sum().alias(prefix + "n_negative"),
        (usable == 0.0).sum().alias(prefix + "n_zero"),
        usable.min().alias(prefix + "min"),
        usable.max().alias(prefix + "max"),
        usable.mean().alias(prefix + "mean"),
        usable.median().alias(prefix + "median"),
        usable.std().alias(prefix + "std"),
    ]


def _se_scale_result(method, comparable, informative, log_matches, raw_matches, policies):
    """Summarize one Z- or p-based SE-scale comparison."""
    comparable_count = int(comparable.sum())
    informative_count = int(informative.sum())
    log_count = int((log_matches & informative).sum())
    raw_count = int((raw_matches & informative).sum())
    log_fraction = log_count / informative_count if informative_count else 0.0
    raw_fraction = raw_count / informative_count if informative_count else 0.0
    minimum = int(policies.effect.se_scale_min_variants)
    required = float(policies.effect.se_scale_min_agreement_fraction)
    margin = float(policies.effect.se_scale_min_agreement_margin)

    decision = None
    if informative_count >= minimum:
        if log_fraction >= required and log_fraction - raw_fraction >= margin:
            decision = "log_odds"
        elif raw_fraction >= required and raw_fraction - log_fraction >= margin:
            decision = "as_given"

    return decision, {
        "comparison_source": method,
        "comparable_variants": comparable_count,
        "informative_variants": informative_count,
        "log_odds_matches": log_count,
        "raw_odds_ratio_matches": raw_count,
        "log_odds_agreement_fraction": round(log_fraction, 6),
        "raw_odds_ratio_agreement_fraction": round(raw_fraction, 6),
        "minimum_informative_variants": minimum,
        "minimum_agreement_fraction": required,
        "minimum_agreement_margin": margin,
        "detected_scale": decision,
        "status": (
            "decided"
            if decision is not None
            else "insufficient_informative_variants"
            if informative_count < minimum
            else "ambiguous_or_inconsistent"
        ),
    }


def detect_or_standard_error_scale(
    df: pl.DataFrame,
    effect_col: str,
    se_col: str,
    *,
    z_col: Optional[str] = None,
    pvalue_col: Optional[str] = None,
    pvalue_type: Optional[str] = None,
    policies: Optional[Any] = None,
) -> Tuple[Optional[str], Dict[str, Any]]:
    """Determine whether an OR study's SE is on the log-OR or raw-OR scale.

    A supplied Z score is the strongest check because it directly tests
    ``Z = log(OR) / SE_log``. When no sufficiently informative Z comparison is
    available, the same two candidate Z scores are compared with the reported
    two-sided p-value on the numerically stable ``-log10(p)`` scale. Variants
    for which the two candidate interpretations are indistinguishable under the
    configured tolerance are excluded from the decision rather than counted for
    both sides. A score-test, likelihood-ratio or one-sided p-value is not a
    Wald identity and therefore requires a supplied Z or an explicit SE scale.
    """
    pol = resolve_policies(policies)
    missing = [
        column for column in (effect_col, se_col)
        if not column or column not in df.columns
    ]
    if missing:
        raise EffectTypeError(
            "OR standard-error scale detection requires effect and SE columns; "
            "missing: %s." % ", ".join(repr(column) for column in missing)
        )

    import numpy as np

    selected = [effect_col, se_col]
    for column in (z_col, pvalue_col):
        if column and column in df.columns and column not in selected:
            selected.append(column)
    numeric = df.select([
        pl.col(column).cast(pl.Float64, strict=False).alias(column)
        for column in selected
    ])
    odds_ratio = numeric.get_column(effect_col).to_numpy()
    standard_error = numeric.get_column(se_col).to_numpy()
    usable = (
        np.isfinite(odds_ratio)
        & (odds_ratio > 0.0)
        & np.isfinite(standard_error)
        & (standard_error > 0.0)
    )
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        beta = np.log(odds_ratio)
        candidate_log = beta / standard_error
        candidate_raw = beta / (standard_error / odds_ratio)

    evidence = {
        "rows": int(df.height),
        "effect_column": effect_col,
        "standard_error_column": se_col,
        "usable_or_se_variants": int(usable.sum()),
        "z_attempt": None,
        "pvalue_attempt": None,
    }
    absolute = float(pol.get("validation.beta_se_z_absolute_tolerance"))
    relative = float(pol.get("validation.beta_se_z_relative_tolerance"))

    if z_col and z_col in numeric.columns:
        reference = numeric.get_column(z_col).to_numpy()
        comparable = usable & np.isfinite(reference)
        candidate_difference = np.abs(candidate_log - candidate_raw)
        discrimination_tolerance = absolute + relative * np.maximum(
            np.abs(candidate_log), np.abs(candidate_raw)
        )
        informative = comparable & np.isfinite(candidate_difference) & (
            candidate_difference > discrimination_tolerance
        )
        log_matches = comparable & (
            np.abs(reference - candidate_log)
            <= absolute + relative * np.abs(candidate_log)
        )
        raw_matches = comparable & (
            np.abs(reference - candidate_raw)
            <= absolute + relative * np.abs(candidate_raw)
        )
        decision, attempt = _se_scale_result(
            "supplied_z", comparable, informative, log_matches, raw_matches, pol,
        )
        attempt["absolute_z_tolerance"] = absolute
        attempt["relative_z_tolerance"] = relative
        evidence["z_attempt"] = attempt
        if attempt["informative_variants"] >= attempt["minimum_informative_variants"]:
            evidence.update(attempt)
            return decision, evidence

    if pvalue_col and pvalue_col in numeric.columns and pvalue_type in (
        "raw", "neglog10", "negln"
    ):
        values = numeric.get_column(pvalue_col).to_numpy()
        if pvalue_type == "raw":
            valid_p = np.isfinite(values) & (values > 0.0) & (values <= 1.0)
            with np.errstate(divide="ignore", invalid="ignore"):
                observed = -np.log10(values)
        elif pvalue_type == "neglog10":
            valid_p = np.isfinite(values) & (values >= 0.0)
            observed = values
        else:
            valid_p = np.isfinite(values) & (values >= 0.0)
            observed = values / math.log(10.0)

        try:
            expected_log = two_sided_negative_log10_p_from_z(candidate_log)
            expected_raw = two_sided_negative_log10_p_from_z(candidate_raw)
        except Exception as exc:  # pragma: no cover - environment dependent
            raise EffectTypeError(
                "OR standard-error scale detection from p-values requires scipy, "
                "but it could not be imported (%s). Supply a Z column, install the "
                "harmonisation dependencies, or set effect.se_scale explicitly."
                % exc
            ) from exc

        comparable = usable & valid_p & np.isfinite(observed)
        p_tolerance = float(pol.get("validation.z_pval_tolerance_log10"))
        candidate_difference = np.abs(expected_log - expected_raw)
        informative = comparable & np.isfinite(candidate_difference) & (
            candidate_difference > p_tolerance
        )
        log_matches = comparable & (
            np.abs(observed - expected_log) <= p_tolerance
        )
        raw_matches = comparable & (
            np.abs(observed - expected_raw) <= p_tolerance
        )
        decision, attempt = _se_scale_result(
            "reported_two_sided_pvalue",
            comparable,
            informative,
            log_matches,
            raw_matches,
            pol,
        )
        attempt["pvalue_type"] = pvalue_type
        attempt["pvalue_tolerance_log10"] = p_tolerance
        evidence["pvalue_attempt"] = attempt
        evidence.update(attempt)
        return decision, evidence

    if evidence["z_attempt"] is not None:
        evidence.update(evidence["z_attempt"])
        return None, evidence

    evidence.update({
        "comparison_source": None,
        "comparable_variants": 0,
        "informative_variants": 0,
        "detected_scale": None,
        "status": "no_usable_z_or_pvalue_comparison",
        "minimum_informative_variants": int(
            pol.effect.se_scale_min_variants
        ),
    })
    return None, evidence


def _decision_sentence(decision, evidence):
    """The one-line evidence string handed to ``logger.decide``."""
    return (
        "%.3f%% of %s usable values are zero or negative (%s zero, %s negative; "
        "%s null or non-finite); median %s, minimum %s, maximum %s"
        % (
            evidence.get("non_positive_fraction", 0.0) * 100.0,
            "{:,}".format(evidence.get("n_usable", 0)),
            "{:,}".format(evidence.get("zero_value_count", 0)),
            "{:,}".format(evidence.get("negative_value_count", 0)),
            "{:,}".format(evidence.get("n_null_or_non_finite", 0)),
            _fmt(evidence.get("median")),
            _fmt(evidence.get("min")),
            _fmt(evidence.get("max")),
        )
    )


def _automatic_inference_warning(effect_type, evidence):
    """Explain the scientific limitation whenever an inferred type is applied."""
    processing = (
        "Values will be converted using BETA = ln(effect)."
        if effect_type == "odds_ratio"
        else "Values will be treated directly as beta coefficients."
    )
    return (
        "AUTOMATIC EFFECT-TYPE INFERENCE: column %r was inferred as %r from "
        "its value distribution (%.3f%% of %s usable values are zero or "
        "negative; median %s). This heuristic cannot prove whether the column "
        "contains beta coefficients or odds ratios; an all-positive beta column "
        "can resemble an odds-ratio column. %s For scientifically reliable "
        "harmonisation, set the sample-sheet effect_type field to 'beta' or "
        "'odds_ratio', or set effect.type explicitly in YAML."
        % (
            evidence.get("effect_col"),
            effect_type,
            float(evidence.get("non_positive_fraction") or 0.0) * 100.0,
            "{:,}".format(int(evidence.get("n_usable") or 0)),
            _fmt(evidence.get("median")),
            processing,
        )
    )


# ---------------------------------------------------------------------------
# 2.  Application  (per chromosome)
# ---------------------------------------------------------------------------
def harmonise_effect_estimates(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    logger: Optional[Any] = None,
    policies: Optional[Any] = None,
    rejects: Optional[Any] = None,
    decision: Optional[Any] = None,
    step_number: Optional[int] = None,
    step_total: Optional[int] = None,
) -> Tuple[pl.DataFrame, dict, dict]:
    """Apply the effect-type decision and convert odds ratios to log-odds.

    Parameters
    ----------
    chromosome, df, sample_column_dict
        Resolved input columns used by downstream steps.
    logger
        A ``pipeline_logging.PipelineLogger``. ``None`` disables logging entirely;
        this module no longer writes a log file of its own.
    policies
        A ``policies.Policies``.  ``None`` means all defaults, which reproduces
        the behaviour this module had before it was migrated.
    rejects
        A ``rejects.RejectCollector``.  Only used when
        ``effect.or_non_positive`` is ``reject``.
    decision
        The study-level answer from the dataset step ('beta' / 'odds_ratio', or
        a dict carrying one).  When given it is authoritative: per-chromosome
        detection can only warn, never override it.  ``None`` falls back to
        detecting from this chromosome alone.

    Returns ``(df, qc_info, sample_column_dict)`` exactly as before.
    """
    policies = resolve_policies(policies)
    number = STEP_NUMBER if step_number is None else step_number
    total_steps = STEP_TOTAL if step_total is None else step_total

    configured_se_scale = str(policies.effect.se_scale)
    resolved_se_scale = _decision_se_scale(decision, configured_se_scale)
    qc_info = {
        "initial_variants": df.height,
        "se_input_scale": resolved_se_scale or "not_decided",
        "se_scale_decision_source": (
            "study-level decision"
            if isinstance(decision, dict) and decision.get("se_scale")
            else "policy"
            if configured_se_scale != "auto"
            else "not_decided"
        ),
    }

    effect_col = sample_column_dict.get("beta_or_col", None)
    if not effect_col or effect_col not in df.columns:
        raise EffectTypeError(
            "Effect size column ('beta_or_col') is missing from the config or "
            "not present in the data frame for chromosome %s." % (chromosome,)
        )

    with step_context(
        logger,
        None,
        number=number,
        total=total_steps,
        title=STEP_TITLE,
        operation="effect_type.harmonise_effect_estimates",
        rows_in=df.height,
        policy_keys=POLICY_KEYS,
    ) as ctx:
        # ------------------------------------------------------------------
        # Cast to a number.  effect.cast_strict is true by default, which is
        # what this module did before: one unparseable cell aborts the
        # chromosome.  false turns those cells into nulls and reports them.
        # ------------------------------------------------------------------
        strict = bool(policies.effect.cast_strict)
        qc_info["effect_cast_strict"] = strict
        if strict:
            df = df.with_columns(pl.col(effect_col).cast(pl.Float64, strict=True))
        else:
            n_unparseable = int(
                df.select(
                    (
                        pl.col(effect_col).is_not_null()
                        & pl.col(effect_col).cast(pl.Float64, strict=False).is_null()
                    )
                    .sum()
                    .alias("n")
                ).item()
                or 0
            )
            before = df.height
            df = df.with_columns(pl.col(effect_col).cast(pl.Float64, strict=False))
            qc_info["effect_unparseable_values"] = n_unparseable
            if n_unparseable:
                ctx.qc(
                    "effect size is a number",
                    "%s values in '%s' are not numbers. Policy 'effect.cast_strict' "
                    "is false, so they were blanked out instead of stopping the "
                    "chromosome." % ("{:,}".format(n_unparseable), effect_col),
                    before,
                    df.height,
                    changed=n_unparseable,
                    warn=True,
                    step=STEP_LABEL,
                )

        # ------------------------------------------------------------------
        # Which decision applies?
        # ------------------------------------------------------------------
        study_decision = _normalise_decision(decision)
        evidence = None
        if study_decision is not None:
            effect_type = study_decision
            decision_source = "study-level decision (made once, before the split)"
            if policies.effect.verify_per_chromosome:
                # Verification may only warn.  A chromosome whose own values
                # cannot be classified must not fail a run whose effect type is
                # already known.
                try:
                    local, evidence = detect_effect_type(
                        df, effect_col,
                        policies.with_overrides({"effect.type": "auto"}),
                    )
                except EffectTypeError as exc:
                    local = None
                    ctx.warn(
                        "Re-checking the effect type on this chromosome alone was "
                        "inconclusive, which changes nothing: the study-level "
                        "decision '%s' was applied. %s" % (effect_type, exc)
                    )
                else:
                    if local != effect_type:
                        ctx.warn(
                            "Re-checking on this chromosome alone suggests '%s', but "
                            "the study-level decision '%s' is authoritative and was "
                            "applied. %s"
                            % (local, effect_type, _decision_sentence(local, evidence))
                        )
                    else:
                        ctx.info(
                            "Re-checked on this chromosome: agrees with the "
                            "study-level decision '%s'." % (effect_type,)
                        )
                qc_info["effect_type_chromosome_check"] = local
        else:
            effect_type, evidence = detect_effect_type(df, effect_col, policies)
            decision_source = evidence.get("decision_source", "automatic detection")
            if decision_source == "automatic detection":
                warning = _automatic_inference_warning(effect_type, evidence)
                qc_info["effect_type_inference_warning"] = warning
                emit_high_visibility_warning(logger, warning)

        if evidence is None:
            evidence = {"effect_col": effect_col}
            qc_info.update({"effect_col": effect_col})
        else:
            qc_info.update(
                {
                    "effect_col": effect_col,
                    "negative_value_count": evidence.get("negative_value_count"),
                    "negative_fraction": evidence.get("negative_fraction"),
                    "non_positive_value_count": evidence.get(
                        "non_positive_value_count"
                    ),
                    "non_positive_fraction": evidence.get("non_positive_fraction"),
                    "negative_fraction_all_rows": evidence.get(
                        "negative_fraction_all_rows"
                    ),
                    "usable_effect_values": evidence.get("n_usable"),
                    "effect_median": evidence.get("median"),
                    "effect_min": evidence.get("min"),
                    "effect_max": evidence.get("max"),
                }
            )

        qc_info["effect_type"] = effect_type
        qc_info["effect_decision_source"] = decision_source
        sample_column_dict["effect_type"] = effect_type

        ctx.decide(
            "Effect column '%s'" % (effect_col,),
            _decision_sentence(effect_type, evidence)
            if evidence.get("n_usable") is not None
            else decision_source,
            "%s (%s)"
            % (
                "beta" if effect_type == "beta" else "odds ratio",
                decision_source,
            ),
        )
        # ------------------------------------------------------------------
        # 3A.  BETA - nothing to convert.
        # ------------------------------------------------------------------
        if effect_type == "beta":
            sample_column_dict["beta_col"] = effect_col
            sample_column_dict["beta_or_col"] = effect_col
            sample_column_dict["harmonised_effect_type"] = "beta"
            beta_col = effect_col

        # ------------------------------------------------------------------
        # 3B.  ODDS RATIO - convert to log-odds.
        # ------------------------------------------------------------------
        else:
            pre_stats = df.select(
                [
                    pl.col(effect_col).min().alias("min_pre"),
                    pl.col(effect_col).max().alias("max_pre"),
                    pl.col(effect_col).mean().alias("mean_pre"),
                    pl.col(effect_col).std().alias("std_pre"),
                    (pl.col(effect_col) == 0.0).sum().alias("n_zero"),
                    (pl.col(effect_col) < 0.0).sum().alias("n_negative"),
                ]
            ).to_dicts()[0]
            qc_info.update(
                {
                    "pre_conversion_min": pre_stats["min_pre"],
                    "pre_conversion_max": pre_stats["max_pre"],
                    "pre_conversion_mean": pre_stats["mean_pre"],
                    "pre_conversion_std": pre_stats["std_pre"],
                }
            )
            ctx.info(
                "Odds ratios before conversion: min %s, max %s, mean %s."
                % (
                    _fmt(pre_stats["min_pre"]),
                    _fmt(pre_stats["max_pre"]),
                    _fmt(pre_stats["mean_pre"]),
                )
            )

            # Non-positive odds ratios have no logarithm. Count zeros and
            # negatives separately so the QC report describes both conditions.
            n_zero = int(pre_stats["n_zero"] or 0)
            n_negative = int(pre_stats["n_negative"] or 0)
            n_non_positive = n_zero + n_negative
            action = policies.effect.or_non_positive
            qc_info.update(
                {
                    "or_zero_count": n_zero,
                    "or_negative_count": n_negative,
                    "or_non_positive_count": n_non_positive,
                    "or_non_positive_action": action,
                }
            )

            detail = "%s exactly zero, %s negative" % (
                "{:,}".format(n_zero),
                "{:,}".format(n_negative),
            )
            plain = (
                "%s odds ratios are zero or below (%s). The logarithm of a "
                "number that is not positive does not exist."
                % ("{:,}".format(n_non_positive), detail)
            )

            if n_non_positive and action == "fail":
                raise EffectTypeError(
                    "%s Policy 'effect.or_non_positive' is 'fail'. Chromosome %s."
                    % (plain, chromosome)
                )

            if n_non_positive and action == "reject":
                if rejects is None:
                    ctx.warn(
                        "Policy 'effect.or_non_positive' is 'reject' but no reject "
                        "collector was supplied, so the %s non-positive odds ratios "
                        "were blanked out instead of being rejected and recorded."
                        % ("{:,}".format(n_non_positive),)
                    )
                    action = "null"
                else:
                    # reject() logs its own before/after through logger.qc().
                    # fill_null(False) on purpose: reject() treats a null mask
                    # value as "reject", and a missing odds ratio is not a
                    # non-positive one - it is dealt with further down the line.
                    df = rejects.reject(
                        df,
                        (pl.col(effect_col) <= 0.0).fill_null(False),
                        STEP_LABEL,
                        "effect_non_positive_or",
                        detail=detail,
                    )

            if action == "keep":
                # No blanking: log(0) is -inf and log(negative) is NaN, and the
                # effect-statistics gate downstream is what removes them.
                df = df.with_columns(pl.col(effect_col).log().alias("beta"))
                if n_non_positive:
                    ctx.qc(
                        "non-positive odds ratios",
                        "%s Policy 'effect.or_non_positive' is 'keep', so they were "
                        "left in place; the logarithm is infinite or not a number "
                        "for those variants." % (plain,),
                        df.height,
                        df.height,
                        changed=n_non_positive,
                        warn=True,
                        step=STEP_LABEL,
                    )
            else:
                # 'null' (the default, and what this module always did) and the
                # remainder of 'reject', where the rows are already gone.
                df = df.with_columns(
                    pl.when(pl.col(effect_col) > 0.0)
                    .then(pl.col(effect_col).log())
                    .otherwise(None)
                    .alias("beta")
                )
                if n_non_positive and action == "null":
                    ctx.qc(
                        "non-positive odds ratios",
                        "%s Policy 'effect.or_non_positive' is 'null', so their "
                        "effect size was blanked out and the variants kept."
                        % (plain,),
                        df.height,
                        df.height,
                        changed=n_non_positive,
                        warn=True,
                        step=STEP_LABEL,
                    )

            sample_column_dict["beta_col"] = "beta"
            sample_column_dict["beta_or_col"] = "beta"
            sample_column_dict["harmonised_effect_type"] = "beta"
            beta_col = "beta"
            qc_info["conversion"] = "OR_to_Beta_log_transform_applied"

            # The standard error, if it arrived on the odds-ratio scale.
            se_col = optional_text(sample_column_dict.get("se_col"))
            rescaled = False
            if se_col is not None and se_col in df.columns and resolved_se_scale is None:
                raise EffectTypeError(
                    "Chromosome %s has odds ratios and a supplied standard-error "
                    "column, but the study-level SE-scale decision is missing. Run "
                    "dataset-level study-property resolution first, or explicitly set "
                    "effect.se_scale to 'log_odds' or 'as_given'."
                    % chromosome
                )
            if resolved_se_scale == "as_given":
                if se_col is not None and se_col in df.columns:
                    df = df.with_columns(
                        pl.when(pl.col(effect_col) > 0.0)
                        .then(
                            pl.col(se_col).cast(pl.Float64, strict=False)
                            / pl.col(effect_col)
                        )
                        .otherwise(None)
                        .alias(se_col)
                    )
                    rescaled = True
                    ctx.info(
                        "The study-level effect.se_scale decision is '%s', so the standard error "
                        "'%s' was divided by the odds ratio to put it on the "
                        "log-odds scale." % (resolved_se_scale, se_col)
                    )
                else:
                    ctx.warn(
                        "The study-level effect.se_scale decision is '%s', which asks for the "
                        "standard error to be divided by the odds ratio, but no "
                        "standard-error column is available."
                        % (resolved_se_scale,)
                    )
            qc_info["se_rescaled_by_or"] = rescaled

            post_stats = df.select(
                [
                    pl.col("beta").min().alias("min_post"),
                    pl.col("beta").max().alias("max_post"),
                    pl.col("beta").mean().alias("mean_post"),
                    pl.col("beta").std().alias("std_post"),
                ]
            ).to_dicts()[0]
            qc_info.update(
                {
                    "post_conversion_min": post_stats["min_post"],
                    "post_conversion_max": post_stats["max_post"],
                    "post_conversion_mean": post_stats["mean_post"],
                    "post_conversion_std": post_stats["std_post"],
                }
            )
            ctx.info(
                "Log-odds after conversion: min %s, max %s, mean %s."
                % (
                    _fmt(post_stats["min_post"]),
                    _fmt(post_stats["max_post"]),
                    _fmt(post_stats["mean_post"]),
                )
            )
            stats = {
                "min_beta": post_stats["min_post"],
                "max_beta": post_stats["max_post"],
                "mean_beta": post_stats["mean_post"],
                "std_beta": post_stats["std_post"],
            }

        # ------------------------------------------------------------------
        # Summary statistics, always on the column that now holds the beta.
        # (They used to be computed on the untransformed odds ratios, so
        # min_beta / max_beta / mean_beta described odds ratios.)
        # ------------------------------------------------------------------
        if effect_type == "beta":
            stats = df.select(
                [
                    pl.col(beta_col).min().alias("min_beta"),
                    pl.col(beta_col).max().alias("max_beta"),
                    pl.col(beta_col).mean().alias("mean_beta"),
                    pl.col(beta_col).std().alias("std_beta"),
                ]
            ).to_dicts()[0]
        qc_info.update(
            {
                "beta_col": beta_col,
                "min_beta": stats["min_beta"],
                "max_beta": stats["max_beta"],
                "mean_beta": stats["mean_beta"],
                "std_beta": stats["std_beta"],
            }
        )
        ctx.info(
            "Effect sizes in '%s': min %s, max %s, mean %s."
            % (
                beta_col,
                _fmt(stats["min_beta"]),
                _fmt(stats["max_beta"]),
                _fmt(stats["mean_beta"]),
            )
        )

        qc_info["final_total"] = df.height
        ctx.set_rows(df.height)

    return df, qc_info, sample_column_dict
