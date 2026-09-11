"""
study_properties.py — decide the study-wide properties once, on the whole file.

Plan v3, section 5.1 (dataset-level step 07, immediately before the chromosome
split).

Effect type, odds-ratio SE scale, p-value scale and whether the frequency column
is MAF-like are properties of the STUDY, not of a chromosome. This step decides
them once, on the full input, records the evidence as DECIDE lines, and hands
the answers to every worker as data. A MAF-like distribution is only a
suspicion; chromosome reference checks are consolidated into the final dataset
frequency type after the workers complete.

**Why the estimator has to be sound.**  Optional per-chromosome re-verification
of the effect and p-value decisions is off by default, and the study-level answer
is always the one applied.  There is therefore no second line of defence: every
statistic here is computed over the NON-NULL subset, and the non-null count is
stated in the log, because a denominator that includes nulls is exactly how a
mostly-empty effect column gets misread as odds ratios and log-transformed
(plan 5.1, defect table).

The classification itself is delegated to the pure detector functions exposed
by the per-step modules: ``detect_effect_type`` in ``effect_type.py`` and
``detect_or_standard_error_scale`` in ``effect_type.py`` and
``detect_p_value_type`` in ``p_values.py``. This module does not carry a second,
potentially divergent implementation of those scientific decisions.

Explicit sample-sheet declarations take precedence, but the corresponding
canonical detector still runs as a cross-check. A disagreement is handled by
``validation.declaration_mismatch_action``. Explicit YAML policy values retain
their existing authoritative, no-cross-check behaviour.

Python 3.8 compatible.  Imports polars and the standard library at module level.
"""

import math
from typing import Any, Dict, Optional

import polars as pl

from postgwas.core.dataframes import numeric_column
from postgwas.core.values import (
    format_count, format_number, format_percentage, optional_text,
)

from .shared.runtime import (
    NullLogger,
    configured_column,
    emit_high_visibility_warning,
    format_inference_warning,
    resolve_policies,
)
from .shared.statistics import valid_frequency_mask

__all__ = [
    "resolve_study_properties",
    "finalise_eaf_decision_from_chromosomes",
    "StudyPropertyError",
    "EFFECT_COLUMN_KEYS",
    "PVALUE_COLUMN_KEYS",
    "EAF_COLUMN_KEYS",
]


#: sample_column_dict keys that may hold the raw effect-size column.
EFFECT_COLUMN_KEYS = ("beta_or_col", "beta_col")

#: sample_column_dict keys that may hold the p-value column.
PVALUE_COLUMN_KEYS = ("pval_col",)

#: sample_column_dict keys that may hold the study's own frequency column.
EAF_COLUMN_KEYS = ("eaf_col",)

_DATASET_STEP_LABEL = "07 resolve_study_properties"

_POLICY_KEYS = [
    "effect.type",
    "effect.or_detection_max_non_positive_fraction",
    "effect.or_detection_require_median_range",
    "effect.verify_per_chromosome",
    "effect_from_z.beta_z_sign_mismatch",
    "effect_from_z.max_beta_z_sign_mismatch_fraction",
    "effect.se_scale",
    "effect.se_scale_min_variants",
    "effect.se_scale_min_agreement_fraction",
    "effect.se_scale_min_agreement_margin",
    "pvalue.type",
    "pvalue.mlogp_detect_threshold",
    "pvalue.mlogp_detect_proportion",
    "pvalue.mlogp_detect_min_count",
    "pvalue.mlogp_expected_median",
    "pvalue.mlogp_median_tolerance",
    "pvalue.max_negative_fraction",
    "pvalue.verify_per_chromosome",
    "eaf.maf_decision_cutoff",
    "validation.beta_se_z_absolute_tolerance",
    "validation.beta_se_z_relative_tolerance",
    "validation.z_pval_tolerance_log10",
    "validation.declaration_mismatch_action",
]


class StudyPropertyError(Exception):
    """A study-wide property could not be decided with confidence.

    Carries a plain-English description of what was ambiguous and which policy
    the user should set explicitly.  A typed exception rather than sys.exit() so
    the dataset-level retry can catch it (plan Part 7).
    """


_column = configured_column
_numeric = numeric_column


def _clean(value):
    # type: (Any) -> Any
    """Make a polars scalar safe to put in a JSON manifest."""
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return value
    if math.isnan(number) or math.isinf(number):
        return None
    return number


_count = format_count
_percent = format_percentage


# =============================================================================
# Evidence — every statistic over the NON-NULL subset
# =============================================================================


def _effect_evidence(df, column, statistics=None):
    # type: (pl.DataFrame, str) -> Dict[str, Any]
    """Statistics for the effect column, all denominators non-null."""
    row = statistics or df.select(_effect_evidence_expressions(column)).to_dicts()[0]

    non_null = int(row["non_null"] or 0)
    negative = int(row["negative"] or 0)
    evidence = {
        "column": column,
        "rows": int(row["rows"] or 0),
        "non_null": non_null,
        "null": int(row["rows"] or 0) - non_null,
        "unparseable": int(row["non_null_before_cast"] or 0) - non_null,
        "negative": negative,
        "negative_fraction_of_non_null": (float(negative) / non_null) if non_null else 0.0,
        "zero": int(row["zero"] or 0),
        "non_positive": int(row["non_positive"] or 0),
        "min": _clean(row["min"]),
        "max": _clean(row["max"]),
        "mean": _clean(row["mean"]),
        "median": _clean(row["median"]),
    }
    return evidence


def _pvalue_evidence(df, column, detect_threshold, statistics=None):
    # type: (pl.DataFrame, str, float) -> Dict[str, Any]
    """Statistics for the p-value column, all denominators non-null."""
    row = statistics or df.select(
        _pvalue_evidence_expressions(column, detect_threshold)
    ).to_dicts()[0]

    non_null = int(row["non_null"] or 0)
    above = int(row["above_threshold"] or 0)
    evidence = {
        "column": column,
        "rows": int(row["rows"] or 0),
        "non_null": non_null,
        "null": int(row["rows"] or 0) - non_null,
        "unparseable": int(row["non_null_before_cast"] or 0) - non_null,
        "detect_threshold": float(detect_threshold),
        "above_threshold": above,
        "above_threshold_fraction_of_non_null": (float(above) / non_null) if non_null else 0.0,
        "negative": int(row["negative"] or 0),
        "min": _clean(row["min"]),
        "max": _clean(row["max"]),
        "mean": _clean(row["mean"]),
        "median": _clean(row["median"]),
    }
    return evidence


def _eaf_evidence(df, column, statistics=None):
    # type: (pl.DataFrame, str) -> Dict[str, Any]
    """Statistics for the frequency column using only usable probabilities.

    The MAF-like fraction is ``usable values <= 0.5 / usable values``. Missing,
    unparseable, non-finite and finite out-of-range values are excluded from
    both sides of that scientific comparison and reported separately.
    """
    row = statistics or df.select(_eaf_evidence_expressions(column)).to_dicts()[0]

    rows = int(row["rows"] or 0)
    raw_non_null = int(row["non_null_before_cast"] or 0)
    non_null = int(row["non_null"] or 0)
    usable = int(row["usable"] or 0)
    low = int(row["le_0_5"] or 0)
    evidence = {
        "column": column,
        "rows": rows,
        "input_missing": rows - raw_non_null,
        "non_null": non_null,
        "null": rows - non_null,
        "unparseable": raw_non_null - non_null,
        "usable": usable,
        "unusable": rows - usable,
        "non_finite": int(row["non_finite"] or 0),
        "out_of_range": int(row["out_of_range"] or 0),
        "at_or_below_0_5": low,
        "low_fraction_of_usable": (float(low) / usable) if usable else 0.0,
        "exactly_zero": int(row["eq_0"] or 0),
        "exactly_one": int(row["eq_1"] or 0),
        "min": _clean(row["min"]),
        "max": _clean(row["max"]),
        "mean": _clean(row["mean"]),
        "median": _clean(row["median"]),
    }
    return evidence


def _effect_evidence_expressions(column, prefix=""):
    value = pl.col(column).cast(pl.Float64, strict=False)
    return [
        pl.len().alias(prefix + "rows"),
        pl.col(column).is_not_null().sum().alias(prefix + "non_null_before_cast"),
        value.is_not_null().sum().alias(prefix + "non_null"),
        (value < 0.0).sum().alias(prefix + "negative"),
        (value == 0.0).sum().alias(prefix + "zero"),
        (value <= 0.0).sum().alias(prefix + "non_positive"),
        value.min().alias(prefix + "min"),
        value.max().alias(prefix + "max"),
        value.mean().alias(prefix + "mean"),
        value.median().alias(prefix + "median"),
    ]


def _pvalue_evidence_expressions(column, detect_threshold, prefix=""):
    value = pl.col(column).cast(pl.Float64, strict=False)
    return [
        pl.len().alias(prefix + "rows"),
        pl.col(column).is_not_null().sum().alias(prefix + "non_null_before_cast"),
        value.is_not_null().sum().alias(prefix + "non_null"),
        (value > float(detect_threshold)).sum().alias(prefix + "above_threshold"),
        (value < 0.0).sum().alias(prefix + "negative"),
        value.min().alias(prefix + "min"),
        value.max().alias(prefix + "max"),
        value.mean().alias(prefix + "mean"),
        value.median().alias(prefix + "median"),
    ]


def _eaf_evidence_expressions(column, prefix=""):
    value = pl.col(column).cast(pl.Float64, strict=False)
    usable = valid_frequency_mask(value)
    return [
        pl.len().alias(prefix + "rows"),
        pl.col(column).is_not_null().sum().alias(prefix + "non_null_before_cast"),
        value.is_not_null().sum().alias(prefix + "non_null"),
        usable.sum().alias(prefix + "usable"),
        (value.is_not_null() & ~value.is_finite()).sum().alias(prefix + "non_finite"),
        (
            value.is_not_null()
            & value.is_finite()
            & ~value.is_between(0.0, 1.0)
        ).sum().alias(prefix + "out_of_range"),
        (usable & (value <= 0.5)).sum().alias(prefix + "le_0_5"),
        (usable & (value == 0.0)).sum().alias(prefix + "eq_0"),
        (usable & (value == 1.0)).sum().alias(prefix + "eq_1"),
        pl.when(usable).then(value).min().alias(prefix + "min"),
        pl.when(usable).then(value).max().alias(prefix + "max"),
        pl.when(usable).then(value).mean().alias(prefix + "mean"),
        pl.when(usable).then(value).median().alias(prefix + "median"),
    ]


def _prefixed_statistics(row, prefix):
    return {
        key[len(prefix):]: value
        for key, value in row.items()
        if key.startswith(prefix)
    }


def _collect_study_statistics(df, sample_column_dict, pol):
    """Collect all study-property evidence in one full-frame aggregation."""
    from .effect_type import _effect_type_statistic_expressions
    from .effect_from_z import _beta_z_sign_statistic_expressions
    from .p_values import _pvalue_type_statistic_expressions

    effect_col = _column(sample_column_dict, EFFECT_COLUMN_KEYS, df)
    pvalue_col = _column(sample_column_dict, PVALUE_COLUMN_KEYS, df)
    eaf_col = _column(sample_column_dict, EAF_COLUMN_KEYS, df)
    z_col = _column(sample_column_dict, ("imp_z_col", "z_col"), df)
    se_col = _column(sample_column_dict, ("se_col",), df)
    threshold = float(pol.get("pvalue.mlogp_detect_threshold"))
    expressions = []
    if effect_col is not None:
        expressions.extend(
            _effect_evidence_expressions(effect_col, "effect_evidence__")
        )
        expressions.extend(
            _effect_type_statistic_expressions(effect_col, "effect_detector__")
        )
        if z_col is not None:
            expressions.extend(
                _beta_z_sign_statistic_expressions(
                    effect_col,
                    z_col,
                    se_col,
                    "beta_z_sign__",
                )
            )
    if pvalue_col is not None:
        expressions.extend(
            _pvalue_evidence_expressions(pvalue_col, threshold, "pvalue_evidence__")
        )
        expressions.extend(
            _pvalue_type_statistic_expressions(
                pvalue_col, threshold, "pvalue_detector__",
            )
        )
    if eaf_col is not None:
        expressions.extend(_eaf_evidence_expressions(eaf_col, "eaf_evidence__"))
    if not expressions:
        return {}
    row = df.select(expressions).to_dicts()[0]
    return {
        "effect_evidence": _prefixed_statistics(row, "effect_evidence__"),
        "effect_detector": _prefixed_statistics(row, "effect_detector__"),
        "pvalue_evidence": _prefixed_statistics(row, "pvalue_evidence__"),
        "pvalue_detector": _prefixed_statistics(row, "pvalue_detector__"),
        "eaf_evidence": _prefixed_statistics(row, "eaf_evidence__"),
        "beta_z_sign": _prefixed_statistics(row, "beta_z_sign__"),
    }


# =============================================================================
# The three decisions
# =============================================================================


def _decide_effect_type(
    df, sample_column_dict, pol, ctx, decisions, statistics, logger=None,
):
    configured = pol.get("effect.type")
    declared = optional_text(sample_column_dict.get("declared_effect_type"))
    declared = str(declared).lower() if declared is not None else None
    declared = declared if declared in ("beta", "odds_ratio") else None
    column = _column(sample_column_dict, EFFECT_COLUMN_KEYS, df)

    decisions["effect_column"] = column

    if configured in ("beta", "odds_ratio") and declared is None:
        decisions["effect_type"] = configured
        decisions["effect_type_source"] = "policy"
        decisions["effect_evidence"] = (
            _effect_evidence(df, column, statistics.get("effect_evidence"))
            if column is not None else {}
        )
        ctx.decide(
            "Effect size scale",
            "policy effect.type is set explicitly, so no detection was run",
            "%s. Applied to all chromosomes."
            % ("odds ratio" if configured == "odds_ratio" else "beta"),
        )
        return

    if column is None:
        decisions["effect_type"] = None
        decisions["effect_type_source"] = "no_effect_column"
        decisions["effect_evidence"] = {}
        ctx.skip(
            "No effect-size column is configured, so the effect-size scale was not "
            "decided here. It is derived from the Z score per chromosome instead."
        )
        return

    evidence = _effect_evidence(df, column, statistics.get("effect_evidence"))
    decisions["effect_evidence"] = evidence

    if evidence["non_null"] == 0:
        raise StudyPropertyError(
            "The effect-size column '%s' has no usable numeric value in any of its %s "
            "rows, so beta cannot be told from odds ratio. Check that the column name "
            "in the config is the right one, or set effect.type to 'beta' or "
            "'odds_ratio' explicitly." % (column, _count(evidence["rows"]))
        )

    from .effect_type import _automatic_inference_warning, detect_effect_type

    detector_policies = (
        pol.with_overrides({"effect.type": "auto"}) if declared is not None else pol
    )
    try:
        verdict, detector_evidence = detect_effect_type(
            df, column, detector_policies,
            statistics=statistics.get("effect_detector"),
        )
    except ValueError as exc:
        if declared is not None:
            decisions["effect_type"] = declared
            decisions["effect_type_source"] = "sample_sheet"
            decisions["effect_type_detected"] = None
            decisions["effect_type_matches_declaration"] = None
            ctx.warn(
                "The effect_type declaration %r remains authoritative, but automatic "
                "detection could not independently classify column %r: %s"
                % (declared, column, exc)
            )
            return
        raise StudyPropertyError(
            "The effect-size detector refused to classify column '%s': %s"
            % (column, exc)
        ) from exc
    if detector_evidence:
        evidence["detector_stats"] = _detector_stats(detector_evidence)

    decisions["effect_type_detected"] = verdict
    if declared is not None:
        decisions["effect_type"] = declared
        decisions["effect_type_source"] = "sample_sheet"
        decisions["effect_type_matches_declaration"] = verdict == declared
        if verdict != declared:
            action = str(pol.get("validation.declaration_mismatch_action"))
            message = (
                "DECLARATION MISMATCH: sample sheet declares effect_type=%r, but "
                "automatic detection found %r for column %r. Policy "
                "validation.declaration_mismatch_action=%r; the declaration remains "
                "authoritative."
                % (declared, verdict, column, action)
            )
            decisions.setdefault("declaration_mismatches", []).append(message)
            emit_high_visibility_warning(logger, message)
            if action == "fail":
                raise StudyPropertyError(message)
    else:
        decisions["effect_type"] = verdict
        decisions["effect_type_source"] = "detector"
        warning = _automatic_inference_warning(verdict, detector_evidence)
        decisions["effect_type_inference_warning"] = warning
        emit_high_visibility_warning(
            logger, warning,
            screen_message=_automatic_inference_warning(verdict, detector_evidence, screen=True),
        )

    ctx.decide(
        "Effect column '%s'" % column,
        {
            "non-null values": _count(evidence["non_null"]),
            "negative": "%s (%s)" % (_count(evidence["negative"]),
                                     _percent(evidence["negative"], evidence["non_null"])),
            "zero": _count(evidence["zero"]),
            "non-positive": "%s (%s)" % (
                _count(evidence["non_positive"]),
                _percent(evidence["non_positive"], evidence["non_null"]),
            ),
            "median": evidence["median"],
            "minimum": evidence["min"],
            "decided by": decisions["effect_type_source"],
        },
        "%s. Applied to all chromosomes."
        % ("odds ratio (log-transformed to beta per chromosome)"
           if decisions["effect_type"] == "odds_ratio" else "beta"),
    )


def _detector_stats(extra):
    # type: (Dict[str, Any]) -> Dict[str, Any]
    """Keep a detector's own statistics, but only the manifest-safe scalars."""
    kept = {}
    for key, value in extra.items():
        if isinstance(value, (int, float, str, bool)) or value is None:
            kept[str(key)] = _clean(value) if isinstance(value, float) else value
    return kept


def _decide_beta_z_sign_consistency(
    df, sample_column_dict, pol, ctx, decisions, statistics, logger=None,
):
    """Validate signed BETA/Z agreement once across the complete study."""
    effect_col = _column(sample_column_dict, EFFECT_COLUMN_KEYS, df)
    z_col = _column(sample_column_dict, ("imp_z_col", "z_col"), df)
    se_col = _column(sample_column_dict, ("se_col",), df)
    if effect_col is None or z_col is None:
        decisions["beta_z_sign_consistency"] = {
            "applies": False,
            "reason": "effect or signed Z column is not supplied",
        }
        return

    from .effect_from_z import (
        BetaZSignMismatchError,
        assess_beta_z_sign_discordance,
    )

    try:
        assessment = assess_beta_z_sign_discordance(
            statistics.get("beta_z_sign") or {},
            effect_type=decisions["effect_type"],
            effect_col=effect_col,
            z_col=z_col,
            se_col=se_col,
            policies=pol,
        )
    except BetaZSignMismatchError as exc:
        message = "BETA/Z SIGN DISCORDANCE LIMIT EXCEEDED: {}".format(exc)
        emit_high_visibility_warning(logger, message)
        raise StudyPropertyError(str(exc)) from exc

    decisions["beta_z_sign_consistency"] = assessment
    discordant = assessment["discordant_rows"]
    comparable = assessment["comparable_missing_se_rows"]
    ctx.decide(
        "BETA/Z direction for missing-SE recovery",
        {
            "comparable missing-SE rows": _count(comparable),
            "discordant rows": "%s (%s)" % (
                _count(discordant),
                _percent(discordant, comparable),
            ),
            "maximum discordant fraction": "{:.2%}".format(
                assessment["maximum_fraction"]
            ),
            "row action": assessment["row_action"],
        },
        (
            "continue; discordant rows are rejected and reported per chromosome"
            if discordant
            else "continue; no directionally incompatible pair was found"
        ),
    )
    if discordant:
        message = (
            "BETA/Z SIGN DISCORDANCE: {:,} of {:,} comparable missing-SE "
            "rows ({:.4%}) have incompatible directions. This does not exceed "
            "effect_from_z.max_beta_z_sign_mismatch_fraction={:g} ({:.2%}), "
            "so the run will continue. Every affected variant that reaches "
            "chromosome step 06 will be removed and reported with reject reason "
            "beta_z_sign_discordant; the configured rejected-variant output "
            "records it when that audit is enabled."
            .format(
                discordant,
                comparable,
                assessment["discordant_fraction"],
                assessment["maximum_fraction"],
                assessment["maximum_fraction"],
            )
        )
        decisions["beta_z_sign_warning"] = message
        emit_high_visibility_warning(logger, message)


def _decide_pvalue_type(
    df, sample_column_dict, pol, ctx, decisions, statistics, logger=None,
):
    configured = pol.get("pvalue.type")
    declared = optional_text(sample_column_dict.get("declared_pvalue_type"))
    declared = str(declared).lower() if declared is not None else None
    declared = declared if declared in ("raw", "neglog10") else None
    column = _column(sample_column_dict, PVALUE_COLUMN_KEYS, df)
    threshold = pol.get("pvalue.mlogp_detect_threshold")

    decisions["pvalue_column"] = column

    from .p_values import (
        ExcessiveNegativePValueError,
        assess_negative_pvalues,
        detect_p_value_type,
    )

    detector_statistics = statistics.get("pvalue_detector") or {}

    def validate_negative_fraction():
        if column is None:
            return {}
        try:
            result = assess_negative_pvalues(detector_statistics, column, pol)
        except ExcessiveNegativePValueError as exc:
            raise StudyPropertyError(str(exc)) from exc
        if result["n_negative"]:
            value_word = "value" if result["n_negative"] == 1 else "values"
            action = str(pol.get("pvalue.out_of_range"))
            handling = (
                "Any such rows that reach chromosome p-value harmonisation will "
                "have their p-values set to missing under pvalue.out_of_range=null; "
                "later completeness checks still apply."
                if action == "null" else
                "Any such rows that reach chromosome p-value harmonisation will "
                "stop that chromosome under pvalue.out_of_range=fail."
                if action == "fail" else
                "Any such rows that reach chromosome p-value harmonisation will "
                "be rejected with reason pval_out_of_range."
            )
            message = (
                "NEGATIVE P-VALUES: column %r contains %s negative usable %s "
                "out of %s (%s). This does not exceed "
                "pvalue.max_negative_fraction=%.6g (%s), so harmonisation will "
                "continue. %s They will never be converted to p=1."
                % (
                    column,
                    _count(result["n_negative"]),
                    value_word,
                    _count(detector_statistics.get("n_usable") or 0),
                    "%.4f%%" % (result["negative_fraction"] * 100.0),
                    result["max_negative_fraction"],
                    "%.2f%%" % (result["max_negative_fraction"] * 100.0),
                    handling,
                )
            )
            decisions["pvalue_negative_warning"] = message
            emit_high_visibility_warning(logger, message)
        return result

    if configured in ("raw", "neglog10") and declared is None:
        decisions["pvalue_type"] = configured
        decisions["pvalue_type_source"] = "policy"
        decisions["pvalue_evidence"] = (
            _pvalue_evidence(
                df, column, threshold, statistics.get("pvalue_evidence"),
            ) if column is not None else {}
        )
        decisions["pvalue_evidence"].update(validate_negative_fraction())
        ctx.decide(
            "P-value scale",
            "policy pvalue.type is set explicitly, so no detection was run",
            "%s. Applied to all chromosomes." % _pvalue_words(configured),
        )
        return

    if column is None:
        decisions["pvalue_type"] = None
        decisions["pvalue_type_source"] = "no_pvalue_column"
        decisions["pvalue_evidence"] = {}
        ctx.skip(
            "No p-value column is configured, so the p-value scale was not decided "
            "here. It is derived from the effect size and standard error instead."
        )
        return

    evidence = _pvalue_evidence(
        df, column, threshold, statistics.get("pvalue_evidence"),
    )
    decisions["pvalue_evidence"] = evidence

    if evidence["non_null"] == 0:
        raise StudyPropertyError(
            "The p-value column '%s' has no usable numeric value in any of its %s rows, "
            "so its scale cannot be determined. Check the column name in the config, or "
            "set pvalue.type to 'raw' or 'neglog10' explicitly."
            % (column, _count(evidence["rows"]))
        )

    evidence.update(validate_negative_fraction())

    detector_policies = (
        pol.with_overrides({"pvalue.type": "auto"}) if declared is not None else pol
    )
    try:
        verdict, detector_evidence = detect_p_value_type(
            df, column, detector_policies,
            statistics=detector_statistics,
        )
    except ValueError as exc:
        if declared is not None:
            decisions["pvalue_type"] = declared
            decisions["pvalue_type_source"] = "sample_sheet"
            decisions["pvalue_type_detected"] = None
            decisions["pvalue_type_matches_declaration"] = None
            ctx.warn(
                "The p_value_type declaration %r remains authoritative, but automatic "
                "detection could not independently classify column %r: %s"
                % (declared, column, exc)
            )
            return
        raise StudyPropertyError(
            "The p-value detector refused to classify column '%s': %s" % (column, exc)
        ) from exc
    if detector_evidence:
        evidence["detector_stats"] = _detector_stats(detector_evidence)

    decisions["pvalue_type_detected"] = verdict
    if declared is not None:
        decisions["pvalue_type"] = declared
        decisions["pvalue_type_source"] = "sample_sheet"
        decisions["pvalue_type_matches_declaration"] = verdict == declared
        if verdict != declared:
            action = str(pol.get("validation.declaration_mismatch_action"))
            message = (
                "DECLARATION MISMATCH: sample sheet declares p_value_type=%r, but "
                "automatic detection found %r for column %r. Policy "
                "validation.declaration_mismatch_action=%r; the declaration remains "
                "authoritative."
                % (declared, verdict, column, action)
            )
            decisions.setdefault("declaration_mismatches", []).append(message)
            emit_high_visibility_warning(logger, message)
            if action == "fail":
                raise StudyPropertyError(message)
    else:
        decisions["pvalue_type"] = verdict
        decisions["pvalue_type_source"] = "detector"
        warning = _automatic_pvalue_inference_warning(verdict, detector_evidence)
        decisions["pvalue_type_inference_warning"] = warning
        emit_high_visibility_warning(
            logger, warning,
            screen_message=_automatic_pvalue_inference_warning(
                verdict, detector_evidence, screen=True,
            ),
        )

    ctx.decide(
        "P-value column '%s'" % column,
        {
            "non-null values": _count(evidence["non_null"]),
            "above %s" % threshold: "%s (%s)" % (
                _count(evidence["above_threshold"]),
                _percent(evidence["above_threshold"], evidence["non_null"]),
            ),
            "median": evidence["median"],
            "maximum": evidence["max"],
            "decided by": decisions["pvalue_type_source"],
        },
        "%s. Applied to all chromosomes."
        % _pvalue_words(decisions["pvalue_type"]),
    )


def _pvalue_words(kind):
    # type: (Optional[str]) -> str
    return {
        "raw": "raw p-values",
        "neglog10": "-log10 p-values",
    }.get(kind, str(kind))


def _automatic_pvalue_inference_warning(kind, evidence, *, screen=False):
    """Explain the applied scale using detector evidence, without another scan."""
    usable = int(evidence["n_usable"])
    above = int(evidence["n_above_threshold"])
    threshold = format_number(evidence["detect_threshold"])
    median = format_number(evidence["median"])
    if kind == "neglog10":
        low, high = evidence["median_acceptance_range"]
        median += "; -log10 rule: %s to %s inclusive" % (
            format_number(low), format_number(high),
        )
        processing = (
            "Convert -log10(P) back to raw P using P = 10^(-input value); "
            "chromosome QC still applies."
        )
        explanation = "The count, percentage and median requirements for -log10(P) were met"
    else:
        median += "; not used to select raw scale"
        processing = "Keep values on the raw probability scale; chromosome QC still applies."
        explanation = "The count and percentage requirements for -log10(P) were not both met"
    fields = [
        ("Usable values", _count(usable) + " finite numeric values"),
        ("Excluded values", _count(int(evidence["n_rows"]) - usable)
         + " missing, non-numeric or non-finite values; excluded from inference, not removed here"),
        ("Observed range", "%s to %s" % (
            format_number(evidence["min"]), format_number(evidence["max"]),
        )),
        ("Above threshold", "%s / %s (%s) exceed %s" % (
            _count(above), _count(usable), _percent(above, usable, places=4), threshold,
        )),
        ("-log10 rule", "At least %s values AND at least %s%% must exceed %s" % (
            _count(evidence["detect_min_count"]),
            format_number(float(evidence["detect_proportion"]) * 100.0), threshold,
        )),
        ("Median value", median),
        ("Why inferred", explanation),
    ]
    if evidence["n_negative"]:
        fields.append(("Negative values", "%s / %s (%s); maximum %s%%; handled separately by chromosome QC" % (
            _count(evidence["n_negative"]), _count(usable),
            _percent(evidence["n_negative"], usable, places=4),
            format_number(float(evidence["max_negative_fraction"]) * 100.0),
        )))
    reason = (
        "A value-based rule cannot confirm the column's meaning or that every "
        "value is valid. It does not determine whether tests are one- or two-sided."
    )
    recommendation = (
        "Check the study documentation. If confirmed, set p_value_type = %s "
        "in the sample sheet, or pvalue.type in YAML." % kind
    )
    if screen:
        return format_inference_warning(
            "P-value type automatically inferred",
            input_column=evidence["pval_col"],
            inferred_type="Raw p-values" if kind == "raw" else "-log10(P)",
            evidence=fields, planned_action=processing,
            reason=reason, recommendation=recommendation,
        )
    return ("AUTOMATIC P-VALUE-TYPE INFERENCE: column %r was inferred as %r. "
            "%s %s %s %s Continuing with the inferred type.") % (
        evidence["pval_col"], kind,
        "; ".join("%s: %s" % field for field in fields) + ".",
        processing, reason, recommendation,
    )


def _decide_or_se_scale(df, sample_column_dict, pol, ctx, decisions, logger=None):
    """Resolve one SE-scale decision for an odds-ratio study."""
    effect_type = decisions.get("effect_type")
    configured = str(pol.get("effect.se_scale"))
    effect_col = _column(sample_column_dict, EFFECT_COLUMN_KEYS, df)
    se_col = _column(sample_column_dict, ("se_col",), df)
    z_col = _column(sample_column_dict, ("imp_z_col", "z_col"), df)
    pvalue_col = _column(sample_column_dict, PVALUE_COLUMN_KEYS, df)
    decisions["se_column"] = se_col

    if effect_type != "odds_ratio":
        decisions.update({
            "se_scale": None,
            "se_scale_source": "not_applicable_to_beta",
            "se_scale_detected": None,
            "se_scale_evidence": {},
        })
        ctx.decide(
            "Standard-error scale",
            "the study effect is not an odds ratio",
            "not applicable; beta-scale SE is left unchanged",
        )
        return

    if se_col is None:
        decisions.update({
            "se_scale": None,
            "se_scale_source": "no_standard_error_column",
            "se_scale_detected": None,
            "se_scale_evidence": {},
        })
        ctx.skip(
            "The study effect is an odds ratio, but no standard-error column is "
            "supplied. There is no SE scale to determine; a later chromosome step "
            "may derive SE after OR has been converted to log-odds beta."
        )
        return

    from .effect_type import EffectTypeError, detect_or_standard_error_scale

    try:
        detected, evidence = detect_or_standard_error_scale(
            df,
            effect_col,
            se_col,
            z_col=z_col,
            pvalue_col=pvalue_col,
            pvalue_type=decisions.get("pvalue_type"),
            policies=pol,
        )
    except EffectTypeError as exc:
        if configured == "auto":
            raise StudyPropertyError(str(exc)) from exc
        detected = None
        evidence = {"status": "cross_check_failed", "error": str(exc)}
        ctx.warn(
            "Configured effect.se_scale=%r remains authoritative, but its dataset-level "
            "Z/p-value cross-check could not run: %s" % (configured, exc)
        )

    decisions["se_scale_detected"] = detected
    decisions["se_scale_evidence"] = evidence
    if configured == "auto":
        if detected is None:
            source = evidence.get("comparison_source") or "no usable Z or p-value comparison"
            raise StudyPropertyError(
                "The study effect is odds_ratio and a standard-error column is supplied, "
                "but PostGWAS could not determine whether SE is on the log-odds or raw "
                "odds-ratio scale. The %s comparison provided %s informative variants; "
                "%s are required, and one interpretation must reach "
                "effect.se_scale_min_agreement_fraction with the configured agreement margin. "
                "Set effect.se_scale explicitly to 'log_odds' (usual GWAS output) or "
                "'as_given' (raw OR-scale SE), or provide consistent Z scores or "
                "non-zero p-values. No chromosome processing was started."
                % (
                    source,
                    _count(evidence.get("informative_variants", 0)),
                    _count(evidence.get("minimum_informative_variants", 0)),
                )
            )
        selected = detected
        decision_source = "dataset_z_pvalue_cross_check"
        matches = None
    else:
        selected = configured
        decision_source = "policy"
        matches = None if detected is None else detected == configured
        if matches is False:
            action = str(pol.get("validation.declaration_mismatch_action"))
            message = (
                "DECLARATION MISMATCH: effect.se_scale=%r, but the dataset-level "
                "%s comparison supports %r. Policy "
                "validation.declaration_mismatch_action=%r; the explicit configuration "
                "remains authoritative."
                % (
                    configured,
                    evidence.get("comparison_source"),
                    detected,
                    action,
                )
            )
            decisions.setdefault("declaration_mismatches", []).append(message)
            emit_high_visibility_warning(logger, message)
            if action == "fail":
                raise StudyPropertyError(message)
        elif detected is None:
            ctx.warn(
                "effect.se_scale=%r is explicit and remains authoritative, but the "
                "dataset did not contain enough discriminating Z/p-value evidence to "
                "validate that declaration." % configured
            )

    decisions["se_scale"] = selected
    decisions["se_scale_source"] = decision_source
    decisions["se_scale_matches_declaration"] = matches
    ctx.decide(
        "Standard-error column '%s'" % se_col,
        {
            "comparison": evidence.get("comparison_source") or "not available",
            "informative variants": _count(evidence.get("informative_variants", 0)),
            "log-odds agreement": _percent(
                evidence.get("log_odds_matches", 0),
                evidence.get("informative_variants", 0),
            ),
            "raw-OR agreement": _percent(
                evidence.get("raw_odds_ratio_matches", 0),
                evidence.get("informative_variants", 0),
            ),
            "decided by": decision_source,
        },
        "%s; applied to all chromosomes before allele orientation" % selected,
    )


def _decide_eaf_is_maf(df, sample_column_dict, pol, ctx, decisions, statistics):
    # type: (pl.DataFrame, Dict[str, Any], Any, Any, Dict[str, Any]) -> None
    column = _column(sample_column_dict, EAF_COLUMN_KEYS, df)
    cutoff = pol.get("eaf.maf_decision_cutoff")

    decisions["eaf_column"] = column

    if column is None:
        decisions["eaf_is_maf"] = None
        decisions["eaf_is_maf_source"] = "no_study_frequency_column"
        decisions["eaf_evidence"] = {}
        ctx.skip(
            "The study provides no frequency column, so there is nothing to test for "
            "being a minor allele frequency. Frequencies come from the reference panel "
            "per chromosome."
        )
        return

    evidence = _eaf_evidence(df, column, statistics.get("eaf_evidence"))
    decisions["eaf_evidence"] = evidence

    if evidence["usable"] == 0:
        raise StudyPropertyError(
            "The frequency column '%s' has no usable finite numeric frequency between "
            "0 and 1 in any of its %s rows, so PostGWAS cannot decide whether it is MAF "
            "or EAF. Missing, unparseable, NaN, infinite and out-of-range values do not "
            "provide frequency evidence. Check the configured column or use a valid "
            "external EAF source."
            % (column, _count(evidence["rows"]))
        )

    low_fraction = evidence["low_fraction_of_usable"]
    suspected = low_fraction > cutoff
    decisions["eaf_is_maf"] = bool(suspected)
    decisions["eaf_is_maf_source"] = "study_level_statistic"
    decisions["eaf_maf_decision_cutoff"] = cutoff

    ctx.decide(
        "Frequency column '%s'" % column,
        {
            "usable frequencies": _count(evidence["usable"]),
            "excluded from MAF screen": _count(evidence["unusable"]),
            "at or below 0.5": "%s (%s)" % (
                _count(evidence["at_or_below_0_5"]),
                _percent(evidence["at_or_below_0_5"], evidence["usable"]),
            ),
            "decision cutoff": cutoff,
            "median": evidence["median"],
        },
        ("MAF-like and requiring reference confirmation - the "
         "per-chromosome cross-check against the reference panel decides what to do "
         "about it. Applied to all chromosomes."
         if suspected else
         "an effect allele frequency, not a minor allele frequency. Applied to all "
         "chromosomes."),
    )


# =============================================================================
# The step
# =============================================================================


def resolve_study_properties(
    df,
    sample_column_dict,
    policies=None,
    logger=None,
    step_number=7,
    step_total=8,
):
    # type: (pl.DataFrame, Dict[str, Any], Any, Any, int, int) -> Dict[str, Any]
    """Dataset step 07 - decide the study-wide properties once, before the split.

    Returns a decisions dictionary for ``main.py`` to hand to every chromosome
    worker::

        {
          "effect_type": "odds_ratio" | "beta" | None,
          "effect_type_source": "policy" | "detector" | "no_effect_column",
          "effect_evidence": {...},          # every denominator non-null
          "beta_z_sign_consistency": {...},  # complete-study missing-SE guard
          "se_scale": "log_odds" | "as_given" | None,
          "se_scale_source": ...,
          "se_scale_evidence": {...},
          "pvalue_type": "raw" | "neglog10" | None,
          "pvalue_type_source": ...,
          "pvalue_evidence": {...},
          "eaf_is_maf": True | False | None,
          "eaf_is_maf_source": ...,
          "eaf_evidence": {...},
          "verify_per_chromosome": {"effect": bool, "pvalue": bool},
          "rows": int,
        }

    The frame is never modified.  Raises :class:`StudyPropertyError` when a
    property cannot be decided with confidence, naming the policy to set.
    """
    pol = resolve_policies(policies)
    log = logger if logger is not None else NullLogger()

    decisions = {
        "step": _DATASET_STEP_LABEL,
        "rows": df.height,
        "columns": df.columns,
        "verify_per_chromosome": {
            "effect": bool(pol.get("effect.verify_per_chromosome")),
            "pvalue": bool(pol.get("pvalue.verify_per_chromosome")),
        },
    }  # type: Dict[str, Any]

    with log.step(step_number, step_total, "Study-wide properties",
                  "study_properties.resolve_study_properties",
                  rows_in=df.height, policy_keys=_POLICY_KEYS) as ctx:

        verify = decisions["verify_per_chromosome"]
        if not any(verify.values()):
            ctx.info(
                "These decisions are made once, on the whole file, and applied unchanged "
                "to every chromosome. Per-chromosome re-verification is off, so every "
                "statistic below is computed over the non-null values only - a "
                "denominator that included the nulls is what makes a mostly-empty "
                "column get misread."
            )
        else:
            ctx.info(
                "These decisions are made once, on the whole file, and applied unchanged "
                "to every chromosome. Re-verification is on for: %s; it can only warn, "
                "never change the answer."
                % ", ".join(sorted(name for name, on in verify.items() if on))
            )

        if df.height == 0:
            raise StudyPropertyError(
                "The input has no rows, so no study-wide property can be decided. "
                "Check the input file and the column configuration."
            )

        statistics = _collect_study_statistics(df, sample_column_dict, pol)

        _decide_effect_type(
            df, sample_column_dict, pol, ctx, decisions, statistics, logger=logger,
        )
        _decide_beta_z_sign_consistency(
            df, sample_column_dict, pol, ctx, decisions, statistics, logger=logger,
        )
        _decide_pvalue_type(
            df, sample_column_dict, pol, ctx, decisions, statistics, logger=logger,
        )
        _decide_or_se_scale(
            df, sample_column_dict, pol, ctx, decisions, logger=logger,
        )
        _decide_eaf_is_maf(
            df, sample_column_dict, pol, ctx, decisions, statistics,
        )

        ctx.set_rows(df.height, removed=0)
        ctx.extra["study_properties"] = {
            "effect_type": decisions.get("effect_type"),
            "beta_z_sign_consistency": decisions.get(
                "beta_z_sign_consistency"
            ),
            "se_scale": decisions.get("se_scale"),
            "pvalue_type": decisions.get("pvalue_type"),
            "eaf_is_maf": decisions.get("eaf_is_maf"),
        }

    return decisions


def finalise_eaf_decision_from_chromosomes(
    study_decisions,
    chromosome_summaries,
    completed_chromosomes=None,
):
    # type: (Dict[str, Any], Dict[str, Any], Optional[list]) -> Dict[str, Any]
    """Resolve the final dataset frequency semantics from chromosome evidence.

    The study-level distribution is only an initial MAF suspicion. When it is
    MAF-like, each completed chromosome performs the reference comparison. A
    final EAF/MAF label is assigned only when every completed chromosome has the
    same conclusive result; otherwise the final type remains unresolved.
    """
    result = dict(study_decisions or {})
    has_initial_decision = (
        "eaf_is_maf_initial" in result or "eaf_is_maf" in result
    )
    initial = result.get("eaf_is_maf_initial", result.get("eaf_is_maf"))
    initial_source = result.get(
        "eaf_is_maf_initial_source", result.get("eaf_is_maf_source")
    )
    result["eaf_is_maf_initial"] = initial
    result["eaf_is_maf_initial_source"] = initial_source

    summaries = chromosome_summaries if isinstance(chromosome_summaries, dict) else {}
    if completed_chromosomes is None:
        completed = [
            chromosome
            for chromosome, summary in summaries.items()
            if not isinstance(summary, dict) or summary.get("status") in (None, "ok")
        ]
    else:
        completed = list(completed_chromosomes)

    counts = {"eaf": 0, "maf": 0, "inconclusive": 0, "missing": 0}
    provenance_counts = {
        "study_supplied": 0,
        "reference_imputed": 0,
        "missing": 0,
    }
    for chromosome in completed:
        summary = summaries.get(chromosome, summaries.get(str(chromosome), {}))
        stage_qc = summary.get("stage_qc", {}) if isinstance(summary, dict) else {}
        eaf_qc = stage_qc.get("eaf_qc", {}) if isinstance(stage_qc, dict) else {}
        provenance = (
            eaf_qc.get("eaf_provenance") if isinstance(eaf_qc, dict) else None
        )
        if provenance in ("study_supplied", "reference_imputed"):
            provenance_counts[provenance] += 1
        else:
            provenance_counts["missing"] += 1
    if initial is True:
        for chromosome in completed:
            summary = summaries.get(chromosome, summaries.get(str(chromosome), {}))
            stage_qc = summary.get("stage_qc", {}) if isinstance(summary, dict) else {}
            eaf_qc = stage_qc.get("eaf_qc", {}) if isinstance(stage_qc, dict) else {}
            decision = (
                eaf_qc.get("maf_reference_decision")
                if isinstance(eaf_qc, dict)
                else None
            )
            if decision in ("eaf", "maf", "inconclusive"):
                counts[decision] += 1
            else:
                counts["missing"] += 1

    completed_count = len(completed)
    if not has_initial_decision:
        final = None
        source = "frequency_decision_missing"
    elif initial is False:
        final = False
        source = initial_source or "study_level_statistic"
    elif initial is None:
        # No internal study-frequency column means the validated external
        # source is effect-allele frequency by contract.
        final = False
        source = "external_effect_allele_frequency"
    elif completed_count and counts["eaf"] == completed_count:
        final = False
        source = "chromosome_reference_consensus"
    elif completed_count and counts["maf"] == completed_count:
        final = True
        source = "chromosome_reference_consensus"
    else:
        final = None
        source = "chromosome_reference_unresolved"

    result["eaf_is_maf"] = final
    result["eaf_is_maf_source"] = source
    result["eaf_reference_decisions"] = counts
    observed_provenance = [
        name
        for name in ("study_supplied", "reference_imputed")
        if provenance_counts[name]
    ]
    result["eaf_provenance"] = (
        observed_provenance[0]
        if len(observed_provenance) == 1
        and provenance_counts["missing"] == 0
        else "mixed"
        if len(observed_provenance) > 1
        else "unavailable"
    )
    result["eaf_provenance_counts"] = provenance_counts
    result["frequency_type"] = (
        "minor_allele_frequency"
        if final is True
        else "effect_allele_frequency"
        if final is False
        else "unresolved"
    )
    return result
