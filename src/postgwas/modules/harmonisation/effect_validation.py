"""
effect_validation.py — the two gates that decide whether a variant is fit to export.

Plan v3, sections 5.5 (per-chromosome step 10) and 5.7 (per-chromosome step 13).

Both functions follow the calling convention every other per-chromosome step in
``main.py`` uses::

    df, qc_info, sample_column_dict = validate_effect_statistics(
        chromosome=chromosome, df=df, sample_column_dict=sample_column_dict,
        policies=policies, logger=logger, rejects=rejects,
    )

``policies``, ``logger`` and ``rejects`` are keyword arguments with a default of
None so that an existing call site keeps working; with all three left out the
module still runs, using the registry defaults and a silent logger, and removes
variants with a plain filter instead of the collector.  That path is for unit
tests only - in the pipeline all three must be passed, because a removal that
does not go through the collector never reaches the reject file.

--------------------------------------------------------------------------
Step 10 - why the order of the SE tests is the whole point
--------------------------------------------------------------------------
Finiteness is tested BEFORE positivity.  ``inf > 0`` is True in Python and in
polars, so a ``> 0`` test alone lets an infinite SE straight through.  Infinite
SEs are not hypothetical: the default EAF policy is ``clip``, which puts an
out-of-range frequency at exactly 0 or 1, and the denominator
``2*p*(1-p)*(Neff + z**2)`` is then exactly zero, which polars evaluates to
``inf`` rather than null.  Testing ``is_finite()`` first is what makes keeping
``eaf.out_of_range = clip`` a safe default.

The second rule has no other line of defence at all. Clipping an invalid p-value
and then deriving SE can produce a positive, finite value that is invisible to
the other tests. The standard-error step therefore consolidates upstream high
and low clipping provenance into ``__se_from_clipped_pval`` only for rows whose
SE was actually derived from that p-value; this step consumes that tag.

--------------------------------------------------------------------------
Step 13 - the final completeness gate
--------------------------------------------------------------------------
Immediately before export, any variant with a missing value in a required field
is removed, where "missing" means null, NaN OR inf - a NaN is as unusable
downstream as a null.  Fields are evaluated ONE AT A TIME, in the configured
order, so a variant missing three fields is attributed to the first one, which
is what keeps "a variant appears exactly once" true.  INFO is deliberately not
in the list; it stays under ``info.on_missing``.

Python 3.8 compatible.  Imports polars and the standard library at module level;
policies / numpy / scipy are imported lazily so the module stays importable in a
partial tree.
"""

from functools import partial
from typing import Any, Dict, List, Optional, Sequence, Tuple

import polars as pl

from postgwas.core.dataframes import count_matching_rows, numeric_column
from postgwas.core.values import format_percentage

from .shared.runtime import NullLogger, configured_column, resolve_policies
from .shared.statistics import two_sided_negative_log10_p_from_z

__all__ = [
    "validate_effect_statistics",
    "final_completeness_check",
    "EffectValidationError",
    "STEP_LABEL",
    "FINAL_STEP_LABEL",
    "FIELD_COLUMN_KEYS",
    "FIELD_REJECT_REASONS",
    "CLIPPED_PVAL_TAGS",
]


#: The label written into the ``reject_step`` column by step 10.
STEP_LABEL = "10 effect_statistics_validation"

#: The label written into the ``reject_step`` column by step 13.
FINAL_STEP_LABEL = "13 final_completeness_check"

#: Boolean working column set only when SE was computed from a p-value that had
#: been clipped. A supplied or usable Z-derived SE must not be rejected merely
#: because the separately reported p-value needed range correction.
CLIPPED_PVAL_TAGS = ("__se_from_clipped_pval",)

#: policies field vocabulary  ->  the sample_column_dict key holding the real
#: column name.  Resolved at run time; never hardcode a column name.
FIELD_COLUMN_KEYS = {
    "chr": ("chr_col",),
    "pos": ("pos_col",),
    "snp": ("snp_id_col",),
    "ea": ("ea_col",),
    "oa": ("oa_col",),
    "eaf": ("eaf_col",),
    "beta": ("beta_col", "beta_or_col"),
    "se": ("se_col",),
    "zscore": ("imp_z_col", "z_col"),
    "pval": ("pval_col",),
    "info": ("imp_info_col",),
    "n": ("neff_col", "ncontrol_col", "ncase_col"),
}

#: policies field vocabulary  ->  the registered reject reason for step 13.
#: Every field accepted by ``final_check.require`` has a reason, so configured
#: completeness checks can never be silently skipped for lack of provenance.
FIELD_REJECT_REASONS = {
    "chr": "final_missing_chr",
    "pos": "final_missing_pos",
    "snp": "final_missing_snp",
    "ea": "final_missing_ea",
    "oa": "final_missing_oa",
    "eaf": "final_missing_eaf",
    "beta": "final_missing_beta",
    "se": "final_missing_se",
    "zscore": "final_missing_zscore",
    "pval": "final_missing_pval",
    "info": "final_missing_info",
    "n": "final_missing_n",
}

#: Fields whose values are numeric, so a text column spelling out "nan"/"inf"
#: is unusable.  chr / snp / ea / oa are excluded on purpose: "inf" is not a
#: plausible allele, and casting them to a number would flag every row.
_NUMERIC_FIELDS = frozenset(("pos", "eaf", "beta", "se", "zscore", "pval", "info", "n"))

_UNUSABLE_TEXT = ("nan", "-nan", "+nan", "inf", "-inf", "+inf", "infinity",
                  "-infinity", "+infinity")

_FLOAT_DTYPES = (pl.Float32, pl.Float64)

class EffectValidationError(Exception):
    """A validation rule set to ``fail`` fired, or too much of the chromosome went.

    A typed exception rather than ``sys.exit()`` so the dataset-level retry can
    catch it (plan Part 7); ``logger.step()`` logs it as FAILED and re-raises.
    """


# =============================================================================
# Small helpers
# =============================================================================


_column = configured_column
_numeric = numeric_column
_count_mask = partial(count_matching_rows, null_is_match=True)
_percent = format_percentage


def _apply_rule(df, mask, reason, action, ctx, rejects, step_label,
                check_name, plain_english, detail=None, warn_on_keep=True):
    # type: (pl.DataFrame, pl.Expr, str, str, Any, Any, str, str, str, Optional[str], bool) -> Tuple[pl.DataFrame, int, int]
    """Run one rule under its policy and return ``(df, rows removed, rows matched)``.

    ``reject`` removes the rows through the collector (which logs the counts
    itself), ``keep`` counts and reports them without removing anything, and
    ``fail`` counts them and raises if there were any.

    The two counts are returned separately and deliberately.  Under ``keep`` and
    ``fail`` some rows *match* the rule but none are *removed*, so a caller that
    conflates the two reports removals that never happened.
    """
    before = df.height

    if action == "reject":
        if rejects is not None:
            df = rejects.reject(df, mask, step_label, reason, detail=detail)
            removed = before - df.height
            return df, removed, removed
        # No collector (unit tests): filter, and log the counts ourselves so
        # rule 5 - never a QC action without before and after - still holds.
        df = df.filter(~mask.fill_null(True))
        removed = before - df.height
        ctx.qc(check_name, plain_english, before, df.height,
               reason=reason, step=step_label, warn=removed > 0)
        return df, removed, removed

    matched = _count_mask(df, mask)

    if action == "keep":
        ctx.qc(
            check_name,
            "%s Policy for this check is 'keep', so they were left in the data." % plain_english,
            before, before, reason=reason, matched=matched, outcome="kept",
            warn=matched > 0 and warn_on_keep,
            step=step_label,
        )
        return df, 0, matched

    if action == "fail":
        ctx.qc(check_name, plain_english, before, before,
               reason=reason, matched=matched, outcome="fail",
               warn=matched > 0, step=step_label)
        if matched:
            raise EffectValidationError(
                "%s %s variants (%s) match, and the policy for this check is 'fail'. "
                "Set the policy to 'reject' to remove them and record them in the "
                "reject file, or to 'keep' to export them unchanged."
                % (plain_english, "{:,}".format(matched), _percent(matched, before))
            )
        return df, 0, matched

    # An unknown action can only come from a policy registry change; say so
    # rather than silently doing nothing.
    ctx.warn(
        "Check '%s' has an unrecognised policy value %r, so it did nothing. "
        "Nothing was removed and nothing was counted." % (check_name, action)
    )
    return df, 0, 0


# =============================================================================
# Step 10 - effect statistics validation
# =============================================================================

_EFFECT_VALIDATION_POLICY_KEYS = [
    "validation.se_invalid",
    "validation.se_from_clipped_pval",
    "validation.beta_invalid",
    "validation.beta_zero",
    "validation.z_invalid",
    "validation.beta_se_z_concordance",
    "validation.beta_se_z_absolute_tolerance",
    "validation.beta_se_z_relative_tolerance",
    "validation.z_pval_concordance",
    "validation.z_pval_tolerance_log10",
    "validation.warn_reject_fraction",
    "validation.max_reject_fraction",
]

_BASIC_CHECK_POLICIES = {
    "se_invalid": "validation.se_invalid",
    "beta_invalid": "validation.beta_invalid",
    "beta_zero": "validation.beta_zero",
    "z_invalid": "validation.z_invalid",
}


def _reusable_basic_checks(state, df, columns, policies):
    # type: (Any, pl.DataFrame, Dict[str, Optional[str]], Any) -> set
    """Return upstream checks that are valid for this exact DataFrame.

    Object identity is deliberate: row-count equality alone cannot prove that
    values or row order were unchanged.  A subset created later in this step is
    still safe because removing rows cannot invalidate a per-row check already
    applied to the parent frame.
    """
    if not isinstance(state, dict):
        return set()
    if state.get("frame") is not df or state.get("row_count") != df.height:
        return set()
    if state.get("columns") != columns:
        return set()
    if float(state.get("se_division_floor", 0.0)) <= 0.0:
        return set()

    completed = set(state.get("completed") or ())
    recorded_policies = state.get("policies") or {}
    return {
        check
        for check, policy_key in _BASIC_CHECK_POLICIES.items()
        if check in completed
        and recorded_policies.get(policy_key) == policies.get(policy_key)
    }


def validate_effect_statistics(
    chromosome,
    df,
    sample_column_dict,
    policies=None,
    logger=None,
    rejects=None,
    step_number=10,
    step_total=16,
    upstream_validation=None,
):
    # type: (str, pl.DataFrame, Dict[str, Any], Any, Any, Any, int, int, Any) -> Tuple[pl.DataFrame, Dict[str, Any], Dict[str, Any]]
    """Pipeline step 10 - remove variants whose effect statistics are unusable.

    Order is load-bearing: SE null, then SE non-finite, then SE non-positive.
    ``inf > 0`` is True, so a positivity test alone would pass an infinite SE.

    ``upstream_validation`` may contain the internal certificate produced by
    the immediately preceding Z step. Standalone callers omit it and receive
    the complete validation gate.

    Returns ``(df, qc_info, sample_column_dict)``.
    """
    pol = resolve_policies(policies)
    log = logger if logger is not None else NullLogger()

    rows_in = df.height
    se_col = _column(sample_column_dict, FIELD_COLUMN_KEYS["se"], df)
    beta_col = _column(sample_column_dict, FIELD_COLUMN_KEYS["beta"], df)
    z_col = _column(sample_column_dict, FIELD_COLUMN_KEYS["zscore"], df)
    pval_col = _column(sample_column_dict, FIELD_COLUMN_KEYS["pval"], df)
    basic_columns = {"se": se_col, "beta": beta_col, "zscore": z_col}
    reused_checks = _reusable_basic_checks(
        upstream_validation, df, basic_columns, pol,
    )

    qc_info = {
        "chromosome": str(chromosome),
        "step": STEP_LABEL,
        "initial_variants": rows_in,
        "columns_used": {
            "se": se_col, "beta": beta_col, "zscore": z_col, "pval": pval_col,
        },
        "removed_by_reason": {},
        "checks_skipped": [],
        "checks_reused": sorted(reused_checks),
    }  # type: Dict[str, Any]

    with log.step(step_number, step_total, "Effect statistics validation",
                  "effect_validation.validate_effect_statistics",
                  rows_in=rows_in, policy_keys=_EFFECT_VALIDATION_POLICY_KEYS) as ctx:

        def record(reason, removed):
            if removed:
                qc_info["removed_by_reason"][reason] = (
                    qc_info["removed_by_reason"].get(reason, 0) + int(removed)
                )

        if reused_checks:
            ctx.info(
                "Reused the immediately preceding Z-step results for: %s. The exact "
                "same DataFrame and policy values were verified, so these basic "
                "columns were not scanned twice."
                % ", ".join(sorted(reused_checks))
            )

        # ------------------------------------------------------------------
        # 1. Standard error.  Null, then non-finite, then non-positive.
        # ------------------------------------------------------------------
        se_action = pol.get("validation.se_invalid")
        if se_col is None:
            qc_info["checks_skipped"].append("se_invalid:no_se_column")
            ctx.warn(
                "No standard-error column is present, so the standard-error checks "
                "could not run. Every variant passed this part of the gate untested."
            )
        elif "se_invalid" not in reused_checks:
            se = _numeric(df, se_col)

            df, removed, _matched = _apply_rule(
                df, se.is_null(), "se_null", se_action, ctx, rejects, STEP_LABEL,
                "standard error missing",
                "The standard error is missing, so no confidence interval or Z score "
                "can be formed for these variants.",
            )
            record("se_null", removed)

            df, removed, _matched = _apply_rule(
                df, ~se.is_finite(), "se_non_finite", se_action, ctx, rejects, STEP_LABEL,
                "standard error not finite",
                "The standard error is infinite or not a number. This is what a "
                "frequency clipped to exactly 0 or 1 produces, because the variance "
                "denominator 2*p*(1-p)*(Neff + z^2) becomes exactly zero. Checked "
                "before the 'greater than zero' test, because infinity passes that test.",
            )
            record("se_non_finite", removed)

            df, removed, _matched = _apply_rule(
                df, se <= 0, "se_non_positive", se_action, ctx, rejects, STEP_LABEL,
                "standard error not positive",
                "The standard error is zero or negative, which is not a possible "
                "value for a standard error.",
            )
            record("se_non_positive", removed)

        # ------------------------------------------------------------------
        # 2. SE fabricated from a clipped p-value.
        #    Positive and finite, so no other rule can see it.
        # ------------------------------------------------------------------
        clipped_action = pol.get("validation.se_from_clipped_pval")
        tags = [name for name in CLIPPED_PVAL_TAGS if name in df.columns]
        qc_info["clipped_pval_tags_present"] = tags
        if not tags:
            qc_info["checks_skipped"].append("se_from_clipped_pval:no_tag_column")
            ctx.skip(
                "No variant carries a 'standard error came from a clipped p-value' tag "
                "(%s), so this check had nothing to act on. If the upstream SE step "
                "stopped consolidating clipped-p-value provenance into this tag, a "
                "positive finite approximation could pass the ordinary SE checks "
                "without scrutiny." % ", ".join(CLIPPED_PVAL_TAGS)
            )
        else:
            mask = None
            for name in tags:
                dtype = df.schema.get(name)
                if dtype == pl.Boolean:
                    term = pl.col(name).fill_null(False)
                else:
                    term = pl.col(name).cast(pl.Boolean, strict=False).fill_null(False)
                mask = term if mask is None else (mask | term)
            df, removed, _matched = _apply_rule(
                df, mask, "se_from_clipped_pval", clipped_action, ctx, rejects, STEP_LABEL,
                "standard error from a clipped p-value",
                "The standard error was computed from a p-value that had been clipped "
                "at the top or bottom of its retained range. The result can look "
                "positive and finite even though the exact significance was lost, so "
                "no ordinary SE validity check can catch it.",
                detail="tagged by " + ", ".join(tags),
            )
            record("se_from_clipped_pval", removed)

        # ------------------------------------------------------------------
        # 3. BETA
        # ------------------------------------------------------------------
        beta_action = pol.get("validation.beta_invalid")
        zero_action = pol.get("validation.beta_zero")
        if beta_col is None:
            qc_info["checks_skipped"].append("beta_invalid:no_beta_column")
            ctx.warn(
                "No effect-size column is present, so the effect-size checks could not "
                "run. Every variant passed this part of the gate untested."
            )
        elif not {"beta_invalid", "beta_zero"}.issubset(reused_checks):
            beta = _numeric(df, beta_col)

            if "beta_invalid" not in reused_checks:
                df, removed, _matched = _apply_rule(
                    df, beta.is_null() | ~beta.is_finite(), "beta_invalid", beta_action,
                    ctx, rejects, STEP_LABEL,
                    "effect size missing or not finite",
                    "The effect size is missing, infinite or not a number, so the variant "
                    "carries no usable estimate.",
                )
                record("beta_invalid", removed)

            if "beta_zero" not in reused_checks:
                df, removed, _matched = _apply_rule(
                    df, beta == 0, "beta_zero", zero_action, ctx, rejects, STEP_LABEL,
                    "effect size exactly zero",
                    "The effect size is exactly zero. With a valid standard error this is a "
                    "legitimate null result giving a Z score of zero, which is why the "
                    "default is to keep it.",
                    warn_on_keep=False,
                )
                record("beta_zero", removed)

        # ------------------------------------------------------------------
        # 4. Z score
        # ------------------------------------------------------------------
        z_action = pol.get("validation.z_invalid")
        if z_col is None:
            qc_info["checks_skipped"].append("z_invalid:no_z_column")
            ctx.warn(
                "No Z-score column is present, so the Z-score checks could not run. "
                "Every variant passed this part of the gate untested."
            )
        elif "z_invalid" not in reused_checks:
            z = _numeric(df, z_col)
            df, removed, _matched = _apply_rule(
                df, z.is_null() | ~z.is_finite(), "z_invalid", z_action,
                ctx, rejects, STEP_LABEL,
                "Z score missing or not finite",
                "The Z score is missing, infinite or not a number, so the variant "
                "cannot be meta-analysed or fine-mapped.",
            )
            record("z_invalid", removed)

        # ------------------------------------------------------------------
        # 5. Optional supplied-Z vs BETA / SE concordance
        # ------------------------------------------------------------------
        df, concordance = _beta_se_z_concordance(
            df, beta_col, se_col, z_col, pol, ctx, rejects, qc_info,
        )
        qc_info["beta_se_z_concordance"] = concordance
        record("beta_se_z_discordant", concordance.get("removed", 0))

        # ------------------------------------------------------------------
        # 6. Optional Z vs p-value concordance
        # ------------------------------------------------------------------
        df, concordance = _z_pval_concordance(
            df, z_col, pval_col, pol, ctx, rejects, qc_info,
        )
        qc_info["z_pval_concordance"] = concordance
        record("z_pval_discordant", concordance.get("removed", 0))

        # ------------------------------------------------------------------
        # 7. How much of the chromosome did this step take?
        # ------------------------------------------------------------------
        rows_out = df.height
        removed_total = rows_in - rows_out
        fraction = (float(removed_total) / float(rows_in)) if rows_in else 0.0
        qc_info["final_variants"] = rows_out
        qc_info["removed_total"] = removed_total
        qc_info["removed_fraction"] = fraction

        warn_at = pol.get("validation.warn_reject_fraction")
        fail_at = pol.get("validation.max_reject_fraction")

        if rows_in and fraction > fail_at:
            raise EffectValidationError(
                "This step removed %s of the %s variants on chromosome %s (%s), which is "
                "above validation.max_reject_fraction (%s). Losing this much of a "
                "chromosome in one gate is usually a wrong effect-size or p-value scale "
                "upstream, not bad data. Check the DECIDE lines for the study-level "
                "effect type and p-value type, then re-run."
                % ("{:,}".format(removed_total), "{:,}".format(rows_in), chromosome,
                   _percent(removed_total, rows_in), fail_at)
            )
        if rows_in and fraction > warn_at:
            ctx.warn(
                "This step removed %s of %s variants (%s), which is above "
                "validation.warn_reject_fraction (%s). The reject file lists the reason "
                "for every one of them."
                % ("{:,}".format(removed_total), "{:,}".format(rows_in),
                   _percent(removed_total, rows_in), warn_at)
            )

        ctx.set_rows(rows_out)
        ctx.extra["removed_by_reason"] = dict(qc_info["removed_by_reason"])

    return df, qc_info, sample_column_dict


def _beta_se_z_concordance(df, beta_col, se_col, z_col, pol, ctx, rejects, qc_info):
    """Compare a supplied Z score with the Wald identity ``Z = BETA / SE``."""
    action = pol.get("validation.beta_se_z_concordance")
    absolute_tolerance = float(pol.get("validation.beta_se_z_absolute_tolerance"))
    relative_tolerance = float(pol.get("validation.beta_se_z_relative_tolerance"))
    result = {
        "action": action,
        "absolute_tolerance": absolute_tolerance,
        "relative_tolerance": relative_tolerance,
        "checked": 0,
        "discordant": 0,
        "removed": 0,
        "status": "off",
    }

    if action == "off":
        return df, result

    missing = [
        label for label, column in (
            ("BETA", beta_col), ("standard-error", se_col), ("Z-score", z_col),
        ) if column is None
    ]
    if missing:
        result["status"] = "skipped_missing_column"
        qc_info["checks_skipped"].append("beta_se_z_concordance:missing_column")
        ctx.warn(
            "The BETA/SE versus Z agreement check is switched on but the %s column "
            "is not present, so it could not run." % ", ".join(missing)
        )
        return df, result

    if df.height == 0:
        result["status"] = "empty_frame"
        return df, result

    beta = _numeric(df, beta_col)
    se = _numeric(df, se_col)
    z = _numeric(df, z_col)
    expected = beta / se
    comparable = (
        beta.is_not_null() & beta.is_finite()
        & se.is_not_null() & se.is_finite() & (se > 0)
        & z.is_not_null() & z.is_finite()
        & expected.is_finite()
    )
    difference = (z - expected).abs()
    allowed = pl.lit(absolute_tolerance) + pl.lit(relative_tolerance) * expected.abs()
    discordant = comparable & (difference > allowed)
    flags = df.select(
        comparable.cast(pl.Int64).sum().alias("checked"),
        discordant.cast(pl.Int64).sum().alias("discordant"),
    ).row(0, named=True)
    result["checked"] = int(flags["checked"] or 0)
    result["discordant"] = int(flags["discordant"] or 0)
    result["status"] = "ran"

    n_discordant = result["discordant"]
    plain = (
        "The supplied Z score disagrees with BETA / SE beyond the configured "
        "absolute tolerance %s plus relative tolerance %s. Small differences can "
        "come from published rounding; larger differences indicate incompatible "
        "effect, standard-error or Z-score columns."
        % (absolute_tolerance, relative_tolerance)
    )
    if n_discordant == 0:
        ctx.qc(
            "BETA/SE versus Z agreement",
            "Every one of the %s comparable variants agrees within the configured "
            "absolute and relative tolerances." % "{:,}".format(result["checked"]),
            df.height, df.height, step=STEP_LABEL,
        )
        return df, result

    if action == "warn":
        ctx.qc(
            "BETA/SE versus Z agreement", plain, df.height, df.height,
            reason="beta_se_z_discordant", matched=n_discordant,
            outcome="kept", warn=True, step=STEP_LABEL,
        )
        return df, result

    if action == "fail":
        ctx.qc(
            "BETA/SE versus Z agreement", plain, df.height, df.height,
            reason="beta_se_z_discordant", matched=n_discordant,
            outcome="fail", warn=True, step=STEP_LABEL,
        )
        raise EffectValidationError(
            "%s variants (%s) fail the BETA/SE versus Z agreement check and "
            "validation.beta_se_z_concordance is 'fail'. Set it to 'warn' to "
            "report them, or to 'reject' to remove them and record them in the "
            "reject file."
            % ("{:,}".format(n_discordant), _percent(n_discordant, df.height))
        )

    if action == "reject":
        mask = df.select(discordant.alias("__beta_se_z_discordant")).to_series()
        before = df.height
        if rejects is not None:
            df = rejects.reject(
                df, mask, STEP_LABEL, "beta_se_z_discordant",
                detail=(
                    "|Z - BETA/SE| > %s + %s*|BETA/SE|"
                    % (absolute_tolerance, relative_tolerance)
                ),
            )
        else:
            df = df.filter(~mask)
            ctx.qc(
                "BETA/SE versus Z agreement", plain, before, df.height,
                reason="beta_se_z_discordant", step=STEP_LABEL, warn=True,
            )
        result["removed"] = before - df.height
        return df, result

    ctx.warn(
        "validation.beta_se_z_concordance has an unrecognised value %r, so the "
        "check did nothing." % action
    )
    return df, result


def _z_pval_concordance(df, z_col, pval_col, pol, ctx, rejects, qc_info):
    # type: (pl.DataFrame, Optional[str], Optional[str], Any, Any, Any, Dict[str, Any]) -> Tuple[pl.DataFrame, Dict[str, Any]]
    """Compare -log10(2 * P(Z > |z|)) with the harmonised -log10 p.

    One check that catches a fabricated standard error, a capped p-value and a
    misdetected effect type at once. It warns without removing rows by default.
    """
    action = pol.get("validation.z_pval_concordance")
    tolerance = pol.get("validation.z_pval_tolerance_log10")
    result = {"action": action, "tolerance_log10": tolerance,
              "checked": 0, "discordant": 0, "removed": 0, "status": "off"}

    if action == "off":
        return df, result

    if z_col is None or pval_col is None:
        result["status"] = "skipped_missing_column"
        qc_info["checks_skipped"].append("z_pval_concordance:missing_column")
        ctx.warn(
            "The Z versus p-value agreement check is switched on but the %s column is "
            "not present, so it could not run."
            % ("Z-score" if z_col is None else "p-value")
        )
        return df, result

    if df.height == 0:
        result["status"] = "empty_frame"
        return df, result

    try:
        import numpy as np
        two_sided_negative_log10_p_from_z([0.0])
    except Exception as exc:                      # pragma: no cover - env dependent
        result["status"] = "skipped_scipy_unavailable"
        raise EffectValidationError(
            "The enabled Z versus p-value agreement check requires numpy and scipy, "
            "but they could not be imported (%s). Install the harmonisation runtime "
            "dependencies, or explicitly set validation.z_pval_concordance to 'off'."
            % exc
        )

    z_values = df.get_column(z_col).cast(pl.Float64, strict=False).to_numpy()
    p_values = df.get_column(pval_col).cast(pl.Float64, strict=False).to_numpy()

    with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
        expected = two_sided_negative_log10_p_from_z(z_values)
        observed = -np.log10(p_values)
        difference = np.abs(expected - observed)

    comparable = np.isfinite(difference)
    discordant = comparable & (difference > float(tolerance))

    result["checked"] = int(comparable.sum())
    result["discordant"] = int(discordant.sum())
    result["status"] = "ran"

    n_discordant = result["discordant"]
    plain = (
        "The p-value implied by the Z score disagrees with the reported p-value by "
        "more than %s orders of magnitude. That happens when the standard error was "
        "fabricated from a clipped p-value, when the p-value column was capped, or "
        "when the effect-size type was detected wrongly." % tolerance
    )

    if n_discordant == 0:
        ctx.qc("Z versus p-value agreement",
               "Every one of the %s comparable variants agrees to within %s orders of "
               "magnitude." % ("{:,}".format(result["checked"]), tolerance),
               df.height, df.height, step=STEP_LABEL)
        return df, result

    if action == "warn":
        ctx.qc("Z versus p-value agreement", plain, df.height, df.height,
               reason="z_pval_discordant", matched=n_discordant,
               outcome="kept", warn=True, step=STEP_LABEL)
        return df, result

    if action == "fail":
        ctx.qc("Z versus p-value agreement", plain, df.height, df.height,
               reason="z_pval_discordant", matched=n_discordant,
               outcome="fail", warn=True, step=STEP_LABEL)
        raise EffectValidationError(
            "%s variants (%s) fail the Z versus p-value agreement check and "
            "validation.z_pval_concordance is 'fail'. Set it to 'warn' to report them, "
            "or to 'reject' to remove them and record them in the reject file."
            % ("{:,}".format(n_discordant), _percent(n_discordant, df.height))
        )

    if action == "reject":
        mask = pl.Series("__z_pval_discordant", discordant)
        before = df.height
        if rejects is not None:
            df = rejects.reject(df, mask, STEP_LABEL, "z_pval_discordant",
                                detail="|expected - observed| > %s on the -log10 scale"
                                       % tolerance)
        else:
            df = df.filter(~mask)
            ctx.qc("Z versus p-value agreement", plain, before, df.height,
                   reason="z_pval_discordant", step=STEP_LABEL, warn=True)
        result["removed"] = before - df.height
        return df, result

    ctx.warn(
        "validation.z_pval_concordance has an unrecognised value %r, so the check "
        "did nothing." % action
    )
    return df, result


# =============================================================================
# Step 13 - final completeness check
# =============================================================================

_FINAL_CHECK_POLICY_KEYS = [
    "final_check.require",
    "final_check.treat_as_missing",
    "final_check.on_missing",
]


def _missing_mask(df, column, field, treat_as_missing):
    # type: (pl.DataFrame, str, str, Sequence[str]) -> Optional[pl.Expr]
    """The 'this value is unusable' mask for one column, per treat_as_missing.

    null applies to every dtype.  NaN and inf only exist in a float column; when
    the column is still text and the field is a numeric one, the text spellings
    ('nan', 'inf', '-inf', ...) are matched instead, because a column that has
    not been cast yet would otherwise slip through the gate.
    """
    wanted = set(str(item).lower() for item in treat_as_missing)
    dtype = df.schema.get(column)
    terms = []  # type: List[pl.Expr]

    if "null" in wanted:
        terms.append(pl.col(column).is_null())

    if dtype in _FLOAT_DTYPES:
        if "nan" in wanted:
            terms.append(pl.col(column).is_nan().fill_null(False))
        if "inf" in wanted:
            terms.append(pl.col(column).is_infinite().fill_null(False))
    elif dtype == pl.Utf8 and field in _NUMERIC_FIELDS:
        literals = []
        if "nan" in wanted:
            literals.extend([t for t in _UNUSABLE_TEXT if "nan" in t])
        if "inf" in wanted:
            literals.extend([t for t in _UNUSABLE_TEXT if "inf" in t])
        if literals:
            terms.append(
                pl.col(column).str.strip_chars().str.to_lowercase()
                .is_in(literals).fill_null(False)
            )

    if not terms:
        return None
    mask = terms[0]
    for term in terms[1:]:
        mask = mask | term
    return mask


def final_completeness_check(
    chromosome,
    df,
    sample_column_dict,
    policies=None,
    logger=None,
    rejects=None,
    step_number=13,
    step_total=16,
):
    # type: (str, pl.DataFrame, Dict[str, Any], Any, Any, Any, int, int) -> Tuple[pl.DataFrame, Dict[str, Any], Dict[str, Any]]
    """Pipeline step 13 - the last gate before export.

    Removes any variant with a missing value in a required field, where missing
    means null, NaN or inf.  One field at a time, in the configured order, so a
    variant missing three fields is attributed to the first - the same rule that
    makes the reject counts reconcile.  Also drops every internal ``__`` working
    column, so nothing derived leaks into the exported file.

    Returns ``(df, qc_info, sample_column_dict)``.
    """
    pol = resolve_policies(policies)
    log = logger if logger is not None else NullLogger()

    rows_in = df.height
    required = list(pol.get("final_check.require"))
    treat_as_missing = list(pol.get("final_check.treat_as_missing"))
    action = pol.get("final_check.on_missing")

    qc_info = {
        "chromosome": str(chromosome),
        "step": FINAL_STEP_LABEL,
        "initial_variants": rows_in,
        "required_fields": required,
        "treat_as_missing": treat_as_missing,
        "on_missing": action,
        "removed_by_reason": {},
        "missing_by_field": {},
        "fields_skipped": {},
        "columns_used": {},
    }  # type: Dict[str, Any]

    with log.step(step_number, step_total, "Final completeness check",
                  "effect_validation.final_completeness_check",
                  rows_in=rows_in, policy_keys=_FINAL_CHECK_POLICY_KEYS) as ctx:

        for field in required:
            keys = FIELD_COLUMN_KEYS.get(field)
            reason = FIELD_REJECT_REASONS.get(field)

            if reason is None:
                raise EffectValidationError(
                    "Field '%s' is listed in final_check.require but has no registered "
                    "rejection reason. This is an internal configuration error; add the "
                    "field to FIELD_REJECT_REASONS and the rejection-reason registry "
                    "before using it." % field
                )

            column = _column(sample_column_dict, keys or (), df) if keys else None
            qc_info["columns_used"][field] = column

            if column is None:
                raise EffectValidationError(
                    "Field '%s' is required by final_check.require, but no mapped column "
                    "is present in the chromosome data. This is a dataset-level schema "
                    "failure: every variant would be missing the required field. Add or "
                    "correct the sample-sheet column mapping, or remove the field from "
                    "final_check.require if it is not required. final_check.on_missing "
                    "only controls missing values inside a column that exists; it cannot "
                    "permit an entirely absent required column." % field
                )

            mask = _missing_mask(df, column, field, treat_as_missing)
            if mask is None:
                qc_info["fields_skipped"][field] = "nothing_treated_as_missing"
                ctx.skip(
                    "Field '%s': final_check.treat_as_missing lists nothing that applies "
                    "to column '%s', so no variant could be judged missing."
                    % (field, column)
                )
                continue

            plain = (
                "Column '%s' holds the %s and its value is missing here - null, not a "
                "number, or infinite. Nothing downstream can use a variant without it."
                % (column, field)
            )
            before = df.height
            df, _removed_unused, matched = _apply_rule(
                df, mask, reason, action, ctx, rejects, FINAL_STEP_LABEL,
                "final check - %s missing" % field, plain,
            )
            qc_info["missing_by_field"][field] = int(matched)
            if action == "reject":
                removed = before - df.height
                if removed:
                    qc_info["removed_by_reason"][reason] = int(removed)

        # ------------------------------------------------------------------
        # Internal working columns never reach the exported file.
        # ------------------------------------------------------------------
        internal = [name for name in df.columns if name.startswith("__")]
        if internal:
            df = df.drop(internal)
            ctx.info(
                "Dropped %d internal working column%s before export: %s."
                % (len(internal), "" if len(internal) == 1 else "s", ", ".join(internal))
            )
        qc_info["internal_columns_dropped"] = internal

        rows_out = df.height
        qc_info["final_variants"] = rows_out
        qc_info["removed_total"] = rows_in - rows_out
        qc_info["removed_fraction"] = (
            (float(rows_in - rows_out) / float(rows_in)) if rows_in else 0.0
        )

        if action == "fail":
            total_missing = sum(qc_info["missing_by_field"].values())
            if total_missing:
                raise EffectValidationError(
                    "%s variants are missing at least one required field and "
                    "final_check.on_missing is 'fail'. Counts by field: %s. Set it to "
                    "'reject' to remove them and record them in the reject file."
                    % ("{:,}".format(total_missing),
                       ", ".join("%s %s" % (k, "{:,}".format(v))
                                 for k, v in qc_info["missing_by_field"].items() if v))
                )

        ctx.set_rows(rows_out)
        ctx.extra["removed_by_reason"] = dict(qc_info["removed_by_reason"])

    return df, qc_info, sample_column_dict
