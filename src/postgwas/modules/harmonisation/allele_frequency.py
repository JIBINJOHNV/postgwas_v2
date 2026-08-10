"""Step 04 (per chromosome) — reference orientation and EAF harmonisation.

Plan v3 sections 5.2 and 5.3.

What this module does, in order:

1.  Join the raw chromosome reference once on chromosome and position, resolve
    forward/swapped/reverse-complement orientations, and transform allele-specific
    statistics consistently.
2.  If the study supplies its own frequency column, describe it, force out-of-range
    values back into range (policy ``eaf.out_of_range``), and decide whether the
    column is really an *effect* allele frequency or a *minor* allele frequency.
    A MAF-like distribution is checked against the configured default EAF table:
    direct matches use ALT frequency and swapped matches use ``1-AF``. This follows
    the GWAS-SSF requirement that frequency is aligned to the listed effect allele
    (https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics).
3.  Otherwise take the frequency from an external panel, matching each variant in
    the two orientations named by ``eaf.match_orientations`` (direct and swap —
    never complement, see plan section 5.3), optionally arbitrating palindromic
    variants with the study's own frequency.
4.  Remove frequencies of exactly 0 or 1 (policy ``eaf.degenerate``); those make
    the standard-error formula divide by zero further down the pipeline.

Logging and rejected variants
-----------------------------
Every public function accepts ``logger=`` (a ``PipelineLogger``), ``ctx=`` (a
``StepContext``) and ``rejects=`` (a ``RejectCollector``) as optional keyword
arguments. Pipeline runs use one structured chromosome log and one unified
reject file.

Policies
--------
Every threshold is a policy. ``policies=None`` means "use the canonical YAML
defaults"; ``eaf.degenerate`` defaults to ``reject`` as documented there.
"""

import contextlib
import os
from typing import Any, Dict, Optional, Tuple

import polars as pl

from postgwas.core.paths import configured_output_path
from postgwas.core.io.tables import read_delimited_table

from postgwas.core.values import (
    format_count,
    format_fraction_percentage,
    optional_text,
)

from .shared.runtime import resolve_policies
from .shared.allele_join import allele_oriented_left_join
from .strand import (
    POLICY_KEYS as STRAND_POLICY_KEYS,
    REFERENCE_AF_COLUMN,
    STRAND_ACTION_COLUMN,
    harmonise_strand_orientation,
)

__all__ = [
    "AlleleFrequencyError",
    "STEP_LABEL",
    "STEP_TITLE",
    "POLICY_KEYS",
    "PALINDROMIC_FLAG_COLUMN",
    "read_reference_table",
    "summarise_allele_frequency",
    "validate_allele_frequency",
    "confirm_eaf_with_reference",
    "resolve_palindromic_frequency",
    "merge_external_allele_frequencies",
    "harmonise_allele_frequency",
]


# -------------------------------------------------------
# 0. Module constants and the typed error
# -------------------------------------------------------
STEP_LABEL = "04 eaf_harmonisation"
STEP_TITLE = "Allele orientation and effect allele frequency"
FUNC_NAME = "allele_frequency.harmonise_allele_frequency"

#: Column added to the frame when ``eaf.palindromic_handling`` is not ``ignore``.
PALINDROMIC_FLAG_COLUMN = "eaf_palindromic_flag"

#: An ordered allele pair whose complement is itself.
PALINDROMIC_PAIRS = ("AT", "TA", "CG", "GC")

#: The settings block printed at the top of the step.
POLICY_KEYS = (
    "eaf.out_of_range",
    "eaf.clip_tolerance",
    "eaf.degenerate",
    "eaf.maf_decision_cutoff",
    "eaf.maf_reference_min_overlap",
    "eaf.on_maf_check_inconclusive",
    "eaf.external_min_match_fraction",
    "eaf.deduplicate_reference",
    "eaf.match_orientations",
    "eaf.palindromic_handling",
    "eaf.palindromic_ambiguous_lower",
    "eaf.palindromic_ambiguous_upper",
    "eaf.palindromic_resolve_tolerance",
    "eaf.invalid_fraction_cutoff",
    "eaf.missing_fraction_cutoff",
)

class AlleleFrequencyError(RuntimeError):
    """A fatal problem in the EAF step.

    Derives from ``RuntimeError`` so the dataset orchestrator can catch it,
    record the chromosome failure, and apply its retry policy.
    """


# -------------------------------------------------------
class _Emitter(object):
    """Route messages to the active StepContext or PipelineLogger."""

    __slots__ = ("logger", "ctx")

    def __init__(self, logger=None, ctx=None):
        self.logger = logger
        self.ctx = ctx

    @property
    def target(self):
        if self.ctx is not None:
            return self.ctx
        return self.logger

    @property
    def wired(self) -> bool:
        return self.target is not None

    def info(self, message, indent=0):
        target = self.target
        if target is not None:
            target.info(message, indent=indent)

    def warn(self, message, indent=0):
        target = self.target
        if target is not None:
            target.warn(message, indent=indent)

    def decide(self, what, evidence, decision):
        target = self.target
        if target is not None:
            target.decide(what, evidence, decision)

    def qc(self, check_name, plain_english, before, after, reason=None, changed=None, warn=False):
        """Record a QC action.  ``before`` and ``after`` are mandatory (rule 5)."""
        if before is None or after is None:
            raise AlleleFrequencyError(
                "QC check '%s' was recorded without a variant count before and after. "
                "Every QC action must carry both." % check_name
            )
        target = self.target
        if target is not None:
            return target.qc(
                check_name,
                plain_english,
                int(before),
                int(after),
                reason=reason,
                changed=changed,
                warn=warn,
                step=STEP_LABEL,
            )
        removed = int(before) - int(after)
        return {
            "step": STEP_LABEL,
            "check": check_name,
            "reason": reason,
            "before": int(before),
            "after": int(after),
            "removed": removed,
            "changed": changed,
        }

_fmt_int = format_count
_fmt_pct = format_fraction_percentage


def _pol(policies, key, override=None):
    """Policy value for ``key``, unless the caller passed an explicit override."""
    if override is not None:
        return override
    return resolve_policies(policies).get(key)


def _study_decision_flag(study_decision) -> Optional[bool]:
    """Read 'is the study frequency column a MAF?' out of a study-level decision.

    Accepts a plain bool, or a dict carrying one of the documented keys.  Returns
    ``None`` when the decision says nothing about the frequency column.
    """
    if study_decision is None:
        return None
    if isinstance(study_decision, bool):
        return study_decision
    if isinstance(study_decision, dict):
        for key in ("eaf_is_maf", "is_maf", "eaf_is_minor_allele_frequency", "maf_like"):
            if key in study_decision and study_decision[key] is not None:
                return bool(study_decision[key])
        return None
    raise AlleleFrequencyError(
        "study_decision must be a bool or a dict such as {'eaf_is_maf': False}; got %r" % type(study_decision)
    )


# -------------------------------------------------------
# 3. Removing variants
# -------------------------------------------------------
def _reject_rows(
    df: pl.DataFrame,
    mask: pl.Expr,
    reason: str,
    emit: "_Emitter",
    rejects=None,
    detail: Optional[str] = None,
    check_name: Optional[str] = None,
    plain_english: Optional[str] = None,
) -> pl.DataFrame:
    """Remove the rows where ``mask`` is true and return the survivors.

    With a ``RejectCollector`` the rows are written to the rejected-variants file
    and the collector logs the QC line itself.  Without one the rows are still
    removed — the policy asked for it — but a warning says they could not be
    written out, and the QC line is emitted here so the before/after counts are
    never lost.
    """
    before = df.height
    if rejects is not None:
        survivors = rejects.reject(df, mask, STEP_LABEL, reason, detail=detail)
        if getattr(rejects, "logger", None) is None:
            emit.qc(
                check_name or reason,
                plain_english or "",
                before,
                survivors.height,
                reason=reason,
                warn=survivors.height < before,
            )
        return survivors
    survivors = df.filter(~mask.fill_null(True))
    removed = before - survivors.height
    emit.qc(
        check_name or reason,
        plain_english or "",
        before,
        survivors.height,
        reason=reason,
        warn=removed > 0,
    )
    if removed > 0:
        emit.warn(
            "No rejected-variants collector was supplied to this step, so those %s variants were "
            "removed without being written to the rejected-variants file." % _fmt_int(removed)
        )
    return survivors


def read_reference_table(
    path: str,
    policies,
    configured_delimiter: str,
    what: str = "reference file",
    emit: Optional["_Emitter"] = None,
) -> pl.DataFrame:
    """Read a headed text table, detecting the separator from its header line.

    Raises :class:`AlleleFrequencyError` when the file cannot be parsed.
    """
    pol = resolve_policies(policies)
    frame, detected = read_delimited_table(
        path,
        configured_delimiter,
        candidates=list(pol.get("input.delimiter_candidates")),
        minimum_columns=int(pol.get("input.delimiter_min_columns")),
        maximum_columns=int(pol.get("input.delimiter_max_columns")),
        sample_lines=int(pol.get("input.delimiter_sample_rows")),
        null_values=list(pol.get("input.null_values")),
        infer_schema_length=int(pol.get("input.schema_inference_rows")),
        error_type=AlleleFrequencyError,
        description=what,
    )
    if emit is not None:
        emit.info(
            "Read the %s %s: %s rows, %d columns, %s separated."
            % (
                what, os.path.basename(path), _fmt_int(frame.height), frame.width,
                detected.method,
            )
        )
    return frame


# -------------------------------------------------------
# 5. Helper: Compute detailed AF/EAF statistics
# -------------------------------------------------------
def summarise_allele_frequency(
    df: pl.DataFrame,
    colname: str,
    clip_tolerance: Optional[float] = None,
    policies=None,
) -> Dict[str, Any]:
    """
    Comprehensive descriptive statistics for AF/EAF-like columns.
    Works on the DataFrame exactly as passed in.

    ``clip_tolerance`` (policy ``eaf.clip_tolerance``, default 1.05) is the value
    above which a frequency is corruption rather than rounding error.
    """
    if colname not in df.columns:
        return {
            "column_present": False,
            "total_rows": df.height,
        }
    tolerance = float(_pol(policies, "eaf.clip_tolerance", clip_tolerance))
    total_rows = df.height
    if total_rows == 0:
        return {
            "column_present": True,
            "total_rows": 0,
            "non_null_rows": 0,
            "null_rows": 0,
            "clip_tolerance": tolerance,
        }
    stats_row = df.select([
        pl.len().alias("total_rows"),
        pl.col(colname).is_not_null().sum().alias("non_null_rows"),
        pl.col(colname).is_null().sum().alias("null_rows"),
        (pl.col(colname) < 0.0).sum().alias("lt_0_rows"),
        (pl.col(colname) > 1.0).sum().alias("gt_1_rows"),
        (pl.col(colname) > tolerance).sum().alias("gt_1_05_rows"),
        (pl.col(colname) < 0.0).sum().alias("invalid_negative_rows"),
        ((pl.col(colname) < 0.0) | (pl.col(colname) > tolerance)).sum().alias("invalid_total_rows"),
        (pl.col(colname) == 0.0).sum().alias("eq_0_rows"),
        (pl.col(colname) == 1.0).sum().alias("eq_1_rows"),
        (pl.col(colname) <= 0.5).sum().alias("le_0_5_rows"),
        pl.col(colname).min().alias("min"),
        pl.col(colname).max().alias("max"),
        pl.col(colname).mean().alias("mean"),
        pl.col(colname).median().alias("median"),
        pl.col(colname).quantile(0.25).alias("q25"),
        pl.col(colname).quantile(0.75).alias("q75"),
    ]).to_dicts()[0]
    non_null_rows = stats_row["non_null_rows"] or 0
    total_rows = stats_row["total_rows"] or 0

    def frac(n, d):
        return float(n) / float(d) if d else 0.0

    stats_row["null_fraction"] = frac(stats_row["null_rows"], total_rows)
    stats_row["invalid_total_fraction"] = frac(stats_row["invalid_total_rows"], total_rows)
    stats_row["le_0_5_fraction_among_total"] = frac(stats_row["le_0_5_rows"], total_rows)
    stats_row["le_0_5_fraction_among_non_null"] = frac(stats_row["le_0_5_rows"], non_null_rows)
    stats_row["degenerate_rows"] = (stats_row["eq_0_rows"] or 0) + (stats_row["eq_1_rows"] or 0)
    stats_row["clip_tolerance"] = tolerance
    stats_row["column_present"] = True
    return stats_row


# -------------------------------------------------------
# 6. Helper: Validate + clean EAF stats
# -------------------------------------------------------
def validate_allele_frequency(
    df: pl.DataFrame,
    colname: str,
    maf_cutoff: float = 0.95,
    invalid_fraction_cutoff: float = 0.05,
    missing_fraction_cutoff: float = 0.05,
    policies=None,
    logger=None,
    ctx=None,
    rejects=None,
    out_of_range: Optional[str] = None,
    clip_tolerance: Optional[float] = None,
) -> Tuple[bool, bool, int, int, pl.DataFrame, float, Dict[str, Any]]:
    """
    Validate and clean AF/EAF column.

    Out-of-range values are counted *before* anything is modified, then handled
    according to ``eaf.out_of_range``:

    ``clip``   force them back into 0..1 (the default, and what the code did
               before this migration)
    ``null``   blank them, keeping the variant
    ``reject`` remove the variant, reason ``eaf_out_of_range``
    ``fail``   raise :class:`AlleleFrequencyError`

    Returns
    -------
    (
        is_valid_range,
        is_suspicious_maf,
        out_of_range_count,
        missing_count,
        cleaned_df,
        low_freq_pct_after_cleaning,
        stats_dict
    )

    ``is_valid_range`` is now reported honestly (it used to be hardcoded
    ``True``): it is ``False`` when the share of values outside 0..tolerance
    exceeds ``invalid_fraction_cutoff``.  ``missing_count`` is the null count of
    the column — it used to be a row-count difference across ``with_columns``,
    which never changes the row count and so was always 0.
    """
    emit = _Emitter(logger=logger, ctx=ctx)
    tolerance = float(_pol(policies, "eaf.clip_tolerance", clip_tolerance))
    action = str(_pol(policies, "eaf.out_of_range", out_of_range))

    if colname not in df.columns:
        stats = {"column_present": False, "total_rows": df.height}
        return False, False, 0, df.height, df, 0.0, stats
    total_before = df.height
    if total_before == 0:
        stats = summarise_allele_frequency(df, colname, clip_tolerance=tolerance)
        return False, False, 0, 0, df, 0.0, stats

    initial_stats = summarise_allele_frequency(df, colname, clip_tolerance=tolerance)

    # ``summarise_allele_frequency`` already calculated these values in one
    # aggregate pass over the untouched column; reuse them rather than scanning
    # the chromosome a second time before applying the configured action.
    out_of_range_count = int(initial_stats["invalid_total_rows"] or 0)
    outside_unit_count = int(
        (initial_stats["lt_0_rows"] or 0) + (initial_stats["gt_1_rows"] or 0)
    )
    n_missing = int(initial_stats["null_rows"] or 0)

    out_of_range_fraction = out_of_range_count / total_before if total_before else 0.0
    is_valid_range = out_of_range_fraction <= invalid_fraction_cutoff

    if out_of_range_count and out_of_range_fraction > invalid_fraction_cutoff:
        emit.warn(
            "Column '%s' has %s values below 0 or above %s, which is %s of all %s variants - far more "
            "than the %s allowed by eaf.invalid_fraction_cutoff. The column may not be a frequency at all."
            % (
                colname,
                _fmt_int(out_of_range_count),
                tolerance,
                _fmt_pct(out_of_range_fraction),
                _fmt_int(total_before),
                _fmt_pct(invalid_fraction_cutoff),
            )
        )

    # ---- apply the out-of-range policy ---------------------------------
    cleaned_df = df
    if out_of_range_count and action == "fail":
        raise AlleleFrequencyError(
            "Column '%s' has %s frequency values below 0 or above %s and policy eaf.out_of_range is "
            "'fail'." % (colname, _fmt_int(out_of_range_count), tolerance)
        )
    if action == "reject":
        # fill_null(False): a null frequency is not out of range, it is missing,
        # and missingness is somebody else's decision
        cleaned_df = _reject_rows(
            cleaned_df,
            ((pl.col(colname) < 0.0) | (pl.col(colname) > tolerance)).fill_null(False),
            "eaf_out_of_range",
            emit,
            rejects=rejects,
            detail="%s outside 0 to %s" % (colname, tolerance),
            check_name="frequency range",
            plain_english=(
                "Frequencies below 0 or above %s cannot be real; policy eaf.out_of_range is 'reject', "
                "so those variants were removed." % tolerance
            ),
        )
    elif action == "null":
        before = cleaned_df.height
        cleaned_df = cleaned_df.with_columns(
            pl.when((pl.col(colname) < 0.0) | (pl.col(colname) > tolerance))
            .then(None)
            .otherwise(pl.col(colname))
            .alias(colname)
        )
        emit.qc(
            "frequency range",
            "%s values are outside 0 to %s; policy eaf.out_of_range is 'null', so those frequencies "
            "were blanked and the variants kept." % (_fmt_int(out_of_range_count), tolerance),
            before,
            cleaned_df.height,
            changed=out_of_range_count,
            warn=out_of_range_count > 0,
        )
    else:  # "clip" - the behaviour before this migration
        before = cleaned_df.height
        emit.qc(
            "frequency range",
            "%s values are outside 0 to 1; policy eaf.out_of_range is 'clip', so they were forced "
            "back into range." % _fmt_int(outside_unit_count),
            before,
            before,
            changed=outside_unit_count,
            warn=outside_unit_count > 0,
        )

    # every surviving value is squeezed into 0..1: this is a no-op for values
    # already in range and turns rounding error such as 1.02 into 1.0
    cleaned_df = cleaned_df.with_columns(
        pl.when(pl.col(colname).is_not_null())
        .then(pl.col(colname).clip(0.0, 1.0))
        .otherwise(None)
        .alias(colname)
    )

    # ---- missingness (null_count, not a row-count difference) ----------
    missing_fraction = n_missing / total_before if total_before else 0.0
    if missing_fraction > missing_fraction_cutoff:
        emit.warn(
            "Column '%s' is empty for %s of the %s variants (%s), more than the %s allowed by "
            "eaf.missing_fraction_cutoff."
            % (
                colname,
                _fmt_int(n_missing),
                _fmt_int(total_before),
                _fmt_pct(missing_fraction),
                _fmt_pct(missing_fraction_cutoff),
            )
        )

    # ---- MAF-likeness, measured over the NON-NULL values ---------------
    total_after = cleaned_df.height
    final_stats = summarise_allele_frequency(cleaned_df, colname, clip_tolerance=tolerance)
    low_freq = int(final_stats.get("le_0_5_rows") or 0)
    low_freq_pct = float(final_stats.get("le_0_5_fraction_among_non_null") or 0.0)
    non_null_rows = int(final_stats.get("non_null_rows") or 0)
    is_suspicious = low_freq_pct > maf_cutoff
    if is_suspicious:
        emit.warn(
            "Column '%s': %s of its %s non-empty values are at or below 0.5, more than the %s set by "
            "eaf.maf_decision_cutoff. That is what a minor allele frequency looks like, not an effect "
            "allele frequency."
            % (colname, _fmt_pct(low_freq_pct), _fmt_int(non_null_rows), _fmt_pct(maf_cutoff))
        )

    combined_stats = {
        "initial": initial_stats,
        "final": final_stats,
        "total_before_cleaning": total_before,
        "total_after_cleaning": total_after,
        "out_of_range_count": out_of_range_count,
        "out_of_range_fraction": out_of_range_fraction,
        "out_of_range_action": action,
        "outside_unit_interval_count": outside_unit_count,
        "is_valid_range": is_valid_range,
        # kept under the old key for the QC summary; it is now the honest null
        # count rather than a row-count difference that was always 0
        "missing_count_removed": n_missing,
        "missing_fraction_removed": missing_fraction,
        "missing_count": n_missing,
        "missing_fraction": missing_fraction,
        "low_freq_count_after_cleaning": low_freq,
        "low_freq_fraction_after_cleaning": low_freq_pct,
        "maf_like_flag": is_suspicious,
    }
    return is_valid_range, is_suspicious, out_of_range_count, n_missing, cleaned_df, low_freq_pct, combined_stats


# -------------------------------------------------------
# 7. Confirm a MAF-like study column against default EAF
# -------------------------------------------------------
def confirm_eaf_with_reference(
    df: pl.DataFrame,
    input_af_col: str,
    reference_file: str,
    population_col: str,
    colmap: Dict[str, str],
    std_cols: Dict[str, str],
    policies=None,
    logger=None,
    ctx=None,
    aligned_reference_col: Optional[str] = None,
    study_columns_canonical: bool = False,
) -> Tuple[str, Dict[str, Any]]:
    """Return ``eaf``, ``maf`` or ``inconclusive`` for a MAF-like column.

    The reference population column is ALT frequency. Direct matches use it as
    supplied; allele-swapped matches use ``1-AF`` so every comparison is against
    the frequency of the study's listed effect allele. Palindromic SNPs are
    excluded because their strand cannot be inferred from the allele pair.
    """
    emit = _Emitter(logger=logger, ctx=ctx)
    resolved = resolve_policies(policies)
    minimum = int(resolved.get("eaf.maf_reference_min_overlap"))
    maf_cutoff = float(resolved.get("eaf.maf_decision_cutoff"))
    keys = [std_cols[name] for name in ("chr", "pos", "ea", "oa")]
    allele_expressions = [pl.col(std_cols[name]) for name in ("ea", "oa")]
    if not study_columns_canonical:
        allele_expressions = [
            expression.cast(pl.Utf8).str.to_uppercase()
            for expression in allele_expressions
        ]
    pair = pl.concat_str(allele_expressions)
    use_aligned = aligned_reference_col in df.columns if aligned_reference_col else False
    selected = keys + [input_af_col]
    if use_aligned:
        selected.append(aligned_reference_col)
        if STRAND_ACTION_COLUMN in df.columns:
            selected.append(STRAND_ACTION_COLUMN)
    study = df.select(selected).filter(
        pl.col(input_af_col).is_not_null()
        & ~pair.is_in(list(PALINDROMIC_PAIRS)).fill_null(False)
    ).rename({input_af_col: "_study_af"})
    if use_aligned:
        reference_col = "__maf_expected_study_eaf"
        if STRAND_ACTION_COLUMN in study.columns:
            swapped = pl.col(STRAND_ACTION_COLUMN).is_in([
                "forward_swapped", "reverse_complement_swapped"
            ])
            merged = study.with_columns(
                pl.when(swapped)
                .then(1.0 - pl.col(aligned_reference_col))
                .otherwise(pl.col(aligned_reference_col))
                .alias(reference_col)
            )
        else:
            merged = study.with_columns(
                pl.col(aligned_reference_col).alias(reference_col)
            )
        matched_rows = merged.filter(pl.col(reference_col).is_not_null())
        matched_count = matched_rows.height
        if STRAND_ACTION_COLUMN in matched_rows.columns:
            actions = matched_rows.select([
                pl.col(STRAND_ACTION_COLUMN).is_in([
                    "forward", "reverse_complement"
                ]).sum().alias("direct"),
                pl.col(STRAND_ACTION_COLUMN).is_in([
                    "forward_swapped", "reverse_complement_swapped"
                ]).sum().alias("swapped"),
            ]).to_dicts()[0]
        else:
            actions = {"direct": matched_count, "swapped": 0}
        merge_stats = {
            "direct_match_rows": int(actions["direct"] or 0),
            "flip_match_rows": int(actions["swapped"] or 0),
        }
        emit.info(
            "Reused the population AF already aligned during reference orientation; "
            "the chromosome reference was not read a second time."
        )
    else:
        merged, reference_col, merge_stats = merge_external_allele_frequencies(
            study,
            reference_file,
            population_col,
            colmap,
            std_cols,
            policies=policies,
            logger=logger,
            ctx=ctx,
            min_match_fraction=0.0,
            palindromic_handling="ignore",
            study_columns_canonical=study_columns_canonical,
        )
    matched = merged.filter(pl.col(reference_col).is_not_null())

    stats = {
        "direct_matches": merge_stats["direct_match_rows"],
        "swapped_matches": merge_stats["flip_match_rows"],
        "comparable_variants": matched.height,
        "minimum_overlap": minimum,
        "reference_effect_allele_minor_fraction": None,
        "mean_eaf_error": None,
        "mean_maf_error": None,
    }
    comparison = matched.select([
        (
            (pl.col(reference_col) < 0.0)
            | (pl.col(reference_col) > 1.0)
        ).sum().alias("invalid_reference"),
        (pl.col("_study_af") - pl.col(reference_col)).abs().mean().alias("eaf"),
        (
            pl.col("_study_af")
            - pl.min_horizontal(
                pl.col(reference_col), 1.0 - pl.col(reference_col)
            )
        ).abs().mean().alias("maf"),
        (pl.col(reference_col) <= 0.5).mean().alias("minor_fraction"),
    ]).to_dicts()[0]
    invalid_reference = int(comparison["invalid_reference"] or 0)
    if invalid_reference:
        raise AlleleFrequencyError(
            "The default EAF table '%s' contains %s matched values outside 0 to 1 in column '%s'."
            % (reference_file, _fmt_int(invalid_reference), population_col)
        )
    if matched.height < minimum:
        emit.warn(
            "The MAF-like column '%s' had only %s non-palindromic matches to the default EAF table; "
            "%s are required by eaf.maf_reference_min_overlap."
            % (input_af_col, _fmt_int(matched.height), _fmt_int(minimum))
        )
        return "inconclusive", stats

    errors = {"eaf": comparison["eaf"], "maf": comparison["maf"]}
    stats["mean_eaf_error"] = float(errors["eaf"])
    stats["mean_maf_error"] = float(errors["maf"])
    reference_minor_fraction = comparison["minor_fraction"]
    stats["reference_effect_allele_minor_fraction"] = float(reference_minor_fraction)
    if reference_minor_fraction > maf_cutoff:
        decision = "eaf"
    else:
        decision = "maf" if errors["maf"] < errors["eaf"] else "eaf"
    emit.decide(
        "MAF-like frequency column '%s'" % input_af_col,
        {
            "non-palindromic reference matches": _fmt_int(matched.height),
            "mean error as EAF": errors["eaf"],
            "mean error as MAF": errors["maf"],
            "reference effect allele is minor": _fmt_pct(reference_minor_fraction),
        },
        "MAF not aligned to the effect allele" if decision == "maf" else "EAF aligned to the effect allele",
    )
    return decision, stats


# -------------------------------------------------------
# 8. Helper: Cross-Check against Default Reference
# -------------------------------------------------------
def resolve_palindromic_frequency(
    df: pl.DataFrame,
    eaf_col: str,
    ea_col: str,
    oa_col: str,
    study_eaf_col: Optional[str] = None,
    mode: str = "resolve_by_frequency",
    lower: float = 0.4,
    upper: float = 0.6,
    tolerance: float = 0.2,
    flag_column: str = PALINDROMIC_FLAG_COLUMN,
    emit: Optional["_Emitter"] = None,
) -> Tuple[pl.DataFrame, Dict[str, Any]]:
    """Flag, and where possible repair, palindromic variants after a panel match.

    An A/T or C/G variant complements to itself, so the ordered join key is
    symmetric and one of the two orientations always matches.  A strand error
    therefore assigns ``1-AF`` instead of ``AF`` silently - the dangerous case
    described in plan section 5.3.

    ``mode``:

    ``null``                 blank the frequency of every matched palindromic
                             variant and flag it ``palindromic_ambiguous``.
    ``reject``               resolve first when the study has its own frequency,
                             then leave the unresolvable ones flagged so the
                             caller can remove them.
    ``resolve_by_frequency`` the plan's arbitration, using the study's own
                             frequency::

                                 d_keep = |EAF_study - AF_matched|
                                 d_flip = |EAF_study - (1 - AF_matched)|

                             both inside [lower, upper]      -> null, ambiguous
                             d_keep <= d_flip, d_keep <= tol -> keep
                             d_flip <  d_keep, d_flip <= tol -> use 1-AF, flagged
                                                                strand_corrected
                             otherwise                       -> null, discordant

    Returns ``(frame_with_flag_column, stats)``.
    """
    flag_values = ("palindromic_kept", "strand_corrected", "palindromic_ambiguous",
                   "af_discordant", "palindromic_unresolved")
    stats: Dict[str, Any] = {
        "mode": mode,
        "rows": df.height,
        "palindromic_rows": 0,
        "resolved_with_study_frequency": False,
    }
    for value in flag_values:
        stats["flag_%s" % value] = 0
    if df.height == 0:
        return df.with_columns(pl.lit(None, dtype=pl.Utf8).alias(flag_column)), stats

    is_palindromic = (
        pl.concat_str([
            pl.col(ea_col).cast(pl.Utf8).str.to_uppercase(),
            pl.col(oa_col).cast(pl.Utf8).str.to_uppercase(),
        ])
        .is_in(list(PALINDROMIC_PAIRS))
        .fill_null(False)
    )

    can_resolve = bool(study_eaf_col) and study_eaf_col in df.columns and study_eaf_col != eaf_col
    stats["resolved_with_study_frequency"] = bool(can_resolve and mode in ("resolve_by_frequency", "reject"))

    if mode == "null" or (mode == "reject" and not can_resolve) or (
        mode == "resolve_by_frequency" and not can_resolve
    ):
        if mode == "resolve_by_frequency":
            # nothing to arbitrate with: leave the values alone but say so loudly
            flag_expr = (
                pl.when(is_palindromic)
                .then(pl.lit("palindromic_unresolved"))
                .otherwise(pl.lit(None, dtype=pl.Utf8))
                .alias(flag_column)
            )
            out = df.with_columns(flag_expr)
            if emit is not None:
                emit.warn(
                    "eaf.palindromic_handling is 'resolve_by_frequency', but this study has no frequency "
                    "column of its own to arbitrate with, so palindromic variants were left exactly as the "
                    "panel matched them. A strand error would be invisible here - see plan section 5.4."
                )
        else:
            flag_expr = (
                pl.when(is_palindromic)
                .then(pl.lit("palindromic_ambiguous"))
                .otherwise(pl.lit(None, dtype=pl.Utf8))
                .alias(flag_column)
            )
            out = df.with_columns(flag_expr).with_columns(
                pl.when(pl.col(flag_column) == "palindromic_ambiguous")
                .then(None)
                .otherwise(pl.col(eaf_col))
                .alias(eaf_col)
            )
    else:
        study = pl.col(study_eaf_col).cast(pl.Float64)
        matched = pl.col(eaf_col).cast(pl.Float64)
        d_keep = (study - matched).abs()
        d_flip = (study - (1.0 - matched)).abs()
        both_in_band = (
            study.is_between(lower, upper) & matched.is_between(lower, upper)
        ).fill_null(False)
        unusable = study.is_null() | matched.is_null()
        flag_expr = (
            pl.when(~is_palindromic).then(pl.lit(None, dtype=pl.Utf8))
            .when(unusable).then(pl.lit("palindromic_ambiguous"))
            .when(both_in_band).then(pl.lit("palindromic_ambiguous"))
            .when((d_keep <= d_flip) & (d_keep <= tolerance)).then(pl.lit("palindromic_kept"))
            .when((d_flip < d_keep) & (d_flip <= tolerance)).then(pl.lit("strand_corrected"))
            .otherwise(pl.lit("af_discordant"))
            .alias(flag_column)
        )
        out = df.with_columns(flag_expr)
        out = out.with_columns(
            pl.when(pl.col(flag_column) == "strand_corrected")
            .then(1.0 - pl.col(eaf_col).cast(pl.Float64))
            .when(pl.col(flag_column).is_in(["palindromic_ambiguous", "af_discordant"]))
            .then(None)
            .otherwise(pl.col(eaf_col))
            .alias(eaf_col)
        )

    counted = out.select([
        is_palindromic.sum().alias("palindromic_rows"),
    ]).to_dicts()[0]
    stats["palindromic_rows"] = int(counted["palindromic_rows"] or 0)
    value_counts = out.get_column(flag_column).value_counts()
    for row in value_counts.to_dicts():
        label = row.get(flag_column)
        if label:
            stats["flag_%s" % label] = int(row.get("count") or 0)
    return out, stats


# -------------------------------------------------------
# 9. Helper: Load Generic External
# ------------------------------------------------------
def merge_external_allele_frequencies(
    df: pl.DataFrame,
    path: str,
    eaf_col: str,
    colmap: dict,
    std_cols: dict,
    chromosome: str = "",
    policies=None,
    logger=None,
    ctx=None,
    rejects=None,
    study_eaf_col: Optional[str] = None,
    min_match_fraction: Optional[float] = None,
    palindromic_handling: Optional[str] = None,
    study_columns_canonical: bool = False,
) -> Tuple[pl.DataFrame, str, Dict[str, Any]]:
    """Take the frequency from an external panel.

    Variants are matched in the orientations listed by ``eaf.match_orientations``
    (``direct`` then ``swap``); complement passes are deliberately not offered,
    see plan section 5.3. Unmatched variants keep a null frequency. The coverage
    gate counts only allele-key matches whose frequency is finite and between 0
    and 1; if that usable share falls below
    ``eaf.external_min_match_fraction`` the step fails instead of quietly
    reporting success.
    """
    emit = _Emitter(logger=logger, ctx=ctx)
    min_match_fraction = float(_pol(policies, "eaf.external_min_match_fraction", min_match_fraction))
    deduplicate_reference = bool(_pol(policies, "eaf.deduplicate_reference"))
    palindromic_handling = str(_pol(policies, "eaf.palindromic_handling", palindromic_handling))
    orientations = [
        str(value).lower()
        for value in resolve_policies(policies).get("eaf.match_orientations")
    ]
    unsupported = [o for o in orientations if o not in ("direct", "swap")]
    if unsupported:
        raise AlleleFrequencyError(
            "eaf.match_orientations asks for %s, but only 'direct' and 'swap' are implemented. Strand is "
            "a property of the study, not of individual variants (plan section 5.3): use dataset-level "
            "strand detection instead of complement passes." % unsupported
        )
    stats: Dict[str, Any] = {
        "external_rows_loaded": 0,
        "input_rows_before_merge": df.height,
        "direct_match_rows": 0,
        "flip_match_rows": 0,
        "key_match_rows": 0,
        "key_match_fraction": 0.0,
        "usable_direct_match_rows": 0,
        "usable_flip_match_rows": 0,
        "usable_match_rows": 0,
        "usable_match_fraction": 0.0,
        "matched_missing_frequency_rows": 0,
        "matched_non_finite_frequency_rows": 0,
        "matched_out_of_range_frequency_rows": 0,
        "matched_unusable_frequency_rows": 0,
        "unusable_frequency_rows": 0,
        "unmatched_rows": 0,
        "total_missing": 0,
        "missing_fraction": 0.0,
        "match_fraction": 0.0,
        "min_match_fraction": min_match_fraction,
        "orientations": orientations,
        "palindromic_handling": palindromic_handling,
    }
    # --------------------------------------------------
    # Load file
    # --------------------------------------------------
    load_path = path
    if not os.path.exists(load_path):
        raise AlleleFrequencyError(
            "The external effect-allele-frequency file was not found: '%s'. Check the 'eaffile' setting "
            "for this dataset." % load_path
        )
    ext_df = read_reference_table(
        load_path,
        policies,
        colmap["delimiter"],
        what="external frequency file",
        emit=emit,
    )
    stats["external_rows_loaded"] = ext_df.height
    # --------------------------------------------------
    # Validate columns
    # --------------------------------------------------
    missing_cols = []
    for key in ["chr", "pos", "a1", "a2"]:
        if colmap[key] not in ext_df.columns:
            missing_cols.append(colmap[key])
    if eaf_col not in ext_df.columns:
        missing_cols.append(eaf_col)
    if missing_cols:
        raise AlleleFrequencyError(
            "The external frequency file '%s' does not contain the columns %s. It has: %s"
            % (load_path, missing_cols, ", ".join(ext_df.columns))
        )
    matched_all, orientation_col, join_stats = allele_oriented_left_join(
        df,
        ext_df,
        study_columns=std_cols,
        reference_columns={
            "chr": colmap["chr"],
            "pos": colmap["pos"],
            "ea": colmap["a1"],
            "oa": colmap["a2"],
        },
        value_column=eaf_col,
        output_column=eaf_col,
        orientations=orientations,
        swapped_value="one_minus",
        deduplicate_reference=deduplicate_reference,
        prefer_non_null_value=True,
        study_columns_canonical=study_columns_canonical,
        error_type=AlleleFrequencyError,
        warn=emit.warn,
        reference_label="external frequency file '%s'" % load_path,
        duplicate_policy_name="eaf.deduplicate_reference",
    )
    duplicate_rows = join_stats["reference_duplicate_rows_removed"]
    preferred_usable = join_stats[
        "reference_duplicate_groups_preferred_usable_value"
    ]
    if duplicate_rows:
        emit.info(
            "External frequency reference duplicates: %s repeated rows removed; "
            "%s allele-specific keys used a later finite value because the first "
            "value was missing or non-finite. Conflicting finite values are not "
            "selected and stop the chromosome."
            % (_fmt_int(duplicate_rows), _fmt_int(preferred_usable))
        )
    stats.update(join_stats)

    # An allele-key match is only evidence that the panel contains the variant;
    # it does not supply an EAF when the matched value is null, non-finite or
    # outside the probability interval. Count all value states in one Polars
    # aggregation so the configured coverage gate measures usable frequencies.
    reference_value = pl.col(eaf_col).cast(pl.Float64, strict=False)
    has_key_match = pl.col(orientation_col).is_not_null()
    finite_value = reference_value.is_finite().fill_null(False)
    in_frequency_range = (
        (reference_value >= 0.0) & (reference_value <= 1.0)
    ).fill_null(False)
    usable_value = has_key_match & finite_value & in_frequency_range
    value_counts = matched_all.select([
        ((pl.col(orientation_col) == 0) & usable_value).sum().alias("direct_usable"),
        ((pl.col(orientation_col) == 1) & usable_value).sum().alias("swapped_usable"),
        usable_value.sum().alias("usable"),
        (has_key_match & reference_value.is_null()).sum().alias("matched_missing"),
        (
            has_key_match & reference_value.is_not_null() & ~finite_value
        ).sum().alias("matched_non_finite"),
        (
            has_key_match & finite_value & ~in_frequency_range
        ).sum().alias("matched_out_of_range"),
    ]).row(0, named=True)

    grouped = matched_all.partition_by(
        orientation_col, as_dict=True, maintain_order=True,
    )
    empty = matched_all.head(0)
    df1 = grouped.get((0,), empty).drop(orientation_col)
    df2 = grouped.get((1,), empty).drop(orientation_col)
    unmatched_df = grouped.get((None,), empty).drop(orientation_col)
    stats["direct_match_rows"] = join_stats["direct_key_matches"]
    stats["flip_match_rows"] = join_stats["swapped_key_matches"]
    stats["key_match_rows"] = stats["direct_match_rows"] + stats["flip_match_rows"]
    stats["key_match_fraction"] = (
        stats["key_match_rows"] / df.height if df.height else 0.0
    )
    stats["usable_direct_match_rows"] = int(value_counts["direct_usable"] or 0)
    stats["usable_flip_match_rows"] = int(value_counts["swapped_usable"] or 0)
    stats["usable_match_rows"] = int(value_counts["usable"] or 0)
    stats["usable_match_fraction"] = (
        stats["usable_match_rows"] / df.height if df.height else 0.0
    )
    stats["matched_missing_frequency_rows"] = int(
        value_counts["matched_missing"] or 0
    )
    stats["matched_non_finite_frequency_rows"] = int(
        value_counts["matched_non_finite"] or 0
    )
    stats["matched_out_of_range_frequency_rows"] = int(
        value_counts["matched_out_of_range"] or 0
    )
    stats["matched_unusable_frequency_rows"] = (
        stats["matched_missing_frequency_rows"]
        + stats["matched_non_finite_frequency_rows"]
        + stats["matched_out_of_range_frequency_rows"]
    )

    # ==================================================
    # STEP 3: UNMATCHED -> NULL EAF from the left join
    # ==================================================
    stats["unmatched_rows"] = unmatched_df.height
    stats["total_missing"] = (
        stats["unmatched_rows"] + stats["matched_missing_frequency_rows"]
    )
    stats["missing_fraction"] = (
        stats["total_missing"] / df.height if df.height else 0.0
    )
    stats["unusable_frequency_rows"] = df.height - stats["usable_match_rows"]
    # Keep the established public name, but make its meaning scientifically
    # correct: it is now the usable EAF share, not merely the allele-key share.
    stats["match_fraction"] = stats["usable_match_fraction"]
    # ==================================================
    # STEP 4: PALINDROMIC VARIANTS (plan section 5.3)
    # ==================================================
    if palindromic_handling != "ignore":
        resolved = resolve_policies(policies)
        lower = float(resolved.get("eaf.palindromic_ambiguous_lower"))
        upper = float(resolved.get("eaf.palindromic_ambiguous_upper"))
        tolerance = float(resolved.get("eaf.palindromic_resolve_tolerance"))
        df1, pal1 = resolve_palindromic_frequency(
            df1, eaf_col, std_cols["ea"], std_cols["oa"],
            study_eaf_col=study_eaf_col, mode=palindromic_handling,
            lower=lower, upper=upper, tolerance=tolerance, emit=emit,
        )
        df2, pal2 = resolve_palindromic_frequency(
            df2, eaf_col, std_cols["ea"], std_cols["oa"],
            study_eaf_col=study_eaf_col, mode=palindromic_handling,
            lower=lower, upper=upper, tolerance=tolerance, emit=None,
        )
        unmatched_df = unmatched_df.with_columns(
            pl.lit(None, dtype=pl.Utf8).alias(PALINDROMIC_FLAG_COLUMN)
        )
        palindromic_stats = {}
        for key in set(list(pal1.keys()) + list(pal2.keys())):
            left, right = pal1.get(key), pal2.get(key)
            if isinstance(left, int) and isinstance(right, int):
                palindromic_stats[key] = left + right
            else:
                palindromic_stats[key] = left if left is not None else right
        stats["palindromic"] = palindromic_stats
    # ==================================================
    # FINAL CONCAT (preserve all input variants)
    # ==================================================
    final_df = pl.concat([df1, df2, unmatched_df], how="vertical")

    if palindromic_handling != "ignore":
        pal = stats.get("palindromic", {})
        emit.info(
            "Palindromic variants (A/T or C/G, which look the same on either strand): %s of the matched "
            "variants. Policy eaf.palindromic_handling is '%s'."
            % (_fmt_int(pal.get("palindromic_rows", 0)), palindromic_handling)
        )
        if pal.get("flag_strand_corrected"):
            emit.warn(
                "%s palindromic variants matched the panel better after flipping the frequency, so "
                "1-AF was used and the rows were flagged '%s = strand_corrected'."
                % (_fmt_int(pal["flag_strand_corrected"]), PALINDROMIC_FLAG_COLUMN)
            )
    if palindromic_handling == "reject":
        # fill_null(False): the flag is null for every non-palindromic variant
        final_df = _reject_rows(
            final_df,
            (pl.col(PALINDROMIC_FLAG_COLUMN) == "palindromic_ambiguous").fill_null(False),
            "palindromic_ambiguous",
            emit,
            rejects=rejects,
            detail="A/T or C/G variant that could not be oriented",
            check_name="palindromic variants",
            plain_english=(
                "An A/T or C/G variant looks identical on either strand. These could not be oriented, "
                "so policy eaf.palindromic_handling='reject' removed them."
            ),
        )
        final_df = _reject_rows(
            final_df,
            (pl.col(PALINDROMIC_FLAG_COLUMN) == "af_discordant").fill_null(False),
            "af_discordant",
            emit,
            rejects=rejects,
            detail="study frequency far from the panel in both orientations",
            check_name="palindromic variants, frequency discordant",
            plain_english=(
                "The study frequency was far from the panel frequency in both orientations, so neither "
                "reading can be trusted and the variants were removed."
            ),
        )

    emit.info("Direct allele-key matches  : %s" % _fmt_int(stats["direct_match_rows"]))
    emit.info("Swapped allele-key matches : %s" % _fmt_int(stats["flip_match_rows"]))
    emit.info("No allele-key match        : %s" % _fmt_int(stats["unmatched_rows"]))
    emit.info(
        "Usable external frequencies : %s of %s variants (%s)."
        % (
            _fmt_int(stats["usable_match_rows"]),
            _fmt_int(df.height),
            _fmt_pct(stats["usable_match_fraction"]),
        )
    )
    if stats["matched_unusable_frequency_rows"]:
        emit.warn(
            "%s allele-key matches did not supply a usable EAF: %s missing, %s "
            "non-finite and %s outside 0 to 1."
            % (
                _fmt_int(stats["matched_unusable_frequency_rows"]),
                _fmt_int(stats["matched_missing_frequency_rows"]),
                _fmt_int(stats["matched_non_finite_frequency_rows"]),
                _fmt_int(stats["matched_out_of_range_frequency_rows"]),
            )
        )

    if df.height and stats["match_fraction"] < min_match_fraction:
        raise AlleleFrequencyError(
            "The external frequency file '%s' matched chromosome, position and alleles for %s of %s "
            "variants on chromosome %s (%s), but supplied a usable EAF for only %s variants (%s). "
            "A usable EAF is finite and between 0 and 1. This is below the %s required by "
            "eaf.external_min_match_fraction. Among allele-key matches, %s frequencies were missing, "
            "%s were non-finite and %s were outside 0 to 1; %s variants had no allele-key match. "
            "Check the reference panel's genome build, chromosome and allele columns, frequency "
            "column, and completeness, or explicitly lower the threshold for an intentionally "
            "limited panel."
            % (
                load_path,
                _fmt_int(stats["key_match_rows"]),
                _fmt_int(df.height),
                chromosome,
                _fmt_pct(stats["key_match_fraction"]),
                _fmt_int(stats["usable_match_rows"]),
                _fmt_pct(stats["usable_match_fraction"]),
                _fmt_pct(min_match_fraction),
                _fmt_int(stats["matched_missing_frequency_rows"]),
                _fmt_int(stats["matched_non_finite_frequency_rows"]),
                _fmt_int(stats["matched_out_of_range_frequency_rows"]),
                _fmt_int(stats["unmatched_rows"]),
            )
        )
    return final_df, eaf_col, stats


# -------------------------------------------------------
# 10. Helper: Write structured statistics into qc_info
# -------------------------------------------------------
def _add_prefixed_statistics(
    qc_info: Dict[str, Any],
    prefix: str,
    stats: Dict[str, Any]
) -> None:
    for key, value in stats.items():
        if isinstance(value, dict):
            for sub_key, sub_val in value.items():
                qc_info["%s_%s_%s" % (prefix, key, sub_key)] = sub_val
        else:
            qc_info["%s_%s" % (prefix, key)] = value


def _out_of_bound_frame(frame: pl.DataFrame, colname: str) -> pl.DataFrame:
    """Rows whose frequency lies outside 0..1, taken BEFORE any clipping."""
    if colname not in frame.columns:
        return frame.head(0)
    # Cast defensively: a frequency column read with schema inference disabled
    # arrives as Utf8, and comparing that against a float raises ComputeError.
    # Every sibling module casts with strict=False first; match them.
    value = pl.col(colname).cast(pl.Float64, strict=False)
    return frame.filter(((value < 0.0) | (value > 1.0)).fill_null(False))


# -------------------------------------------------------
# 11. MAIN FUNCTION
# -------------------------------------------------------
def harmonise_allele_frequency(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    output_layout: Dict[str, str],
    output_delimiter: str,
    eaffile: Optional[str] = None,
    external_eaf_colmap: Optional[Dict] = None,
    default_eaf_file: Optional[str] = None,
    default_eaf_column: Optional[str] = None,
    default_eaf_colmap: Optional[Dict] = None,
    policies=None,
    logger=None,
    rejects=None,
    step_context=None,
    step_number: int = 4,
    step_total: int = 16,
    study_decision: Optional[Any] = None,
) -> Tuple[pl.DataFrame, Dict, dict]:
    """Give every variant an effect allele frequency.

    Parameters added by the plan v3 migration, all optional:

    policies       a ``Policies`` object; ``None`` means the registry defaults.
    logger         a ``PipelineLogger``.  When given, this function opens its own
                   ``step()`` and nothing is written to stdout.
    rejects        a ``RejectCollector``.  Removed variants are written to the
                   per-chromosome rejected-variants file with their reason.
    step_context   an already-open ``StepContext``, when the caller prefers to
                   own the step.
    study_decision the dataset-level decision from step 07, e.g.
                   ``{"eaf_is_maf": False}``. The decision is made once from
                   the whole study and then applied consistently to each chromosome.
    default_eaf_*  the raw chromosome REF/ALT/population-AF reference used for
                   strand orientation, and reused for a MAF-like study column.
    """
    gwas_outputname = sample_column_dict["gwas_outputname"]
    output_dir = sample_column_dict["output_folder"]
    merge_stats = {}

    pol = resolve_policies(policies)
    maf_cutoff = float(pol.get("eaf.maf_decision_cutoff"))
    invalid_fraction_cutoff = float(pol.get("eaf.invalid_fraction_cutoff"))
    missing_fraction_cutoff = float(pol.get("eaf.missing_fraction_cutoff"))
    degenerate_action = str(pol.get("eaf.degenerate"))
    inconclusive_action = str(pol.get("eaf.on_maf_check_inconclusive"))

    # Create subdirectory for EAF QC outputs
    eaf_qc_dir = configured_output_path(
        output_dir, output_layout["frequency_qc_directory"],
    )
    eaf_qc_dir.mkdir(parents=True, exist_ok=True)

    qc_info: Dict[str, Any] = {
        "chromosome": chromosome,
        "gwas_outputname": gwas_outputname,
        "Total_number_of_variants": df.height,
        "final_status": "started",
        "decision_source": None,
        "final_eaf_col": None,
    }

    final_df = None
    with contextlib.ExitStack() as stack:
        ctx = step_context
        if ctx is None and logger is not None:
            ctx = stack.enter_context(
                logger.step(
                    step_number,
                    step_total,
                    STEP_TITLE,
                    FUNC_NAME,
                    rows_in=df.height,
                    policy_keys=list(STRAND_POLICY_KEYS) + list(POLICY_KEYS),
                )
            )
        emit = _Emitter(logger=logger, ctx=ctx)
        emit.info(
            "Giving every variant on chromosome %s an effect allele frequency. Starting with %s variants."
            % (chromosome, _fmt_int(df.height))
        )

        try:
            std_cols = {
                "chr": sample_column_dict["chr_col"],
                "pos": sample_column_dict["pos_col"],
                "ea": sample_column_dict["ea_col"],
                "oa": sample_column_dict["oa_col"],
            }

            df, strand_qc, sample_column_dict = harmonise_strand_orientation(
                chromosome=chromosome,
                df=df,
                sample_column_dict=sample_column_dict,
                reference_file=default_eaf_file,
                population_col=default_eaf_column,
                reference_colmap=default_eaf_colmap,
                study_decision=study_decision,
                policies=pol,
                logger=logger,
                ctx=ctx,
                rejects=rejects,
            )
            qc_info["strand_orientation"] = strand_qc

            internal_eaf = optional_text(sample_column_dict.get("eaf_col"))

            if external_eaf_colmap is None:
                raise AlleleFrequencyError(
                    "External EAF column mapping is missing from the resolved YAML configuration."
                )

            final_eaf_col = None
            out_of_bound_df = df.head(0)
            study_says_maf = _study_decision_flag(study_decision)
            qc_info["study_decision_eaf_is_maf"] = study_says_maf

            # ------------------------------
            # STEP 1: Internal EAF
            # ------------------------------
            if internal_eaf and internal_eaf in df.columns:
                emit.info("The study supplies its own frequency column '%s'." % internal_eaf)

                # captured BEFORE anything is clipped, so the count is real
                out_of_bound_df = _out_of_bound_frame(df, internal_eaf)

                valid_range, is_suspicious, out_of_range, n_missing, cleaned_df, low_freq_pct, stats = \
                    validate_allele_frequency(
                        df,
                        internal_eaf,
                        maf_cutoff,
                        invalid_fraction_cutoff,
                        missing_fraction_cutoff,
                        policies=pol,
                        logger=logger,
                        ctx=ctx,
                        rejects=rejects,
                    )

                qc_info["Internal_AF"] = "present"
                _add_prefixed_statistics(qc_info, "internal_af", stats)
                qc_info["internal_af_valid_range"] = valid_range

                # NOTE: `valid_range` is reported, not used as the gate for
                # adopting the column.  eaf.out_of_range has already decided
                # what happens to out-of-range values; diverting an otherwise
                # usable column to the reference panel because a few values were
                # bad would change behaviour silently.
                if cleaned_df.height == 0:
                    raise AlleleFrequencyError(
                        "No variants remain after validating effect-allele-frequency column "
                        "'%s'." % internal_eaf
                    )
                if study_says_maf is not None:
                    emit.decide(
                        "Frequency column '%s'" % internal_eaf,
                        "study-level decision, made once on the whole file before the split",
                        "a minor allele frequency" if study_says_maf else "an effect allele frequency",
                    )
                    if study_says_maf:
                        reference_file = optional_text(default_eaf_file)
                        reference_column = optional_text(default_eaf_column)
                        if reference_file is None or reference_column is None or default_eaf_colmap is None:
                            raise AlleleFrequencyError(
                                "Column '%s' is MAF-like, but the resolved default EAF file, population "
                                "column, or column mapping is missing." % internal_eaf
                            )
                        reference_decision, reference_stats = confirm_eaf_with_reference(
                            cleaned_df,
                            internal_eaf,
                            reference_file,
                            reference_column,
                            default_eaf_colmap,
                            std_cols,
                            policies=pol,
                            logger=logger,
                            ctx=ctx,
                            aligned_reference_col=REFERENCE_AF_COLUMN,
                            study_columns_canonical=True,
                        )
                        _add_prefixed_statistics(qc_info, "maf_reference", reference_stats)
                        qc_info["maf_reference_decision"] = reference_decision
                        if reference_decision == "maf":
                            raise AlleleFrequencyError(
                                "Column '%s' was confirmed as minor allele frequency rather than frequency "
                                "aligned to the listed effect allele. Provide an effect-allele-frequency "
                                "column or a dataset-supplied external EAF file." % internal_eaf
                            )
                        if reference_decision == "inconclusive":
                            qc_info["maf_check_inconclusive_action"] = inconclusive_action
                            if inconclusive_action == "fail":
                                raise AlleleFrequencyError(
                                    "Column '%s' is MAF-like, but too few variants matched the default EAF "
                                    "table to determine whether it is aligned to the listed effect allele. "
                                    "Policy eaf.on_maf_check_inconclusive is 'fail'." % internal_eaf
                                )
                            emit.warn(
                                "The default EAF comparison was inconclusive, but policy "
                                "eaf.on_maf_check_inconclusive is 'accept', so column '%s' is being used."
                                % internal_eaf
                            )
                            qc_info["decision_source"] = "internal_eaf_reference_inconclusive"
                        else:
                            if STRAND_ACTION_COLUMN in cleaned_df.columns:
                                swapped = pl.col(STRAND_ACTION_COLUMN).is_in([
                                    "forward_swapped", "reverse_complement_swapped"
                                ])
                                cleaned_df = cleaned_df.with_columns(
                                    pl.when(swapped)
                                    .then(1.0 - pl.col(internal_eaf))
                                    .otherwise(pl.col(internal_eaf))
                                    .alias(internal_eaf)
                                )
                            qc_info["decision_source"] = "internal_eaf_reference_validation"
                    else:
                        qc_info["decision_source"] = "internal_eaf_study_decision"
                    final_df = cleaned_df
                    final_eaf_col = internal_eaf
                elif not is_suspicious:
                    final_df = cleaned_df
                    final_eaf_col = internal_eaf
                    qc_info["decision_source"] = "internal_eaf_direct"
                else:
                    qc_info["maf_check_inconclusive_action"] = inconclusive_action
                    if inconclusive_action == "fail":
                        qc_info["final_status"] = "failed"
                        raise AlleleFrequencyError(
                            "Column '%s' has a frequency distribution consistent with a minor "
                            "allele frequency, but no study-wide decision is available. Policy "
                            "eaf.on_maf_check_inconclusive is %r."
                            % (internal_eaf, inconclusive_action)
                        )
                    emit.warn(
                        "Column '%s' has a frequency distribution consistent with a minor "
                        "allele frequency, but policy eaf.on_maf_check_inconclusive is "
                        "'accept', so it is being used as effect-allele frequency."
                        % internal_eaf
                    )
                    final_df = cleaned_df
                    final_eaf_col = internal_eaf
                    qc_info["decision_source"] = "internal_eaf_unverified"
            else:
                qc_info["Internal_AF"] = "Absent"
                emit.info("The study has no frequency column of its own.")

            # ------------------------------
            # STEP 2: Dataset-supplied external source
            # ------------------------------
            if final_df is None:
                emit.info(
                    "Taking the effect allele frequency from the dataset's configured "
                    "external frequency file."
                )

                external_file = optional_text(eaffile)
                if external_file is None or not os.path.exists(external_file):
                    raise AlleleFrequencyError(
                        "No usable effect-allele-frequency source remains. Provide either "
                        "effect_allele_frequency_column or external_eaf_file plus "
                        "external_eaf_column in the sample sheet."
                    )
                target_file = external_file
                target_col = optional_text(sample_column_dict.get("eafcolumn"))
                if target_col is None:
                    raise AlleleFrequencyError("external_eaf_column is required with external_eaf_file.")
                target_map = external_eaf_colmap
                qc_info["external_eaf_file"] = os.path.basename(target_file)
                qc_info["external_eaf_column"] = target_col

                study_frequency_col = internal_eaf if (internal_eaf and internal_eaf in df.columns) else None

                merged_df, new_col_name, merge_stats = merge_external_allele_frequencies(
                    df,
                    target_file,
                    target_col,
                    target_map,
                    std_cols,
                    chromosome,
                    policies=pol,
                    logger=logger,
                    ctx=ctx,
                    rejects=rejects,
                    study_eaf_col=study_frequency_col,
                    study_columns_canonical=True,
                )

                _add_prefixed_statistics(qc_info, "external_merge", merge_stats)

                # captured BEFORE the clip inside validate_allele_frequency
                out_of_bound_df = _out_of_bound_frame(merged_df, new_col_name)

                _, _, _, _, cleaned_df, _, stats = validate_allele_frequency(
                    merged_df,
                    new_col_name,
                    maf_cutoff,
                    invalid_fraction_cutoff,
                    missing_fraction_cutoff,
                    policies=pol,
                    logger=logger,
                    ctx=ctx,
                    rejects=rejects,
                )

                _add_prefixed_statistics(qc_info, "external_af", stats)

                final_df = cleaned_df
                final_eaf_col = new_col_name
                if qc_info["decision_source"] != "internal_eaf_rejected_inconclusive":
                    qc_info["decision_source"] = "external_eaf"

            # ------------------------------
            # FINALIZE
            # ------------------------------
            if final_df is None or final_eaf_col is None:
                raise AlleleFrequencyError("Failed to determine which column holds the effect allele frequency.")

            missing_eaf_df = final_df.filter(pl.col(final_eaf_col).is_null())
            missing_eaf_count = missing_eaf_df.height
            out_of_bound_count = out_of_bound_df.height

            # Save QC files inside subdirectory
            missing_file = configured_output_path(
                output_dir,
                output_layout["missing_eaf"],
                dataset_id=gwas_outputname,
                chromosome=chromosome,
            )
            outofbound_file = configured_output_path(
                output_dir,
                output_layout["out_of_range_eaf"],
                dataset_id=gwas_outputname,
                chromosome=chromosome,
            )

            if missing_eaf_count > 0:
                missing_eaf_df.write_csv(missing_file, separator=output_delimiter)
                emit.info("Variants with no frequency were written to %s" % missing_file)

            if out_of_bound_count > 0:
                out_of_bound_df.write_csv(outofbound_file, separator=output_delimiter)
                emit.info("Variants whose frequency was outside 0 to 1 were written to %s" % outofbound_file)

            qc_info["missing_eaf_count"] = missing_eaf_count
            qc_info["out_of_bound_eaf_count"] = out_of_bound_count

            # Reference AF is supporting concordance evidence after alleles and
            # the final EAF source are aligned. It must not arbitrate a column
            # that remains only provisionally accepted as possible MAF.
            frequency_is_aligned = (
                not study_says_maf
                or qc_info.get("maf_reference_decision") == "eaf"
            )
            if REFERENCE_AF_COLUMN in final_df.columns and frequency_is_aligned:
                af_difference_col = "strand_af_difference"
                tolerance = float(pol.get("strand.af_tolerance"))
                discordance_action = str(pol.get("strand.af_discordance_action"))
                final_df = final_df.with_columns(
                    (
                        pl.col(final_eaf_col).cast(pl.Float64, strict=False)
                        - pl.col(REFERENCE_AF_COLUMN).cast(pl.Float64, strict=False)
                    ).abs().alias(af_difference_col)
                )
                comparable = (
                    pl.col(final_eaf_col).is_not_null()
                    & pl.col(REFERENCE_AF_COLUMN).is_not_null()
                )
                discordant = (comparable & (pl.col(af_difference_col) > tolerance)).fill_null(False)
                counts = final_df.select([
                    comparable.sum().alias("comparable"),
                    discordant.sum().alias("discordant"),
                ]).to_dicts()[0]
                n_comparable = int(counts["comparable"] or 0)
                n_discordant = int(counts["discordant"] or 0)
                qc_info["strand_af_comparable"] = n_comparable
                qc_info["strand_af_discordant"] = n_discordant
                qc_info["strand_af_tolerance"] = tolerance
                qc_info["strand_af_discordance_action"] = discordance_action
                if n_discordant and discordance_action == "fail":
                    raise AlleleFrequencyError(
                        "%s variants differ from the aligned population reference AF by more "
                        "than %s and strand.af_discordance_action is 'fail'."
                        % (_fmt_int(n_discordant), tolerance)
                    )
                if discordance_action == "reject":
                    final_df = _reject_rows(
                        final_df,
                        discordant,
                        "af_discordant",
                        emit,
                        rejects=rejects,
                        detail="absolute EAF difference above %s" % tolerance,
                        check_name="study/reference frequency concordance",
                        plain_english=(
                            "The aligned study EAF differs from the configured population "
                            "reference by more than %s." % tolerance
                        ),
                    )
                else:
                    emit.qc(
                        "study/reference frequency concordance",
                        "%s of %s comparable variants differ from the aligned population "
                        "reference AF by more than %s; policy action is '%s'."
                        % (
                            _fmt_int(n_discordant), _fmt_int(n_comparable),
                            tolerance, discordance_action,
                        ),
                        final_df.height,
                        final_df.height,
                        changed=n_discordant,
                        warn=n_discordant > 0,
                    )
            else:
                qc_info["strand_af_comparable"] = 0
                qc_info["strand_af_discordant"] = 0

            final_df = final_df.with_columns([
                pl.when(pl.col(final_eaf_col).is_null())
                .then(None)
                .when(pl.col(final_eaf_col) <= 0.5)
                .then(pl.col(final_eaf_col))
                .otherwise(1 - pl.col(final_eaf_col))
                .alias("zmaf"),

                pl.col(final_eaf_col).clip(0.0, 1.0).alias(final_eaf_col),
            ])

            # ---- frequencies of exactly 0 or 1 ----------------------------
            # These make the standard-error denominator 2*p*(1-p)*(Neff+z^2)
            # exactly zero, which yields inf rather than null and slips through
            # every later gate.  Clipping is what creates them, so this has to
            # run after the clip above.
            # fill_null(False): a missing frequency is not a degenerate one
            degenerate_expr = (
                (pl.col(final_eaf_col) == 0.0) | (pl.col(final_eaf_col) == 1.0)
            ).fill_null(False)
            n_degenerate = int(final_df.select(degenerate_expr.sum()).item() or 0)
            qc_info["eaf_degenerate_count"] = n_degenerate
            qc_info["eaf_degenerate_action"] = degenerate_action
            if degenerate_action == "fail" and n_degenerate:
                raise AlleleFrequencyError(
                    "%s variants have a frequency of exactly 0 or 1 and policy eaf.degenerate is 'fail'."
                    % _fmt_int(n_degenerate)
                )
            if degenerate_action == "reject":
                before_degenerate = final_df.height
                final_df = _reject_rows(
                    final_df,
                    degenerate_expr,
                    "eaf_degenerate",
                    emit,
                    rejects=rejects,
                    detail="%s is exactly 0 or 1" % final_eaf_col,
                    check_name="degenerate frequencies",
                    plain_english=(
                        "A frequency of exactly 0 or 1 makes the standard-error formula divide by zero, "
                        "so those variants were removed."
                    ),
                )
                qc_info["eaf_degenerate_removed"] = before_degenerate - final_df.height
            elif degenerate_action == "null":
                before_degenerate = final_df.height
                final_df = final_df.with_columns([
                    pl.when(degenerate_expr).then(None).otherwise(pl.col(final_eaf_col)).alias(final_eaf_col),
                    pl.when(degenerate_expr).then(None).otherwise(pl.col("zmaf")).alias("zmaf"),
                ])
                emit.qc(
                    "degenerate frequencies",
                    "%s variants have a frequency of exactly 0 or 1, which makes the standard-error "
                    "formula divide by zero; policy eaf.degenerate is 'null', so the frequency was "
                    "blanked and the variants kept." % _fmt_int(n_degenerate),
                    before_degenerate,
                    final_df.height,
                    changed=n_degenerate,
                    warn=n_degenerate > 0,
                )
                qc_info["eaf_degenerate_removed"] = 0
            else:  # "keep"
                emit.qc(
                    "degenerate frequencies",
                    "%s variants have a frequency of exactly 0 or 1; policy eaf.degenerate is 'keep', so "
                    "they were left in place." % _fmt_int(n_degenerate),
                    final_df.height,
                    final_df.height,
                    warn=n_degenerate > 0,
                )
                qc_info["eaf_degenerate_removed"] = 0

            final_stats = summarise_allele_frequency(final_df, final_eaf_col, policies=pol)
            _add_prefixed_statistics(qc_info, "final_eaf", final_stats)

            sample_column_dict["eaf_col"] = final_eaf_col
            qc_info["final_eaf_col"] = final_eaf_col
            qc_info["final_variants"] = final_df.height
            qc_info["variants_removed"] = qc_info["Total_number_of_variants"] - final_df.height
            qc_info["final_status"] = "success"

            emit.info("Frequency column used          : %s" % final_eaf_col)
            emit.info("Initial variants               : %s" % _fmt_int(qc_info["Total_number_of_variants"]))
            emit.info("Variants with no frequency     : %s" % _fmt_int(missing_eaf_count))
            emit.info("Frequencies outside 0 to 1     : %s" % _fmt_int(out_of_bound_count))
            emit.info("Final variants                 : %s" % _fmt_int(qc_info["final_variants"]))
            emit.info("Variants removed               : %s" % _fmt_int(qc_info["variants_removed"]))
            if ctx is not None:
                ctx.set_rows(final_df.height)

        except Exception as exc:
            qc_info["final_status"] = "failed"
            qc_info["error_message"] = str(exc)
            raise

    return final_df, qc_info, sample_column_dict
