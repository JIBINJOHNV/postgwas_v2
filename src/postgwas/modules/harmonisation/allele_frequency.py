"""Step 04 (per chromosome) — reference orientation and EAF harmonisation.

Plan v3 sections 5.2 and 5.3.

What this module does, in order:

1.  Join the raw chromosome reference once on chromosome and position, resolve
    forward/swapped/reverse-complement orientations, and transform allele-specific
    statistics consistently.
2.  If the study supplies its own frequency column, describe it, apply the configured
    out-of-range action without silently clipping under other policies, and decide whether the
    column is really an *effect* allele frequency or a *minor* allele frequency.
    A MAF-like distribution is checked against the configured default EAF table:
    direct matches use ALT frequency and swapped matches use ``1-AF``. This follows
    the GWAS-SSF requirement that frequency is aligned to the listed effect allele
    (https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics).
3.  Otherwise take the frequency from an external panel, matching each variant in
    the two orientations named by ``eaf.match_orientations`` (direct and swap —
    never complement). Before swap alignment, screen the deduplicated raw
    frequency column for a MAF-like distribution. A non-MAF-like result retains
    the user's declared ALT/EAF contract. A MAF-like result must be confirmed
    against an independent default EAF panel; using the same file, confirmed
    MAF, or inconclusive evidence stops rather than guessing or auto-converting.
    Strand orientation has already been validated against reference REF/ALT. The
    external file's configured allele mapping determines whether its frequency
    is used as listed or as ``1-AF``; it never re-decides strand or changes study
    beta, Z, or ``strand_action``.
4.  Re-run the configured duplicate validation after REF/ALT and EAF alignment,
    catching opposite-strand representations that now identify the same
    physical variant.
5.  Apply policy ``eaf.degenerate`` to frequencies of exactly 0 or 1; the
    default removes them, and Z-only reconstruction always rejects any endpoint
    retained by an explicit alternative because its denominator would be zero.

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

import math
import os
from typing import Any, Dict, Optional, Tuple

import polars as pl

from postgwas.core.paths import configured_output_path

from postgwas.core.values import (
    format_count,
    format_fraction_percentage,
    optional_text,
)

from .shared.runtime import (
    active_context,
    reject_rows,
    resolve_policies,
    step_context as shared_step_context,
)
from .shared.allele_join import allele_oriented_left_join
from .shared.statistics import valid_frequency_mask
from .shared.variant_columns import (
    palindromic_snp_expression,
    read_reference_variant_table,
)
from .summary_statistics_io import resolve_duplicate_variants
from .strand import (
    POLICY_KEYS as STRAND_POLICY_KEYS,
    REFERENCE_AF_COLUMN,
    RESOLVED_STRAND_ACTIONS,
    STEP_LABEL as STRAND_STEP_LABEL,
    STRAND_ACTION_COLUMN,
    harmonise_strand_orientation,
)

__all__ = [
    "AlleleFrequencyError",
    "STEP_LABEL",
    "STEP_TITLE",
    "POLICY_KEYS",
    "summarise_allele_frequency",
    "validate_allele_frequency",
    "confirm_eaf_with_reference",
    "merge_external_allele_frequencies",
    "harmonise_allele_frequency",
]


# -------------------------------------------------------
# 0. Module constants and the typed error
# -------------------------------------------------------
STEP_LABEL = "04 eaf_harmonisation"
STEP_TITLE = "Allele orientation and effect allele frequency"
FUNC_NAME = "allele_frequency.harmonise_allele_frequency"

#: The settings block printed at the top of the step.
POLICY_KEYS = (
    "eaf.out_of_range",
    "eaf.clip_tolerance",
    "eaf.degenerate",
    "eaf.maf_decision_cutoff",
    "eaf.maf_reference_min_overlap",
    "eaf.maf_reference_correlation_method",
    "eaf.maf_reference_min_correlation",
    "eaf.maf_reference_max_mean_absolute_difference",
    "eaf.reference_minor_fraction_cutoff",
    "eaf.maf_reference_error_margin",
    "eaf.external_min_match_fraction",
    "external_reference.exact_duplicate_action",
    "external_reference.non_identical_duplicate_action",
    "eaf.match_orientations",
    "chromosome.strip_chr_prefix",
    "chromosome.strip_leading_zero",
    "chromosome.rename_map",
    "eaf.invalid_fraction_cutoff",
    "eaf.missing_fraction_cutoff",
    "duplicates.key",
    "duplicates.consistency_fields",
    "duplicates.post_orientation_relative_tolerance",
    "duplicates.quality_fields",
    "duplicates.selection_order",
    "duplicates.conflicting_action",
)

class AlleleFrequencyError(RuntimeError):
    """A fatal problem in the EAF step.

    Derives from ``RuntimeError`` so the dataset orchestrator can catch it,
    record the chromosome failure, and apply its retry policy.
    """


_fmt_int = format_count
_fmt_pct = format_fraction_percentage


def _fmt_metric(value: Any) -> str:
    """Render one optional finite QC metric without hiding missing evidence."""
    if value is None:
        return "not available"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    return "%.6g" % number if math.isfinite(number) else "not available"


def _maf_reference_suitability_text(evidence: Dict[str, Any]) -> str:
    """Summarise the YAML-controlled folded-MAF suitability evidence."""
    return (
        "folded-MAF %s correlation: %s (minimum %s); mean absolute "
        "folded-MAF difference: %s (maximum %s)"
        % (
            evidence.get("maf_correlation_method", "unknown"),
            _fmt_metric(evidence.get("maf_correlation")),
            _fmt_metric(evidence.get("minimum_maf_correlation")),
            _fmt_metric(evidence.get("mean_absolute_maf_difference")),
            _fmt_metric(evidence.get(
                "maximum_mean_absolute_maf_difference"
            )),
        )
    )


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

    ``clip_tolerance`` (policy ``eaf.clip_tolerance``) is retained as a
    diagnostic boundary for overshoots. It never makes a value above 1 valid.
    """
    if colname not in df.columns:
        return {
            "column_present": False,
            "total_rows": df.height,
        }
    tolerance = float(
        clip_tolerance
        if clip_tolerance is not None
        else resolve_policies(policies).get("eaf.clip_tolerance")
    )
    total_rows = df.height
    if total_rows == 0:
        return {
            "column_present": True,
            "total_rows": 0,
            "non_null_rows": 0,
            "null_rows": 0,
            "clip_tolerance": tolerance,
        }
    value = pl.col(colname).cast(pl.Float64, strict=False)
    usable = valid_frequency_mask(value)
    stats_row = df.select([
        pl.len().alias("total_rows"),
        value.is_not_null().sum().alias("non_null_rows"),
        value.is_null().sum().alias("null_rows"),
        usable.sum().alias("usable_rows"),
        (value < 0.0).sum().alias("lt_0_rows"),
        (value > 1.0).sum().alias("gt_1_rows"),
        (value > tolerance).sum().alias("gt_1_05_rows"),
        (
            value.is_not_null() & ~value.is_finite()
        ).sum().alias("non_finite_rows"),
        (value < 0.0).sum().alias("invalid_negative_rows"),
        (
            value.is_not_null()
            & (
                ~value.is_finite()
                | (value < 0.0)
                | (value > 1.0)
            )
        ).sum().alias("invalid_total_rows"),
        (usable & (value == 0.0)).sum().alias("eq_0_rows"),
        (usable & (value == 1.0)).sum().alias("eq_1_rows"),
        (usable & (value <= 0.5)).sum().alias("le_0_5_rows"),
        pl.when(usable).then(value).min().alias("min"),
        pl.when(usable).then(value).max().alias("max"),
        pl.when(usable).then(value).mean().alias("mean"),
        pl.when(usable).then(value).median().alias("median"),
        pl.when(usable).then(value).quantile(0.25).alias("q25"),
        pl.when(usable).then(value).quantile(0.75).alias("q75"),
    ]).to_dicts()[0]
    usable_rows = stats_row["usable_rows"] or 0
    total_rows = stats_row["total_rows"] or 0

    def frac(n, d):
        return float(n) / float(d) if d else 0.0

    stats_row["null_fraction"] = frac(stats_row["null_rows"], total_rows)
    stats_row["invalid_total_fraction"] = frac(stats_row["invalid_total_rows"], total_rows)
    stats_row["unusable_rows"] = total_rows - usable_rows
    stats_row["usable_fraction"] = frac(usable_rows, total_rows)
    stats_row["le_0_5_fraction_among_total"] = frac(stats_row["le_0_5_rows"], total_rows)
    stats_row["le_0_5_fraction_among_usable"] = frac(stats_row["le_0_5_rows"], usable_rows)
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

    ``clip``   explicitly force them back into 0..1
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

    ``is_valid_range`` is reported honestly: it is ``False`` when the share of
    non-finite values or values outside the closed interval 0..1
    exceeds ``invalid_fraction_cutoff``.  ``missing_count`` is the null count of
    the column — it used to be a row-count difference across ``with_columns``,
    which never changes the row count and so was always 0.
    """
    emit = active_context(logger, ctx, df.height)
    resolved = resolve_policies(policies)
    tolerance = float(
        clip_tolerance
        if clip_tolerance is not None
        else resolved.get("eaf.clip_tolerance")
    )
    action = str(
        out_of_range
        if out_of_range is not None
        else resolved.get("eaf.out_of_range")
    )

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
    outside_unit_count = out_of_range_count
    n_missing = int(initial_stats["null_rows"] or 0)

    out_of_range_fraction = out_of_range_count / total_before if total_before else 0.0
    is_valid_range = out_of_range_fraction <= invalid_fraction_cutoff

    if out_of_range_count and out_of_range_fraction > invalid_fraction_cutoff:
        emit.warn(
            "Column '%s' has %s non-finite values or values outside 0 to 1, which is %s of all %s variants - far more "
            "than the %s allowed by eaf.invalid_fraction_cutoff. The column may not be a frequency at all."
            % (
                colname,
                _fmt_int(out_of_range_count),
                _fmt_pct(out_of_range_fraction),
                _fmt_int(total_before),
                _fmt_pct(invalid_fraction_cutoff),
            )
        )

    # ---- apply the out-of-range policy ---------------------------------
    cleaned_df = df
    invalid_frequency = (
        pl.col(colname).is_not_null()
        & (
            ~pl.col(colname).is_finite()
            | (pl.col(colname) < 0.0)
            | (pl.col(colname) > 1.0)
        )
    ).fill_null(False)
    if out_of_range_count and action == "fail":
        raise AlleleFrequencyError(
            "Column '%s' has %s non-finite frequency values or values outside 0 to 1 and policy "
            "eaf.out_of_range is 'fail'." % (colname, _fmt_int(out_of_range_count))
        )
    if action == "reject":
        # fill_null(False): a null frequency is not out of range, it is missing,
        # and missingness is somebody else's decision
        cleaned_df, _ = reject_rows(
            cleaned_df,
            invalid_frequency,
            step_label=STEP_LABEL,
            reason="eaf_out_of_range",
            context=emit,
            collector=rejects,
            detail="%s non-finite or outside 0 to 1" % colname,
            check_name="frequency range",
            description=(
                "Non-finite frequencies and frequencies outside 0 to 1 cannot be probabilities; "
                "policy eaf.out_of_range is 'reject', so those variants were removed."
            ),
            warn_on_remove=True,
            warn_without_collector=True,
            record_empty=True,
        )
    elif action == "null":
        before = cleaned_df.height
        cleaned_df = cleaned_df.with_columns(
            pl.when(invalid_frequency)
            .then(None)
            .otherwise(pl.col(colname))
            .alias(colname)
        )
        emit.qc(
            "frequency range",
            "%s values are non-finite or outside 0 to 1; policy eaf.out_of_range is 'null', so those "
            "frequencies were blanked and the variants kept." % _fmt_int(out_of_range_count),
            before,
            cleaned_df.height,
            changed=out_of_range_count,
            warn=out_of_range_count > 0,
            step=STEP_LABEL,
        )
    elif action == "clip":  # explicit opt-in correction
        before = cleaned_df.height
        cleaned_df = cleaned_df.with_columns(
            pl.when(pl.col(colname).is_not_null() & pl.col(colname).is_finite())
            .then(pl.col(colname).clip(0.0, 1.0))
            .otherwise(None)
            .alias(colname)
        )
        emit.qc(
            "frequency range",
            "%s values are non-finite or outside 0 to 1; policy eaf.out_of_range is 'clip', so "
            "finite values were forced into range and non-finite values were blanked."
            % _fmt_int(outside_unit_count),
            before,
            before,
            changed=outside_unit_count,
            warn=outside_unit_count > 0,
            step=STEP_LABEL,
        )
    else:  # "fail" with no invalid values
        emit.qc(
            "frequency range",
            "No non-finite frequencies or values outside 0 to 1 were found; "
            "policy eaf.out_of_range is 'fail', so no action was needed.",
            cleaned_df.height,
            cleaned_df.height,
            changed=0,
            step=STEP_LABEL,
        )

    # ---- missingness (null_count, not a row-count difference) ----------
    final_stats = summarise_allele_frequency(cleaned_df, colname, clip_tolerance=tolerance)
    final_missing = int(final_stats.get("null_rows") or 0)
    missing_fraction = final_missing / total_before if total_before else 0.0
    if missing_fraction > missing_fraction_cutoff:
        emit.warn(
            "Column '%s' is empty for %s of the %s variants (%s), more than the %s allowed by "
            "eaf.missing_fraction_cutoff."
            % (
                colname,
                _fmt_int(final_missing),
                _fmt_int(total_before),
                _fmt_pct(missing_fraction),
                _fmt_pct(missing_fraction_cutoff),
            )
        )

    # ---- MAF-likeness, measured over usable finite probabilities --------
    total_after = cleaned_df.height
    low_freq = int(final_stats.get("le_0_5_rows") or 0)
    low_freq_pct = float(final_stats.get("le_0_5_fraction_among_usable") or 0.0)
    usable_rows = int(final_stats.get("usable_rows") or 0)
    is_suspicious = low_freq_pct > maf_cutoff
    if is_suspicious:
        emit.warn(
            "Column '%s': %s of its %s usable finite frequencies are at or below 0.5, more than the %s set by "
            "eaf.maf_decision_cutoff. That is what a minor allele frequency looks like, not an effect "
            "allele frequency."
            % (colname, _fmt_pct(low_freq_pct), _fmt_int(usable_rows), _fmt_pct(maf_cutoff))
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
        "initial_missing_count": n_missing,
        "missing_count_removed": final_missing,
        "missing_fraction_removed": missing_fraction,
        "missing_count": final_missing,
        "missing_fraction": missing_fraction,
        "low_freq_count_after_cleaning": low_freq,
        "low_freq_fraction_after_cleaning": low_freq_pct,
        "maf_like_flag": is_suspicious,
    }
    return is_valid_range, is_suspicious, out_of_range_count, final_missing, cleaned_df, low_freq_pct, combined_stats


# -------------------------------------------------------
# 7. Confirm a MAF-like frequency column against default EAF
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
    excluded because their strand cannot be inferred from the allele pair. Before
    testing direction, both frequencies are folded as ``min(AF, 1-AF)`` and the
    configured correlation plus mean-absolute-difference gates must show that
    the independent default panel meets the selected suitability criteria for
    these variants. Folding is a scientific invariant that removes allele
    direction; it cannot itself prove EAF. A reference set dominated by minor
    effect alleles is inconclusive because aligned EAF and folded MAF are then
    nearly identical; otherwise the two mean errors must differ by the
    configured minimum margin.
    """
    emit = active_context(logger, ctx, df.height)
    resolved = resolve_policies(policies)
    minimum = int(resolved.get("eaf.maf_reference_min_overlap"))
    correlation_method = str(
        resolved.get("eaf.maf_reference_correlation_method")
    )
    minimum_maf_correlation = float(
        resolved.get("eaf.maf_reference_min_correlation")
    )
    maximum_maf_difference = float(
        resolved.get("eaf.maf_reference_max_mean_absolute_difference")
    )
    reference_minor_cutoff = float(
        resolved.get("eaf.reference_minor_fraction_cutoff")
    )
    minimum_error_margin = float(
        resolved.get("eaf.maf_reference_error_margin")
    )
    keys = [std_cols[name] for name in ("chr", "pos", "ea", "oa")]
    palindromic = palindromic_snp_expression(
        pl.col(std_cols["ea"]), pl.col(std_cols["oa"])
    )
    use_aligned = aligned_reference_col in df.columns if aligned_reference_col else False
    selected = keys + [input_af_col]
    if use_aligned:
        selected.append(aligned_reference_col)
        if STRAND_ACTION_COLUMN in df.columns:
            selected.append(STRAND_ACTION_COLUMN)
    study = df.select(selected).filter(
        valid_frequency_mask(pl.col(input_af_col))
        & ~palindromic
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
            study_columns_canonical=study_columns_canonical,
        )
    matched = merged.filter(pl.col(reference_col).is_not_null())

    stats = {
        "direct_matches": merge_stats["direct_match_rows"],
        "swapped_matches": merge_stats["flip_match_rows"],
        "comparable_variants": matched.height,
        "minimum_overlap": minimum,
        "maf_correlation_method": correlation_method,
        "maf_correlation": None,
        "minimum_maf_correlation": minimum_maf_correlation,
        "mean_absolute_maf_difference": None,
        "maximum_mean_absolute_maf_difference": maximum_maf_difference,
        "reference_suitability_status": "not_evaluated",
        "reference_suitability_reason": None,
        "reference_effect_allele_minor_fraction": None,
        "reference_minor_fraction_cutoff": reference_minor_cutoff,
        "mean_eaf_error": None,
        "mean_maf_error": None,
        "absolute_error_difference": None,
        "minimum_error_margin": minimum_error_margin,
        "decision_reason": None,
    }
    reference_value = pl.col(reference_col).cast(pl.Float64, strict=False)
    invalid_reference = int(matched.select(
        (
            ~reference_value.is_finite()
            | (reference_value < 0.0)
            | (reference_value > 1.0)
        ).sum().alias("invalid_reference")
    ).item() or 0)
    if invalid_reference:
        raise AlleleFrequencyError(
            "The default EAF table '%s' contains %s matched non-finite values or values outside "
            "0 to 1 in column '%s'."
            % (reference_file, _fmt_int(invalid_reference), population_col)
        )
    if matched.height < minimum:
        stats["decision_reason"] = "insufficient_overlap"
        stats["reference_suitability_status"] = "inconclusive"
        stats["reference_suitability_reason"] = "insufficient_overlap"
        emit.warn(
            "The MAF-like column '%s' had only %s non-palindromic matches to the default EAF table; "
            "%s are required by eaf.maf_reference_min_overlap."
            % (input_af_col, _fmt_int(matched.height), _fmt_int(minimum))
        )
        return "inconclusive", stats

    study_value = pl.col("_study_af").cast(pl.Float64, strict=False)
    study_maf = pl.min_horizontal(study_value, 1.0 - study_value)
    reference_maf = pl.min_horizontal(
        reference_value, 1.0 - reference_value
    )
    comparison = matched.select([
        pl.corr(
            study_maf, reference_maf, method=correlation_method
        ).alias("maf_correlation"),
        (study_maf - reference_maf).abs().mean().alias("maf_difference"),
        (study_value - reference_value).abs().mean().alias("eaf"),
        (
            study_value - reference_maf
        ).abs().mean().alias("maf"),
        (reference_value <= 0.5).mean().alias("minor_fraction"),
    ]).to_dicts()[0]

    eaf_error = float(comparison["eaf"])
    maf_error = float(comparison["maf"])
    maf_difference = float(comparison["maf_difference"])
    correlation_value = comparison["maf_correlation"]
    maf_correlation = None
    if correlation_value is not None:
        candidate_correlation = float(correlation_value)
        if math.isfinite(candidate_correlation):
            maf_correlation = candidate_correlation
    error_difference = abs(eaf_error - maf_error)
    reference_minor_fraction = float(comparison["minor_fraction"])
    stats["maf_correlation"] = maf_correlation
    stats["mean_absolute_maf_difference"] = maf_difference
    stats["mean_eaf_error"] = eaf_error
    stats["mean_maf_error"] = maf_error
    stats["absolute_error_difference"] = error_difference
    stats["reference_effect_allele_minor_fraction"] = reference_minor_fraction

    if maf_correlation is None:
        suitability_reason = "maf_correlation_unavailable"
    elif maf_correlation < minimum_maf_correlation:
        suitability_reason = "maf_correlation_below_minimum"
    elif maf_difference > maximum_maf_difference:
        suitability_reason = "maf_difference_above_maximum"
    else:
        suitability_reason = "maf_correlation_and_difference_passed"
    reference_suitable = (
        suitability_reason == "maf_correlation_and_difference_passed"
    )
    stats["reference_suitability_status"] = (
        "suitable" if reference_suitable else "inconclusive"
    )
    stats["reference_suitability_reason"] = suitability_reason
    emit.decide(
        "Default EAF reference suitability",
        {
            "non-palindromic matches": _fmt_int(matched.height),
            "folded-MAF correlation method": correlation_method,
            "folded-MAF correlation": _fmt_metric(maf_correlation),
            "minimum correlation": minimum_maf_correlation,
            "mean absolute folded-MAF difference": maf_difference,
            "maximum mean absolute difference": maximum_maf_difference,
        },
        (
            "suitable for directional MAF/EAF confirmation"
            if reference_suitable
            else "inconclusive; the default panel is not proven suitable"
        ),
    )
    if not reference_suitable:
        stats["decision_reason"] = suitability_reason
        return "inconclusive", stats

    if reference_minor_fraction > reference_minor_cutoff:
        decision = "inconclusive"
        reason = "reference_effect_allele_mostly_minor"
    elif eaf_error + minimum_error_margin < maf_error:
        decision = "eaf"
        reason = "eaf_error_lower_by_required_margin"
    elif maf_error + minimum_error_margin < eaf_error:
        decision = "maf"
        reason = "maf_error_lower_by_required_margin"
    else:
        decision = "inconclusive"
        reason = "error_difference_below_required_margin"
    stats["decision_reason"] = reason
    emit.decide(
        "MAF-like frequency column '%s'" % input_af_col,
        {
            "non-palindromic reference matches": _fmt_int(matched.height),
            "mean error as EAF": eaf_error,
            "mean error as MAF": maf_error,
            "absolute error difference": error_difference,
            "required error margin": minimum_error_margin,
            "reference effect allele is minor": _fmt_pct(reference_minor_fraction),
            "uninformative-panel cutoff": _fmt_pct(reference_minor_cutoff),
        },
        {
            "eaf": "EAF aligned to the effect allele",
            "maf": "MAF not aligned to the effect allele",
            "inconclusive": "inconclusive; PostGWAS will not guess",
        }[decision],
    )
    return decision, stats


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
    min_match_fraction: Optional[float] = None,
    study_columns_canonical: bool = False,
    validation_reference_col: Optional[str] = None,
    validation_reference_file: Optional[str] = None,
    validation_reference_column: Optional[str] = None,
) -> Tuple[pl.DataFrame, str, Dict[str, Any]]:
    """Take the frequency from an external panel.

    Variants are matched in the orientations listed by ``eaf.match_orientations``
    (``direct`` then ``swap``); complement passes are deliberately not offered,
    see plan section 5.3. Unmatched variants keep a null frequency. The coverage
    gate counts only allele-key matches whose frequency is finite and between 0
    and 1; if that usable share falls below
    ``eaf.external_min_match_fraction`` the step fails instead of quietly
    reporting success. The raw external column is screened before swap-based
    alignment. When that screen is MAF-like and ``validation_reference_col``
    is available on the study frame, the existing EAF-versus-folded-MAF
    classifier validates the external column against that independent,
    allele-aligned reference before the caller adopts it.
    """
    resolved_policies = resolve_policies(policies)
    emit = active_context(logger, ctx, df.height)
    min_match_fraction = float(
        min_match_fraction
        if min_match_fraction is not None
        else resolved_policies.get("eaf.external_min_match_fraction")
    )
    duplicate_exact_action = str(
        resolved_policies.get("external_reference.exact_duplicate_action")
    )
    duplicate_non_identical_action = str(
        resolved_policies.get(
            "external_reference.non_identical_duplicate_action"
        )
    )
    orientations = [
        str(value).lower()
        for value in resolved_policies.get("eaf.match_orientations")
    ]
    unsupported = [o for o in orientations if o not in ("direct", "swap")]
    if unsupported:
        raise AlleleFrequencyError(
            "eaf.match_orientations asks for %s, but only 'direct' and 'swap' are implemented. "
            "Study rows have already received a validated per-variant strand_action; the external "
            "frequency join must not perform a second complement pass or re-decide effect direction."
            % unsupported
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
        "palindromic_evaluated_rows": 0,
        "palindromic_orientation_resolved_rows": 0,
        "palindromic_orientation_rejected_rows": 0,
        "palindromic_orientation_state_counts": {},
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
    ext_df, detected = read_reference_variant_table(
        load_path,
        colmap,
        policies,
        value_columns=[eaf_col],
        error_type=AlleleFrequencyError,
        description="external frequency file",
    )
    input_format = (
        detected.method
        if detected.kind == "unknown"
        else "%s separated" % detected.method
    )
    emit.info(
        "Read the external frequency file %s: %s rows, %d required columns, %s."
        % (
            os.path.basename(load_path),
            _fmt_int(ext_df.height),
            ext_df.width,
            input_format,
        )
    )
    stats["external_rows_loaded"] = ext_df.height
    row_col = "__postgwas_external_eaf_row"
    while row_col in df.columns:
        row_col += "_"
    study_for_join = df.with_row_index(row_col)
    matched_all, orientation_col, join_stats = allele_oriented_left_join(
        study_for_join,
        ext_df,
        study_columns=std_cols,
        reference_columns={
            "chr": colmap["chr"],
            "pos": colmap["pos"],
            "ea": colmap["a1"],
            "oa": colmap["a2"],
        },
        policies=resolved_policies,
        value_column=eaf_col,
        output_column=eaf_col,
        duplicate_exact_action=duplicate_exact_action,
        duplicate_non_identical_action=duplicate_non_identical_action,
        orientations=orientations,
        swapped_value="one_minus",
        prefer_non_null_value=True,
        study_columns_canonical=study_columns_canonical,
        error_type=AlleleFrequencyError,
        warn=emit.warn,
        reference_label="external frequency file '%s'" % load_path,
    )
    if join_stats["reference_duplicate_groups"]:
        emit.info(
            "External frequency reference duplicates: %s exact key/value groups "
            "(action '%s'); %s non-identical key groups covering %s rows "
            "(action '%s'); %s reference rows removed in total."
            % (
                _fmt_int(join_stats["reference_exact_duplicate_groups"]),
                duplicate_exact_action,
                _fmt_int(
                    join_stats["reference_non_identical_duplicate_groups"]
                ),
                _fmt_int(join_stats["reference_non_identical_duplicate_rows"]),
                duplicate_non_identical_action,
                _fmt_int(join_stats["reference_duplicate_rows_removed"]),
            )
        )
    stats.update(join_stats)

    raw_frequency_screen = dict(
        join_stats.get("raw_frequency_screen") or {}
    )
    stats["raw_frequency_screen"] = raw_frequency_screen
    stats["maf_reference_decision"] = "not_required"
    stats["maf_reference_validation"] = {}
    if (
        raw_frequency_screen.get("maf_like") is True
        and validation_reference_col is not None
    ):
        if validation_reference_col not in matched_all.columns:
            raise AlleleFrequencyError(
                "Internal PostGWAS execution-order error: external AF was "
                "MAF-like, but the aligned default reference column '%s' was "
                "not available for confirmation."
                % validation_reference_col
            )
        candidate_col = "__postgwas_external_raw_af"
        while candidate_col in matched_all.columns:
            candidate_col += "_"
        expected_col = "__postgwas_external_validation_eaf"
        while expected_col in matched_all.columns or expected_col == candidate_col:
            expected_col += "_"
        swapped = pl.col(orientation_col) == 1
        semantic_frame = matched_all.select([
            *[pl.col(std_cols[name]) for name in ("chr", "pos", "ea", "oa")],
            pl.when(swapped)
            .then(1.0 - pl.col(eaf_col).cast(pl.Float64, strict=False))
            .otherwise(pl.col(eaf_col).cast(pl.Float64, strict=False))
            .alias(candidate_col),
            pl.when(swapped)
            .then(
                1.0
                - pl.col(validation_reference_col).cast(
                    pl.Float64, strict=False,
                )
            )
            .otherwise(
                pl.col(validation_reference_col).cast(
                    pl.Float64, strict=False,
                )
            )
            .alias(expected_col),
        ])
        decision, validation = confirm_eaf_with_reference(
            semantic_frame,
            candidate_col,
            validation_reference_file or "aligned default EAF reference",
            validation_reference_column or validation_reference_col,
            colmap,
            std_cols,
            policies=resolved_policies,
            logger=logger,
            ctx=ctx,
            aligned_reference_col=expected_col,
            study_columns_canonical=True,
        )
        stats["maf_reference_decision"] = decision
        stats["maf_reference_validation"] = validation

    # The study row is already reference-oriented. The external mapping declares
    # which columns are REF and ALT/effect allele, so the shared direct/swap join
    # has already put its frequency on the final study EA (= reference ALT).
    # A palindrome is retained only when upstream strand_action proves that this
    # orientation was resolved; a reference-panel frequency must never be used
    # to infer the missing study strand.
    is_palindromic = palindromic_snp_expression(
        pl.col(std_cols["ea"]), pl.col(std_cols["oa"])
    )
    has_key_match = pl.col(orientation_col).is_not_null()
    matched_palindromic = is_palindromic & has_key_match
    orientation_resolved = (
        pl.col(STRAND_ACTION_COLUMN).is_in(list(RESOLVED_STRAND_ACTIONS))
        if STRAND_ACTION_COLUMN in matched_all.columns
        else pl.lit(False)
    ).fill_null(False)
    palindromic_orientation = (
        matched_all.filter(matched_palindromic)
        .select(row_col, orientation_resolved.alias("__orientation_resolved"))
        .with_columns(
            pl.when(pl.col("__orientation_resolved"))
            .then(pl.lit("resolved_by_study_orientation"))
            .otherwise(pl.lit("orientation_unavailable"))
            .alias("palindromic_orientation_state")
        )
        .drop("__orientation_resolved")
    )
    matched_all = matched_all.join(
        palindromic_orientation, on=row_col, how="left", coalesce=True
    ).with_columns(
        pl.when(
            matched_palindromic
            & (
                pl.col("palindromic_orientation_state")
                != "resolved_by_study_orientation"
            ).fill_null(True)
        )
        .then(None)
        .otherwise(pl.col(eaf_col))
        .alias(eaf_col)
    )
    palindromic_state_counts = {
        str(row["palindromic_orientation_state"]): int(row["len"])
        for row in palindromic_orientation.group_by(
            "palindromic_orientation_state"
        ).len().to_dicts()
    }
    stats["palindromic_evaluated_rows"] = palindromic_orientation.height
    stats["palindromic_orientation_resolved_rows"] = int(
        palindromic_state_counts.get("resolved_by_study_orientation", 0)
    )
    stats["palindromic_orientation_rejected_rows"] = (
        palindromic_orientation.height
        - stats["palindromic_orientation_resolved_rows"]
    )
    stats["palindromic_orientation_state_counts"] = palindromic_state_counts

    # An allele-key match is only evidence that the panel contains the variant;
    # it does not supply an EAF when the matched value is null, non-finite or
    # outside the probability interval. Count all value states in one Polars
    # aggregation so the configured coverage gate measures usable frequencies.
    reference_value = pl.col(eaf_col).cast(pl.Float64, strict=False)
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
    stats["unmatched_rows"] = join_stats["unmatched_rows"]
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
    final_df = matched_all
    if stats["palindromic_evaluated_rows"]:
        emit.info(
            "External EAF for palindromic variants: %s evaluated after study "
            "orientation, %s retained, and %s rejected because orientation "
            "evidence was unavailable."
            % (
                _fmt_int(stats["palindromic_evaluated_rows"]),
                _fmt_int(stats["palindromic_orientation_resolved_rows"]),
                _fmt_int(stats["palindromic_orientation_rejected_rows"]),
            )
        )
    final_df, _ = reject_rows(
        final_df,
        matched_palindromic
        & (
            pl.col("palindromic_orientation_state")
            != "resolved_by_study_orientation"
        ).fill_null(True),
        step_label=STEP_LABEL,
        reason="palindromic_orientation_unavailable",
        context=emit,
        collector=rejects,
        detail=(
            "external reference AF was available, but no resolved upstream "
            "strand_action proved the study palindrome's orientation"
        ),
        check_name="external-EAF palindromic orientation",
        description=(
            "A matched A/T or C/G external frequency is used only after automatic "
            "consensus/internal-study-EAF resolution or an explicit "
            "reference_aligned declaration has oriented the study row."
        ),
        warn_on_remove=True,
        warn_without_collector=True,
        record_empty=True,
    )

    final_df = final_df.drop([
        orientation_col,
        row_col,
        "palindromic_orientation_state",
    ])

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
    """Rows whose frequency is non-finite or outside 0..1 before correction."""
    if colname not in frame.columns:
        return frame.head(0)
    # Cast defensively: a frequency column read with schema inference disabled
    # arrives as Utf8, and comparing that against a float raises ComputeError.
    # Every sibling module casts with strict=False first; match them.
    value = pl.col(colname).cast(pl.Float64, strict=False)
    return frame.filter(
        (
            value.is_not_null()
            & (~value.is_finite() | (value < 0.0) | (value > 1.0))
        ).fill_null(False)
    )


def _resolve_reference_aligned_duplicates(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    output_layout: Dict[str, str],
    output_delimiter: str,
    policies,
    logger=None,
    ctx=None,
    rejects=None,
) -> Tuple[pl.DataFrame, Dict[str, Any]]:
    """Resolve rows that become identical only after reference orientation.

    Strand orientation has already made the study effect allele reference ALT
    and the other allele reference REF, and the final effect-frequency source
    has been aligned to ALT. Running the shared duplicate resolver at this
    boundary catches opposite-strand representations without applying a
    complement-based key to unnormalised indels. With the canonical key,
    distinct alternate alleles at one coordinate remain distinct because both
    REF and ALT participate in the key.
    """
    retained, report, statistics = resolve_duplicate_variants(
        df,
        sample_column_dict,
        policies=policies,
        logger=logger,
        rejects=rejects,
        step_label=STRAND_STEP_LABEL,
        reference_aligned=True,
    )
    report_path = configured_output_path(
        sample_column_dict["output_folder"],
        output_layout["post_orientation_duplicates"],
        error_type=AlleleFrequencyError,
        dataset_id=sample_column_dict["gwas_outputname"],
        chromosome=chromosome,
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report.write_csv(report_path, separator=output_delimiter)
    statistics["report"] = str(report_path)
    if ctx is not None:
        ctx.extra["post_orientation_duplicates"] = dict(statistics)
        if report.height:
            ctx.info(
                "Reference-aligned duplicate evidence for %s row%s was written "
                "to %s."
                % (
                    "{:,}".format(report.height),
                    "" if report.height == 1 else "s",
                    report_path,
                )
            )

    unsafe_groups = (
        int(statistics["conflicting_groups"])
        + int(statistics["swapped_orientation_groups"])
    )
    if (
        unsafe_groups
        and str(policies.get("duplicates.conflicting_action"))
        == "fail_dataset"
    ):
        raise AlleleFrequencyError(
            "Chromosome %s contains %s unsafe reference-aligned duplicate "
            "group%s. Their rows were recorded in %s and the rejected-variant "
            "output; processing stopped because duplicates.conflicting_action "
            "is 'fail_dataset'."
            % (
                chromosome,
                "{:,}".format(unsafe_groups),
                "" if unsafe_groups == 1 else "s",
                report_path,
            )
        )
    return retained, statistics


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
                   strand orientation and independent MAF/EAF confirmation.
    """
    gwas_outputname = sample_column_dict["gwas_outputname"]
    output_dir = sample_column_dict["output_folder"]
    merge_stats = {}

    pol = resolve_policies(policies)
    maf_cutoff = float(pol.get("eaf.maf_decision_cutoff"))
    invalid_fraction_cutoff = float(pol.get("eaf.invalid_fraction_cutoff"))
    missing_fraction_cutoff = float(pol.get("eaf.missing_fraction_cutoff"))
    degenerate_action = str(pol.get("eaf.degenerate"))

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
        "eaf_provenance": None,
        "final_eaf_col": None,
    }

    final_df = None
    with shared_step_context(
        logger,
        step_context,
        number=step_number,
        total=step_total,
        title=STEP_TITLE,
        operation=FUNC_NAME,
        rows_in=df.height,
        policy_keys=list(STRAND_POLICY_KEYS) + list(POLICY_KEYS),
    ) as ctx:
        emit = ctx
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
                qc_info["eaf_provenance"] = "study_supplied"

                # Capture the original invalid values before the configured
                # reject, null, clip or fail action changes the working frame.
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
                            minor_fraction = reference_stats.get(
                                "reference_effect_allele_minor_fraction"
                            )
                            minor_text = (
                                "not available"
                                if minor_fraction is None
                                else _fmt_pct(minor_fraction)
                            )
                            raise AlleleFrequencyError(
                                "Column '%s' is MAF-like, but the reference comparison was "
                                "inconclusive (%s). Comparable non-palindromic variants: %s "
                                "(minimum %s); %s; mean error as EAF: %s; mean error as MAF: %s; "
                                "required error separation: %s; reference-aligned EAF values at "
                                "or below 0.5: %s (uninformative above %s). PostGWAS will not "
                                "guess or apply the deferred allele-swap frequency flip. Provide "
                                "a frequency aligned to the listed effect allele, or configure a "
                                "dataset-supplied external EAF file."
                                % (
                                    internal_eaf,
                                    reference_stats.get("decision_reason", "unknown reason"),
                                    _fmt_int(reference_stats.get("comparable_variants", 0)),
                                    _fmt_int(reference_stats.get("minimum_overlap", 0)),
                                    _maf_reference_suitability_text(
                                        reference_stats
                                    ),
                                    reference_stats.get("mean_eaf_error"),
                                    reference_stats.get("mean_maf_error"),
                                    reference_stats.get("minimum_error_margin"),
                                    minor_text,
                                    _fmt_pct(reference_stats.get(
                                        "reference_minor_fraction_cutoff", 0.0
                                    )),
                                )
                            )
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
                    qc_info["final_status"] = "failed"
                    raise AlleleFrequencyError(
                        "Column '%s' has a frequency distribution consistent with a minor "
                        "allele frequency, but no study-wide MAF/EAF decision is available. "
                        "PostGWAS will not use an unoriented frequency column. Run the normal "
                        "dataset-level study-property step, or provide a dataset-supplied "
                        "external EAF file."
                        % internal_eaf
                    )
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

                mapping_keys = ("chr", "pos", "a1", "a2")
                target_delimiter = target_map.get("delimiter")
                default_delimiter = (
                    default_eaf_colmap.get("delimiter")
                    if default_eaf_colmap
                    else None
                )
                compatible_delimiter = (
                    target_delimiter == default_delimiter
                    or target_delimiter == "auto"
                    or default_delimiter == "auto"
                )
                same_reference_file = (
                    optional_text(default_eaf_file) is not None
                    and os.path.realpath(target_file)
                    == os.path.realpath(str(default_eaf_file))
                )
                same_reference_source = (
                    same_reference_file
                    and optional_text(default_eaf_column) == target_col
                    and default_eaf_colmap is not None
                    and compatible_delimiter
                    and all(
                        target_map.get(key) == default_eaf_colmap.get(key)
                        for key in mapping_keys
                    )
                )
                if same_reference_source:
                    if REFERENCE_AF_COLUMN not in df.columns:
                        raise AlleleFrequencyError(
                            "Internal PostGWAS execution-order error: the explicitly selected "
                            "external EAF source is identical to the strand reference, but its "
                            "validated aligned frequency column is unavailable."
                        )
                    raw_frequency_screen = strand_qc.get(
                        "raw_frequency_screen"
                    )
                    if not isinstance(raw_frequency_screen, dict):
                        raise AlleleFrequencyError(
                            "Internal PostGWAS execution-order error: the raw "
                            "external-frequency screen was not retained from "
                            "the identical strand reference."
                        )
                    reference_value = pl.col(REFERENCE_AF_COLUMN).cast(
                        pl.Float64, strict=False
                    )
                    finite_value = reference_value.is_finite().fill_null(False)
                    in_range = (
                        (reference_value >= 0.0) & (reference_value <= 1.0)
                    ).fill_null(False)
                    usable_value = finite_value & in_range
                    value_counts = df.select([
                        usable_value.sum().alias("usable"),
                        reference_value.is_null().sum().alias("missing"),
                        (
                            reference_value.is_not_null() & ~finite_value
                        ).sum().alias("non_finite"),
                        (
                            reference_value.is_not_null()
                            & finite_value
                            & ~in_range
                        ).sum().alias("out_of_range"),
                        palindromic_snp_expression(
                            pl.col(std_cols["ea"]), pl.col(std_cols["oa"])
                        ).sum().alias("palindromic"),
                    ]).to_dicts()[0]
                    usable_rows = int(value_counts["usable"] or 0)
                    usable_fraction = usable_rows / df.height if df.height else 0.0
                    minimum_fraction = float(
                        pol.get("eaf.external_min_match_fraction")
                    )
                    if df.height and usable_fraction < minimum_fraction:
                        raise AlleleFrequencyError(
                            "The external frequency file '%s' is identical to the validated "
                            "strand reference, but its aligned column '%s' supplied a usable EAF "
                            "for only %s of %s retained variants on chromosome %s (%s). This is "
                            "below the %s required by eaf.external_min_match_fraction."
                            % (
                                target_file,
                                target_col,
                                _fmt_int(usable_rows),
                                _fmt_int(df.height),
                                chromosome,
                                _fmt_pct(usable_fraction),
                                _fmt_pct(minimum_fraction),
                            )
                        )
                    actions = strand_qc.get("actions", {})
                    direct_rows = int(actions.get("forward", 0) or 0) + int(
                        actions.get("reverse_complement", 0) or 0
                    )
                    swapped_rows = int(
                        actions.get("forward_swapped", 0) or 0
                    ) + int(actions.get("reverse_complement_swapped", 0) or 0)
                    palindromic_rows = int(value_counts["palindromic"] or 0)
                    missing_rows = int(value_counts["missing"] or 0)
                    non_finite_rows = int(value_counts["non_finite"] or 0)
                    out_of_range_rows = int(value_counts["out_of_range"] or 0)
                    merge_stats = {
                        "external_rows_loaded": 0,
                        "input_rows_before_merge": df.height,
                        "direct_match_rows": direct_rows,
                        "flip_match_rows": swapped_rows,
                        "key_match_rows": df.height,
                        "key_match_fraction": 1.0 if df.height else 0.0,
                        "usable_match_rows": usable_rows,
                        "usable_match_fraction": usable_fraction,
                        "match_fraction": usable_fraction,
                        "matched_missing_frequency_rows": missing_rows,
                        "matched_non_finite_frequency_rows": non_finite_rows,
                        "matched_out_of_range_frequency_rows": out_of_range_rows,
                        "matched_unusable_frequency_rows": (
                            missing_rows + non_finite_rows + out_of_range_rows
                        ),
                        "unusable_frequency_rows": df.height - usable_rows,
                        "unmatched_rows": 0,
                        "total_missing": missing_rows,
                        "missing_fraction": (
                            missing_rows / df.height if df.height else 0.0
                        ),
                        "min_match_fraction": minimum_fraction,
                        "palindromic_evaluated_rows": palindromic_rows,
                        "palindromic_orientation_resolved_rows": palindromic_rows,
                        "palindromic_orientation_rejected_rows": 0,
                        "palindromic_orientation_state_counts": (
                            {"resolved_by_study_orientation": palindromic_rows}
                            if palindromic_rows
                            else {}
                        ),
                        "reused_aligned_strand_reference": True,
                        "raw_frequency_screen": raw_frequency_screen,
                        "maf_reference_decision": "not_required",
                        "maf_reference_validation": {},
                    }
                    merged_df = df
                    new_col_name = REFERENCE_AF_COLUMN
                    emit.info(
                        "Reused the allele-aligned population AF already retained from "
                        "strand orientation because the explicitly selected external EAF "
                        "file, column and allele mapping are identical; the reference "
                        "table was not read a second time."
                    )
                else:
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
                        study_columns_canonical=True,
                        validation_reference_col=(
                            None
                            if same_reference_file
                            else REFERENCE_AF_COLUMN
                        ),
                        validation_reference_file=default_eaf_file,
                        validation_reference_column=default_eaf_column,
                    )

                raw_frequency_screen = dict(
                    merge_stats.get("raw_frequency_screen") or {}
                )
                _add_prefixed_statistics(
                    qc_info,
                    "external_raw_af",
                    raw_frequency_screen,
                )
                if raw_frequency_screen:
                    emit.decide(
                        "Raw external frequency column '%s'" % target_col,
                        {
                            "usable frequencies": _fmt_int(
                                raw_frequency_screen.get("usable", 0)
                            ),
                            "at or below 0.5": "%s (%s)"
                            % (
                                _fmt_int(raw_frequency_screen.get(
                                    "at_or_below_0_5", 0
                                )),
                                _fmt_pct(raw_frequency_screen.get(
                                    "low_fraction_of_usable", 0.0
                                )),
                            ),
                            "decision cutoff": raw_frequency_screen.get(
                                "decision_cutoff"
                            ),
                        },
                        (
                            "MAF-like and requiring independent reference confirmation"
                            if raw_frequency_screen.get("maf_like") is True
                            else (
                                "not MAF-like; retaining the declared external "
                                "ALT/EAF contract"
                            )
                        ),
                    )
                if raw_frequency_screen.get("maf_like") is True:
                    if same_reference_file:
                        raise AlleleFrequencyError(
                            "The external AF column '%s' is MAF-like: %s of %s "
                            "usable values are at or below 0.5 (%s), exceeding "
                            "eaf.maf_decision_cutoff=%s. external_eaf_file and "
                            "the resource selected by "
                            "modules.harmonisation.default_eaf.source resolve "
                            "to the same file '%s', so PostGWAS cannot independently "
                            "determine whether the column is ALT/EAF or MAF. "
                            "Configure a different modules.harmonisation.default_eaf.source "
                            "and rerun."
                            % (
                                target_col,
                                _fmt_int(raw_frequency_screen.get(
                                    "at_or_below_0_5", 0
                                )),
                                _fmt_int(raw_frequency_screen.get("usable", 0)),
                                _fmt_pct(raw_frequency_screen.get(
                                    "low_fraction_of_usable", 0.0
                                )),
                                raw_frequency_screen.get("decision_cutoff"),
                                target_file,
                            )
                        )
                    external_decision = merge_stats.get(
                        "maf_reference_decision"
                    )
                    external_validation = dict(
                        merge_stats.get("maf_reference_validation") or {}
                    )
                    _add_prefixed_statistics(
                        qc_info,
                        "external_maf_reference",
                        external_validation,
                    )
                    qc_info["external_maf_reference_decision"] = (
                        external_decision
                    )
                    if external_decision == "maf":
                        raise AlleleFrequencyError(
                            "The external AF column '%s' was confirmed as minor "
                            "allele frequency rather than ALT/effect-allele "
                            "frequency. Configure an external_eaf_column containing "
                            "the frequency of the file's declared ALT allele."
                            % target_col
                        )
                    if external_decision != "eaf":
                        raise AlleleFrequencyError(
                            "The external AF column '%s' is MAF-like, but comparison "
                            "with the independent default EAF reference was "
                            "inconclusive (%s). Comparable non-palindromic variants: "
                            "%s (minimum %s); %s; mean error as EAF: %s; mean error as "
                            "MAF: %s; required error separation: %s. PostGWAS will "
                            "not guess. Configure a different external ALT/EAF "
                            "column or a more suitable independent "
                            "modules.harmonisation.default_eaf reference."
                            % (
                                target_col,
                                external_validation.get(
                                    "decision_reason", "unknown reason"
                                ),
                                _fmt_int(external_validation.get(
                                    "comparable_variants", 0
                                )),
                                _fmt_int(external_validation.get(
                                    "minimum_overlap", 0
                                )),
                                _maf_reference_suitability_text(
                                    external_validation
                                ),
                                external_validation.get("mean_eaf_error"),
                                external_validation.get("mean_maf_error"),
                                external_validation.get("minimum_error_margin"),
                            )
                        )

                _add_prefixed_statistics(qc_info, "external_merge", merge_stats)

                # Capture invalid reference values before the configured
                # frequency-range action changes the working frame.
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
                qc_info["eaf_provenance"] = "reference_imputed"

            # ------------------------------
            # FINALIZE
            # ------------------------------
            if final_df is None or final_eaf_col is None:
                raise AlleleFrequencyError("Failed to determine which column holds the effect allele frequency.")

            # The reference has now fixed ALT/REF orientation, and any deferred
            # MAF-like-column proof has also applied the required 1-EAF change
            # to swapped rows. Comparing duplicates before this boundary could
            # falsely call an EAF disagreement that step 04 was still required
            # to reconcile.
            sample_column_dict["eaf_col"] = final_eaf_col
            final_df, post_orientation_duplicate_qc = (
                _resolve_reference_aligned_duplicates(
                    chromosome=chromosome,
                    df=final_df,
                    sample_column_dict=sample_column_dict,
                    output_layout=output_layout,
                    output_delimiter=output_delimiter,
                    policies=pol,
                    logger=logger,
                    ctx=ctx,
                    rejects=rejects,
                )
            )
            qc_info["post_orientation_duplicates"] = (
                post_orientation_duplicate_qc
            )

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
                non_palindromic_action = str(
                    pol.get("strand.af_discordance_action")
                )
                palindromic_action = str(
                    pol.get("strand.palindromic_af_discordance_action")
                )
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
                palindromic = palindromic_snp_expression(
                    pl.col(std_cols["ea"]), pl.col(std_cols["oa"])
                )
                discordant = (
                    comparable & (pl.col(af_difference_col) > tolerance)
                ).fill_null(False)
                palindromic_discordant = discordant & palindromic
                non_palindromic_discordant = discordant & ~palindromic
                counts = final_df.select([
                    comparable.sum().alias("comparable"),
                    discordant.sum().alias("discordant"),
                    (comparable & palindromic).sum().alias("pal_comparable"),
                    palindromic_discordant.sum().alias("pal_discordant"),
                    (comparable & ~palindromic).sum().alias("non_pal_comparable"),
                    non_palindromic_discordant.sum().alias("non_pal_discordant"),
                ]).to_dicts()[0]
                n_comparable = int(counts["comparable"] or 0)
                n_discordant = int(counts["discordant"] or 0)
                n_pal_comparable = int(counts["pal_comparable"] or 0)
                n_pal_discordant = int(counts["pal_discordant"] or 0)
                n_non_pal_comparable = int(counts["non_pal_comparable"] or 0)
                n_non_pal_discordant = int(
                    counts["non_pal_discordant"] or 0
                )
                qc_info["strand_af_comparable"] = n_comparable
                qc_info["strand_af_discordant"] = n_discordant
                qc_info["strand_af_tolerance"] = tolerance
                qc_info["strand_af_discordance_action"] = non_palindromic_action
                qc_info["strand_palindromic_af_comparable"] = n_pal_comparable
                qc_info["strand_palindromic_af_discordant"] = n_pal_discordant
                qc_info["strand_palindromic_af_discordance_action"] = (
                    palindromic_action
                )
                qc_info["strand_non_palindromic_af_comparable"] = (
                    n_non_pal_comparable
                )
                qc_info["strand_non_palindromic_af_discordant"] = (
                    n_non_pal_discordant
                )
                if n_pal_discordant and palindromic_action == "fail":
                    raise AlleleFrequencyError(
                        "%s palindromic variants differ from the aligned population "
                        "reference AF by more than %s after resolved study orientation, "
                        "and strand.palindromic_af_discordance_action is 'fail'."
                        % (_fmt_int(n_pal_discordant), tolerance)
                    )
                if n_non_pal_discordant and non_palindromic_action == "fail":
                    raise AlleleFrequencyError(
                        "%s non-palindromic variants differ from the aligned population "
                        "reference AF by more than %s and "
                        "strand.af_discordance_action is 'fail'."
                        % (_fmt_int(n_non_pal_discordant), tolerance)
                    )
                if palindromic_action == "reject":
                    final_df, _ = reject_rows(
                        final_df,
                        palindromic_discordant,
                        step_label=STEP_LABEL,
                        reason="af_discordant",
                        context=emit,
                        collector=rejects,
                        detail=(
                            "palindromic variant has absolute EAF difference above %s "
                            "after resolved study orientation" % tolerance
                        ),
                        check_name="palindromic study/reference frequency concordance",
                        description=(
                            "After automatic consensus/internal study EAF or an explicit "
                            "reference-aligned declaration selected the orientation, "
                            "the palindromic variant's aligned study EAF differs from the "
                            "configured population reference by more than %s." % tolerance
                        ),
                        warn_on_remove=True,
                        warn_without_collector=True,
                        record_empty=True,
                    )
                else:
                    emit.qc(
                        "palindromic study/reference frequency concordance",
                        "%s of %s comparable variants differ from the aligned population "
                        "reference AF by more than %s; policy action is '%s'."
                        % (
                            _fmt_int(n_pal_discordant), _fmt_int(n_pal_comparable),
                            tolerance, palindromic_action,
                        ),
                        final_df.height,
                        final_df.height,
                        matched=n_pal_discordant,
                        outcome="retained_without_frequency_change",
                        warn=n_pal_discordant > 0,
                        step=STEP_LABEL,
                    )
                if non_palindromic_action == "reject":
                    final_df, _ = reject_rows(
                        final_df,
                        non_palindromic_discordant,
                        step_label=STEP_LABEL,
                        reason="af_discordant",
                        context=emit,
                        collector=rejects,
                        detail="absolute EAF difference above %s" % tolerance,
                        check_name="non-palindromic study/reference frequency concordance",
                        description=(
                            "The aligned study EAF differs from the configured population "
                            "reference by more than %s." % tolerance
                        ),
                        warn_on_remove=True,
                        warn_without_collector=True,
                        record_empty=True,
                    )
                else:
                    emit.qc(
                        "non-palindromic study/reference frequency concordance",
                        "%s of %s comparable variants differ from the aligned population "
                        "reference AF by more than %s; policy action is '%s'."
                        % (
                            _fmt_int(n_non_pal_discordant),
                            _fmt_int(n_non_pal_comparable), tolerance,
                            non_palindromic_action,
                        ),
                        final_df.height,
                        final_df.height,
                        matched=n_non_pal_discordant,
                        outcome="retained_without_frequency_change",
                        warn=n_non_pal_discordant > 0,
                        step=STEP_LABEL,
                    )
            else:
                qc_info["strand_af_comparable"] = 0
                qc_info["strand_af_discordant"] = 0
                qc_info["strand_palindromic_af_comparable"] = 0
                qc_info["strand_palindromic_af_discordant"] = 0
                qc_info["strand_non_palindromic_af_comparable"] = 0
                qc_info["strand_non_palindromic_af_discordant"] = 0

            # ---- frequencies of exactly 0 or 1 ----------------------------
            # These make either configured step-06 reconstruction denominator
            # exactly zero, which yields inf rather than null without its
            # reconstruction-specific guard.
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
                final_df, _ = reject_rows(
                    final_df,
                    degenerate_expr,
                    step_label=STEP_LABEL,
                    reason="eaf_degenerate",
                    context=emit,
                    collector=rejects,
                    detail="%s is exactly 0 or 1" % final_eaf_col,
                    check_name="degenerate frequencies",
                    description=(
                        "A frequency of exactly 0 or 1 makes the standard-error formula divide by zero, "
                        "so those variants were removed."
                    ),
                    warn_on_remove=True,
                    warn_without_collector=True,
                    record_empty=True,
                )
                qc_info["eaf_degenerate_removed"] = before_degenerate - final_df.height
            elif degenerate_action == "null":
                before_degenerate = final_df.height
                final_df = final_df.with_columns(
                    pl.when(degenerate_expr)
                    .then(None)
                    .otherwise(pl.col(final_eaf_col))
                    .alias(final_eaf_col)
                )
                emit.qc(
                    "degenerate frequencies",
                    "%s variants have a frequency of exactly 0 or 1, which makes the standard-error "
                    "formula divide by zero; policy eaf.degenerate is 'null', so the frequency was "
                    "blanked and the variants kept." % _fmt_int(n_degenerate),
                    before_degenerate,
                    final_df.height,
                    changed=n_degenerate,
                    warn=n_degenerate > 0,
                    step=STEP_LABEL,
                )
                qc_info["eaf_degenerate_removed"] = 0
            elif degenerate_action == "keep":
                emit.qc(
                    "degenerate frequencies",
                    "%s variants have a frequency of exactly 0 or 1; policy eaf.degenerate is 'keep', so "
                    "they were left in place." % _fmt_int(n_degenerate),
                    final_df.height,
                    final_df.height,
                    warn=n_degenerate > 0,
                    step=STEP_LABEL,
                )
                qc_info["eaf_degenerate_removed"] = 0
            else:  # "fail" with no degenerate values
                emit.qc(
                    "degenerate frequencies",
                    "No frequencies of exactly 0 or 1 were found; policy "
                    "eaf.degenerate is 'fail', so no action was needed.",
                    final_df.height,
                    final_df.height,
                    changed=0,
                    step=STEP_LABEL,
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
