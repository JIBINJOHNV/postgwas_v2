"""Step 12 - imputation quality (INFO) harmonisation.

Implements plan v3 section 5.6.  The INFO score is taken from the study's own
column when it has one, otherwise from a per-chromosome or whole-genome
reference file, and is then brought inside the configured window in exactly the
way the effect-allele-frequency step does.

There is deliberately no INFO equivalent of ``eaf_degenerate``: an INFO of 0 is
scientifically poor but not mathematically fatal. Harmonisation applies its
configured range action here, and the QC module independently assesses the
retained value against ``modules.qc_summary.rules`` after VCF creation.

All messages go through the PipelineLogger supplied by the orchestrator; the
module never writes to stdout.
"""

import os
from typing import Dict, Tuple

import polars as pl

from postgwas.core.dataframes import (
    chromosome_expression,
    count_non_null,
    position_expression,
    validate_cast_retention,
)
from postgwas.core.io.tables import read_delimited_table
from postgwas.core.values import optional_text

from .shared.runtime import resolve_policies, step_context
from .shared.allele_join import allele_oriented_left_join
from .shared.variant_columns import has_canonical_variant_columns

__all__ = ["harmonise_imputation_quality", "ImputationQualityError", "POLICY_KEYS", "STEP_LABEL"]


# The free-text step label recorded in the reject file.
STEP_LABEL = "11 info_harmonisation"

# Every policy this step reads.  Passed to logger.step() so the log shows the
# value and the plain-English explanation of each one.
POLICY_KEYS = [
    "info.source",
    "info.out_of_range",
    "info.clip_min",
    "info.clip_max",
    "info.clip_tolerance",
    "info.low_quality_threshold",
    "info.on_missing",
    "external_reference.exact_duplicate_action",
    "external_reference.non_identical_duplicate_action",
]

_REASON_TEXT = {
    "info_out_of_range": (
        "The imputation quality is outside the configured window and policy "
        "'info.out_of_range' is 'reject'."
    ),
    "info_missing": (
        "The variant has no imputation quality and policy 'info.on_missing' is "
        "'reject'."
    ),
}


class ImputationQualityError(RuntimeError):
    """INFO harmonisation cannot continue.

    A typed exception rather than sys.exit(), so the dataset-level retry can
    catch it.
    """


def harmonise_imputation_quality(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    info_file: str | None = None,
    info_column: str = "info",
    external_info_colmap=None,
    logger=None,
    policies=None,
    rejects=None,
    ctx=None,
    step_number: int = 11,
    step_total: int = 16,
) -> Tuple[pl.DataFrame, Dict, dict]:
    """Attach an imputation quality score to every variant and bring it in range.

    logger   : a PipelineLogger.  When given and ``ctx`` is not, this function
               opens its own step on it.
    policies : a Policies object.  When None the registry defaults are used, so
               behaviour is unchanged.
    rejects  : a RejectCollector.  Required for the 'reject' actions to record
               the removed variants.
    ctx      : a StepContext, when the caller opened the step itself.
    """
    policies = resolve_policies(policies)
    out_of_range = policies.get("info.out_of_range")
    clip_min = float(policies.get("info.clip_min"))
    clip_max = float(policies.get("info.clip_max"))
    clip_tolerance = float(policies.get("info.clip_tolerance"))
    # The tolerance is expressed as an absolute value ("values up to 1.05 are
    # rounding error"), so it can never sit below the ceiling itself.
    tolerance_limit = max(clip_tolerance, clip_max)
    low_quality = float(policies.get("info.low_quality_threshold"))
    on_missing = policies.get("info.on_missing")
    duplicate_exact_action = str(
        policies.get("external_reference.exact_duplicate_action")
    )
    duplicate_non_identical_action = str(
        policies.get("external_reference.non_identical_duplicate_action")
    )
    source = policies.get("info.source")

    qc_info = {"initial_variants": df.height}
    rows_in = df.height

    with step_context(
        logger, ctx,
        number=step_number, total=step_total,
        title="Imputation quality (INFO)",
        operation="imputation_quality.harmonise_imputation_quality",
        rows_in=rows_in, policy_keys=POLICY_KEYS,
    ) as step:
        chr_col = sample_column_dict["chr_col"]
        pos_col = sample_column_dict["pos_col"]
        ea_col = sample_column_dict["ea_col"]
        oa_col = sample_column_dict["oa_col"]
        imp_info_col = optional_text(sample_column_dict.get("imp_info_col"))
        fixed_info = sample_column_dict.get("fixed_info")
        fixed_info_column = optional_text(
            sample_column_dict.get("fixed_info_column")
        )
        resolved_source = optional_text(sample_column_dict.get("info_source"))
        if resolved_source not in (None, "internal", "external", "fixed_cli"):
            raise ImputationQualityError(
                "Unknown resolved INFO source %r on chromosome %s."
                % (resolved_source, chromosome)
            )

        # ------------------------------------------------------------------
        # NORMALISE THE STUDY COORDINATES
        # ------------------------------------------------------------------
        if not has_canonical_variant_columns(df, sample_column_dict):
            schema = dict(df.schema)
            pos_before = count_non_null(df, pos_col)
            df = df.with_columns([
                chromosome_expression(chr_col, schema[chr_col]),
                position_expression(pos_col, schema[pos_col]),
                pl.col(ea_col).cast(pl.String).str.to_uppercase().str.strip_chars(),
                pl.col(oa_col).cast(pl.String).str.to_uppercase().str.strip_chars(),
            ])
            validate_cast_retention(
                pos_before,
                count_non_null(df, pos_col),
                "study position",
                error_type=ImputationQualityError,
                warn=step.warn,
            )

        has_study_column = imp_info_col is not None and imp_info_col in df.columns
        require_study = (
            resolved_source == "internal"
            or (resolved_source is None and source == "study")
        )
        if require_study and not has_study_column:
            raise ImputationQualityError(
                "The resolved INFO source requires study imputation-quality "
                "column %r, but that column is not present on chromosome %s."
                % (imp_info_col, chromosome)
            )
        use_study_column = has_study_column and (
            resolved_source == "internal"
            or (resolved_source is None and source in ("auto", "study"))
        )

        # ------------------------------------------------------------------
        # STEP 1: THE STUDY'S OWN INFO COLUMN
        # ------------------------------------------------------------------
        external_file = optional_text(info_file)
        external_column = optional_text(info_column)
        has_external_source = (
            external_file is not None
            and external_column is not None
            and os.path.exists(external_file)
        )

        if use_study_column:
            step.info(
                "The study supplies its own imputation quality column %r, so no "
                "reference file was needed." % imp_info_col
            )
            df = df.with_columns(pl.col(imp_info_col).cast(pl.Float64, strict=False))
            info_col = imp_info_col
            qc_info["source"] = "study_column:%s" % imp_info_col
        elif has_external_source and (
            resolved_source == "external"
            or (
                resolved_source is None
                and source in ("auto", "reference", "external")
            )
        ):
            info_file_to_use = external_file
            info_col_to_use = external_column

            step.info(
                "Imputation quality comes from %s, column %r."
                % (info_file_to_use, info_col_to_use)
            )

            # --------------------------------------------------------------
            # STEP 2: LOAD THE REFERENCE
            # --------------------------------------------------------------
            if not isinstance(external_info_colmap, dict):
                raise ImputationQualityError(
                    "External imputation-quality column mapping is not configured."
                )
            required_mapping = ("chr", "pos", "a1", "a2", "delimiter")
            missing_mapping = [
                key for key in required_mapping
                if not str(external_info_colmap.get(key, "")).strip()
            ]
            if missing_mapping:
                raise ImputationQualityError(
                    "External imputation-quality mapping is missing: %s"
                    % ", ".join(missing_mapping)
                )
            ref_chr = external_info_colmap["chr"]
            ref_pos = external_info_colmap["pos"]
            ref_effect = external_info_colmap["a1"]
            ref_other = external_info_colmap["a2"]
            info_raw, _detected = read_delimited_table(
                info_file_to_use,
                external_info_colmap["delimiter"],
                candidates=list(policies.get("input.delimiter_candidates")),
                minimum_columns=int(policies.get("input.delimiter_min_columns")),
                maximum_columns=int(policies.get("input.delimiter_max_columns")),
                sample_lines=int(policies.get("input.delimiter_sample_rows")),
                null_values=list(policies.get("input.null_values")),
                infer_schema_length=int(policies.get("input.schema_inference_rows")),
                error_type=ImputationQualityError,
                description="external imputation-quality table",
            )
            if info_col_to_use in df.columns:
                step.warn(
                    "The study frame already has a column named %r; it is being "
                    "replaced by the imputation quality taken from the reference."
                    % info_col_to_use
                )
            df, orientation_col, join_stats = allele_oriented_left_join(
                df,
                info_raw,
                study_columns={
                    "chr": chr_col, "pos": pos_col, "ea": ea_col, "oa": oa_col,
                },
                reference_columns={
                    "chr": ref_chr, "pos": ref_pos,
                    "ea": ref_effect, "oa": ref_other,
                },
                value_column=info_col_to_use,
                output_column=info_col_to_use,
                duplicate_exact_action=duplicate_exact_action,
                duplicate_non_identical_action=duplicate_non_identical_action,
                swapped_value="same",
                prefer_non_null_value=True,
                study_columns_canonical=True,
                error_type=ImputationQualityError,
                warn=step.warn,
                reference_label="external imputation-quality table",
            )
            df = df.drop(orientation_col)
            reference_before = join_stats["reference_rows_before_deduplication"]
            reference_after = join_stats["reference_rows_after_deduplication"]
            step.info("Imputation quality reference rows: {:,}".format(reference_before))
            if join_stats["reference_duplicate_groups"]:
                step.info(
                    "External INFO reference duplicates: {:,} exact key/value "
                    "groups (action {!r}); {:,} non-identical key groups covering "
                    "{:,} rows (action {!r}); {:,} reference rows removed "
                    "({:,} before, {:,} after)."
                    .format(
                        join_stats["reference_exact_duplicate_groups"],
                        duplicate_exact_action,
                        join_stats[
                            "reference_non_identical_duplicate_groups"
                        ],
                        join_stats["reference_non_identical_duplicate_rows"],
                        duplicate_non_identical_action,
                        join_stats["reference_duplicate_rows_removed"],
                        reference_before,
                        reference_after,
                    )
                )
            matched_direct = join_stats["direct_value_matches"]
            matched_flipped = join_stats["swapped_value_matches"]

            step.info(
                "Matched {:,} variants with the alleles as given and {:,} more "
                "with the alleles the other way round, out of {:,}.".format(
                    matched_direct, matched_flipped, df.height
                )
            )
            qc_info["matched_direct"] = matched_direct
            qc_info["matched_flipped"] = matched_flipped
            qc_info.update(join_stats)

            info_col = info_col_to_use
            sample_column_dict["imp_info_col"] = info_col
            qc_info["source"] = "external_file:%s" % info_file_to_use
        elif fixed_info is not None and (
            resolved_source == "fixed_cli"
            or (resolved_source is None and source == "auto")
        ):
            if fixed_info_column is None:
                raise ImputationQualityError(
                    "The fixed INFO working column is not configured."
                )
            if fixed_info_column in df.columns:
                raise ImputationQualityError(
                    "Cannot assign fixed INFO on chromosome %s because the reserved "
                    "working column %r already exists in the study data. Rename that "
                    "unconfigured input column before retrying."
                    % (chromosome, fixed_info_column)
                )
            fixed_value = float(fixed_info)
            if not 0 <= fixed_value <= 1:
                raise ImputationQualityError(
                    "Fixed INFO must be between 0 and 1; received %r."
                    % fixed_info
                )
            df = df.with_columns(
                pl.lit(fixed_value, dtype=pl.Float64).alias(fixed_info_column)
            )
            info_col = fixed_info_column
            sample_column_dict["imp_info_col"] = info_col
            qc_info["source"] = "fixed_cli:%g" % fixed_value
            qc_info["fixed_info"] = fixed_value
            qc_info["fixed_info_user_assigned"] = True
            step.info(
                "Assigned INFO={:g} to all {:,} variants from the explicit "
                "--fixed-info option. This constant was supplied by the user; "
                "it was not measured per variant.".format(fixed_value, df.height)
            )
        else:
            requested = " Policy 'info.source' is %r." % source
            raise ImputationQualityError(
                "No usable imputation-quality source remains for chromosome %s.%s "
                "Provide imputation_info_column, provide external_info_file plus "
                "external_info_column, or explicitly supply --fixed-info VALUE "
                "with the default source resolution."
                % (chromosome, requested)
            )

        qc_info["info_column"] = info_col

        # ------------------------------------------------------------------
        # QC METRICS - measured before anything is changed
        # ------------------------------------------------------------------
        total = df.height
        stats = df.select([
            pl.col(info_col).is_null().sum().alias("missing"),
            (pl.col(info_col) < clip_min).sum().alias("lt_zero"),
            (pl.col(info_col) > clip_max).sum().alias("gt_one"),
            (pl.col(info_col) > tolerance_limit).sum().alias("gt_tolerance"),
            (
                (pl.col(info_col) > clip_max)
                & (pl.col(info_col) <= tolerance_limit)
            ).sum().alias("within_tolerance"),
            (
                pl.col(info_col).is_not_null()
                & (pl.col(info_col) <= low_quality)
            ).sum().alias("low_info"),
            (
                (pl.col(info_col) > low_quality)
                & (pl.col(info_col) <= clip_max)
            ).sum().alias("good_info"),
            pl.col(info_col).min().alias("min"),
            pl.col(info_col).max().alias("max"),
            pl.col(info_col).mean().alias("mean"),
            pl.col(info_col).median().alias("median"),
        ]).to_dicts()[0]
        missing = int(stats["missing"] or 0)
        lt_zero = int(stats["lt_zero"] or 0)
        gt_one = int(stats["gt_one"] or 0)
        gt_1_05 = int(stats["gt_tolerance"] or 0)
        within_tolerance = int(stats["within_tolerance"] or 0)
        invalid = lt_zero + gt_one
        low_info = int(stats["low_info"] or 0)
        good_info = int(stats["good_info"] or 0)

        qc_info["total_variant_infile"] = total
        qc_info["missing_info"] = missing
        qc_info["lt_zero"] = lt_zero
        qc_info["gt_one"] = gt_one
        qc_info["gt_1_05"] = gt_1_05
        qc_info["invalid"] = invalid
        qc_info["low_info"] = low_info
        qc_info["good_info"] = good_info
        qc_info["min_info"] = stats["min"]
        qc_info["max_info"] = stats["max"]
        qc_info["mean_info"] = stats["mean"]
        qc_info["median_info"] = stats["median"]
        qc_info["info_clip_min"] = clip_min
        qc_info["info_clip_max"] = clip_max
        qc_info["info_clip_tolerance"] = tolerance_limit
        qc_info["info_out_of_range_action"] = out_of_range
        qc_info["info_on_missing_action"] = on_missing

        step.info(
            "Imputation quality ranges from {} to {} (mean {}, median {}). "
            "{:,} variants are at or below {:g}, which is poor imputation, and "
            "{:,} are above it.".format(
                stats["min"], stats["max"], stats["mean"], stats["median"],
                low_info, low_quality, good_info,
            )
        )

        # ------------------------------------------------------------------
        # ROUNDING ERROR JUST ABOVE THE CEILING
        # ------------------------------------------------------------------
        if within_tolerance > 0:
            before = df.height
            df = df.with_columns(
                pl.when(
                    (pl.col(info_col) > clip_max)
                    & (pl.col(info_col) <= tolerance_limit)
                )
                .then(pl.lit(clip_max))
                .otherwise(pl.col(info_col))
                .alias(info_col)
            )
            step.qc(
                "imputation quality rounding error",
                "{:,} variants have a quality just above {:g}, no higher than "
                "{:g}. That is rounding error, so they were set to {:g}.".format(
                    within_tolerance, clip_max, tolerance_limit, clip_max
                ),
                before, df.height, changed=within_tolerance, step=STEP_LABEL,
            )
        qc_info["rescaled_within_tolerance"] = within_tolerance

        # ------------------------------------------------------------------
        # OUT OF RANGE
        # ------------------------------------------------------------------
        out_of_range_mask = (
            (pl.col(info_col) < clip_min) | (pl.col(info_col) > tolerance_limit)
        ).fill_null(False)
        n_out_of_range = lt_zero + gt_1_05
        qc_info["out_of_range"] = n_out_of_range
        qc_info["rejected_out_of_range"] = 0

        if n_out_of_range > 0:
            explanation = (
                "{:,} variants have an imputation quality below {:g} or above "
                "{:g}.".format(n_out_of_range, clip_min, tolerance_limit)
            )
            if out_of_range == "fail":
                raise ImputationQualityError(
                    explanation
                    + " Policy 'info.out_of_range' is 'fail', so the run stopped."
                )
            if out_of_range == "clip":
                before = df.height
                df = df.with_columns(
                    pl.col(info_col).clip(clip_min, clip_max).alias(info_col)
                )
                step.qc(
                    "imputation quality range",
                    explanation
                    + " Policy 'info.out_of_range' is 'clip', so they were "
                      "forced back to between {:g} and {:g}.".format(
                          clip_min, clip_max
                      ),
                    before, df.height, changed=n_out_of_range, warn=True,
                    step=STEP_LABEL,
                )
            elif out_of_range == "null":
                before = df.height
                df = df.with_columns(
                    pl.when(out_of_range_mask)
                    .then(pl.lit(None, dtype=pl.Float64))
                    .otherwise(pl.col(info_col))
                    .alias(info_col)
                )
                step.qc(
                    "imputation quality range",
                    explanation
                    + " Policy 'info.out_of_range' is 'null', so their quality "
                      "was blanked out.",
                    before, df.height, changed=n_out_of_range, warn=True,
                    step=STEP_LABEL,
                )
            elif out_of_range == "reject":
                df, removed = _remove(
                    df, out_of_range_mask, step, rejects, "info_out_of_range",
                    "quality outside {:g} to {:g}".format(clip_min, tolerance_limit),
                )
                qc_info["rejected_out_of_range"] = removed

        # ------------------------------------------------------------------
        # MISSING
        # ------------------------------------------------------------------
        missing_now = missing + (
            n_out_of_range if out_of_range == "null" else 0
        )
        qc_info["missing_info_after"] = missing_now
        qc_info["rejected_missing"] = 0
        if missing_now > 0:
            explanation = (
                "{:,} variants have no imputation quality.".format(missing_now)
            )
            if on_missing == "fail":
                raise ImputationQualityError(
                    explanation
                    + " Policy 'info.on_missing' is 'fail', so the run stopped."
                )
            if on_missing == "reject":
                df, removed = _remove(
                    df, pl.col(info_col).is_null(), step, rejects, "info_missing",
                    None,
                )
                qc_info["rejected_missing"] = removed
            else:
                step.qc(
                    "missing imputation quality",
                    explanation
                    + " Policy 'info.on_missing' is 'keep', so they were left "
                      "unchanged for the independently configured post-merge "
                      "modules.qc_summary rules.",
                    df.height, df.height, step=STEP_LABEL,
                )
        else:
            step.qc(
                "missing imputation quality",
                "Every variant has an imputation quality.",
                df.height, df.height, step=STEP_LABEL,
            )

        qc_info["final_variants"] = df.height
        step.set_rows(df.height, removed=rows_in - df.height)

    return df, qc_info, sample_column_dict


def _remove(df, mask, step, rejects, reason, detail):
    """Remove the rows the mask selects, recording them.

    The reject collector logs its own before/after counts, so nothing is
    logged here on that path.  Without a collector the variants are still
    removed - the policy asked for it - but the log says they were not written
    to the rejected-variants file.
    """
    mask = mask.fill_null(False)
    before = df.height
    if rejects is not None:
        df = rejects.reject(df, mask, STEP_LABEL, reason, detail=detail)
        return df, before - df.height

    kept = df.filter(~mask)
    removed = before - kept.height
    step.qc(
        reason.replace("_", " "),
        _REASON_TEXT[reason],
        before, kept.height, reason=reason, warn=removed > 0, step=STEP_LABEL,
    )
    if removed > 0:
        step.warn(
            "No rejected-variant collector was wired in, so these {:,} variants "
            "were removed without being written to the rejected-variants "
            "file.".format(removed)
        )
    return kept, removed
