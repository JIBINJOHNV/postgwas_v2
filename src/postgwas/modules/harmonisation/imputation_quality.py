"""Dataset-wide INFO scale resolution and chromosome step 11 harmonisation.

The INFO score is taken from the study's own column, a per-chromosome or
whole-genome reference file, or an explicit fixed value. Its numerical type is
resolved once from the complete selected source before chromosome fan-out.
Chromosome workers then apply the corresponding configured standard-INFO or
MaCH-Rsq range without independently re-detecting the type.

There is deliberately no INFO equivalent of ``eaf_degenerate``: an INFO of 0 is
scientifically poor but not mathematically fatal. Harmonisation applies its
configured range action here, and the QC module independently assesses the
retained value against ``modules.qc_summary.rules`` after VCF creation.

All messages go through the PipelineLogger supplied by the orchestrator; the
module never writes to stdout.
"""

import math
import os
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Sequence, Tuple

import polars as pl

from postgwas.core.io.delimiters import resolve_delimiter
from postgwas.core.values import optional_text

from .external_reference_staging import (
    ExternalReferenceStagingError,
    projected_reference_batches,
)
from .shared.runtime import reject_rows, resolve_policies, step_context
from .shared.allele_join import allele_oriented_left_join
from .shared.variant_columns import (
    canonical_chromosome_expression,
    canonicalize_variant_frame,
    has_canonical_variant_columns,
    read_reference_variant_table,
)

__all__ = [
    "harmonise_imputation_quality",
    "resolve_dataset_info_score_type",
    "resolve_info_score_type_from_batches",
    "ImputationQualityError",
    "POLICY_KEYS",
    "STEP_LABEL",
]


# The free-text step label recorded in the reject file.
STEP_LABEL = "11 info_harmonisation"

# Every policy this step reads.  Passed to logger.step() so the log shows the
# value and the plain-English explanation of each one.
POLICY_KEYS = [
    "info.source",
    "info.score_type",
    "info.out_of_range",
    "info.clip_min",
    "info.clip_max",
    "info.clip_tolerance",
    "info.mach_rsq_max",
    "info.auto_mach_rsq_fraction",
    "info.maximum_invalid_fraction",
    "info.low_quality_threshold",
    "info.on_missing",
    "external_reference.exact_duplicate_action",
    "external_reference.non_identical_duplicate_action",
    "chromosome.strip_chr_prefix",
    "chromosome.strip_leading_zero",
    "chromosome.rename_map",
]

_REASON_TEXT = {
    "info_out_of_range": (
        "The imputation quality is outside the configured window and policy "
        "'info.out_of_range' is 'reject'."
    ),
    "info_out_of_range_after_null": (
        "The imputation quality was outside the configured window, was set to "
        "missing by policy 'info.out_of_range' = 'null', and was subsequently "
        "removed by policy 'info.on_missing' = 'reject'."
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


_DETECTION_VALUE_COLUMN = "__postgwas_info_score_detection_value"


def _empty_score_statistics() -> dict[str, Any]:
    return {
        "source_rows": 0,
        "finite_values": 0,
        "unusable_values": 0,
        "below_minimum": 0,
        "standard_range_values": 0,
        "mach_range_values": 0,
        "above_mach_maximum": 0,
        "minimum": None,
        "maximum": None,
    }


def _summarise_info_score_batches(
    batches: Iterable[pl.DataFrame],
    info_column: str,
    *,
    policies,
) -> dict[str, Any]:
    """Aggregate INFO-scale evidence without materialising source values."""
    clip_min = float(policies.get("info.clip_min"))
    tolerance = float(policies.get("info.clip_tolerance"))
    mach_max = float(policies.get("info.mach_rsq_max"))
    combined = _empty_score_statistics()

    for batch in batches:
        if info_column not in batch.columns:
            raise ImputationQualityError(
                "The selected INFO source is missing configured score column %r."
                % info_column
            )
        values = batch.select(
            pl.col(info_column).cast(pl.Float64, strict=False).alias(
                _DETECTION_VALUE_COLUMN
            )
        )
        value = pl.col(_DETECTION_VALUE_COLUMN)
        finite = (value.is_not_null() & value.is_finite()).fill_null(False)
        finite_value = pl.when(finite).then(value).otherwise(
            pl.lit(None, dtype=pl.Float64)
        )
        current = values.select([
            pl.len().alias("source_rows"),
            finite.sum().alias("finite_values"),
            (~finite).sum().alias("unusable_values"),
            (finite & (value < clip_min)).sum().alias("below_minimum"),
            (
                finite & (value >= clip_min) & (value <= tolerance)
            ).sum().alias("standard_range_values"),
            (
                finite & (value > tolerance) & (value <= mach_max)
            ).sum().alias("mach_range_values"),
            (finite & (value > mach_max)).sum().alias("above_mach_maximum"),
            finite_value.min().alias("minimum"),
            finite_value.max().alias("maximum"),
        ]).to_dicts()[0]
        for key in (
            "source_rows", "finite_values", "unusable_values",
            "below_minimum", "standard_range_values", "mach_range_values",
            "above_mach_maximum",
        ):
            combined[key] += int(current[key] or 0)
        if current["minimum"] is not None:
            combined["minimum"] = (
                current["minimum"]
                if combined["minimum"] is None
                else min(combined["minimum"], current["minimum"])
            )
        if current["maximum"] is not None:
            combined["maximum"] = (
                current["maximum"]
                if combined["maximum"] is None
                else max(combined["maximum"], current["maximum"])
            )
    return combined


def _resolve_info_score_statistics(
    statistics: Mapping[str, Any],
    *,
    policies,
    source_description: str,
) -> dict[str, Any]:
    """Resolve one dataset-wide INFO type from configured fraction guards."""
    finite_count = int(statistics["finite_values"])
    if finite_count == 0:
        raise ImputationQualityError(
            "INFO score-type detection found no finite numeric values in %s. "
            "Verify the selected imputation-quality column and its missing-value "
            "encoding; no chromosome worker was started."
            % source_description
        )

    configured_type = str(policies.get("info.score_type"))
    mach_count = int(statistics["mach_range_values"])
    invalid_count = int(statistics["above_mach_maximum"])
    mach_fraction = mach_count / finite_count
    invalid_fraction = invalid_count / finite_count
    auto_fraction = float(policies.get("info.auto_mach_rsq_fraction"))
    maximum_invalid_fraction = float(
        policies.get("info.maximum_invalid_fraction")
    )
    mach_max = float(policies.get("info.mach_rsq_max"))

    if invalid_fraction > maximum_invalid_fraction:
        raise ImputationQualityError(
            "INFO validation stopped before chromosome processing: {:,} of "
            "{:,} finite values ({:.6%}) in {} exceed "
            "info.mach_rsq_max={:g}. This is above "
            "info.maximum_invalid_fraction={:.6%}. Verify "
            "imputation_info_column/external_info_column and info.score_type; "
            "the selected field may not be an imputation-quality score."
            .format(
                invalid_count,
                finite_count,
                invalid_fraction,
                source_description,
                mach_max,
                maximum_invalid_fraction,
            )
        )

    if configured_type == "auto":
        resolved_type = (
            "mach_rsq" if mach_fraction >= auto_fraction else "standard_info"
        )
        decision_source = "automatic_full_dataset"
    else:
        resolved_type = configured_type
        decision_source = "explicit_configuration"

    result = dict(statistics)
    result.update({
        "configured_type": configured_type,
        "resolved_type": resolved_type,
        "decision_source": decision_source,
        "source": source_description,
        "mach_range_fraction": mach_fraction,
        "invalid_fraction": invalid_fraction,
        "auto_mach_rsq_fraction": auto_fraction,
        "maximum_invalid_fraction": maximum_invalid_fraction,
        "standard_info_maximum": float(policies.get("info.clip_tolerance")),
        "mach_rsq_maximum": mach_max,
    })
    return result


def resolve_info_score_type_from_batches(
    batches: Iterable[pl.DataFrame],
    info_column: str,
    *,
    policies=None,
    source_description: str = "the selected INFO source",
) -> dict[str, Any]:
    """Resolve INFO versus MaCH Rsq once from a bounded dataset-wide scan."""
    resolved_policies = resolve_policies(policies)
    statistics = _summarise_info_score_batches(
        batches, info_column, policies=resolved_policies,
    )
    return _resolve_info_score_statistics(
        statistics,
        policies=resolved_policies,
        source_description=source_description,
    )


def _external_info_batches(
    resource_maps: Mapping[str, Mapping[str, Any]],
    chromosomes: Sequence[str],
    *,
    info_column: str,
    external_info_colmap: Mapping[str, Any],
    external_reference_staging: Mapping[str, Any],
    policies,
) -> Iterable[pl.DataFrame]:
    chromosome_column = optional_text(external_info_colmap.get("chr"))
    configured_delimiter = optional_text(external_info_colmap.get("delimiter"))
    if chromosome_column is None or configured_delimiter is None:
        raise ImputationQualityError(
            "External INFO mapping requires configured 'chr' and 'delimiter' fields."
        )

    paths: dict[Path, set[str]] = {}
    for chromosome in chromosomes:
        path_text = optional_text(
            (resource_maps.get(str(chromosome)) or {}).get("user_info_file")
        )
        if path_text is None:
            raise ImputationQualityError(
                "External INFO source has no preflighted file for chromosome %s."
                % chromosome
            )
        path = Path(path_text).expanduser().resolve()
        paths.setdefault(path, set()).add(str(chromosome))

    selected_columns = tuple(dict.fromkeys([chromosome_column, info_column]))
    for path, path_chromosomes in paths.items():
        separator = None
        if path.suffix.lower() != ".parquet":
            detected = resolve_delimiter(
                path,
                configured_delimiter,
                candidates=list(policies.get("input.delimiter_candidates")),
                minimum_columns=int(
                    policies.get("input.delimiter_min_columns")
                ),
                maximum_columns=int(
                    policies.get("input.delimiter_max_columns")
                ),
                sample_lines=int(policies.get("input.delimiter_sample_rows")),
            )
            separator = detected.value
            if separator is None:
                raise ImputationQualityError(
                    "Could not resolve the delimiter for external INFO source %s."
                    % path
                )
        try:
            for batch in projected_reference_batches(
                path,
                columns=selected_columns,
                separator=separator,
                null_values=list(policies.get("input.null_values")),
                infer_schema_length=int(
                    policies.get("input.schema_inference_rows")
                ),
                batch_rows=int(external_reference_staging["batch_rows"]),
                compressed_suffixes=tuple(
                    external_reference_staging["compressed_suffixes"]
                ),
            ):
                yield batch.filter(
                    canonical_chromosome_expression(
                        pl.col(chromosome_column), policies,
                    ).is_in(sorted(path_chromosomes))
                )
        except ExternalReferenceStagingError as exc:
            raise ImputationQualityError(str(exc)) from exc


def resolve_dataset_info_score_type(
    *,
    df: pl.DataFrame | None,
    sample_column_dict: Mapping[str, Any],
    resource_maps: Mapping[str, Mapping[str, Any]] | None,
    chromosomes: Sequence[str],
    external_info_column: str | None,
    external_info_colmap: Mapping[str, Any],
    external_reference_staging: Mapping[str, Any],
    policies=None,
) -> dict[str, Any]:
    """Resolve the selected INFO source's score type before chromosome fan-out."""
    resolved_policies = resolve_policies(policies)
    source = optional_text(sample_column_dict.get("info_source"))
    if source == "internal":
        info_column = optional_text(sample_column_dict.get("imp_info_col"))
        if df is None or info_column is None or info_column not in df.columns:
            raise ImputationQualityError(
                "Dataset-wide INFO detection requires the selected internal "
                "imputation-quality column before chromosome partitioning."
            )
        return resolve_info_score_type_from_batches(
            (df.select(info_column),),
            info_column,
            policies=resolved_policies,
            source_description="study column %r" % info_column,
        )

    if source == "fixed_cli":
        fixed_value = sample_column_dict.get("fixed_info")
        if fixed_value is None:
            raise ImputationQualityError(
                "Dataset-wide INFO detection requires the resolved --fixed-info value."
            )
        numeric = float(fixed_value)
        row_count = int(df.height if df is not None else 1)
        finite = math.isfinite(numeric)
        clip_min = float(resolved_policies.get("info.clip_min"))
        tolerance = float(resolved_policies.get("info.clip_tolerance"))
        mach_max = float(resolved_policies.get("info.mach_rsq_max"))
        statistics = {
            "source_rows": row_count,
            "finite_values": row_count if finite else 0,
            "unusable_values": 0 if finite else row_count,
            "below_minimum": row_count if finite and numeric < clip_min else 0,
            "standard_range_values": (
                row_count if finite and clip_min <= numeric <= tolerance else 0
            ),
            "mach_range_values": (
                row_count if finite and tolerance < numeric <= mach_max else 0
            ),
            "above_mach_maximum": (
                row_count if finite and numeric > mach_max else 0
            ),
            "minimum": numeric if finite else None,
            "maximum": numeric if finite else None,
        }
        return _resolve_info_score_statistics(
            statistics,
            policies=resolved_policies,
            source_description="the explicit --fixed-info value",
        )

    if source == "external":
        info_column = optional_text(external_info_column)
        if info_column is None or resource_maps is None:
            raise ImputationQualityError(
                "Dataset-wide INFO detection requires the preflighted external "
                "INFO files and configured external_info_column."
            )
        batches = _external_info_batches(
            resource_maps,
            chromosomes,
            info_column=info_column,
            external_info_colmap=external_info_colmap,
            external_reference_staging=external_reference_staging,
            policies=resolved_policies,
        )
        return resolve_info_score_type_from_batches(
            batches,
            info_column,
            policies=resolved_policies,
            source_description="external INFO column %r" % info_column,
        )

    raise ImputationQualityError(
        "Dataset-wide INFO detection requires a resolved source of 'internal', "
        "'external', or 'fixed_cli'; received %r." % source
    )


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
    policies : a Policies object. When None the schema-validated registry
               defaults are used.
    rejects  : a RejectCollector.  Required for the 'reject' actions to record
               the removed variants.
    ctx      : a StepContext, when the caller opened the step itself.
    """
    policies = resolve_policies(policies)
    out_of_range = policies.get("info.out_of_range")
    clip_min = float(policies.get("info.clip_min"))
    clip_max = float(policies.get("info.clip_max"))
    tolerance_limit = float(policies.get("info.clip_tolerance"))
    mach_rsq_max = float(policies.get("info.mach_rsq_max"))
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
        configured_score_type = str(policies.get("info.score_type"))
        resolved_score_type = optional_text(
            sample_column_dict.get("resolved_info_score_type")
        )
        if resolved_score_type is None and configured_score_type != "auto":
            resolved_score_type = configured_score_type
        if resolved_score_type is None and (
            resolved_source == "fixed_cli" or fixed_info is not None
        ):
            resolved_score_type = "standard_info"
        if resolved_score_type not in ("standard_info", "mach_rsq"):
            raise ImputationQualityError(
                "INFO score type was not resolved at the complete-dataset level "
                "before chromosome %s entered the worker. Set info.score_type "
                "explicitly or run the dataset preparation stage; per-chromosome "
                "automatic detection is not permitted."
                % chromosome
            )
        upper_limit = (
            tolerance_limit
            if resolved_score_type == "standard_info"
            else mach_rsq_max
        )

        # ------------------------------------------------------------------
        # NORMALISE THE STUDY COORDINATES
        # ------------------------------------------------------------------
        if not has_canonical_variant_columns(df, sample_column_dict):
            df, _normalization = canonicalize_variant_frame(
                df,
                {"chr": chr_col, "pos": pos_col, "ea": ea_col, "oa": oa_col},
                policies,
                error_type=ImputationQualityError,
                warn=step.warn,
                label="study",
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
            info_raw, _detected = read_reference_variant_table(
                info_file_to_use,
                external_info_colmap,
                policies,
                value_columns=[info_col_to_use],
                error_type=ImputationQualityError,
                description="external imputation-quality table",
            )
            ref_chr = external_info_colmap["chr"]
            ref_pos = external_info_colmap["pos"]
            ref_effect = external_info_colmap["a1"]
            ref_other = external_info_colmap["a2"]
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
                policies=policies,
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
            if not clip_min <= fixed_value <= clip_max:
                raise ImputationQualityError(
                    "Fixed INFO must be between info.clip_min={:g} and "
                    "info.clip_max={:g}; received {!r}."
                    .format(clip_min, clip_max, fixed_info)
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
        qc_info["info_score_type"] = resolved_score_type
        qc_info["info_score_type_configured"] = configured_score_type
        qc_info["info_score_type_source"] = sample_column_dict.get(
            "info_score_type_source",
            "explicit_configuration"
            if configured_score_type != "auto"
            else "fixed_info",
        )
        qc_info["info_valid_maximum"] = upper_limit
        numeric_info = pl.col(info_col).cast(pl.Float64, strict=False)
        df = df.with_columns(
            pl.when(numeric_info.is_finite())
            .then(numeric_info)
            .otherwise(pl.lit(None, dtype=pl.Float64))
            .alias(info_col)
        )
        step.info(
            "The dataset-wide INFO score type is {}; this chromosome uses the "
            "configured valid window {:g} to {:g}."
            .format(resolved_score_type, clip_min, upper_limit)
        )

        # ------------------------------------------------------------------
        # QC METRICS - measured before anything is changed
        # ------------------------------------------------------------------
        total = df.height
        stats = df.select([
            pl.col(info_col).is_null().sum().alias("missing"),
            (pl.col(info_col) < clip_min).sum().alias("lt_zero"),
            (pl.col(info_col) > clip_max).sum().alias("gt_one"),
            (pl.col(info_col) > tolerance_limit).sum().alias(
                "gt_standard_tolerance"
            ),
            (pl.col(info_col) > upper_limit).sum().alias("gt_valid_maximum"),
            (
                (pl.col(info_col) > clip_max)
                & (pl.col(info_col) <= tolerance_limit)
                & pl.lit(resolved_score_type == "standard_info")
            ).sum().alias("within_tolerance"),
            (
                pl.col(info_col).is_not_null()
                & (pl.col(info_col) <= low_quality)
            ).sum().alias("low_info"),
            (
                (pl.col(info_col) > low_quality)
                & (pl.col(info_col) <= upper_limit)
            ).sum().alias("good_info"),
            pl.col(info_col).min().alias("min"),
            pl.col(info_col).max().alias("max"),
            pl.col(info_col).mean().alias("mean"),
            pl.col(info_col).median().alias("median"),
        ]).to_dicts()[0]
        missing = int(stats["missing"] or 0)
        lt_zero = int(stats["lt_zero"] or 0)
        gt_one = int(stats["gt_one"] or 0)
        gt_standard_tolerance = int(stats["gt_standard_tolerance"] or 0)
        gt_valid_maximum = int(stats["gt_valid_maximum"] or 0)
        within_tolerance = int(stats["within_tolerance"] or 0)
        invalid = lt_zero + gt_valid_maximum
        low_info = int(stats["low_info"] or 0)
        good_info = int(stats["good_info"] or 0)

        qc_info["total_variant_infile"] = total
        qc_info["missing_info"] = missing
        qc_info["lt_zero"] = lt_zero
        qc_info["gt_one"] = gt_one
        qc_info["gt_standard_info_tolerance"] = gt_standard_tolerance
        qc_info["gt_valid_maximum"] = gt_valid_maximum
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
        qc_info["info_mach_rsq_max"] = mach_rsq_max
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
        if resolved_score_type == "standard_info" and within_tolerance > 0:
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
            (pl.col(info_col) < clip_min) | (pl.col(info_col) > upper_limit)
        ).fill_null(False)
        n_out_of_range = lt_zero + gt_valid_maximum
        # Preserve the scientific origin only for the policy combination that
        # needs it. Otherwise the default direct-rejection path should not pay
        # for an additional materialised chromosome-length mask.
        out_of_range_rows = (
            df.select(out_of_range_mask).to_series()
            if n_out_of_range > 0
            and out_of_range == "null"
            and on_missing == "reject"
            else None
        )
        qc_info["out_of_range"] = n_out_of_range
        qc_info["rejected_out_of_range"] = 0
        qc_info["nullified_out_of_range"] = 0

        if n_out_of_range > 0:
            explanation = (
                "{:,} variants have an imputation quality below {:g} or above "
                "{:g} for the resolved {} score type."
                .format(
                    n_out_of_range,
                    clip_min,
                    upper_limit,
                    resolved_score_type,
                )
            )
            if out_of_range == "fail":
                raise ImputationQualityError(
                    explanation
                    + " Policy 'info.out_of_range' is 'fail', so the run stopped."
                )
            if out_of_range == "null":
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
                qc_info["nullified_out_of_range"] = n_out_of_range
            elif out_of_range == "reject":
                df, removed = reject_rows(
                    df,
                    out_of_range_mask,
                    step_label=STEP_LABEL,
                    reason="info_out_of_range",
                    context=step,
                    collector=rejects,
                    detail="quality outside {:g} to {:g}".format(
                        clip_min, upper_limit
                    ),
                    check_name="info out of range",
                    description=_REASON_TEXT["info_out_of_range"],
                    warn_on_remove=True,
                    warn_without_collector=True,
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
            if out_of_range == "null" and n_out_of_range > 0:
                explanation = (
                    "{:,} variants have no imputation quality after range "
                    "handling: {:,} were already missing or non-finite and "
                    "{:,} were set to missing because their INFO was outside "
                    "the valid window."
                    .format(missing_now, missing, n_out_of_range)
                )
            else:
                explanation = (
                    "{:,} variants have no imputation quality."
                    .format(missing_now)
                )
            if on_missing == "fail":
                raise ImputationQualityError(
                    explanation
                    + " Policy 'info.on_missing' is 'fail', so the run stopped."
                )
            if on_missing == "reject":
                if out_of_range == "null" and n_out_of_range > 0:
                    if out_of_range_rows is None:
                        raise AssertionError(
                            "The preserved out-of-range INFO mask is missing."
                        )
                    df, removed = reject_rows(
                        df,
                        out_of_range_rows,
                        step_label=STEP_LABEL,
                        reason="info_out_of_range",
                        context=step,
                        collector=rejects,
                        detail=(
                            "quality outside {:g} to {:g}; set to missing by "
                            "info.out_of_range='null' and removed by "
                            "info.on_missing='reject'"
                        ).format(clip_min, upper_limit),
                        check_name="nullified out-of-range INFO",
                        description=_REASON_TEXT[
                            "info_out_of_range_after_null"
                        ],
                        warn_on_remove=True,
                        warn_without_collector=True,
                    )
                    qc_info["rejected_out_of_range"] += removed
                df, removed = reject_rows(
                    df,
                    pl.col(info_col).is_null(),
                    step_label=STEP_LABEL,
                    reason="info_missing",
                    context=step,
                    collector=rejects,
                    detail=None,
                    check_name="info missing",
                    description=_REASON_TEXT["info_missing"],
                    warn_on_remove=True,
                    warn_without_collector=True,
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
