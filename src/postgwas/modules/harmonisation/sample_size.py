"""Step 05 - sample sizes and the effective sample size.

Reads the case and control counts, from a column or from a single configured
number, handles per-variant missing counts under an explicit dataset-wide
policy, and works out the effective sample size Neff.

The main formula is

    Neff = 4 / (1/Ncase + 1/Ncontrol)

which is correct and is not changed here.

Policies read here
------------------
sample_size.trait_type       case/control study or quantitative trait
sample_size.controls_only    what a lone control count means
sample_size.cases_only       what a lone case count means
sample_size.min_value        smallest sample size accepted
sample_size.missing_action   remove, impute or fail on missing counts
sample_size.max_missing_fraction  largest fraction eligible for that action
"""

from typing import Any, Dict, Optional, Tuple

import polars as pl

from postgwas.core.values import optional_text, parse_integer

from .policies import default_policies
from .shared.runtime import reject_rows, step_context


__all__ = [
    "effective_sample_size_expression",
    "sample_count_expression",
    "prepare_missing_sample_sizes",
    "harmonise_sample_sizes",
]


STEP_LABEL = "05 sample_size"

POLICY_KEYS = [
    "sample_size.trait_type",
    "sample_size.controls_only",
    "sample_size.cases_only",
    "sample_size.min_value",
    "sample_size.missing_action",
    "sample_size.max_missing_fraction",
]


# =============================================================================
# small helpers
# =============================================================================

_clean_name = optional_text
_parse_count = parse_integer


def sample_count_expression(df, column):
    """Read one sample-count column exactly as the chromosome step does."""
    expr = pl.col(column)
    if df.schema[column] == pl.String:
        expr = expr.str.replace_all(r"[,\s_]", "")
    return expr.cast(pl.Float64, strict=False).cast(pl.Int64, strict=False)


def effective_sample_size_expression(ncase, ncontrol):
    """Polars expression for the standard case-control effective sample size."""
    cases = ncase.cast(pl.Float64, strict=False)
    controls = ncontrol.cast(pl.Float64, strict=False)
    return 4 / (1 / cases + 1 / controls)


def _to_counts(df, column, ctx, what):
    # type: (pl.DataFrame, str, Any, str) -> Tuple[pl.DataFrame, int]
    """Cast a count column to whole numbers, and say how many could not be read.

    Text is cleaned of thousands separators first and goes through Float64, so
    '50,000', '50000.0' and '5e4' all survive; a direct cast to Int64 turns all
    three into nulls.
    """
    before_null = int(df.select(pl.col(column).null_count()).item() or 0)
    df = df.with_columns(
        sample_count_expression(df, column).alias(column)
    )
    after_null = int(df.select(pl.col(column).null_count()).item() or 0)
    unreadable = after_null - before_null
    if unreadable > 0:
        ctx.warn(
            "%s of the %s values in '%s' could not be read as a number and are now empty."
            % ("{:,}".format(unreadable), what, column)
        )
    if after_null == df.height and df.height:
        ctx.warn(
            "Every value in the %s column '%s' is empty. Nothing that depends on it can be "
            "worked out." % (what, column)
        )
    return df, unreadable


def _configured_count_columns(df, sample_column_dict):
    """Return configured per-variant count columns in case/control order.

    A fixed configured count cannot be missing per variant, so it is not part
    of the dataset missingness denominator or the imputation plan.
    """
    columns = []
    for column_key in ("ncase_col", "ncontrol_col"):
        column = _clean_name(sample_column_dict.get(column_key))
        if column is not None and column in df.columns:
            if column not in columns:
                columns.append(column)
    return columns


def prepare_missing_sample_sizes(
    df,
    sample_column_dict,
    policies=None,
    ctx=None,
    scope="dataset",
):
    # type: (pl.DataFrame, dict, Any, Any, str) -> Dict[str, Any]
    """Build the one missing-sample-size plan used by every chromosome.

    The fraction and optional fill values are calculated over the complete
    dataset before it is split.  This avoids chromosome-specific imputation
    values.  The function does not alter or remove rows: chromosome workers
    apply the plan so removed variants retain their reject-file provenance.

    Missing-N imputation is deliberately opt-in. MungeSumstats likewise drops
    missing N by default and describes sample-size imputation as a last resort:
    https://www.bioconductor.org/packages/release/bioc/vignettes/MungeSumstats/inst/doc/MungeSumstats.html
    """
    policies = default_policies() if policies is None else policies
    action = str(policies.get("sample_size.missing_action"))
    max_fraction = float(policies.get("sample_size.max_missing_fraction"))
    min_value = int(policies.get("sample_size.min_value"))
    columns = _configured_count_columns(df, sample_column_dict)
    total = df.height

    if not columns:
        return {
            "scope": scope,
            "action": action,
            "max_fraction": max_fraction,
            "total_variants": total,
            "missing_variants": 0,
            "missing_fraction": 0.0,
            "columns": [],
            "missing_by_column": {},
            "fill_values": {},
        }

    expressions = []
    missing_terms = []
    for index, column in enumerate(columns):
        parsed = sample_count_expression(df, column)
        missing = parsed.is_null()
        missing_terms.append(missing)
        expressions.extend([
            missing.sum().alias("__missing_%d" % index),
        ])
        if action in ("median", "mean"):
            usable = parsed.filter(parsed >= min_value)
            expressions.append(
                getattr(usable, action)().alias("__fill_%d" % index)
            )
    expressions.append(
        pl.any_horizontal(missing_terms).sum().alias("__missing_variants")
    )
    summary = df.select(expressions).to_dicts()[0]
    missing_variants = int(summary["__missing_variants"] or 0)
    missing_fraction = (float(missing_variants) / total) if total else 0.0
    missing_by_column = {
        column: int(summary["__missing_%d" % index] or 0)
        for index, column in enumerate(columns)
    }

    if missing_fraction > max_fraction:
        raise ValueError(
            "Sample-size validation failed for %s: %s of %s variants (%.2f%%) "
            "have a missing case, control or total sample size, exceeding "
            "sample_size.max_missing_fraction = %.2f%%. Provide a more complete "
            "sample-size column or increase the configured maximum."
            % (
                scope,
                "{:,}".format(missing_variants),
                "{:,}".format(total),
                100.0 * missing_fraction,
                100.0 * max_fraction,
            )
        )
    if action == "fail" and missing_variants:
        raise ValueError(
            "Sample-size validation failed for %s: %s of %s variants (%.2f%%) "
            "have a missing case, control or total sample size and "
            "sample_size.missing_action is 'fail'."
            % (
                scope,
                "{:,}".format(missing_variants),
                "{:,}".format(total),
                100.0 * missing_fraction,
            )
        )

    fill_values = {}
    if missing_variants and action in ("median", "mean"):
        for index, column in enumerate(columns):
            if not missing_by_column[column]:
                continue
            value = summary["__fill_%d" % index]
            if value is None:
                raise ValueError(
                    "Sample-size validation failed for %s: column '%s' has no usable "
                    "values, so its missing values cannot be filled with the %s."
                    % (scope, column, action)
                )
            # Counts must remain whole individuals. This matches the existing
            # sample-size normalisation and the GWAS-VCF integer NEF contract.
            fill_values[column] = int(round(float(value)))

    plan = {
        "scope": scope,
        "action": action,
        "max_fraction": max_fraction,
        "total_variants": total,
        "missing_variants": missing_variants,
        "missing_fraction": missing_fraction,
        "columns": columns,
        "missing_by_column": missing_by_column,
        "fill_values": fill_values,
    }
    if ctx is not None:
        ctx.info(
            "Sample-size missingness across %s: %s of %s variants (%.2f%%); "
            "configured maximum %.2f%%; action '%s'."
            % (
                scope,
                "{:,}".format(missing_variants),
                "{:,}".format(total),
                100.0 * missing_fraction,
                100.0 * max_fraction,
                action,
            )
        )
        if fill_values:
            ctx.decide(
                "Missing sample-size imputation",
                {"statistic": action, "dataset fill values": fill_values},
                "missing values are filled in chromosome workers using these dataset-wide values",
            )
    return plan


def _column_stats(df, columns):
    # type: (pl.DataFrame, list[str]) -> Dict[str, Dict[str, Any]]
    """Summarize all sample-size columns in one aggregate pass.

    Without null_count a column that is entirely empty reports a total equal to
    the number of variants, min and max of None, and looks healthy.
    """
    unique_columns = list(dict.fromkeys(columns))
    expressions = []
    for index, column in enumerate(unique_columns):
        prefix = "__sample_size_%d" % index
        expressions.extend([
            pl.col(column).min().alias(prefix + "_min"),
            pl.col(column).max().alias(prefix + "_max"),
            pl.col(column).mean().alias(prefix + "_mean"),
            pl.col(column).median().alias(prefix + "_median"),
            pl.col(column).null_count().alias(prefix + "_null_count"),
            pl.col(column).is_not_null().sum().alias(prefix + "_non_null"),
        ])
    row = df.select(expressions).to_dicts()[0] if expressions else {}
    total = df.height
    result = {}
    for index, column in enumerate(unique_columns):
        prefix = "__sample_size_%d" % index
        nulls = int(row[prefix + "_null_count"] or 0)
        result[column] = {
            "min": row[prefix + "_min"],
            "max": row[prefix + "_max"],
            "mean": row[prefix + "_mean"],
            "median": row[prefix + "_median"],
            "total": total,
            "null_count": nulls,
            "non_null": int(row[prefix + "_non_null"] or 0),
            "null_fraction": (float(nulls) / total) if total else 0.0,
        }
    return result


# =============================================================================
# the step
# =============================================================================

def harmonise_sample_sizes(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    logger: Any = None,
    policies: Any = None,
    rejects: Any = None,
    ctx: Any = None,
    missing_plan: Optional[dict] = None,
    step_number: int = 5,
    step_total: int = 16,
) -> Tuple[pl.DataFrame, Dict, dict]:
    """Harmonise the case and control counts and compute Neff.

    Returns ``(df, qc_info, sample_column_dict)``.

    Parameters
    ----------
    logger : the pipeline logger; direct library calls may omit it.
    policies : a Policies object; the defaults are used when it is None.
    rejects : a RejectCollector.  Without one the variants are still removed and
        still counted, they are just not written to the reject file.
    ctx : an open StepContext, when the caller has already opened the step.
    """
    policies = default_policies() if policies is None else policies

    rows_in = df.height
    counters = {}                                   # type: Dict[str, int]
    qc_info = {"initial_variants": rows_in}         # type: Dict[str, Any]

    with step_context(
        logger, ctx,
        number=step_number, total=step_total,
        title="Sample size", operation="sample_size.harmonise_sample_sizes",
        rows_in=rows_in, policy_keys=POLICY_KEYS,
    ) as step_ctx:
        df, sample_column_dict = _harmonise(
            df=df,
            sample_column_dict=sample_column_dict,
            policies=policies,
            ctx=step_ctx,
            rejects=rejects,
            counters=counters,
            qc_info=qc_info,
            chromosome=chromosome,
            missing_plan=missing_plan,
        )
        qc_info["removed_by_reason"] = dict(counters)
        qc_info["final_variants"] = df.height
        step_ctx.set_rows(df.height)
        step_ctx.extra.update(qc_info)

    return df, qc_info, sample_column_dict


def _harmonise(
    df, sample_column_dict, policies, ctx, rejects, counters, qc_info,
    chromosome, missing_plan,
):
    # type: (pl.DataFrame, dict, Any, Any, Any, Dict[str, int], Dict[str, Any]) -> Tuple[pl.DataFrame, dict]
    ncontrol_col = _clean_name(sample_column_dict.get("ncontrol_col"))
    ncase_col = _clean_name(sample_column_dict.get("ncase_col"))
    ncontrol_val = sample_column_dict.get("ncontrol")
    ncase_val = sample_column_dict.get("ncase")

    trait_type = str(policies.sample_size.trait_type)
    min_value = policies.sample_size.min_value
    min_value = None if min_value is None else int(min_value)
    qc_info["trait_type"] = trait_type
    qc_info["min_value"] = min_value

    if trait_type == "quantitative" and (
        ncase_col is not None or _clean_name(ncase_val) is not None
    ):
        raise ValueError(
            "sample_size.trait_type is 'quantitative', so case-count inputs are not allowed. "
            "Provide the quantitative study's total sample size through the control-count "
            "column or fixed control-count field."
        )

    # ------------------------------------------------------------------
    # 1. controls
    # ------------------------------------------------------------------
    if ncontrol_col is not None and ncontrol_col in df.columns:
        ctx.info("Control counts come from the column '%s'." % ncontrol_col)
        df, unreadable = _to_counts(df, ncontrol_col, ctx, "control count")
        qc_info["ncontrol_source"] = "column:%s" % ncontrol_col
        qc_info["ncontrol_unreadable"] = unreadable
    elif _clean_name(ncontrol_val) is not None:
        parsed = _parse_count(ncontrol_val)
        if parsed is None:
            raise ValueError(
                "The control count in the configuration, %r, is not a number. Write it as a "
                "plain integer such as 50000; 50,000 and 5e4 are also accepted."
                % (ncontrol_val,)
            )
        ctx.info(
            "Every variant is given the same control count, %s, from the configuration."
            % "{:,}".format(parsed)
        )
        df = df.with_columns(pl.lit(parsed, dtype=pl.Int64).alias("ncontrol"))
        ncontrol_col = "ncontrol"
        sample_column_dict["ncontrol_col"] = "ncontrol"
        qc_info["ncontrol_source"] = "fixed_value:%d" % parsed
    else:
        ctx.skip("No control count column or value is configured.")
        qc_info["ncontrol_source"] = "missing"

    # ------------------------------------------------------------------
    # 2. cases
    # ------------------------------------------------------------------
    if ncase_col is not None and ncase_col in df.columns:
        ctx.info("Case counts come from the column '%s'." % ncase_col)
        df, unreadable = _to_counts(df, ncase_col, ctx, "case count")
        qc_info["ncase_source"] = "column:%s" % ncase_col
        qc_info["ncase_unreadable"] = unreadable
    elif _clean_name(ncase_val) is not None:
        parsed = _parse_count(ncase_val)
        if parsed is None:
            raise ValueError(
                "The case count in the configuration, %r, is not a number. Write it as a plain "
                "integer such as 12000; 12,000 and 1.2e4 are also accepted." % (ncase_val,)
            )
        ctx.info(
            "Every variant is given the same case count, %s, from the configuration."
            % "{:,}".format(parsed)
        )
        df = df.with_columns(pl.lit(parsed, dtype=pl.Int64).alias("ncase"))
        ncase_col = "ncase"
        sample_column_dict["ncase_col"] = "ncase"
        qc_info["ncase_source"] = "fixed_value:%d" % parsed
    else:
        ctx.skip("No case count column or value is configured.")
        qc_info["ncase_source"] = "missing"

    has_controls = ncontrol_col is not None and ncontrol_col in df.columns
    has_cases = ncase_col is not None and ncase_col in df.columns

    # ------------------------------------------------------------------
    # 3. missing per-variant counts
    # ------------------------------------------------------------------
    if missing_plan is None:
        missing_plan = prepare_missing_sample_sizes(
            df,
            sample_column_dict,
            policies=policies,
            ctx=ctx,
            scope="chromosome %s" % chromosome,
        )
    missing_columns = [
        column for column in missing_plan.get("columns", [])
        if column in df.columns
    ]
    chromosome_missing_summary = (
        df.select([
            pl.col(column).null_count().alias("__missing_%d" % index)
            for index, column in enumerate(missing_columns)
        ]).to_dicts()[0]
        if missing_columns else {}
    )
    chromosome_missing_by_column = {
        column: int(chromosome_missing_summary["__missing_%d" % index] or 0)
        for index, column in enumerate(missing_columns)
    }
    columns_missing_here = [
        column for column in missing_columns
        if chromosome_missing_by_column[column] > 0
    ]
    missing_mask = (
        pl.any_horizontal([pl.col(column).is_null() for column in columns_missing_here])
        if columns_missing_here else None
    )
    missing_rows = (
        int(df.select(missing_mask.sum()).item() or 0)
        if missing_mask is not None else 0
    )
    missing_action = str(missing_plan.get("action", policies.sample_size.missing_action))
    missing_qc = {
        "action": missing_action,
        "dataset_missing_fraction": missing_plan.get("missing_fraction"),
        "dataset_max_missing_fraction": missing_plan.get("max_fraction"),
        "variants_affected": missing_rows,
        "missing_by_column": chromosome_missing_by_column,
        "fill_values": {},
    }
    if missing_rows and missing_action == "remove":
        df, removed = reject_rows(
            df,
            missing_mask,
            step_label=STEP_LABEL,
            reason="sample_size_invalid",
            context=ctx,
            collector=rejects,
            counters=counters,
            detail="missing sample size; sample_size.missing_action = remove",
        )
        missing_qc["variants_removed"] = removed
    elif missing_rows and missing_action in ("median", "mean"):
        fill_values = dict(missing_plan.get("fill_values") or {})
        absent = [column for column in columns_missing_here if column not in fill_values]
        if absent:
            raise ValueError(
                "Chromosome %s has missing sample-size values in %s, but the "
                "dataset-wide %s plan has no fill value for those columns."
                % (chromosome, ", ".join(absent), missing_action)
            )
        df = df.with_columns([
            pl.col(column).fill_null(fill_values[column]).alias(column)
            for column in columns_missing_here
        ])
        missing_qc["variants_imputed"] = missing_rows
        missing_qc["fill_values"] = {
            column: fill_values[column] for column in columns_missing_here
        }
        ctx.info(
            "Filled %s variants with missing sample size using the dataset-wide %s: %s."
            % (
                "{:,}".format(missing_rows),
                missing_action,
                ", ".join(
                    "%s=%s" % (column, fill_values[column])
                    for column in columns_missing_here
                ),
            )
        )
    elif missing_rows:
        raise ValueError(
            "Chromosome %s still has %s variants with missing sample size under "
            "sample_size.missing_action = '%s'."
            % (chromosome, "{:,}".format(missing_rows), missing_action)
        )
    qc_info["missing_sample_size"] = missing_qc

    # ------------------------------------------------------------------
    # 4. remove counts that cannot be right
    # ------------------------------------------------------------------
    if min_value is not None:
        for present, column, what in ((has_cases, ncase_col, "case"),
                                      (has_controls, ncontrol_col, "control")):
            if not present:
                continue
            # Missing values were already handled by the explicit policy above;
            # this separate rule handles present values below the numeric floor.
            df, _ = reject_rows(
                df,
                pl.col(column).is_not_null() & (pl.col(column) < min_value),
                step_label=STEP_LABEL, reason="sample_size_invalid",
                context=ctx, collector=rejects, counters=counters,
                detail="%s count below sample_size.min_value (%d)" % (what, min_value),
            )

    # ------------------------------------------------------------------
    # 5. the effective sample size
    # ------------------------------------------------------------------
    neff_created = False

    if trait_type == "binary" and not (has_cases and has_controls):
        raise ValueError(
            "sample_size.trait_type is 'binary', which needs both a case count and a control "
            "count, but %s. Either configure the missing count or set sample_size.trait_type to "
            "'auto'." % ("neither is configured" if not (has_cases or has_controls)
                         else ("only the case count is configured" if has_cases
                               else "only the control count is configured"))
        )

    if has_controls and has_cases:
        ctx.info("Neff = 4 / (1/Ncase + 1/Ncontrol)")
        df = df.with_columns(
            effective_sample_size_expression(
                pl.col(ncase_col), pl.col(ncontrol_col),
            )
            .round(0)
            .cast(pl.Int64, strict=False)
            .alias("Neff")
        )
        neff_created = True
        qc_info["Neff_status"] = "calculated_from_case_control"

    elif has_controls:
        action = str(policies.sample_size.controls_only)
        if action == "fail":
            raise ValueError(
                "Only a control count is configured and sample_size.controls_only is 'fail'. "
                "Set it to 'use_ncontrol_as_total' if the lone count really is the total sample "
                "size, which is what it means for a quantitative trait."
            )
        ctx.decide(
            "Sample size with only a control count",
            {"column": ncontrol_col, "policy sample_size.controls_only": action},
            "the lone control count is the total sample size, so Neff := that count",
        )
        df = df.with_columns(pl.col(ncontrol_col).cast(pl.Int64, strict=False).alias("Neff"))
        neff_created = True
        qc_info["Neff_status"] = "fallback_ncontrol_only:%s" % ncontrol_col

    elif has_cases:
        raise ValueError(
            "Only a case count is configured and sample_size.cases_only is 'fail'. A case "
            "count alone cannot determine case/control effective sample size: using it as "
            "Neff would understate the balanced-design value by two times and the value with "
            "many controls by up to four times. Configure the control count."
        )

    else:
        ctx.warn(
            "Neither a case count nor a control count is available, so the effective sample size "
            "cannot be worked out and no Neff column is created."
        )
        qc_info["Neff_status"] = "not_computed_no_inputs"

    # The name is recorded only when the column was actually created. It used to
    # be set unconditionally, so every later step was told there was an Neff
    # column even when no branch had made one.
    if neff_created:
        sample_column_dict["neff_col"] = "Neff"
    else:
        sample_column_dict["neff_col"] = None

    # ------------------------------------------------------------------
    # 6. remove effective sample sizes that cannot be right
    # ------------------------------------------------------------------
    if neff_created and min_value is not None:
        df, _ = reject_rows(
            df,
            pl.col("Neff").is_not_null()
            & (~pl.col("Neff").is_finite() | (pl.col("Neff") < min_value)),
            step_label=STEP_LABEL, reason="neff_invalid",
            context=ctx, collector=rejects, counters=counters,
            detail="Neff below sample_size.min_value (%d), or not a finite number" % min_value,
        )

    # ------------------------------------------------------------------
    # 7. summary statistics
    # ------------------------------------------------------------------
    sample_stats = {}       # type: Dict[str, Any]
    configured_columns = [
        sample_column_dict.get(key)
        for key in ("ncase_col", "ncontrol_col", "neff_col")
    ]
    present_columns = [
        column for column in configured_columns
        if optional_text(column) is not None and column in df.columns
    ]
    stats_by_column = _column_stats(df, present_columns)
    for key in ("ncase_col", "ncontrol_col", "neff_col"):
        real_col = sample_column_dict.get(key)
        if optional_text(real_col) is None or real_col not in df.columns:
            sample_stats[key] = "missing"
            continue
        stats = stats_by_column[real_col]
        sample_stats[real_col] = stats
        ctx.info(
            "%s ('%s'): %s values, %s empty (%.1f%%); smallest %s, largest %s, median %s."
            % (key, real_col,
               "{:,}".format(int(stats["total"])),
               "{:,}".format(int(stats["null_count"])),
               100.0 * stats["null_fraction"],
               stats["min"], stats["max"], stats["median"])
        )
        if stats["null_count"] and stats["null_count"] == stats["total"]:
            ctx.warn(
                "Every value in '%s' is empty, so its summary statistics are meaningless."
                % real_col
            )

    qc_info["sample_size_stats"] = sample_stats
    return df, sample_column_dict
