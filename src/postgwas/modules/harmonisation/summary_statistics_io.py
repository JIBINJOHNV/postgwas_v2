import os
import re
import sys
import zipfile
from functools import partial
from pathlib import Path

import polars as pl

from postgwas.core.paths import configured_output_path
from postgwas.core.io.delimiters import open_text, resolve_delimiter

from postgwas.core.ui.screen import screen_field
from postgwas.core.values import optional_text

from .coordinates import harmonise_coordinates_and_alleles
from .rejects import RejectCollector, SOURCE_INPUT_ROW_COLUMN
from .sample_size import effective_sample_size_expression
from .shared.runtime import emit_message, reject_rows, resolve_policies
from .shared.variant_columns import (
    mark_canonical_variant_columns,
    study_string_schema,
)


# ======================================================================
# Policy / logging / reject plumbing shared by everything in this module
# ======================================================================
#: Which config key supplies each of the field names used by the list-valued
#: policies (columns.mandatory, duplicates.key).  'zscore', not 'z'.
FIELD_TO_CONFIG_KEY = {
    "chr": "chr_col",
    "pos": "pos_col",
    "snp": "snp_id_col",
    "ea": "ea_col",
    "oa": "oa_col",
    "eaf": "eaf_col",
    "beta": "beta_or_col",
    "se": "se_col",
    "zscore": "imp_z_col",
    "pval": "pval_col",
    "info": "imp_info_col",
    "n": "ncontrol_col",
    "ncontrol": "ncontrol_col",
    "ncase": "ncase_col",
}

#: Config keys that name columns read out of the study file.
COLUMN_CONFIG_KEYS = (
    "chr_col", "pos_col", "chr_pos_col", "snp_id_col",
    "ea_col", "oa_col", "eaf_col",
    "beta_or_col", "se_col", "imp_z_col", "pval_col",
    "ncontrol_col", "ncase_col", "imp_info_col",
)

#: Configured scientific measurements that may be converted to Float64 after
#: coordinate/allele normalization. Unknown input columns retain the type that
#: Polars inferred instead of being reinterpreted merely because they look numeric.
NUMERIC_COLUMN_CONFIG_KEYS = (
    "eaf_col", "beta_or_col", "se_col", "imp_z_col", "pval_col",
    "ncontrol_col", "ncase_col", "imp_info_col",
)

_emit = partial(emit_message, screen_when_unlogged=True)

DUPLICATE_REPORT_COLUMNS = (
    "duplicate_input_row",
    "duplicate_class",
    "duplicate_conflicting_fields",
    "duplicate_action",
)


def _field_column_pairs(fields, sample_column_dict, df_columns):
    """Resolved ``(policy field, study column)`` pairs in policy order."""
    pairs = []
    for field in fields or []:
        field_name = str(field).strip().lower()
        key = FIELD_TO_CONFIG_KEY.get(field_name)
        if not key:
            continue
        column = optional_text(sample_column_dict.get(key)) or ""
        if column and column in df_columns:
            pairs.append((field_name, column))
    return pairs


def _fields_to_columns(fields, sample_column_dict, df_columns):
    """Policy field names -> the study columns that carry them, de-duplicated."""
    columns = []
    for _field, column in _field_column_pairs(
        fields, sample_column_dict, df_columns,
    ):
        if column not in columns:
            columns.append(column)
    return columns


def _materialise_mask(df, expr):
    """Evaluate a boolean expression once so keep and reject are exact complements."""
    return df.select(expr.fill_null(True).alias("__mask")).get_column("__mask")


# ----------------------------------------------------------------------
# 4. Count data lines – **safe** shell command (list form)
# ----------------------------------------------------------------------
def count_data_lines(path: str, skip_hash: bool = True) -> int:
    """Count text records without shell commands or temporary files."""
    with open_text(path) as handle:
        return sum(
            1 for line in handle
            if not (skip_hash and line.startswith("##"))
        )

def load_summary_statistics_table(file_path: str, sample_column_dict: dict, policies=None, logger=None):
    policies = resolve_policies(policies)
    null_values = list(policies.get("input.null_values"))
    infer_rows = int(policies.get("input.schema_inference_rows"))
    comment_prefix = "##" if bool(policies.get("input.strip_double_hash_lines")) else None
    delimiter = resolve_delimiter(
        file_path,
        str(policies.get("input.delimiter")),
        candidates=list(policies.get("input.delimiter_candidates")),
        minimum_columns=int(policies.get("input.delimiter_min_columns")),
        maximum_columns=int(policies.get("input.delimiter_max_columns")),
        sample_lines=int(policies.get("input.delimiter_sample_rows")),
    )
    sep = delimiter.value
    string_schema = study_string_schema(sample_column_dict)
    if string_schema:
        _emit(
            logger,
            "Variant identity columns read as text to preserve labels exactly: %s."
            % ", ".join(string_schema),
        )
    _emit(
        logger,
        "Separator resolved by %s: %r (%s sampled columns)."
        % (delimiter.method, sep, delimiter.column_count or "configured"),
    )
    _emit(logger, "Detected delimiter: %r" % sep)
    try:
        df = pl.read_csv(
            file_path,
            separator=sep,
            infer_schema_length=infer_rows,
            ignore_errors=False,
            truncate_ragged_lines=False,
            has_header=True,
            null_values=null_values,
            comment_prefix=comment_prefix,
            schema_overrides=string_schema,
        )
    except Exception as e:
        _emit(logger, "Primary read failed: %s" % e, warn=True)
        _emit(logger, "Retrying with low_memory=True.", warn=True)
        df = pl.read_csv(
            file_path,
            separator=sep,
            infer_schema_length=infer_rows,
            low_memory=True,
            truncate_ragged_lines=False,
            has_header=True,
            null_values=null_values,
            comment_prefix=comment_prefix,
            schema_overrides=string_schema,
        )
    # ----------------------------
    # 4. Required columns  (de-duplicated: two config keys may name one column)
    # ----------------------------
    wanted = []
    for key in COLUMN_CONFIG_KEYS:
        column = optional_text(sample_column_dict.get(key)) or ""
        if column:
            wanted.append((key, column))
    required_cols = list(dict.fromkeys(column for _key, column in wanted))
    # ----------------------------
    # 5. Missing check
    # ----------------------------
    missing_cols = [c for c in required_cols if c not in df.columns]
    total_cols = len(df.columns)
    # Structural checks.  Polars names a repeated header 'NAME_duplicated_0' —
    # a SUFFIX — and never produces pandas' 'Unnamed: 0'; the detectors that
    # looked for a '_duplicated' prefix and an 'unnamed' prefix could not fire.
    empty_columns = [c for c in df.columns if str(c).strip() == ""]
    duplicated_columns = [c for c in df.columns if re.search(r"_duplicated_\d+$", str(c))]
    duplicate_action = policies.get("validation.on_duplicate_header")
    if duplicated_columns:
        _emit(
            logger,
            "The input file repeats %d column name(s); polars renamed the later "
            "copies to %s. Policy validation.on_duplicate_header is '%s'."
            % (
                len(duplicated_columns),
                ", ".join(duplicated_columns),
                duplicate_action,
            ),
            warn=True,
        )
    # ----------------------------
    # 6. STRUCTURE VALIDATION
    # ----------------------------
    structure_broken = (
        bool(missing_cols)
        or total_cols == 1
        or len(empty_columns) > 0
        or (bool(duplicated_columns) and duplicate_action == "fail")
    )
    if structure_broken:
        raise RuntimeError(
            "The summary-statistics table structure is invalid. Missing configured "
            "columns: %s; parsed columns: %d; empty column names: %s; repeated "
            "column names: %s. Check the delimiter and sample-sheet column names. "
            "PostGWAS will not rewrite the input because collapsing whitespace can "
            "shift empty fields and corrupt variant values."
            % (missing_cols, total_cols, empty_columns, duplicated_columns)
        )
    # ----------------------------
    # 7. FINAL SANITY
    # ----------------------------
    if df.height == 0:
        raise RuntimeError("❌ Empty dataframe")
    _emit(logger, "Loaded %d rows x %d columns" % (df.height, len(df.columns)))
    return df, sample_column_dict


def normalise_string_columns(
    df: pl.DataFrame,
    preserved_columns=None,
) -> pl.DataFrame:
    """
    1. Trim whitespace from all string columns
    2. Replace empty strings with null
    3. Preserve numeric precision (no unsafe downcasting)
    """
    preserved = set(preserved_columns or ())
    string_columns = [
        column
        for column, dtype in df.schema.items()
        if dtype == pl.String and column not in preserved
    ]
    if not string_columns:
        return df
    df = (
        df.lazy()
        # Step 1: trim whitespace
        .with_columns(
            pl.col(string_columns).str.strip_chars()
        )
        # Step 2: empty string → null
        .with_columns(
            pl.col(string_columns).replace("", None)
        )
        .collect()
    )
    return df


# ======================================================================
# HELPER 1: DataFrame Cleanup (Headers, Types, Chromosome)
# ======================================================================
def normalise_summary_statistics_values(
    df: pl.DataFrame,
    chr_col: str,
    policies=None,
    logger=None,
    strings_normalised=False,
    preserved_columns=None,
    numeric_columns=None,
) -> pl.DataFrame:
    """
    Normalize configured scientific numeric columns after the first read.

    Unconfigured columns retain the type inferred by Polars. Variant identity
    columns are forced to strings by :func:`load_summary_statistics_table`, so
    numeric-looking IDs never pass through Float64.

    A configured numeric column is converted only when every non-null value is
    parseable. This lossless gate is an input-integrity invariant: mixed fields
    such as comma-separated INFO values remain intact for their field-specific
    parser instead of becoming null here.
    """
    if df.height == 0:
        return df
    resolve_policies(policies)
    preserved = set(preserved_columns or ())
    numeric = set(numeric_columns or ())
    # 1. Clean Headers
    df = df.rename({c: c.strip() for c in df.columns})
    # 2. Clean String Values — one pass over every string column, not one
    #    frame rebuild per column. The dataset reader has already applied the
    #    stronger trim-and-empty-to-null normalization immediately beforehand.
    if not strings_normalised:
        df = df.with_columns(pl.col(pl.String).str.strip_chars())
    # 3. Type Conversion
    valid_cols = []
    for col in df.columns:
        if col == chr_col or col in preserved or col not in numeric:
            valid_cols.append(pl.col(col))
            continue
        original = df[col]
        numeric_series = original.cast(pl.Float64, strict=False)
        new_nulls = numeric_series.null_count() - original.null_count()
        if new_nulls <= 0:
            valid_cols.append(numeric_series)
            continue
        non_null = original.len() - original.null_count()
        # Preserve every original value when even one value is not parseable.
        # The field-specific step later decides whether those values are invalid
        # or use an intentional compound representation such as INFO lists.
        valid_cols.append(pl.col(col))
        if new_nulls < non_null:
            _emit(
                logger,
                "Configured numeric column '%s' stayed as text because %d of %d "
                "non-null values do not parse as numbers. Its field-specific "
                "validation step will handle those values."
                % (col, new_nulls, non_null),
                warn=True,
            )
    df = df.with_columns(valid_cols)
    # 4. Chromosome Fix (12.0 -> 12)
    if chr_col and chr_col in df.columns and chr_col not in preserved:
        df = df.with_columns(
            pl.coalesce(
                pl.col(chr_col).cast(pl.Float64, strict=False).cast(pl.Int64, strict=False).cast(pl.String),
                pl.col(chr_col)
            ).alias(chr_col)
        )
    return df


# =========================================================
# 4. READ SUMSTATS (CRITICAL FIXES)
# =========================================================

def normalise_imputation_quality_column(
    df: pl.DataFrame,
    sample_column_dict: dict,
    policies,
):
    """
    Process IMPINFO column:
    - If imp_info_col is missing/NA -> do nothing
    - If IMPINFO is single-valued -> do nothing
    - If multi-valued and external infofile+infocolumn are provided -> do nothing
    - Otherwise compute row-wise median and update sample_column_dict['imp_info_col']
    """

    imp_col = optional_text(sample_column_dict.get("imp_info_col"))
    infofile = optional_text(sample_column_dict.get("infofile"))
    infocolumn = optional_text(sample_column_dict.get("infocolumn"))

    if not imp_col or imp_col not in df.columns:
        if infofile and infocolumn:
            print(screen_field(
                "info", "Imputation quality",
                "study column not provided; external file %s, column %s will be used"
                % (infofile, infocolumn),
                indent=4, label_width=20,
            ))
        return df, sample_column_dict

    # detect whether any row has comma-separated multi-values
    has_multi = df.select(
        pl.col(imp_col)
        .cast(pl.Utf8, strict=False)
        .str.contains(",")
        .any()
    ).item()

    if not has_multi:
        return df, sample_column_dict

    # if external info resource is provided, prefer that
    if infofile and infocolumn:
        print(screen_field(
            "info", "Imputation quality",
            "column %s contains multiple values; external file %s, column %s will be used"
            % (imp_col, infofile, infocolumn),
            indent=4, label_width=20,
        ))
        return df, sample_column_dict

    print(screen_field(
        "analysis", "Imputation quality",
        "column %s contains multiple values; calculating one row-wise median"
        % imp_col,
        indent=4, label_width=20,
    ))
    new_col = f"{imp_col}_median"

    missing_tokens = [
        str(value).strip().upper()
        for value in policies.get("input.null_values")
    ]
    df = df.with_columns(
        pl.col(imp_col)
        .cast(pl.Utf8, strict=False)
        .str.split(",")
        .list.eval(
            pl.when(
                pl.element()
                .cast(pl.Utf8, strict=False)
                .str.strip_chars()
                .str.to_uppercase()
                .is_in(missing_tokens)
            )
            .then(None)
            .otherwise(pl.element().cast(pl.Float64, strict=False))
        )
        .list.drop_nulls()
        .list.median()
        .alias(new_col)
    )

    sample_column_dict["imp_info_col"] = new_col

    missing_count = df.select(pl.col(new_col).is_null().sum()).item()
    print(screen_field(
        "success", "INFO column", "%s created" % new_col,
        indent=4, label_width=20,
    ))
    print(screen_field(
        "warning" if missing_count else "info", "Missing INFO medians",
        "{:,}".format(missing_count), indent=4, label_width=20,
    ))
    return df, sample_column_dict


class AmbiguousColumnMappingError(ValueError):
    """Two config keys name the same column in the input file.

    Raised in place of the polars ``DuplicateError``, which named the *column*
    rather than the config keys that caused it.
    """


def _resolve_required_columns(sample_column_dict, df_columns, policies, logger):
    """The configured columns present in the frame, each exactly once.

    ``required_cols`` used to be built from 13 config keys with no
    de-duplication, so a config with (say) ``beta_or_col`` and ``imp_z_col``
    both set to ``Z`` produced two selection expressions with the same output
    name and died with ``DuplicateError: projections contained duplicate output
    name 'Z'`` — a polars error naming a column, not the config key.
    """
    by_column = {}
    order = []
    for key in COLUMN_CONFIG_KEYS:
        column = optional_text(sample_column_dict.get(key)) or ""
        if not column:
            continue
        if column not in by_column:
            by_column[column] = []
            order.append(column)
        by_column[column].append(key)

    collisions = [(c, by_column[c]) for c in order if len(by_column[c]) > 1]
    for column, keys in collisions:
        _emit(
            logger,
            "Config keys %s name the same column '%s'. It is read once."
            % (" and ".join(keys), column),
            warn=True,
        )
    if collisions and str(policies.get("validation.on_ambiguous_column_mapping")).lower() == "fail":
        raise AmbiguousColumnMappingError(
            "Two or more config keys name the same input column: "
            + "; ".join(
                "%s -> '%s'" % (" and ".join(keys), column) for column, keys in collisions
            )
            + ". Point each key at its own column, or set policy "
            "validation.on_ambiguous_column_mapping to 'warn'."
        )

    return [c for c in order if c in df_columns]


def _missing_data_mask(df, sample_column_dict, policies, logger):
    """Return the mask for variants missing a scientifically required field."""
    # A nested entry is an alternative: the variant needs at least ONE of its
    # members. That is how the effect estimate is expressed as [beta, zscore]:
    # either value is sufficient because the downstream steps derive the other.
    spec = policies.get("columns.mandatory")
    plain, alternatives = [], []
    for entry in spec or []:
        if isinstance(entry, (list, tuple)):
            alternatives.append(list(entry))
        else:
            plain.append(entry)
    mandatory = _fields_to_columns(plain, sample_column_dict, set(df.columns))
    alt_columns = [
        _fields_to_columns(group, sample_column_dict, set(df.columns))
        for group in alternatives
    ]
    alt_columns = [cols for cols in alt_columns if cols]
    if not mandatory and not alt_columns:
        _emit(
            logger,
            "None of the mandatory fields (%s) map to a column in this file; "
            "the input schema validation should have rejected this mapping."
            % ", ".join(str(f) for f in spec),
            warn=True,
        )
        return pl.lit(False), []
    conditions = [pl.col(c).is_null() for c in mandatory]
    for cols in alt_columns:
        conditions.append(pl.all_horizontal([pl.col(c).is_null() for c in cols]))
    used = list(mandatory)
    for cols in alt_columns:
        used.extend(cols)
    return pl.any_horizontal(conditions), used


def _temporary_column(existing, base):
    """Return a private column name that cannot overwrite an input column."""
    candidate = base
    suffix = 1
    while candidate in existing:
        candidate = "%s_%d" % (base, suffix)
        suffix += 1
    existing.add(candidate)
    return candidate


def _numeric_ranking_expression(df, column):
    """Read a ranking value without changing the original study column."""
    expression = pl.col(column)
    if df.schema[column] == pl.String:
        expression = expression.str.replace_all(r"[,\s_]", "")
    return expression.cast(pl.Float64, strict=False)


def _empty_duplicate_report(df):
    """The stable report schema used when no duplicate group exists."""
    return (
        df.head(0)
        .with_columns([
            pl.lit(None, dtype=pl.UInt32).alias("duplicate_input_row"),
            pl.lit(None, dtype=pl.String).alias("duplicate_class"),
            pl.lit(None, dtype=pl.String).alias("duplicate_conflicting_fields"),
            pl.lit(None, dtype=pl.String).alias("duplicate_action"),
        ])
        .select(list(df.columns) + list(DUPLICATE_REPORT_COLUMNS))
    )


def resolve_duplicate_variants(
    df,
    sample_column_dict,
    policies=None,
    logger=None,
    rejects=None,
):
    """Validate duplicate groups and retain one scientifically consistent row.

    Rows are grouped by ``duplicates.key``. Within a group, non-null values in
    every available ``duplicates.consistency_fields`` column must agree exactly.
    Conflicting groups are removed in full. Consistent groups retain one row by
    the configured quality ranking, ending with original input order so the
    result is reproducible. The original row order is restored before return.

    Returns ``(retained_frame, duplicate_report, statistics)``.
    """
    policies = resolve_policies(policies)
    working_columns = list(df.columns)
    report_columns = [
        column for column in df.columns if column != SOURCE_INPUT_ROW_COLUMN
    ]
    collisions = sorted(set(report_columns) & set(DUPLICATE_REPORT_COLUMNS))
    if collisions:
        raise ValueError(
            "The input uses reserved duplicate-report column name(s): %s. Rename "
            "those input columns before harmonisation."
            % ", ".join(collisions)
        )

    key_fields = [str(field).strip().lower() for field in policies.get("duplicates.key")]
    key_pairs = _field_column_pairs(key_fields, sample_column_dict, set(df.columns))
    resolved_key_fields = [field for field, _column in key_pairs]
    unresolved = [field for field in key_fields if field not in resolved_key_fields]
    key_source_columns = [column for _field, column in key_pairs]
    repeated_sources = sorted({
        column for column in key_source_columns
        if key_source_columns.count(column) > 1
    })
    if unresolved or repeated_sources:
        details = []
        if unresolved:
            details.append("unresolved fields: %s" % ", ".join(unresolved))
        if repeated_sources:
            details.append(
                "multiple key fields resolve to the same column: %s"
                % ", ".join(repeated_sources)
            )
        raise ValueError(
            "duplicates.key cannot be applied (%s). Configure every key field to "
            "one distinct input column; duplicate validation is not skipped."
            % "; ".join(details)
        )

    used_names = set(working_columns)
    input_order_col = _temporary_column(used_names, "__postgwas_duplicate_input_order")
    key_columns = []
    key_expressions = []
    for index, (field, column) in enumerate(key_pairs):
        key_column = _temporary_column(
            used_names, "__postgwas_duplicate_key_%d" % index,
        )
        expression = pl.col(column)
        if field in {"chr", "ea", "oa"}:
            expression = (
                expression.cast(pl.String, strict=False)
                .str.strip_chars()
                .str.to_uppercase()
            )
        key_columns.append(key_column)
        key_expressions.append(expression.alias(key_column))

    duplicate_col = _temporary_column(used_names, "__postgwas_is_duplicate")
    working = (
        df.with_row_index(input_order_col)
        .with_columns(key_expressions)
        .with_columns(
            pl.struct(key_columns).is_duplicated().alias(duplicate_col)
        )
    )
    duplicate_rows = working.filter(pl.col(duplicate_col))
    if duplicate_rows.is_empty():
        statistics = {
            "duplicate_groups": 0,
            "duplicate_rows": 0,
            "consistent_groups": 0,
            "conflicting_groups": 0,
            "consistent_rows_removed": 0,
            "conflicting_rows_removed": 0,
            "rows_removed": 0,
        }
        return (
            df.select(working_columns),
            _empty_duplicate_report(df.select(report_columns)),
            statistics,
        )

    consistency_fields = [
        str(field).strip().lower()
        for field in policies.get("duplicates.consistency_fields")
    ]
    consistency_pairs = _field_column_pairs(
        consistency_fields, sample_column_dict, set(df.columns),
    )
    consistency_columns = []
    resolved_consistency_fields = []
    for field, column in consistency_pairs:
        if column not in consistency_columns:
            resolved_consistency_fields.append(field)
            consistency_columns.append(column)
    if not consistency_columns:
        raise ValueError(
            "None of duplicates.consistency_fields resolves to an available input "
            "column. Configure at least one scientific field before duplicate "
            "validation."
        )

    conflict_columns = []
    conflict_expressions = []
    for index, (field, column) in enumerate(zip(
        resolved_consistency_fields, consistency_columns,
    )):
        conflict_column = _temporary_column(
            used_names, "__postgwas_duplicate_conflict_%d" % index,
        )
        conflict_columns.append(conflict_column)
        conflict_expressions.append(
            (
                pl.col(column).drop_nulls().n_unique().over(key_columns) > 1
            ).alias(conflict_column)
        )

    conflict_col = _temporary_column(used_names, "__postgwas_duplicate_conflicting")
    duplicate_rows = (
        duplicate_rows.with_columns(conflict_expressions)
        .with_columns(
            pl.any_horizontal([pl.col(column) for column in conflict_columns])
            .alias(conflict_col)
        )
    )

    quality_columns = _fields_to_columns(
        policies.get("duplicates.quality_fields"),
        sample_column_dict,
        set(df.columns),
    )
    selection_order = list(policies.get("duplicates.selection_order"))
    if "completeness" in selection_order and not quality_columns:
        raise ValueError(
            "duplicates.selection_order uses completeness, but none of "
            "duplicates.quality_fields resolves to an available input column."
        )
    completeness_col = _temporary_column(
        used_names, "__postgwas_duplicate_completeness",
    )
    completeness_expression = pl.sum_horizontal([
        pl.col(column).is_not_null().cast(pl.Int64)
        for column in quality_columns
    ])
    duplicate_rows = duplicate_rows.with_columns(
        completeness_expression.alias(completeness_col)
    )

    sample_size_col = _temporary_column(
        used_names, "__postgwas_duplicate_sample_size",
    )
    sample_size_columns = _fields_to_columns(
        ["ncontrol", "ncase"], sample_column_dict, set(df.columns),
    )
    if len(sample_size_columns) == 2:
        control_column, case_column = sample_size_columns
        sample_size_expression = effective_sample_size_expression(
            _numeric_ranking_expression(df, case_column),
            _numeric_ranking_expression(df, control_column),
        )
    elif sample_size_columns:
        sample_size_expression = _numeric_ranking_expression(
            df, sample_size_columns[0],
        )
    else:
        sample_size_expression = pl.lit(None, dtype=pl.Float64)
    duplicate_rows = duplicate_rows.with_columns(
        sample_size_expression.alias(sample_size_col)
    )

    info_rank_col = _temporary_column(used_names, "__postgwas_duplicate_info")
    info_columns = _fields_to_columns(
        ["info"], sample_column_dict, set(df.columns),
    )
    info_expression = (
        _numeric_ranking_expression(df, info_columns[0])
        if info_columns
        else pl.lit(None, dtype=pl.Float64)
    )
    duplicate_rows = duplicate_rows.with_columns(
        info_expression.alias(info_rank_col)
    )

    ranking_columns = {
        "completeness": (completeness_col, True),
        "sample_size": (sample_size_col, True),
        "info": (info_rank_col, True),
        "input_order": (input_order_col, False),
    }
    sort_columns = list(key_columns)
    descending = [False] * len(key_columns)
    for criterion in selection_order:
        column, high_first = ranking_columns[str(criterion)]
        sort_columns.append(column)
        descending.append(high_first)

    group_head_col = _temporary_column(used_names, "__postgwas_duplicate_group_head")
    ordered = (
        duplicate_rows.sort(
            sort_columns, descending=descending, nulls_last=True,
        )
        .with_columns(
            pl.struct(key_columns).is_first_distinct().alias(group_head_col)
        )
    )

    conflict_fields_expression = pl.concat_str([
        pl.when(pl.col(conflict_column))
        .then(pl.lit(field))
        .otherwise(None)
        for field, conflict_column in zip(
            resolved_consistency_fields, conflict_columns,
        )
    ], separator=",", ignore_nulls=True)
    conflict_action = str(policies.get("duplicates.conflicting_action"))
    report = (
        ordered.with_columns([
            (
                pl.col(SOURCE_INPUT_ROW_COLUMN)
                if SOURCE_INPUT_ROW_COLUMN in ordered.columns
                else pl.col(input_order_col) + 1
            ).alias("duplicate_input_row"),
            pl.when(pl.col(conflict_col))
            .then(pl.lit("conflicting"))
            .otherwise(pl.lit("consistent"))
            .alias("duplicate_class"),
            conflict_fields_expression.alias("duplicate_conflicting_fields"),
            pl.when(pl.col(conflict_col))
            .then(pl.lit(conflict_action))
            .when(pl.col(group_head_col))
            .then(pl.lit("kept"))
            .otherwise(pl.lit("removed"))
            .alias("duplicate_action"),
        ])
        .sort(input_order_col)
        .select(report_columns + list(DUPLICATE_REPORT_COLUMNS))
    )

    group_heads = ordered.filter(pl.col(group_head_col))
    conflicting_groups = group_heads.filter(pl.col(conflict_col)).height
    consistent_groups = group_heads.height - conflicting_groups

    retained_duplicates, conflicting_removed = reject_rows(
        ordered,
        pl.col(conflict_col),
        step_label="02 fix_chr_pos_allele",
        reason="conflicting_duplicate",
        collector=rejects,
        detail="non-empty values disagree in one or more of %s"
        % ",".join(resolved_consistency_fields),
    )
    retained_duplicates, consistent_removed = reject_rows(
        retained_duplicates,
        ~pl.col(group_head_col),
        step_label="02 fix_chr_pos_allele",
        reason="duplicate_variant",
        collector=rejects,
        detail="scientific values agree; selected by %s"
        % ",".join(selection_order),
    )

    selected_input_rows = retained_duplicates.get_column(input_order_col).implode()
    retained = working.filter(
        ~pl.col(duplicate_col)
        | pl.col(input_order_col).is_in(selected_input_rows)
    ).select(
        working_columns
    )
    statistics = {
        "duplicate_groups": group_heads.height,
        "duplicate_rows": ordered.height,
        "consistent_groups": consistent_groups,
        "conflicting_groups": conflicting_groups,
        "consistent_rows_removed": consistent_removed,
        "conflicting_rows_removed": conflicting_removed,
        "rows_removed": consistent_removed + conflicting_removed,
    }
    _emit(
        logger,
        "Duplicate validation: %d groups (%d consistent, %d conflicting); "
        "%d rows removed. Consistency fields: %s. Selection order: %s."
        % (
            statistics["duplicate_groups"],
            consistent_groups,
            conflicting_groups,
            statistics["rows_removed"],
            ", ".join(resolved_consistency_fields),
            ", ".join(selection_order),
        ),
        warn=bool(conflicting_groups),
    )
    return retained, report, statistics


def read_summary_statistics(
    sumstat_file: str,
    output_dir: str,
    sample_column_dict: dict,
    output_layout: dict[str, str],
    table_delimiter: str,
    policies=None,
    logger=None,
    rejects=None,
    input_line_count=None,
):
    """Read the study file and apply the input-stage QC.

    Returns the retained frame, its immutable pre-harmonisation source snapshot,
    row counts, the resolved column mapping, and input-stage QC counts needed by
    the orchestrator. Removed variants are recorded through the supplied reject
    collector.
    """
    policies = resolve_policies(policies)
    os.makedirs(output_dir, exist_ok=True)

    if input_line_count is None:
        input_line_count = count_data_lines(sumstat_file, skip_hash=True) - 1
    df, sample_column_dict = load_summary_statistics_table(
        sumstat_file, sample_column_dict, policies=policies, logger=logger
    )
    python_read_count = df.height

    removed_coords = 0
    non_standard_allele_count = 0
    removed_duplicates_count = 0
    removed_missing_count = 0

    gwas_name = sample_column_dict["gwas_outputname"]

    # Capture the parsed study values before any scientific normalization or
    # transformation. Polars frames are immutable, so clone() is a cheap,
    # stable snapshot until it is partitioned for chromosome workers.
    original_columns = list(df.columns)
    if SOURCE_INPUT_ROW_COLUMN in original_columns:
        raise ValueError(
            "The input uses reserved internal column name %r. Rename that column "
            "before harmonisation." % SOURCE_INPUT_ROW_COLUMN
        )
    df = df.with_row_index(SOURCE_INPUT_ROW_COLUMN, offset=1)
    source_snapshot = df.clone()
    own_collector = False
    if rejects is None and bool(policies.get("rejects.enabled")):
        reject_dir = configured_output_path(
            output_dir, output_layout["rejected_directory"],
        )
        reject_dir.mkdir(parents=True, exist_ok=True)
        rejects = RejectCollector(
            source_snapshot=source_snapshot,
            logger=logger,
            out_path=str(configured_output_path(
                output_dir,
                output_layout["input_reject"],
                dataset_id=gwas_name,
            )),
            delimiter=table_delimiter,
            compress=bool(policies.get("rejects.compress")),
        )
        own_collector = True

    try:
        _emit(
            logger,
            "Sumstat reading successful: rows loaded %d, columns %s"
            % (df.height, df.columns),
        )

        # =========================================================
        # EARLY MISSINGNESS FILTER
        # =========================================================
        required_cols = _resolve_required_columns(
            sample_column_dict, set(df.columns), policies, logger
        )

        if required_cols:
            total_rows = df.height

            missing_stats = df.select([
                pl.col(c).is_null().sum().alias(c) for c in required_cols
            ])

            _emit(logger, "Missing values per configured column (early QC):")
            for c in required_cols:
                missing_count = missing_stats[c][0]
                missing_pct = (missing_count / total_rows) * 100 if total_rows else 0.0
                _emit(logger, "   %s: %d (%.2f%%)" % (c, missing_count, missing_pct))

            before_rows = total_rows
            remove_expr, rule_cols = _missing_data_mask(
                df, sample_column_dict, policies, logger
            )
            remove_mask = _materialise_mask(df, remove_expr)

            df, _ = reject_rows(
                df,
                remove_mask,
                step_label="01 read_summary_statistics",
                reason="missing_required_columns",
                collector=rejects,
                detail="missing one or more mandatory fields across %d column(s)"
                % len(rule_cols),
            )

            after_rows = df.height
            removed_missing_count = before_rows - after_rows

            _emit(logger, "Mandatory-field validation:")
            _emit(logger, "   Before : %d" % before_rows)
            _emit(logger, "   After  : %d" % after_rows)
            _emit(logger, "   Removed: %d" % removed_missing_count)

        # =========================================================
        # COLUMN FIX
        # =========================================================
        # policies / logger / rejects are forwarded deliberately.  Without them
        # this sub-step falls back to a plain df.filter(), so variants dropped
        # here leave no reject-file row and no QC line - the exact silent loss
        # the reject collector exists to eliminate.  Forwarding also makes the
        # dataset-level read answer chromosome.rename_map, chromosome.drop_mt,
        # chromosome.allowed_after_split, position.extraction and
        # position.min_value the SAME way the per-chromosome steps do; with the
        # registry defaults used here and the user's values used there, one run
        # would otherwise resolve the same policy two different ways.
        df, sample_column_dict = harmonise_coordinates_and_alleles(
            chromosome="All_Chrs",
            df=df,
            sample_column_dict=sample_column_dict,
            drop_mt=None,
            policies=policies,
            logger=logger,
            rejects=rejects,
        )

        chr_col = sample_column_dict.get("chr_col")
        pos_col = sample_column_dict.get("pos_col")
        ea_col = sample_column_dict.get("ea_col")
        oa_col = sample_column_dict.get("oa_col")

        numeric_columns = [
            column
            for key in NUMERIC_COLUMN_CONFIG_KEYS
            if (column := optional_text(sample_column_dict.get(key))) is not None
        ]

        canonical_columns = [
            column for column in (chr_col, pos_col, ea_col, oa_col) if column
        ]
        df = normalise_string_columns(
            df,
            preserved_columns=canonical_columns,
        )
        df = normalise_summary_statistics_values(
            df,
            chr_col,
            policies=policies,
            logger=logger,
            strings_normalised=True,
            preserved_columns=canonical_columns,
            numeric_columns=numeric_columns,
        )

        _emit(logger, "Recovery successful: %d rows loaded." % df.height)

        # =========================================================
        # CHR FIX
        # =========================================================
        valid_chr = [str(c).upper() for c in policies.get("chromosome.allowed")]

        if chr_col in df.columns:
            before_chr = df.height

            # Same rows removed as before, but each one now carries the reason
            # that explains it rather than a single "invalid coordinate" count.
            df, _ = reject_rows(
                df,
                ~pl.col(chr_col).is_in(valid_chr),
                step_label="02 fix_chr_pos_allele",
                reason="invalid_chromosome",
                collector=rejects,
                detail="not one of %s" % ",".join(valid_chr),
            )
            if pos_col in df.columns:
                df, _ = reject_rows(
                    df,
                    pl.col(pos_col).is_null(),
                    step_label="02 fix_chr_pos_allele",
                    reason="invalid_position",
                    collector=rejects,
                )

            removed_coords = before_chr - df.height
            _emit(logger, "Removed %d invalid coordinate rows" % removed_coords)

        # =========================================================
        # ALLELE FILTER
        # =========================================================
        if ea_col in df.columns and oa_col in df.columns:
            pattern = str(policies.get("allele.pattern"))

            # str.contains returns NULL for a null input, and `~null` is null,
            # which filter treats as false — so a variant with a null allele
            # used to be excluded by BOTH the keep filter and the reject filter
            # and vanish without trace.  fill_null makes the two halves exact
            # complements, and null alleles get their own reason.
            valid_expr = (
                pl.col(ea_col).str.contains(pattern)
                & pl.col(oa_col).str.contains(pattern)
            ).fill_null(False)
            null_allele_expr = (
                pl.col(ea_col).is_null() | pl.col(oa_col).is_null()
            ).fill_null(True)

            before_alleles = df.height
            df, _ = reject_rows(
                df,
                null_allele_expr,
                step_label="02 fix_chr_pos_allele",
                reason="null_allele",
                collector=rejects,
            )
            df, _ = reject_rows(
                df,
                ~valid_expr,
                step_label="02 fix_chr_pos_allele",
                reason="non_standard_allele",
                collector=rejects,
                detail="allele does not match %s" % pattern,
            )
            non_standard_allele_count = before_alleles - df.height

        _emit(logger, "Removed %d invalid allele rows" % non_standard_allele_count)

        # This private contract is set only after both coordinate and allele
        # filtering have completed. Content validation can then reuse the
        # result without repeating a full-frame scientific validity scan.
        mark_canonical_variant_columns(sample_column_dict)

        # =========================================================
        # DUPLICATES
        # =========================================================
        df, dup_df, duplicate_statistics = resolve_duplicate_variants(
            df,
            sample_column_dict,
            policies=policies,
            logger=logger,
            rejects=rejects,
        )
        removed_duplicates_count = duplicate_statistics["rows_removed"]
        dup_file = configured_output_path(
            output_dir,
            output_layout["duplicates"],
            dataset_id=gwas_name,
        )
        if dup_df.height > 0:
            dup_df.write_csv(dup_file, separator=table_delimiter)
            _emit(logger, "Duplicate assessment saved to %s" % dup_file)

        if (
            duplicate_statistics["conflicting_groups"]
            and policies.get("duplicates.conflicting_action") == "fail_dataset"
        ):
            raise ValueError(
                "%d conflicting duplicate group(s) were found. Their rows were "
                "recorded as conflicting_duplicate in %s and the dataset stopped "
                "because duplicates.conflicting_action is 'fail_dataset'."
                % (duplicate_statistics["conflicting_groups"], dup_file)
            )

        _emit(logger, "Removed %d duplicate variants" % removed_duplicates_count)

        # =========================================================
        # IMPINFO PROCESSING
        # =========================================================
        df, sample_column_dict = normalise_imputation_quality_column(
            df, sample_column_dict, policies,
        )

        if abs(input_line_count - python_read_count) > 1:
            _emit(
                logger,
                "Input line count (%d) and Polars row count (%d) disagree."
                % (input_line_count, python_read_count),
                warn=True,
            )

        if rejects is None and SOURCE_INPUT_ROW_COLUMN in df.columns:
            # Provenance is disabled and no external collector was supplied,
            # so do not carry an unused internal column through every
            # chromosome file and scientific transformation.
            df = df.drop(SOURCE_INPUT_ROW_COLUMN)

        return (
            df,
            source_snapshot,
            input_line_count,
            python_read_count,
            sample_column_dict,
            removed_coords,
            non_standard_allele_count,
            removed_duplicates_count,
            removed_missing_count,
        )
    finally:
        # Flushed from a finally, so a failure half way through still leaves
        # everything rejected up to that point on disk.  Do not replace an
        # analysis exception with a secondary provenance error, but a write
        # failure after otherwise successful input QC must fail the dataset.
        if own_collector and rejects is not None:
            analysis_failed = sys.exc_info()[0] is not None
            try:
                rejects.flush()
            except Exception:
                if not analysis_failed:
                    raise


def _inspect_text_file(file_path):
    """Validate one text stream while retaining its last line and record count."""
    last_line = ""
    record_count = 0
    with open_text(file_path) as handle:
        for line in handle:
            last_line = line
            if not line.startswith("##"):
                record_count += 1
    return last_line.strip(), record_count


def inspect_summary_statistics_file(file_path, policies=None, logger=None):
    """Inspect the input stream once for truncation and physical data-line count.

    The count follows :func:`count_data_lines` exactly: double-hash metadata
    records are excluded and one header record is subtracted.  Returning that
    count lets the normal pipeline avoid decompressing the complete study file
    a second time immediately before the Polars parse.
    """
    policies = resolve_policies(policies)
    source = Path(file_path)
    if not source.is_file():
        raise FileNotFoundError("Input file not found: %s" % source)
    if source.stat().st_size == 0:
        raise ValueError("Input file is empty: %s" % source)

    warnings = []
    record_count = None
    try:
        last_line, record_count = _inspect_text_file(source)
    except (OSError, EOFError, ValueError, zipfile.BadZipFile) as exc:
        if bool(policies.get("input.check_truncation")):
            warnings.append("The file could not be read completely: %s" % exc)
    else:
        if bool(policies.get("input.check_truncation")):
            if not last_line:
                warnings.append("The file has no readable data line")
            elif last_line.endswith("."):
                warnings.append("The last line ends with '.', so the file may be truncated")
        else:
            _emit(logger, "Truncation check skipped: policy input.check_truncation is off.")

    if warnings:
        _emit(logger, "Possible file truncation detected:", warn=True)
        for warning in warnings:
            _emit(logger, "   - %s" % warning, warn=True)

    data_line_count = None if record_count is None else record_count - 1
    return warnings, data_line_count


def check_file_truncation(file_path, policies=None, logger=None):
    """Compatibility wrapper returning only truncation warnings."""
    warnings, _line_count = inspect_summary_statistics_file(
        file_path, policies=policies, logger=logger,
    )
    return warnings


def resolve_resource_file(
    input_file, grch_version, chromosome, *, must_exist=True,
):
    """Resolve an optional user-file template and optionally verify the path."""
    input_file = optional_text(input_file)
    if input_file is None:
        return None
    try:
        resolved = str(input_file).format(
            build=grch_version, chromosome=chromosome
        )
    except (KeyError, ValueError) as exc:
        raise ValueError("Invalid external-file path template %r: %s" % (input_file, exc))
    if must_exist and not Path(resolved).is_file():
        raise FileNotFoundError(
            "External file not found for chromosome %s: %s" % (chromosome, resolved)
        )
    return resolved
