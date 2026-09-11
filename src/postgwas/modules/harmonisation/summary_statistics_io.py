import copy
import json
import os
import re
import sys
import zipfile
from functools import partial
from pathlib import Path

import polars as pl

from postgwas.core.paths import configured_output_path
from postgwas.core.io.delimiters import open_binary, open_text, resolve_delimiter

from postgwas.core.values import missing_tokens, optional_text

from .coordinates import harmonise_coordinates_and_alleles
from .policies import FIELD_LIFECYCLE
from .p_values import preserve_pvalue_source_text
from .rejects import (
    RejectCollector,
    ReconciliationError,
    SOURCE_INPUT_ROW_COLUMN,
)
from .sample_size import (
    effective_sample_size_expression,
    sample_count_expression,
)
from .shared.runtime import emit_message, reject_rows, resolve_policies
from .shared.variant_columns import (
    mark_canonical_variant_columns,
    minimal_allele_representation_series,
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
    null_values = list(missing_tokens(
        policies.get("input.null_values"),
        include_standard=False,
        exact=True,
    ))
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
    pvalue_column = optional_text(sample_column_dict.get("pval_col"))
    if pvalue_column:
        # The source token is required to distinguish a literal zero from a
        # positive probability below Float64 (for example 1e-400).
        string_schema[pvalue_column] = pl.String
    if string_schema:
        _emit(
            logger,
            "Variant identity and p-value source columns read as text to "
            "preserve their exact tokens: %s."
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
    sample_count_columns=None,
    compound_numeric_columns=None,
    conversion_counts=None,
) -> pl.DataFrame:
    """
    Normalize configured scientific numeric columns after the first read.

    Unconfigured columns retain the type inferred by Polars. Variant identity
    columns are forced to strings by :func:`load_summary_statistics_table`, so
    numeric-looking IDs never pass through Float64.

    Scalar numeric values that cannot be parsed become null and are reported;
    the immutable source snapshot retains their original text for rejection
    provenance. Sample-count columns reuse their whole-number parser, while
    configured compound numeric columns such as delimited INFO lists remain
    text for their field-specific parser.

    When supplied, conversion_counts receives existing per-column conversion
    counters without another data scan or any change to the returned values.
    """
    if df.height == 0:
        return df
    policies = resolve_policies(policies)
    preserved = set(preserved_columns or ())
    numeric = set(numeric_columns or ())
    sample_counts = set(sample_count_columns or ())
    compound_numeric = set(compound_numeric_columns or ())
    compound_delimiter = str(policies.get("info.multi_value_delimiter"))
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
        original = df.get_column(col)
        if (
            col in compound_numeric
            and original.dtype == pl.String
            and bool(
                original.str.contains(
                    compound_delimiter, literal=True,
                ).any()
            )
        ):
            if conversion_counts is not None:
                conversion_counts[col] = {"status": "compound_values_deferred", "rows": df.height}
            valid_cols.append(pl.col(col))
            continue
        numeric_series = (
            df.select(sample_count_expression(df, col).alias(col)).get_column(col)
            if col in sample_counts
            else original.cast(pl.Float64, strict=False)
        )
        new_nulls = numeric_series.null_count() - original.null_count()
        if conversion_counts is not None:
            conversion_counts[col] = {
                "status": "assessed", "rows": df.height,
                "non_missing": original.len() - original.null_count(),
                "invalid": new_nulls,
            }
        valid_cols.append(numeric_series)
        if new_nulls > 0:
            non_null = original.len() - original.null_count()
            _emit(
                logger,
                "Configured numeric column '%s' had %d of %d non-missing values "
                "that could not be parsed; those cells were set to missing. "
                "Original input values remain available in rejection provenance."
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
    policies=None,
    logger=None,
):
    """Reduce a delimited study INFO list under the resolved YAML policy.

    A scalar study INFO column is returned unchanged. When at least one row
    contains multiple values separated by ``info.multi_value_delimiter``, every
    row is parsed to the
    collision-checked working column configured by
    ``info.multi_value_output_column``. Finite non-missing values are reduced
    by the unweighted row-wise median or mean, or the dataset stops when the
    policy is ``fail``. The original input text remains in the frame and in the
    immutable source snapshot; the column mapping is updated so all later
    scientific operations consume the scalar value.

    Sample-size weighting is deliberately unavailable here: the input contract
    has no ordered cohort-level sample-size vector that can be paired safely
    with the ordered INFO values.
    """
    policies = resolve_policies(policies)

    imp_col = optional_text(sample_column_dict.get("imp_info_col"))
    infofile = optional_text(sample_column_dict.get("infofile"))
    infocolumn = optional_text(sample_column_dict.get("infocolumn"))
    info_source = optional_text(sample_column_dict.get("info_source"))

    if not imp_col or imp_col not in df.columns:
        if infofile and infocolumn:
            _emit(
                logger,
                "Study INFO column was not selected; external file %s, column "
                "%s will be used." % (infofile, infocolumn),
            )
        return df, sample_column_dict

    # The validated sample-sheet contract gives the internal source priority
    # and clears infofile/infocolumn in that case. Direct engine callers may
    # explicitly select the external source while retaining an unused study
    # column, so source selection—not mere file presence—controls this branch.
    if info_source == "external":
        _emit(
            logger,
            "Study INFO column %s was not aggregated because the resolved "
            "INFO source is the external file %s, column %s."
            % (imp_col, infofile or "not recorded", infocolumn or "not recorded"),
        )
        return df, sample_column_dict

    multi_value_delimiter = str(policies.get("info.multi_value_delimiter"))
    multi_mask = (
        pl.col(imp_col)
        .cast(pl.String, strict=False)
        .str.contains(multi_value_delimiter, literal=True)
        .fill_null(False)
    )
    multi_value_rows = int(df.select(multi_mask.sum()).item() or 0)

    if not multi_value_rows:
        return df, sample_column_dict

    aggregation = str(policies.get("info.multi_value_aggregation"))
    if aggregation == "fail":
        raise ValueError(
            "Study INFO column %r contains values separated by %r in %d row(s), "
            "and info.multi_value_aggregation is 'fail'. Supply one scalar INFO "
            "value per variant, or explicitly select the unweighted 'median' or "
            "'mean' policy in the run configuration."
            % (imp_col, multi_value_delimiter, multi_value_rows)
        )

    output_col = str(policies.get("info.multi_value_output_column"))
    if output_col in df.columns:
        raise ValueError(
            "Configured INFO aggregation output column %r already exists in the "
            "study table. Rename the input column or change "
            "info.multi_value_output_column; PostGWAS will not overwrite it."
            % output_col
        )

    normalized_missing_tokens = missing_tokens(
        policies.get("input.null_values"),
        include_standard=True,
    )
    used_names = set(df.columns)
    invalid_count_col = _temporary_column(
        used_names, "__postgwas_info_invalid_token_count",
    )
    token = pl.element().cast(pl.String, strict=False).str.strip_chars()
    numeric = token.cast(pl.Float64, strict=False)
    token_is_missing = token.str.to_uppercase().is_in(
        normalized_missing_tokens
    )
    numeric_is_finite = numeric.is_finite().fill_null(False)
    token_lists = (
        pl.col(imp_col)
        .cast(pl.String, strict=False)
        .str.split(multi_value_delimiter)
    )
    finite_values = (
        token_lists.list.eval(
            pl.when(token_is_missing)
            .then(None)
            .when(numeric_is_finite)
            .then(numeric)
            .otherwise(None)
        )
        .list.drop_nulls()
    )
    aggregate_expression = (
        finite_values.list.median()
        if aggregation == "median"
        else finite_values.list.mean()
    )
    df = df.with_columns([
        aggregate_expression.cast(pl.Float64).alias(output_col),
        token_lists.list.eval(
            pl.when(token_is_missing | numeric_is_finite)
            .then(pl.lit(0, dtype=pl.UInt32))
            .otherwise(pl.lit(1, dtype=pl.UInt32))
        )
        .list.sum()
        .alias(invalid_count_col),
    ])

    statistics = df.select([
        pl.col(invalid_count_col).fill_null(0).sum().alias("invalid_tokens"),
        (pl.col(invalid_count_col).fill_null(0) > 0)
        .sum()
        .alias("rows_with_invalid_tokens"),
        pl.col(output_col).is_null().sum().alias("missing_aggregates"),
    ]).row(0, named=True)
    invalid_tokens = int(statistics["invalid_tokens"] or 0)
    invalid_rows = int(statistics["rows_with_invalid_tokens"] or 0)
    invalid_token_action = str(
        policies.get("info.multi_value_invalid_token_action")
    )
    if invalid_tokens and invalid_token_action == "fail":
        raise ValueError(
            "Study INFO column %r contains %d non-missing token(s) that are not "
            "finite numbers across %d row(s), and "
            "info.multi_value_invalid_token_action is 'fail'. Correct those "
            "tokens or explicitly select 'ignore' in the run configuration."
            % (imp_col, invalid_tokens, invalid_rows)
        )
    df = df.drop(invalid_count_col)

    sample_column_dict["info_multi_value_source_column"] = imp_col
    sample_column_dict["info_multi_value_aggregation"] = aggregation
    sample_column_dict["info_multi_value_output_column"] = output_col
    sample_column_dict["imp_info_col"] = output_col

    _emit(
        logger,
        "INFO aggregation decision: %d row(s) in source column %r contained "
        "values separated by %r; info.multi_value_aggregation=%r created "
        "scalar column %r before duplicate validation. The aggregation is "
        "unweighted because no aligned cohort-level sample-size vector was "
        "provided."
        % (
            multi_value_rows,
            imp_col,
            multi_value_delimiter,
            aggregation,
            output_col,
        ),
    )
    if invalid_tokens:
        _emit(
            logger,
            "INFO aggregation applied "
            "info.multi_value_invalid_token_action='ignore' to %d non-missing "
            "token(s) that were not "
            "finite numbers across %d row(s); valid finite tokens in those rows "
            "were still aggregated."
            % (invalid_tokens, invalid_rows),
            warn=True,
        )
    missing_aggregates = int(statistics["missing_aggregates"] or 0)
    _emit(
        logger,
        "INFO aggregation produced %d missing scalar value(s) among %d study "
        "row(s); chromosome INFO missing-value policy will handle them."
        % (missing_aggregates, df.height),
        warn=bool(missing_aggregates),
    )
    return df, sample_column_dict


class AmbiguousColumnMappingError(ValueError):
    """Two config keys name the same column in the input file.

    Raised in place of the polars ``DuplicateError``, which named the *column*
    rather than the config keys that caused it.
    """


def _resolve_mapped_columns(sample_column_dict, df_columns, policies, logger):
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


def _effective_field_requirements(policies):
    """Return effective read requirements from ``columns.mandatory``.

    Plain fields are independently mandatory. Nested lists are alternative
    groups for which at least one member must be present on each row.
    """
    plain = []
    alternatives = []
    for entry in policies.get("columns.mandatory") or []:
        if isinstance(entry, (list, tuple)):
            alternatives.append([
                str(field).strip().lower() for field in entry
            ])
        else:
            plain.append(str(entry).strip().lower())
    return plain, alternatives


def _selected_required_input_alternative(field, sample_column_dict, policies):
    """Find the configured source alternative associated with one field.

    The canonical registry owns the input alternatives. This function searches
    that resolved registry rather than encoding special cases for EAF, INFO,
    effect estimates, or sample size in reporting code.
    """
    config_key = FIELD_TO_CONFIG_KEY.get(field)
    if config_key is None:
        return None, []
    for group, specification in policies.required_inputs.items():
        alternatives = specification.get("any_of") or []
        if not any(config_key in alternative for alternative in alternatives):
            continue
        for alternative in alternatives:
            if all(
                optional_text(sample_column_dict.get(key)) is not None
                for key in alternative
            ):
                return str(group), [str(key) for key in alternative]
    return None, []


def _input_field_completeness(df, sample_column_dict, policies):
    """Build the structured post-parse completeness and recovery plan.

    These counts are intentionally measured immediately after parsing and
    before numeric normalization. A later unparseable numeric token is reported
    by numeric normalization and is not mislabelled here as an input null.
    """
    plain, alternatives = _effective_field_requirements(policies)
    alternative_by_field = {}
    for group in alternatives:
        for field in group:
            alternative_by_field[field] = list(group)

    configured_columns = {}
    for field in FIELD_LIFECYCLE:
        config_key = FIELD_TO_CONFIG_KEY.get(field)
        column = optional_text(sample_column_dict.get(config_key)) if config_key else None
        if column is not None and column in df.columns:
            configured_columns[field] = column

    unique_columns = list(dict.fromkeys(
        column
        for key in COLUMN_CONFIG_KEYS
        if (
            (column := optional_text(sample_column_dict.get(key))) is not None
            and column in df.columns
        )
    ))
    missing_by_column = {}
    if unique_columns:
        counts = df.select([
            pl.col(column).is_null().sum().alias("__missing_%d" % index)
            for index, column in enumerate(unique_columns)
        ]).row(0)
        missing_by_column = dict(zip(unique_columns, (int(value) for value in counts)))

    effective_fields = list(FIELD_LIFECYCLE)
    for field in plain + [member for group in alternatives for member in group]:
        if field not in effective_fields:
            effective_fields.append(field)
    for field in policies.get("final_check.require") or []:
        field = str(field).strip().lower()
        if field not in effective_fields:
            effective_fields.append(field)

    final_order = {
        str(field).strip().lower(): index
        for index, field in enumerate(
            policies.get("final_check.require") or [],
            start=1,
        )
    }
    final_required = set(final_order)
    fields = []
    for field in effective_fields:
        lifecycle = dict(FIELD_LIFECYCLE.get(field) or {})
        config_key = FIELD_TO_CONFIG_KEY.get(field)
        configured_column = (
            optional_text(sample_column_dict.get(config_key))
            if config_key else None
        )
        direct_input_column = configured_columns.get(field)
        source_group, source_keys = _selected_required_input_alternative(
            field, sample_column_dict, policies,
        )
        source_input_columns = [
            str(sample_column_dict[key])
            for key in source_keys
            if (
                key in COLUMN_CONFIG_KEYS
                and optional_text(sample_column_dict.get(key)) in df.columns
            )
        ]
        input_column = direct_input_column
        if input_column is None and source_group == "coordinates" and len(source_input_columns) == 1:
            input_column = source_input_columns[0]
        input_missing = (
            missing_by_column[input_column]
            if input_column in missing_by_column else None
        )
        recovery = optional_text(lifecycle.get("recovered"))

        if direct_input_column is not None:
            source_kind = "study_column"
            source_label = "study column %s" % direct_input_column
            source_keys = [config_key] if config_key else []
            source_input_columns = [direct_input_column]
        elif configured_column is not None:
            source_kind = "configured_column_absent"
            source_label = "configured column %s is absent" % configured_column
        elif source_keys:
            source_kind = "configured_alternative"
            source_label = "configured %s" % " + ".join(source_keys)
        elif recovery is not None:
            source_kind = "recovery_step"
            source_label = "not supplied; recover at %s" % recovery
        else:
            source_kind = "not_configured"
            source_label = "not configured"

        alternative_group = alternative_by_field.get(field)
        if (
            field in plain
            and configured_column is None
            and direct_input_column is None
            and source_keys
        ):
            read_requirement = "required_via_configured_alternative"
            read_action = "validate_configured_alternative_during_normalisation"
        elif field in plain:
            read_requirement = "required"
            read_action = "reject_row_if_missing"
        elif alternative_group is not None:
            read_requirement = "alternative"
            read_action = "reject_row_if_all_alternatives_missing"
        elif recovery is not None and direct_input_column is not None:
            read_requirement = "downstream_validated"
            read_action = "retain_for_field_specific_downstream_policy"
        elif recovery is not None:
            read_requirement = "recoverable"
            read_action = "retain_for_configured_recovery_and_downstream_policy"
        else:
            read_requirement = "optional"
            read_action = "not_required_at_read"

        is_final_required = field in final_required
        fields.append({
            "field": field,
            "configured_column": configured_column,
            "input_column": input_column,
            "source_kind": source_kind,
            "source_label": source_label,
            "source_requirement_group": source_group,
            "source_config_keys": source_keys,
            "source_values": {key: sample_column_dict[key] for key in source_keys},
            "source_input_columns": source_input_columns,
            "read_requirement": read_requirement,
            "read_alternative_group": alternative_group or [],
            "read_action": read_action,
            "recovery_step": recovery,
            "input_rows": int(df.height),
            "post_parse_missing": input_missing,
            "post_parse_missing_fraction": (
                float(input_missing) / float(df.height)
                if input_missing is not None and df.height else None
            ),
            "final_required": is_final_required,
            "final_order": final_order.get(field),
            "final_missing": None,
            "final_action": (
                str(policies.get("final_check.on_missing"))
                if is_final_required else "not_checked_by_final_gate"
            ),
            "final_status": "pending" if is_final_required else "not_required",
            "lifecycle_note": optional_text(lifecycle.get("note")),
        })

    alternative_groups = []
    for group in alternatives:
        columns = [
            configured_columns[field]
            for field in group
            if field in configured_columns
        ]
        missing_all = None
        if columns:
            missing_all = int(df.select(
                pl.all_horizontal([
                    pl.col(column).is_null() for column in columns
                ]).sum().alias("missing_all")
            ).item() or 0)
        alternative_groups.append({
            "fields": list(group),
            "input_columns": columns,
            "input_rows": int(df.height),
            "post_parse_missing_all": missing_all,
            "post_parse_missing_all_fraction": (
                float(missing_all) / float(df.height)
                if missing_all is not None and df.height else None
            ),
            "read_action": "reject_row_if_all_alternatives_missing",
        })

    return {
        "input_measurement_stage": "post_parse_before_numeric_normalisation",
        "input_rows": int(df.height),
        "numeric_conversion_by_column": {},
        "info_missing_action": str(policies.get("info.on_missing")),
        "post_parse_missing_by_column": missing_by_column,
        "read_mandatory_policy": copy.deepcopy(
            policies.get("columns.mandatory")
        ),
        "rows_removed_at_read_mandatory_gate": 0,
        "fields": fields,
        "read_alternative_groups": alternative_groups,
        "final_check": {
            "required_fields": list(policies.get("final_check.require")),
            "on_missing": str(policies.get("final_check.on_missing")),
            "measurement": "pending post-recovery assessment",
            "expected_chromosomes": [],
            "completed_chromosomes": [],
            "unassessed_chromosomes": [],
            "rows_entering_gate": None,
            "rows_retained_after_gate": None,
            "status": "pending",
        },
    }


def finalise_field_completeness_report(
    report, chromosome_summaries, completed_chromosomes,
    expected_chromosomes=None,
):
    """Add post-recovery final-gate counts from completed chromosomes."""
    completed = [str(chromosome) for chromosome in completed_chromosomes or []]
    expected = [
        str(chromosome)
        for chromosome in (expected_chromosomes or completed)
    ]
    final = copy.deepcopy(report)
    required = set(final["final_check"]["required_fields"])
    action = str(final["final_check"].get("on_missing") or "")
    if action == "reject":
        measurement = "sequential rule attribution in final_check.require order"
    elif action == "keep":
        measurement = "independent per-field matches; rows may appear in multiple counts"
    else:
        measurement = "fail-fast field assessment in final_check.require order"
    final_qc_by_chromosome = {}
    for chromosome in completed:
        summary = chromosome_summaries.get(chromosome) or {}
        final_qc = (
            (summary.get("stage_qc") or {}).get("final_completeness_qc")
            or {}
        )
        if final_qc:
            final_qc_by_chromosome[chromosome] = final_qc

    for field_report in final["fields"]:
        field = field_report["field"]
        if field not in required:
            continue
        counts = []
        skipped = []
        complete = bool(completed)
        for chromosome in completed:
            final_qc = final_qc_by_chromosome.get(chromosome)
            if final_qc is None:
                complete = False
                break
            missing_by_field = final_qc.get("missing_by_field") or {}
            if field in missing_by_field:
                counts.append(int(missing_by_field[field]))
            elif field in (final_qc.get("fields_skipped") or {}):
                skipped.append(chromosome)
            else:
                complete = False
                break
        chromosome_scope_complete = set(completed) == set(expected) and bool(expected)
        if complete and counts and not skipped:
            field_report["final_missing"] = sum(counts)
            field_report["final_status"] = (
                "assessed" if chromosome_scope_complete
                else "partial_chromosome_coverage"
            )
        elif complete and skipped and not counts:
            field_report["final_missing"] = None
            field_report["final_status"] = "not_assessed_by_missing_value_policy"
        else:
            field_report["final_missing"] = None
            field_report["final_status"] = "incomplete"

    assessed = [
        chromosome for chromosome in completed
        if chromosome in final_qc_by_chromosome
    ]
    unassessed = [
        chromosome for chromosome in expected
        if chromosome not in final_qc_by_chromosome
    ]
    complete_qc = not unassessed and bool(expected)
    final["final_check"].update({
        "measurement": measurement,
        "expected_chromosomes": expected,
        "completed_chromosomes": assessed,
        "unassessed_chromosomes": unassessed,
        "rows_entering_gate": (
            sum(
                int(final_qc_by_chromosome[chromosome].get("initial_variants") or 0)
                for chromosome in assessed
            )
            if assessed else None
        ),
        "rows_retained_after_gate": (
            sum(
                int(final_qc_by_chromosome[chromosome].get("final_variants") or 0)
                for chromosome in assessed
            )
            if assessed else None
        ),
        "status": (
            "complete" if complete_qc
            else "partial" if assessed
            else "incomplete"
        ),
    })
    return final


def write_field_completeness_report(path, report, delimiter):
    """Write the manifest's field-completeness evidence as a flat QC table."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    removed = report.get("rows_removed_at_read_mandatory_gate")
    final_check = report.get("final_check") or {}
    shared = {
        "input_measurement_stage": report.get("input_measurement_stage"),
        "read_gate_rows_rejected_total": removed,
        "final_measurement": final_check.get("measurement"),
        "final_expected_chromosomes": json.dumps(
            final_check.get("expected_chromosomes") or []
        ),
        "final_assessed_chromosomes": json.dumps(
            final_check.get("completed_chromosomes") or []
        ),
        "final_unassessed_chromosomes": json.dumps(
            final_check.get("unassessed_chromosomes") or []
        ),
        "final_report_status": final_check.get("status"),
    }
    for field in report.get("fields") or []:
        rows.append({
            **shared,
            "record_type": "field",
            "field": field["field"],
            "configured_column": field.get("configured_column"),
            "source": field.get("source_label"),
            "source_config_keys": json.dumps(
                field.get("source_config_keys") or []
            ),
            "source_input_columns": json.dumps(
                field.get("source_input_columns") or []
            ),
            "read_requirement": field.get("read_requirement"),
            "read_action": field.get("read_action"),
            "read_alternative_group": json.dumps(
                field.get("read_alternative_group") or []
            ),
            "recovery_step": field.get("recovery_step"),
            "input_rows": field.get("input_rows"),
            "post_parse_missing": field.get("post_parse_missing"),
            "post_parse_missing_fraction": field.get(
                "post_parse_missing_fraction"
            ),
            "numeric_conversion": json.dumps(
                (report.get("numeric_conversion_by_column") or {}).get(field.get("input_column"))
            ),
            "final_required": field.get("final_required"),
            "final_order": field.get("final_order"),
            "final_missing": field.get("final_missing"),
            "final_action": field.get("final_action"),
            "final_status": field.get("final_status"),
            "lifecycle_note": field.get("lifecycle_note"),
        })
    for group in report.get("read_alternative_groups") or []:
        rows.append({
            **shared,
            "record_type": "read_alternative_group",
            "field": json.dumps(group["fields"]),
            "configured_column": json.dumps(group["input_columns"]),
            "source": "at least one alternative is required",
            "source_config_keys": None,
            "source_input_columns": json.dumps(group["input_columns"]),
            "read_requirement": "alternative_group",
            "read_action": group.get("read_action"),
            "read_alternative_group": json.dumps(group["fields"]),
            "recovery_step": None,
            "input_rows": group.get("input_rows"),
            "post_parse_missing": group.get("post_parse_missing_all"),
            "post_parse_missing_fraction": group.get(
                "post_parse_missing_all_fraction"
            ),
            "final_required": None,
            "final_order": None,
            "final_missing": None,
            "final_action": group.get("read_action"),
            "final_status": "read_gate_assessed",
            "lifecycle_note": None,
        })
    pl.DataFrame(rows).write_csv(path, separator=delimiter)
    report["report_path"] = str(path)
    return str(path)


def _missing_data_mask(df, sample_column_dict, policies, logger):
    """Return the mask for variants missing a scientifically required field."""
    # A nested entry is an alternative: the variant needs at least ONE of its
    # members. That is how the effect estimate is expressed as [beta, zscore]:
    # either value is sufficient because the downstream steps derive the other.
    spec = policies.get("columns.mandatory")
    plain, alternatives = _effective_field_requirements(policies)
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
    *,
    step_label="02 fix_chr_pos_allele",
    reference_aligned=False,
):
    """Validate duplicate groups and retain one scientifically consistent row.

    Rows are grouped by ``duplicates.key``. When position, effect allele and
    other allele participate in that key, shared VCF padding is minimally
    trimmed in private key columns before the upper-case lexicographically
    normalized allele pair is formed. Original study coordinates and alleles
    are never changed, and no reference-based left-alignment occurs. Thus A/G
    and G/A, as well as equivalent padded indel representations, identify the
    same physical allele pair. A group that contains both ordered
    representations is removed in full because effect, Z, and frequency values
    cannot be reconciled before their semantics and reference orientation are
    resolved. Among same-orientation rows, non-null values in every available
    ``duplicates.consistency_fields`` column must agree. Dataset-stage values
    use exact equality; the reference-aligned pass permits only the configured
    relative floating-point roundoff. Conflicting groups are removed in full.
    Consistent groups retain one row by the configured quality ranking, ending
    with original input order so the result is reproducible. The original row
    order is restored before return.

    ``reference_aligned=True`` identifies the second, chromosome-level pass:
    strand step 04 has already placed REF in the other-allele column and ALT in
    the effect-allele column, transformed allele-dependent statistics, and
    finished any reference-backed effect-frequency alignment. This lets
    reverse-complement input representations such as A/G and T/C converge on
    one physical variant before the same conservative consistency and ranking
    rules are applied. With the canonical key, both alleles remain key
    components, so distinct alternate alleles at a multiallelic coordinate are
    never collapsed merely because their chromosome and position agree.

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
    key_columns_by_field = dict(key_pairs)
    unordered_allele_key = {"ea", "oa"}.issubset(key_columns_by_field)
    minimally_trimmed_key = {"pos", "ea", "oa"}.issubset(
        key_columns_by_field
    )
    ordered_allele_col = None
    allele_key_expressions = {}
    working = df.with_row_index(input_order_col)
    if unordered_allele_key:
        if minimally_trimmed_key:
            minimal_position_col = _temporary_column(
                used_names, "__postgwas_duplicate_minimal_position",
            )
            minimal_effect_col = _temporary_column(
                used_names, "__postgwas_duplicate_minimal_effect_allele",
            )
            minimal_other_col = _temporary_column(
                used_names, "__postgwas_duplicate_minimal_other_allele",
            )
            working = working.with_columns(
                minimal_allele_representation_series(
                    df,
                    key_columns_by_field["pos"],
                    key_columns_by_field["ea"],
                    key_columns_by_field["oa"],
                    position_output=minimal_position_col,
                    first_allele_output=minimal_effect_col,
                    second_allele_output=minimal_other_col,
                )
            )
            effect_allele = pl.col(minimal_effect_col)
            other_allele = pl.col(minimal_other_col)
        else:
            effect_allele = (
                pl.col(key_columns_by_field["ea"])
                .cast(pl.String, strict=False)
                .str.strip_chars()
                .str.to_uppercase()
            )
            other_allele = (
                pl.col(key_columns_by_field["oa"])
                .cast(pl.String, strict=False)
                .str.strip_chars()
                .str.to_uppercase()
            )
        effect_first = effect_allele <= other_allele
        allele_key_expressions = {
            "ea": pl.when(effect_first).then(effect_allele).otherwise(other_allele),
            "oa": pl.when(effect_first).then(other_allele).otherwise(effect_allele),
        }
        if minimally_trimmed_key:
            allele_key_expressions["pos"] = pl.col(minimal_position_col)
        ordered_allele_col = _temporary_column(
            used_names, "__postgwas_duplicate_ordered_alleles",
        )
        key_expressions.append(
            pl.concat_str([effect_allele, other_allele], separator=">")
            .alias(ordered_allele_col)
        )
    for index, (field, column) in enumerate(key_pairs):
        key_column = _temporary_column(
            used_names, "__postgwas_duplicate_key_%d" % index,
        )
        expression = allele_key_expressions.get(field, pl.col(column))
        if field == "chr" or (field in {"ea", "oa"} and not unordered_allele_key):
            expression = (
                expression.cast(pl.String, strict=False)
                .str.strip_chars()
                .str.to_uppercase()
            )
        key_columns.append(key_column)
        key_expressions.append(expression.alias(key_column))

    duplicate_col = _temporary_column(used_names, "__postgwas_is_duplicate")
    swapped_orientation_col = _temporary_column(
        used_names, "__postgwas_duplicate_swapped_orientation",
    )
    working = (
        working.with_columns(key_expressions)
        .with_columns([
            pl.struct(key_columns).is_duplicated().alias(duplicate_col),
            (
                pl.col(ordered_allele_col).n_unique().over(key_columns) > 1
                if ordered_allele_col is not None
                else pl.lit(False)
            ).alias(swapped_orientation_col),
        ])
    )
    duplicate_rows = working.filter(pl.col(duplicate_col))
    if duplicate_rows.is_empty():
        statistics = {
            "duplicate_groups": 0,
            "duplicate_rows": 0,
            "consistent_groups": 0,
            "conflicting_groups": 0,
            "swapped_orientation_groups": 0,
            "swapped_orientation_rows": 0,
            "consistent_rows_removed": 0,
            "conflicting_rows_removed": 0,
            "swapped_orientation_rows_removed": 0,
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
    post_orientation_relative_tolerance = (
        float(policies.get("duplicates.post_orientation_relative_tolerance"))
        if reference_aligned
        else 0.0
    )
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
        exact_conflict = (
            pl.col(column).drop_nulls().n_unique().over(key_columns) > 1
        )
        if (
            post_orientation_relative_tolerance > 0.0
            and field in {"beta", "se", "zscore", "pval", "eaf", "info"}
        ):
            numeric = pl.col(column).cast(pl.Float64, strict=False)
            non_null_count = numeric.is_not_null().sum().over(key_columns)
            all_finite = numeric.drop_nulls().is_finite().all().over(key_columns)
            group_min = numeric.drop_nulls().min().over(key_columns)
            group_max = numeric.drop_nulls().max().over(key_columns)
            scale = pl.max_horizontal(group_min.abs(), group_max.abs())
            conflict = (
                pl.when(non_null_count <= 1)
                .then(pl.lit(False))
                .when(all_finite)
                .then(
                    (group_max - group_min).abs()
                    > post_orientation_relative_tolerance * scale
                )
                .otherwise(exact_conflict)
            )
        else:
            conflict = exact_conflict
        conflict_expressions.append(conflict.alias(conflict_column))

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
            pl.when(pl.col(swapped_orientation_col))
            .then(pl.lit("swapped_orientation"))
            .when(pl.col(conflict_col))
            .then(pl.lit("conflicting"))
            .otherwise(pl.lit("consistent"))
            .alias("duplicate_class"),
            conflict_fields_expression.alias("duplicate_conflicting_fields"),
            pl.when(pl.col(swapped_orientation_col))
            .then(pl.lit(conflict_action))
            .when(pl.col(conflict_col))
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
    swapped_orientation_groups = group_heads.filter(
        pl.col(swapped_orientation_col)
    ).height
    swapped_orientation_rows = ordered.filter(
        pl.col(swapped_orientation_col)
    ).height
    conflicting_groups = group_heads.filter(
        pl.col(conflict_col) & ~pl.col(swapped_orientation_col)
    ).height
    consistent_groups = (
        group_heads.height - conflicting_groups - swapped_orientation_groups
    )

    retained_duplicates, swapped_orientation_removed = reject_rows(
        ordered,
        pl.col(swapped_orientation_col),
        step_label=step_label,
        reason="swapped_orientation_duplicate",
        collector=rejects,
        detail=(
            "effect and other allele order differs within the reference-aligned "
            "allele-pair group"
            if reference_aligned else
            "effect and other allele order differs within the normalized "
            "allele-pair group"
        ),
    )

    retained_duplicates, conflicting_removed = reject_rows(
        retained_duplicates,
        pl.col(conflict_col),
        step_label=step_label,
        reason="conflicting_duplicate",
        collector=rejects,
        detail=(
            "%snon-empty values disagree in one or more of %s"
            % (
                "reference-aligned " if reference_aligned else "",
                ",".join(resolved_consistency_fields),
            )
        ),
    )
    retained_duplicates, consistent_removed = reject_rows(
        retained_duplicates,
        ~pl.col(group_head_col),
        step_label=step_label,
        reason="duplicate_variant",
        collector=rejects,
        detail=(
            "%sscientific values agree; selected by %s"
            % (
                "reference-aligned " if reference_aligned else "",
                ",".join(selection_order),
            )
        ),
    )

    selected_input_rows = retained_duplicates.get_column(input_order_col).to_list()
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
        "swapped_orientation_groups": swapped_orientation_groups,
        "swapped_orientation_rows": swapped_orientation_rows,
        "consistent_rows_removed": consistent_removed,
        "conflicting_rows_removed": conflicting_removed,
        "swapped_orientation_rows_removed": swapped_orientation_removed,
        "rows_removed": (
            consistent_removed
            + conflicting_removed
            + swapped_orientation_removed
        ),
    }
    _emit(
        logger,
        "%sduplicate validation: %d groups (%d consistent, %d conflicting, "
        "%d swapped orientation); %d rows removed (%d swapped orientation). "
        "Consistency fields: %s. Selection order: %s."
        % (
            "Post-orientation " if reference_aligned else "",
            statistics["duplicate_groups"],
            consistent_groups,
            conflicting_groups,
            swapped_orientation_groups,
            statistics["rows_removed"],
            swapped_orientation_removed,
            ", ".join(resolved_consistency_fields),
            ", ".join(selection_order),
        ),
        warn=bool(conflicting_groups or swapped_orientation_groups),
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
    row counts, the resolved column mapping, and disjoint input-stage QC counts
    for invalid coordinates, unsupported chromosomes, unusable alleles,
    duplicates, and missing mandatory values, followed by the structured
    field-completeness and recovery plan. Removed variants are recorded through
    the supplied reject collector.
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
    removed_unsupported_chromosomes = 0
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
    df, sample_column_dict = preserve_pvalue_source_text(
        df, sample_column_dict
    )
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
    reject_counts_before = rejects.counts() if rejects is not None else {}

    try:
        _emit(
            logger,
            "Sumstat reading successful: rows loaded %d, columns %s"
            % (df.height, df.columns),
        )

        # =========================================================
        # EARLY MISSINGNESS FILTER
        # =========================================================
        mapped_cols = _resolve_mapped_columns(
            sample_column_dict, set(df.columns), policies, logger
        )
        field_completeness = _input_field_completeness(
            df, sample_column_dict, policies,
        )

        if mapped_cols:
            total_rows = df.height
            missing_by_column = field_completeness[
                "post_parse_missing_by_column"
            ]

            _emit(
                logger,
                "Post-parse missing values by configured study column "
                "(before numeric normalization):",
            )
            for c in mapped_cols:
                missing_count = missing_by_column[c]
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
                reason="missing_read_mandatory_value",
                collector=rejects,
                detail=(
                    "missing a read-stage mandatory field, or every member of "
                    "one mandatory alternative group, across %d mapped column(s)"
                    % len(rule_cols)
                ),
            )

            after_rows = df.height
            removed_missing_count = before_rows - after_rows
            field_completeness[
                "rows_removed_at_read_mandatory_gate"
            ] = removed_missing_count

            _emit(logger, "Read-stage mandatory-field validation:")
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
        df, sample_column_dict, coordinate_qc = harmonise_coordinates_and_alleles(
            chromosome="All_Chrs",
            df=df,
            sample_column_dict=sample_column_dict,
            drop_mt=None,
            policies=policies,
            logger=logger,
            rejects=rejects,
            return_qc_info=True,
        )
        coordinate_reasons = coordinate_qc["removed_by_reason"]
        removed_coords = int(
            coordinate_reasons.get("invalid_chromosome", 0)
            + coordinate_reasons.get("invalid_position", 0)
        )
        removed_unsupported_chromosomes = int(
            coordinate_reasons.get("unsupported_chromosome", 0)
        )
        if (
            removed_coords + removed_unsupported_chromosomes
            != int(coordinate_qc["removed_total"])
        ):
            raise ReconciliationError(
                "Coordinate validation counts do not balance: %d invalid "
                "coordinate row(s) + %d unsupported-chromosome row(s) != %d "
                "total coordinate-step removal(s)."
                % (
                    removed_coords,
                    removed_unsupported_chromosomes,
                    int(coordinate_qc["removed_total"]),
                )
            )
        _emit(logger, "Removed %d invalid coordinate rows" % removed_coords)
        _emit(
            logger,
            "Removed %d rows outside the configured chromosome scope"
            % removed_unsupported_chromosomes,
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
        sample_count_columns = [
            column
            for key in ("ncontrol_col", "ncase_col")
            if (column := optional_text(sample_column_dict.get(key))) is not None
        ]
        compound_numeric_columns = [
            column
            for key in ("imp_info_col",)
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
            sample_count_columns=sample_count_columns,
            compound_numeric_columns=compound_numeric_columns,
            conversion_counts=field_completeness["numeric_conversion_by_column"],
        )

        _emit(logger, "Recovery successful: %d rows loaded." % df.height)

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
        # MULTI-VALUE STUDY INFO
        # =========================================================
        # Scalarise a configured delimited study INFO column before
        # duplicate validation. This makes the configured INFO tie-breaker
        # operational and ensures every later dataset/chromosome step sees the
        # same schema-validated aggregation decision.
        df, sample_column_dict = normalise_imputation_quality_column(
            df, sample_column_dict, policies, logger=logger,
        )

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
            dup_file.parent.mkdir(parents=True, exist_ok=True)
            dup_df.write_csv(dup_file, separator=table_delimiter)
            _emit(logger, "Duplicate assessment saved to %s" % dup_file)

        unsafe_duplicate_groups = (
            duplicate_statistics["conflicting_groups"]
            + duplicate_statistics["swapped_orientation_groups"]
        )
        if (
            unsafe_duplicate_groups
            and policies.get("duplicates.conflicting_action") == "fail_dataset"
        ):
            categories = []
            if duplicate_statistics["conflicting_groups"]:
                categories.append(
                    "%d conflicting duplicate group(s)"
                    % duplicate_statistics["conflicting_groups"]
                )
            if duplicate_statistics["swapped_orientation_groups"]:
                categories.append(
                    "%d swapped-orientation duplicate group(s)"
                    % duplicate_statistics["swapped_orientation_groups"]
                )
            raise ValueError(
                "%s were found. Their rows were recorded with distinct duplicate "
                "reasons in %s and the dataset stopped because "
                "duplicates.conflicting_action is 'fail_dataset'."
                % (" and ".join(categories), dup_file)
            )

        _emit(logger, "Removed %d duplicate variants" % removed_duplicates_count)
        measured_ledger = {
            "missing_read_mandatory_value": int(removed_missing_count),
            "invalid_coordinates": int(removed_coords),
            "unsupported_chromosomes": int(removed_unsupported_chromosomes),
            "non_standard_alleles": int(non_standard_allele_count),
            "duplicate_variants": int(removed_duplicates_count),
        }
        if rejects is not None:
            reject_counts_after = rejects.counts()
            stage_reject_counts = {
                reason: int(count) - int(reject_counts_before.get(reason, 0))
                for reason, count in reject_counts_after.items()
                if int(count) - int(reject_counts_before.get(reason, 0))
            }
            collector_ledger = {
                "missing_read_mandatory_value": int(
                    stage_reject_counts.get("missing_read_mandatory_value", 0)
                ),
                "invalid_coordinates": int(
                    stage_reject_counts.get("invalid_chromosome", 0)
                    + stage_reject_counts.get("invalid_position", 0)
                ),
                "unsupported_chromosomes": int(
                    stage_reject_counts.get("unsupported_chromosome", 0)
                ),
                "non_standard_alleles": int(
                    stage_reject_counts.get("null_allele", 0)
                    + stage_reject_counts.get("non_standard_allele", 0)
                ),
                "duplicate_variants": int(
                    stage_reject_counts.get("duplicate_variant", 0)
                    + stage_reject_counts.get("conflicting_duplicate", 0)
                    + stage_reject_counts.get(
                        "swapped_orientation_duplicate", 0,
                    )
                ),
            }
            if sum(stage_reject_counts.values()) != sum(
                collector_ledger.values()
            ):
                raise ReconciliationError(
                    "Input validation produced a rejected-variant reason that "
                    "is not assigned to the input ledger: %r."
                    % stage_reject_counts
                )
            if collector_ledger != measured_ledger:
                raise ReconciliationError(
                    "Input validation counts disagree with rejected-variant "
                    "provenance. Measured=%r; rejected-file counts=%r."
                    % (measured_ledger, collector_ledger)
                )

        accounted_rows = df.height + sum(measured_ledger.values())
        if accounted_rows != python_read_count:
            raise ReconciliationError(
                "Input validation counts do not balance: %d variants read, "
                "%d retained, and %d attributed to input rejection categories."
                % (
                    python_read_count,
                    df.height,
                    sum(measured_ledger.values()),
                )
            )
        _emit(
            logger,
            "Input validation accounting balanced: %d read = %d rejected + "
            "%d ready for harmonisation."
            % (python_read_count, sum(measured_ledger.values()), df.height),
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
            removed_unsupported_chromosomes,
            non_standard_allele_count,
            removed_duplicates_count,
            removed_missing_count,
            field_completeness,
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


def _previous_line_break(data, end):
    """Return the preceding universal-newline span in ``data[:end]``."""
    lf_index = data.rfind(b"\n", 0, end)
    cr_index = data.rfind(b"\r", 0, end)
    index = max(lf_index, cr_index)
    if index < 0:
        return None
    if data[index] == 10 and index > 0 and data[index - 1] == 13:
        return index - 1, index + 1
    return index, index + 1


def _last_two_line_suffix(data):
    """Retain enough bytes to reconstruct the final logical text line."""
    last_break = _previous_line_break(data, len(data))
    if last_break is None:
        return data
    previous_break = _previous_line_break(data, last_break[0])
    if previous_break is None:
        return data
    return data[previous_break[1]:]


def _final_text_line(data):
    """Decode the final universal-newline-delimited line like ``open_text``."""
    normalized = data.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    if normalized.endswith(b"\n"):
        normalized = normalized[:-1]
    line = normalized.rsplit(b"\n", 1)[-1]
    return line.decode("utf-8", errors="replace").strip()


def _inspect_text_file(file_path, *, io_buffer_bytes):
    """Validate and count a text stream in bounded binary chunks.

    Counting newline bytes in C avoids one Python iteration per GWAS row. The
    bookkeeping preserves ``open_text`` semantics: CR, LF and CRLF are logical
    line endings, ``##`` metadata lines are excluded, an unterminated final
    line is counted, and only the final logical line is decoded.
    """
    if io_buffer_bytes <= 0:
        raise ValueError("io_buffer_bytes must be greater than zero")

    line_terminators = 0
    metadata_lines = 0
    first_bytes = b""
    pattern_overlap = b""
    final_lines = b""
    previous_byte = b""
    last_byte = b""
    has_bytes = False

    with open_binary(file_path) as handle:
        while chunk := handle.read(io_buffer_bytes):
            has_bytes = True
            if len(first_bytes) < 2:
                first_bytes += chunk[:2 - len(first_bytes)]

            # TextIOWrapper's universal-newline mode treats CRLF as one line
            # ending. A pair crossing a chunk boundary therefore needs the
            # same one-count adjustment as a pair contained in the chunk.
            line_terminators += (
                chunk.count(b"\n")
                + chunk.count(b"\r")
                - chunk.count(b"\r\n")
            )
            if previous_byte == b"\r" and chunk.startswith(b"\n"):
                line_terminators -= 1

            # The two-byte overlap is shorter than either search pattern, so
            # it detects chunk-boundary line starts without double counting.
            metadata_window = pattern_overlap + chunk
            metadata_lines += metadata_window.count(b"\n##")
            metadata_lines += metadata_window.count(b"\r##")
            pattern_overlap = metadata_window[-2:]

            final_lines = _last_two_line_suffix(final_lines + chunk)
            previous_byte = chunk[-1:]
            last_byte = previous_byte

    if first_bytes.startswith(b"##"):
        metadata_lines += 1
    physical_lines = line_terminators + int(
        has_bytes and last_byte not in (b"\r", b"\n")
    )
    record_count = physical_lines - metadata_lines
    return _final_text_line(final_lines), record_count


def inspect_summary_statistics_file(
    file_path,
    *,
    io_buffer_bytes,
    policies=None,
    logger=None,
):
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
        last_line, record_count = _inspect_text_file(
            source,
            io_buffer_bytes=io_buffer_bytes,
        )
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


def check_file_truncation(
    file_path,
    *,
    io_buffer_bytes,
    policies=None,
    logger=None,
):
    """Compatibility wrapper returning only truncation warnings."""
    warnings, _line_count = inspect_summary_statistics_file(
        file_path,
        io_buffer_bytes=io_buffer_bytes,
        policies=policies,
        logger=logger,
    )
    return warnings
