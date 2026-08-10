"""One row-preserving direct/swapped allele join for reference annotations."""

from __future__ import annotations

from typing import Callable, Mapping, Sequence, Type

import polars as pl

from postgwas.core.dataframes import (
    chromosome_expression,
    count_non_null,
    position_expression,
    validate_cast_retention,
)


def _unused_name(prefix: str, *frames: pl.DataFrame) -> str:
    name = prefix
    occupied = {column for frame in frames for column in frame.columns}
    while name in occupied:
        name += "_"
    return name


def _deduplicate_reference_values(
    reference: pl.DataFrame,
    keys: Sequence[str],
    value_column: str,
    *,
    source_count_column: str,
    reference_label: str,
    error_type: Type[Exception],
) -> tuple[pl.DataFrame, dict[str, int]]:
    """Collapse safe duplicates and reject conflicting scientific values.

    A finite annotation is preferred over an earlier null, NaN, infinity or
    unparseable value. Equivalent finite values are harmless duplicates. Two
    distinct finite values for the same exact allele-specific key are
    scientifically ambiguous, so the caller must repair the reference rather
    than having this function select one arbitrarily.
    """
    count_name = _unused_name("__postgwas_duplicate_count__", reference)
    first_name = _unused_name("__postgwas_first_value__", reference)
    usable_name = _unused_name("__postgwas_first_usable_value__", reference)
    distinct_name = _unused_name("__postgwas_distinct_usable_values__", reference)
    value = pl.col(value_column)
    usable = value.is_not_null() & value.is_finite()
    grouped = reference.group_by(list(keys), maintain_order=True).agg([
        pl.len().cast(pl.UInt32).alias(count_name),
        value.first().alias(first_name),
        value.filter(usable).first().alias(usable_name),
        value.filter(usable).n_unique().alias(distinct_name),
    ])

    conflicts = grouped.filter(pl.col(distinct_name) > 1)
    if conflicts.height:
        example = conflicts.row(0, named=True)
        key_mask = pl.lit(True)
        for key in keys:
            key_value = example[key]
            condition = (
                pl.col(key).is_null()
                if key_value is None
                else pl.col(key) == pl.lit(key_value)
            )
            key_mask = key_mask & condition
        conflicting_values = (
            reference.filter(key_mask & usable)
            .get_column(value_column)
            .unique(maintain_order=True)
            .head(10)
            .to_list()
        )
        key_text = ", ".join(
            "%s=%r" % (key, example[key]) for key in keys
        )
        raise error_type(
            "%s contains conflicting finite values for duplicate allele-specific "
            "reference key %s: %r. PostGWAS cannot choose one without changing "
            "scientific data. Make the duplicate values identical or remove the "
            "incorrect rows."
            % (reference_label, key_text, conflicting_values)
        )

    duplicate_groups = int(
        grouped.select((pl.col(count_name) > 1).sum()).item() or 0
    )
    preferred_usable = int(
        grouped.select(
            (
                (pl.col(count_name) > 1)
                & ~(
                    pl.col(first_name).is_not_null()
                    & pl.col(first_name).is_finite()
                )
                & pl.col(usable_name).is_not_null()
            ).sum()
        ).item()
        or 0
    )
    resolved = grouped.select([
        *[pl.col(key) for key in keys],
        pl.coalesce([pl.col(usable_name), pl.col(first_name)]).alias(value_column),
        pl.col(count_name).alias(source_count_column),
    ])
    return resolved, {
        "reference_duplicate_groups": duplicate_groups,
        "reference_duplicate_groups_preferred_usable_value": preferred_usable,
        "reference_conflicting_duplicate_groups": 0,
    }


def allele_oriented_left_join(
    study: pl.DataFrame,
    reference: pl.DataFrame,
    *,
    study_columns: Mapping[str, str],
    reference_columns: Mapping[str, str],
    value_column: str,
    output_column: str,
    orientations: Sequence[str] = ("direct", "swap"),
    swapped_value: str = "same",
    deduplicate_reference: bool = True,
    prefer_non_null_value: bool = False,
    study_columns_canonical: bool = True,
    error_type: Type[Exception] = RuntimeError,
    warn: Callable[[str], None] | None = None,
    reference_label: str = "reference table",
    duplicate_policy_name: str | None = None,
) -> tuple[pl.DataFrame, str, dict[str, int]]:
    """Attach one numeric annotation using direct and/or swapped alleles.

    The function performs exactly one study/reference join and never changes
    the number or order of study rows. When reference deduplication is enabled,
    a finite value is preferred over an earlier missing or non-finite duplicate;
    conflicting finite values fail rather than being selected arbitrarily.
    Direct orientation has precedence when both allele orders provide a value;
    ``prefer_non_null_value=True`` permits a populated swapped record to fill
    an empty direct record.
    ``swapped_value='one_minus'`` is appropriate only for an allele frequency;
    allele-independent annotations use ``'same'``.
    """
    requested = tuple(str(value).lower() for value in orientations)
    unsupported = [
        value for value in requested if value not in ("direct", "swap")
    ]
    if not requested or unsupported:
        raise error_type(
            "Allele-oriented matching requires 'direct' and/or 'swap'; received %r."
            % (list(orientations),)
        )
    if swapped_value not in ("same", "one_minus"):
        raise error_type(
            "swapped_value must be 'same' or 'one_minus'; received %r."
            % swapped_value
        )

    required_keys = ("chr", "pos", "ea", "oa")
    missing_study = [key for key in required_keys if not study_columns.get(key)]
    missing_reference = [
        key for key in required_keys if not reference_columns.get(key)
    ]
    required_reference_columns = [
        reference_columns[key] for key in required_keys if reference_columns.get(key)
    ] + [value_column]
    absent_reference = [
        column for column in required_reference_columns if column not in reference.columns
    ]
    if missing_study or missing_reference or absent_reference:
        raise error_type(
            "%s cannot be allele-matched: missing study mappings %s, reference "
            "mappings %s, or reference columns %s."
            % (reference_label, missing_study, missing_reference, absent_reference)
        )

    study_keys = [study_columns[key] for key in required_keys]
    absent_study = [column for column in study_keys if column not in study.columns]
    if absent_study:
        raise error_type(
            "%s cannot be allele-matched because the study is missing columns %s."
            % (reference_label, absent_study)
        )

    if not study_columns_canonical:
        schema = dict(study.schema)
        position_before = count_non_null(study, study_columns["pos"])
        study = study.with_columns([
            chromosome_expression(
                study_columns["chr"], schema[study_columns["chr"]]
            ),
            position_expression(
                study_columns["pos"], schema[study_columns["pos"]]
            ),
            pl.col(study_columns["ea"])
            .cast(pl.String).str.to_uppercase().str.strip_chars(),
            pl.col(study_columns["oa"])
            .cast(pl.String).str.to_uppercase().str.strip_chars(),
        ])
        validate_cast_retention(
            position_before,
            count_non_null(study, study_columns["pos"]),
            "study position",
            error_type=error_type,
            warn=warn,
        )

    reference = reference.select(required_reference_columns)
    reference_schema = dict(reference.schema)
    reference_position_before = count_non_null(
        reference, reference_columns["pos"]
    )
    reference = reference.with_columns([
        chromosome_expression(
            reference_columns["chr"], reference_schema[reference_columns["chr"]]
        ),
        position_expression(
            reference_columns["pos"], reference_schema[reference_columns["pos"]]
        ),
        pl.col(reference_columns["ea"])
        .cast(pl.String).str.to_uppercase().str.strip_chars(),
        pl.col(reference_columns["oa"])
        .cast(pl.String).str.to_uppercase().str.strip_chars(),
        pl.col(value_column).cast(pl.Float64, strict=False),
    ]).rename({
        reference_columns[key]: study_columns[key] for key in required_keys
    })
    validate_cast_retention(
        reference_position_before,
        count_non_null(reference, study_columns["pos"]),
        "%s position" % reference_label,
        error_type=error_type,
        warn=warn,
    )

    source_duplicates = _unused_name(
        "__postgwas_source_duplicate_count__", study, reference
    )
    rows_before_deduplication = reference.height
    duplicate_stats = {
        "reference_duplicate_groups": 0,
        "reference_duplicate_groups_preferred_usable_value": 0,
        "reference_conflicting_duplicate_groups": 0,
    }
    if deduplicate_reference:
        reference, duplicate_stats = _deduplicate_reference_values(
            reference,
            study_keys,
            value_column,
            source_count_column=source_duplicates,
            reference_label=reference_label,
            error_type=error_type,
        )
    else:
        reference = reference.with_columns(
            pl.len().over(study_keys).cast(pl.UInt32).alias(source_duplicates)
        )
    rows_after_deduplication = reference.height

    value_name = _unused_name("__postgwas_reference_value__", study, reference)
    orientation_name = _unused_name(
        "__postgwas_allele_match_orientation__", study, reference
    )
    priority_name = _unused_name("__postgwas_match_priority__", study, reference)
    matched_duplicates = _unused_name(
        "__postgwas_matched_duplicate_count__", study, reference
    )
    reference = reference.rename({value_column: value_name})

    frames = []
    if "direct" in requested:
        frames.append(reference.with_columns(
            pl.lit(0, dtype=pl.Int8).alias(orientation_name)
        ))
    if "swap" in requested:
        swapped_expression = pl.col(value_name)
        if swapped_value == "one_minus":
            swapped_expression = 1.0 - swapped_expression
        frames.append(reference.select([
            pl.col(study_columns["chr"]),
            pl.col(study_columns["pos"]),
            pl.col(study_columns["oa"]).alias(study_columns["ea"]),
            pl.col(study_columns["ea"]).alias(study_columns["oa"]),
            swapped_expression.alias(value_name),
            pl.col(source_duplicates),
            pl.lit(1, dtype=pl.Int8).alias(orientation_name),
        ]))

    oriented = pl.concat(frames, how="vertical")
    duplicate_expression = pl.col(source_duplicates).max().over(study_keys)
    oriented = oriented.with_columns(
        duplicate_expression.alias(matched_duplicates)
    )
    if prefer_non_null_value:
        oriented = (
            oriented
            .with_columns(
                (
                    pl.col(orientation_name)
                    + pl.when(
                        pl.col(value_name).is_null()
                        | ~pl.col(value_name).is_finite()
                    ).then(2).otherwise(0)
                ).alias(priority_name)
            )
            .sort(priority_name, maintain_order=True)
            .unique(subset=study_keys, keep="first", maintain_order=True)
            .drop(priority_name)
        )
    else:
        oriented = oriented.unique(
            subset=study_keys, keep="first", maintain_order=True,
        )

    rows_before_join = study.height
    joined = study.join(
        oriented, on=study_keys, how="left", maintain_order="left",
    )
    if joined.height != rows_before_join:
        raise error_type(
            "Matching %s changed the study row count from %d to %d; the "
            "allele-oriented join is not row preserving."
            % (reference_label, rows_before_join, joined.height)
        )

    duplicated_matches = int(joined.select(
        (
            pl.col(orientation_name).is_not_null()
            & (pl.col(matched_duplicates) > 1)
        ).sum()
    ).item() or 0)
    if duplicated_matches and not deduplicate_reference:
        policy = (
            " Set %s to true, or deduplicate the reference file."
            % duplicate_policy_name
            if duplicate_policy_name else " Deduplicate the reference file."
        )
        raise error_type(
            "%s contains repeated allele-specific rows matching %d study variants.%s"
            % (reference_label, duplicated_matches, policy)
        )

    counts = joined.select([
        (pl.col(orientation_name) == 0).sum().alias("direct_key_matches"),
        (pl.col(orientation_name) == 1).sum().alias("swapped_key_matches"),
        pl.col(orientation_name).is_null().sum().alias("unmatched_rows"),
        (
            (pl.col(orientation_name) == 0) & pl.col(value_name).is_not_null()
        ).sum().alias("direct_value_matches"),
        (
            (pl.col(orientation_name) == 1) & pl.col(value_name).is_not_null()
        ).sum().alias("swapped_value_matches"),
        pl.col(value_name).is_null().sum().alias("missing_values"),
    ]).row(0, named=True)
    stats = {key: int(value or 0) for key, value in counts.items()}
    stats.update({
        "reference_rows_before_deduplication": rows_before_deduplication,
        "reference_rows_after_deduplication": rows_after_deduplication,
        "reference_duplicate_rows_removed": (
            rows_before_deduplication - rows_after_deduplication
        ),
    })
    stats.update(duplicate_stats)

    joined = joined.with_columns(
        pl.col(value_name).alias(output_column)
    ).drop([value_name, source_duplicates, matched_duplicates])
    return joined, orientation_name, stats


__all__ = ["allele_oriented_left_join"]
