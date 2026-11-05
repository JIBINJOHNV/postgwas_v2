"""One row-preserving direct/swapped allele join for reference annotations."""

from __future__ import annotations

from typing import Callable, Mapping, Sequence, Type

import polars as pl

from .statistics import frequency_maf_screen
from .variant_columns import canonicalize_variant_frame


def unused_column_name(
    prefix: str, *frames: pl.DataFrame, reserved: Sequence[str] = (),
) -> str:
    """Preserve existing columns and destinations planned for later expressions."""
    name = prefix
    occupied = {column for frame in frames for column in frame.columns}
    occupied.update(reserved)
    while name in occupied:
        name += "_"
    return name


def biallelic_reference_site_expression(columns):
    """Recognize one DNA allele pair across all listed records at a site.

    Frequency inversion requires two alleles. A reference table's declared
    complete site membership remains an input assumption; this rejects known
    extra alleles before matching can hide them, while retaining sequence indels.
    """
    first, second = pl.col(columns["ea"]), pl.col(columns["oa"])
    return (
        (pl.struct([columns["ea"], columns["oa"]]).n_unique().over(
            [columns["chr"], columns["pos"]],
        ) == 1)
        & first.str.contains(r"^[ACGT]+$")
        & second.str.contains(r"^[ACGT]+$")
        & (first != second)
    ).fill_null(False)


def _reference_frequency_sites(
    grouped, columns, value_column, retained_component, output_columns,
    relevant_sites=None,
):
    """Resolve REF AF from complete, mutually exclusive population alleles.

    The supplied reference must retain all population alleles at each site.
    A single DNA pair retains the existing sequence-indel contract. Multiple
    pairs are recoverable only as distinct SNV alternatives to one common REF;
    collocated indels need original haplotype membership and are not summed.
    VCF AF is defined per ALT: https://samtools.github.io/hts-specs/VCFv4.5.pdf.
    Column roles remain configurable; REF means the declared other allele.
    """
    if relevant_sites is not None:
        # Coordinate-only projection retains EVERY reference allele at each
        # needed site. Duplicate policy was already applied to the full panel.
        grouped = grouped.join(
            relevant_sites, on=[columns["chr"], columns["pos"]], how="semi",
        )
    frequency, reason, multiallelic = output_columns
    names = {
        key: unused_column_name(
            "__postgwas_site_" + key, grouped, reserved=output_columns,
        )
        for key in ("alleles", "refs", "dna", "snv", "retained", "valid", "sum")
    }
    alt, ref = pl.col(columns["ea"]), pl.col(columns["oa"])
    value = pl.col(value_column)
    dna_pair = (
        alt.str.contains(r"^[ACGT]+$") & ref.str.contains(r"^[ACGT]+$")
        & (alt != ref)
    ).fill_null(False)
    valid_value = (
        value.is_not_null() & value.is_finite()
        & (value >= 0.0) & (value <= 1.0)
    ).fill_null(False)
    site_keys = [columns["chr"], columns["pos"]]
    multiple = unused_column_name("__postgwas_multiple_allele_keys", grouped)
    prepared = grouped.with_columns(
        pl.struct(site_keys).is_duplicated().alias(multiple),
        dna_pair.alias(names["dna"]),
        ((alt.str.len_chars() == 1) & (ref.str.len_chars() == 1))
        .fill_null(False).alias(names["snv"]),
        retained_component.alias(names["retained"]),
        valid_value.alias(names["valid"]),
    )
    # Most sites have one allele key. Avoid per-group aggregate state for their
    # existing 1-AF calculation; full-set aggregation is needed only for repeats.
    single = prepared.filter(~pl.col(multiple)).select(
        *site_keys,
        pl.lit(1, dtype=pl.UInt32).alias(names["alleles"]),
        pl.lit(1, dtype=pl.UInt32).alias(names["refs"]),
        *[pl.col(names[key]) for key in ("dna", "snv", "retained", "valid")],
        value.alias(names["sum"]),
    )
    combined = prepared.filter(pl.col(multiple)).group_by(site_keys).agg([
        pl.len().alias(names["alleles"]),
        ref.n_unique().alias(names["refs"]),
        *[pl.col(names[key]).all().alias(names[key])
          for key in ("dna", "snv", "retained", "valid")],
        # sum() skips nulls; a missing ALT must invalidate the derived value.
        value.sum().alias(names["sum"]),
    ])
    sites = pl.concat([single, combined], how="vertical").with_columns(
        pl.when(
            (pl.col(names["refs"]) != 1) | ~pl.col(names["dna"])
            | ((pl.col(names["alleles"]) > 1) & ~pl.col(names["snv"]))
        ).then(pl.lit("unsupported_reference_allele_set"))
        .when(~pl.col(names["retained"]))
        .then(pl.lit("discarded_or_conflicting_reference_component"))
        .when(~pl.col(names["valid"]))
        .then(pl.lit("missing_or_invalid_reference_frequency_component"))
        .when((pl.col(names["sum"]) < 0.0) | (pl.col(names["sum"]) > 1.0))
        .then(pl.lit("reference_alt_frequency_sum_out_of_range"))
        .otherwise(pl.lit(None, dtype=pl.String)).alias(reason),
        (pl.col(names["alleles"]) > 1).alias(multiallelic),
    ).with_columns(
        pl.when(pl.col(reason).is_null())
        .then(1.0 - pl.col(names["sum"]))
        .otherwise(None).alias(frequency),
    )
    return sites.select(columns["chr"], columns["pos"], *output_columns)


def deduplicate_reference_values(
    reference: pl.DataFrame,
    keys: Sequence[str],
    value_column: str,
    *,
    exact_action: str,
    non_identical_action: str,
    reference_label: str,
    error_type: Type[Exception],
    exact_policy_key: str = "external_reference.exact_duplicate_action",
    non_identical_policy_key: str = (
        "external_reference.non_identical_duplicate_action"
    ),
    frequency_site_columns: Mapping[str, str] | None = None,
    frequency_output_columns: tuple[str, str, str] | None = None,
    frequency_sites: pl.DataFrame | None = None,
    frequency_site_metadata: dict | None = None,
) -> tuple[pl.DataFrame, dict[str, int]]:
    """Resolve repeated allele keys without selecting among different values."""
    if exact_action not in ("keep_one", "discard_all", "fail"):
        raise error_type(
            "%s must be 'keep_one', 'discard_all' or 'fail'; received %r."
            % (exact_policy_key, exact_action)
        )
    if non_identical_action not in ("discard_all", "fail"):
        raise error_type(
            "%s must be 'discard_all' or 'fail'; received %r."
            % (non_identical_policy_key, non_identical_action)
        )
    count_name = unused_column_name("__postgwas_duplicate_count__", reference)
    first_name = unused_column_name("__postgwas_first_value__", reference)
    distinct_name = unused_column_name("__postgwas_distinct_values__", reference)
    value = pl.col(value_column)
    grouped = reference.group_by(list(keys), maintain_order=True).agg([
        pl.len().cast(pl.UInt32).alias(count_name),
        value.first().alias(first_name),
        # n_unique includes null and distinguishes NaN and infinities. Thus a
        # missing value beside a finite annotation is non-identical rather
        # than an invitation to select the finite row silently.
        value.n_unique().alias(distinct_name),
    ])

    duplicate = pl.col(count_name) > 1
    exact = duplicate & (pl.col(distinct_name) == 1)
    non_identical = duplicate & (pl.col(distinct_name) > 1)

    def fail_for(mask: pl.Expr, kind: str, policy_key: str) -> None:
        conflicts = grouped.filter(mask)
        if not conflicts.height:
            return
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
        values = (
            reference.filter(key_mask)
            .get_column(value_column)
            .unique(maintain_order=True)
            .head(10)
            .to_list()
        )
        key_text = ", ".join(
            "%s=%r" % (key, example[key]) for key in keys
        )
        raise error_type(
            "%s contains %s duplicate values for allele-specific reference key "
            "%s: %r. Policy '%s' is 'fail'. Make the values identical, repair "
            "the reference, or choose a permitted discard policy."
            % (reference_label, kind, key_text, values, policy_key)
        )

    if exact_action == "fail":
        fail_for(
            exact,
            "exact",
            exact_policy_key,
        )
    if non_identical_action == "fail":
        fail_for(
            non_identical,
            "non-identical",
            non_identical_policy_key,
        )

    discard = pl.lit(False)
    if exact_action == "discard_all":
        discard = discard | exact
    if non_identical_action == "discard_all":
        discard = discard | non_identical
    counts = grouped.select([
        exact.sum().alias("exact_groups"),
        non_identical.sum().alias("non_identical_groups"),
        pl.when(exact).then(pl.col(count_name)).otherwise(0).sum().alias("exact_rows"),
        pl.when(non_identical).then(pl.col(count_name)).otherwise(0).sum().alias("non_identical_rows"),
        discard.sum().alias("discarded_groups"),
    ]).row(0, named=True)
    counts = {key: int(value or 0) for key, value in counts.items()}

    resolved = grouped.filter(~discard).select([
        *[pl.col(key) for key in keys],
        pl.col(first_name).alias(value_column),
    ])
    statistics = {
        "reference_duplicate_groups": counts["exact_groups"] + counts["non_identical_groups"],
        "reference_exact_duplicate_groups": counts["exact_groups"],
        "reference_non_identical_duplicate_groups": counts["non_identical_groups"],
        "reference_duplicate_groups_discarded": counts["discarded_groups"],
        "reference_exact_duplicate_rows": counts["exact_rows"],
        "reference_non_identical_duplicate_rows": counts["non_identical_rows"],
    }
    if frequency_output_columns is not None:
        sites = _reference_frequency_sites(
            grouped, frequency_site_columns, first_name,
            ~discard & (pl.col(distinct_name) == 1), frequency_output_columns,
            relevant_sites=frequency_sites,
        )
        if frequency_site_metadata is not None:
            frequency_site_metadata["sites"] = sites
        frequency, reason, multiallelic = frequency_output_columns
        statistics.update(sites.select(
            (pl.col(multiallelic) & pl.col(frequency).is_not_null()).sum()
            .alias("reference_multiallelic_sites"),
            pl.col(reason).is_not_null().sum()
            .alias("reference_frequency_unresolved_sites"),
        ).row(0, named=True))
        resolved = resolved.join(
            sites, on=[frequency_site_columns["chr"], frequency_site_columns["pos"]],
            how="left", coalesce=True,
        )
    return resolved, statistics


def prepare_reference_values(
    reference, columns, value_column, policies, *, exact_action,
    non_identical_action, error_type=RuntimeError, warn=None, label="reference table",
    output_columns=None, frequency_output_columns=None, frequency_sites=None,
    frequency_site_metadata=None,
):
    """Normalize allele keys and resolve duplicates before AF screening/joining."""
    keys = [columns[key] for key in ("chr", "pos", "ea", "oa")]
    reference, _normalization = canonicalize_variant_frame(
        reference.select([*keys, value_column]), columns, policies,
        error_type=error_type, warn=warn, label=label,
    )
    reference = reference.with_columns(pl.col(value_column).cast(pl.Float64, strict=False))
    if output_columns is not None:
        reference = reference.rename({columns[key]: output_columns[key] for key in ("chr", "pos", "ea", "oa")})
        keys = [output_columns[key] for key in ("chr", "pos", "ea", "oa")]
    site_columns = output_columns if output_columns is not None else columns
    rows_before = reference.height
    reference, statistics = deduplicate_reference_values(
        reference, keys, value_column, exact_action=exact_action,
        non_identical_action=non_identical_action, reference_label=label,
        error_type=error_type,
        frequency_site_columns=site_columns,
        frequency_output_columns=frequency_output_columns,
        frequency_sites=frequency_sites,
        frequency_site_metadata=frequency_site_metadata,
    )
    return reference, statistics, rows_before


def summarise_allele_assignments(joined, orientation_name, value_name, frequency_columns=None):
    """Recount raw key/value assignments on exactly the caller's surviving rows."""
    count_expressions = [
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
    ]
    if frequency_columns is not None:
        count_expressions.extend([
            ((pl.col(orientation_name) == 1) & pl.col(frequency_columns[2])
             & pl.col(value_name).is_not_null()).sum()
            .alias("reference_frequency_reconstructed_rows"),
            ((pl.col(orientation_name) == 1) & pl.col(value_name).is_null()).sum()
            .alias("reference_frequency_unresolved_rows"),
        ])
    counts = joined.select(count_expressions).row(0, named=True)
    stats = {key: int(value or 0) for key, value in counts.items()}
    return stats


def allele_oriented_left_join(
    study: pl.DataFrame,
    reference: pl.DataFrame,
    *,
    study_columns: Mapping[str, str],
    reference_columns: Mapping[str, str],
    policies,
    value_column: str,
    output_column: str,
    duplicate_exact_action: str,
    duplicate_non_identical_action: str,
    orientations: Sequence[str] = ("direct", "swap"),
    swapped_value: str = "same",
    prefer_non_null_value: bool = False,
    study_columns_canonical: bool = True,
    error_type: Type[Exception] = RuntimeError,
    warn: Callable[[str], None] | None = None,
    reference_label: str = "reference table",
    screen_frequency: bool = True,
    screen_chromosomes: Sequence[str] | None = None,
    require_complete_frequency_site: bool = False,
    frequency_reason_column: str | None = None,
    frequency_biallelic_column: str | None = None,
    retain_assignment_metadata: bool = False,
) -> tuple[pl.DataFrame, str, dict[str, object]]:
    """Attach one numeric annotation using direct and/or swapped alleles.

    The function performs exactly one study/reference join and never changes
    the number or order of study rows. Exact allele-key plus annotation-value
    duplicates follow ``duplicate_exact_action``. Different values for one key
    follow ``duplicate_non_identical_action`` and are never selected by row
    order.
    Direct orientation has precedence when both allele orders provide a value;
    ``prefer_non_null_value=True`` permits a populated swapped record to fill
    an empty direct record.
    ``swapped_value='one_minus'`` is appropriate only for an allele frequency;
    allele-independent annotations use ``'same'``.
    ``require_complete_frequency_site`` resolves the other-allele frequency
    once per complete reference site, including common-REF multiallelic SNPs.
    Unresolved swaps remain null for the caller's configured missing policy.
    Private ``retain_assignment_metadata`` keeps raw assignment/site columns
    and returns their descriptor in statistics for delayed shared finalization.
    It does not select strand or apply final frequency/palindrome policies.
    Study and reference chromosome keys use the same resolved chromosome
    policies before deduplication and matching.
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
    if require_complete_frequency_site and swapped_value != "one_minus":
        raise error_type("Complete frequency-site validation requires an allele-frequency swap.")
    reserved_outputs = tuple(
        name for name in (
            output_column, frequency_reason_column, frequency_biallelic_column,
        ) if name is not None
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
        study, _study_normalization = canonicalize_variant_frame(
            study,
            study_columns,
            policies,
            error_type=error_type,
            warn=warn,
            label="study",
        )

    frequency_columns = None
    if require_complete_frequency_site and (
        "swap" in requested or frequency_biallelic_column is not None
    ):
        planned = list(reserved_outputs)
        frequency_columns = []
        for role in ("frequency", "reason", "multiallelic"):
            name = unused_column_name(
                "__postgwas_reference_" + role, study, reference, reserved=planned,
            )
            planned.append(name)
            frequency_columns.append(name)
        frequency_columns = tuple(frequency_columns)
    site_metadata = {} if retain_assignment_metadata else None
    reference, duplicate_stats, rows_before_deduplication = prepare_reference_values(
        reference, reference_columns, value_column, policies,
        exact_action=duplicate_exact_action, non_identical_action=duplicate_non_identical_action,
        error_type=error_type, warn=warn, label=reference_label, output_columns=study_columns,
        frequency_output_columns=frequency_columns,
        frequency_site_metadata=site_metadata,
        frequency_sites=(
            study.select(study_columns["chr"], study_columns["pos"]).unique()
            if frequency_columns is not None else None
        ),
    )
    rows_after_deduplication = reference.height
    if swapped_value == "one_minus" and screen_frequency:
        duplicate_stats["raw_frequency_screen"] = frequency_maf_screen(
            (reference.filter(pl.col(study_columns["chr"]).is_in(screen_chromosomes))
             if screen_chromosomes is not None else reference),
            value_column,
            float(policies.get("eaf.maf_decision_cutoff")),
        )

    value_name = unused_column_name(
        "__postgwas_reference_value__", study, reference, reserved=reserved_outputs,
    )
    orientation_name = unused_column_name(
        "__postgwas_allele_match_orientation__", study, reference,
        reserved=(*reserved_outputs, value_name),
    )
    priority_name = unused_column_name(
        "__postgwas_match_priority__", study, reference,
        reserved=(*reserved_outputs, value_name, orientation_name),
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
            swapped_expression = (
                pl.col(frequency_columns[0]) if frequency_columns is not None
                else 1.0 - swapped_expression
            )
        frames.append(reference.select([
            pl.col(study_columns["chr"]),
            pl.col(study_columns["pos"]),
            pl.col(study_columns["oa"]).alias(study_columns["ea"]),
            pl.col(study_columns["ea"]).alias(study_columns["oa"]),
            swapped_expression.alias(value_name),
            *([pl.col(column) for column in frequency_columns]
              if frequency_columns is not None else []),
            pl.lit(1, dtype=pl.Int8).alias(orientation_name),
        ]))

    oriented = pl.concat(frames, how="vertical")
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
    # Polars 0.20, which is the package-pinned runtime, has no join-level
    # ``maintain_order`` argument. Preserve the scientific input order
    # explicitly with a collision-safe row index and remove it immediately
    # after the one-to-one left join.
    row_order_name = unused_column_name(
        "__postgwas_study_row_order__", study, oriented, reserved=reserved_outputs,
    )
    joined = (
        study.with_row_index(row_order_name)
        .join(oriented, on=study_keys, how="left", coalesce=True)
        .sort(row_order_name)
        .drop(row_order_name)
    )
    if joined.height != rows_before_join:
        raise error_type(
            "Matching %s changed the study row count from %d to %d; the "
            "allele-oriented join is not row preserving."
            % (reference_label, rows_before_join, joined.height)
        )

    stats = {
        "reference_rows_before_deduplication": rows_before_deduplication,
        "reference_rows_after_deduplication": rows_after_deduplication,
        "reference_duplicate_rows_removed": (
            rows_before_deduplication - rows_after_deduplication
        ),
    }
    stats.update(duplicate_stats)
    if retain_assignment_metadata:
        # Coordinate facts cover unmatched keys and discarded allele groups too.
        # They must be recounted on the actual worker coordinates, not copied
        # from the larger preparation candidate set.
        site_columns = []
        if frequency_columns is not None:
            sites = site_metadata["sites"]
            for role in ("multiallelic_resolved_site", "unresolved_site"):
                site_columns.append(unused_column_name("__postgwas_" + role, joined, reserved=site_columns))
            site_keys = [study_columns["chr"], study_columns["pos"]]
            site_facts = sites.select(
                *site_keys,
                (pl.col(frequency_columns[2]) & pl.col(frequency_columns[0]).is_not_null()).alias(site_columns[0]),
                pl.col(frequency_columns[1]).is_not_null().alias(site_columns[1]),
            )
            joined = (joined.with_row_index(row_order_name)
                .join(site_facts, on=site_keys, how="left", coalesce=True)
                .sort(row_order_name).drop(row_order_name))
        stats["assignment_metadata"] = {
            "value": value_name, "orientation": orientation_name,
            "frequency_columns": frequency_columns, "site_columns": site_columns,
        }
        return joined, orientation_name, stats
    return finish_allele_assignments(
        joined, orientation_name, value_name, stats, output_column=output_column,
        orientations=requested, frequency_columns=frequency_columns, frequency_reason_column=frequency_reason_column,
        frequency_biallelic_column=frequency_biallelic_column, warn=warn,
        reference_label=reference_label,
    )


def finish_allele_assignments(
    joined, orientation_name, value_name, stats, *, output_column, orientations,
    frequency_columns=None, frequency_reason_column=None, frequency_biallelic_column=None,
    site_columns=(), site_keys=(), warn=None, reference_label="reference table",
):
    """Finalize either fresh or cached assignments with identical row/subset QC."""
    stats = dict(stats)
    if frequency_columns is not None and "swap" not in orientations and frequency_biallelic_column is None:
        # Preparation retains complete-site evidence for diagnostics even when
        # the worker's direct-only attachment has no site/reconstruction QC.
        # Apply the same metadata condition as fresh preparation above.
        joined = joined.drop([*frequency_columns, *site_columns])
        frequency_columns, site_columns = None, ()
        for key in ("reference_multiallelic_sites", "reference_frequency_unresolved_sites"):
            stats.pop(key, None)
    stats.update(summarise_allele_assignments(joined, orientation_name, value_name, frequency_columns))
    if site_columns:
        sites = joined.select(*site_keys, *site_columns).unique()
        stats.update(sites.select(
            pl.col(site_columns[0]).sum().alias("reference_multiallelic_sites"),
            pl.col(site_columns[1]).sum().alias("reference_frequency_unresolved_sites"),
        ).row(0, named=True))
        joined = joined.drop(list(site_columns))
    if frequency_columns is not None:
        unresolved = (pl.col(orientation_name) == 1) & pl.col(value_name).is_null()
        if stats["reference_frequency_unresolved_rows"]:
            reasons = joined.filter(unresolved).group_by(frequency_columns[1]).len()
            stats["reference_frequency_unresolved_reasons"] = {
                row[frequency_columns[1]]: int(row["len"])
                for row in reasons.to_dicts()
            }
            if warn is not None:
                warn(
                    "%s: %s REF-frequency assignments remain missing because the "
                    "complete allele set is unsupported or has invalid/discarded "
                    "frequency components; the configured missing-EAF policy applies. "
                    "Reasons: %s."
                    % (reference_label, stats["reference_frequency_unresolved_rows"],
                       stats["reference_frequency_unresolved_reasons"])
                )
        if frequency_reason_column is not None and stats["reference_frequency_unresolved_rows"]:
            joined = joined.with_columns(
                pl.when(unresolved).then(pl.col(frequency_columns[1]))
                .otherwise(pl.lit(None, dtype=pl.String)).alias(frequency_reason_column)
            )
        if frequency_biallelic_column is not None:
            joined = joined.with_columns(
                (~pl.col(frequency_columns[2]) & pl.col(frequency_columns[1]).is_null())
                .fill_null(False).alias(frequency_biallelic_column)
            )
        joined = joined.drop(list(frequency_columns))

    joined = joined.with_columns(
        pl.col(value_name).alias(output_column)
    ).drop(value_name)
    return joined, orientation_name, stats


__all__ = [
    "allele_oriented_left_join",
    "deduplicate_reference_values",
    "unused_column_name",
]
