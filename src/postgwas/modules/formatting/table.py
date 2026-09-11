"""Canonical in-memory representation of one harmonised GWAS-VCF."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import get_args, Iterable, Literal

import polars as pl

from postgwas.config.models.modules.formatting import FormattingDuplicatePolicy
from postgwas.core.dataframes import position_expression
from postgwas.core.errors import FormattingError
from postgwas.core.paths import configured_output_path
from postgwas.core.statistics import negative_log10_to_raw_p
from postgwas.core.variant_identifiers import (
    IDENTIFIER_TEMPLATE_FIELDS,
    identifier_template_parts,
)
from postgwas.core.vcf import (
    extract_vcf_table,
    read_vcf_header,
    validate_harmonised_vcf_header,
)


@dataclass(frozen=True)
class StudyDesign:
    """Study type inferred from the harmonised per-variant sample counts."""

    trait_type: Literal["binary", "quantitative"]
    rows: int
    case_counts_present: int
    control_counts_present: int


@dataclass(frozen=True)
class FormattingColumnSpec:
    """One ordered output expression, including duplicate canonical sources."""

    source: str
    destination: str
    transformation: str | None = None
    integer: bool = False


def load_harmonised_vcf(
    vcf_path: str | Path,
    work_table: str | Path,
    bcftools: str,
    config,
    *,
    logger=None,
    return_header_evidence: bool = False,
    validated_header_evidence: dict[str, str | list[str]] | None = None,
    validated_sample: str | None = None,
) -> pl.DataFrame | tuple[pl.DataFrame, dict[str, str | list[str]]]:
    """Validate and parse one PostGWAS-harmonised VCF without re-harmonising."""
    columns = dict(config.vcf_fields.root)
    canonical = config.canonical_columns
    header_evidence = validated_header_evidence
    if header_evidence is None:
        header = read_vcf_header(
            vcf_path,
            bcftools,
            logger=logger,
            error_type=FormattingError,
        )
        header_evidence = validate_harmonised_vcf_header(
            header, config, vcf_path=vcf_path, logger=logger,
        )
    if logger is not None:
        logger.record(
            "VALIDATE",
            "formatter_input_vcf_contract",
            status="PASSED",
            **header_evidence,
        )
    extract_vcf_table(
        vcf_path,
        work_table,
        str(header_evidence["postgwas_dataset_id"]),
        columns,
        bcftools,
        delimiter=config.runtime.table_delimiter,
        io_buffer_bytes=config.runtime.io_buffer_bytes,
        include_expression=config.vcf_include_expression,
        allow_undefined_tags=True,
        validated_sample=validated_sample,
        logger=logger,
        error_type=FormattingError,
        purpose="Extracting harmonised GWAS-VCF fields for formatting",
    )
    try:
        string_schema = {name: pl.String for name in columns}
        frame = pl.read_csv(
            work_table,
            separator=config.runtime.table_delimiter,
            null_values=config.runtime.input_null_values,
            schema=string_schema,
            low_memory=True,
            ignore_errors=False,
        )
    except (OSError, pl.exceptions.PolarsError) as exc:
        raise FormattingError("Cannot read the extracted VCF table: %s" % exc) from exc
    if frame.is_empty():
        raise FormattingError("The VCF contains no biallelic variants.")

    numeric = list(config.numeric_columns)
    frame = frame.with_columns(
        position_expression(canonical.position, frame.schema[canonical.position]),
        *[pl.col(column).cast(pl.Float64, strict=False) for column in numeric],
    )
    chromosome = pl.col(canonical.chromosome)
    noncanonical_chromosome = (
        chromosome.is_null()
        | (chromosome == "")
        | (chromosome != chromosome.str.strip_chars())
        | ~chromosome.str.contains(config.input_contract.chromosome_pattern)
    )

    def invalid_allele(column: str) -> pl.Expr:
        allele = pl.col(column)
        return (
            allele.is_null()
            | (allele == "")
            | (allele != allele.str.strip_chars())
            | ~allele.str.contains(config.input_contract.allele_pattern)
        )

    invalid_position = (
        pl.col(canonical.position).is_null()
        | (pl.col(canonical.position) < 1)
    )
    invalid_alleles = invalid_allele(
        canonical.reference_allele
    ) | invalid_allele(canonical.alternate_allele)
    violations = frame.select(
        noncanonical_chromosome.fill_null(True).sum().alias("chromosome"),
        invalid_position.fill_null(True).sum().alias("position"),
        invalid_alleles.fill_null(True).sum().alias("alleles"),
    ).row(0, named=True)
    if any(int(value) for value in violations.values()):
        raise FormattingError(
            "Input violates the configured PostGWAS-harmonised GWAS-VCF row "
            "contract: %s noncanonical chromosome labels; %s invalid "
            "positions; %s noncanonical or invalid REF/ALT allele records. "
            "Re-run PostGWAS harmonisation; formatter does not repair "
            "chromosome or allele values."
            % tuple(f"{int(value):,}" for value in violations.values())
        )
    if logger is not None:
        logger.record(
            "VALIDATE",
            "formatter_input_row_contract",
            status="PASSED",
            variants=frame.height,
            chromosome_pattern=config.input_contract.chromosome_pattern,
            allele_pattern=config.input_contract.allele_pattern,
            values_rewritten=0,
        )
    if return_header_evidence:
        return frame, header_evidence
    return frame


def variant_identifier_expression(config, identifier_type: str) -> pl.Expr:
    """Build the configured formatter identifier expression for one target."""
    canonical = config.canonical_columns
    policy = config.variant_identifiers
    if identifier_type == "rsid":
        return (
            pl.col(canonical.variant_id)
            .cast(pl.String)
            .str.strip_chars()
            .str.extract(policy.rsid_extraction_pattern, group_index=1)
        )
    if identifier_type == "unique":
        fields = {
            field: getattr(canonical, field)
            for field in IDENTIFIER_TEMPLATE_FIELDS
        }
        identifier = pl.lit("")
        for literal, field, _, _ in identifier_template_parts(
            policy.unique_id_template,
        ):
            if literal:
                identifier += pl.lit(literal)
            if field is not None:
                identifier += pl.col(fields[field]).cast(pl.String)
        return identifier
    # Pydantic rejects this; retain a clear boundary error for artifact consumers.
    raise FormattingError("Unknown variant identifier type: %s" % identifier_type)


def select_variant_identifiers(
    frame: pl.DataFrame,
    config,
    identifier_type: str,
    *,
    duplicate_policy: str | None = None,
) -> tuple[pl.DataFrame, dict[str, int | str]]:
    """Select one configured identifier convention for a formatter target."""
    canonical = config.canonical_columns
    identifier = variant_identifier_expression(config, identifier_type)

    selected = frame.with_columns(
        identifier.alias(canonical.resolved_variant_id)
    )
    usable = selected.filter(
        pl.col(canonical.resolved_variant_id).is_not_null()
        & (pl.col(canonical.resolved_variant_id) != "")
    )
    missing_excluded = selected.height - usable.height
    if usable.is_empty():
        raise FormattingError(
            "No variants have a usable %s identifier. Check the VCF ID field or "
            "select the configured unique-ID convention." % identifier_type
        )
    duplicated = usable.filter(
        pl.col(canonical.resolved_variant_id).is_duplicated()
    )
    duplicate_rows = duplicated.height
    duplicate_groups = (
        int(duplicated.select(
            pl.col(canonical.resolved_variant_id).n_unique()
        ).item())
        if duplicate_rows
        else 0
    )
    resolved_policy = (
        config.variant_identifiers.default_duplicate_policy
        if duplicate_policy is None
        else duplicate_policy
    )
    usable, resolution_qc = resolve_duplicate_identifiers(
        usable,
        config,
        identifier_type,
        duplicate_policy=resolved_policy,
    )
    duplicate_rows_excluded = int(
        resolution_qc["identifier_duplicate_rows_excluded"]
    )
    total_excluded = missing_excluded + duplicate_rows_excluded
    return usable, {
        "variant_id_type": identifier_type,
        "identifier_rows_in": selected.height,
        "identifier_rows_out": usable.height,
        "identifier_rows_excluded": total_excluded,
        "identifier_missing_rows_excluded": missing_excluded,
        "identifier_duplicate_groups": duplicate_groups,
        "identifier_duplicate_rows": duplicate_rows,
        **resolution_qc,
    }


def _temporary_column(columns: Iterable[str], stem: str) -> str:
    """Return an internal column name that cannot replace user data."""
    occupied = set(columns)
    candidate = stem
    while candidate in occupied:
        candidate += "_"
    return candidate


def _duplicate_ranking_expression(config, duplicate_policy: str) -> tuple[str, pl.Expr]:
    """Return the configured source and valid score for one ranked policy."""
    canonical = config.canonical_columns
    if duplicate_policy in {"lowest_p", "most_significant"}:
        source = canonical.negative_log10_p_value
        value = pl.col(source).cast(pl.Float64, strict=False)
        valid = value.is_not_null() & value.is_finite() & (value >= 0)
        score = pl.when(valid).then(value).otherwise(None)
    elif duplicate_policy == "highest_maf":
        source = canonical.effect_allele_frequency
        value = pl.col(source).cast(pl.Float64, strict=False)
        valid = (
            value.is_not_null()
            & value.is_finite()
            & (value >= 0)
            & (value <= 1)
        )
        score = pl.when(valid).then(
            pl.min_horizontal(value, 1.0 - value)
        ).otherwise(None)
    elif duplicate_policy == "highest_info":
        source = canonical.imputation_quality
        value = pl.col(source).cast(pl.Float64, strict=False)
        valid = (
            value.is_not_null()
            & value.is_finite()
            & (value >= 0)
            & (value <= 1)
        )
        score = pl.when(valid).then(value).otherwise(None)
    else:
        raise FormattingError(
            "Duplicate-ID policy %s has no ranking definition."
            % duplicate_policy
        )
    return source, score


def resolve_duplicate_identifiers(
    frame: pl.DataFrame,
    config,
    identifier_type: str,
    *,
    duplicate_policy: str,
) -> tuple[pl.DataFrame, dict[str, int | str]]:
    """Resolve already-selected IDs without choosing a row by input order.

    Exact repeated records are safe to collapse because all extracted values
    agree. Ranked policies keep a row only when one valid score is strictly
    greater than every other score in the conflicting ID group. Missing ranks
    and tied maxima exclude the complete group.
    """
    identifier = config.canonical_columns.resolved_variant_id
    if identifier not in frame.columns:
        raise FormattingError(
            "Resolved variant identifier column is missing: %s" % identifier
        )
    accepted = {"allow", *get_args(FormattingDuplicatePolicy)}
    if duplicate_policy not in accepted:
        raise FormattingError(
            "Unknown duplicate-ID policy: %s" % duplicate_policy
        )

    if duplicate_policy == "allow":
        return frame, {
            "identifier_duplicate_policy": duplicate_policy,
            "identifier_exact_duplicate_rows_collapsed": 0,
            "identifier_conflicting_duplicate_groups": 0,
            "identifier_conflicting_duplicate_rows": 0,
            "identifier_duplicate_groups_resolved_by_policy": 0,
            "identifier_duplicate_groups_unresolved": 0,
            "identifier_ranked_winners_retained": 0,
            "identifier_duplicate_rows_excluded": 0,
        }

    duplicate_marker = _temporary_column(
        frame.columns, "__postgwas_duplicate_identifier",
    )
    initial = frame.with_columns(
        pl.col(identifier).is_duplicated().alias(duplicate_marker)
    )
    duplicated = initial.filter(pl.col(duplicate_marker))
    if duplicated.is_empty():
        collapsed = frame
        exact_rows_collapsed = 0
    else:
        row_column = _temporary_column(
            [*frame.columns, duplicate_marker], "__postgwas_duplicate_input_row",
        )
        working = initial.with_row_index(row_column)
        nonduplicated = working.filter(~pl.col(duplicate_marker))
        duplicated = working.filter(pl.col(duplicate_marker))
        collapsed_duplicates = duplicated.unique(
            subset=list(frame.columns),
            maintain_order=True,
        )
        exact_rows_collapsed = duplicated.height - collapsed_duplicates.height
        collapsed = pl.concat(
            [nonduplicated, collapsed_duplicates],
            how="vertical",
        ).sort(row_column).select(frame.columns)
    duplicate_mask = pl.col(identifier).is_duplicated()
    conflicting = collapsed.filter(duplicate_mask)
    conflict_rows = conflicting.height
    conflict_groups = (
        conflicting.select(pl.col(identifier).n_unique()).item()
        if conflict_rows
        else 0
    )

    if conflict_rows and duplicate_policy == "error":
        raise FormattingError(
            "%s conflicting VCF records across %s duplicated %s identifier "
            "groups remain after collapsing %s exact repeated records. Set "
            "--duplicate-id-policy exclude_all to remove every ambiguous row, "
            "or provide an appropriate target reference."
            % (
                f"{conflict_rows:,}",
                f"{int(conflict_groups):,}",
                identifier_type,
                f"{exact_rows_collapsed:,}",
            )
        )

    resolved_groups = 0
    unresolved_groups = 0
    ranked_winners = 0
    ranking_column = None
    if not conflict_rows:
        output = collapsed
    elif duplicate_policy == "lowest_p":
        ranking_column, score = _duplicate_ranking_expression(
            config, duplicate_policy,
        )
        row_column = _temporary_column(
            collapsed.columns, "__postgwas_duplicate_row",
        )
        score_column = _temporary_column(
            [*collapsed.columns, row_column], "__postgwas_duplicate_score",
        )
        working = collapsed.with_row_index(row_column).with_columns(
            score.alias(score_column)
        )
        nonduplicated = working.filter(~pl.col(identifier).is_duplicated())
        winners = (
            working.filter(pl.col(identifier).is_duplicated())
            .sort(
                [identifier, score_column, row_column],
                descending=[False, True, False],
                nulls_last=True,
            )
            .unique(subset=[identifier], keep="first", maintain_order=True)
        )
        resolved_groups = int(conflict_groups)
        ranked_winners = winners.height
        output = pl.concat(
            [
                nonduplicated.select(working.columns),
                winners.select(working.columns),
            ],
            how="vertical",
        ).sort(row_column).drop(row_column, score_column)
    elif duplicate_policy == "exclude_all":
        output = collapsed.filter(~duplicate_mask)
        unresolved_groups = int(conflict_groups)
    else:
        ranking_column, score = _duplicate_ranking_expression(
            config, duplicate_policy,
        )
        if ranking_column not in collapsed.columns:
            raise FormattingError(
                "Duplicate-ID policy %s requires configured canonical column %s."
                % (duplicate_policy, ranking_column)
            )
        row_column = _temporary_column(collapsed.columns, "__postgwas_duplicate_row")
        score_column = _temporary_column(
            [*collapsed.columns, row_column], "__postgwas_duplicate_score",
        )
        maximum_column = _temporary_column(
            [*collapsed.columns, row_column, score_column],
            "__postgwas_duplicate_maximum",
        )
        working = collapsed.with_row_index(row_column).with_columns(
            score.alias(score_column)
        )
        nonduplicated = working.filter(~pl.col(identifier).is_duplicated())
        candidates = working.filter(
            pl.col(identifier).is_duplicated()
            & pl.col(score_column).is_not_null()
        )
        maxima = candidates.group_by(identifier).agg(
            pl.col(score_column).max().alias(maximum_column)
        )
        top = candidates.join(maxima, on=identifier, how="inner").filter(
            (pl.col(maximum_column) - pl.col(score_column)).abs()
            <= config.variant_identifiers.duplicate_rank_tolerance
        )
        unique_top_ids = top.group_by(identifier).len().filter(
            pl.col("len") == 1
        ).select(identifier)
        winners = top.join(unique_top_ids, on=identifier, how="inner")
        resolved_groups = unique_top_ids.height
        unresolved_groups = int(conflict_groups) - resolved_groups
        ranked_winners = winners.height
        output = pl.concat(
            [
                nonduplicated.select(working.columns),
                winners.select(working.columns),
            ],
            how="vertical",
        ).sort(row_column).drop(row_column, score_column)

    duplicate_rows_excluded = frame.height - output.height
    if output.is_empty():
        raise FormattingError(
            "No variants remain after duplicate-ID policy %s handled %s rows "
            "across %s duplicated %s identifier groups. Provide a compatible "
            "reference, choose unique identifiers only when they match the "
            "downstream reference, or use policy error to inspect the conflict."
            % (
                duplicate_policy,
                f"{conflict_rows:,}",
                f"{int(conflict_groups):,}",
                identifier_type,
            )
        )
    metrics: dict[str, int | str] = {
        "identifier_duplicate_policy": duplicate_policy,
        "identifier_exact_duplicate_rows_collapsed": exact_rows_collapsed,
        "identifier_conflicting_duplicate_groups": int(conflict_groups),
        "identifier_conflicting_duplicate_rows": conflict_rows,
        "identifier_duplicate_groups_resolved_by_policy": resolved_groups,
        "identifier_duplicate_groups_unresolved": unresolved_groups,
        "identifier_ranked_winners_retained": ranked_winners,
        "identifier_duplicate_rows_excluded": duplicate_rows_excluded,
    }
    if ranking_column is not None:
        metrics["identifier_duplicate_ranking_column"] = ranking_column
    return output, metrics


def validate_unique_identifiers(
    frame: pl.DataFrame,
    config,
    identifier_type: str,
    duplicate_policy: str,
) -> dict[str, int | bool]:
    """Prove that a final formatter candidate has one valid row per ID."""
    identifier = config.canonical_columns.resolved_variant_id
    if identifier not in frame.columns:
        raise FormattingError(
            "Resolved variant identifier column is missing: %s" % identifier
        )
    if frame.is_empty():
        raise FormattingError(
            "Formatter identifier validation cannot accept an empty candidate table."
        )
    observed = frame.select(
        (
            pl.col(identifier).is_null()
            | (pl.col(identifier).cast(pl.String).str.strip_chars() == "")
        ).sum().alias("invalid"),
        pl.col(identifier).is_duplicated().sum().alias("duplicated"),
    ).row(0, named=True)
    invalid = int(observed["invalid"])
    duplicated = int(observed["duplicated"])
    if invalid or duplicated:
        raise FormattingError(
            "Formatter identifier validation failed after duplicate policy %s: "
            "%s invalid and %s duplicated %s candidate rows remain. Refusing "
            "to write downstream artifacts."
            % (
                duplicate_policy,
                f"{invalid:,}",
                f"{duplicated:,}",
                identifier_type,
            )
        )
    return {
        "identifier_uniqueness_validated": True,
        "identifier_validation_rows": frame.height,
        "identifier_unique_values": frame.height,
    }


def complete_rows(
    frame: pl.DataFrame,
    required: Iterable[str],
    condition: pl.Expr | None = None,
) -> tuple[pl.DataFrame, int]:
    """Retain finite, complete required values and report the exact loss."""
    required = list(required)
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise FormattingError("Required canonical fields are missing: %s" % ", ".join(missing))
    mask = pl.all_horizontal(pl.col(required).is_not_null())
    numeric = [
        column for column in required
        if frame.schema.get(column) in (pl.Float32, pl.Float64)
    ]
    if numeric:
        mask &= pl.all_horizontal(pl.col(numeric).is_finite())
    if condition is not None:
        mask &= condition
    retained = frame.filter(mask.fill_null(False))
    if retained.is_empty():
        raise FormattingError(
            "No variants have all values required by this output format."
        )
    return retained, frame.height - retained.height


def validation_expression(validation, *, positive_columns=()) -> pl.Expr:
    """Build reusable row-validity checks from a typed formatter schema."""
    condition = pl.lit(True)
    for column in [*validation.positive_columns, *positive_columns]:
        condition &= pl.col(column) > 0
    for column in validation.nonnegative_columns:
        condition &= pl.col(column) >= 0
    for column in validation.open_unit_interval_columns:
        condition &= (pl.col(column) > 0) & (pl.col(column) < 1)
    for column in validation.closed_unit_interval_columns:
        condition &= (pl.col(column) >= 0) & (pl.col(column) <= 1)
    return condition


def mapped_table(
    frame: pl.DataFrame,
    schema,
    minimum_p_value: float,
    *,
    columns: dict[str, str] | None = None,
    column_specs: Iterable[FormattingColumnSpec] | None = None,
) -> pl.DataFrame:
    """Apply one ordered canonical-column -> tool-column mapping."""
    if columns is not None and column_specs is not None:
        raise FormattingError("Provide columns or column_specs, not both.")
    if column_specs is None:
        mapping = schema.columns if columns is None else columns
        specs = [
            FormattingColumnSpec(
                source=source,
                destination=destination,
                transformation=schema.transformations.get(source),
                integer=source in schema.integer_columns,
            )
            for source, destination in mapping.items()
        ]
    else:
        specs = list(column_specs)
    missing = sorted({
        spec.source for spec in specs if spec.source not in frame.columns
    })
    if missing:
        raise FormattingError(
            "Configured formatter source columns are missing: %s" % ", ".join(missing)
        )
    expressions = []
    for spec in specs:
        source = spec.source
        destination = spec.destination
        transform = spec.transformation
        if transform == "negative_log10_to_raw_p":
            expression = negative_log10_to_raw_p(
                source, minimum_p_value, output_name=destination,
            )
        elif transform == "effect_frequency_to_minor_frequency":
            value = pl.col(source).cast(pl.Float64, strict=False)
            expression = pl.min_horizontal(value, 1.0 - value).alias(destination)
        elif transform is None:
            expression = pl.col(source).alias(destination)
        else:  # Pydantic rejects this; retain a clear boundary error for callers.
            raise FormattingError("Unknown formatter transformation: %s" % transform)
        if spec.integer:
            expression = expression.cast(pl.Int64, strict=False).alias(destination)
        expressions.append(expression)
    return frame.select(*expressions)


def infer_study_design(
    frame: pl.DataFrame,
    case_count_column: str,
    control_count_column: str,
) -> StudyDesign:
    """Infer binary versus quantitative from NC/NCO values, never metadata."""
    required = [
        column for column in (case_count_column, control_count_column)
        if column not in frame.columns
    ]
    if required:
        raise FormattingError(
            "Cannot infer the study type because canonical sample-count fields are "
            "missing: %s." % ", ".join(required)
        )
    counts = frame.select(
        pl.col(case_count_column).is_not_null().sum().alias("cases"),
        pl.col(control_count_column).is_not_null().sum().alias("controls"),
    ).row(0, named=True)
    cases = int(counts["cases"])
    controls = int(counts["controls"])
    if controls == 0:
        raise FormattingError(
            "Cannot infer the study type because configured control/sample-size "
            "column %s is missing for every variant. Harmonisation must supply it."
            % control_count_column
        )
    return StudyDesign(
        trait_type="quantitative" if cases == 0 else "binary",
        rows=frame.height,
        case_counts_present=cases,
        control_counts_present=controls,
    )


def count_bounded_negative_log10_values(
    frame: pl.DataFrame,
    source_column: str,
    minimum_p_value: float,
) -> int:
    """Count -log10(p) values bounded by a configured raw-p minimum."""
    maximum_lp = -math.log10(float(minimum_p_value))
    return int(
        frame.select((pl.col(source_column) > maximum_lp).sum()).item() or 0
    )


def transformation_source(schema, transformation: str) -> str:
    """Resolve the unique canonical source for a configured transformation."""
    sources = [
        source for source, configured in schema.transformations.items()
        if configured == transformation
    ]
    if len(sources) != 1:
        raise FormattingError(
            "Formatter schema must configure exactly one %s transformation; found %s."
            % (transformation, len(sources))
        )
    return sources[0]


__all__ = [
    "FormattingColumnSpec",
    "FormattingError",
    "StudyDesign",
    "complete_rows",
    "count_bounded_negative_log10_values",
    "infer_study_design",
    "load_harmonised_vcf",
    "mapped_table",
    "resolve_duplicate_identifiers",
    "select_variant_identifiers",
    "transformation_source",
    "validate_harmonised_vcf_header",
    "validate_unique_identifiers",
    "validation_expression",
    "variant_identifier_expression",
]
