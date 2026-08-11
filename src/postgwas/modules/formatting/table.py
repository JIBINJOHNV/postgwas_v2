"""Canonical in-memory representation of one harmonised GWAS-VCF."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal

import polars as pl

from postgwas.core.dataframes import chromosome_expression, position_expression
from postgwas.core.paths import configured_output_path
from postgwas.core.statistics import negative_log10_to_raw_p
from postgwas.core.variant_identifiers import (
    IDENTIFIER_TEMPLATE_FIELDS,
    identifier_template_parts,
)
from postgwas.core.vcf import extract_vcf_table


class FormattingError(RuntimeError):
    """A downstream input cannot be produced without corrupting meaning."""


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
    dataset_id: str,
    bcftools: str,
    config,
    *,
    logger=None,
) -> pl.DataFrame:
    """Extract all reusable fields once and return a typed canonical frame."""
    columns = dict(config.vcf_fields.root)
    canonical = config.canonical_columns
    extract_vcf_table(
        vcf_path,
        work_table,
        dataset_id,
        columns,
        bcftools,
        delimiter=config.runtime.table_delimiter,
        io_buffer_bytes=config.runtime.io_buffer_bytes,
        include_expression=config.vcf_include_expression,
        allow_undefined_tags=True,
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
    chromosome = chromosome_expression(
        canonical.chromosome,
        frame.schema[canonical.chromosome],
        strip_chr_prefix=False,
    ).str.replace(config.chromosome_labels.prefix_pattern, "")
    if config.chromosome_labels.aliases:
        chromosome = chromosome.replace(config.chromosome_labels.aliases)
    frame = frame.with_columns(
        chromosome.alias(canonical.chromosome),
        position_expression(canonical.position, frame.schema[canonical.position]),
        pl.col(canonical.reference_allele)
        .cast(pl.String).str.to_uppercase().str.strip_chars(),
        pl.col(canonical.alternate_allele)
        .cast(pl.String).str.to_uppercase().str.strip_chars(),
        *[pl.col(column).cast(pl.Float64, strict=False) for column in numeric],
    )
    invalid_coordinates = frame.filter(
        pl.col(canonical.chromosome).is_null()
        | pl.col(canonical.position).is_null()
        | (pl.col(canonical.position) < 1)
        | pl.col(canonical.reference_allele).is_null()
        | pl.col(canonical.alternate_allele).is_null()
        | (pl.col(canonical.reference_allele) == "")
        | (pl.col(canonical.alternate_allele) == "")
    ).height
    if invalid_coordinates:
        raise FormattingError(
            "%s biallelic VCF records have invalid coordinates or alleles."
            % f"{invalid_coordinates:,}"
        )
    return frame


def select_variant_identifiers(
    frame: pl.DataFrame,
    config,
    identifier_type: str,
) -> tuple[pl.DataFrame, dict[str, int | str]]:
    """Select one configured identifier convention for a formatter target."""
    canonical = config.canonical_columns
    policy = config.variant_identifiers
    if identifier_type == "rsid":
        identifier = (
            pl.col(canonical.variant_id)
            .cast(pl.String)
            .str.strip_chars()
            .str.extract(policy.rsid_extraction_pattern, group_index=1)
        )
    elif identifier_type == "unique":
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
    else:  # Pydantic rejects this; retain a clear boundary error for callers.
        raise FormattingError("Unknown variant identifier type: %s" % identifier_type)

    selected = frame.with_columns(
        identifier.alias(canonical.resolved_variant_id)
    )
    usable = selected.filter(
        pl.col(canonical.resolved_variant_id).is_not_null()
        & (pl.col(canonical.resolved_variant_id) != "")
    )
    excluded = selected.height - usable.height
    if usable.is_empty():
        raise FormattingError(
            "No variants have a usable %s identifier. Check the VCF ID field or "
            "select the configured unique-ID convention." % identifier_type
        )
    duplicate_rows = int(
        usable.select(
            pl.col(canonical.resolved_variant_id).is_duplicated().sum()
        ).item()
        or 0
    )
    if duplicate_rows:
        raise FormattingError(
            "%s VCF records produce duplicate %s identifiers; downstream tools "
            "require unique identifiers." % (f"{duplicate_rows:,}", identifier_type)
        )
    return usable, {
        "variant_id_type": identifier_type,
        "identifier_rows_in": selected.height,
        "identifier_rows_out": usable.height,
        "identifier_rows_excluded": excluded,
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
    "select_variant_identifiers",
    "transformation_source",
    "validation_expression",
]
