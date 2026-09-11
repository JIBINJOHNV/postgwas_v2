"""Lossless input and canonical internal types for study identity/provenance.

Reference tables are separate scientific inputs and must still be normalized
when read. String overrides apply from the first study-file read; the canonical
contract applies after dataset-level coordinate and allele validation.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Mapping, Sequence, Type

import polars as pl

from postgwas.core.dataframes import count_non_null, validate_cast_retention
from postgwas.core.io.delimiters import DelimiterDetectionResult
from postgwas.core.io.tables import read_delimited_table, read_parquet_table


CANONICAL_VARIANT_COLUMNS_KEY = "_variant_columns_canonical"

# These four ordered pairs are the complete set of single-base DNA allele
# pairs that equal their own reverse complement.  This is a biological
# invariant, not a user-selectable analysis threshold.
PALINDROMIC_SNP_PAIRS = ("AT", "TA", "CG", "GC")

# A combined variant identifier may carry alleles or an interval after its
# position. Capture the complete leading numeric token so exact integer-valued
# forms such as ``100.0`` and ``1e5`` retain their mathematical value before the
# shared integer-coordinate validation runs.
_LEADING_POSITION_NUMBER = (
    r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
)

_STRING_VARIANT_COLUMN_KEYS = (
    "chr_col",
    "chr_pos_col",
    "snp_id_col",
    "ea_col",
    "oa_col",
    "pvalue_source_text_col",
)


def canonical_chromosome_expression(
    expression: pl.Expr,
    policies,
    *,
    apply_rename_map: bool = True,
) -> pl.Expr:
    """Return one policy-driven chromosome representation for every input.

    Study tables, reference tables and extracted VCF tables must use this same
    expression before chromosome values are compared. ``apply_rename_map`` is
    false only while the authoritative coordinate step counts which aliases
    the subsequent canonicalization will change for its audit report.
    """
    value = expression.cast(pl.String, strict=False).str.strip_chars()
    if bool(policies.get("chromosome.strip_chr_prefix")):
        value = value.str.replace(r"(?i)^chr", "")
    # Integral-looking chromosome labels commonly arise when a table reader
    # inferred a numeric column. A non-integral label remains unchanged and
    # therefore cannot acquire a valid chromosome silently.
    value = value.str.replace(r"\.0+$", "")
    if bool(policies.get("chromosome.strip_leading_zero")):
        value = value.str.replace(r"^0+([0-9]+)$", "${1}")
    value = value.str.to_uppercase()

    if apply_rename_map:
        rename_map = {
            str(source).strip().upper(): str(target).strip().upper()
            for source, target in dict(
                policies.get("chromosome.rename_map") or {}
            ).items()
        }
        if rename_map:
            value = value.replace(rename_map)
    return value


def canonical_allele_expression(expression: pl.Expr) -> pl.Expr:
    """Return the canonical uppercase DNA-sequence representation."""
    return (
        expression.cast(pl.String, strict=False)
        .str.strip_chars()
        .replace("", None)
        .str.to_uppercase()
    )


def reverse_complement_expression(expression: pl.Expr) -> pl.Expr:
    """Reverse-complement an A/C/G/T sequence without changing row order."""
    allele = canonical_allele_expression(expression)
    return pl.when(allele.str.contains(r"^[ACGT]+$")).then(
        allele
        .str.replace_many(["A", "C", "G", "T"], ["T", "G", "C", "A"])
        .str.reverse()
    ).otherwise(None)


def palindromic_snp_expression(first: pl.Expr, second: pl.Expr) -> pl.Expr:
    """Identify A/T, T/A, C/G or G/C SNP pairs after canonicalization."""
    first = canonical_allele_expression(first)
    second = canonical_allele_expression(second)
    pair = pl.concat_str([first, second])
    return (
        (first.str.len_chars() == 1)
        & (second.str.len_chars() == 1)
        & pair.is_in(PALINDROMIC_SNP_PAIRS)
    ).fill_null(False)


def _trim_shared_allele_padding(
    position: int,
    first: str,
    second: str,
) -> tuple[int, str, str]:
    """Return the VCF-minimal representation without reference left-alignment.

    VCF requires each simple-indel allele to retain an anchor base.  Removing
    identical trailing bases first, then identical leading bases while both
    alleles remain non-empty, is the representation-only part of normalization
    described by the VCF specification and ``bcftools norm``:
    https://samtools.github.io/hts-specs/VCFv4.5.pdf and
    https://samtools.github.io/bcftools/bcftools#norm. It does not use a FASTA
    and therefore cannot left-align an indel in a repeat.
    """
    while len(first) > 1 and len(second) > 1 and first[-1] == second[-1]:
        first = first[:-1]
        second = second[:-1]
    while len(first) > 1 and len(second) > 1 and first[0] == second[0]:
        first = first[1:]
        second = second[1:]
        position += 1
    return position, first, second


def minimal_allele_representation_series(
    frame: pl.DataFrame,
    position_column: str,
    first_allele_column: str,
    second_allele_column: str,
    *,
    position_output: str,
    first_allele_output: str,
    second_allele_output: str,
) -> list[pl.Series]:
    """Build minimally trimmed key columns while preserving source columns.

    Only rows containing a multi-base allele enter the Python trimming loop;
    the usual SNP majority remains in native Polars buffers.  The returned
    columns are for identity keys and concordance only.  They must never replace
    study coordinates or alleles because reference-backed left-alignment is a
    separate downstream operation.
    """
    required = [position_column, first_allele_column, second_allele_column]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(
            "Minimal allele representation is missing column(s): %s."
            % ", ".join(missing)
        )

    canonical = frame.select(
        canonical_position_expression(pl.col(position_column)).alias(
            position_output
        ),
        canonical_allele_expression(pl.col(first_allele_column)).alias(
            first_allele_output
        ),
        canonical_allele_expression(pl.col(second_allele_column)).alias(
            second_allele_output
        ),
    )
    candidate = (
        pl.col(position_output).is_not_null()
        & pl.col(first_allele_output).is_not_null()
        & pl.col(second_allele_output).is_not_null()
        & (pl.col(first_allele_output) != pl.col(second_allele_output))
        & (
            (pl.col(first_allele_output).str.len_chars() > 1)
            | (pl.col(second_allele_output).str.len_chars() > 1)
        )
    )
    indices = canonical.select(
        pl.arg_where(candidate).alias("__minimal_allele_row")
    ).get_column("__minimal_allele_row")
    if indices.is_empty():
        return [
            canonical.get_column(column)
            for column in (
                position_output,
                first_allele_output,
                second_allele_output,
            )
        ]

    positions: list[int] = []
    first_alleles: list[str] = []
    second_alleles: list[str] = []
    for position, first, second in canonical[indices].iter_rows():
        trimmed_position, trimmed_first, trimmed_second = (
            _trim_shared_allele_padding(position, first, second)
        )
        positions.append(trimmed_position)
        first_alleles.append(trimmed_first)
        second_alleles.append(trimmed_second)

    return [
        canonical.get_column(position_output).clone().scatter(indices, positions),
        canonical.get_column(first_allele_output).clone().scatter(
            indices, first_alleles
        ),
        canonical.get_column(second_allele_output).clone().scatter(
            indices, second_alleles
        ),
    ]


def position_text_expression(expression: pl.Expr, extraction: str) -> pl.Expr:
    """Apply the configured combined-coordinate position extraction policy."""
    text = expression.cast(pl.String, strict=False).str.strip_chars()
    if extraction == "leading_digits":
        return text.str.extract(_LEADING_POSITION_NUMBER, 1)
    if extraction == "none":
        return text
    if extraction == "strip_non_digits":
        # Explicit compatibility policy. The canonical YAML warns that this
        # can join distinct numeric components into a different coordinate.
        return text.str.replace_all(r"\D", "")
    raise ValueError("Unknown position.extraction value %r." % extraction)


def _position_components(
    expression: pl.Expr,
) -> tuple[pl.Expr, pl.Expr, pl.Expr, pl.Expr]:
    """Return text, parsed Float64, direct Int64 and validated Int64 forms."""
    text = expression.cast(pl.String, strict=False).str.strip_chars()
    direct_integer = text.cast(pl.Int64, strict=False)
    numeric = text.cast(pl.Float64, strict=False)
    finite_integer = (
        numeric.is_not_null()
        & numeric.is_finite()
        & (numeric == numeric.floor())
    ).fill_null(False)
    converted_integer = numeric.cast(pl.Int64, strict=False)
    canonical = (
        pl.when(direct_integer.is_not_null())
        .then(direct_integer)
        .when(finite_integer)
        .then(converted_integer)
        .otherwise(None)
    )
    return text, numeric, direct_integer, canonical


def canonical_position_expression(expression: pl.Expr) -> pl.Expr:
    """Return a finite, mathematically integral Int64 coordinate or null.

    Human genome coordinates are 1-based integers. Positivity is enforced by
    the authoritative coordinate policy after this representation check, while
    this shared expression prevents non-integral values such as ``100.7`` from
    being silently truncated. Exact integer-valued source forms such as
    ``100``, ``100.0`` and ``1e5`` are retained.
    """
    _text, _numeric, _direct_integer, canonical = _position_components(expression)
    return canonical


def canonicalize_variant_frame(
    frame: pl.DataFrame,
    columns: Mapping[str, object],
    policies,
    *,
    error_type: Type[Exception] = ValueError,
    warn: Callable[[str], None] | None = None,
    label: str = "variant",
    require_position_retention: bool = True,
) -> tuple[pl.DataFrame, dict[str, int]]:
    """Canonicalize configured chromosome, position and allele columns once.

    ``columns`` may contain any of ``chr``, ``pos``, ``ea`` and ``oa``.  Every
    supplied mapping must exist. Positions must be finite and mathematically
    integral; decimal-looking integer values such as ``100.0`` and exact
    scientific notation such as ``1e5`` are accepted, while ``100.7`` becomes
    null rather than being truncated. The returned counts let the authoritative
    coordinate step report malformed values without rescanning the frame.
    ``require_position_retention=False`` is reserved for that study-coordinate
    step so an entirely invalid study column can still be attributed row by row
    through the reject collector; reference callers remain fail-closed.
    """
    supported = ("chr", "pos", "ea", "oa")
    unknown = sorted(set(columns).difference(supported))
    if unknown:
        raise error_type(
            "%s normalization received unsupported variant roles: %s."
            % (label, ", ".join(unknown))
        )
    resolved = {
        role: str(columns[role])
        for role in supported
        if columns.get(role) is not None and str(columns[role]).strip()
    }
    absent = [name for name in resolved.values() if name not in frame.columns]
    if absent:
        raise error_type(
            "%s normalization is missing configured column(s): %s."
            % (label, ", ".join(absent))
        )
    if not resolved:
        return frame, {
            "positions_present": 0,
            "positions_retained": 0,
            "positions_unreadable": 0,
            "positions_non_finite": 0,
            "positions_non_integral": 0,
            "positions_outside_int64": 0,
        }

    expressions: list[pl.Expr] = []
    if "chr" in resolved:
        expressions.append(
            canonical_chromosome_expression(
                pl.col(resolved["chr"]), policies,
            ).alias(resolved["chr"])
        )

    positions_present = 0
    positions_unreadable = 0
    positions_non_finite = 0
    positions_non_integral = 0
    positions_outside_int64 = 0
    if "pos" in resolved:
        position = resolved["pos"]
        source_position = pl.col(position)
        (
            _position_text,
            parsed_position,
            direct_integer,
            canonical_position,
        ) = _position_components(source_position)
        parsed_finite = parsed_position.is_not_null() & parsed_position.is_finite()
        parsed_integral = parsed_finite & (parsed_position == parsed_position.floor())
        position_counts = frame.select([
            source_position.is_not_null().sum().alias("present"),
            (
                source_position.is_not_null() & canonical_position.is_null()
            ).sum().alias("unreadable"),
            (
                source_position.is_not_null()
                & parsed_position.is_not_null()
                & ~parsed_position.is_finite()
            ).sum().alias("non_finite"),
            (
                source_position.is_not_null()
                & parsed_finite
                & (parsed_position != parsed_position.floor())
            ).sum().alias("non_integral"),
            (
                source_position.is_not_null()
                & direct_integer.is_null()
                & parsed_integral
                & canonical_position.is_null()
            ).sum().alias("outside_int64"),
        ]).row(0, named=True)
        positions_present = int(position_counts["present"] or 0)
        positions_unreadable = int(position_counts["unreadable"] or 0)
        positions_non_finite = int(position_counts["non_finite"] or 0)
        positions_non_integral = int(position_counts["non_integral"] or 0)
        positions_outside_int64 = int(position_counts["outside_int64"] or 0)
        expressions.append(
            canonical_position.alias(position)
        )

    for role in ("ea", "oa"):
        if role in resolved:
            expressions.append(
                canonical_allele_expression(
                    pl.col(resolved[role])
                ).alias(resolved[role])
            )

    normalized = frame.with_columns(expressions)
    positions_retained = 0
    if "pos" in resolved:
        positions_retained = count_non_null(normalized, resolved["pos"])
        if require_position_retention:
            validate_cast_retention(
                positions_present,
                positions_retained,
                "%s position" % label,
                error_type=error_type,
                warn=warn,
            )
    return normalized, {
        "positions_present": positions_present,
        "positions_retained": positions_retained,
        "positions_unreadable": positions_unreadable,
        "positions_non_finite": positions_non_finite,
        "positions_non_integral": positions_non_integral,
        "positions_outside_int64": positions_outside_int64,
    }


def read_reference_variant_table(
    path,
    column_mapping: Mapping[str, object],
    policies,
    *,
    value_columns: Sequence[str] = (),
    error_type: Type[Exception] = RuntimeError,
    description: str = "variant reference table",
):
    """Read and validate one allele-keyed harmonisation reference table.

    File parsing and required CHR/POS/A1/A2 schema validation are identical for
    build, strand, EAF and INFO references.  Source-specific value validation,
    allele semantics and duplicate policy remain with their scientific caller.
    Only columns needed by that caller are retained.
    """
    mapping_keys = ("chr", "pos", "a1", "a2", "delimiter")
    missing_mapping = [
        key
        for key in mapping_keys
        if column_mapping.get(key) is None
        or not str(column_mapping.get(key)).strip()
    ]
    if missing_mapping:
        raise error_type(
            "%s column mapping is missing: %s."
            % (description, ", ".join(missing_mapping))
        )
    required_columns = list(dict.fromkeys([
        *[str(column_mapping[key]) for key in ("chr", "pos", "a1", "a2")],
        *[str(column) for column in value_columns],
    ]))
    if Path(path).suffix.lower() == ".parquet":
        frame = read_parquet_table(
            path,
            columns=required_columns,
            error_type=error_type,
            description=description,
        )
        detected = DelimiterDetectionResult(
            "unknown", method="dataset-level projected Parquet partition",
        )
    else:
        frame, detected = read_delimited_table(
            path,
            str(column_mapping["delimiter"]),
            candidates=list(policies.get("input.delimiter_candidates")),
            minimum_columns=int(policies.get("input.delimiter_min_columns")),
            maximum_columns=int(policies.get("input.delimiter_max_columns")),
            sample_lines=int(policies.get("input.delimiter_sample_rows")),
            null_values=list(policies.get("input.null_values")),
            infer_schema_length=int(policies.get("input.schema_inference_rows")),
            columns=required_columns,
            error_type=error_type,
            description=description,
        )
    return frame, detected


def study_string_schema(
    column_mapping: Mapping[str, object],
) -> dict[str, pl.DataType]:
    """Keep variant labels and exact source provenance as strings."""
    return {
        str(column_mapping[key]): pl.String
        for key in _STRING_VARIANT_COLUMN_KEYS
        if column_mapping.get(key)
    }


def canonical_variant_schema(
    column_mapping: Mapping[str, object],
) -> dict[str, pl.DataType]:
    """Return lossless partition-read types for variant identity and position."""
    schema = study_string_schema(column_mapping)
    position = column_mapping.get("pos_col")
    if position:
        schema[str(position)] = pl.Int64
    return schema


def has_canonical_variant_columns(
    frame: pl.DataFrame,
    column_mapping: Mapping[str, object],
) -> bool:
    """Check the pipeline marker and exact dtypes without scanning any rows."""
    if column_mapping.get(CANONICAL_VARIANT_COLUMNS_KEY) is not True:
        return False
    expected = canonical_variant_schema(column_mapping)
    return bool(expected) and all(
        frame.schema.get(name) == dtype for name, dtype in expected.items()
    )


def mark_canonical_variant_columns(column_mapping: dict) -> None:
    """Record that normalization and scientific input validation succeeded."""
    column_mapping[CANONICAL_VARIANT_COLUMNS_KEY] = True


__all__ = [
    "CANONICAL_VARIANT_COLUMNS_KEY",
    "PALINDROMIC_SNP_PAIRS",
    "canonical_allele_expression",
    "canonical_chromosome_expression",
    "canonical_position_expression",
    "canonicalize_variant_frame",
    "canonical_variant_schema",
    "has_canonical_variant_columns",
    "mark_canonical_variant_columns",
    "minimal_allele_representation_series",
    "palindromic_snp_expression",
    "position_text_expression",
    "read_reference_variant_table",
    "reverse_complement_expression",
    "study_string_schema",
]
