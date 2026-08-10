"""Lossless input and canonical internal types for study variant columns.

Reference tables are separate scientific inputs and must still be normalized
when read. String overrides apply from the first study-file read; the canonical
contract applies after dataset-level coordinate and allele validation.
"""

from __future__ import annotations

from typing import Mapping

import polars as pl


CANONICAL_VARIANT_COLUMNS_KEY = "_variant_columns_canonical"

_STRING_VARIANT_COLUMN_KEYS = (
    "chr_col",
    "chr_pos_col",
    "snp_id_col",
    "ea_col",
    "oa_col",
)


def study_string_schema(
    column_mapping: Mapping[str, object],
) -> dict[str, pl.DataType]:
    """Keep configured variant labels losslessly as strings from first read."""
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
    "canonical_variant_schema",
    "has_canonical_variant_columns",
    "mark_canonical_variant_columns",
    "study_string_schema",
]
