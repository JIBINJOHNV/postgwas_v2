"""Consistent normalization and display of optional scalar values."""

from __future__ import annotations

import math
from typing import Any


MISSING_TEXT = frozenset({
    "", "NA", "N/A", "NAN", "NONE", "NULL", "-", ".", "NOT_APPLICABLE",
})


def optional_text(value: Any) -> str | None:
    """Return stripped text, or ``None`` for standard missing placeholders."""
    if value is None:
        return None
    text = str(value).strip()
    return None if text.upper() in MISSING_TEXT else text


def parse_integer(value: Any) -> int | None:
    """Parse integers, decimals, scientific notation, and grouped numbers."""
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return None if not math.isfinite(value) else int(round(value))
    text = str(value).strip().replace(",", "").replace(" ", "").replace("_", "")
    if not text:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return None if not math.isfinite(number) else int(round(number))


def format_count(value: Any, *, missing: str = "unknown") -> str:
    """Format an integer-like value with thousands separators."""
    if value is None:
        return missing
    try:
        return "{:,}".format(int(value))
    except (TypeError, ValueError):
        return str(value)


def format_number(value: Any, pattern: str = "%.6g", *, missing: str = "not available") -> str:
    """Format a finite statistic without failing on null or NaN."""
    if value is None:
        return missing
    try:
        if not math.isfinite(float(value)):
            return missing
        return pattern % value
    except (TypeError, ValueError):
        return str(value)


def format_percentage(
    part: Any,
    whole: Any,
    *,
    places: int = 2,
    missing: str = "not available",
) -> str:
    """Format ``part / whole`` as a percentage, including the zero denominator."""
    try:
        denominator = float(whole)
        numerator = float(part)
    except (TypeError, ValueError):
        return missing
    if not math.isfinite(denominator) or not math.isfinite(numerator):
        return missing
    value = 0.0 if denominator == 0 else 100.0 * numerator / denominator
    return ("%%.%df%%%%" % places) % value


def format_fraction_percentage(
    fraction: Any,
    *,
    places: int = 2,
    missing: str = "not available",
) -> str:
    """Format a fraction on the 0-to-1 scale as a percentage."""
    return format_percentage(fraction, 1, places=places, missing=missing)


__all__ = [
    "MISSING_TEXT",
    "format_count",
    "format_fraction_percentage",
    "format_number",
    "format_percentage",
    "optional_text",
    "parse_integer",
]
