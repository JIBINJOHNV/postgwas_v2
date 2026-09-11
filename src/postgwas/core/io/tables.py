"""Configuration-driven readers for headed delimited tables."""

from __future__ import annotations

import csv
import gzip
import re
import shutil
from pathlib import Path
from typing import Sequence, Type

import polars as pl

from .delimiters import DelimiterDetectionResult, open_text, resolve_delimiter


def read_pandas_table(
    path, delimiter: str | None, label: str, *, engine=None, error_type=ValueError,
):
    """Read a nonempty table using the existing pandas consumer semantics.

    Do not cache arbitrary full tables: outputs may legitimately change at a
    later stage. Read-only resource validators cache their compact evidence.
    """
    import pandas as pd

    try:
        table = pd.read_csv(path, sep=delimiter, engine=engine)
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise error_type("Cannot read %s %s: %s" % (label, path, exc)) from exc
    if table.empty:
        raise error_type("%s contains no data rows: %s" % (label, path))
    return table


def require_table_columns(table, columns, label: str, *, error_type=ValueError) -> None:
    """Require named columns without changing the consumer's parsed table."""
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise error_type("%s is missing required columns: %s" % (label, ", ".join(missing)))


def _projected_columns(
    path: str | Path,
    available_columns: Sequence[str],
    columns: Sequence[str] | None,
    *,
    error_type: Type[Exception],
    description: str,
) -> list[str] | None:
    """Validate a requested projection without reading any table rows."""
    if columns is None:
        return None
    requested = [str(column) for column in columns]
    if not requested:
        raise error_type("%s column projection must not be empty" % description)
    available = [str(column) for column in available_columns]
    absent = [column for column in requested if column not in available]
    if absent:
        raise error_type(
            "%s '%s' is missing configured column(s): %s. Available columns: %s."
            % (
                description,
                path,
                ", ".join(absent),
                ", ".join(available),
            )
        )
    return requested


def _delimited_header(path: str | Path, separator: str) -> list[str]:
    """Read only the first non-comment header row from a delimited table."""
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not line.strip() or line.startswith("##"):
                continue
            if separator == " ":
                columns = re.split(r"\s+", line.strip())
            else:
                columns = next(csv.reader([line], delimiter=separator))
            if columns:
                columns[0] = columns[0].lstrip("\ufeff")
            return columns
    raise ValueError("table contains no header row")


def read_delimited_table(
    path: str | Path,
    configured_delimiter: str,
    *,
    candidates: Sequence[str],
    minimum_columns: int,
    maximum_columns: int,
    sample_lines: int,
    null_values: Sequence[str],
    infer_schema_length: int,
    columns: Sequence[str] | None = None,
    has_header: bool = True,
    column_names: Sequence[str] | None = None,
    error_type: Type[Exception] = RuntimeError,
    description: str = "table",
) -> tuple[pl.DataFrame, DelimiterDetectionResult]:
    """Read one configured table and physically project requested columns."""
    try:
        if not has_header and columns is not None:
            raise error_type(
                "%s cannot project named columns from a headerless table"
                % description
            )
        configured_names = (
            None if column_names is None else [str(name) for name in column_names]
        )
        if has_header and configured_names is not None:
            raise error_type(
                "%s column_names is valid only for a headerless table" % description
            )
        if configured_names is not None and (
            not configured_names or len(configured_names) != len(set(configured_names))
        ):
            raise error_type(
                "%s column_names must contain one or more unique names" % description
            )
        detected = resolve_delimiter(
            path,
            configured_delimiter,
            candidates=candidates,
            minimum_columns=minimum_columns,
            maximum_columns=maximum_columns,
            sample_lines=sample_lines,
        )
        projected = (
            None
            if columns is None or not has_header
            else _projected_columns(
                path,
                _delimited_header(path, detected.value),
                columns,
                error_type=error_type,
                description=description,
            )
        )
        if detected.value == " ":
            try:
                import pandas as pd
            except ImportError as exc:
                raise error_type(
                    "%s uses whitespace runs, but pandas is not installed" % description
                ) from exc
            frame = pl.from_pandas(
                pd.read_csv(
                    path,
                    sep=r"\s+",
                    usecols=projected,
                    header=0 if has_header else None,
                    names=configured_names,
                    na_values=list(null_values),
                    keep_default_na=False,
                )
            )
        else:
            frame = pl.read_csv(
                path,
                separator=detected.value,
                columns=projected,
                has_header=has_header,
                new_columns=configured_names,
                null_values=list(null_values),
                infer_schema_length=infer_schema_length,
            )
        if projected is not None:
            frame = frame.select(projected)
    except error_type:
        raise
    except (
        OSError,
        UnicodeDecodeError,
        ValueError,
        csv.Error,
        pl.exceptions.PolarsError,
    ) as exc:
        raise error_type("Could not read %s '%s': %s" % (description, path, exc)) from exc
    return frame, detected


def read_parquet_table(
    path: str | Path,
    *,
    columns: Sequence[str] | None = None,
    error_type: Type[Exception] = RuntimeError,
    description: str = "table",
) -> pl.DataFrame:
    """Read Parquet metadata first, then physically project requested columns."""
    try:
        projected = (
            None
            if columns is None
            else _projected_columns(
                path,
                list(pl.read_parquet_schema(path)),
                columns,
                error_type=error_type,
                description=description,
            )
        )
        frame = pl.read_parquet(path, columns=projected)
    except error_type:
        raise
    except (OSError, ValueError, pl.exceptions.PolarsError) as exc:
        raise error_type(
            "Could not read %s '%s': %s" % (description, path, exc)
        ) from exc
    if projected is not None:
        frame = frame.select(projected)
    return frame


def write_dataframe_table(
    frame: pl.DataFrame,
    path: str | Path,
    *,
    overwrite: bool,
    runtime,
    include_header: bool = True,
    error_type: Type[Exception] = RuntimeError,
) -> str:
    """Write one configured dataframe atomically, with optional gzip output."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not overwrite:
        raise error_type(
            "Output already exists: %s. Enable overwrite or choose another "
            "output directory." % destination
        )
    temporary = destination.with_name(
        ".%s%s" % (destination.name, runtime.atomic_output_suffix)
    )
    uncompressed = destination.with_name(
        ".%s%s" % (destination.name, runtime.atomic_uncompressed_suffix)
    )
    try:
        if destination.suffix.lower() in runtime.compressed_suffixes:
            frame.write_csv(
                uncompressed,
                separator=runtime.table_delimiter,
                null_value=runtime.output_null_value,
                include_header=include_header,
            )
            with uncompressed.open("rb") as source, gzip.open(temporary, "wb") as target:
                shutil.copyfileobj(source, target, length=runtime.io_buffer_bytes)
        else:
            frame.write_csv(
                temporary,
                separator=runtime.table_delimiter,
                null_value=runtime.output_null_value,
                include_header=include_header,
            )
        temporary.replace(destination)
    except (OSError, pl.exceptions.PolarsError) as exc:
        raise error_type("Cannot write %s: %s" % (destination, exc)) from exc
    finally:
        temporary.unlink(missing_ok=True)
        uncompressed.unlink(missing_ok=True)
    return str(destination)


__all__ = [
    "read_pandas_table",
    "require_table_columns",
    "read_delimited_table",
    "read_parquet_table",
    "write_dataframe_table",
]
