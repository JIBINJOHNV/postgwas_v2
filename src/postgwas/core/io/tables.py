"""Configuration-driven readers for headed delimited tables."""

from __future__ import annotations

import gzip
from pathlib import Path
import shutil
from typing import Sequence, Type

import polars as pl

from .delimiters import DelimiterDetectionResult, resolve_delimiter


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
    error_type: Type[Exception] = RuntimeError,
    description: str = "table",
) -> tuple[pl.DataFrame, DelimiterDetectionResult]:
    """Read one configured table, including files separated by whitespace runs."""
    try:
        detected = resolve_delimiter(
            path,
            configured_delimiter,
            candidates=candidates,
            minimum_columns=minimum_columns,
            maximum_columns=maximum_columns,
            sample_lines=sample_lines,
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
                    na_values=list(null_values),
                    keep_default_na=False,
                )
            )
        else:
            frame = pl.read_csv(
                path,
                separator=detected.value,
                has_header=True,
                null_values=list(null_values),
                infer_schema_length=infer_schema_length,
            )
    except error_type:
        raise
    except (OSError, UnicodeDecodeError, ValueError, pl.exceptions.PolarsError) as exc:
        raise error_type("Could not read %s '%s': %s" % (description, path, exc)) from exc
    return frame, detected


def write_dataframe_table(
    frame: pl.DataFrame,
    path: str | Path,
    *,
    overwrite: bool,
    runtime,
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
            )
            with uncompressed.open("rb") as source, gzip.open(temporary, "wb") as target:
                shutil.copyfileobj(source, target, length=runtime.io_buffer_bytes)
        else:
            frame.write_csv(
                temporary,
                separator=runtime.table_delimiter,
                null_value=runtime.output_null_value,
            )
        temporary.replace(destination)
    except (OSError, pl.exceptions.PolarsError) as exc:
        raise error_type("Cannot write %s: %s" % (destination, exc)) from exc
    finally:
        temporary.unlink(missing_ok=True)
        uncompressed.unlink(missing_ok=True)
    return str(destination)


__all__ = ["read_delimited_table", "write_dataframe_table"]
