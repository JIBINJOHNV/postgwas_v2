"""Deterministic, compressed-file-aware delimiter resolution."""

from __future__ import annotations

import bz2
import csv
import gzip
import io
import lzma
import zipfile
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Literal, Sequence, TextIO


NAMED_DELIMITERS = {
    "tab": "\t",
    "comma": ",",
    "semicolon": ";",
    "space": " ",
    "whitespace": " ",
    "pipe": "|",
}


@dataclass(frozen=True)
class DelimiterDetectionResult:
    """Resolved separator and the evidence used to select it."""

    kind: Literal["delimiter", "unknown"]
    value: str | None = None
    column_count: int = 0
    method: str = ""


def delimiter_character(name: str) -> str:
    """Resolve one validated configuration name to a parser character."""
    try:
        return NAMED_DELIMITERS[str(name)]
    except KeyError as exc:
        raise ValueError("Unknown configured delimiter name: %s" % name) from exc


@contextmanager
def open_text(path: str | Path) -> Iterator[TextIO]:
    """Open one plain or supported compressed text file without loading it."""
    source = Path(path)
    lowered = source.name.lower()
    if lowered.endswith(".zip"):
        with zipfile.ZipFile(source) as archive:
            members = [
                name for name in archive.namelist()
                if not name.endswith("/") and not name.startswith("__MACOSX/")
            ]
            if len(members) != 1:
                raise ValueError(
                    "ZIP archive must contain exactly one data file; found %d in %s"
                    % (len(members), source)
                )
            with archive.open(members[0]) as raw:
                with io.TextIOWrapper(
                    raw, encoding="utf-8", errors="replace",
                ) as text:
                    yield text
        return
    opener = open
    if lowered.endswith((".gz", ".bgz")):
        opener = gzip.open
    elif lowered.endswith(".bz2"):
        opener = bz2.open
    elif lowered.endswith((".xz", ".lzma")):
        opener = lzma.open
    with opener(source, "rt", encoding="utf-8", errors="replace") as handle:
        yield handle


def _sample_lines(
    path: str | Path,
    *,
    limit: int,
    comment_prefix: str | None,
) -> list[str]:
    lines = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not line.strip() or (
                comment_prefix and line.startswith(comment_prefix)
            ):
                continue
            lines.append(line)
            if len(lines) >= limit:
                break
    return lines


def detect_delimiter(
    path: str | Path,
    *,
    candidates: Sequence[str],
    minimum_columns: int,
    maximum_columns: int,
    sample_lines: int,
    comment_prefix: str | None = "##",
    fail_policy: Literal["raise", "unknown"] = "raise",
) -> DelimiterDetectionResult:
    """Select the candidate producing the widest consistent valid table."""
    if minimum_columns < 1 or maximum_columns < minimum_columns:
        raise ValueError("Invalid delimiter column-count bounds")
    if sample_lines < 1:
        raise ValueError("sample_lines must be greater than zero")
    separators = list(dict.fromkeys(
        delimiter_character(name) for name in candidates
    ))
    if not separators:
        raise ValueError("delimiter candidates must not be empty")
    lines = _sample_lines(
        path, limit=sample_lines, comment_prefix=comment_prefix,
    )
    best = None
    for separator in separators:
        counts = [
            len(next(csv.reader([line], delimiter=separator)))
            for line in lines
        ]
        if not counts or len(set(counts)) != 1:
            continue
        columns = counts[0]
        if minimum_columns <= columns <= maximum_columns:
            if best is None or columns > best.column_count:
                best = DelimiterDetectionResult(
                    "delimiter", separator, columns, "configured_candidates",
                )
    if best is not None:
        return best
    if fail_policy == "unknown":
        return DelimiterDetectionResult("unknown")
    raise ValueError(
        "No configured delimiter produced a consistent table with %d to %d columns: %s"
        % (minimum_columns, maximum_columns, path)
    )


def resolve_delimiter(
    path: str | Path,
    configured: str,
    *,
    candidates: Sequence[str],
    minimum_columns: int,
    maximum_columns: int,
    sample_lines: int,
    comment_prefix: str | None = "##",
) -> DelimiterDetectionResult:
    """Use an explicit configured delimiter or run deterministic detection."""
    if configured != "auto":
        return DelimiterDetectionResult(
            "delimiter", delimiter_character(configured), 0, "configured",
        )
    return detect_delimiter(
        path,
        candidates=candidates,
        minimum_columns=minimum_columns,
        maximum_columns=maximum_columns,
        sample_lines=sample_lines,
        comment_prefix=comment_prefix,
    )


__all__ = [
    "DelimiterDetectionResult",
    "NAMED_DELIMITERS",
    "delimiter_character",
    "detect_delimiter",
    "open_text",
    "resolve_delimiter",
]
