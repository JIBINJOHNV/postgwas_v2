"""Read-only validation of GCTA's set/SNP/END block-list protocol.

https://yanglab.westlake.edu.cn/software/gcta/#fastBAT
Counts and exact duplicate checks use a temporary disk-backed SQLite table,
not a genome-wide Python SNP set. Membership filtering is a consumer task.
"""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
import sqlite3

from postgwas.core.input_validation import (
    capture_preflight_file_identities,
    record_file_validation,
    require_unchanged_preflight_files,
    validate_once,
)
from postgwas.core.paths import require_nonempty_file


def validate_fastbat_set_list(path: str | Path, *, error_type=ValueError) -> dict[str, int]:
    """Validate the static block structure once, retaining only compact counts."""
    source = require_nonempty_file(path, "fastBAT set list", error_type=error_type)

    def inspect():
        names: set[str] = set()
        memberships = 0

        def rows():
            nonlocal memberships
            current = None
            current_count = 0
            with source.open("r", encoding="utf-8") as handle:
                for number, raw in enumerate(handle, 1):
                    value = raw.strip()
                    if not value:
                        continue
                    if current is None:
                        if value.upper() == "END":
                            raise ValueError("fastBAT set list has END without a set ID at line %d." % number)
                        if len(value.split()) != 1:
                            raise ValueError("fastBAT set ID at line %d must not contain whitespace." % number)
                        if value in names:
                            raise ValueError("fastBAT set list repeats set ID: %s" % value)
                        names.add(value)
                        current = value
                        current_count = 0
                    elif value.upper() == "END":
                        if not current_count:
                            raise ValueError("fastBAT set %s contains no variant IDs." % current)
                        current = None
                    else:
                        if len(value.split()) != 1:
                            raise ValueError("fastBAT set-list line %d must contain one variant ID." % number)
                        current_count += 1
                        memberships += 1
                        yield current, value
                if current is not None:
                    raise ValueError("fastBAT set %s is missing its terminating END line." % current)
            if not names:
                raise ValueError("fastBAT set list contains no sets: %s" % source)

        connection = sqlite3.connect("")
        try:
            connection.execute("CREATE TABLE membership (set_id TEXT, variant_id TEXT)")
            # Bulk ingestion avoids per-row index maintenance. The empty SQLite
            # filename is temporary disk-backed storage, removed on close.
            connection.executemany("INSERT INTO membership VALUES (?, ?)", rows())
            duplicate = connection.execute(
                "SELECT set_id, variant_id FROM membership GROUP BY set_id, variant_id "
                "HAVING COUNT(*) > 1 LIMIT 1"
            ).fetchone()
            if duplicate is not None:
                raise ValueError("fastBAT set %s repeats variant ID %s." % duplicate)
            unique_variants = connection.execute(
                "SELECT COUNT(DISTINCT variant_id) FROM membership"
            ).fetchone()[0]
        except (OSError, UnicodeError, ValueError, sqlite3.Error) as exc:
            record_file_validation(source, "fastBAT SNP sets", status="failed", message=str(exc))
            raise error_type(str(exc)) from exc
        finally:
            connection.close()
        result = {
            "input_sets": len(names),
            "requested_set_variants": memberships,
            "unique_requested_set_variants": unique_variants,
        }
        record_file_validation(
            source, "fastBAT SNP sets",
            checks=("one token per line", "unique nonempty sets", "END terminators", "no duplicate membership within a set"),
            metrics=result,
            message="GWAS/reference membership and configured set-size policies remain separate checks.",
        )
        return result

    return dict(validate_once(
        (source,), {"validator": "fastbat_set_list", "version": 1},
        inspect, error_type=error_type,
    ))


@contextmanager
def open_fastbat_set_memberships(path: str | Path, *, error_type=ValueError):
    """Yield validated counts and streaming membership/end events for analysis.

    The file is read again to consume its data, not to repeat the static checks.
    An event's second value is ``None`` at the end of its set. Identity checks
    protect the interval between validation and consumption; consumers should
    publish scientific outputs only after this context exits successfully.
    """
    source = Path(path).expanduser().resolve()
    identity = capture_preflight_file_identities((source,), error_type=error_type)
    metrics = validate_fastbat_set_list(source, error_type=error_type)

    def events(handle):
        current = None
        for raw in handle:
            value = raw.strip()
            if not value:
                continue
            if current is None:
                current = value
            elif value.upper() == "END":
                yield current, None
                current = None
            else:
                yield current, value

    with source.open("r", encoding="utf-8") as handle:
        yield metrics, events(handle)
    require_unchanged_preflight_files(identity, error_type=error_type)


__all__ = ["open_fastbat_set_memberships", "validate_fastbat_set_list"]
