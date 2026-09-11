"""Shared variant-identifier templates and PLINK BIM convention inspection."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from contextlib import closing
from pathlib import Path
import re
import sqlite3
from string import Formatter
from typing import Literal, Type

from postgwas.core.input_validation import record_file_validation, validate_once
from postgwas.core.values import missing_tokens


VariantIdentifierType = Literal["rsid", "unique"]
# These names define the scientific components of PostGWAS unique variant IDs.
IDENTIFIER_TEMPLATE_FIELDS = (
    "chromosome",
    "position",
    "reference_allele",
    "alternate_allele",
)


@dataclass(frozen=True)
class BimIdentifierSummary:
    """Complete counts supporting one homogeneous BIM identifier decision."""

    identifier_type: VariantIdentifierType
    variants: int
    rsids: int
    unique_ids: int
    other_ids: int
    bim_file: Path
    missing_ids: int = 0
    duplicate_ids: int = 0
    duplicate_rows: int = 0


def identifier_template_parts(template: str):
    """Return validated templates as reusable literal/field components."""
    return tuple(Formatter().parse(template))


def render_identifier(template, values: dict[str, object]) -> str:
    """Render one configured identifier without embedding a naming convention."""
    pieces = []
    parts = (
        identifier_template_parts(template)
        if isinstance(template, str)
        else template
    )
    for literal, field, _, _ in parts:
        pieces.append(literal)
        if field is not None:
            pieces.append(str(values[field]))
    return "".join(pieces)


def inspect_bim_identifier_type(
    bim_file: str | Path,
    *,
    column_roles: list[str],
    delimiter_pattern: str,
    rsid_pattern: str,
    unique_id_template: str,
    chromosome_prefix_pattern: str,
    chromosome_aliases: dict[str, str],
    error_type: Type[Exception] = ValueError,
) -> BimIdentifierSummary:
    """Inspect once per unchanged file and exact parsing/identifier contract.

    Duplicate counts do not impose a new universal rejection policy: consumers
    retain their scientifically appropriate full-panel or matched-subset rule.
    BIM A1/A2 are not guaranteed REF/ALT; both token orders are recognised.
    """
    path = Path(bim_file).expanduser().resolve()
    options = {
        "column_roles": column_roles,
        "delimiter_pattern": delimiter_pattern,
        "rsid_pattern": rsid_pattern,
        "unique_id_template": unique_id_template,
        "chromosome_prefix_pattern": chromosome_prefix_pattern,
        "chromosome_aliases": chromosome_aliases,
    }

    def inspect():
        try:
            return _scan_bim_identifier_type(path, **options, error_type=error_type)
        except (OSError, UnicodeError) as exc:
            record_file_validation(path, "PLINK BIM", status="failed", message=str(exc))
            raise error_type("Cannot read PLINK BIM %s: %s" % (path, exc)) from exc

    return validate_once(
        (path,), {"validator": "bim_identifier_convention", "version": 2, **options}, inspect,
    )


def _scan_bim_identifier_type(
    path: Path,
    *,
    column_roles: list[str],
    delimiter_pattern: str,
    rsid_pattern: str,
    unique_id_template: str,
    chromosome_prefix_pattern: str,
    chromosome_aliases: dict[str, str],
    error_type: Type[Exception],
) -> BimIdentifierSummary:
    """Stream BIM rows into a temporary on-disk ID tally with bounded memory."""
    if not path.is_file() or path.stat().st_size <= 0:
        raise error_type("PLINK BIM file does not exist or is empty: %s" % path)
    required = {
        "chromosome", "variant_id", "position", "allele1", "allele2",
    }
    if not required.issubset(column_roles):
        raise error_type(
            "BIM column roles are missing: %s"
            % ", ".join(sorted(required - set(column_roles)))
        )
    indexes = {role: column_roles.index(role) for role in required}
    try:
        delimiter = re.compile(delimiter_pattern)
        rsid = re.compile(rsid_pattern)
        chromosome_prefix = re.compile(chromosome_prefix_pattern)
    except re.error as exc:
        raise error_type(
            "Invalid configured BIM delimiter, rsID, or chromosome pattern: %s"
            % exc
        ) from exc

    template = identifier_template_parts(unique_id_template)
    counts = {"variants": 0, "rsids": 0, "unique_ids": 0, "other_ids": 0, "missing_ids": 0}
    null_identifiers = frozenset(missing_tokens())

    def identifiers(handle):
        for line_number, raw in enumerate(handle, 1):
            text = raw.strip()
            if not text:
                raise error_type(
                    "PLINK BIM contains a blank row at line %d" % line_number
                )
            fields = delimiter.split(text)
            if len(fields) != len(column_roles):
                raise error_type(
                    "PLINK BIM line %d has %d fields; expected %d"
                    % (line_number, len(fields), len(column_roles))
                )
            identifier = fields[indexes["variant_id"]].strip()
            position_text = fields[indexes["position"]]
            # BIM coordinates are decimal text, not Python numeric literals.
            # int() alone also accepts underscores and non-ASCII digits, which
            # can disagree with the coordinate parsed by a PLINK consumer.
            if not position_text.isascii() or not position_text.isdecimal():
                raise error_type(
                    "PLINK BIM contains an invalid position at line %d; "
                    "expected ASCII decimal digits" % line_number
                )
            try:
                position = int(position_text)
            except ValueError as exc:
                raise error_type(
                    "PLINK BIM contains an invalid position at line %d" % line_number
                ) from exc
            if position < 1:
                raise error_type(
                    "PLINK BIM contains a non-positive position at line %d"
                    % line_number
                )
            chromosome = fields[indexes["chromosome"]].strip()
            allele1 = fields[indexes["allele1"]].strip().upper()
            allele2 = fields[indexes["allele2"]].strip().upper()
            if not chromosome or not allele1 or not allele2:
                raise error_type(
                    "PLINK BIM contains an empty chromosome or allele at line %d"
                    % line_number
                )

            counts["variants"] += 1
            if identifier.upper() in null_identifiers:
                counts["missing_ids"] += 1
                continue
            yield (identifier,)
            if rsid.fullmatch(identifier):
                counts["rsids"] += 1
                continue
            normalized_chromosome = chromosome_prefix.sub("", chromosome)
            normalized_chromosome = chromosome_aliases.get(
                normalized_chromosome, normalized_chromosome,
            )
            common = {
                "chromosome": normalized_chromosome,
                "position": position,
            }
            first = render_identifier(
                template,
                {
                    **common,
                    "reference_allele": allele1,
                    "alternate_allele": allele2,
                },
            )
            second = render_identifier(
                template,
                {
                    **common,
                    "reference_allele": allele2,
                    "alternate_allele": allele1,
                },
            )
            if identifier in {first, second}:
                counts["unique_ids"] += 1
            else:
                counts["other_ids"] += 1

    # An empty SQLite database name creates an automatically removed temporary
    # disk database. Bulk insertion followed by grouping avoids an in-memory ID
    # set and per-row primary-key/index maintenance on genome-wide references.
    with closing(sqlite3.connect("")) as connection:
        connection.execute("CREATE TABLE identifiers (identifier TEXT NOT NULL)")
        with path.open("r", encoding="utf-8") as handle:
            connection.executemany("INSERT INTO identifiers VALUES (?)", identifiers(handle))
        duplicate_ids, duplicate_rows = connection.execute(
            "SELECT COUNT(*), COALESCE(SUM(n - 1), 0) FROM "
            "(SELECT COUNT(*) AS n FROM identifiers GROUP BY identifier HAVING COUNT(*) > 1)"
        ).fetchone()
    counts.update(duplicate_ids=duplicate_ids, duplicate_rows=duplicate_rows)
    if counts["variants"] == 0:
        raise error_type("PLINK BIM contains no variants: %s" % path)
    if counts["rsids"] == counts["variants"]:
        identifier_type: VariantIdentifierType = "rsid"
    elif counts["unique_ids"] == counts["variants"]:
        identifier_type = "unique"
    else:
        message = (
            "PLINK BIM identifiers are mixed or unsupported: rsID=%d, unique=%d, "
            "other=%d, missing=%d, duplicated IDs=%d. Automatic formatter selection requires one homogeneous "
            "rsID or configured unique-ID convention."
            % (counts["rsids"], counts["unique_ids"], counts["other_ids"],
               counts["missing_ids"], counts["duplicate_ids"])
        )
        record_file_validation(
            path, "PLINK BIM", checks=("all-row identifier convention",
                                      "coordinate-style IDs agree with BIM coordinates/alleles where applicable"),
            metrics=counts, status="failed", message=message,
        )
        raise error_type(message)
    summary = BimIdentifierSummary(
        identifier_type=identifier_type,
        variants=counts["variants"],
        rsids=counts["rsids"],
        unique_ids=counts["unique_ids"],
        other_ids=counts["other_ids"],
        bim_file=path,
        missing_ids=counts["missing_ids"],
        duplicate_ids=counts["duplicate_ids"],
        duplicate_rows=counts["duplicate_rows"],
    )
    record_file_validation(
        path, "PLINK BIM",
        checks=("all-row field counts", "positive integer positions", "identifier convention",
                "coordinate-style IDs agree with BIM coordinates/alleles where applicable",
                "missing and duplicate ID counts"),
        metrics={key: value for key, value in asdict(summary).items() if key != "bim_file"},
        status="warning" if duplicate_ids else "passed",
        message=("Duplicate IDs require the consuming method's compatibility policy. " if duplicate_ids else "")
        + "BIM allele order does not establish reference-genome or effect orientation.",
    )
    return summary


__all__ = [
    "BimIdentifierSummary",
    "IDENTIFIER_TEMPLATE_FIELDS",
    "VariantIdentifierType",
    "identifier_template_parts",
    "inspect_bim_identifier_type",
    "render_identifier",
]
