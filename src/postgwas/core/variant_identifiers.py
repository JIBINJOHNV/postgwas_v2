"""Shared variant-identifier templates and PLINK BIM convention inspection."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from string import Formatter
from typing import Literal, Type


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
    """Scan every BIM row and require one formatter-supported ID convention."""
    path = Path(bim_file).expanduser().resolve()
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
    counts = {"variants": 0, "rsids": 0, "unique_ids": 0, "other_ids": 0}
    with path.open("r", encoding="utf-8") as handle:
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
            if not identifier:
                raise error_type(
                    "PLINK BIM contains an empty variant identifier at line %d"
                    % line_number
                )
            try:
                position = int(fields[indexes["position"]])
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

    if counts["variants"] == 0:
        raise error_type("PLINK BIM contains no variants: %s" % path)
    if counts["rsids"] == counts["variants"]:
        identifier_type: VariantIdentifierType = "rsid"
    elif counts["unique_ids"] == counts["variants"]:
        identifier_type = "unique"
    else:
        raise error_type(
            "PLINK BIM identifiers are mixed or unsupported: rsID=%d, unique=%d, "
            "other=%d. Automatic formatter selection requires one homogeneous "
            "rsID or configured unique-ID convention."
            % (counts["rsids"], counts["unique_ids"], counts["other_ids"])
        )
    return BimIdentifierSummary(
        identifier_type=identifier_type,
        variants=counts["variants"],
        rsids=counts["rsids"],
        unique_ids=counts["unique_ids"],
        other_ids=counts["other_ids"],
        bim_file=path,
    )


__all__ = [
    "BimIdentifierSummary",
    "IDENTIFIER_TEMPLATE_FIELDS",
    "VariantIdentifierType",
    "identifier_template_parts",
    "inspect_bim_identifier_type",
    "render_identifier",
]
