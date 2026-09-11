"""Reusable one-based gene-location validation, without annotation policies.

GCTA gene lists contain chromosome/start/end/gene fields; PostGWAS MAGMA
locations additionally require strand and permit an alternate identifier.
This reader preserves those existing contracts, not BED's zero-based format.
Gene windows, chromosome normalization, aliases and reference compatibility
remain the responsibility of their scientific consumers.

https://yanglab.westlake.edu.cn/software/gcta/#fastBAT
https://ctg.cncr.nl/software/magma
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
import re

from postgwas.core.input_validation import record_file_validation, validate_once
from postgwas.core.io.delimiters import open_text
from postgwas.core.paths import require_nonempty_file


@dataclass(frozen=True)
class GeneCoordinate:
    gene: str
    chromosome: str
    start: int
    end: int
    strand: str | None = None
    alternate: str | None = None


def read_gene_coordinates(
    path: str | Path,
    *,
    column_roles: Sequence[str],
    delimiter_pattern: str | None = None,
    has_header: bool = False,
    comment_prefix: str | None = None,
    compressed: bool = False,
    label: str = "Gene list",
    error_type: type[Exception] = ValueError,
) -> tuple[GeneCoordinate, ...]:
    """Read immutable gene rows once per exact format and unchanged file.

    Only gene-sized data are cached. No SNP identifiers or GWAS frames are
    retained here. The optional alias is a trailing protocol field; it is not
    used to infer, replace or combine primary identifiers.
    """
    source = require_nonempty_file(path, label, error_type=error_type)
    roles = tuple(column_roles)
    required = {"gene", "chromosome", "start", "end"}
    optional = {"strand", "alternate_gene_id"}
    if len(set(roles)) != len(roles) or not required <= set(roles) or set(roles) - required - optional:
        raise error_type("Gene-coordinate column roles must contain gene, chromosome, start and end.")
    indexes = {role: roles.index(role) for role in roles}
    alias_index = indexes.get("alternate_gene_id")
    if alias_index is not None and alias_index != len(roles) - 1:
        raise error_type("The optional alternate gene identifier must be the last configured column.")
    counts = {len(roles)}
    if alias_index is not None:
        counts.add(len(roles) - 1)
    # Historical MAGMA input permits an unused sixth field when only the five
    # primary roles are configured; do not reinterpret that field as an alias.
    if "strand" in indexes and alias_index is None:
        counts.add(len(roles) + 1)

    def inspect():
        genes: set[str] = set()
        rows: list[GeneCoordinate] = []
        delimiter = re.compile(delimiter_pattern) if delimiter_pattern is not None else None
        header_pending = has_header
        try:
            opener = open_text(source) if compressed else source.open("r", encoding="utf-8")
            with opener as handle:
                for number, raw in enumerate(handle, 1):
                    text = raw.strip()
                    if not text or (comment_prefix and text.startswith(comment_prefix)):
                        continue
                    if header_pending:
                        header_pending = False
                        continue
                    fields = delimiter.split(text) if delimiter is not None else text.split()
                    if len(fields) not in counts:
                        expected = (
                            "five columns plus an optional alternate gene identifier"
                            if "strand" in indexes else "%d columns" % len(roles)
                        )
                        raise ValueError("%s line %d must contain %s; found %d" % (
                            label, number, expected, len(fields),
                        ))
                    gene = fields[indexes["gene"]].strip()
                    chromosome = fields[indexes["chromosome"]].strip()
                    try:
                        coordinates = [fields[indexes[role]].strip() for role in ("start", "end")]
                        # External integer fields must not accept Python-only
                        # literals such as 1_000 or non-ASCII decimal digits.
                        if any(re.fullmatch(r"[+-]?[0-9]+", value) is None for value in coordinates):
                            raise ValueError("coordinates must use decimal integer syntax")
                        start, end = map(int, coordinates)
                    except ValueError as exc:
                        raise ValueError("%s line %d requires integer start and end coordinates" % (
                            label, number,
                        )) from exc
                    strand = fields[indexes["strand"]].strip() if "strand" in indexes else None
                    if not gene or not chromosome or start < 1 or end < start or (
                        strand is not None and strand not in {"+", "-"}
                    ):
                        raise ValueError("%s line %d has an invalid gene ID, chromosome, interval, or strand" % (
                            label, number,
                        ))
                    if gene in genes:
                        raise ValueError("%s contains duplicate gene ID %s at line %d" % (label, gene, number))
                    alternate = (
                        fields[alias_index].strip() or None
                        if alias_index is not None and alias_index < len(fields) else None
                    )
                    genes.add(gene)
                    rows.append(GeneCoordinate(gene, chromosome, start, end, strand, alternate))
            if not rows:
                raise ValueError("%s contains no genes: %s" % (label, source))
        except (OSError, UnicodeError, ValueError) as exc:
            record_file_validation(source, "Gene coordinates", status="failed", message=str(exc))
            raise error_type(str(exc)) from exc
        record_file_validation(
            source, "Gene coordinates",
            checks=("all-row configured fields", "integer one-based intervals", "unique primary gene IDs"),
            metrics={"genes": len(rows), "chromosomes": sorted({row.chromosome for row in rows})},
            message="Coordinates preserved; reference-build and annotation policies are separate checks.",
        )
        return tuple(rows)

    return validate_once((source,), {
        "validator": "gene_coordinate_rows", "version": 1,
        "column_roles": roles, "delimiter_pattern": delimiter_pattern,
        "has_header": has_header, "comment_prefix": comment_prefix,
        "compressed": compressed,
    }, inspect, error_type=error_type)


__all__ = ["GeneCoordinate", "read_gene_coordinates"]
