"""Streaming validation of gzip BED4 annotation resources.

BED's first three columns use integer, zero-based, half-open coordinates;
zero-length intervals are legal. These are format invariants, not tunable
scientific filters (https://genome.ucsc.edu/FAQ/FAQformat.html#format1).
The fourth column is required here because consumers request an annotation
label. No chromosome normalization or interval filtering is performed.
"""

from __future__ import annotations

import gzip
from collections.abc import Collection, Sequence
from dataclasses import dataclass
from pathlib import Path
import re
import zlib

from postgwas.core.input_validation import record_file_validation, validate_once
from postgwas.core.paths import require_nonempty_file


@dataclass(frozen=True)
class Bed4Validation:
    """Compact all-record inventory, without retaining every interval."""

    path: Path
    contigs: frozenset[str]
    labels: frozenset[str]
    block_count: int


def validate_bed4_file(path: str | Path) -> Bed4Validation:
    """Validate every gzip BED4 row once per unchanged pipeline resource."""
    bed = Path(path).expanduser().resolve()

    def inspect() -> Bed4Validation:
        try:
            require_nonempty_file(bed, "BED4 annotation reference")
            contigs: set[str] = set()
            labels: set[str] = set()
            count = 0
            with gzip.open(bed, "rt", encoding="utf-8") as handle:
                for line_number, line in enumerate(handle, start=1):
                    if not line.strip() or line.startswith("#"):
                        continue
                    fields = line.rstrip("\r\n").split("\t")
                    if len(fields) < 4 or not fields[0] or not fields[3]:
                        raise ValueError(
                            "LD-block BED row %d in %s must contain CHROM, START, "
                            "END, and a nonempty annotation label" % (line_number, bed)
                        )
                    try:
                        if any(re.fullmatch(r"[+-]?[0-9]+", value) is None for value in fields[1:3]):
                            raise ValueError("coordinates must use decimal integer syntax")
                        start, end = int(fields[1]), int(fields[2])
                    except ValueError as exc:
                        raise ValueError(
                            "BED row %d in %s has non-integer START or END"
                            % (line_number, bed)
                        ) from exc
                    if start < 0 or end < start:
                        raise ValueError(
                            "BED row %d in %s must have 0 <= START <= END "
                            "(zero-based, half-open coordinates)" % (line_number, bed)
                        )
                    contigs.add(fields[0])
                    labels.add(fields[3])
                    count += 1
            if not count:
                raise ValueError("LD-block BED contains no annotation rows: %s" % bed)
        except (OSError, EOFError, UnicodeError, ValueError, zlib.error) as exc:
            record_file_validation(
                bed, "BED4 annotation reference", status="failed", message=str(exc),
            )
            raise ValueError("Invalid gzip BED4 annotation reference %s: %s" % (bed, exc)) from exc
        result = Bed4Validation(bed, frozenset(contigs), frozenset(labels), count)
        record_file_validation(
            bed,
            "BED4 annotation reference",
            checks=(
                "gzip text read through EOF", "all rows: required BED4 fields",
                "all rows: integer zero-based half-open coordinates",
                "nonempty annotation labels",
            ),
            metrics={
                "rows": count, "contigs": sorted(contigs),
                "unique_labels": len(labels), "duplicate_labels": count - len(labels),
            },
        )
        return result

    return validate_once((bed,), ("gzip_bed4", "zero_based_half_open"), inspect)


def _contig_aliases(contigs: Collection[str]) -> dict[str, set[str]]:
    """Group exact labels only to diagnose optional ``chr``-prefix mismatches."""
    aliases: dict[str, set[str]] = {}
    for contig in contigs:
        identity = contig[3:] if contig.lower().startswith("chr") else contig
        aliases.setdefault(identity, set()).add(contig)
    return aliases


def validate_vcf_bed_contigs(
    vcf_contigs: Sequence[str],
    bed_contigs: Collection[str],
    bed: Path,
) -> None:
    """Require exact VCF/BED labels; never silently add or remove ``chr``."""
    vcf_aliases = _contig_aliases(vcf_contigs)
    bed_aliases = _contig_aliases(bed_contigs)
    shared_identities = set(vcf_aliases).intersection(bed_aliases)
    mismatched = sorted(
        identity
        for identity in shared_identities
        if vcf_aliases[identity] != bed_aliases[identity]
    )
    if set(vcf_contigs).intersection(bed_contigs) and not mismatched:
        return

    examples = ", ".join(
        "%s versus %s"
        % (
            "/".join(sorted(vcf_aliases[identity])),
            "/".join(sorted(bed_aliases[identity])),
        )
        for identity in mismatched[:3]
    )
    detail = " Mismatched examples: %s." % examples if examples else ""
    raise ValueError(
        "Chromosome names are incompatible between the input VCF and LD-block "
        "BED %s.%s VCF contigs include: %s. BED contigs include: %s. Names "
        "must match exactly; PostGWAS does not add or remove the 'chr' prefix."
        % (
            bed,
            detail,
            ", ".join(vcf_contigs[:5]),
            ", ".join(sorted(bed_contigs)[:5]),
        )
    )


__all__ = ["Bed4Validation", "validate_bed4_file", "validate_vcf_bed_contigs"]
