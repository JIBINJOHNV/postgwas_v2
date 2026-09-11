"""Streaming GMT/native gene-set validation and immutable file evidence.

GMT's tab-separated name, description and gene fields are format invariants:
https://docs.gsea-msigdb.org/GSEA/Data_Formats/#gmt-gene-matrix-transposed-file-format-gmt
Consumers explicitly supply their existing missing/duplicate/name policies.
Only the gene universe and counts are retained, never all set memberships.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator
from dataclasses import asdict, dataclass, replace
from pathlib import Path
import re

from postgwas.core.input_validation import (
    PreflightFileIdentity,
    capture_preflight_file_identities,
    record_file_validation,
    require_unchanged_preflight_files,
    validate_once,
)
from postgwas.core.io.delimiters import open_text
from postgwas.core.paths import require_nonempty_file


@dataclass(frozen=True)
class GeneSetFormat:
    input_format: str
    delimiter_pattern: str | None
    comment_prefix: str | None
    compressed: bool
    empty_gene_policy: str
    duplicate_gene_policy: str
    reject_name_whitespace: bool
    reserved_names: tuple[str, ...] = ()
    membership_has_header: bool = False
    membership_set_column: int | None = None
    membership_gene_column: int | None = None


@dataclass(frozen=True)
class GeneSetRow:
    index: int
    name: str
    description: str
    genes: tuple[str, ...]
    duplicate_gene_entries: int
    input_format: str


@dataclass(frozen=True)
class GeneSetSource:
    path: Path
    options: GeneSetFormat
    identities: tuple[PreflightFileIdentity, ...]
    gene_sets: int
    genes: frozenset[str]
    detected_format: str

    def rows(self, *, error_type=ValueError) -> Iterator[GeneSetRow]:
        """Consume memberships in source order without repeating static checks."""
        require_unchanged_preflight_files(self.identities, error_type=error_type)
        options = replace(self.options, input_format=self.detected_format)
        try:
            yield from _iter_gene_set_rows(self.path, options, validate=False)
        except (OSError, UnicodeError, ValueError) as exc:
            raise error_type(str(exc)) from exc
        require_unchanged_preflight_files(self.identities, error_type=error_type)


def _iter_gene_set_rows(
    source: Path, options: GeneSetFormat, *, validate: bool,
) -> Iterator[GeneSetRow]:
    """Shared decoder; validation-only guards are disabled for proven files."""
    if options.input_format == "membership":
        yield from _iter_gene_membership_rows(source, options, validate=validate)
        return
    detected_format = None
    delimiter = re.compile(options.delimiter_pattern) if options.delimiter_pattern else None
    seen: set[str] = set()
    index = 0
    opener = open_text(source) if options.compressed else source.open("r", encoding="utf-8")
    with opener as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or (
                options.comment_prefix and raw.lstrip().startswith(options.comment_prefix)
            ):
                continue
            text = raw.rstrip("\r\n")
            fields = text.split("\t")
            line_format = "gmt" if len(fields) >= 3 else "magma"
            if detected_format is None:
                detected_format = line_format if options.input_format == "auto" else options.input_format
            if validate and options.input_format == "auto" and line_format != detected_format:
                raise ValueError("Gene-set file mixes GMT and native MAGMA records at line %d" % number)
            if detected_format == "gmt":
                if validate and (len(fields) < 3 or not fields[0].strip()):
                    raise ValueError("Malformed GMT gene-set line %d: expected name, description and genes" % number)
                name, description = fields[0].strip(), fields[1].strip()
                values = fields[2:]
            else:
                fields = delimiter.split(text) if delimiter is not None else text.split()
                if validate and (len(fields) < 2 or not fields[0].strip()):
                    raise ValueError("Malformed native MAGMA gene-set line %d" % number)
                name, description = fields[0].strip(), ""
                values = fields[1:]
            genes = [value.strip() for value in values]
            if validate:
                if options.reject_name_whitespace and any(character.isspace() for character in name):
                    raise ValueError("GMT line %d has a whitespace-containing pathway ID" % number)
                if name in options.reserved_names:
                    raise ValueError("GMT pathway ID %s is reserved by the consumer" % name)
                if name in seen:
                    raise ValueError("Duplicate gene-set name: %s" % name)
                if options.empty_gene_policy == "reject" and any(not gene for gene in genes):
                    raise ValueError("GMT pathway %r contains an empty gene identifier" % name)
                seen.add(name)
            if options.empty_gene_policy == "drop":
                genes = [gene for gene in genes if gene]
            unique = tuple(dict.fromkeys(genes))
            duplicate_count = len(genes) - len(unique)
            if validate and not unique:
                raise ValueError("Gene set %s contains no valid gene IDs" % name)
            if validate and duplicate_count and options.duplicate_gene_policy == "error":
                raise ValueError("GMT pathway %r contains %d duplicate gene entries" % (name, duplicate_count))
            yield GeneSetRow(index, name, description, unique, duplicate_count, detected_format)
            index += 1
    if validate and not index:
        raise ValueError("Gene-set file contains no gene sets: %s" % source)


def _iter_gene_membership_rows(
    source: Path, options: GeneSetFormat, *, validate: bool,
) -> Iterator[GeneSetRow]:
    """Group the configured column-based gene memberships, not SNP sets.

    Grouping is required to create the consumer's gene-set rows. It is not
    retained in validation evidence, and repeated copies of each membership
    are not accumulated when deduplicating an unsorted input table.
    """
    delimiter = re.compile(options.delimiter_pattern) if options.delimiter_pattern else None
    grouped: dict[str, dict[str, None]] = {}
    duplicate_counts: dict[str, int] = {}
    maximum = max(options.membership_set_column, options.membership_gene_column)
    opener = open_text(source) if options.compressed else source.open("r", encoding="utf-8")
    with opener as handle:
        for number, raw in enumerate(handle, 1):
            if options.membership_has_header and number == 1:
                continue
            text = raw.strip()
            if not text or (options.comment_prefix and text.startswith(options.comment_prefix)):
                continue
            fields = delimiter.split(text) if delimiter is not None else text.split()
            if validate and len(fields) <= maximum:
                raise ValueError("Malformed two-column gene-set membership at line %d" % number)
            name = fields[options.membership_set_column].strip()
            gene = fields[options.membership_gene_column].strip()
            if validate and (not name or not gene):
                raise ValueError("Empty set or gene identifier in membership line %d" % number)
            genes = grouped.setdefault(name, {})
            if gene in genes:
                duplicate_counts[name] = duplicate_counts.get(name, 0) + 1
            genes[gene] = None
    if validate and not grouped:
        raise ValueError("Gene-set file contains no gene sets: %s" % source)
    for index, (name, genes) in enumerate(grouped.items()):
        yield GeneSetRow(index, name, "", tuple(genes), duplicate_counts.get(name, 0), "membership")


def validate_gene_set_source(
    path: str | Path, options: GeneSetFormat, *, error_type=ValueError,
    evidence: GeneSetSource | None = None,
    on_row: Callable[[GeneSetRow], None] | None = None,
) -> GeneSetSource:
    """Scan once for an exact contract and cache compact, read-only evidence.

    ``on_row`` may construct the caller's already-required in-memory output
    during first validation. Cache hits stream data without rechecking static
    rules. The callback must not publish outputs before this function returns:
    the complete read and final file-identity check must succeed first.
    """
    source = require_nonempty_file(path, "gene-set file", error_type=error_type)
    if options.input_format not in {"auto", "gmt", "magma", "membership"}:
        raise error_type("Unsupported gene-set input format: %s" % options.input_format)
    if options.empty_gene_policy not in {"reject", "drop"} or options.duplicate_gene_policy not in {"error", "deduplicate"}:
        raise error_type("Invalid gene-set missing/duplicate policy")
    if options.input_format == "membership":
        indexes = (options.membership_set_column, options.membership_gene_column)
        if any(not isinstance(index, int) or index < 0 for index in indexes) or indexes[0] == indexes[1]:
            raise error_type("Gene membership requires distinct nonnegative configured set/gene columns")
        if (options.empty_gene_policy != "reject" or options.duplicate_gene_policy != "deduplicate"
                or options.reject_name_whitespace or options.reserved_names):
            raise error_type("Column-based membership requires its native nonempty/deduplication/name contract")
    if evidence is not None:
        if evidence.path != source or evidence.options != options:
            raise error_type("Validated gene-set evidence does not match the requested file and format policies")
        require_unchanged_preflight_files(evidence.identities, error_type=error_type)
    consumed = False

    def inspect():
        nonlocal consumed
        identities = capture_preflight_file_identities((source,), error_type=error_type)
        genes: set[str] = set()
        count = 0
        duplicates = 0
        detected_format = None
        try:
            for row in _iter_gene_set_rows(source, options, validate=True):
                genes.update(row.genes)
                count += 1
                duplicates += row.duplicate_gene_entries
                detected_format = row.input_format
                if on_row is not None:
                    on_row(row)
            require_unchanged_preflight_files(identities, error_type=error_type)
        except (OSError, UnicodeError, ValueError) as exc:
            record_file_validation(source, "Gene-set reference", status="failed", message=str(exc))
            raise error_type(str(exc)) from exc
        record_file_validation(
            source, "Gene-set reference",
            checks=("all-row format",
                    "set-name grouping" if detected_format == "membership" else "unique set names",
                    "configured empty/duplicate gene policies"),
            metrics={"gene_sets": count, "unique_genes": len(genes), "detected_format": detected_format,
                     "duplicate_gene_entries": duplicates},
            message="Gene/reference overlap is a separate scientific compatibility check.",
        )
        consumed = on_row is not None
        return GeneSetSource(source, options, identities, count, frozenset(genes), detected_format)

    result = validate_once(
        (source,), {"validator": "gene_set_source", "version": 1, **asdict(options)},
        inspect if evidence is None else lambda: evidence, error_type=error_type,
    )
    if on_row is not None and not consumed:
        for row in result.rows(error_type=error_type):
            on_row(row)
    return result


__all__ = ["GeneSetFormat", "GeneSetRow", "GeneSetSource", "validate_gene_set_source"]
