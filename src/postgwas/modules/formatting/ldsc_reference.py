"""HapMap3 rsID-and-allele selection for the LDSC formatter target."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from postgwas.core.io.tables import read_delimited_table

from .table import FormattingError


# These are the strand-unambiguous biallelic SNV pairs accepted by the official
# LDSC munge_sumstats implementation. Keeping the protocol invariant here lets
# reference selection agree with the later LDSC allele merge.
_LDSC_VALID_ALLELE_PAIRS = frozenset({
    "AC", "AG", "CA", "CT", "GA", "GT", "TC", "TG",
})
_DNA_COMPLEMENTS = {"A": "T", "C": "G", "G": "C", "T": "A"}


def resolve_ldsc_merge_alleles_file(config, selected: list[str]) -> Path | None:
    """Validate and canonicalize the optional LDSC reference before extraction."""
    value = config.ldsc_reference.merge_alleles_file
    if "ldsc" not in selected or value is None:
        return None
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise FormattingError(
            "LDSC --merge-alleles file does not exist or is not a file: %s" % path
        )
    if path.stat().st_size <= 0:
        raise FormattingError("LDSC --merge-alleles file is empty: %s" % path)
    return path


def _read_reference(path: Path, config, rsid_pattern: str) -> pl.DataFrame:
    reference_config = config.ldsc_reference
    columns = reference_config.columns
    reference, _delimiter = read_delimited_table(
        path,
        reference_config.delimiter,
        candidates=reference_config.delimiter_candidates,
        minimum_columns=3,
        maximum_columns=reference_config.maximum_columns,
        sample_lines=reference_config.sample_lines,
        null_values=reference_config.null_values,
        infer_schema_length=reference_config.infer_schema_length,
        error_type=FormattingError,
        description="LDSC --merge-alleles reference",
    )
    required = [columns.variant_id, columns.allele_1, columns.allele_2]
    missing = [column for column in required if column not in reference.columns]
    if missing:
        raise FormattingError(
            "LDSC --merge-alleles reference is missing configured columns: %s"
            % ", ".join(missing)
        )
    reference = reference.select(required).rename({
        columns.variant_id: "__ldsc_reference_id",
        columns.allele_1: "__ldsc_reference_a1",
        columns.allele_2: "__ldsc_reference_a2",
    }).with_columns(
        pl.col("__ldsc_reference_id")
        .cast(pl.String).str.strip_chars().str.to_lowercase(),
        pl.col("__ldsc_reference_a1")
        .cast(pl.String).str.strip_chars().str.to_uppercase(),
        pl.col("__ldsc_reference_a2")
        .cast(pl.String).str.strip_chars().str.to_uppercase(),
    )

    invalid_id = (
        pl.col("__ldsc_reference_id").is_null()
        | ~pl.col("__ldsc_reference_id").str.contains(rsid_pattern)
    )
    allele_pair = (
        pl.col("__ldsc_reference_a1") + pl.col("__ldsc_reference_a2")
    )
    invalid_alleles = (
        pl.col("__ldsc_reference_a1").is_null()
        | pl.col("__ldsc_reference_a2").is_null()
        | ~allele_pair.is_in(_LDSC_VALID_ALLELE_PAIRS)
    )
    invalid_ids = int(reference.select(invalid_id.sum()).item() or 0)
    invalid_pairs = int(reference.select(invalid_alleles.sum()).item() or 0)
    if invalid_ids or invalid_pairs:
        raise FormattingError(
            "LDSC --merge-alleles reference contains %s invalid rsIDs and %s "
            "invalid or strand-ambiguous allele pairs. Expected unique rsIDs "
            "and non-palindromic A/C/G/T SNP alleles."
            % (f"{invalid_ids:,}", f"{invalid_pairs:,}")
        )

    duplicate_rows = int(
        reference.select(
            pl.col("__ldsc_reference_id").is_duplicated().sum()
        ).item()
        or 0
    )
    if duplicate_rows:
        duplicate_groups = reference.group_by("__ldsc_reference_id").len().filter(
            pl.col("len") > 1
        ).height
        raise FormattingError(
            "LDSC --merge-alleles reference contains %s duplicated rsIDs across "
            "%s groups; each reference rsID must occur exactly once."
            % (f"{duplicate_rows:,}", f"{duplicate_groups:,}")
        )
    return reference


def _complement(column: str) -> pl.Expr:
    # Unknown alleles remain unchanged and therefore cannot match the validated
    # A/C/G/T LDSC reference. This is equivalent to a null default for the
    # compatibility predicate and supports the repository's Polars baseline.
    return pl.col(column).replace(_DNA_COMPLEMENTS)


def select_ldsc_reference_variants(
    frame: pl.DataFrame,
    config,
    merge_alleles_file: str | Path,
) -> tuple[pl.DataFrame, dict[str, int | bool]]:
    """Retain one HapMap3 allele-compatible record for each selected rsID."""
    canonical = config.canonical_columns
    identifier = canonical.resolved_variant_id
    reference = _read_reference(
        Path(merge_alleles_file),
        config,
        config.variant_identifiers.rsid_pattern,
    )
    original_columns = list(frame.columns)
    occupied = set(original_columns)
    row_column = "__postgwas_ldsc_row"
    while row_column in occupied:
        row_column += "_"
    key_column = "__postgwas_ldsc_rsid"
    while key_column in occupied or key_column == row_column:
        key_column += "_"

    study = frame.with_row_index(row_column).with_columns(
        pl.col(identifier)
        .cast(pl.String).str.strip_chars().str.to_lowercase().alias(key_column)
    )
    reference = reference.rename({"__ldsc_reference_id": key_column})
    joined = study.join(reference, on=key_column, how="inner")

    alternate = canonical.alternate_allele
    reference_allele = canonical.reference_allele
    a1 = "__ldsc_reference_a1"
    a2 = "__ldsc_reference_a2"
    compatible = (
        ((pl.col(alternate) == pl.col(a1)) & (pl.col(reference_allele) == pl.col(a2)))
        | ((pl.col(alternate) == pl.col(a2)) & (pl.col(reference_allele) == pl.col(a1)))
        | ((_complement(alternate) == pl.col(a1)) & (_complement(reference_allele) == pl.col(a2)))
        | ((_complement(alternate) == pl.col(a2)) & (_complement(reference_allele) == pl.col(a1)))
    ).fill_null(False)
    match_column = "__postgwas_ldsc_allele_match"
    joined = joined.with_columns(compatible.alias(match_column))

    input_groups = study.group_by(key_column).len().rename({"len": "input_rows"})
    duplicate_groups = input_groups.filter(pl.col("input_rows") > 1)
    duplicate_rows = int(
        duplicate_groups.select(pl.col("input_rows").sum()).item() or 0
    )
    candidate_groups = joined.group_by(key_column).agg(
        pl.len().alias("reference_candidate_rows"),
        pl.col(match_column).sum().alias("compatible_rows"),
    ).join(input_groups, on=key_column, how="left", coalesce=True)
    ambiguous = candidate_groups.filter(pl.col("compatible_rows") > 1)

    matched = joined.filter(pl.col(match_column))
    if matched.is_empty():
        raise FormattingError(
            "No VCF records match both an rsID and an allele pair in the supplied "
            "LDSC --merge-alleles reference."
        )
    resolved_duplicate_groups = candidate_groups.filter(
        (pl.col("input_rows") > 1) & (pl.col("compatible_rows") == 1)
    ).height
    rows_not_in_reference = frame.height - joined.height
    allele_mismatch_rows = joined.height - matched.height
    output = matched.sort(row_column).select(original_columns)
    return output, {
        "reference_selection_applied": True,
        "reference_variants": reference.height,
        "reference_rows_in": frame.height,
        "reference_rows_out": output.height,
        "rows_excluded_not_in_reference": rows_not_in_reference,
        "rows_excluded_reference_allele_mismatch": allele_mismatch_rows,
        "identifier_duplicate_groups": duplicate_groups.height,
        "identifier_duplicate_rows": duplicate_rows,
        "identifier_duplicate_groups_resolved_by_reference": (
            resolved_duplicate_groups
        ),
        "identifier_duplicate_groups_unresolved_after_reference": ambiguous.height,
    }


__all__ = [
    "resolve_ldsc_merge_alleles_file",
    "select_ldsc_reference_variants",
]
