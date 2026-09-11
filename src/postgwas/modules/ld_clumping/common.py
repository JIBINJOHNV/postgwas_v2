"""Errors, identifier rules, and genomic filters shared by both clumping methods.

Both LD-clumping analyses normalise chromosomes, build canonical variant
identifiers, and exclude the MHC. Keeping one implementation here is what stops
the two methods from disagreeing about what ``chr01`` means.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping

import polars as pl

from postgwas.core.errors import ConfigurationError, PostGWASError


class LDClumpingError(ConfigurationError):
    """LD-clumping configuration, reference, or input contract is invalid.

    Deriving from :class:`ConfigurationError` is deliberate: every condition
    raised as this type is detected before any analysis runs, so module CLIs
    report it as a configuration failure rather than a traceback.
    """


class LDRegionClumpingError(PostGWASError):
    """Annotated-region clumping cannot produce scientifically valid output."""


class PipelineStageError(PostGWASError):
    """Analysis error carrying a stable stage, function, and execution context."""

    def __init__(self, stage, function, message, **context):
        details = " | ".join(
            f"{key}={value}" for key, value in context.items() if value is not None
        )
        text = f"[STAGE: {stage}] [FUNCTION: {function}] {message}"
        if details:
            text = f"{text} | {details}"
        super().__init__(text)
        self.stage = stage
        self.function = function
        self.context = context


def function_error(stage, function, error, **context):
    """Create a contextual error while retaining the original exception."""
    return PipelineStageError(
        stage,
        function,
        f"{type(error).__name__}: {error}",
        **context,
    )


MISSING_ALLELES = frozenset({"", ".", "NA"})


def normalise_chromosome(value) -> str:
    """Return one canonical chromosome label: upper case, no prefix, no zeros.

    ``chr01``, ``CHR1`` and ``1`` all normalise to ``1``; ``chrX`` to ``X``.
    """
    text = str(value).strip().upper()
    if text.startswith("CHR"):
        text = text[3:]
    return str(int(text)) if text.isdigit() else text


def chromosome_expression(column: str) -> pl.Expr:
    """Return the column-wise equivalent of :func:`normalise_chromosome`."""
    text = pl.col(column).cast(pl.Utf8).str.strip_chars().str.to_uppercase()
    text = pl.when(text.str.starts_with("CHR")).then(
        text.str.slice(3)
    ).otherwise(text)
    return (
        pl.when(text.str.contains(r"^\d+$"))
        .then(text.cast(pl.Int64, strict=False).cast(pl.Utf8))
        .otherwise(text)
    )


def chromosome_sort_key(value):
    """Order chromosomes numerically first, then by label."""
    chromosome = normalise_chromosome(value)
    return (0, int(chromosome), "") if chromosome.isdigit() else (1, 0, chromosome)


def canonical_variant_expression(
    chromosome_column: str,
    position_column: str,
    allele_1_column: str,
    allele_2_column: str,
) -> pl.Expr:
    """Return the column-wise equivalent of :func:`canonical_variant_id`."""
    chromosome = chromosome_expression(chromosome_column)
    position = pl.col(position_column).cast(pl.Int64, strict=False)
    allele_1 = pl.col(allele_1_column).cast(pl.Utf8).str.to_uppercase()
    allele_2 = pl.col(allele_2_column).cast(pl.Utf8).str.to_uppercase()
    low = pl.when(allele_1 <= allele_2).then(allele_1).otherwise(allele_2)
    high = pl.when(allele_1 <= allele_2).then(allele_2).otherwise(allele_1)
    invalid = (
        chromosome.is_null()
        | position.is_null()
        | allele_1.is_null()
        | allele_2.is_null()
        | allele_1.is_in(list(MISSING_ALLELES))
        | allele_2.is_in(list(MISSING_ALLELES))
    )
    return pl.when(invalid).then(None).otherwise(
        pl.concat_str(
            [chromosome, position.cast(pl.Utf8), low, high],
            separator="_",
        )
    )


def canonical_variant_id(chrom, pos, allele_1, allele_2):
    """Create one allele-order-independent ID for summary and LD variants."""
    if any(value is None for value in (chrom, pos, allele_1, allele_2)):
        return None
    try:
        pos = int(pos)
    except (TypeError, ValueError):
        return None
    alleles = sorted((str(allele_1).upper(), str(allele_2).upper()))
    if any(allele in MISSING_ALLELES for allele in alleles):
        return None
    return f"{normalise_chromosome(chrom)}_{pos}_{alleles[0]}_{alleles[1]}"


def normalize_ld_id(ld_id):
    """Normalize chr:pos:a1:a2 or chr_pos_a1_a2 LD identifiers."""
    if ld_id is None:
        return None
    parts = str(ld_id).replace(":", "_").split("_")
    if len(parts) != 4:
        return None
    return canonical_variant_id(parts[0], parts[1], parts[2], parts[3])


_REFERENCE_EXCLUSION_DISPLAY = {
    "missing_chromosome_reference": (
        "Chromosome LD reference unavailable",
        "on a chromosome with an unavailable LD reference",
        "on chromosomes with unavailable LD references",
        "chromosome reference unavailable",
    ),
    "missing_from_reference": (
        "Position not found in LD reference",
        "not found at its exact chromosome/position in the LD reference",
        "not found at their exact chromosome/position in the LD reference",
        "not found",
    ),
    "allele_mismatch": (
        "Alleles do not match LD reference",
        "allele mismatch with the LD reference",
        "allele mismatches with the LD reference",
        "allele mismatch",
    ),
    "below_reference_maf": (
        "Below reference MAF threshold",
        "below the reference MAF threshold",
        "below the reference MAF threshold",
        "below MAF",
    ),
}


def reference_exclusion_reason_items(reason_counts: Mapping[str, int]):
    """Return positive exclusion counts in one stable, meaningful order."""
    positive = {
        str(reason): int(count)
        for reason, count in reason_counts.items()
        if int(count) > 0
    }
    ordered_reasons = [
        reason for reason in _REFERENCE_EXCLUSION_DISPLAY if reason in positive
    ]
    ordered_reasons.extend(sorted(set(positive).difference(ordered_reasons)))
    return [(reason, positive[reason]) for reason in ordered_reasons]


def reference_exclusion_reason_counts(
    exclusions: Iterable[Mapping[str, object]],
) -> dict[str, int]:
    """Count the validated reason attached to every excluded index variant."""
    counts: dict[str, int] = {}
    for exclusion in exclusions:
        reason = str(exclusion.get("reason") or "unspecified_reference_exclusion")
        counts[reason] = counts.get(reason, 0) + 1
    return dict(reference_exclusion_reason_items(counts))


def reference_exclusion_reason_label(reason: str) -> str:
    """Return a user-facing label without changing the machine reason code."""
    display = _REFERENCE_EXCLUSION_DISPLAY.get(str(reason))
    return display[0] if display is not None else str(reason).replace("_", " ")


def format_reference_exclusion_reasons(
    reason_counts: Mapping[str, int],
) -> str:
    """Format skipped-index counts so a generic warning total is unnecessary."""
    parts = []
    for reason, count in reference_exclusion_reason_items(reason_counts):
        display = _REFERENCE_EXCLUSION_DISPLAY.get(reason)
        if display is None:
            phrase = reference_exclusion_reason_label(reason)
        else:
            phrase = display[1] if count == 1 else display[2]
        parts.append(f"{count:,} {phrase}")
    return " · ".join(parts) if parts else "none"


def format_compact_reference_exclusion_counts(
    reason_counts: Mapping[str, int],
) -> str:
    """Format a short label-first reason breakdown for live progress rows."""
    parts = []
    for reason, count in reference_exclusion_reason_items(reason_counts):
        display = _REFERENCE_EXCLUSION_DISPLAY.get(reason)
        label = (
            display[3]
            if display is not None
            else reference_exclusion_reason_label(reason)
        )
        parts.append(f"{label} {count:,}")
    return " · ".join(parts) if parts else "none"


def reference_exclusion_analysis_limitation(count: int) -> str:
    """Explain the result limitation caused by excluded GWS variants."""
    count = int(count)
    subject = "variant was" if count == 1 else "variants were"
    return (
        f"{count:,} GWS {subject} excluded from LD clumping; additional "
        "independent signals or loci may be missing"
    )


def exclude_mhc(frame, configuration, *, chromosome_column, position_column):
    """Drop the build-specific MHC interval and report what was removed.

    Returns ``(frame, statistics)`` so every caller audits the exclusion with
    identical fields regardless of which log it writes to.
    """
    rows_in = frame.height
    statistics = {
        "enabled": configuration.remove_mhc,
        "genome_build": configuration.genome_build.value,
        "rows_in": rows_in,
        "rows_out": rows_in,
        "removed": 0,
    }
    if not configuration.remove_mhc:
        return frame, statistics
    mhc = configuration.mhc_regions[configuration.genome_build]
    frame = frame.filter(
        ~(
            (chromosome_expression(chromosome_column) == normalise_chromosome(mhc.chromosome))
            & pl.col(position_column).is_between(mhc.start, mhc.end, closed="both")
        )
    )
    statistics["rows_out"] = frame.height
    statistics["removed"] = rows_in - frame.height
    statistics["region"] = "%s:%d-%d" % (mhc.chromosome, mhc.start, mhc.end)
    return frame, statistics


__all__ = [
    "LDClumpingError",
    "LDRegionClumpingError",
    "PipelineStageError",
    "canonical_variant_expression",
    "canonical_variant_id",
    "chromosome_expression",
    "chromosome_sort_key",
    "exclude_mhc",
    "format_compact_reference_exclusion_counts",
    "format_reference_exclusion_reasons",
    "function_error",
    "normalise_chromosome",
    "normalize_ld_id",
    "reference_exclusion_analysis_limitation",
    "reference_exclusion_reason_counts",
    "reference_exclusion_reason_items",
    "reference_exclusion_reason_label",
]
