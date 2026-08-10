"""Shared dense-LD resource accounting for fine-mapping engines."""

from __future__ import annotations


VARIANT_LIMIT_FAILURE_REASON = "maximum_variants_per_locus_exceeded"
MEMORY_WARNING_REASON = "estimated_ld_peak_memory_exceeds_reserved_worker_memory"

# Both engines materialize IEEE-754 double-precision dense LD matrices.
_FLOAT64_BYTES = 8
_BYTES_PER_GIBIBYTE = 1024**3


def dense_ld_resource_qc(
    variant_count: int,
    maximum_variants: int,
    peak_matrix_multiplier: float,
    reserved_worker_memory_gb: float,
) -> dict:
    """Return auditable limit and approximate dense-LD memory information."""
    variant_count = int(variant_count)
    maximum_variants = int(maximum_variants)
    peak_matrix_multiplier = float(peak_matrix_multiplier)
    reserved_worker_memory_gb = float(reserved_worker_memory_gb)
    if variant_count < 0:
        raise ValueError("variant_count cannot be negative")
    if maximum_variants < 1:
        raise ValueError("maximum_variants_per_locus must be at least 1")
    if peak_matrix_multiplier < 1:
        raise ValueError("LD peak matrix multiplier must be at least 1")
    if reserved_worker_memory_gb <= 0:
        raise ValueError("reserved worker memory must be greater than zero")
    matrix_gb = (
        variant_count * variant_count * _FLOAT64_BYTES / _BYTES_PER_GIBIBYTE
    )
    estimated_peak_gb = matrix_gb * peak_matrix_multiplier
    exceeds_limit = variant_count > maximum_variants
    return {
        "input_variant_count": variant_count,
        "maximum_variants_per_locus": maximum_variants,
        "dense_ld_matrix_gb": matrix_gb,
        "ld_peak_matrix_multiplier": peak_matrix_multiplier,
        "estimated_peak_ld_memory_gb": estimated_peak_gb,
        "reserved_worker_memory_gb": reserved_worker_memory_gb,
        "resource_warning_reason": (
            MEMORY_WARNING_REASON
            if not exceeds_limit and estimated_peak_gb > reserved_worker_memory_gb
            else None
        ),
        "resource_failure_reason": (
            VARIANT_LIMIT_FAILURE_REASON if exceeds_limit else None
        ),
    }


def variant_limit_failure_detail(qc: dict) -> str:
    """Build one actionable YAML/CLI override message for an oversized locus."""
    return (
        f"Locus contains {qc['input_variant_count']} harmonized variants, exceeding "
        f"maximum_variants_per_locus={qc['maximum_variants_per_locus']}. "
        "Override with --maximum-variants-per-locus COUNT or edit "
        "modules.fine_mapping.ld_resource_guard.maximum_variants_per_locus in "
        "the run YAML. Increasing the limit can require substantially more RAM."
    )


def combine_reasons(*reasons: str | None) -> str | None:
    """Combine non-empty reason codes once while preserving their order."""
    unique = []
    for reason in reasons:
        if reason and reason not in unique:
            unique.append(reason)
    return ";".join(unique) if unique else None
