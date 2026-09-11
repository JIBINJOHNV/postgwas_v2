"""Shared resolution of build-specific MHC and chromosome analysis scope."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any


def resolve_genomic_analysis_scope(
    *,
    genome_build: object,
    genomes: Mapping[str, Any],
    mhc,
    chromosomes,
    normalize_chromosome: Callable[[object], str],
    analysis_name: str,
    error_type: type[Exception],
) -> dict:
    """Resolve one validated scope without inferring a genome build."""
    build = getattr(genome_build, "value", genome_build)
    build_name = str(build)
    genome = genomes.get(build_name)
    if genome is None:
        raise error_type(
            "No genome resources are configured for %s" % build_name
        )
    region = mhc.region_override or genome.regions.get("mhc")
    if (mhc.excludes_snps or mhc.excludes_genes) and region is None:
        raise error_type(
            "%s MHC exclusion requires resources.genomes.%s.regions.mhc or "
            "the complete --mhc-chrom/--mhc-start/--mhc-end override."
            % (analysis_name, build_name)
        )
    if region is not None and region.start < 1:
        raise error_type("%s requires a one-based MHC start position" % analysis_name)
    normalized_exclusions = [
        normalize_chromosome(value) for value in chromosomes.exclude
    ]
    if len(normalized_exclusions) != len(set(normalized_exclusions)):
        raise error_type(
            "%s excluded chromosome labels are duplicated after normalization"
            % analysis_name
        )
    return {
        "mhc_policy": mhc.policy,
        "exclude_mhc_snps": mhc.excludes_snps,
        "exclude_mhc_genes": mhc.excludes_genes,
        "mhc_region": (
            {
                "chromosome": normalize_chromosome(region.chromosome),
                "start": region.start,
                "end": region.end,
            }
            if region is not None else None
        ),
        "mhc_source": "override" if mhc.region_override else "genome_resources",
        "exclude_chromosomes": normalized_exclusions,
    }


__all__ = ["resolve_genomic_analysis_scope"]
