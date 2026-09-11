"""LD clumping: annotated regions, FUMA-style r² clumping, and GCTA-COJO.

``region`` selects the strongest association in each annotated population LD
block. ``standard`` performs FUMA-style two-stage r² clumping against a
manifest-validated, bidirectionally indexed population LD reference and defines
genomic risk loci.
``cojo-slct`` reuses GCTA stepwise selection and groups the unchanged selected
signals into configured physical-distance loci.
"""

from postgwas.modules.ld_clumping.common import (
    LDClumpingError,
    LDRegionClumpingError,
    PipelineStageError,
)

__all__ = [
    "LDClumpingError",
    "LDRegionClumpingError",
    "PipelineStageError",
    "ld_clump_by_regions",
    "ld_clump_standard",
    "preflight_ld_clumping",
    "run_ld_clump_direct",
    "validate_ld_clumping_configuration",
]


def __getattr__(name):
    """Expose analysis entry points without importing polars at package import."""
    if name == "ld_clump_by_regions":
        from postgwas.modules.ld_clumping.ld_prune_region import ld_clump_by_regions

        return ld_clump_by_regions
    if name == "ld_clump_standard":
        from postgwas.modules.ld_clumping.ld_prune_standard import ld_clump_standard

        return ld_clump_standard
    if name in (
        "preflight_ld_clumping",
        "run_ld_clump_direct",
        "validate_ld_clumping_configuration",
    ):
        from postgwas.modules.ld_clumping import service

        return getattr(service, name)
    raise AttributeError("module %r has no attribute %r" % (__name__, name))
