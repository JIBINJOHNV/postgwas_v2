"""Named defaults for the fine-mapping package.

All analysis, resource, tolerance, timeout, and output-control values used by
the Python workflow are defined here. Public functions expose the applicable
values as default arguments so callers can inspect or override operational
settings without editing function bodies. The credible-set coverage is
intentionally fixed at 0.95 by scientific requirement.
"""

DEFAULT_CREDIBLE_SET_COVERAGE = 0.95
DEFAULT_COVERAGE_TOLERANCE = 1e-3
DEFAULT_MAX_MAF = 0.5
DEFAULT_MIN_POSITION = 1
DEFAULT_BASES_PER_KILOBASE = 1_000
DEFAULT_SCHEMA_INFERENCE_LENGTH = 10_000

DEFAULT_MHC_CHROMOSOME = "6"
DEFAULT_MHC_START = 25_000_000
DEFAULT_MHC_END = 35_000_000

DEFAULT_TOOL_VERSION_TIMEOUT_SECONDS = 15
DEFAULT_SOFTWARE_VERSION_TIMEOUT_SECONDS = 20
DEFAULT_FINEMAP_RAM_PER_WORKER_GB = 14
DEFAULT_FALLBACK_MEMORY_GB = 16
DEFAULT_PLINK_MEMORY_MB = 2_000
DEFAULT_BGEN_BITS = 8
DEFAULT_EXTERNAL_TOOL_THREADS = 1
DEFAULT_SUSIE_LP_THRESHOLD = 7.3
DEFAULT_SUSIE_MAX_CAUSAL_COMPONENTS = 10
DEFAULT_SUSIE_LD_TIMEOUT_SECONDS = 180
DEFAULT_SUSIE_TIMEOUT_SECONDS = 180
DEFAULT_SUSIE_MIN_RAM_PER_WORKER_GB = 4
DEFAULT_SUSIE_RAM_PER_WORKER_GB = 14

DEFAULT_LD_MAX_CORRELATION = 1.0
DEFAULT_LD_CORRELATION_TOLERANCE = 1e-6
DEFAULT_LD_SYMMETRY_TOLERANCE = 1e-5
DEFAULT_LD_DIAGONAL_TOLERANCE = 1e-4
DEFAULT_LD_EIGENVALUE_TOLERANCE = 1e-4

DEFAULT_PROGRESS_DECIMAL_PLACES = 1
DEFAULT_FINEMAP_PIPELINE_STAGE_TOTAL = 7
DEFAULT_SUSIE_PIPELINE_STAGE_TOTAL = 7
DEFAULT_GENOME_BUILD = "GRCh37"
SUPPORTED_GENOME_BUILDS = ("GRCh37", "GRCh38")


def get_finemap_defaults():
    """Return a documented copy of all Python workflow defaults."""
    return {
        "credible_set_coverage": DEFAULT_CREDIBLE_SET_COVERAGE,
        "coverage_tolerance": DEFAULT_COVERAGE_TOLERANCE,
        "max_maf": DEFAULT_MAX_MAF,
        "minimum_position": DEFAULT_MIN_POSITION,
        "bases_per_kilobase": DEFAULT_BASES_PER_KILOBASE,
        "schema_inference_length": DEFAULT_SCHEMA_INFERENCE_LENGTH,
        "mhc_chromosome": DEFAULT_MHC_CHROMOSOME,
        "mhc_start": DEFAULT_MHC_START,
        "mhc_end": DEFAULT_MHC_END,
        "tool_version_timeout_seconds": DEFAULT_TOOL_VERSION_TIMEOUT_SECONDS,
        "software_version_timeout_seconds": DEFAULT_SOFTWARE_VERSION_TIMEOUT_SECONDS,
        "finemap_ram_per_worker_gb": DEFAULT_FINEMAP_RAM_PER_WORKER_GB,
        "fallback_memory_gb": DEFAULT_FALLBACK_MEMORY_GB,
        "plink_memory_mb": DEFAULT_PLINK_MEMORY_MB,
        "bgen_bits": DEFAULT_BGEN_BITS,
        "external_tool_threads": DEFAULT_EXTERNAL_TOOL_THREADS,
        "susie_lp_threshold": DEFAULT_SUSIE_LP_THRESHOLD,
        "susie_max_causal_components": DEFAULT_SUSIE_MAX_CAUSAL_COMPONENTS,
        "susie_ld_timeout_seconds": DEFAULT_SUSIE_LD_TIMEOUT_SECONDS,
        "susie_timeout_seconds": DEFAULT_SUSIE_TIMEOUT_SECONDS,
        "susie_min_ram_per_worker_gb": DEFAULT_SUSIE_MIN_RAM_PER_WORKER_GB,
        "susie_ram_per_worker_gb": DEFAULT_SUSIE_RAM_PER_WORKER_GB,
        "ld_max_correlation": DEFAULT_LD_MAX_CORRELATION,
        "ld_correlation_tolerance": DEFAULT_LD_CORRELATION_TOLERANCE,
        "ld_symmetry_tolerance": DEFAULT_LD_SYMMETRY_TOLERANCE,
        "ld_diagonal_tolerance": DEFAULT_LD_DIAGONAL_TOLERANCE,
        "ld_eigenvalue_tolerance": DEFAULT_LD_EIGENVALUE_TOLERANCE,
        "progress_decimal_places": DEFAULT_PROGRESS_DECIMAL_PLACES,
        "finemap_pipeline_stage_total": DEFAULT_FINEMAP_PIPELINE_STAGE_TOTAL,
        "susie_pipeline_stage_total": DEFAULT_SUSIE_PIPELINE_STAGE_TOTAL,
        "genome_build": DEFAULT_GENOME_BUILD,
        "supported_genome_builds": list(SUPPORTED_GENOME_BUILDS),
    }
