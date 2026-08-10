"""Fine-mapping protocol invariants.

User-selectable scientific, execution, resource, and output controls belong to
the schema-validated canonical YAML. Values kept here are mathematical units,
domain bounds, supported protocol identifiers, or fixed progress-stage counts.
"""

MAXIMUM_MINOR_ALLELE_FREQUENCY = 0.5
MINIMUM_GENOMIC_POSITION = 1
BASE_PAIRS_PER_KILOBASE = 1_000
MAXIMUM_ABSOLUTE_CORRELATION = 1.0

PROGRESS_DECIMAL_PLACES = 1
FINEMAP_PIPELINE_STAGE_TOTAL = 7
SUSIE_PIPELINE_STAGE_TOTAL = 7
SUPPORTED_GENOME_BUILDS = ("GRCh37", "GRCh38")


def get_finemap_defaults(configuration=None):
    """Return runtime values derived from a validated fine-mapping model."""
    if configuration is None:
        from postgwas.config import load_configuration

        configuration = load_configuration().modules.fine_mapping
    validation = configuration.validation
    runtime = configuration.runtime
    preparation = configuration.summary_statistics_preparation
    finemap = configuration.engines.finemap
    susie = configuration.engines.susie
    return {
        "credible_set_coverage": configuration.credible_set_coverage,
        "coverage_tolerance": validation.credible_set_coverage_tolerance,
        "max_maf": MAXIMUM_MINOR_ALLELE_FREQUENCY,
        "minimum_position": MINIMUM_GENOMIC_POSITION,
        "bases_per_kilobase": BASE_PAIRS_PER_KILOBASE,
        "schema_inference_length": preparation.schema_inference_length,
        "mhc_chromosome": str(configuration.mhc_chromosome),
        "mhc_start": configuration.mhc_start,
        "mhc_end": configuration.mhc_end,
        "tool_version_timeout_seconds": runtime.tool_version_timeout_seconds,
        "software_version_timeout_seconds": runtime.software_version_timeout_seconds,
        "finemap_ram_per_worker_gb": configuration.memory_per_worker_gb,
        "fallback_memory_gb": runtime.fallback_memory_gb,
        "plink_memory_mb": finemap.plink_memory_mb,
        "bgen_bits": finemap.bgen_bits,
        "external_tool_threads": finemap.external_tool_threads,
        "susie_lp_threshold": configuration.locus_lp_threshold,
        "susie_max_causal_components": susie.max_causal_components,
        "susie_ld_timeout_seconds": susie.execution.ld_timeout_seconds,
        "susie_timeout_seconds": susie.execution.susie_timeout_seconds,
        "susie_min_ram_per_worker_gb": configuration.memory_per_worker_gb,
        "susie_ram_per_worker_gb": configuration.memory_per_worker_gb,
        "ld_max_correlation": MAXIMUM_ABSOLUTE_CORRELATION,
        "ld_correlation_tolerance": validation.ld_correlation_tolerance,
        "ld_symmetry_tolerance": validation.ld_symmetry_tolerance,
        "ld_diagonal_tolerance": validation.ld_diagonal_tolerance,
        "ld_eigenvalue_tolerance": validation.ld_eigenvalue_tolerance,
        "progress_decimal_places": PROGRESS_DECIMAL_PLACES,
        "finemap_pipeline_stage_total": FINEMAP_PIPELINE_STAGE_TOTAL,
        "susie_pipeline_stage_total": SUSIE_PIPELINE_STAGE_TOTAL,
        "genome_build": configuration.genome_build.value,
        "supported_genome_builds": list(SUPPORTED_GENOME_BUILDS),
    }


# Compatibility names are derived from the canonical YAML rather than carrying
# independent values. Internal callers should prefer a resolved run model.
_PACKAGED = get_finemap_defaults()
DEFAULT_CREDIBLE_SET_COVERAGE = _PACKAGED["credible_set_coverage"]
DEFAULT_COVERAGE_TOLERANCE = _PACKAGED["coverage_tolerance"]
DEFAULT_MAX_MAF = MAXIMUM_MINOR_ALLELE_FREQUENCY
DEFAULT_MIN_POSITION = MINIMUM_GENOMIC_POSITION
DEFAULT_BASES_PER_KILOBASE = BASE_PAIRS_PER_KILOBASE
DEFAULT_SCHEMA_INFERENCE_LENGTH = _PACKAGED["schema_inference_length"]
DEFAULT_MHC_CHROMOSOME = _PACKAGED["mhc_chromosome"]
DEFAULT_MHC_START = _PACKAGED["mhc_start"]
DEFAULT_MHC_END = _PACKAGED["mhc_end"]
DEFAULT_TOOL_VERSION_TIMEOUT_SECONDS = _PACKAGED["tool_version_timeout_seconds"]
DEFAULT_SOFTWARE_VERSION_TIMEOUT_SECONDS = _PACKAGED[
    "software_version_timeout_seconds"
]
DEFAULT_FINEMAP_RAM_PER_WORKER_GB = _PACKAGED["finemap_ram_per_worker_gb"]
DEFAULT_FALLBACK_MEMORY_GB = _PACKAGED["fallback_memory_gb"]
DEFAULT_PLINK_MEMORY_MB = _PACKAGED["plink_memory_mb"]
DEFAULT_BGEN_BITS = _PACKAGED["bgen_bits"]
DEFAULT_EXTERNAL_TOOL_THREADS = _PACKAGED["external_tool_threads"]
DEFAULT_SUSIE_LP_THRESHOLD = _PACKAGED["susie_lp_threshold"]
DEFAULT_SUSIE_MAX_CAUSAL_COMPONENTS = _PACKAGED["susie_max_causal_components"]
DEFAULT_SUSIE_LD_TIMEOUT_SECONDS = _PACKAGED["susie_ld_timeout_seconds"]
DEFAULT_SUSIE_TIMEOUT_SECONDS = _PACKAGED["susie_timeout_seconds"]
DEFAULT_SUSIE_MIN_RAM_PER_WORKER_GB = _PACKAGED["susie_min_ram_per_worker_gb"]
DEFAULT_SUSIE_RAM_PER_WORKER_GB = _PACKAGED["susie_ram_per_worker_gb"]
DEFAULT_LD_MAX_CORRELATION = MAXIMUM_ABSOLUTE_CORRELATION
DEFAULT_LD_CORRELATION_TOLERANCE = _PACKAGED["ld_correlation_tolerance"]
DEFAULT_LD_SYMMETRY_TOLERANCE = _PACKAGED["ld_symmetry_tolerance"]
DEFAULT_LD_DIAGONAL_TOLERANCE = _PACKAGED["ld_diagonal_tolerance"]
DEFAULT_LD_EIGENVALUE_TOLERANCE = _PACKAGED["ld_eigenvalue_tolerance"]
DEFAULT_PROGRESS_DECIMAL_PLACES = PROGRESS_DECIMAL_PLACES
DEFAULT_FINEMAP_PIPELINE_STAGE_TOTAL = FINEMAP_PIPELINE_STAGE_TOTAL
DEFAULT_SUSIE_PIPELINE_STAGE_TOTAL = SUSIE_PIPELINE_STAGE_TOTAL
DEFAULT_GENOME_BUILD = _PACKAGED["genome_build"]
