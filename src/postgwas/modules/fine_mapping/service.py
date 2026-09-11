"""Shared service boundary for fine-mapping engines."""

from __future__ import annotations

import logging

from postgwas.config import load_module_configuration, load_run_configuration_for_module
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.core.execution.runtime import safe_thread_count
from postgwas.core.io.reports import write_delimited_report
from postgwas.core.preflight import (
    PipelinePreflightEvidence,
    pipeline_preflight_evidence,
    require_pipeline_input_vcf,
)
from postgwas.core.plink import (
    PLINK_BIM_COLUMN_ROLES,
    PLINK_TABLE_DELIMITER_PATTERN,
    validate_plink_bundle_dimensions,
)
from postgwas.core.required_arguments import (
    RequiredArgument,
    require_resolved_arguments,
)
from postgwas.modules.fine_mapping.logging_utils import detailed_file_logging
from postgwas.modules.fine_mapping.output_layout import resolve_output_paths
from postgwas.modules.fine_mapping.overlap_resolution import resolve_overlapping_results
from postgwas.modules.fine_mapping.presentation import FineMappingScreen
from postgwas.modules.fine_mapping.preflight import (
    FineMappingResourcePreflight,
    run_fine_mapping_preflight,
    run_fine_mapping_resource_preflight,
)
from postgwas.modules.fine_mapping.reporting import write_fine_mapping_html_report
from postgwas.modules.formatting.reference_identifiers import (
    BimIdentifierRequirement,
    configure_reference_variant_identifiers,
)


logger = logging.getLogger("postgwas.modules.fine_mapping")

_COMMON_OVERRIDES = {
    "locus_file": "input.locus_file",
    "susie_input_file": "input.susie_summary_statistics_file",
    "finemap_in_files": "input.finemap_summary_statistics_file",
    "finemap_ld_reference": "input.ld_reference_prefix",
    "locus_type": "locus_type",
    "window_kb": "locus_window_kb",
    "lp_threshold": "locus_lp_threshold",
    "finemap_skip_mhc": "skip_mhc",
    "finemap_mhc_chromosome": "mhc_chromosome",
    "finemap_mhc_start": "mhc_start",
    "finemap_mhc_end": "mhc_end",
    "genome_build": "genome_build",
    "minimum_memory_per_worker_gb": "memory_per_worker_gb",
    "maximum_variants_per_locus": (
        "ld_resource_guard.maximum_variants_per_locus"
    ),
}

_SUSIE_OVERRIDES = {
    "L": "engines.susie.max_causal_components",
    "minimum_purity": "engines.susie.minimum_purity",
    "ld_timeout_seconds": "engines.susie.execution.ld_timeout_seconds",
    "susie_timeout_seconds": "engines.susie.execution.susie_timeout_seconds",
    "susie_main_max_iter": "engines.susie.fitting.main_max_iter",
    "susie_recovery_max_iter": "engines.susie.recovery.max_iter",
    "susie_repaired_max_iter": "engines.susie.recovery.repaired_max_iter",
    "susie_reduced_l_max": "engines.susie.recovery.reduced_l_max",
    "susie_ld_eigenvalue_tolerance": (
        "engines.susie.ld_validation.eigenvalue_tolerance"
    ),
    "susie_ld_mismatch_warning_threshold": (
        "engines.susie.ld_validation.mismatch_warning_threshold"
    ),
    "susie_ld_mismatch_failure_threshold": (
        "engines.susie.ld_validation.mismatch_failure_threshold"
    ),
    "susie_ld_repair_maximum_change": (
        "engines.susie.ld_validation.repair_maximum_change"
    ),
    "susie_plink_timeout_seconds": (
        "engines.susie.execution.plink_timeout_seconds"
    ),
    "susie_memory_used_threshold_percent": (
        "engines.susie.memory.used_threshold_percent"
    ),
    "susie_memory_maximum_wait_seconds": (
        "engines.susie.memory.maximum_wait_seconds"
    ),
}

_FINEMAP_OVERRIDES = {
    "ldstore_timeout_seconds": "engines.finemap.ldstore_timeout_seconds",
    "finemap_timeout_seconds": "engines.finemap.finemap_timeout_seconds",
    "finemap_termination_grace_seconds": (
        "engines.finemap.termination_grace_seconds"
    ),
    "n_causal_snps": "engines.finemap.n_causal_snps",
    "n_iter": "engines.finemap.n_iter",
    "n_conv_sss": "engines.finemap.n_conv_sss",
    "prob_conv_sss_tol": "engines.finemap.prob_conv_sss_tol",
    "n_configs_top": "engines.finemap.n_configs_top",
    "corr_config": "engines.finemap.corr_config",
    "pvalue_snps": "engines.finemap.pvalue_snps",
    "cond_pvalue": "engines.finemap.cond_pvalue",
    "prior_std": "engines.finemap.prior_std",
    "prior_k": "engines.finemap.prior_k",
    "force_n_samples": "engines.finemap.force_n_samples",
    "std_effects": "engines.finemap.std_effects",
    "prob_cred_set": "credible_set_coverage",
}


def _resolve_fine_mapping_arguments(args):
    """Resolve YAML defaults and explicit CLI overrides exactly once."""
    overrides = explicit_overrides(
        args,
        {**_COMMON_OVERRIDES, **_SUSIE_OVERRIDES, **_FINEMAP_OVERRIDES},
    )
    if hasattr(args, "finemap_method"):
        overrides["engine"] = args.finemap_method
    if hasattr(args, "sss"):
        overrides["engines.finemap.algorithm"] = "sss"
    elif hasattr(args, "cond"):
        overrides["engines.finemap.algorithm"] = "cond"

    selection_overrides = {
        path: value
        for path, value in overrides.items()
        if path in set(_COMMON_OVERRIDES.values()) | {"engine"}
    }
    plink_was_explicit = hasattr(args, "plink")
    global_overrides = explicit_overrides(
        args,
        {
            "plink": "resources.executables.plink",
            "rscript": "resources.executables.rscript",
            "threads": "execution.threads",
            "memory_gb": "execution.memory_gb",
            "seed": "execution.random_seed",
        },
    )
    selected_configuration = load_run_configuration_for_module(
        "fine_mapping",
        getattr(args, "run_config", None),
        module_overrides=selection_overrides,
        global_overrides=global_overrides,
    )
    selected = selected_configuration.modules.fine_mapping
    engine_specific = (
        _FINEMAP_OVERRIDES if selected.engine == "susie" else _SUSIE_OVERRIDES
    )
    inapplicable = sorted(
        name for name in engine_specific if hasattr(args, name)
    )
    if selected.engine == "susie" and (hasattr(args, "sss") or hasattr(args, "cond")):
        inapplicable.append("algorithm")
    if inapplicable:
        raise ValueError(
            "%s controls do not apply when engine is %s: %s"
            % (
                "FINEMAP" if selected.engine == "susie" else "SuSiE",
                selected.engine,
                ", ".join(sorted(set(inapplicable))),
            )
        )
    configuration = load_run_configuration_for_module(
        "fine_mapping",
        getattr(args, "run_config", None),
        module_overrides=overrides,
        global_overrides=global_overrides,
    )
    module = configuration.modules.fine_mapping
    args.threads = configuration.execution.threads
    args.memory_gb = configuration.execution.memory_gb
    args.seed = configuration.execution.random_seed
    args.plink = (
        configuration.resources.executables.plink
        if module.engine == "susie" or plink_was_explicit
        else configuration.resources.executables.plink2
    )
    args.rscript = configuration.resources.executables.rscript
    args.bgenix = configuration.resources.executables.bgenix
    args.ldstore = configuration.resources.executables.ldstore
    args.finemap_executable = configuration.resources.executables.finemap
    args.finemap_method = module.engine
    args.locus_file = module.input.locus_file
    args.susie_input_file = module.input.susie_summary_statistics_file
    args.finemap_in_files = module.input.finemap_summary_statistics_file
    args.finemap_ld_reference = module.input.ld_reference_prefix
    args.locus_type = module.locus_type
    args.window_kb = module.locus_window_kb
    args.lp_threshold = module.locus_lp_threshold
    args.finemap_skip_mhc = module.skip_mhc
    args.finemap_include_mhc = not module.skip_mhc
    args.finemap_mhc_chromosome = module.mhc_chromosome
    args.finemap_mhc_start = module.mhc_start
    args.finemap_mhc_end = module.mhc_end
    args.genome_build = module.genome_build.value
    args.sample_size_policy = module.sample_size.policy
    args.sample_size_summary_statistic = module.sample_size.summary_statistic
    args.sample_size_relative_range_warning_threshold = (
        module.sample_size.relative_range_warning_threshold
    )
    args.maximum_variants_per_locus = (
        module.ld_resource_guard.maximum_variants_per_locus
    )
    args.fine_mapping_validation = module.validation.model_dump()
    args.fine_mapping_runtime = module.runtime.model_dump()
    args.fine_mapping_html_report = module.html_report.model_dump()
    args.finemap_ld_peak_matrix_multiplier = (
        module.ld_resource_guard.finemap_peak_matrix_multiplier
    )
    args.susie_ld_peak_matrix_multiplier = (
        module.ld_resource_guard.susie_peak_matrix_multiplier
    )
    args.overlap_resolution = module.overlap_resolution.model_dump()
    args.fine_mapping_output_layout = module.output_layout.model_dump()
    args.fine_mapping_logging = configuration.logging.model_dump()
    args.summary_statistics_preparation = (
        module.summary_statistics_preparation.model_dump()
    )
    args.resolved_fine_mapping_configuration = module.model_dump(mode="json")
    from postgwas.modules.fine_mapping.defaults import get_finemap_defaults

    args.fine_mapping_runtime_defaults = get_finemap_defaults(module)
    if module.engine == "finemap":
        settings = module.engines.finemap
        if settings.prior_k:
            raise ValueError(
                "FINEMAP prior_k requires a K file, but the current FINEMAP "
                "master-file generator does not create a K-file column. Disable "
                "prior_k until K-file input is implemented."
            )
        for name, value in settings.model_dump().items():
            if name != "algorithm":
                setattr(args, name, value)
        args.minimum_memory_per_worker_gb = module.memory_per_worker_gb
        args.prob_cred_set = module.credible_set_coverage
        args.finemap_termination_grace_seconds = (
            settings.termination_grace_seconds
        )
        args.algorithm = settings.algorithm
        args.sss = settings.algorithm == "sss"
        args.cond = settings.algorithm == "cond"
    else:
        settings = module.engines.susie
        args.L = settings.max_causal_components
        args.minimum_purity = settings.minimum_purity
        args.ld_timeout_seconds = settings.execution.ld_timeout_seconds
        args.susie_timeout_seconds = settings.execution.susie_timeout_seconds
        args.recovery_audit_filename = settings.recovery_audit_filename
        args.minimum_memory_per_worker_gb = module.memory_per_worker_gb
    return args


def _require_fine_mapping_inputs(args) -> None:
    """Validate all engine-specific required values after YAML/CLI resolution."""
    summary_requirement = (
        RequiredArgument(
            "--susie-input-file",
            "modules.fine_mapping.input.susie_summary_statistics_file",
            args.susie_input_file,
        )
        if args.finemap_method == "susie"
        else RequiredArgument(
            "--finemap-in-files",
            "modules.fine_mapping.input.finemap_summary_statistics_file",
            args.finemap_in_files,
        )
    )
    require_resolved_arguments(
        (
            RequiredArgument(
                "--locus-file",
                "modules.fine_mapping.input.locus_file",
                args.locus_file,
            ),
            summary_requirement,
            RequiredArgument(
                "--finemap-ld-reference",
                "modules.fine_mapping.input.ld_reference_prefix",
                args.finemap_ld_reference,
            ),
        )
    )


def preflight_fine_mapping_pipeline(
    args,
    *,
    preflight_evidence=None,
) -> PipelinePreflightEvidence:
    """Validate fine-mapping references and tools before generated inputs exist."""
    require_pipeline_input_vcf(preflight_evidence)
    args = _resolve_fine_mapping_arguments(args)
    require_resolved_arguments((
        RequiredArgument(
            "--finemap-ld-reference",
            "modules.fine_mapping.input.ld_reference_prefix",
            args.finemap_ld_reference,
        ),
    ))
    safe_thread_count(
        args.threads,
        args.minimum_memory_per_worker_gb,
        available_ram_gb=args.memory_gb,
        enforce_memory_budget=True,
        reporter=None,
    )
    resources = run_fine_mapping_resource_preflight(args)
    # Both engines extract reference variants by the exact SNP identifier;
    # the column named rsid in FINEMAP is not an rsID-only scientific contract.
    reference_files = dict(resources.reference_files)
    configure_reference_variant_identifiers(
        args,
        load_module_configuration("formatting", getattr(args, "run_config", None)),
        [BimIdentifierRequirement(
            consumer="Fine-mapping (%s)" % args.finemap_method,
            formatter_target=args.finemap_method,
            bim_file=reference_files[".bim"],
            column_roles=PLINK_BIM_COLUMN_ROLES,
            delimiter_pattern=PLINK_TABLE_DELIMITER_PATTERN,
        )],
    )
    validate_plink_bundle_dimensions(
        reference_files,
        variants=args.variant_id_observations[args.finemap_method]["variants"],
    )
    return pipeline_preflight_evidence(
        "finemap",
        preflight_evidence,
        resources=resources,
        deferred_checks=(
            "Validate formatter summary statistics and LD-clumping loci.",
            "Validate locus-level variant and coordinate overlap with the PLINK panel.",
        ),
    )


def run_fine_mapping(
    args,
    *,
    resource_preflight: FineMappingResourcePreflight | None = None,
    upstream_results=None,
):
    """Dispatch a resolved fine-mapping request to its selected engine."""
    args = _resolve_fine_mapping_arguments(args)
    _require_fine_mapping_inputs(args)
    safe_thread_count(
        args.threads,
        args.minimum_memory_per_worker_gb,
        available_ram_gb=args.memory_gb,
        enforce_memory_budget=True,
        reporter=None,
    )
    method = args.finemap_method
    if method == "susie":
        from postgwas.modules.fine_mapping.engines.susie.adapter import (
            run_parallel_susie,
        )

        run_engine = run_parallel_susie
        engine_label = "SuSiE-RSS"
    elif method == "finemap":
        from postgwas.modules.fine_mapping.engines.finemap.adapter import (
            run_finemap_pipeline,
        )

        run_engine = run_finemap_pipeline
        engine_label = "FINEMAP"
    else:
        raise ValueError(
            "Unknown fine-mapping engine %r. Expected 'susie' or 'finemap'."
            % method
        )

    screen = FineMappingScreen(
        engine_label,
        show_progress=args.fine_mapping_logging["show_progress"],
        label_width=args.fine_mapping_logging["terminal_label_width"],
    )
    screen.start()
    try:
        initial_output_paths = resolve_output_paths(
            args.output_directory,
            args.fine_mapping_output_layout,
            args.dataset_id,
        )
        with detailed_file_logging(
            "postgwas.modules.fine_mapping",
            initial_output_paths["pipeline_log_file"],
            args.fine_mapping_logging["file_level"],
        ):
            logger.info(
                "Fine-mapping preflight started: engine=%s dataset=%s",
                method,
                args.dataset_id,
            )
            try:
                preflight = run_fine_mapping_preflight(
                    args,
                    resource_preflight=resource_preflight,
                )
            except BaseException as exc:
                write_delimited_report(
                    (
                        {
                            "category": "preflight",
                            "check": "input_and_resource_validation",
                            "status": "failed",
                            "value": "",
                            "path": "",
                            "detail": str(exc) or type(exc).__name__,
                        },
                    ),
                    initial_output_paths["preflight_validation_file"],
                    fieldnames=(
                        "category",
                        "check",
                        "status",
                        "value",
                        "path",
                        "detail",
                    ),
                    delimiter="\t",
                    null_value="",
                )
                logger.exception("Fine-mapping preflight failed: %s", exc)
                raise
            preflight_report = write_delimited_report(
                preflight.audit_rows(),
                initial_output_paths["preflight_validation_file"],
                fieldnames=(
                    "category",
                    "check",
                    "status",
                    "value",
                    "path",
                    "detail",
                ),
                delimiter="\t",
                null_value="",
            )
            logger.info(
                "Fine-mapping preflight completed: summary_variants=%d "
                "eligible_loci=%d reference_variants=%d reference_samples=%d "
                "report=%s",
                preflight.inputs.summary_statistics_rows,
                preflight.inputs.eligible_loci,
                preflight.reference.variants,
                preflight.reference.samples,
                preflight_report,
            )
            if (
                preflight.reference.loci_with_reference_variants
                < preflight.reference.loci_requested
            ):
                logger.warning(
                    "Fine-mapping preflight found eligible loci without a "
                    "coordinate-concordant PLINK reference variant: covered=%d "
                    "eligible=%d; locus processing will retain locus-scoped "
                    "failure handling",
                    preflight.reference.loci_with_reference_variants,
                    preflight.reference.loci_requested,
                )
            for tool in preflight.tools:
                logger.info(
                    "Fine-mapping runtime validated: tool=%s path=%s version=%s",
                    tool.name,
                    tool.path,
                    tool.version,
                )
        args._fine_mapping_preflight = preflight
        screen.complete(1, preflight.input_screen_fields())
        screen.complete(2, preflight.resource_screen_fields())
        primary_result = run_engine(args, screen=screen)
        output_paths = resolve_output_paths(
            primary_result["output_dir"],
            args.fine_mapping_output_layout,
            args.dataset_id,
        )
        with detailed_file_logging(
            "postgwas.modules.fine_mapping",
            output_paths["pipeline_log_file"],
            args.fine_mapping_logging["file_level"],
            mode="a",
        ):
            result = resolve_overlapping_results(
                args,
                primary_result,
                run_engine,
            )
            result["preflight_validation"] = str(preflight_report)
            html_report = write_fine_mapping_html_report(
                args,
                preflight,
                result,
                output_paths,
                upstream_results=upstream_results,
            )
            result["html_report"] = str(html_report)
            logger.info(
                "Fine-mapping scientific HTML report completed: %s",
                html_report,
            )
        screen.complete(
            8,
            [
                (
                    "analysis",
                    "Connected overlap groups",
                    result.get("n_overlap_groups", 0),
                ),
                (
                    "success",
                    "Successful joint reruns",
                    result.get("n_joint_rerun_successful", 0),
                ),
                (
                    "success",
                    "Authoritative credible sets",
                    result.get(
                        "n_final_credible_sets",
                        result.get("n_credible_sets", 0),
                    ),
                ),
                (
                    "success" if result["status"] == "success" else "warning",
                    "Overall status",
                    result["status"].replace("_", " "),
                ),
                ("success", "Scientific HTML report", html_report.name),
            ],
        )
        screen.print_final_summary(result, output_paths["pipeline_log_file"])
        return result
    except BaseException:
        screen.fail()
        raise
