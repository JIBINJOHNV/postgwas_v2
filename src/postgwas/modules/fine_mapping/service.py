"""Shared service boundary for fine-mapping engines."""

from __future__ import annotations

from postgwas.config import load_run_configuration_for_module
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.modules.fine_mapping.logging_utils import detailed_file_logging
from postgwas.modules.fine_mapping.output_layout import resolve_output_paths
from postgwas.modules.fine_mapping.overlap_resolution import resolve_overlapping_results
from postgwas.modules.fine_mapping.presentation import FineMappingScreen

_COMMON_OVERRIDES = {
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
    global_overrides = explicit_overrides(
        args,
        {
            "plink": "resources.executables.plink",
            "rscript": "resources.executables.rscript",
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
    args.plink = configuration.resources.executables.plink
    args.rscript = configuration.resources.executables.rscript
    args.finemap_method = module.engine
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


def run_fine_mapping(args):
    """Dispatch a resolved fine-mapping request to its selected engine."""
    args = _resolve_fine_mapping_arguments(args)
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
            ],
        )
        screen.print_final_summary(result, output_paths["pipeline_log_file"])
        return result
    except BaseException:
        screen.fail()
        raise
