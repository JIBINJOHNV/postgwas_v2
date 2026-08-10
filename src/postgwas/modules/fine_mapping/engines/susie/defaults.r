# SuSiE run controls are loaded from the schema-validated configuration emitted
# by Python. Only protocol invariants required while sourcing function
# definitions are defined in R.

SUSIE_PROTOCOL_INVARIANTS <- list(
  minimum_position = 1L,
  supported_genome_builds = c("GRCh37", "GRCh38"),
  workers = 1L,
  progress_digits = 1L,
  pipeline_stage_total = 5L
)

require_susie_config_value <- function(value, path, length_required = 1L) {
  invalid_length <- !is.null(length_required) && length(value) != length_required
  invalid_missing <- !is.list(value) && anyNA(value)
  if (is.null(value) || invalid_length || invalid_missing) {
    stop(
      "Resolved SuSiE configuration is missing or invalid at '", path, "'",
      call. = FALSE
    )
  }
  value
}

load_susie_defaults <- function(configuration_file) {
  if (length(configuration_file) != 1L || !file.exists(configuration_file)) {
    stop("Resolved SuSiE configuration file does not exist", call. = FALSE)
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The R package 'jsonlite' is required to load SuSiE configuration")
  }
  configuration <- jsonlite::fromJSON(
    configuration_file, simplifyVector = TRUE
  )
  engine <- require_susie_config_value(
    configuration$engines$susie, "engines.susie", NULL
  )
  fitting <- require_susie_config_value(
    engine$fitting, "engines.susie.fitting", NULL
  )
  recovery <- require_susie_config_value(
    engine$recovery, "engines.susie.recovery", NULL
  )
  ld_validation <- require_susie_config_value(
    engine$ld_validation, "engines.susie.ld_validation", NULL
  )
  execution <- require_susie_config_value(
    engine$execution, "engines.susie.execution", NULL
  )
  memory <- require_susie_config_value(
    engine$memory, "engines.susie.memory", NULL
  )
  plotting <- require_susie_config_value(
    engine$plotting, "engines.susie.plotting", NULL
  )
  sample_size <- require_susie_config_value(
    configuration$sample_size, "sample_size", NULL
  )
  resource_guard <- require_susie_config_value(
    configuration$ld_resource_guard, "ld_resource_guard", NULL
  )

  resolved <- list(
    credible_set_coverage = require_susie_config_value(
      configuration$credible_set_coverage, "credible_set_coverage"
    ),
    credible_set_tolerance = require_susie_config_value(
      engine$credible_set_tolerance, "engines.susie.credible_set_tolerance"
    ),
    min_abs_corr = require_susie_config_value(
      engine$minimum_purity, "engines.susie.minimum_purity"
    ),
    ld_eigenvalue_tolerance = require_susie_config_value(
      ld_validation$eigenvalue_tolerance,
      "engines.susie.ld_validation.eigenvalue_tolerance"
    ),
    ld_eigenvalue_floor = require_susie_config_value(
      ld_validation$eigenvalue_floor,
      "engines.susie.ld_validation.eigenvalue_floor"
    ),
    ld_lambda_warning = require_susie_config_value(
      ld_validation$mismatch_warning_threshold,
      "engines.susie.ld_validation.mismatch_warning_threshold"
    ),
    ld_lambda_failure = require_susie_config_value(
      ld_validation$mismatch_failure_threshold,
      "engines.susie.ld_validation.mismatch_failure_threshold"
    ),
    ld_repair_max_change = require_susie_config_value(
      ld_validation$repair_maximum_change,
      "engines.susie.ld_validation.repair_maximum_change"
    ),
    lp_threshold = require_susie_config_value(
      configuration$locus_lp_threshold, "locus_lp_threshold"
    ),
    max_causal_components = require_susie_config_value(
      engine$max_causal_components, "engines.susie.max_causal_components"
    ),
    main_max_iter = require_susie_config_value(
      fitting$main_max_iter, "engines.susie.fitting.main_max_iter"
    ),
    recovery_max_iter = require_susie_config_value(
      recovery$max_iter, "engines.susie.recovery.max_iter"
    ),
    repaired_max_iter = require_susie_config_value(
      recovery$repaired_max_iter, "engines.susie.recovery.repaired_max_iter"
    ),
    reduced_l_max = require_susie_config_value(
      recovery$reduced_l_max, "engines.susie.recovery.reduced_l_max"
    ),
    ld_retry_multiplier = require_susie_config_value(
      recovery$ld_retry_multiplier,
      "engines.susie.recovery.ld_retry_multiplier"
    ),
    ld_retry_add_seconds = require_susie_config_value(
      recovery$ld_retry_add_seconds,
      "engines.susie.recovery.ld_retry_add_seconds"
    ),
    recovery_timeout_multiplier = require_susie_config_value(
      recovery$timeout_multiplier,
      "engines.susie.recovery.timeout_multiplier"
    ),
    ld_repair_timeout_multiplier = require_susie_config_value(
      recovery$ld_repair_timeout_multiplier,
      "engines.susie.recovery.ld_repair_timeout_multiplier"
    ),
    repaired_timeout_multiplier = require_susie_config_value(
      recovery$repaired_timeout_multiplier,
      "engines.susie.recovery.repaired_timeout_multiplier"
    ),
    ld_timeout_seconds = require_susie_config_value(
      execution$ld_timeout_seconds,
      "engines.susie.execution.ld_timeout_seconds"
    ),
    susie_timeout_seconds = require_susie_config_value(
      execution$susie_timeout_seconds,
      "engines.susie.execution.susie_timeout_seconds"
    ),
    plink_timeout_seconds = require_susie_config_value(
      execution$plink_timeout_seconds,
      "engines.susie.execution.plink_timeout_seconds"
    ),
    process_poll_seconds = require_susie_config_value(
      execution$process_poll_seconds,
      "engines.susie.execution.process_poll_seconds"
    ),
    process_terminate_grace_seconds = require_susie_config_value(
      execution$termination_grace_seconds,
      "engines.susie.execution.termination_grace_seconds"
    ),
    stderr_tail_lines = require_susie_config_value(
      execution$stderr_tail_lines, "engines.susie.execution.stderr_tail_lines"
    ),
    verbose = require_susie_config_value(
      execution$verbose, "engines.susie.execution.verbose"
    ),
    memory_used_threshold_pct = require_susie_config_value(
      memory$used_threshold_percent,
      "engines.susie.memory.used_threshold_percent"
    ),
    memory_check_interval_seconds = require_susie_config_value(
      memory$check_interval_seconds,
      "engines.susie.memory.check_interval_seconds"
    ),
    memory_max_wait_seconds = require_susie_config_value(
      memory$maximum_wait_seconds,
      "engines.susie.memory.maximum_wait_seconds"
    ),
    min_ram_per_worker_gb = require_susie_config_value(
      configuration$memory_per_worker_gb, "memory_per_worker_gb"
    ),
    meminfo_read_lines = require_susie_config_value(
      memory$meminfo_read_lines, "engines.susie.memory.meminfo_read_lines"
    ),
    skip_mhc = require_susie_config_value(configuration$skip_mhc, "skip_mhc"),
    mhc_chromosome = as.character(require_susie_config_value(
      configuration$mhc_chromosome, "mhc_chromosome"
    )),
    mhc_start = require_susie_config_value(configuration$mhc_start, "mhc_start"),
    mhc_end = require_susie_config_value(configuration$mhc_end, "mhc_end"),
    genome_build = require_susie_config_value(
      configuration$genome_build, "genome_build"
    ),
    recovery_audit_filename = require_susie_config_value(
      engine$recovery_audit_filename, "engines.susie.recovery_audit_filename"
    ),
    sample_size_policy = require_susie_config_value(
      sample_size$policy, "sample_size.policy"
    ),
    sample_size_summary_statistic = require_susie_config_value(
      sample_size$summary_statistic, "sample_size.summary_statistic"
    ),
    sample_size_relative_range_warning_threshold = require_susie_config_value(
      sample_size$relative_range_warning_threshold,
      "sample_size.relative_range_warning_threshold"
    ),
    maximum_variants_per_locus = require_susie_config_value(
      resource_guard$maximum_variants_per_locus,
      "ld_resource_guard.maximum_variants_per_locus"
    ),
    ld_peak_matrix_multiplier = require_susie_config_value(
      resource_guard$susie_peak_matrix_multiplier,
      "ld_resource_guard.susie_peak_matrix_multiplier"
    ),
    reserved_worker_memory_gb = require_susie_config_value(
      configuration$memory_per_worker_gb, "memory_per_worker_gb"
    ),
    plot_min_pip_label = plotting$minimum_pip_to_label,
    plot_background_size = plotting$background_point_size,
    plot_background_alpha = plotting$background_alpha,
    plot_cs_size = plotting$credible_set_point_size,
    plot_cs_alpha = plotting$credible_set_alpha,
    plot_max_ld_dimension = plotting$maximum_ld_dimension,
    plot_width = plotting$width,
    plot_height = plotting$height,
    plot_dpi = plotting$dpi,
    plot_ld_palette_colors = plotting$ld_palette_colors,
    plot_legend_gradient_points = plotting$legend_gradient_points,
    plot_legend_title_size = plotting$legend_title_size,
    plot_legend_text_size = plotting$legend_text_size,
    plot_legend_key_size_cm = plotting$legend_key_size_cm,
    plot_legend_scale = plotting$legend_scale,
    plot_credible_set_legend_size = plotting$credible_set_legend_size,
    plot_labels_per_credible_set = plotting$labels_per_credible_set,
    plot_large_locus_threshold = plotting$large_locus_threshold,
    plot_medium_locus_threshold = plotting$medium_locus_threshold,
    plot_large_locus_scale = plotting$large_locus_scale,
    plot_medium_locus_scale = plotting$medium_locus_scale,
    plot_label_lift = plotting$label_lift,
    plot_label_size = plotting$label_size,
    plot_label_segment_alpha = plotting$label_segment_alpha,
    plot_label_segment_width = plotting$label_segment_width,
    plot_panel_heights = require_susie_config_value(
      plotting$panel_heights, "engines.susie.plotting.panel_heights", 3L
    ),
    plot_panel_widths = require_susie_config_value(
      plotting$panel_widths, "engines.susie.plotting.panel_widths", 2L
    )
  )
  utils::modifyList(SUSIE_PROTOCOL_INVARIANTS, resolved)
}

# Populated by run_susie_parallel_cli.R before any analysis function is called.
SUSIE_DEFAULTS <- SUSIE_PROTOCOL_INVARIANTS
