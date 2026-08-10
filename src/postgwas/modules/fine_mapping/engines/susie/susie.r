#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(glue)
  library(susieR)
  library(Matrix)
  library(R.utils)
  library(ggplot2)
  library(ggrepel)
  library(parallel)
})

# =============================================================================
# 0. Locate and source utilities from same folder
# =============================================================================
get_script_dir <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  fileArg <- grep("^--file=", cmdArgs, value = TRUE)
  if (length(fileArg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", fileArg))))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  getwd()
}

this_dir <- get_script_dir()
DEFAULTS_R <- glue("{this_dir}/defaults.r")
if (!file.exists(DEFAULTS_R)) {
  stop(glue("Defaults file not found: {DEFAULTS_R}"))
}
source(DEFAULTS_R, local = FALSE)
UTILITIES_R <- glue("{this_dir}/utlities.r")  # keep your filename as-is
if (!file.exists(UTILITIES_R)) {
  stop(glue("Utilities file not found: {UTILITIES_R}"))
}
source(UTILITIES_R, local = FALSE)

# =============================================================================
# Safe operator
# =============================================================================
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

# =============================================================================
# Cross-platform per-locus logger (NO open-connection leak)
# =============================================================================
make_locus_logger <- function(logfile) {
  dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
  function(...) {
    txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""), "\n")
    cat(txt, file = logfile, append = TRUE)
    cat(txt)  # also to console
  }
}

initialize_susie_progress <- function(progress_file) {
  data.table::fwrite(
    data.table::data.table(
      timestamp_utc = character(), scope = character(), stage = character(),
      status = character(), completed = integer(), total = integer(),
      percentage = numeric(), remaining = integer(), message = character()
    ),
    progress_file,
    sep = "\t"
  )
  invisible(progress_file)
}

record_susie_progress <- function(
  progress_file, scope, stage, completed, total, status,
  message_text = "", digits = SUSIE_DEFAULTS$progress_digits
) {
  completed <- as.integer(completed)
  total <- as.integer(total)
  if (!is.finite(completed) || !is.finite(total) ||
      completed < 0L || total < 0L || completed > total) {
    stop("Invalid SuSiE progress values")
  }
  percentage <- if (total == 0L) 100 else round(100 * completed / total, digits)
  remaining <- total - completed
  row <- data.table::data.table(
    timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    scope = as.character(scope), stage = as.character(stage),
    status = as.character(status), completed = completed, total = total,
    percentage = percentage, remaining = remaining,
    message = as.character(message_text)
  )
  data.table::fwrite(row, progress_file, sep = "\t", append = TRUE, col.names = FALSE)
  percentage_text <- formatC(percentage, format = "f", digits = digits)
  message(sprintf(
    "[PROGRESS] scope=%s stage=%s status=%s completed=%d/%d percentage=%s%% remaining=%d message=%s",
    scope, stage, status, completed, total, percentage_text, remaining, message_text
  ))
  invisible(row)
}

# =============================================================================
# 3. Memory throttle (Linux only) with max-wait safeguard
# =============================================================================
wait_for_memory <- function(threshold_used_pct = SUSIE_DEFAULTS$memory_used_threshold_pct,
                            check_interval_sec = SUSIE_DEFAULTS$memory_check_interval_seconds,
                            verbose = TRUE,
                            max_wait_sec = SUSIE_DEFAULTS$memory_max_wait_seconds) {
  if (Sys.info()[["sysname"]] %in% c("Darwin", "Windows")) return(invisible(TRUE))

  t0 <- Sys.time()

  repeat {
    mem <- tryCatch(
      readLines("/proc/meminfo", n = SUSIE_DEFAULTS$meminfo_read_lines),
      error = function(e) character(0)
    )
    if (length(mem) == 0) {
      Sys.sleep(check_interval_sec)
      next
    }

    get_kb <- function(k) {
      line <- mem[grepl(paste0("^", k, ":"), mem)]
      if (length(line) == 0) return(0)
      val <- suppressWarnings(as.numeric(sub(".* ([0-9]+) kB.*", "\\1", line[1])))
      val %||% 0
    }

    total <- get_kb("MemTotal")
    # Better available approximation
    available <- get_kb("MemAvailable")
    if (!is.finite(available) || available <= 0) {
      available <- get_kb("MemFree") + get_kb("Buffers") + get_kb("Cached") + get_kb("SReclaimable")
    }

    if (!is.finite(total) || total <= 0) {
      Sys.sleep(check_interval_sec)
      next
    }

    used_pct <- 100 * (1 - available / total)

    if (used_pct < threshold_used_pct) break

    if (verbose) message("High memory: ", round(used_pct, 1), "% used — waiting ", check_interval_sec, "s...")

    if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > max_wait_sec) {
      if (verbose) message("wait_for_memory(): max_wait_sec reached → proceeding anyway")
      break
    }

    Sys.sleep(check_interval_sec)
  }

  invisible(TRUE)
}

# =============================================================================
# 4. Timeout wrapper (FIXED: expr evaluated inside withTimeout)
# =============================================================================
run_with_timeout <- function(expr, timeout, step, tag, log_msg, verbose = TRUE) {
  expr_sub <- substitute(expr)

  if (!is.finite(timeout) || is.na(timeout) || timeout <= 0) {
    return(tryCatch(eval(expr_sub, envir = parent.frame()),
                    error = function(e) {
                      log_msg("Error during ", step, " for ", tag, ": ", e$message)
                      NULL
                    }))
  }

  tryCatch(
    R.utils::withTimeout(eval(expr_sub, envir = parent.frame()),
                         timeout = timeout,
                         onTimeout = "error"),
    TimeoutException = function(e) {
      log_msg("Timeout (>", timeout, "s) during ", step, " for ", tag)
      structure(list(timeout = TRUE), class = "timeout")
    },
    error = function(e) {
      log_msg("Error during ", step, " for ", tag, ": ", e$message)
      NULL
    }
  )
}

# =============================================================================
# 5. Fast LD repair (same logic)
# =============================================================================
LD_fix_fast <- function(
  ld.mat, tag = "", log_msg,
  max_allowed_change = SUSIE_DEFAULTS$ld_repair_max_change,
  eigenvalue_floor = SUSIE_DEFAULTS$ld_eigenvalue_floor
) {
  ld.mat <- (ld.mat + t(ld.mat)) / 2
  original <- ld.mat
  start_t <- Sys.time()
  log_msg("Fast LD repair started for ", tag)

  e <- eigen(ld.mat, symmetric = TRUE)
  e$values[e$values < eigenvalue_floor] <- eigenvalue_floor
  repaired <- e$vectors %*% diag(e$values) %*% t(e$vectors)
  repaired <- (repaired + t(repaired)) / 2
  d <- sqrt(diag(repaired))
  if (any(!is.finite(d)) || any(d <= 0)) stop("LD repair produced an invalid diagonal")
  repaired <- repaired / tcrossprod(d)
  diag(repaired) <- 1
  max_change <- max(abs(repaired - original), na.rm = TRUE)
  log_msg("LD repair maximum absolute change = ", signif(max_change, 5))
  if (!is.finite(max_change) || max_change > max_allowed_change) {
    stop("LD repair changed the correlation matrix excessively")
  }

  ok <- isTRUE(tryCatch({ chol(repaired); TRUE }, error = function(e) FALSE))
  if (ok) {
    log_msg("Fast LD repair OK (", round(difftime(Sys.time(), start_t, units = "secs"), 2), "s)")
    return(repaired)
  }

  log_msg("Fast repair failed → falling back to nearPD")
  npd <- Matrix::nearPD(ld.mat, corr = TRUE, keepDiag = TRUE)
  repaired <- as.matrix(npd$mat)
  repaired <- (repaired + t(repaired)) / 2
  diag(repaired) <- 1
  max_change <- max(abs(repaired - original), na.rm = TRUE)
  log_msg("nearPD maximum absolute change = ", signif(max_change, 5))
  if (!is.finite(max_change) || max_change > max_allowed_change) {
    stop("nearPD changed the correlation matrix excessively")
  }
  repaired
}

# =============================================================================
# Helpers
# =============================================================================
make_safe_tag <- function(tag) gsub("[^A-Za-z0-9_.-]", "_", tag)

make_susie_audit_row <- function(
  locus, sequence, stage, action, status, reason,
  ld_status = NA_character_, repairable = NA,
  L = NA_integer_, max_iter = NA_integer_, timeout_seconds = NA_real_,
  ld_z_mismatch_lambda = NA_real_, ld_min_eigenvalue = NA_real_,
  details = NA_character_
) {
  data.table::data.table(
    timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    locus_chr = as.character(locus$chr),
    locus_start = as.integer(locus$start),
    locus_end = as.integer(locus$end),
    genomic_locus = as.character(locus$genomic_locus),
    sequence = as.integer(sequence),
    stage = as.character(stage),
    action = as.character(action),
    status = as.character(status),
    reason = as.character(reason),
    ld_status = as.character(ld_status),
    repairable = as.logical(repairable),
    L = as.integer(L),
    max_iter = as.integer(max_iter),
    timeout_seconds = as.numeric(timeout_seconds),
    ld_z_mismatch_lambda = as.numeric(ld_z_mismatch_lambda),
    ld_min_eigenvalue = as.numeric(ld_min_eigenvalue),
    details = as.character(details)
  )
}


summarize_locus_nef <- function(values, policy, summary_statistic, warning_threshold) {
  values <- as.numeric(values)
  if (summary_statistic != "median" || policy != "warn") {
    stop("Unsupported sample-size policy or summary statistic")
  }
  if (!is.finite(warning_threshold) || warning_threshold < 0) {
    stop("Sample-size relative-range warning threshold must be non-negative")
  }
  if (!length(values) || any(!is.finite(values)) || any(values <= 0)) {
    stop("Effective sample sizes must be finite and greater than zero")
  }
  nef_median <- stats::median(values)
  nef_min <- min(values)
  nef_max <- max(values)
  relative_range <- (nef_max - nef_min) / nef_median
  warning_reason <- if (relative_range > warning_threshold) {
    "nef_relative_range_exceeds_threshold"
  } else {
    NA_character_
  }
  list(
    nef_min = nef_min, nef_max = nef_max, nef_median = nef_median,
    nef_relative_range = relative_range, selected_nef = nef_median,
    nef_policy = policy, nef_summary_statistic = summary_statistic,
    nef_warning_threshold = warning_threshold,
    warning_reason = warning_reason
  )
}


annotate_sample_size_qc <- function(table, sample_size_qc) {
  table <- data.table::copy(data.table::as.data.table(table))
  for (name in names(sample_size_qc)) {
    data.table::set(table, j = name, value = sample_size_qc[[name]])
  }
  table
}


combine_warning_reasons <- function(...) {
  reasons <- unique(as.character(unlist(list(...), use.names = FALSE)))
  reasons <- reasons[!is.na(reasons) & nzchar(reasons)]
  if (length(reasons)) paste(reasons, collapse = ";") else NA_character_
}


summarize_locus_ld_resources <- function(
  variant_count, maximum_variants, peak_matrix_multiplier,
  reserved_worker_memory_gb
) {
  variant_count <- as.integer(variant_count)
  maximum_variants <- as.integer(maximum_variants)
  peak_matrix_multiplier <- as.numeric(peak_matrix_multiplier)
  reserved_worker_memory_gb <- as.numeric(reserved_worker_memory_gb)
  if (!is.finite(variant_count) || variant_count < 0L) {
    stop("variant_count must be a non-negative integer")
  }
  if (!is.finite(maximum_variants) || maximum_variants < 1L) {
    stop("maximum_variants_per_locus must be at least 1")
  }
  if (!is.finite(peak_matrix_multiplier) || peak_matrix_multiplier < 1) {
    stop("LD peak matrix multiplier must be at least 1")
  }
  if (!is.finite(reserved_worker_memory_gb) || reserved_worker_memory_gb <= 0) {
    stop("reserved worker memory must be greater than zero")
  }
  # Protocol invariant: PLINK/SuSiE LD is materialized as float64 p-by-p data.
  dense_ld_matrix_gb <- variant_count^2 * 8 / 1024^3
  estimated_peak_ld_memory_gb <- dense_ld_matrix_gb * peak_matrix_multiplier
  exceeds_limit <- variant_count > maximum_variants
  list(
    input_variant_count = variant_count,
    maximum_variants_per_locus = maximum_variants,
    dense_ld_matrix_gb = dense_ld_matrix_gb,
    ld_peak_matrix_multiplier = peak_matrix_multiplier,
    estimated_peak_ld_memory_gb = estimated_peak_ld_memory_gb,
    reserved_worker_memory_gb = reserved_worker_memory_gb,
    resource_warning_reason = if (!exceeds_limit &&
      estimated_peak_ld_memory_gb > reserved_worker_memory_gb) {
      "estimated_ld_peak_memory_exceeds_reserved_worker_memory"
    } else {
      NA_character_
    },
    resource_failure_reason = if (exceeds_limit) {
      "maximum_variants_per_locus_exceeded"
    } else {
      NA_character_
    }
  )
}


variant_limit_failure_detail <- function(resource_qc) {
  paste0(
    "Locus contains ", resource_qc$input_variant_count,
    " harmonized variants, exceeding maximum_variants_per_locus=",
    resource_qc$maximum_variants_per_locus,
    ". Override with --maximum-variants-per-locus COUNT or edit ",
    "modules.fine_mapping.ld_resource_guard.maximum_variants_per_locus ",
    "in the run YAML. Increasing the limit can require substantially more RAM."
  )
}


annotate_resource_qc <- function(table, resource_qc) {
  table <- data.table::copy(data.table::as.data.table(table))
  for (name in names(resource_qc)) {
    data.table::set(table, j = name, value = resource_qc[[name]])
  }
  table
}







# =============================================================================
# 6. Main locus processor (PASS 1)
# =============================================================================
process_locus <- function(
  locus, df, sample_id, ld_ref, plink, analysis_folder,
  lp_threshold, L, timeout_ld, timeout_susie, skip_mhc,
  mhc_chr, mhc_start, mhc_end, verbose = TRUE,
  genome_build = SUSIE_DEFAULTS$genome_build,
  defaults = SUSIE_DEFAULTS
) {
  analysis_folder <- normalizePath(analysis_folder, mustWork = TRUE)

  tag <- glue("{locus$chr}_{locus$start}_{locus$end}")
  safe_tag <- make_safe_tag(tag)

  logfile <- file.path(analysis_folder, "logs",
                       glue("{sample_id}_{safe_tag}.log"))
  msg <- make_locus_logger(logfile)

  if (skip_mhc && locus$chr == mhc_chr &&
      locus$start < mhc_end && locus$end > mhc_start) {
    msg("Skipping MHC region")
    return(list(status = "SKIP", reason = "mhc", locus = locus))
  }

  if (!nrow(df))
    return(list(status = "SKIP", reason = "no_snps", locus = locus))
  outside_boundary <- df$CHR != as.character(locus$chr) |
    df$BP < locus$start | df$BP > locus$end
  if (any(outside_boundary)) {
    stop("Exact locus summary-statistics file contains variants outside its boundary")
  }
  selected <- data.table::copy(data.table::as.data.table(df))

  if ("LP" %in% names(selected)) {
    maxlp <- suppressWarnings(max(selected$LP, na.rm = TRUE))
    if (is.finite(maxlp) && maxlp < lp_threshold)
      return(list(status = "SKIP", reason = "low_signal", locus = locus))
  }

  resource_qc <- summarize_locus_ld_resources(
    nrow(selected), defaults$maximum_variants_per_locus,
    defaults$ld_peak_matrix_multiplier,
    defaults$reserved_worker_memory_gb
  )
  if (!is.na(resource_qc$resource_failure_reason)) {
    failure_detail <- variant_limit_failure_detail(resource_qc)
    msg("RESOURCE GUARD FAILED [", resource_qc$resource_failure_reason,
        "]: ", failure_detail)
    return(list(
      status = "FAILED",
      reason = resource_qc$resource_failure_reason,
      failure_detail = failure_detail,
      locus = locus,
      resource_qc = resource_qc,
      audit = list(make_susie_audit_row(
        locus, 1L, "resource_guard", "validate_variant_count", "failed",
        resource_qc$resource_failure_reason,
        ld_status = "not_constructed", repairable = FALSE,
        details = failure_detail
      ))
    ))
  }
  if (!is.na(resource_qc$resource_warning_reason)) {
    msg("WARNING [", resource_qc$resource_warning_reason,
        "]: estimated peak dense-LD memory ",
        signif(resource_qc$estimated_peak_ld_memory_gb, 5),
        " GiB exceeds reserved worker memory ",
        signif(resource_qc$reserved_worker_memory_gb, 5), " GiB")
  }

  wait_for_memory(verbose = verbose)

  # ---------------- LD MATRIX ----------------
  msg("LD extraction START")
  ld <- tryCatch(
    postgwas_ld_matrix(
      variants = selected$SNP,
      bfile = ld_ref,
      plink_bin = plink,
      tag = gsub(":", "-", tag),
      with_alleles = FALSE,
      logfile = glue("{analysis_folder}/ld_matrix_related/{sample_id}_{safe_tag}.log"),
      output_folder = analysis_folder,
      timeout_seconds = timeout_ld
    ),
    error = function(e) {
      msg("❌ LD extraction ERROR (attempt 1): ", conditionMessage(e))
      NULL
    }
  )

  # ---- Retry with longer timeout ----
  if (is.null(ld)) {
    retry_timeout <- max(
      timeout_ld * defaults$ld_retry_multiplier,
      timeout_ld + defaults$ld_retry_add_seconds
    )
    msg("🔁 LD retry with longer timeout (", retry_timeout, "s)")

    ld <- tryCatch(
      postgwas_ld_matrix(
        variants = selected$SNP,
        bfile = ld_ref,
        plink_bin = plink,
        tag = gsub(":", "-", tag),
        with_alleles = FALSE,
        logfile = glue("{analysis_folder}/ld_matrix_related/{sample_id}_{safe_tag}.log"),
        output_folder = analysis_folder,
        timeout_seconds = retry_timeout
      ),
      error = function(e) {
        msg("❌ LD extraction ERROR (attempt 2): ", conditionMessage(e))
        NULL
      }
    )
  }


  if (is.null(ld)) {
    msg("LD extraction FAILED after retry")
    return(list(
      status = "FAILED", reason = "ld_extraction_failed", locus = locus,
      resource_qc = resource_qc
    ))
  }

  aligned <- tryCatch(
    align_sumstats_to_ld(selected, ld, msg),
    error = function(e) {
      msg("❌ LD allele alignment failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(aligned)) {
    return(list(
      status = "FAILED", reason = "allele_alignment_failed", locus = locus,
      resource_qc = resource_qc
    ))
  }
  selected <- aligned$selected
  ld <- aligned$ld
  data.table::set(selected, j = "variable_index", value = seq_len(nrow(selected)))

  z <- selected$EZ
  sample_size_qc <- tryCatch(
    summarize_locus_nef(
      selected$NEF, defaults$sample_size_policy,
      defaults$sample_size_summary_statistic,
      defaults$sample_size_relative_range_warning_threshold
    ),
    error = function(e) NULL
  )
  if (is.null(sample_size_qc)) {
    return(list(
      status = "FAILED", reason = "invalid_sample_size", locus = locus,
      resource_qc = resource_qc
    ))
  }
  n_eff <- sample_size_qc$selected_nef
  sample_size_qc$warning_reason <- combine_warning_reasons(
    sample_size_qc$warning_reason, resource_qc$resource_warning_reason
  )
  if (!is.na(sample_size_qc$warning_reason)) {
    msg("WARNING [", sample_size_qc$warning_reason, "]: NEF min=",
        signif(sample_size_qc$nef_min, 7), ", median=",
        signif(sample_size_qc$nef_median, 7), ", max=",
        signif(sample_size_qc$nef_max, 7), ", relative range=",
        signif(sample_size_qc$nef_relative_range, 5))
  }

  ld_validation <- tryCatch(
    list(
      ok = TRUE,
      value = validate_ld_for_susie(
        ld, z, n_eff, msg,
        eigenvalue_tolerance = defaults$ld_eigenvalue_tolerance,
        lambda_warning = defaults$ld_lambda_warning,
        lambda_failure = defaults$ld_lambda_failure
      )
    ),
    susie_ld_validation_error = function(e) list(ok = FALSE, error = e),
    error = function(e) list(
      ok = FALSE,
      error = structure(
        list(
          message = conditionMessage(e), reason = "ld_validation_internal_error",
          repairable = FALSE, lambda = NA_real_, min_eigenvalue = NA_real_
        ),
        class = c("susie_ld_validation_error", "error", "condition")
      )
    )
  )
  if (!isTRUE(ld_validation$ok)) {
    error <- ld_validation$error
    msg("❌ LD validation failed [", error$reason, "]: ", error$message)
    audit <- list(make_susie_audit_row(
      locus, 1L, "ld_validation", "validate_original_ld", "failed",
      error$reason,
      ld_status = "invalid",
      repairable = error$repairable,
      ld_z_mismatch_lambda = error$lambda,
      ld_min_eigenvalue = error$min_eigenvalue,
      details = error$message
    ))
    if (!isTRUE(error$repairable)) {
      return(list(
        status = "FAILED", reason = error$reason, locus = locus,
        ld_qc_lambda = error$lambda,
        ld_min_eigenvalue = error$min_eigenvalue,
        sample_size_qc = sample_size_qc,
        resource_qc = resource_qc,
        audit = audit
      ))
    }
    return(list(
      status = "RECOVERY", reason = error$reason, recovery_route = "ld_repair",
      locus = locus, selected = selected, ld.mat = ld, n_eff = n_eff,
      sample_size_qc = sample_size_qc,
      resource_qc = resource_qc,
      ld_qc_lambda = error$lambda, ld_min_eigenvalue = error$min_eigenvalue,
      audit = audit
    ))
  }
  ld_qc <- ld_validation$value
  ld <- ld_qc$ld
  sample_size_qc$warning_reason <- combine_warning_reasons(
    sample_size_qc$warning_reason, ld_qc$warning_reason
  )

  fit_input <- NULL
  fit_input_transferred <- FALSE
  on.exit({
    if (!fit_input_transferred) cleanup_susie_fit_input(fit_input, msg)
  }, add = TRUE)
  fit_input_error <- NULL
  fit_input <- tryCatch(
    create_susie_fit_input(
      z, ld, n_eff, "valid_original", safe_tag, msg
    ),
    error = function(e) {
      fit_input_error <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(fit_input)) {
    msg("❌ SuSiE fit input serialization failed: ", fit_input_error)
    return(list(
      status = "FAILED", reason = "susie_fit_input_serialization_failed",
      failure_detail = fit_input_error, locus = locus,
      ld_qc_lambda = ld_qc$lambda,
      ld_min_eigenvalue = ld_qc$min_eigenvalue,
      sample_size_qc = sample_size_qc,
      resource_qc = resource_qc,
      fit_input_qc = summarize_susie_fit_input(NULL),
      audit = list(make_susie_audit_row(
        locus, 1L, "fit_input", "serialize_valid_ld", "failed",
        "susie_fit_input_serialization_failed",
        ld_status = "valid", repairable = FALSE,
        ld_z_mismatch_lambda = ld_qc$lambda,
        ld_min_eigenvalue = ld_qc$min_eigenvalue,
        details = fit_input_error
      ))
    ))
  }

  # ---------------- SuSiE (HARD TIMEOUT) ----------------
  fitted <- run_susie_hard_timeout(
    fit_input = fit_input,
    L = L,
    max_iter = defaults$main_max_iter,
    timeout_s = timeout_susie,
    tag = safe_tag,
    log_msg = msg,
    susie_verbose = TRUE,
    coverage = defaults$credible_set_coverage,
    min_abs_corr = defaults$min_abs_corr
  )

  if (!is_susie_fit_success(fitted)) {
    fit_reason <- if (inherits(fitted, "susie_fit_failure")) {
      fitted$reason
    } else {
      "susie_invalid_result"
    }
    fit_message <- if (inherits(fitted, "susie_fit_failure")) {
      fitted$message
    } else {
      "SuSiE returned an invalid result"
    }
    msg("SuSiE FAILED [", fit_reason, "] → fitting recovery")
    recovery_result <- list(
      status = "RECOVERY", reason = fit_reason, recovery_route = "fit_retry",
      locus = locus, selected = selected, ld.mat = ld, n_eff = n_eff,
      fit_input = fit_input,
      sample_size_qc = sample_size_qc,
      resource_qc = resource_qc,
      ld_qc_lambda = ld_qc$lambda,
      ld_min_eigenvalue = ld_qc$min_eigenvalue,
      audit = list(make_susie_audit_row(
        locus, 1L, "primary_fit", "fit_valid_ld", "failed", fit_reason,
        ld_status = "valid", repairable = FALSE,
        L = L, max_iter = defaults$main_max_iter,
        timeout_seconds = timeout_susie,
        ld_z_mismatch_lambda = ld_qc$lambda,
        ld_min_eigenvalue = ld_qc$min_eigenvalue,
        details = paste(
          fit_message, susie_fit_input_audit_detail(fit_input), sep = "; "
        )
      ))
    )
    fit_input_transferred <- TRUE
    return(recovery_result)
  }

  # ---------------- SUCCESS PATH ----------------
  selected_sus <- annotate_susie(
    fitted, selected, locus, min_abs_corr = defaults$min_abs_corr
  )
  selected_sus <- annotate_sample_size_qc(selected_sus, sample_size_qc)
  selected_sus <- annotate_resource_qc(selected_sus, resource_qc)
  selected_sus[, `:=`(
    analysis_status = "success",
    outcome_reason = "primary_fit_converged"
  )]
  cred_df <- selected_sus[!is.na(cs) & cs != -1]

  plot_susie_loci_ld(
    df = selected_sus,
    ld = ld,
    sample_id = sample_id,
    outdir = glue("{analysis_folder}/plots/"),
    min_pip_label = defaults$plot_min_pip_label,
    bg_size = defaults$plot_background_size,
    bg_alpha = defaults$plot_background_alpha,
    cs_size = defaults$plot_cs_size,
    cs_alpha = defaults$plot_cs_alpha
  )

  generate_flames_files(
    fitted = fitted,
    ld = ld,
    snp_df = selected,
    outfile = glue("{analysis_folder}/flames_input/{sample_id}_{safe_tag}"),
    genomic_locus = locus$genomic_locus,
    genome_build = genome_build,
    target_coverage = defaults$credible_set_coverage,
    resource_qc = resource_qc
  )

  dir.create(glue("{analysis_folder}/locus_files"), recursive = TRUE, showWarnings = FALSE)
  dir.create(glue("{analysis_folder}/rds_files"),   recursive = TRUE, showWarnings = FALSE)

  cred_file <- glue("{analysis_folder}/locus_files/{sample_id}_{safe_tag}_SUSIE_credible_sets.cred1")
  data.table::fwrite(cred_df, cred_file, sep = " ", quote = FALSE)
  saveRDS(fitted, glue("{analysis_folder}/rds_files/{sample_id}_{safe_tag}_SUSIE.rds"))

  msg("Locus COMPLETED successfully")

  list(
    status = "OK",
    vars = selected_sus,
    cs = cred_df,
    fitted = fitted,
    locus = locus,
    cred_file = cred_file,
    ld_qc_lambda = ld_qc$lambda,
    ld_min_eigenvalue = ld_qc$min_eigenvalue,
    sample_size_qc = sample_size_qc,
    resource_qc = resource_qc,
    fit_input_qc = summarize_susie_fit_input(fit_input),
    audit = list(make_susie_audit_row(
      locus, 1L, "primary_fit", "fit_valid_ld", "success", "primary_fit_converged",
      ld_status = "valid", repairable = FALSE,
      L = L, max_iter = defaults$main_max_iter,
      timeout_seconds = timeout_susie,
      ld_z_mismatch_lambda = ld_qc$lambda,
      ld_min_eigenvalue = ld_qc$min_eigenvalue,
      details = susie_fit_input_audit_detail(fit_input)
    ))
  )
}

# =============================================================================
# 7. Recovery (PASS 2)
# =============================================================================
recover_locus <- function(
  job, sample_id, analysis_folder,
  L, timeout_ld, timeout_susie,
  verbose = TRUE,
  genome_build = SUSIE_DEFAULTS$genome_build,
  defaults = SUSIE_DEFAULTS
) {
  locus <- job$locus
  selected <- job$selected
  ld.mat <- job$ld.mat
  n_eff <- job$n_eff
  ld_qc_lambda <- job$ld_qc_lambda
  ld_min_eigenvalue <- job$ld_min_eigenvalue
  sample_size_qc <- job$sample_size_qc
  resource_qc <- job$resource_qc
  fit_input <- job$fit_input %||% NULL
  route <- job$recovery_route
  audit <- job$audit %||% list()

  tag <- glue("{locus$chr}_{locus$start}_{locus$end}")
  safe_tag <- make_safe_tag(tag)
  logfile <- glue("{analysis_folder}/logs/{sample_id}_{safe_tag}_RECOVERY.log")
  msg <- make_locus_logger(logfile)
  on.exit(cleanup_susie_fit_input(fit_input, msg), add = TRUE)

  add_audit <- function(...) {
    audit[[length(audit) + 1L]] <<- make_susie_audit_row(
      locus, length(audit) + 1L, ...
    )
  }
  failure_reason <- function(fit) {
    if (inherits(fit, "susie_fit_failure")) fit$reason else "susie_invalid_result"
  }
  failure_message <- function(fit) {
    if (inherits(fit, "susie_fit_failure")) fit$message else "SuSiE returned an invalid result"
  }
  failed_result <- function(reason) {
    list(
      status = "FAILED", locus = locus, reason = reason,
      ld_qc_lambda = ld_qc_lambda,
      ld_min_eigenvalue = ld_min_eigenvalue,
      sample_size_qc = sample_size_qc,
      resource_qc = resource_qc,
      fit_input_qc = summarize_susie_fit_input(fit_input),
      audit = audit
    )
  }
  sus_to_out <- function(fitted, ld_matrix, reason, suffix = "") {
    selected_sus <- annotate_susie(
      fitted = fitted, selected = selected, locus = locus,
      min_abs_corr = defaults$min_abs_corr
    )
    selected_sus <- annotate_sample_size_qc(selected_sus, sample_size_qc)
    selected_sus <- annotate_resource_qc(selected_sus, resource_qc)
    selected_sus[, `:=`(
      analysis_status = "recovered",
      outcome_reason = as.character(reason)
    )]
    plot_susie_loci_ld(
      df = selected_sus, ld = ld_matrix, sample_id = sample_id,
      outdir = glue("{analysis_folder}/plots/"),
      min_pip_label = defaults$plot_min_pip_label,
      bg_size = defaults$plot_background_size,
      bg_alpha = defaults$plot_background_alpha,
      cs_size = defaults$plot_cs_size,
      cs_alpha = defaults$plot_cs_alpha
    )
    generate_flames_files(
      fitted = fitted, ld = ld_matrix, snp_df = selected,
      outfile = glue("{analysis_folder}/flames_input/{sample_id}_{safe_tag}"),
      genomic_locus = locus$genomic_locus, genome_build = genome_build,
      target_coverage = defaults$credible_set_coverage,
      resource_qc = resource_qc
    )
    cred_df <- selected_sus[!is.na(cs) & cs != -1]
    dir.create(glue("{analysis_folder}/locus_files"), recursive = TRUE, showWarnings = FALSE)
    dir.create(glue("{analysis_folder}/rds_files"), recursive = TRUE, showWarnings = FALSE)
    cred_file <- glue("{analysis_folder}/locus_files/{sample_id}_{safe_tag}_SUSIE_credible_sets{suffix}.cred1")
    data.table::fwrite(cred_df, cred_file, sep = " ", quote = FALSE)
    saveRDS(fitted, glue("{analysis_folder}/rds_files/{sample_id}_{safe_tag}_SUSIE{suffix}.rds"))
    list(
      status = "RECOVERED", reason = reason, vars = selected_sus, cs = cred_df,
      fitted = fitted, locus = locus, cred_file = cred_file,
      ld_qc_lambda = ld_qc_lambda,
      ld_min_eigenvalue = ld_min_eigenvalue,
      sample_size_qc = sample_size_qc,
      resource_qc = resource_qc,
      fit_input_qc = summarize_susie_fit_input(fit_input),
      audit = audit
    )
  }
  run_fit_attempt <- function(attempt_name, attempt_l, max_iter, timeout, ld_status) {
    fit <- run_susie_hard_timeout(
      fit_input = fit_input, L = attempt_l,
      max_iter = max_iter, timeout_s = timeout,
      tag = paste0(safe_tag, "_", attempt_name), log_msg = msg,
      susie_verbose = TRUE, coverage = defaults$credible_set_coverage,
      min_abs_corr = defaults$min_abs_corr
    )
    succeeded <- is_susie_fit_success(fit)
    add_audit(
      stage = "recovery_fit", action = attempt_name,
      status = if (succeeded) "success" else "failed",
      reason = if (succeeded) "susie_converged" else failure_reason(fit),
      ld_status = ld_status, repairable = FALSE,
      L = attempt_l, max_iter = max_iter, timeout_seconds = timeout,
      ld_z_mismatch_lambda = ld_qc_lambda,
      ld_min_eigenvalue = ld_min_eigenvalue,
      details = paste(
        if (succeeded) "SuSiE converged" else failure_message(fit),
        susie_fit_input_audit_detail(fit_input), sep = "; "
      )
    )
    fit
  }

  if (is.null(ld.mat)) {
    add_audit(
      "recovery_validation", "validate_recovery_input", "failed", "ld_mat_null",
      ld_status = "missing", repairable = FALSE, details = "Recovery LD matrix is NULL"
    )
    return(failed_result("ld_mat_null"))
  }
  if (!is.finite(n_eff) || n_eff <= 0) {
    add_audit(
      "recovery_validation", "validate_recovery_input", "failed", "invalid_sample_size",
      ld_status = "unknown", repairable = FALSE, details = "Recovery sample size is invalid"
    )
    return(failed_result("invalid_sample_size"))
  }

  if (identical(route, "fit_retry")) {
    fit_input_valid <- tryCatch(
      {
        validate_susie_fit_input(
          fit_input,
          expected_matrix_state = "valid_original",
          expected_variant_count = nrow(selected),
          expected_n_eff = n_eff
        )
        TRUE
      },
      error = function(e) {
        add_audit(
          "recovery_validation", "validate_fit_input", "failed",
          "susie_fit_input_unavailable", ld_status = "valid",
          repairable = FALSE, details = conditionMessage(e)
        )
        FALSE
      }
    )
    if (!fit_input_valid) return(failed_result("susie_fit_input_unavailable"))
    retry_timeout <- timeout_susie * defaults$recovery_timeout_multiplier
    msg("RECOVERY: valid LD; retrying only the SuSiE fit")
    fit <- run_fit_attempt(
      "higher_iterations", L, defaults$recovery_max_iter,
      retry_timeout, "valid"
    )
    if (is_susie_fit_success(fit)) {
      msg("RECOVERY SUCCESS: higher iterations on validated LD")
      return(sus_to_out(fit, ld.mat, "higher_iterations"))
    }

    reduced_l <- min(defaults$reduced_l_max, as.integer(L))
    fit <- run_fit_attempt(
      "reduced_L", reduced_l, defaults$recovery_max_iter,
      retry_timeout, "valid"
    )
    if (is_susie_fit_success(fit)) {
      msg("RECOVERY SUCCESS: reduced L on validated LD")
      return(sus_to_out(fit, ld.mat, glue("reduced_L{reduced_l}")))
    }
    msg("FITTING RECOVERY FAILED; validated LD was not altered")
    return(failed_result("fit_recovery_exhausted"))
  }

  if (!identical(route, "ld_repair")) {
    add_audit(
      "recovery_validation", "route_recovery", "failed", "invalid_recovery_route",
      ld_status = "unknown", repairable = FALSE, details = as.character(route)
    )
    return(failed_result("invalid_recovery_route"))
  }
  if (!is.null(fit_input)) {
    add_audit(
      "recovery_validation", "validate_fit_input_state", "failed",
      "unexpected_fit_input_for_invalid_ld",
      ld_status = "invalid", repairable = FALSE,
      details = "LD repair recovery received a fit payload before LD repair"
    )
    return(failed_result("unexpected_fit_input_for_invalid_ld"))
  }

  msg("RECOVERY: repairing non-PSD LD before any further SuSiE fit")
  wait_for_memory(verbose = verbose)
  repair_error <- NULL
  ld_repaired <- tryCatch(
    R.utils::withTimeout(
      LD_fix_fast(ld.mat, safe_tag, msg),
      timeout = timeout_ld * defaults$ld_repair_timeout_multiplier,
      onTimeout = "error"
    ),
    error = function(e) {
      repair_error <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(ld_repaired)) {
    add_audit(
      "ld_repair", "repair_non_psd_ld", "failed", "ld_repair_failed",
      ld_status = "invalid", repairable = TRUE,
      timeout_seconds = timeout_ld * defaults$ld_repair_timeout_multiplier,
      ld_z_mismatch_lambda = ld_qc_lambda,
      ld_min_eigenvalue = ld_min_eigenvalue,
      details = repair_error
    )
    msg("LD repair FAILED: ", repair_error)
    return(failed_result("ld_repair_failed"))
  }
  add_audit(
    "ld_repair", "repair_non_psd_ld", "success", "ld_repaired",
    ld_status = "repaired_unvalidated", repairable = TRUE,
    timeout_seconds = timeout_ld * defaults$ld_repair_timeout_multiplier,
    ld_z_mismatch_lambda = ld_qc_lambda,
    ld_min_eigenvalue = ld_min_eigenvalue
  )

  repaired_validation <- tryCatch(
    list(
      ok = TRUE,
      value = validate_ld_for_susie(
        ld_repaired, selected$EZ, n_eff, msg,
        eigenvalue_tolerance = defaults$ld_eigenvalue_tolerance,
        lambda_warning = defaults$ld_lambda_warning,
        lambda_failure = defaults$ld_lambda_failure
      )
    ),
    susie_ld_validation_error = function(e) list(ok = FALSE, error = e),
    error = function(e) list(
      ok = FALSE,
      error = structure(
        list(
          message = conditionMessage(e), reason = "repaired_ld_validation_internal_error",
          repairable = FALSE, lambda = NA_real_, min_eigenvalue = NA_real_
        ),
        class = c("susie_ld_validation_error", "error", "condition")
      )
    )
  )
  if (!isTRUE(repaired_validation$ok)) {
    error <- repaired_validation$error
    add_audit(
      "ld_revalidation", "validate_repaired_ld", "failed", error$reason,
      ld_status = "invalid", repairable = error$repairable,
      ld_z_mismatch_lambda = error$lambda,
      ld_min_eigenvalue = error$min_eigenvalue,
      details = error$message
    )
    msg("Repaired LD validation FAILED [", error$reason, "]: ", error$message)
    ld_qc_lambda <- error$lambda
    ld_min_eigenvalue <- error$min_eigenvalue
    return(failed_result(error$reason))
  }

  repaired_qc <- repaired_validation$value
  ld_repaired <- repaired_qc$ld
  ld_qc_lambda <- repaired_qc$lambda
  ld_min_eigenvalue <- repaired_qc$min_eigenvalue
  sample_size_qc$warning_reason <- combine_warning_reasons(
    sample_size_qc$warning_reason, repaired_qc$warning_reason
  )
  add_audit(
    "ld_revalidation", "validate_repaired_ld", "success", "repaired_ld_valid",
    ld_status = "valid_repaired", repairable = FALSE,
    ld_z_mismatch_lambda = ld_qc_lambda,
    ld_min_eigenvalue = ld_min_eigenvalue
  )

  fit_input_error <- NULL
  fit_input <- tryCatch(
    create_susie_fit_input(
      selected$EZ, ld_repaired, n_eff, "valid_repaired", safe_tag, msg
    ),
    error = function(e) {
      fit_input_error <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(fit_input)) {
    add_audit(
      "fit_input", "serialize_repaired_ld", "failed",
      "susie_fit_input_serialization_failed",
      ld_status = "valid_repaired", repairable = FALSE,
      ld_z_mismatch_lambda = ld_qc_lambda,
      ld_min_eigenvalue = ld_min_eigenvalue,
      details = fit_input_error
    )
    msg("Repaired-LD fit input serialization FAILED: ", fit_input_error)
    return(failed_result("susie_fit_input_serialization_failed"))
  }

  repaired_timeout <- timeout_susie * defaults$repaired_timeout_multiplier
  fit <- run_fit_attempt(
    "fit_repaired_ld", L, defaults$repaired_max_iter,
    repaired_timeout, "valid_repaired"
  )
  if (is_susie_fit_success(fit)) {
    msg("RECOVERY SUCCESS: repaired and revalidated LD")
    return(sus_to_out(fit, ld_repaired, "repaired_ld"))
  }

  reduced_l <- max(
    1L,
    min(
      as.integer(defaults$reduced_l_max),
      as.integer(L),
      floor(as.integer(L) / 2)
    )
  )
  retry_timeout <- timeout_susie * defaults$recovery_timeout_multiplier
  fit <- run_fit_attempt(
    "fit_repaired_ld_reduced_L", reduced_l,
    defaults$recovery_max_iter, retry_timeout, "valid_repaired"
  )
  if (is_susie_fit_success(fit)) {
    msg("RECOVERY SUCCESS: repaired LD with reduced L")
    return(sus_to_out(
      fit, ld_repaired, glue("repaired_ld_L{reduced_l}"),
      suffix = glue("_L{reduced_l}")
    ))
  }

  msg("ALL REPAIRED-LD FITTING ATTEMPTS FAILED")
  failed_result("repaired_ld_fit_recovery_exhausted")
}


# =============================================================================
# 8. MAIN DRIVER (SEQUENTIAL; your structure kept)
# =============================================================================
run_susie_finemap_parallel <- function(
  locus_file, sumstat_manifest, sample_id, ld_ref, plink, SUSIE_Analysis_folder,
  configuration_file
) {
  defaults <- load_susie_defaults(configuration_file)
  SUSIE_DEFAULTS <<- defaults
  lp_threshold <- defaults$lp_threshold
  L <- defaults$max_causal_components
  minimum_purity <- defaults$min_abs_corr
  workers <- defaults$workers
  min_ram_per_worker_gb <- defaults$min_ram_per_worker_gb
  verbose <- defaults$verbose
  timeout_ld_seconds <- defaults$ld_timeout_seconds
  timeout_susie_seconds <- defaults$susie_timeout_seconds
  recovery_audit_filename <- defaults$recovery_audit_filename
  sample_size_policy <- defaults$sample_size_policy
  sample_size_summary_statistic <- defaults$sample_size_summary_statistic
  sample_size_relative_range_warning_threshold <-
    defaults$sample_size_relative_range_warning_threshold
  maximum_variants_per_locus <- defaults$maximum_variants_per_locus
  ld_peak_matrix_multiplier <- defaults$ld_peak_matrix_multiplier
  skip_mhc <- defaults$skip_mhc
  mhc_chrom <- defaults$mhc_chromosome
  mhc_start <- defaults$mhc_start
  mhc_end <- defaults$mhc_end
  genome_build <- defaults$genome_build

  if (length(configuration_file) != 1L || !file.exists(configuration_file)) {
    stop("Resolved SuSiE configuration file is missing")
  }
  if (!is.finite(minimum_purity) || minimum_purity < 0 || minimum_purity > 1) {
    stop("minimum_purity must be between 0 and 1")
  }
  if (!is.finite(maximum_variants_per_locus) ||
      maximum_variants_per_locus < 1) {
    stop("maximum_variants_per_locus must be a positive integer")
  }
  if (!is.finite(ld_peak_matrix_multiplier) ||
      ld_peak_matrix_multiplier < 1) {
    stop("ld_peak_matrix_multiplier must be at least 1")
  }
  if (!is.finite(min_ram_per_worker_gb) || min_ram_per_worker_gb <= 0) {
    stop("min_ram_per_worker_gb must be greater than zero")
  }
  if (length(recovery_audit_filename) != 1L ||
      !grepl("^[^/\\\\]+\\.tsv$", recovery_audit_filename)) {
    stop("recovery_audit_filename must be a basename ending in .tsv")
  }
  summarize_locus_nef(
    1, defaults$sample_size_policy, defaults$sample_size_summary_statistic,
    defaults$sample_size_relative_range_warning_threshold
  )
  if (!is.finite(defaults$credible_set_coverage) ||
      defaults$credible_set_coverage <= 0 ||
      defaults$credible_set_coverage > 1) {
    stop("SuSiE credible-set coverage must be in (0, 1]")
  }
  if (!genome_build %in% defaults$supported_genome_builds) {
    stop("genome_build must be GRCh37 or GRCh38")
  }
  if (!is.finite(L) || L < 1) {
    stop("L must be a positive integer")
  }
  L <- as.integer(L)

  dir.create(SUSIE_Analysis_folder, recursive = TRUE, showWarnings = FALSE)
  for (subdir in c(
    "logs", "locus_files", "rds_files", "plots", "flames_input",
    "ld_matrix_related", "output"
  )) {
    dir.create(
      file.path(SUSIE_Analysis_folder, subdir),
      recursive = TRUE,
      showWarnings = FALSE
    )
  }

  progress_file <- file.path(SUSIE_Analysis_folder, "susie_locus_progress.tsv")
  initialize_susie_progress(progress_file)
  pipeline_stage_total <- defaults$pipeline_stage_total
  record_susie_progress(
    progress_file, "pipeline", "initialization", 1L, pipeline_stage_total,
    "completed", "Output directories and progress log initialized"
  )

  configuration <- data.table::data.table(
    parameter = c(
      "locus_file", "sumstat_manifest", "sample_id", "ld_ref", "plink",
      "output_folder", "lp_threshold", "L", "workers",
      "min_ram_per_worker_gb", "verbose", "timeout_ld_seconds",
      "timeout_susie_seconds", "skip_mhc", "mhc_chrom", "mhc_start",
      "mhc_end", "genome_build", "credible_set_coverage",
      "min_abs_corr", "main_max_iter", "recovery_max_iter",
      "repaired_max_iter", "ld_repair_max_change",
      "maximum_variants_per_locus", "ld_peak_matrix_multiplier",
      "resolved_configuration_file", "resolved_configuration_md5"
    ),
    value = as.character(c(
      locus_file, sumstat_manifest, sample_id, ld_ref, plink,
      SUSIE_Analysis_folder, lp_threshold, L, workers,
      min_ram_per_worker_gb, verbose, timeout_ld_seconds,
      timeout_susie_seconds, skip_mhc, mhc_chrom, mhc_start,
      mhc_end, genome_build, defaults$credible_set_coverage,
      defaults$min_abs_corr, defaults$main_max_iter,
      defaults$recovery_max_iter, defaults$repaired_max_iter,
      defaults$ld_repair_max_change,
      defaults$maximum_variants_per_locus,
      defaults$ld_peak_matrix_multiplier,
      normalizePath(configuration_file, mustWork = TRUE),
      unname(tools::md5sum(configuration_file))
    ))
  )
  resolved_configuration <- data.table::data.table(
    parameter = paste0("resolved.", names(defaults)),
    value = vapply(
      defaults,
      function(value) paste(as.character(value), collapse = ","),
      character(1)
    )
  )
  configuration <- data.table::rbindlist(
    list(configuration, resolved_configuration),
    use.names = TRUE
  )
  data.table::fwrite(
    configuration,
    file.path(SUSIE_Analysis_folder, "run_configuration_r.tsv"),
    sep = "\t"
  )

  # This file is an internal worker fragment.  It must describe only the
  # current run; otherwise rerunning into the same directory duplicates loci.
  index_rows_file <- file.path(
    SUSIE_Analysis_folder, "flames_input", "indexfile_rows.tsv"
  )
  if (file.exists(index_rows_file)) {
    unlink(index_rows_file)
  }

  locus_raw <- data.table::fread(locus_file)
  required_locus <- c("CHR", "START", "END")
  missing_locus <- setdiff(required_locus, names(locus_raw))
  if (length(missing_locus)) {
    stop("Locus file is missing required columns: ", paste(missing_locus, collapse = ", "))
  }
  if (!nrow(locus_raw)) {
    stop("Locus file contains no loci")
  }

  normalize_chr <- function(x) {
    value <- toupper(sub("^chr", "", as.character(x), ignore.case = TRUE))
    value[value == "X"] <- "23"
    value[value == "Y"] <- "24"
    value[value == "XY"] <- "25"
    value[value %in% c("M", "MT")] <- "26"
    sub("\\.0$", "", value)
  }
  mhc_chrom <- normalize_chr(mhc_chrom)
  locus_raw[, CHR := normalize_chr(CHR)]
  locus_raw[, `:=`(START = as.integer(START), END = as.integer(END))]
  if (anyNA(locus_raw$CHR) || any(!nzchar(locus_raw$CHR)) ||
      anyNA(locus_raw$START) || anyNA(locus_raw$END) ||
      any(locus_raw$START < defaults$minimum_position) ||
      any(locus_raw$END < locus_raw$START)) {
    stop("Locus file has invalid CHR/START/END values")
  }
  if (!"GenomicLocus" %in% names(locus_raw)) {
    locus_raw[, GenomicLocus := sprintf("chr%s:%d-%d", CHR, START, END)]
  }
  if (anyNA(locus_raw$GenomicLocus) ||
      any(!nzchar(trimws(as.character(locus_raw$GenomicLocus))))) {
    stop("GenomicLocus values must be non-empty")
  }
  if (anyDuplicated(locus_raw[, .(CHR, START, END)])) {
    stop("Locus file contains duplicate CHR/START/END rows")
  }
  locus_df <- locus_raw[, .(
    chr = as.character(CHR),
    start = as.integer(START),
    end = as.integer(END),
    genomic_locus = as.character(GenomicLocus)
  )]

  sumstat_index <- data.table::fread(sumstat_manifest)
  required_manifest <- c(
    "CHR", "START", "END", "GenomicLocus", "Filename", "n_variants"
  )
  missing_manifest <- setdiff(required_manifest, names(sumstat_index))
  if (length(missing_manifest)) {
    stop("Locus summary-statistics manifest is missing required columns: ",
         paste(missing_manifest, collapse = ", "))
  }
  sumstat_index[, CHR := normalize_chr(CHR)]
  sumstat_index[, `:=`(START = as.integer(START), END = as.integer(END))]
  if (anyNA(sumstat_index$CHR) || anyNA(sumstat_index$START) ||
      anyNA(sumstat_index$END) || anyNA(sumstat_index$Filename) ||
      any(!nzchar(sumstat_index$Filename)) ||
      anyDuplicated(sumstat_index[, .(CHR, START, END)])) {
    stop("Locus summary-statistics manifest has invalid or duplicate interval rows")
  }
  if (any(!file.exists(sumstat_index$Filename))) {
    stop("Locus summary-statistics manifest references missing file(s)")
  }
  missing_locus_inputs <- locus_df[
    !sumstat_index, on = .(chr = CHR, start = START, end = END)
  ]
  if (nrow(missing_locus_inputs)) {
    stop("Locus summary-statistics manifest does not cover every worker locus")
  }

  required_sumstat <- c("CHR", "BP", "SNP", "REF", "ALT", "EZ", "NEF", "LP")

  record_susie_progress(
    progress_file, "pipeline", "input_validation", 2L, pipeline_stage_total,
    "completed",
    sprintf(
      "Validated %d loci and an exact per-locus summary-statistics manifest",
      nrow(locus_df)
    )
  )

  data.table::fwrite(
    data.table::data.table(
      component = c("R", "susieR"),
      version = c(
        paste(R.version$major, R.version$minor, sep = "."),
        as.character(utils::packageVersion("susieR"))
      ),
      genome_build = genome_build,
      credible_set_coverage = defaults$credible_set_coverage
    ),
    file.path(SUSIE_Analysis_folder, "software_versions_r.tsv"),
    sep = "\t"
  )

  message("\t\t\tStarting SEQUENTIAL PASS 1 on ", nrow(locus_df), " loci")

  pass1 <- vector("list", nrow(locus_df))

  for (i in seq_len(nrow(locus_df))) {
    record_susie_progress(
      progress_file, "primary_loci", "fine_mapping", i - 1L, nrow(locus_df),
      "running", sprintf("Starting %s", locus_df$genomic_locus[[i]])
    )
    manifest_row <- sumstat_index[
      CHR == locus_df$chr[[i]] &
      START == locus_df$start[[i]] &
      END == locus_df$end[[i]]
    ]
    if (nrow(manifest_row) != 1L) {
      stop("Expected exactly one summary-statistics manifest row for locus ",
           locus_df$genomic_locus[[i]])
    }
    if (!identical(
      as.character(manifest_row$GenomicLocus[[1]]),
      as.character(locus_df$genomic_locus[[i]])
    )) {
      stop("Summary-statistics manifest locus provenance differs for ",
           locus_df$genomic_locus[[i]])
    }
    df <- data.table::fread(manifest_row$Filename[[1]])
    missing_sumstat <- setdiff(required_sumstat, names(df))
    if (length(missing_sumstat)) {
      stop("Exact locus summary statistics are missing required columns: ",
           paste(missing_sumstat, collapse = ", "))
    }
    df[, CHR := normalize_chr(CHR)]
    df[, `:=`(
      BP = as.integer(BP), SNP = as.character(SNP),
      REF = toupper(as.character(REF)), ALT = toupper(as.character(ALT)),
      EZ = as.numeric(EZ), NEF = as.numeric(NEF), LP = as.numeric(LP)
    )]
    invalid_sumstat <- is.na(df$CHR) | is.na(df$BP) |
      df$BP < defaults$minimum_position |
      is.na(df$SNP) | !nzchar(df$SNP) |
      is.na(df$REF) | !nzchar(df$REF) |
      is.na(df$ALT) | !nzchar(df$ALT) |
      !is.finite(df$EZ) | !is.finite(df$NEF) | df$NEF <= 0 |
      !is.finite(df$LP)
    if (any(invalid_sumstat) || anyDuplicated(df$SNP) ||
        anyDuplicated(df[, .(CHR, BP, REF, ALT)])) {
      stop("Exact locus summary-statistics file failed validation for ",
           locus_df$genomic_locus[[i]])
    }
    if (nrow(df) != as.integer(manifest_row$n_variants[[1]])) {
      stop("Exact locus summary-statistics row count differs from manifest for ",
           locus_df$genomic_locus[[i]])
    }

    pass1[[i]] <- process_locus(
      locus = locus_df[i, ],
      df = df,
      sample_id = sample_id,
      ld_ref = ld_ref,
      plink = plink,
      analysis_folder = SUSIE_Analysis_folder,
      lp_threshold = lp_threshold,
      L = L,
      timeout_ld = timeout_ld_seconds,
      timeout_susie = timeout_susie_seconds,
      skip_mhc = skip_mhc,
      mhc_chr = mhc_chrom,
      mhc_start = mhc_start,
      mhc_end = mhc_end,
      verbose = verbose,
      genome_build = genome_build,
      defaults = defaults
    )
    record_susie_progress(
      progress_file, "primary_loci", "fine_mapping", i, nrow(locus_df),
      pass1[[i]]$status,
      sprintf("Finished %s", locus_df$genomic_locus[[i]])
    )
  }

  record_susie_progress(
    progress_file, "pipeline", "primary_pass", 3L, pipeline_stage_total,
    "completed", sprintf("Primary pass completed for %d loci", nrow(locus_df))
  )

  rec_indices <- which(vapply(
    pass1,
    function(x) identical(x$status, "RECOVERY"),
    logical(1)
  ))
  rec_jobs <- pass1[rec_indices]
  message("\t\t\t", length(rec_jobs), " loci need recovery")

  rec_results <- list()
  if (length(rec_jobs)) {
    message("\t\t\tStarting SEQUENTIAL RECOVERY PASS")
    rec_results <- vector("list", length(rec_jobs))

    for (j in seq_along(rec_jobs)) {
      record_susie_progress(
        progress_file, "recovery_loci", "recovery", j - 1L, length(rec_jobs),
        "running", sprintf("Starting recovery for %s", rec_jobs[[j]]$locus$genomic_locus)
      )
      rec_results[[j]] <- recover_locus(
        job = rec_jobs[[j]],
        sample_id = sample_id,
        analysis_folder = SUSIE_Analysis_folder,
        L = L,
        timeout_ld = timeout_ld_seconds,
        timeout_susie = timeout_susie_seconds,
        verbose = verbose,
        genome_build = genome_build,
        defaults = defaults
      )
      record_susie_progress(
        progress_file, "recovery_loci", "recovery", j, length(rec_jobs),
        rec_results[[j]]$status,
        sprintf("Finished recovery for %s", rec_jobs[[j]]$locus$genomic_locus)
      )
    }
  } else {
    record_susie_progress(
      progress_file, "recovery_loci", "recovery", 0L, 0L,
      "not_required", "No loci required recovery"
    )
  }

  record_susie_progress(
    progress_file, "pipeline", "recovery_pass", 4L, pipeline_stage_total,
    "completed", sprintf("Recovery pass completed for %d loci", length(rec_jobs))
  )

  # Replace recovery placeholders in place.  Appending recovery output would
  # count each recovered locus twice and would leave the original RECOVERY row
  # in QC files.
  final_results <- pass1
  if (length(rec_indices)) {
    for (j in seq_along(rec_indices)) {
      final_results[[rec_indices[[j]]]] <- rec_results[[j]]
    }
  }

  credible_set_df <- aggregate_and_write_results(
    final_results,
    sample_id,
    SUSIE_Analysis_folder,
    verbose,
    target_coverage = defaults$credible_set_coverage,
    recovery_audit_file = file.path(
      SUSIE_Analysis_folder, defaults$recovery_audit_filename
    )
  )

  successful <- sum(vapply(
    final_results,
    function(x) x$status %in% c("OK", "RECOVERED"),
    logical(1)
  ))
  recovered <- sum(vapply(
    final_results, function(x) identical(x$status, "RECOVERED"), logical(1)
  ))
  warned <- sum(vapply(
    final_results,
    function(x) !is.null(x$sample_size_qc) &&
      !is.na(x$sample_size_qc$warning_reason),
    logical(1)
  ))
  failed <- length(final_results) - successful
  variant_limit_failures <- sum(vapply(
    final_results,
    function(x) identical(x$reason, "maximum_variants_per_locus_exceeded"),
    logical(1)
  ))
  fit_input_serializations <- sum(vapply(
    final_results,
    function(x) if (!is.null(x$fit_input_qc)) {
      as.integer(x$fit_input_qc$fit_input_serializations)
    } else {
      0L
    },
    integer(1)
  ))
  fit_child_processes <- sum(vapply(
    final_results,
    function(x) if (!is.null(x$fit_input_qc)) {
      as.integer(x$fit_input_qc$fit_child_processes)
    } else {
      0L
    },
    integer(1)
  ))
  fit_input_reuses <- sum(vapply(
    final_results,
    function(x) if (!is.null(x$fit_input_qc)) {
      as.integer(x$fit_input_qc$fit_input_reuses)
    } else {
      0L
    },
    integer(1)
  ))
  message("\nSuSiE final summary")
  message("  Loci attempted: ", length(final_results))
  message("  Primary successes: ", successful - recovered)
  message("  Recovery successes: ", recovered)
  message("  Loci with warnings: ", warned)
  message("  Failed or skipped: ", failed)
  message("  Variant-limit failures: ", variant_limit_failures)
  message("  Fit payload serializations: ", fit_input_serializations)
  message("  SuSiE child fits: ", fit_child_processes)
  message("  Fit payload reuses: ", fit_input_reuses)
  if (variant_limit_failures) {
    message(
      "  Override: --maximum-variants-per-locus COUNT or ",
      "modules.fine_mapping.ld_resource_guard.maximum_variants_per_locus in YAML"
    )
  }
  message("  Recovery audit TSV: ", file.path(
    SUSIE_Analysis_folder, defaults$recovery_audit_filename
  ))
  message("  Overall status: ", if (successful) "completed" else "failed")
  record_susie_progress(
    progress_file, "pipeline", "aggregation", pipeline_stage_total,
    pipeline_stage_total, if (successful) "completed" else "failed",
    sprintf("Finished with %d successful and %d unsuccessful loci",
            successful, length(final_results) - successful)
  )
  invisible(TRUE)
}
