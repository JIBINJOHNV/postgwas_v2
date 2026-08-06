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
# 1. RAM detection
# =============================================================================
get_free_memory_gb <- function() {
  os <- Sys.info()[["sysname"]]

  if (os == "Darwin") {
    # macOS: true "free" is not reliable; keep your conservative heuristic
    out <- suppressWarnings(system("sysctl hw.memsize", intern = TRUE, ignore.stderr = TRUE))
    if (length(out) == 0) return(Inf)
    total <- as.numeric(sub("hw.memsize: ([0-9]+)", "\\1", out))
    if (is.na(total)) return(Inf)
    return(SUSIE_DEFAULTS$macos_available_memory_fraction * total / 1024^3)
  }

  if (os == "Linux") {
    mem <- tryCatch(readLines("/proc/meminfo"), error = function(e) character(0))
    avail <- grep("^MemAvailable:", mem, value = TRUE)
    if (length(avail) == 0) return(Inf)
    kb <- suppressWarnings(as.numeric(sub(".* ([0-9]+) kB.*", "\\1", avail)))
    if (!is.finite(kb)) return(Inf)
    return(kb / 1024^2)
  }

  Inf
}

auto_detect_workers <- function(
  min_ram_per_worker_gb = SUSIE_DEFAULTS$min_ram_per_worker_gb,
  reserve_cores = SUSIE_DEFAULTS$reserve_cores,
  verbose = TRUE
) {
  phys_cores <- parallel::detectCores(logical = FALSE) %||% parallel::detectCores()
  max_by_cores <- max(1L, as.integer(phys_cores - reserve_cores))

  free_gb <- tryCatch(get_free_memory_gb(), error = function(e) Inf)
  if (!is.finite(free_gb) || free_gb < 1) {
    if (verbose) message("Free RAM unknown → using ", max_by_cores, " workers (core-based)")
    return(max_by_cores)
  }

  max_by_ram <- floor(free_gb / min_ram_per_worker_gb)
  workers <- max(1L, min(max_by_cores, max_by_ram))
  if (verbose) message("Auto workers: ~", round(free_gb, 1), " GB free → ", workers, " workers")
  workers
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

  wait_for_memory(verbose = verbose)

  idx <- which(
    df$CHR == as.character(locus$chr) &
    df$BP  >= locus$start &
    df$BP  <= locus$end
  )
  if (!length(idx))
    return(list(status = "SKIP", reason = "no_snps", locus = locus))

  selected <- data.table::as.data.table(df[idx])

  if ("LP" %in% names(selected)) {
    maxlp <- suppressWarnings(max(selected$LP, na.rm = TRUE))
  if (is.finite(maxlp) && maxlp < lp_threshold)
      return(list(status = "SKIP", reason = "low_signal", locus = locus))
  }

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
    return(list(status = "FAILED", reason = "ld_extraction_failed", locus = locus))
  }

  aligned <- tryCatch(
    align_sumstats_to_ld(selected, ld, msg),
    error = function(e) {
      msg("❌ LD allele alignment failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(aligned)) {
    return(list(status = "FAILED", reason = "allele_alignment_failed", locus = locus))
  }
  selected <- aligned$selected
  ld <- aligned$ld
  data.table::set(selected, j = "variable_index", value = seq_len(nrow(selected)))

  z <- selected$EZ
  n_eff <- suppressWarnings(median(selected$NEF, na.rm = TRUE))
  if (!is.finite(n_eff) || n_eff <= 0) {
    return(list(status = "FAILED", reason = "invalid_sample_size", locus = locus))
  }

  ld_qc <- tryCatch(
    validate_ld_for_susie(
      ld, z, n_eff, msg,
      eigenvalue_tolerance = defaults$ld_eigenvalue_tolerance,
      lambda_warning = defaults$ld_lambda_warning,
      lambda_failure = defaults$ld_lambda_failure
    ),
    error = function(e) {
      msg("❌ LD validation failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(ld_qc)) {
    return(list(
      status = "RECOVERY", reason = "ld_validation_failed", locus = locus,
      selected = selected, ld.mat = ld, n_eff = n_eff, ld_qc_lambda = NA_real_
    ))
  }
  ld <- ld_qc$ld

  # ---------------- SuSiE (HARD TIMEOUT) ----------------
  fitted <- run_susie_hard_timeout(
    z = z,
    R = ld,
    n = n_eff,
    L = L,
    max_iter = defaults$main_max_iter,
    timeout_s = timeout_susie,
    tag = safe_tag,
    log_msg = msg,
    susie_verbose = TRUE,
    coverage = defaults$credible_set_coverage,
    min_abs_corr = defaults$min_abs_corr
  )

  if (is.null(fitted) ||
      isTRUE(fitted$timeout) ||
      !inherits(fitted, "susie") ||
      !isTRUE(fitted$converged)) {
    msg("SuSiE FAILED → recovery")
    return(list(
      status = "RECOVERY",
      locus = locus,
      selected = selected,
      ld.mat = ld,
      n_eff = n_eff,
      ld_qc_lambda = ld_qc$lambda
    ))
  }

  # ---------------- SUCCESS PATH ----------------
  selected_sus <- annotate_susie(
    fitted, selected, locus, min_abs_corr = defaults$min_abs_corr
  )
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
    target_coverage = defaults$credible_set_coverage
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
    ld_qc_lambda = ld_qc$lambda
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
      locus    <- job$locus
      selected <- job$selected
      ld.mat   <- job$ld.mat
      n_eff    <- job$n_eff
      ld_qc_lambda <- job$ld_qc_lambda

      tag      <- glue("{locus$chr}_{locus$start}_{locus$end}")
      safe_tag <- make_safe_tag(tag)

      logfile <- glue("{analysis_folder}/logs/{sample_id}_{safe_tag}_RECOVERY.log")
      msg     <- make_locus_logger(logfile)

      sus_to_out <- function(fitted, ld_matrix, status, reason = " ", suffix = "") {

        selected_sus <- annotate_susie(
          fitted   = fitted,
          selected = selected,
          locus    = locus,
          min_abs_corr = defaults$min_abs_corr
        )

        plot_susie_loci_ld(
          df = selected_sus,
          ld = ld_matrix,
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
          ld     = ld_matrix,
          snp_df = selected,
          outfile = glue("{analysis_folder}/flames_input/{sample_id}_{safe_tag}"),
          genomic_locus = locus$genomic_locus,
          genome_build = genome_build,
          target_coverage = defaults$credible_set_coverage
        )

        cred_df <- selected_sus[!is.na(cs) & cs != -1]

        dir.create(glue("{analysis_folder}/locus_files"), recursive = TRUE, showWarnings = FALSE)
        dir.create(glue("{analysis_folder}/rds_files"),   recursive = TRUE, showWarnings = FALSE)

        cred_file <- glue("{analysis_folder}/locus_files/{sample_id}_{safe_tag}_SUSIE_credible_sets{suffix}.cred1")
        data.table::fwrite(cred_df, cred_file, sep = " ", quote = FALSE)

        saveRDS(
          fitted,
          glue("{analysis_folder}/rds_files/{sample_id}_{safe_tag}_SUSIE{suffix}.rds")
        )

        list(
          status    = status,
          reason    = reason,
          vars      = selected_sus,
          cs        = cred_df,
          fitted    = fitted,
          locus     = locus,
          cred_file = cred_file,
          ld_qc_lambda = ld_qc_lambda
        )
      }

      if (is.null(ld.mat)) {
        msg("RECOVERY FAILED: ld.mat is NULL")
        return(list(status = "FAILED", locus = locus, reason = "ld_mat_null"))
      }

      if (!is.finite(n_eff) || n_eff <= 0) {
        return(list(status = "FAILED", locus = locus, reason = "invalid_sample_size"))
      }

      # ============================================================
      # STEP 0 — Higher iterations
      # ============================================================
      msg("RECOVERY Step 0/4: SuSiE max_iter=", defaults$recovery_max_iter, " (L=", L, ")")

      fit0 <- run_susie_hard_timeout(
        z = selected$EZ,
        R = ld.mat,
        n = n_eff,
        L = L,
        max_iter = defaults$recovery_max_iter,
        timeout_s = timeout_susie * defaults$recovery_timeout_multiplier,
        tag = paste0(safe_tag, "_rec0"),
        log_msg = msg,
        susie_verbose = TRUE,
        coverage = defaults$credible_set_coverage,
        min_abs_corr = defaults$min_abs_corr
      )

      if (!is.null(fit0) &&
          inherits(fit0, "susie") &&
          isTRUE(fit0$converged)) {
        msg("RECOVERY SUCCESS: higher iterations")
        return(sus_to_out(fit0, ld.mat, "RECOVERED", "higher_iter"))
      }

      # ============================================================
      # STEP 1 — Reduce L
      # ============================================================
      L_step1 <- min(defaults$reduced_l_max, as.integer(L))
      msg("RECOVERY Step 1/4: SuSiE L=", L_step1)

      fit1 <- run_susie_hard_timeout(
        z = selected$EZ,
        R = ld.mat,
        n = n_eff,
        L = L_step1,
        max_iter = defaults$recovery_max_iter,
        timeout_s = timeout_susie * defaults$recovery_timeout_multiplier,
        tag = paste0(safe_tag, "_rec1"),
        log_msg = msg,
        susie_verbose = TRUE,
        coverage = defaults$credible_set_coverage,
        min_abs_corr = defaults$min_abs_corr
      )

      if (!is.null(fit1) &&
          inherits(fit1, "susie") &&
          isTRUE(fit1$converged)) {
        msg("RECOVERY SUCCESS: L=", L_step1)
        return(sus_to_out(fit1, ld.mat, "RECOVERED", paste0("L=", L_step1)))
      }

      # ============================================================
      # STEP 2 — LD repair
      # ============================================================
      msg("RECOVERY Step 2/4: LD repair")
      wait_for_memory(verbose = verbose)

      ld_repaired <- tryCatch(
        R.utils::withTimeout(
          LD_fix_fast(ld.mat, safe_tag, msg),
          timeout = timeout_ld * defaults$ld_repair_timeout_multiplier,
          onTimeout = "error"
        ),
        error = function(e) {
          msg("LD repair FAILED or TIMEOUT: ", e$message)
          NULL
        }
      )

      if (is.null(ld_repaired)) {
        msg("RECOVERY FAILED: LD repair")
        return(list(status = "FAILED", locus = locus, reason = "LD_repair_failed"))
      }

      repaired_qc <- tryCatch(
        validate_ld_for_susie(
          ld_repaired, selected$EZ, n_eff, msg,
          eigenvalue_tolerance = defaults$ld_eigenvalue_tolerance,
          lambda_warning = defaults$ld_lambda_warning,
          lambda_failure = defaults$ld_lambda_failure
        ),
        error = function(e) {
          msg("Repaired LD validation FAILED: ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(repaired_qc)) {
        return(list(status = "FAILED", locus = locus, reason = "repaired_ld_validation_failed"))
      }
      ld_repaired <- repaired_qc$ld
      ld_qc_lambda <- repaired_qc$lambda

      # ============================================================
      # STEP 3 — Repaired LD, original L
      # ============================================================
      msg("RECOVERY Step 3/4: SuSiE on repaired LD (L=", L, ")")

      fit3 <- run_susie_hard_timeout(
        z = selected$EZ,
        R = ld_repaired,
        n = n_eff,
        L = L,
        max_iter = defaults$repaired_max_iter,
        timeout_s = timeout_susie * defaults$repaired_timeout_multiplier,
        tag = paste0(safe_tag, "_rec3"),
        log_msg = msg,
        susie_verbose = TRUE,
        coverage = defaults$credible_set_coverage,
        min_abs_corr = defaults$min_abs_corr
      )

      if (!is.null(fit3) &&
          inherits(fit3, "susie") &&
          isTRUE(fit3$converged)) {
        msg("RECOVERY SUCCESS: repaired LD")
        return(sus_to_out(fit3, ld_repaired, "RECOVERED", "repaired_ld"))
      }

      # ============================================================
      # STEP 4 — Repaired LD, reduced L
      # ============================================================
      L_new <- max(1L, min(as.integer(L), floor(as.integer(L) / 2)))
      msg("RECOVERY Step 4/4: SuSiE repaired LD with L=", L_new)

      fit4 <- run_susie_hard_timeout(
        z = selected$EZ,
        R = ld_repaired,
        n = n_eff,
        L = L_new,
        max_iter = defaults$recovery_max_iter,
        timeout_s = timeout_susie * defaults$recovery_timeout_multiplier,
        tag = paste0(safe_tag, "_rec4"),
        log_msg = msg,
        susie_verbose = TRUE,
        coverage = defaults$credible_set_coverage,
        min_abs_corr = defaults$min_abs_corr
      )

      if (!is.null(fit4) &&
          inherits(fit4, "susie") &&
          isTRUE(fit4$converged)) {
        msg("RECOVERY SUCCESS: repaired LD L=", L_new)
        return(
          sus_to_out(
            fit4, ld_repaired,
            "RECOVERED",
            glue("repaired_ld_L{L_new}"),
            suffix = glue("_L{L_new}")
          )
        )
      }

      msg("ALL RECOVERY ATTEMPTS FAILED")
      list(status = "FAILED", locus = locus, reason = "all_recovery_attempts_failed")
}


# =============================================================================
# 8. MAIN DRIVER (SEQUENTIAL; your structure kept)
# =============================================================================
run_susie_finemap_parallel <- function(
  locus_file, sumstat_file, sample_id, ld_ref, plink, SUSIE_Analysis_folder,
  lp_threshold = SUSIE_DEFAULTS$lp_threshold,
  L = SUSIE_DEFAULTS$max_causal_components,
  workers = SUSIE_DEFAULTS$workers,
  min_ram_per_worker_gb = SUSIE_DEFAULTS$min_ram_per_worker_gb,
  verbose = SUSIE_DEFAULTS$verbose,
  timeout_ld_seconds = SUSIE_DEFAULTS$ld_timeout_seconds,
  timeout_susie_seconds = SUSIE_DEFAULTS$susie_timeout_seconds,
  skip_mhc = SUSIE_DEFAULTS$skip_mhc,
  mhc_chrom = SUSIE_DEFAULTS$mhc_chromosome,
  mhc_start = SUSIE_DEFAULTS$mhc_start,
  mhc_end = SUSIE_DEFAULTS$mhc_end,
  genome_build = SUSIE_DEFAULTS$genome_build,
  defaults = SUSIE_DEFAULTS
) {
  if (as.numeric(defaults$credible_set_coverage) !=
      SUSIE_DEFAULTS$credible_set_coverage) {
    stop("SuSiE credible-set coverage must be exactly 0.95")
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
      "locus_file", "sumstat_file", "sample_id", "ld_ref", "plink",
      "output_folder", "lp_threshold", "L", "workers",
      "min_ram_per_worker_gb", "verbose", "timeout_ld_seconds",
      "timeout_susie_seconds", "skip_mhc", "mhc_chrom", "mhc_start",
      "mhc_end", "genome_build", "credible_set_coverage",
      "min_abs_corr", "main_max_iter", "recovery_max_iter",
      "repaired_max_iter", "ld_repair_max_change"
    ),
    value = as.character(c(
      locus_file, sumstat_file, sample_id, ld_ref, plink,
      SUSIE_Analysis_folder, lp_threshold, L, workers,
      min_ram_per_worker_gb, verbose, timeout_ld_seconds,
      timeout_susie_seconds, skip_mhc, mhc_chrom, mhc_start,
      mhc_end, genome_build, defaults$credible_set_coverage,
      defaults$min_abs_corr, defaults$main_max_iter,
      defaults$recovery_max_iter, defaults$repaired_max_iter,
      defaults$ld_repair_max_change
    ))
  )
  default_configuration <- data.table::data.table(
    parameter = paste0("default.", names(defaults)),
    value = vapply(
      defaults,
      function(value) paste(as.character(value), collapse = ","),
      character(1)
    )
  )
  configuration <- data.table::rbindlist(
    list(configuration, default_configuration),
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

  df <- data.table::fread(sumstat_file)
  required_sumstat <- c("CHR", "BP", "SNP", "REF", "ALT", "EZ", "NEF", "LP")
  missing_sumstat <- setdiff(required_sumstat, names(df))
  if (length(missing_sumstat)) {
    stop("Summary statistics are missing required columns: ",
         paste(missing_sumstat, collapse = ", "))
  }
  if (!nrow(df)) {
    stop("Summary-statistics file contains no variants")
  }
  df[, CHR := normalize_chr(CHR)]
  df[, `:=`(
    BP = as.integer(BP),
    SNP = as.character(SNP),
    REF = toupper(as.character(REF)),
    ALT = toupper(as.character(ALT)),
    EZ = as.numeric(EZ),
    NEF = as.numeric(NEF),
    LP = as.numeric(LP)
  )]
  invalid_sumstat <- is.na(df$CHR) | is.na(df$BP) |
    df$BP < defaults$minimum_position |
    is.na(df$SNP) | !nzchar(df$SNP) |
    is.na(df$REF) | !nzchar(df$REF) |
    is.na(df$ALT) | !nzchar(df$ALT) |
    !is.finite(df$EZ) | !is.finite(df$NEF) | df$NEF <= 0 |
    !is.finite(df$LP)
  if (any(invalid_sumstat)) {
    stop("Summary statistics contain ", sum(invalid_sumstat),
         " invalid row(s); fix missing/non-finite coordinates, alleles, EZ, NEF, or LP")
  }
  if (anyDuplicated(df$SNP)) {
    stop("Summary statistics contain duplicate SNP identifiers")
  }
  if (anyDuplicated(df[, .(CHR, BP, REF, ALT)])) {
    stop("Summary statistics contain duplicate chromosome/position/allele variants")
  }

  record_susie_progress(
    progress_file, "pipeline", "input_validation", 2L, pipeline_stage_total,
    "completed",
    sprintf("Validated %d loci and %d summary-statistic variants", nrow(locus_df), nrow(df))
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
    source(UTILITIES_R, local = FALSE)

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
      source(UTILITIES_R, local = FALSE)

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
    target_coverage = defaults$credible_set_coverage
  )

  successful <- sum(vapply(
    final_results,
    function(x) x$status %in% c("OK", "RECOVERED"),
    logical(1)
  ))
  if (!successful) {
    stop("SuSiE completed without any successful loci; inspect the SuSiE QC summary and logs")
  }

  message("\t\t\tSuSiE fine-mapping finished for ", sample_id,
          " (", successful, "/", length(final_results), " successful loci)")
  record_susie_progress(
    progress_file, "pipeline", "aggregation", pipeline_stage_total,
    pipeline_stage_total, "completed",
    sprintf("Finished with %d successful and %d unsuccessful loci",
            successful, length(final_results) - successful)
  )
  invisible(TRUE)
}
