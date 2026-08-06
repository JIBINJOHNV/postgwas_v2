#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(glue)
  library(susieR)
  library(Matrix)
  library(future)
  library(future.apply)
  library(R.utils)
  library(ggplot2)
  library(ggrepel)
  library(parallel)
  library(purrr)
  library(progressr)
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
    return(0.25 * total / 1024^3)  # keep your original conservative value
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

auto_detect_workers <- function(min_ram_per_worker_gb = 4, reserve_cores = 1, verbose = TRUE) {
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
wait_for_memory <- function(threshold_used_pct = 96,
                            check_interval_sec = 10,
                            verbose = TRUE,
                            max_wait_sec = 6 * 3600) {  # 6 hours safeguard
  if (Sys.info()[["sysname"]] %in% c("Darwin", "Windows")) return(invisible(TRUE))

  t0 <- Sys.time()

  repeat {
    mem <- tryCatch(readLines("/proc/meminfo", n = 50), error = function(e) character(0))
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
LD_fix_fast <- function(ld.mat, tag = "", log_msg) {
  ld.mat <- (ld.mat + t(ld.mat)) / 2
  start_t <- Sys.time()
  log_msg("Fast LD repair started for ", tag)

  e <- eigen(ld.mat, symmetric = TRUE)
  e$values[e$values < 1e-8] <- 1e-8
  repaired <- e$vectors %*% diag(e$values) %*% t(e$vectors)

  ok <- isTRUE(tryCatch({ chol(repaired); TRUE }, error = function(e) FALSE))
  if (ok) {
    log_msg("Fast LD repair OK (", round(difftime(Sys.time(), start_t, units = "secs"), 2), "s)")
    return(repaired)
  }

  log_msg("Fast repair failed → falling back to nearPD")
  npd <- Matrix::nearPD(ld.mat, corr = TRUE, keepDiag = TRUE)
  as.matrix(npd$mat)
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
  mhc_chr, mhc_start, mhc_end, verbose = TRUE
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
    if (is.finite(maxlp) && maxlp <= lp_threshold)
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
      timeout_seconds = 600
    ),
    error = function(e) {
      msg("❌ LD extraction ERROR (attempt 1): ", conditionMessage(e))
      NULL
    }
  )

  # ---- Retry with longer timeout ----
  if (is.null(ld)) {
    msg("🔁 LD retry with longer timeout (1500s)")

    ld <- tryCatch(
      postgwas_ld_matrix(
        variants = selected$SNP,
        bfile = ld_ref,
        plink_bin = plink,
        tag = gsub(":", "-", tag),
        with_alleles = FALSE,
        logfile = glue("{analysis_folder}/ld_matrix_related/{sample_id}_{safe_tag}.log"),
        output_folder = analysis_folder,
        timeout_seconds = 1500
      ),
      error = function(e) {
        msg("❌ LD extraction ERROR (attempt 2): ", conditionMessage(e))
        NULL
      }
    )
  }


  if (is.null(ld)) {
    msg("LD extraction FAILED → recovery")
    return(list(status = "RECOVERY", locus = locus,
                selected = selected, ld.mat = NULL, n_eff = NA_real_))
  }

  m <- match(rownames(ld), selected$SNP)
  if (anyNA(m)) stop("LD SNP mismatch")

  selected <- selected[m]
  data.table::set(selected, j = "variable_index", value = seq_len(nrow(selected)))

  z <- selected$EZ
  n_eff <- suppressWarnings(median(selected$NEF, na.rm = TRUE))
  if (!is.finite(n_eff) || n_eff <= 0) n_eff <- 1e5

  # ---------------- SuSiE (HARD TIMEOUT) ----------------
  fitted <- run_susie_hard_timeout(
    z = z,
    R = ld,
    n = n_eff,
    L = L,
    max_iter = 100,
    timeout_s = timeout_susie,
    tag = safe_tag,
    log_msg = msg,
    susie_verbose = TRUE
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
      n_eff = n_eff
    ))
  }

  # ---------------- SUCCESS PATH ----------------
  selected_sus <- annotate_susie(fitted, selected, locus)
  cred_df <- selected_sus[!is.na(cs) & cs != -1]

  plot_susie_loci_ld(
    df = selected_sus,
    ld = ld,
    sample_id = sample_id,
    outdir = glue("{analysis_folder}/plots/"),
    min_pip_label = 0.01,
    bg_size = 0.01,
    bg_alpha = 0.15,
    cs_size = 1,
    cs_alpha = 0.2
  )

  generate_flames_files(
    fitted = fitted,
    ld = ld,
    snp_df = selected,
    outfile = glue("{analysis_folder}/flames_input/{sample_id}_{safe_tag}")
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
    cred_file = cred_file
  )
}

# =============================================================================
# 7. Recovery (PASS 2)
# =============================================================================
recover_locus <- function(
  job, sample_id, analysis_folder,
  L, timeout_ld, timeout_susie,
  verbose = TRUE
) {
      locus    <- job$locus
      selected <- job$selected
      ld.mat   <- job$ld.mat
      n_eff    <- job$n_eff

      tag      <- glue("{locus$chr}_{locus$start}_{locus$end}")
      safe_tag <- make_safe_tag(tag)

      logfile <- glue("{analysis_folder}/logs/{sample_id}_{safe_tag}_RECOVERY.log")
      msg     <- make_locus_logger(logfile)

      sus_to_out <- function(fitted, ld_matrix, status, reason = " ", suffix = "") {

        selected_sus <- annotate_susie(
          fitted   = fitted,
          selected = selected,
          locus    = locus
        )

        plot_susie_loci_ld(
          df = selected_sus,
          ld = ld_matrix,
          sample_id = sample_id,
          outdir = glue("{analysis_folder}/plots/"),
          min_pip_label = 0.01,
          bg_size = 0.01,
          bg_alpha = 0.15,
          cs_size = 1,
          cs_alpha = 0.2
        )

        generate_flames_files(
          fitted = fitted,
          ld     = ld_matrix,
          snp_df = selected,
          outfile = glue("{analysis_folder}/flames_input/{sample_id}_{safe_tag}")
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
          cred_file = cred_file
        )
      }

      if (is.null(ld.mat)) {
        msg("RECOVERY FAILED: ld.mat is NULL")
        return(list(status = "FAILED", locus = locus, reason = "ld_mat_null"))
      }

      if (!is.finite(n_eff) || n_eff <= 0)
        n_eff <- 1e5

      # ============================================================
      # STEP 0 — Higher iterations
      # ============================================================
      msg("RECOVERY Step 0/4: SuSiE max_iter=600 (L=", L, ")")

      fit0 <- run_susie_hard_timeout(
        z = selected$EZ,
        R = ld.mat,
        n = n_eff,
        L = L,
        max_iter = 600,
        timeout_s = timeout_susie * 4,
        tag = paste0(safe_tag, "_rec0"),
        log_msg = msg,
        susie_verbose = TRUE
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
      msg("RECOVERY Step 1/4: SuSiE L=5")

      fit1 <- run_susie_hard_timeout(
        z = selected$EZ,
        R = ld.mat,
        n = n_eff,
        L = 5,
        max_iter = 600,
        timeout_s = timeout_susie * 4,
        tag = paste0(safe_tag, "_rec1"),
        log_msg = msg,
        susie_verbose = TRUE
      )

      if (!is.null(fit1) &&
          inherits(fit1, "susie") &&
          isTRUE(fit1$converged)) {
        msg("RECOVERY SUCCESS: L=5")
        return(sus_to_out(fit1, ld.mat, "RECOVERED", "L=5"))
      }

      # ============================================================
      # STEP 2 — LD repair
      # ============================================================
      msg("RECOVERY Step 2/4: LD repair")
      wait_for_memory(verbose = verbose)

      ld_repaired <- tryCatch(
        R.utils::withTimeout(
          LD_fix_fast(ld.mat, safe_tag, msg),
          timeout = timeout_ld * 3,
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

      # ============================================================
      # STEP 3 — Repaired LD, original L
      # ============================================================
      msg("RECOVERY Step 3/4: SuSiE on repaired LD (L=", L, ")")

      fit3 <- run_susie_hard_timeout(
        z = selected$EZ,
        R = ld_repaired,
        n = n_eff,
        L = L,
        max_iter = 300,
        timeout_s = timeout_susie * 3,
        tag = paste0(safe_tag, "_rec3"),
        log_msg = msg,
        susie_verbose = TRUE
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
      L_new <- max(3L, floor(L / 2))
      msg("RECOVERY Step 4/4: SuSiE repaired LD with L=", L_new)

      fit4 <- run_susie_hard_timeout(
        z = selected$EZ,
        R = ld_repaired,
        n = n_eff,
        L = L_new,
        max_iter = 600,
        timeout_s = timeout_susie * 4,
        tag = paste0(safe_tag, "_rec4"),
        log_msg = msg,
        susie_verbose = TRUE
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
  lp_threshold = 7.3, L = 10, workers = "auto", min_ram_per_worker_gb = 4,
  verbose = TRUE, timeout_ld_seconds = 600, timeout_susie_seconds = 900,
  skip_mhc = TRUE, mhc_chrom = "6", mhc_start = 25e6, mhc_end = 35e6
) {
  dir.create(SUSIE_Analysis_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(SUSIE_Analysis_folder, "logs"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(SUSIE_Analysis_folder, "locus_files"), recursive = TRUE, showWarnings = FALSE)

  locus_df <- data.table::fread(locus_file)[, .(chr = as.character(CHR), start = as.integer(START), end = as.integer(END))]
  df <- data.table::fread(sumstat_file)
  df[, CHR := as.character(CHR)]

  message("\t\t\tStarting SEQUENTIAL PASS 1 on ", nrow(locus_df), " loci")

  pass1 <- vector("list", nrow(locus_df))

  for (i in seq_len(nrow(locus_df))) {
    message(sprintf("\t\t\tProcessing locus %d / %d", i, nrow(locus_df)))
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
      verbose = verbose
    )
  }

  rec_jobs <- purrr::keep(pass1, ~ identical(.x$status, "RECOVERY"))
  message("\t\t\t", length(rec_jobs), " loci need recovery")

  rec_results <- list()
  if (length(rec_jobs)) {
    message("\t\t\tStarting SEQUENTIAL RECOVERY PASS")
    rec_results <- vector("list", length(rec_jobs))

    for (j in seq_along(rec_jobs)) {
      message(sprintf("\t\t\tRecovering locus %d / %d", j, length(rec_jobs)))
      source(UTILITIES_R, local = FALSE)

      rec_results[[j]] <- recover_locus(
        job = rec_jobs[[j]],
        sample_id = sample_id,
        analysis_folder = SUSIE_Analysis_folder,
        L = L,
        timeout_ld = timeout_ld_seconds,
        timeout_susie = timeout_susie_seconds,
        verbose = verbose
      )
    }
  }

  credible_set_df <- aggregate_and_write_results(
    c(pass1, rec_results),
    sample_id,
    SUSIE_Analysis_folder,
    verbose
  )

  message("\t\t\tSuSiE fine-mapping finished for ", sample_id)
  invisible(TRUE)
}
