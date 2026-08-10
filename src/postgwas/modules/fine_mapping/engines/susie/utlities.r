
# You need to install once: install.packages("processx")
library(processx)

# Optional: if you want to use this version with timeout, install once:
# install.packages("processx")
# library(processx)  # Load only when using the timeout feature
postgwas_ld_matrix <- function(
  variants,
  bfile,
  plink_bin,
  output_folder,
  tag,
  with_alleles = TRUE,
  logfile = NULL,
  timeout_seconds = SUSIE_DEFAULTS$plink_timeout_seconds
) {
  if (length(variants) == 0) {
    warning("LD skipped: no variants provided")
    return(NULL)
  }

  # -------------------------------------------------------
  # Create output subfolder for LD
  # -------------------------------------------------------
  ld_dir <- file.path(output_folder, "ld_matrix_related")
  dir.create(ld_dir, recursive = TRUE, showWarnings = FALSE)

  safe_tag <- gsub(":", "-", tag)
  prefix <- file.path(ld_dir, paste0("ld_run_", safe_tag))

  if (is.null(logfile)) {
    logfile <- file.path(ld_dir, paste0("ld_run_", safe_tag, ".log"))
  }
  file.create(logfile)

  # -------------------------------------------------------
  # Write SNP list file
  # -------------------------------------------------------
  snp_list_file <- paste0(prefix, "_snplist.txt")
  write.table(
    data.frame(variants),
    file = snp_list_file,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  version_text <- tryCatch(
    paste(system2(plink_bin, "--version", stdout = TRUE, stderr = TRUE), collapse = " "),
    error = function(e) basename(plink_bin)
  )
  is_plink2 <- grepl("PLINK v?2|plink2", version_text, ignore.case = TRUE)

  run_plink_safe <- function(args, step_name) {
    res <- tryCatch(
      processx::run(
        command = plink_bin,
        args = args,
        timeout = timeout_seconds,
        error_on_status = FALSE
      ),
      error = function(e) {
        cat("\n", step_name, " failed: ", conditionMessage(e), "\n",
            file = logfile, append = TRUE, sep = "")
        NULL
      }
    )
    if (is.null(res)) return(FALSE)
    cat("\n=== ", step_name, " ===\n", res$stdout, "\n", res$stderr, "\n",
        file = logfile, append = TRUE, sep = "")
    if (!identical(res$status, 0L)) {
      warning(step_name, " failed (exit=", res$status, ") | see log: ", logfile)
      return(FALSE)
    }
    TRUE
  }

  # -------------------------------------------------------
  # PLINK call 1 → generate .bim
  # -------------------------------------------------------
  args1 <- c(
    "--bfile", bfile,
    "--extract", snp_list_file,
    "--make-just-bim",
    if (!is_plink2) "--keep-allele-order" else character(0),
    "--out", prefix
  )

  if (!run_plink_safe(args1, "PLINK .bim generation")) {
    return(NULL)
  }

  bim_file <- paste0(prefix, ".bim")
  if (!file.exists(bim_file) || file.size(bim_file) == 0) {
    warning("PLINK produced no .bim file | see log: ", logfile)
    return(NULL)
  }

  bim <- tryCatch(
    read.table(bim_file, stringsAsFactors = FALSE),
    error = function(e) {
      warning("Failed to read .bim: ", e$message)
      NULL
    }
  )
  if (is.null(bim) || nrow(bim) == 0) {
    warning("Empty .bim after PLINK | see log: ", logfile)
    return(NULL)
  }

  # -------------------------------------------------------
  # PLINK call 2 → LD matrix
  # -------------------------------------------------------
  args2 <- if (is_plink2) {
    c(
      "--bfile", bfile,
      "--extract", snp_list_file,
      "--r-unphased", "square", "ref-based",
      "--out", prefix
    )
  } else {
    c(
      "--bfile", bfile,
      "--extract", snp_list_file,
      "--r", "square",
      "--keep-allele-order",
      "--out", prefix
    )
  }

  if (!run_plink_safe(args2, "PLINK LD calculation")) {
    return(NULL)
  }

  if (is_plink2) {
    ld_candidates <- Sys.glob(paste0(prefix, "*.vcor1"))
    ld_candidates <- ld_candidates[!grepl("\\.vars$", ld_candidates)]
    ld_file <- if (length(ld_candidates)) ld_candidates[[1L]] else ""
  } else {
    ld_file <- paste0(prefix, ".ld")
  }
  if (!file.exists(ld_file) || file.size(ld_file) == 0) {
    warning("PLINK produced no LD file | see log: ", logfile)
    return(NULL)
  }

  # -------------------------------------------------------
  # Read LD matrix
  # -------------------------------------------------------
  res <- tryCatch(
    as.matrix(read.table(ld_file, header = FALSE)),
    error = function(e) {
      warning("Failed to read LD matrix: ", e$message)
      NULL
    }
  )
  if (is.null(res)) {
    return(NULL)
  }

  if (nrow(res) != ncol(res)) {
    warning("PLINK LD matrix is not square | see log: ", logfile)
    return(NULL)
  }

  if (is_plink2) {
    vars_file <- paste0(ld_file, ".vars")
    if (file.exists(vars_file)) {
      matrix_ids <- scan(vars_file, what = character(), quiet = TRUE)
      ord <- match(matrix_ids, bim$V2)
      if (anyNA(ord)) {
        warning("PLINK2 LD .vars file does not match the extracted BIM")
        return(NULL)
      }
      bim <- bim[ord, , drop = FALSE]
    }
  }

  if (nrow(res) != nrow(bim)) {
    warning("LD dimension does not equal extracted BIM variant count")
    return(NULL)
  }

  counted_allele <- if (is_plink2) bim$V6 else bim$V5
  other_allele   <- if (is_plink2) bim$V5 else bim$V6

  if (with_alleles) {
    rownames(res) <- colnames(res) <- paste(bim$V2, counted_allele, other_allele, sep = "_")
  } else {
    rownames(res) <- colnames(res) <- bim$V2
  }

  attr(res, "variant_map") <- data.frame(
    SNP = as.character(bim$V2),
    CHR = as.character(bim$V1),
    BP = as.integer(bim$V4),
    counted_allele = toupper(as.character(counted_allele)),
    other_allele = toupper(as.character(other_allele)),
    correlation_coding = if (is_plink2) "PLINK2_REF" else "PLINK1_A1",
    stringsAsFactors = FALSE
  )
  attr(res, "plink_version") <- version_text

  # -------------------------------------------------------
  # Cleanup
  # -------------------------------------------------------
  unlink(Sys.glob(paste0(prefix, "*")), force = TRUE)

  res
}


align_sumstats_to_ld <- function(selected, ld, log_msg) {
  selected <- data.table::copy(data.table::as.data.table(selected))
  variant_map <- attr(ld, "variant_map")
  if (is.null(variant_map)) stop("LD matrix is missing allele metadata")

  ord <- match(rownames(ld), selected$SNP)
  if (anyNA(ord)) stop("LD variants cannot be matched to summary statistics")
  selected <- selected[ord]
  variant_map <- variant_map[match(selected$SNP, variant_map$SNP), , drop = FALSE]
  if (anyNA(variant_map$SNP)) stop("LD allele map is incomplete")

  ea_col <- intersect(c("EA", "effect_allele", "A1", "ALT"), names(selected))
  oa_col <- intersect(c("OA", "non_effect_allele", "A2", "REF"), names(selected))
  if (!length(ea_col) || !length(oa_col)) {
    stop("Summary statistics must contain effect and non-effect allele columns")
  }
  ea_col <- ea_col[[1L]]
  oa_col <- oa_col[[1L]]

  ea <- toupper(as.character(selected[[ea_col]]))
  oa <- toupper(as.character(selected[[oa_col]]))
  counted <- toupper(as.character(variant_map$counted_allele))
  other <- toupper(as.character(variant_map$other_allele))

  complement <- function(x) {
    chartr("ACGT", "TGCA", x)
  }
  is_snv <- nchar(ea) == 1L & nchar(oa) == 1L &
    nchar(counted) == 1L & nchar(other) == 1L &
    ea %in% c("A", "C", "G", "T") & oa %in% c("A", "C", "G", "T")
  palindromic <- paste0(ea, oa) %in% c("AT", "TA", "CG", "GC")

  same <- ea == counted & oa == other
  swapped <- ea == other & oa == counted
  complement_same <- is_snv & !palindromic & complement(ea) == counted & complement(oa) == other
  complement_swapped <- is_snv & !palindromic & complement(ea) == other & complement(oa) == counted

  normalize_chr <- function(x) {
    value <- toupper(sub("^chr", "", as.character(x), ignore.case = TRUE))
    value[value == "X"] <- "23"
    value[value == "Y"] <- "24"
    value[value == "XY"] <- "25"
    value[value %in% c("M", "MT")] <- "26"
    sub("\\.0$", "", value)
  }
  position_ok <- normalize_chr(selected$CHR) == normalize_chr(variant_map$CHR) &
    as.integer(selected$BP) == as.integer(variant_map$BP)
  # Without effect-allele frequency in both datasets, A/T and C/G strand
  # orientation is not identifiable.  Dropping these SNPs is safer than a
  # potentially reversed z-score/LD sign.
  keep <- position_ok & !palindromic &
    (same | swapped | complement_same | complement_swapped)

  n_drop <- sum(!keep)
  if (n_drop > 0L) {
    log_msg("Allele/position QC removed ", n_drop,
            " of ", length(keep), " variants before SuSiE")
  }
  if (!any(keep)) stop("No variants remain after LD allele alignment")

  flip <- swapped | complement_swapped
  flip <- flip & keep
  selected[, EZ := as.numeric(EZ)]
  selected[flip, EZ := -EZ]

  selected <- selected[keep]
  variant_map <- variant_map[keep, , drop = FALSE]
  ld_aligned <- ld[keep, keep, drop = FALSE]

  # After sign alignment, ALT is the effect/counted allele and REF is the other allele.
  # Use data.table setters throughout. Base `[[<-` can replace data.table's
  # internal attributes after row filtering, making the next in-place update
  # fail on current data.table releases.
  data.table::set(selected, j = "ALT", value = variant_map$counted_allele)
  data.table::set(selected, j = "REF", value = variant_map$other_allele)
  data.table::set(selected, j = ea_col, value = variant_map$counted_allele)
  data.table::set(selected, j = oa_col, value = variant_map$other_allele)

  rownames(ld_aligned) <- colnames(ld_aligned) <- selected$SNP
  attr(ld_aligned, "variant_map") <- variant_map
  attr(ld_aligned, "plink_version") <- attr(ld, "plink_version")
  log_msg("LD allele alignment retained ", nrow(selected),
          " variants; flipped ", sum(flip), " z-score signs")

  list(selected = selected, ld = ld_aligned, n_flipped = sum(flip), n_dropped = n_drop)
}


stop_susie_ld_validation <- function(
  reason, message, repairable = FALSE,
  lambda = NA_real_, min_eigenvalue = NA_real_
) {
  condition <- structure(
    list(
      message = as.character(message), call = NULL,
      reason = as.character(reason), repairable = isTRUE(repairable),
      lambda = as.numeric(lambda), min_eigenvalue = as.numeric(min_eigenvalue)
    ),
    class = c("susie_ld_validation_error", "error", "condition")
  )
  stop(condition)
}


validate_ld_for_susie <- function(
  ld, z, n, log_msg,
  eigenvalue_tolerance = SUSIE_DEFAULTS$ld_eigenvalue_tolerance,
  lambda_warning = SUSIE_DEFAULTS$ld_lambda_warning,
  lambda_failure = SUSIE_DEFAULTS$ld_lambda_failure
) {
  if (!is.matrix(ld) || nrow(ld) != ncol(ld)) {
    stop_susie_ld_validation("ld_not_square", "LD must be a square matrix")
  }
  if (nrow(ld) != length(z)) {
    stop_susie_ld_validation(
      "ld_dimension_mismatch", "LD dimension does not match z-score length"
    )
  }
  if (any(!is.finite(ld))) {
    stop_susie_ld_validation("ld_nonfinite", "LD matrix contains non-finite values")
  }
  if (any(!is.finite(z))) {
    stop_susie_ld_validation("z_nonfinite", "z-scores contain non-finite values")
  }

  ld <- (ld + t(ld)) / 2
  diag_values <- diag(ld)
  if (any(!is.finite(diag_values)) || any(diag_values <= 0)) {
    stop_susie_ld_validation(
      "ld_invalid_diagonal", "LD matrix has invalid diagonal values"
    )
  }
  scale_values <- sqrt(diag_values)
  ld <- ld / tcrossprod(scale_values)
  diag(ld) <- 1

  min_eigenvalue <- min(eigen(ld, symmetric = TRUE, only.values = TRUE)$values)
  log_msg("LD minimum eigenvalue = ", signif(min_eigenvalue, 5))
  if (!is.finite(min_eigenvalue) || min_eigenvalue < -eigenvalue_tolerance) {
    stop_susie_ld_validation(
      "ld_non_psd", "LD matrix is not positive semidefinite",
      repairable = is.finite(min_eigenvalue),
      min_eigenvalue = min_eigenvalue
    )
  }

  lambda <- NA_real_
  warning_reason <- NA_character_
  if (requireNamespace("susieR", quietly = TRUE) &&
      "estimate_s_rss" %in% getNamespaceExports("susieR")) {
    lambda <- tryCatch(
      susieR::estimate_s_rss(z, ld, n = n),
      error = function(e) {
        log_msg("LD/z consistency diagnostic unavailable: ", conditionMessage(e))
        NA_real_
      }
    )
  }
  if (is.finite(lambda)) {
    log_msg("LD/z consistency lambda = ", signif(lambda, 5))
    if (lambda > lambda_failure) {
      stop_susie_ld_validation(
        "ld_z_mismatch",
        paste0("Severe LD/z mismatch detected (lambda > ", lambda_failure, ")"),
        lambda = lambda,
        min_eigenvalue = min_eigenvalue
      )
    } else if (lambda > lambda_warning) {
      warning_reason <- "ld_z_mismatch_warning"
      log_msg("WARNING: material LD/z mismatch detected (lambda > ",
              lambda_warning, ")")
    }
  }

  list(
    ld = ld, lambda = lambda, min_eigenvalue = min_eigenvalue,
    warning_reason = warning_reason
  )
}


susie_fit_failure <- function(reason, message) {
  structure(
    list(reason = as.character(reason), message = as.character(message)),
    class = "susie_fit_failure"
  )
}


is_susie_fit_success <- function(value) {
  inherits(value, "susie") && isTRUE(value$converged)
}


create_susie_fit_input <- function(z, R, n, matrix_state, tag, log_msg) {
  z <- as.numeric(z)
  n <- as.numeric(n)
  if (!is.matrix(R) || nrow(R) != ncol(R)) {
    stop("SuSiE fit input LD must be a square matrix")
  }
  if (!length(z) || nrow(R) != length(z)) {
    stop("SuSiE fit input LD dimension must match the z-score length")
  }
  if (length(n) != 1L || !is.finite(n) || n <= 0) {
    stop("SuSiE fit input sample size must be one finite positive value")
  }
  if (any(!is.finite(z)) || any(!is.finite(R))) {
    stop("SuSiE fit input contains non-finite values")
  }
  if (length(matrix_state) != 1L || is.na(matrix_state) || !nzchar(matrix_state)) {
    stop("SuSiE fit input matrix state must be non-empty")
  }

  input_path <- tempfile(
    paste0("susie_fit_", make.names(tag), "_"), fileext = ".rds"
  )
  write_ok <- FALSE
  on.exit({
    if (!write_ok && file.exists(input_path)) unlink(input_path, force = TRUE)
  }, add = TRUE)
  started <- Sys.time()
  saveRDS(
    list(z = z, R = R, n = n, matrix_state = as.character(matrix_state)),
    input_path
  )
  info <- file.info(input_path)
  if (!file.exists(input_path) || is.na(info$size) || info$size <= 0) {
    stop("SuSiE fit input serialization produced no readable payload")
  }
  if (!isTRUE(Sys.chmod(input_path, mode = "0444"))) {
    stop("SuSiE fit input could not be made read-only")
  }
  write_ok <- TRUE

  handle <- new.env(parent = emptyenv())
  handle$path <- normalizePath(input_path, mustWork = TRUE)
  handle$matrix_state <- as.character(matrix_state)
  handle$variant_count <- as.integer(length(z))
  handle$n_eff <- n
  handle$file_size_bytes <- as.numeric(info$size)
  handle$serialization_seconds <- as.numeric(difftime(
    Sys.time(), started, units = "secs"
  ))
  handle$serializations <- 1L
  handle$child_processes <- 0L
  handle$cleaned <- FALSE
  class(handle) <- c("susie_fit_input", "environment")
  log_msg(
    "SuSiE fit payload serialized once | LD state=", handle$matrix_state,
    " p=", handle$variant_count,
    " bytes=", format(handle$file_size_bytes, scientific = FALSE),
    " elapsed=", signif(handle$serialization_seconds, 5), "s"
  )
  handle
}


validate_susie_fit_input <- function(
  fit_input, expected_matrix_state = NULL,
  expected_variant_count = NULL, expected_n_eff = NULL
) {
  if (!inherits(fit_input, "susie_fit_input")) {
    stop("SuSiE fit input handle is missing or invalid")
  }
  if (isTRUE(fit_input$cleaned) || !file.exists(fit_input$path)) {
    stop("SuSiE fit input payload is unavailable")
  }
  info <- file.info(fit_input$path)
  if (is.na(info$size) || info$size <= 0) {
    stop("SuSiE fit input payload is empty")
  }
  if (!is.finite(fit_input$variant_count) || fit_input$variant_count < 1L ||
      !is.finite(fit_input$n_eff) || fit_input$n_eff <= 0) {
    stop("SuSiE fit input handle metadata is invalid")
  }
  if (!is.null(expected_matrix_state) &&
      !identical(fit_input$matrix_state, as.character(expected_matrix_state))) {
    stop("SuSiE fit input LD state does not match the requested recovery route")
  }
  if (!is.null(expected_variant_count) &&
      !identical(fit_input$variant_count, as.integer(expected_variant_count))) {
    stop("SuSiE fit input variant count does not match the recovery data")
  }
  if (!is.null(expected_n_eff) &&
      !identical(fit_input$n_eff, as.numeric(expected_n_eff))) {
    stop("SuSiE fit input sample size does not match the recovery data")
  }
  invisible(TRUE)
}


summarize_susie_fit_input <- function(fit_input) {
  if (!inherits(fit_input, "susie_fit_input")) {
    return(list(
      fit_input_matrix_state = NA_character_,
      fit_input_serializations = 0L,
      fit_child_processes = 0L,
      fit_input_reuses = 0L,
      fit_input_file_size_bytes = NA_real_,
      fit_input_serialization_seconds = NA_real_
    ))
  }
  list(
    fit_input_matrix_state = fit_input$matrix_state,
    fit_input_serializations = as.integer(fit_input$serializations),
    fit_child_processes = as.integer(fit_input$child_processes),
    fit_input_reuses = as.integer(max(0L, fit_input$child_processes - 1L)),
    fit_input_file_size_bytes = as.numeric(fit_input$file_size_bytes),
    fit_input_serialization_seconds = as.numeric(fit_input$serialization_seconds)
  )
}


susie_fit_input_audit_detail <- function(fit_input) {
  qc <- summarize_susie_fit_input(fit_input)
  paste0(
    "fit_input_matrix_state=", qc$fit_input_matrix_state,
    "; fit_input_serializations=", qc$fit_input_serializations,
    "; fit_child_processes=", qc$fit_child_processes,
    "; fit_input_reuses=", qc$fit_input_reuses,
    "; fit_input_file_size_bytes=", qc$fit_input_file_size_bytes
  )
}


cleanup_susie_fit_input <- function(fit_input, log_msg) {
  if (!inherits(fit_input, "susie_fit_input") || isTRUE(fit_input$cleaned)) {
    return(invisible(TRUE))
  }
  if (file.exists(fit_input$path)) {
    try(Sys.chmod(fit_input$path, mode = "0600"), silent = TRUE)
    unlink(fit_input$path, force = TRUE)
  }
  fit_input$cleaned <- TRUE
  removed <- !file.exists(fit_input$path)
  log_msg(
    "SuSiE fit payload cleanup | LD state=", fit_input$matrix_state,
    " child processes=", fit_input$child_processes,
    " reuses=", max(0L, fit_input$child_processes - 1L),
    " status=", if (removed) "removed" else "failed"
  )
  if (!removed) warning("Unable to remove SuSiE fit input payload: ", fit_input$path)
  invisible(removed)
}


run_susie_hard_timeout <- function(
  fit_input, L, max_iter,
  timeout_s,
  tag,
  log_msg,
  susie_verbose = TRUE,
  coverage = SUSIE_DEFAULTS$credible_set_coverage,
  min_abs_corr = SUSIE_DEFAULTS$min_abs_corr,
  poll_seconds = SUSIE_DEFAULTS$process_poll_seconds,
  terminate_grace_seconds = SUSIE_DEFAULTS$process_terminate_grace_seconds,
  stderr_tail_lines = SUSIE_DEFAULTS$stderr_tail_lines
) {
  if (as.numeric(coverage) != SUSIE_DEFAULTS$credible_set_coverage) {
    stop("SuSiE fit coverage differs from the resolved YAML configuration")
  }
  # Absolute safety net — nothing escapes this function
  tryCatch({

    stopifnot(requireNamespace("processx", quietly = TRUE))

    validate_susie_fit_input(fit_input)
    attempt_number <- fit_input$child_processes + 1L

    tmp_control <- tempfile(paste0("susie_control_", tag, "_"), fileext = ".rds")
    tmp_output <- tempfile(paste0("susie_out_",   tag, "_"), fileext = ".rds")
    tmp_script <- tempfile(paste0("susie_run_",   tag, "_"), fileext = ".R")

    tmp_stdout <- tempfile(paste0("susie_stdout_", tag, "_"), fileext = ".log")
    tmp_stderr <- tempfile(paste0("susie_stderr_", tag, "_"), fileext = ".log")

    on.exit({
      unlink(
        c(tmp_control, tmp_output, tmp_script, tmp_stdout, tmp_stderr),
        force = TRUE
      )
    }, add = TRUE)

    saveRDS(list(
      L = L,
      max_iter = max_iter,
      susie_verbose = susie_verbose,
      coverage = coverage,
      min_abs_corr = min_abs_corr,
      expected_matrix_state = fit_input$matrix_state,
      expected_variant_count = fit_input$variant_count,
      expected_n_eff = fit_input$n_eff
    ), tmp_control)

    writeLines(c(
      "suppressPackageStartupMessages(library(susieR))",
      sprintf("inp <- readRDS(%s)", deparse(fit_input$path)),
      sprintf("ctl <- readRDS(%s)", deparse(tmp_control)),
      "if (!is.list(inp) || !is.matrix(inp$R) || nrow(inp$R) != ncol(inp$R) ||",
      "    nrow(inp$R) != length(inp$z) || length(inp$n) != 1L ||",
      "    any(!is.finite(inp$R)) || any(!is.finite(inp$z)) ||",
      "    !is.finite(inp$n) || inp$n <= 0) {",
      "  stop('Serialized SuSiE fit input failed child-process validation')",
      "}",
      "if (!identical(inp$matrix_state, ctl$expected_matrix_state) ||",
      "    !identical(as.integer(length(inp$z)), ctl$expected_variant_count) ||",
      "    !identical(as.numeric(inp$n), ctl$expected_n_eff)) {",
      "  stop('Serialized SuSiE fit input does not match its attempt controls')",
      "}",
      "fit <- susieR::susie_rss(",
      "  z = inp$z,",
      "  R = inp$R,",
      "  n = inp$n,",
      "  L = ctl$L,",
      "  max_iter = ctl$max_iter,",
      "  coverage = ctl$coverage,",
      "  min_abs_corr = ctl$min_abs_corr,",
      "  return_correlation = TRUE,",
      "  verbose = ctl$susie_verbose",
      ")",
      sprintf("saveRDS(fit, %s)", deparse(tmp_output))
    ), tmp_script)

    log_msg(
      "SuSiE START | p=", fit_input$variant_count,
      " L=", L,
      " coverage=", coverage,
      " max_iter=", max_iter,
      " n_eff=", round(fit_input$n_eff, 2),
      " timeout=", timeout_s, "s",
      " LD state=", fit_input$matrix_state,
      " payload attempt=", attempt_number,
      if (attempt_number > 1L) " (reused)" else " (initial read)"
    )

    rscript <- file.path(R.home("bin"), "Rscript")
    if (!file.exists(rscript)) {
      log_msg("❌ SuSiE FAILED: current R installation has no Rscript executable")
      return(susie_fit_failure(
        "rscript_not_found",
        "Current R installation has no Rscript executable"
      ))
    }

    t0 <- Sys.time()

    px <- processx::process$new(
      command = rscript,
      args = c("--vanilla", tmp_script),
      stdout = tmp_stdout,
      stderr = tmp_stderr,
      cleanup = TRUE
    )
    fit_input$child_processes <- attempt_number

    repeat {
      if (!px$is_alive()) break

      if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > timeout_s) {
        log_msg("❌ SuSiE HARD TIMEOUT → terminating PID ", px$get_pid())

        try(px$terminate(), silent = TRUE)
        Sys.sleep(terminate_grace_seconds)

        if (px$is_alive()) {
          log_msg("❌ SuSiE still alive → killing PID ", px$get_pid())
          try(px$kill(), silent = TRUE)
        }

        return(susie_fit_failure(
          "susie_timeout", paste0("SuSiE exceeded ", timeout_s, " seconds")
        ))
      }

      Sys.sleep(poll_seconds)
    }

    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
    log_msg("SuSiE FINISHED in ", elapsed, "s")

    if (!file.exists(tmp_output)) {
      log_msg("❌ SuSiE FAILED: no output produced")

      if (file.size(tmp_stderr) > 0) {
        err <- readLines(tmp_stderr, warn = FALSE)
        log_msg("SuSiE stderr (tail):\n",
                paste(tail(err, stderr_tail_lines), collapse = "\n"))
      }

      return(susie_fit_failure("susie_no_output", "SuSiE produced no output"))
    }

    fit <- readRDS(tmp_output)

    if (!inherits(fit, "susie") || !isTRUE(fit$converged)) {
      log_msg("❌ SuSiE FAILED: did not converge")
      return(susie_fit_failure("susie_nonconvergence", "SuSiE did not converge"))
    }

    fit

  }, error = function(e) {
    # This guarantees the worker NEVER dies
    log_msg("❌ SuSiE INTERNAL ERROR: ", conditionMessage(e))
    susie_fit_failure("susie_internal_error", conditionMessage(e))
  })
}

annotate_susie <- function(
  fitted,
  selected,
  locus,
  min_abs_corr = SUSIE_DEFAULTS$min_abs_corr
) {
    suppressPackageStartupMessages({ require(data.table) }) 

    # ======================================================
    # 0. Convert to data.table
    # ======================================================
    selected <- as.data.table(selected)
    J <- nrow(selected)

    # ======================================================
    # 1. Extract SuSiE summary (vars + cs)
    # ======================================================
    sus <- summary(fitted)

    vars <- if (!is.null(sus$vars)) as.data.table(sus$vars) else data.table()
    cs   <- if (!is.null(sus$cs))   as.data.table(sus$cs)   else data.table()

    # ======================================================
    # 2. Link SuSiE vars back to GWAS SNPs
    # ======================================================
    if (nrow(vars) > 0 && "variable" %in% names(vars)) {

        setnames(vars, "variable", "variable_index")
        vars[, variable_index := as.integer(variable_index)]
        if (anyNA(vars$variable_index) ||
            any(vars$variable_index < 1L | vars$variable_index > J)) {
          stop("SuSiE summary contains invalid variable indices")
        }

        vars[, `:=`(
        SNP = selected$SNP[variable_index],
        CHR = selected$CHR[variable_index],
        BP  = selected$BP[variable_index],
        REF = selected$REF[variable_index],
        ALT = selected$ALT[variable_index]
        )]

        if ("effect" %in% colnames(vars)) {
        vars[, effect_component := effect]
        }
    } else if (nrow(vars) > 0) {
        stop("SuSiE summary variables are missing the variable index column")
    }

    # ======================================================
    # 3. Merge selected SNPs with SuSiE vars
    # ======================================================
    if (nrow(vars) > 0) {
        selected_sus <- merge(
            selected, vars,
            by = c("variable_index", "SNP", "CHR", "BP", "REF", "ALT"),
            all.x = TRUE,
            sort = FALSE
        )
    } else {
        selected_sus <- data.table::copy(selected)
        selected_sus[, `:=`(
            variable_prob = if (!is.null(fitted$pip)) fitted$pip[variable_index] else 0,
            cs = NA_integer_
        )]
    }
    data.table::setorder(selected_sus, variable_index)
    if (!"variable_prob" %in% names(selected_sus)) {
      selected_sus[, variable_prob := if (!is.null(fitted$pip)) {
        fitted$pip[variable_index]
      } else {
        0
      }]
    }
    if (!"cs" %in% names(selected_sus)) {
      selected_sus[, cs := NA_integer_]
    }
    # ======================================================
    # 4. Add locus metadata
    # ======================================================
    selected_sus[, `:=`(
        locus_chr   = locus$chr[1L],
        locus_start = locus$start[1L],
        locus_end   = locus$end[1L]
    )]
    # ======================================================
    # 5. Add credible set assignment
    # ======================================================
    cs_index <- rep(NA_character_, J)

    if (!is.null(fitted$sets$cs)) {
        for (cs_pos in seq_along(fitted$sets$cs)) {
        cs_name <- names(fitted$sets$cs)[cs_pos]
        if (is.null(cs_name) || is.na(cs_name) || !nzchar(cs_name)) {
          cs_name <- paste0("L", cs_pos)
        }
        snp_idxs <- fitted$sets$cs[[cs_pos]]
        cs_index[snp_idxs] <- cs_name
        }
    }
    selected_sus[, credible_set := cs_index[variable_index]]

    # ======================================================
    # 6. Add credible set purity (min/mean/median LD)
    # ======================================================
    if (!is.null(fitted$sets$purity)) {

        pur <- as.data.frame(fitted$sets$purity)

        purity_dt <- data.table(
        cs            = rownames(pur),
        purity_min    = pur$min.abs.corr,
        purity_mean   = pur$mean.abs.corr,
        purity_median = pur$median.abs.corr,
        is_pure       = pur$min.abs.corr > min_abs_corr
        )

        selected_sus <- merge(
        selected_sus,
        purity_dt,
        by.x = "credible_set",
        by.y = "cs",
        all.x = TRUE,
        sort = FALSE
        )
        data.table::setorder(selected_sus, variable_index)
    }

    # ======================================================
    # 7. Add SuSiE model-level diagnostics
    # ======================================================
    selected_sus[, susie_converged := fitted$converged]
    selected_sus[, susie_niter     := fitted$niter]
    selected_sus[, susie_sigma2    := fitted$sigma2]
    selected_sus[, susie_log10BF1  := if (!is.null(fitted$lbf) && length(fitted$lbf)) {
      fitted$lbf[1] / log(10)
    } else {
      NA_real_
    }]
    if (!is.null(fitted$lbf)) {
        for (k in seq_along(fitted$lbf)) {
            data.table::set(
              selected_sus, j = paste0("susie_log10BF_L", k),
              value = fitted$lbf[k] / log(10)
            )
        }
    }

    # ======================================================
    # 8. Add SNP-level SuSiE signals: PIP, mu, mu2, alpha, lbf
    # ======================================================
    ## 8A — Add PIP
    model_index <- selected_sus$variable_index
    if (!is.null(fitted$pip)) {
        data.table::set(
          selected_sus, j = "global_pip", value = fitted$pip[model_index]
        )
    }
    ## 8B — Add lbf_variable (L x J)
    if (!is.null(fitted$lbf_variable)) {
        L <- nrow(fitted$lbf_variable)
        for (k in 1:L) {
          data.table::set(
            selected_sus, j = paste0("lbf_L", k),
            value = fitted$lbf_variable[k, model_index]
          )
        }
    }
    ## 8C — Add mu
    if (!is.null(fitted$mu)) {
        L <- nrow(fitted$mu)
        for (k in 1:L) {
          data.table::set(
            selected_sus, j = paste0("mu_L", k),
            value = fitted$mu[k, model_index]
          )
        }
    }
    ## 8D — Add mu2
    if (!is.null(fitted$mu2)) {
        L <- nrow(fitted$mu2)
        for (k in 1:L) {
          data.table::set(
            selected_sus, j = paste0("mu2_L", k),
            value = fitted$mu2[k, model_index]
          )
        }
    }
    ## 8E — Add alpha
    if (!is.null(fitted$alpha)) {
        L <- nrow(fitted$alpha)
        for (k in 1:L) {
          data.table::set(
            selected_sus, j = paste0("alpha_L", k),
            value = fitted$alpha[k, model_index]
          )
        }
    }
    # ======================================================
    # RETURN
    # ======================================================
    return(selected_sus)
}





plot_susie_loci_ld <- function(
    df,
    ld,
    sample_id      = "sample",
    outdir         = "susie_plots_ld",
    min_pip_label  = SUSIE_DEFAULTS$plot_min_pip_label,
    bg_size        = SUSIE_DEFAULTS$plot_background_size,
    bg_alpha       = SUSIE_DEFAULTS$plot_background_alpha,
    cs_size        = SUSIE_DEFAULTS$plot_cs_size,
    cs_alpha       = SUSIE_DEFAULTS$plot_cs_alpha,
    max_ld_dim     = SUSIE_DEFAULTS$plot_max_ld_dimension,
    plot_width     = SUSIE_DEFAULTS$plot_width,
    plot_height    = SUSIE_DEFAULTS$plot_height,
    plot_dpi       = SUSIE_DEFAULTS$plot_dpi
) {

    # ------------------------------------------------------------
    # PACKAGES
    # ------------------------------------------------------------
    suppressPackageStartupMessages(
      { require(data.table)
        require(ggplot2)
        require(ggrepel)
        require(patchwork)
        require(glue)
        require(RColorBrewer)
        require(cowplot)
        require(grid)
        require(grDevices)
        }
    )

      # ------------------------------------------------------------
      # INPUT CHECKS
      # ------------------------------------------------------------
      DT <- as.data.table(df)

      required_cols <- c(
          "BP","LP","SNP","variable_prob","cs",
          "locus_chr","locus_start","locus_end","variable_index"
      )
      missing <- setdiff(required_cols, names(DT))
      if (length(missing) > 0)
          stop("Missing columns: ", paste(missing, collapse = ", "))

      if (!is.matrix(ld)) stop("'ld' must be a matrix.")

      # ------------------------------------------------------------
      # NORMALIZE CS
      # ------------------------------------------------------------
      DT[, cs := trimws(as.character(cs))]
      suppressWarnings(DT[, cs := as.integer(cs)])
      DT[variable_prob > 0 & is.na(cs), cs := -2L]       # credible but not assigned
      DT[, cs_plot := ifelse(variable_prob > 0, cs, -1L)] # non-credible = -1
      DT[, cs_plot := factor(cs_plot, levels = sort(unique(cs_plot)))]

      # ------------------------------------------------------------
      # OUTPUT DIR
      # ------------------------------------------------------------
      if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

      # ------------------------------------------------------------
      # SINGLE LOCUS INFO
      # ------------------------------------------------------------
      chr   <- unique(DT$locus_chr)
      start <- unique(DT$locus_start)
      end   <- unique(DT$locus_end)
      locus_id <- glue("chr{chr}_{start}_{end}")

      locus_df    <- DT
      credible_df <- DT[variable_prob > 0]
      if (nrow(credible_df) == 0) {
          message("⚠️ No credible variants. Skipping.")
          return(NULL)
      }

      # ------------------------------------------------------------
      # COLORS FOR CS
      # ------------------------------------------------------------
      cs_levels <- levels(DT$cs_plot)
      pal <- rep("grey80", length(cs_levels)); names(pal) <- cs_levels
      if ("-2" %in% cs_levels) pal["-2"] <- "grey60"   # credible, unassigned

      real_cs <- cs_levels[as.numeric(as.character(cs_levels)) >= 0]
      if (length(real_cs) > 0)
          pal[real_cs] <- brewer.pal(8, "Dark2")[seq_along(real_cs)]

      # Shaded CS regions
      cs_regions <- credible_df[cs >= 0, .(
          xmin = min(BP),
          xmax = max(BP),
          cs_plot = factor(unique(cs), levels = cs_levels)
      ), by = cs]

      # ------------------------------------------------------------
      # LD SUBSET + OPTIONAL DOWNSAMPLING (FAST)
      # ------------------------------------------------------------
      idx <- locus_df$variable_index
      ld_sub <- ld[idx, idx, drop = FALSE]
      diag(ld_sub) <- 1

      n_ld <- nrow(ld_sub)
      if (n_ld > max_ld_dim) {
          # downsample indices to max_ld_dim in each dimension
          keep_idx <- unique(round(seq(1, n_ld, length.out = max_ld_dim)))
          ld_sub   <- ld_sub[keep_idx, keep_idx, drop = FALSE]
          #message("📉 Downsampled LD matrix from ", n_ld, " to ", nrow(ld_sub), " per dimension.")
      }

      # ------------------------------------------------------------
      # LD HEATMAP AS RASTER (NO LONG DATA.FRAME)
      # ------------------------------------------------------------
      # values in [0,1]
      z <- abs(ld_sub)
      z[z < 0] <- 0
      z[z > 1] <- 1

      # color palette and mapping
      n_col <- SUSIE_DEFAULTS$plot_ld_palette_colors
      pal_ld <- colorRampPalette(c("white", "red"))(n_col)
      z_idx <- floor(z * (n_col - 1L)) + 1L
      col_mat <- matrix(pal_ld[z_idx], nrow = nrow(z), ncol = ncol(z))

      ras <- as.raster(col_mat)

      ld_grob <- rasterGrob(
          image        = ras,
          x            = 0.5,
          y            = 0.5,
          width        = 1,
          height       = 1,
          interpolate  = FALSE
      )

      p_ld_panel <- patchwork::wrap_elements(ld_grob)

      # ------------------------------------------------------------
      # LD LEGEND ONLY (TINY DUMMY GGPLOT)
      # ------------------------------------------------------------
      df_leg <- data.frame(
          x    = 1,
          y    = seq(0, 1, length.out = SUSIE_DEFAULTS$plot_legend_gradient_points),
          fill = seq(0, 1, length.out = SUSIE_DEFAULTS$plot_legend_gradient_points)
      )

      p_ld_leg <- ggplot(df_leg, aes(x, y, fill = fill)) +
          geom_raster() +
          scale_fill_gradient(
              low    = "white",
              high   = "red",
              limits = c(0, 1),
              name   = "|r|"
          ) +
          theme_void() +
          theme(
              legend.position  = "right",
              legend.title     = element_text(size = SUSIE_DEFAULTS$plot_legend_title_size),
              legend.text      = element_text(size = SUSIE_DEFAULTS$plot_legend_text_size),
              legend.key.size  = unit(SUSIE_DEFAULTS$plot_legend_key_size_cm, "cm")
          )

      ld_leg <- cowplot::get_legend(p_ld_leg)
      ld_leg$widths  <- ld_leg$widths * SUSIE_DEFAULTS$plot_legend_scale
      ld_leg$heights <- ld_leg$heights * SUSIE_DEFAULTS$plot_legend_scale

      # ------------------------------------------------------------
      # CS LEGEND
      # ------------------------------------------------------------
      p_cs_leg_source <- ggplot(
          data.frame(cs_plot = factor(cs_levels, levels = cs_levels))
      ) +
          geom_point(
            aes(x = 1, y = cs_plot, color = cs_plot),
            size = SUSIE_DEFAULTS$plot_credible_set_legend_size
          ) +
          scale_color_manual(values = pal, name = "Credible Set") +
          theme_void() +
          theme(
              legend.position  = "right",
              legend.title     = element_text(size = SUSIE_DEFAULTS$plot_legend_title_size),
              legend.text      = element_text(size = SUSIE_DEFAULTS$plot_legend_text_size),
              legend.key.size  = unit(SUSIE_DEFAULTS$plot_legend_key_size_cm, "cm")
          )

      cs_leg <- cowplot::get_legend(p_cs_leg_source)
      cs_leg$widths  <- cs_leg$widths * SUSIE_DEFAULTS$plot_legend_scale
      cs_leg$heights <- cs_leg$heights * SUSIE_DEFAULTS$plot_legend_scale

      # ------------------------------------------------------------
      # COMBINE ONLY TWO LEGENDS (CS + LD)
      # ------------------------------------------------------------
      legend_column <- patchwork::wrap_elements(
          cowplot::plot_grid(
              cs_leg,
              ld_leg,
              ncol        = 1,
              rel_heights = c(1, 1)
          )
      )

      # ------------------------------------------------------------
      # LABEL VARIANTS
      # ------------------------------------------------------------
      label_df <- credible_df[
          cs >= 0 & variable_prob >= min_pip_label
      ][order(-variable_prob)][
        , head(.SD, SUSIE_DEFAULTS$plot_labels_per_credible_set), by = cs
      ]

      if (nrow(label_df) > 0)
          label_df[, cs_plot := factor(cs, levels = cs_levels)]

      n_snps <- nrow(locus_df)
      bg_size2 <- if (n_snps > SUSIE_DEFAULTS$plot_large_locus_threshold) {
        bg_size * SUSIE_DEFAULTS$plot_large_locus_scale
      } else if (n_snps > SUSIE_DEFAULTS$plot_medium_locus_threshold) {
        bg_size * SUSIE_DEFAULTS$plot_medium_locus_scale
      } else {
        bg_size
      }

      lift_amt <- SUSIE_DEFAULTS$plot_label_lift

      # ------------------------------------------------------------
      # LP PANEL
      # ------------------------------------------------------------
      p_lp <- ggplot() +
          geom_rect(
              data = cs_regions,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = cs_plot),
              alpha = cs_alpha
          ) +
          geom_point(
              data = locus_df,
              aes(x = BP, y = LP),
              alpha = bg_alpha, size = bg_size2, color = "grey75"
          ) +
          geom_point(
              data = credible_df,
              aes(x = BP, y = LP, color = cs_plot),
              size = cs_size
          ) +
          {
              if (nrow(label_df) > 0) geom_text_repel(
                  data = label_df,
                  aes(x = BP, y = LP, label = SNP, color = cs_plot),
                  size = SUSIE_DEFAULTS$plot_label_size,
                  segment.alpha = SUSIE_DEFAULTS$plot_label_segment_alpha,
                  max.overlaps = Inf
              )
          } +
          scale_color_manual(values = pal) +
          scale_fill_manual(values = pal) +
          theme_bw() +
          theme(legend.position = "none") +
          labs(
              y     = "LP (-log10 p)",
              title = glue("{sample_id}: chr{chr}:{start}-{end}")
          )

      # ------------------------------------------------------------
      # PIP PANEL
      # ------------------------------------------------------------
      p_pip <- ggplot() +
          geom_rect(
              data = cs_regions,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = cs_plot),
              alpha = cs_alpha
          ) +
          geom_point(
              data = locus_df,
              aes(x = BP, y = variable_prob),
              alpha = bg_alpha, size = bg_size2, color = "grey75"
          ) +
          geom_point(
              data = credible_df,
              aes(x = BP, y = variable_prob, color = cs_plot),
              size = cs_size
          ) +
          {
              if (nrow(label_df) > 0)
                  list(
                      geom_segment(
                          data = label_df,
                          aes(
                              x    = BP, y    = variable_prob,
                              xend = BP, yend = variable_prob + lift_amt,
                              color = cs_plot
                          ),
                          linewidth = SUSIE_DEFAULTS$plot_label_segment_width
                      ),
                      geom_text_repel(
                          data = transform(label_df, variable_prob = variable_prob + lift_amt),
                          aes(x = BP, y = variable_prob, label = SNP, color = cs_plot),
                          size = SUSIE_DEFAULTS$plot_label_size,
                          segment.alpha = SUSIE_DEFAULTS$plot_label_segment_alpha,
                          max.overlaps = Inf
                      )
                  )
          } +
          scale_color_manual(values = pal) +
          scale_fill_manual(values = pal) +
          theme_bw() +
          theme(legend.position = "none") +
          labs(x = "BP", y = "PIP")

      # ------------------------------------------------------------
      # COMBINE PANELS + LEGENDS
      # ------------------------------------------------------------
      left_panels <- (
          p_lp / p_pip / p_ld_panel
      ) + patchwork::plot_layout(heights = SUSIE_DEFAULTS$plot_panel_heights)

      combined <- left_panels | legend_column
      combined <- combined + patchwork::plot_layout(
        widths = SUSIE_DEFAULTS$plot_panel_widths
      )

      # ------------------------------------------------------------
      # SAVE PNG + PDF
      # ------------------------------------------------------------
      file_base <- glue("{outdir}/{sample_id}_{locus_id}_SUSIE_LD")

      ggsave(glue("{file_base}.png"), combined,
             width = plot_width, height = plot_height, dpi = plot_dpi)
      ggsave(glue("{file_base}.pdf"), combined,
             width = plot_width, height = plot_height)

      invisible(TRUE)
  }








generate_flames_files <- function(
  fitted,
  ld,
  snp_df,
  outfile,
  genomic_locus = NULL,
  genome_build = SUSIE_DEFAULTS$genome_build,
  target_coverage = SUSIE_DEFAULTS$credible_set_coverage,
  coverage_tolerance = SUSIE_DEFAULTS$credible_set_tolerance,
  resource_qc
) {
  suppressPackageStartupMessages({ require(data.table) })
  if (as.numeric(target_coverage) != SUSIE_DEFAULTS$credible_set_coverage) {
    stop("FLAMES export coverage differs from the resolved YAML configuration")
  }
  cs_list <- fitted$sets$cs

  if (length(cs_list) == 0) {
      #message("❌ No credible sets found — skipping FLAMES output...")
      return(invisible(NULL))
  }
  outdir <- dirname(outfile)
  prefix <- tools::file_path_sans_ext(outfile)

  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
    #message("📁 Created folder: ", outdir)
  }

  snp_df <- data.table::as.data.table(snp_df)
  manifest_rows <- list()
  index_rows <- list()

  for (cs_pos in seq_along(cs_list)) {
    cs_name <- names(cs_list)[cs_pos]
    if (is.null(cs_name) || is.na(cs_name) || !nzchar(cs_name)) cs_name <- paste0("L", cs_pos)
    cs <- cs_list[[cs_pos]]
    if (is.null(cs) || length(cs) == 0) next

    # FLAMES explicitly requires the fine-mapping PIP for each SNP. SuSiE
    # credible-set membership is component-specific (alpha), whereas the SNP
    # probability consumed by FLAMES is the model-wide PIP. Do not substitute
    # alpha for prob1; record both definitions in the manifest below.
    pip <- fitted$pip[cs]
    if (any(!is.finite(pip)) || any(pip < 0 | pip > 1)) {
      stop("Invalid SuSiE PIP encountered during FLAMES export")
    }

    # Sort by PIP
    ord <- order(pip, decreasing = TRUE)
    cs  <- cs[ord]
    pip <- pip[ord]

    chr_value <- as.character(snp_df$CHR[cs])
    chr_value <- sub("^chr", "", chr_value, ignore.case = TRUE)
    chr_value[toupper(chr_value) == "X"] <- "23"
    chr_value[toupper(chr_value) == "Y"] <- "24"
    snps <- paste0(
      chr_value, ":", snp_df$BP[cs], ":",
      toupper(snp_df$ALT[cs]), "_", toupper(snp_df$REF[cs])
    )

    # --------------------------
    # LD Statistics
    # --------------------------
    if (length(cs) == 1) {
      min_ld    <- NA
      mean_ld   <- NA
      median_ld <- NA
    } else {
      ld_sub  <- abs(ld[cs, cs, drop = FALSE])
      ld_vals <- ld_sub[upper.tri(ld_sub)]
      min_ld    <- min(ld_vals)
      mean_ld   <- mean(ld_vals)
      median_ld <- median(ld_vals)
    }

    cs_filename <- sprintf("%s_CS_%s.txt", prefix, cs_name)
    flames_df <- data.table::data.table(
      index = seq_along(cs),
      cred1 = snps,
      prob1 = signif(pip, 12)
    )
    data.table::fwrite(flames_df, cs_filename, sep = " ", quote = FALSE)

    component_index <- if (!is.null(fitted$sets$cs_index) &&
                           length(fitted$sets$cs_index) >= cs_pos) {
      as.integer(fitted$sets$cs_index[[cs_pos]])
    } else {
      suppressWarnings(as.integer(sub("^L", "", cs_name)))
    }
    achieved_coverage <- NA_real_
    component_log10bf <- NA_real_
    if (is.finite(component_index) && component_index >= 1L) {
      if (!is.null(fitted$alpha) && component_index <= nrow(fitted$alpha)) {
        achieved_coverage <- sum(fitted$alpha[component_index, cs], na.rm = TRUE)
      }
      if (!is.null(fitted$lbf) && component_index <= length(fitted$lbf)) {
        component_log10bf <- fitted$lbf[component_index] / log(10)
      }
    }
    if (is.finite(achieved_coverage) &&
        achieved_coverage < target_coverage - coverage_tolerance) {
      stop("SuSiE credible set did not achieve the configured component coverage")
    }
    global_pip_sum <- sum(pip)
    if (!is.finite(global_pip_sum) ||
        global_pip_sum < target_coverage - coverage_tolerance) {
      stop("SuSiE FLAMES probabilities did not achieve the configured PIP mass")
    }

    locus_value <- if (!is.null(genomic_locus) && length(genomic_locus) &&
                       !is.na(genomic_locus[[1L]]) && nzchar(as.character(genomic_locus[[1L]]))) {
      as.character(genomic_locus[[1L]])
    } else {
      basename(prefix)
    }
    index_rows[[length(index_rows) + 1L]] <- data.table::data.table(
      Filename = basename(cs_filename),
      GenomicLocus = locus_value
    )
    manifest_rows[[length(manifest_rows) + 1L]] <- data.table::data.table(
      Filename = basename(cs_filename),
      GenomicLocus = locus_value,
      genome_build = genome_build,
      credible_set = cs_name,
      component_index = component_index,
      target_coverage = target_coverage,
      achieved_component_coverage = achieved_coverage,
      global_pip_sum = global_pip_sum,
      credible_set_membership_definition =
        "susie_component_alpha_configured_coverage",
      probability_definition = "susie_model_wide_pip",
      n_variants = length(cs),
      component_log10bf = component_log10bf,
      min_abs_ld = min_ld,
      mean_abs_ld = mean_ld,
      median_abs_ld = median_ld,
      input_variant_count = resource_qc$input_variant_count,
      maximum_variants_per_locus = resource_qc$maximum_variants_per_locus,
      dense_ld_matrix_gb = resource_qc$dense_ld_matrix_gb,
      ld_peak_matrix_multiplier = resource_qc$ld_peak_matrix_multiplier,
      estimated_peak_ld_memory_gb = resource_qc$estimated_peak_ld_memory_gb,
      reserved_worker_memory_gb = resource_qc$reserved_worker_memory_gb,
      resource_warning_reason = resource_qc$resource_warning_reason,
      resource_failure_reason = resource_qc$resource_failure_reason
    )
  }

  if (length(index_rows)) {
    index_file <- file.path(outdir, "indexfile_rows.tsv")
    data.table::fwrite(
      data.table::rbindlist(index_rows), index_file, sep = "\t",
      append = file.exists(index_file), col.names = !file.exists(index_file)
    )
    manifest_file <- sprintf("%s_FLAMES_manifest.tsv", prefix)
    data.table::fwrite(data.table::rbindlist(manifest_rows), manifest_file, sep = "\t")
  }

  invisible(index_rows)
}


aggregate_and_write_results <- function(
  results_all,
  sample_id,
  analysis_folder,
  verbose = TRUE,
  target_coverage = SUSIE_DEFAULTS$credible_set_coverage,
  recovery_audit_file
) {
  if (as.numeric(target_coverage) != SUSIE_DEFAULTS$credible_set_coverage) {
    stop("SuSiE aggregate coverage differs from the resolved YAML configuration")
  }

  # ------------------------------------------------------------------
  # ROBUST msg() wrapper: use parent's msg() if available,
  # otherwise use a simple internal message() fallback
  # ------------------------------------------------------------------
  msg <- NULL
  if (exists("msg", envir = parent.frame(), inherits = FALSE)) {
    msg <- get("msg", envir = parent.frame())
  } else {
    msg <- function(verbose, ...) {
      if (isTRUE(verbose)) message(...)
    }
  }
  # ------------------------------------------------------------------

  combined_vars <- list()
  combined_cs   <- list()
  qc_rows       <- list()
  audit_rows    <- list()

  for (res in results_all) {
    locus <- res$locus

    is_success <- res$status %in% c("OK", "RECOVERED")
    locus_name <- if (!is.null(locus$genomic_locus)) as.character(locus$genomic_locus) else NA_character_
    ld_lambda <- if (!is.null(res$ld_qc_lambda)) as.numeric(res$ld_qc_lambda) else NA_real_
    ld_min_eigenvalue <- if (!is.null(res$ld_min_eigenvalue)) {
      as.numeric(res$ld_min_eigenvalue)
    } else {
      NA_real_
    }
    sample_size_qc <- if (!is.null(res$sample_size_qc)) {
      res$sample_size_qc
    } else {
      list(
        nef_min = NA_real_, nef_max = NA_real_, nef_median = NA_real_,
        nef_relative_range = NA_real_, selected_nef = NA_real_,
        nef_policy = NA_character_, nef_summary_statistic = NA_character_,
        nef_warning_threshold = NA_real_, warning_reason = NA_character_
      )
    }
    resource_qc <- if (!is.null(res$resource_qc)) {
      res$resource_qc
    } else {
      list(
        input_variant_count = NA_integer_,
        maximum_variants_per_locus = NA_integer_,
        dense_ld_matrix_gb = NA_real_,
        ld_peak_matrix_multiplier = NA_real_,
        estimated_peak_ld_memory_gb = NA_real_,
        reserved_worker_memory_gb = NA_real_,
        resource_warning_reason = NA_character_,
        resource_failure_reason = NA_character_
      )
    }
    fit_input_qc <- if (!is.null(res$fit_input_qc)) {
      res$fit_input_qc
    } else {
      summarize_susie_fit_input(NULL)
    }
    if (!is.null(res$audit) && length(res$audit)) {
      audit_rows <- c(audit_rows, res$audit)
    }

    if (is_success) {
      if (!is.null(res$vars)) combined_vars[[length(combined_vars) + 1]] <- res$vars
      if (!is.null(res$cs))   combined_cs[[length(combined_cs) + 1]]   <- res$cs

      qc_rows[[length(qc_rows) + 1]] <- data.frame(
        locus_chr   = locus$chr,
        locus_start = locus$start,
        locus_end   = locus$end,
        genomic_locus = locus_name,
        stage       = if (res$status == "RECOVERED") "recovered" else "final",
        converged   = TRUE,
        note        = if (!is.null(res$reason)) res$reason else res$status,
        recovery_method = if (res$status == "RECOVERED") {
          if (!is.null(res$reason)) res$reason else "recovered"
        } else NA_character_,
        ld_z_mismatch_lambda = ld_lambda,
        ld_min_eigenvalue = ld_min_eigenvalue,
        nef_min = sample_size_qc$nef_min,
        nef_max = sample_size_qc$nef_max,
        nef_median = sample_size_qc$nef_median,
        nef_relative_range = sample_size_qc$nef_relative_range,
        selected_nef = sample_size_qc$selected_nef,
        nef_policy = sample_size_qc$nef_policy,
        nef_summary_statistic = sample_size_qc$nef_summary_statistic,
        nef_warning_threshold = sample_size_qc$nef_warning_threshold,
        warning_reason = sample_size_qc$warning_reason,
        input_variant_count = resource_qc$input_variant_count,
        maximum_variants_per_locus = resource_qc$maximum_variants_per_locus,
        dense_ld_matrix_gb = resource_qc$dense_ld_matrix_gb,
        ld_peak_matrix_multiplier = resource_qc$ld_peak_matrix_multiplier,
        estimated_peak_ld_memory_gb = resource_qc$estimated_peak_ld_memory_gb,
        reserved_worker_memory_gb = resource_qc$reserved_worker_memory_gb,
        resource_warning_reason = resource_qc$resource_warning_reason,
        resource_failure_reason = resource_qc$resource_failure_reason,
        fit_input_matrix_state = fit_input_qc$fit_input_matrix_state,
        fit_input_serializations = fit_input_qc$fit_input_serializations,
        fit_child_processes = fit_input_qc$fit_child_processes,
        fit_input_reuses = fit_input_qc$fit_input_reuses,
        fit_input_file_size_bytes = fit_input_qc$fit_input_file_size_bytes,
        fit_input_serialization_seconds =
          fit_input_qc$fit_input_serialization_seconds,
        failure_reason = NA_character_,
        failure_detail = NA_character_,
        target_coverage = target_coverage,
        stringsAsFactors = FALSE
      )

    } else {
      qc_rows[[length(qc_rows) + 1]] <- data.frame(
        locus_chr   = locus$chr,
        locus_start = locus$start,
        locus_end   = locus$end,
        genomic_locus = locus_name,
        stage       = res$status,
        converged   = FALSE,
        note        = if (!is.null(res$reason)) res$reason else res$status,
        recovery_method = NA_character_,
        ld_z_mismatch_lambda = ld_lambda,
        ld_min_eigenvalue = ld_min_eigenvalue,
        nef_min = sample_size_qc$nef_min,
        nef_max = sample_size_qc$nef_max,
        nef_median = sample_size_qc$nef_median,
        nef_relative_range = sample_size_qc$nef_relative_range,
        selected_nef = sample_size_qc$selected_nef,
        nef_policy = sample_size_qc$nef_policy,
        nef_summary_statistic = sample_size_qc$nef_summary_statistic,
        nef_warning_threshold = sample_size_qc$nef_warning_threshold,
        warning_reason = sample_size_qc$warning_reason,
        input_variant_count = resource_qc$input_variant_count,
        maximum_variants_per_locus = resource_qc$maximum_variants_per_locus,
        dense_ld_matrix_gb = resource_qc$dense_ld_matrix_gb,
        ld_peak_matrix_multiplier = resource_qc$ld_peak_matrix_multiplier,
        estimated_peak_ld_memory_gb = resource_qc$estimated_peak_ld_memory_gb,
        reserved_worker_memory_gb = resource_qc$reserved_worker_memory_gb,
        resource_warning_reason = resource_qc$resource_warning_reason,
        resource_failure_reason = resource_qc$resource_failure_reason,
        fit_input_matrix_state = fit_input_qc$fit_input_matrix_state,
        fit_input_serializations = fit_input_qc$fit_input_serializations,
        fit_child_processes = fit_input_qc$fit_child_processes,
        fit_input_reuses = fit_input_qc$fit_input_reuses,
        fit_input_file_size_bytes = fit_input_qc$fit_input_file_size_bytes,
        fit_input_serialization_seconds =
          fit_input_qc$fit_input_serialization_seconds,
        failure_reason = if (!is.null(res$reason)) res$reason else res$status,
        failure_detail = if (!is.null(res$failure_detail)) {
          res$failure_detail
        } else {
          NA_character_
        },
        target_coverage = target_coverage,
        stringsAsFactors = FALSE
      )
    }
  }

  df_vars <- if (length(combined_vars)) data.table::rbindlist(combined_vars, fill = TRUE) else data.table::data.table()
  df_cs   <- if (length(combined_cs))   data.table::rbindlist(combined_cs,   fill = TRUE) else data.table::data.table()

  out_combined <- glue::glue("{analysis_folder}/{sample_id}_SUSIE_combined_results.csv")
  out_credible <- glue::glue("{analysis_folder}/{sample_id}_SUSIE_combined_credibleset.csv")
  qc_file      <- glue::glue("{analysis_folder}/{sample_id}_SuSiE_QC_summary.tsv")

  if (nrow(df_vars) > 0) {
    data.table::fwrite(df_vars, paste0(out_combined, ".gz"), compress = "gzip")
  } else {
    msg(verbose, glue::glue("⚠ No SuSiE results to write for {sample_id}; combined table is empty."))
  }

  if (length(qc_rows)) { 
    qc_df <- do.call(rbind, qc_rows)
    data.table::fwrite(qc_df, qc_file, sep = "\t")
  } else {
    qc_df <- data.frame(
      locus_chr = character(), locus_start = integer(), locus_end = integer(),
      genomic_locus = character(), stage = character(), converged = logical(),
      note = character(), recovery_method = character(),
      ld_z_mismatch_lambda = numeric(), ld_min_eigenvalue = numeric(),
      nef_min = numeric(), nef_max = numeric(), nef_median = numeric(),
      nef_relative_range = numeric(), selected_nef = numeric(),
      nef_policy = character(), nef_summary_statistic = character(),
      nef_warning_threshold = numeric(), warning_reason = character(),
      input_variant_count = integer(), maximum_variants_per_locus = integer(),
      dense_ld_matrix_gb = numeric(), ld_peak_matrix_multiplier = numeric(),
      estimated_peak_ld_memory_gb = numeric(),
      reserved_worker_memory_gb = numeric(),
      resource_warning_reason = character(),
      resource_failure_reason = character(),
      fit_input_matrix_state = character(),
      fit_input_serializations = integer(), fit_child_processes = integer(),
      fit_input_reuses = integer(), fit_input_file_size_bytes = numeric(),
      fit_input_serialization_seconds = numeric(),
      failure_reason = character(), failure_detail = character(),
      target_coverage = numeric()
    )
    data.table::fwrite(qc_df, qc_file, sep = "\t")
  }

  if (nrow(df_cs) > 0) {
    data.table::fwrite(df_cs, out_credible) 
  } else {
    msg(verbose, glue::glue("⚠ No SuSiE credible set to write for {sample_id}; combined CS is empty."))
  }
  
  data.table::setDT(qc_df)
  failed_df <- qc_df[converged == FALSE]
  fail_file <- glue::glue("{analysis_folder}/{sample_id}_SuSiE_failed_loci.tsv")
  data.table::fwrite(failed_df, fail_file, sep = "\t")
  audit_df <- if (length(audit_rows)) {
    data.table::rbindlist(audit_rows, fill = TRUE)
  } else {
    data.table::data.table(
      timestamp_utc = character(), locus_chr = character(),
      locus_start = integer(), locus_end = integer(), genomic_locus = character(),
      sequence = integer(), stage = character(), action = character(),
      status = character(), reason = character(), ld_status = character(),
      repairable = logical(), L = integer(), max_iter = integer(),
      timeout_seconds = numeric(), ld_z_mismatch_lambda = numeric(),
      ld_min_eigenvalue = numeric(), details = character()
    )
  }
  data.table::fwrite(audit_df, recovery_audit_file, sep = "\t")
  if (nrow(failed_df) > 0) {
    failed_loci <- nrow(failed_df)
    msg(verbose, glue::glue(
      "⚠ {failed_loci} loci failed or were skipped; see {fail_file} for details."
    ))
  } else {
    msg(verbose, "✅ No failed or MHC-skipped loci found.")
  }

  msg(verbose, glue::glue(
    "SuSiE output finalization completed.\n",
    "  • Combined results: {out_combined}.gz\n",
    "  • QC summary:      {qc_file}\n",
    "  • Recovery audit:  {recovery_audit_file}"
  ))

  return(df_vars)
}












# plot_susie_loci_advanced <- function(
#     df,
#     sample_id     = "sample",
#     outdir        = "susie_plots_advanced",
#     n_cores       = 4,
#     min_pip_label = 0.1
# ) {
#     require(data.table)
#     require(ggplot2)
#     require(ggrepel)
#     require(parallel)
#     require(patchwork)
#     require(glue)
#     require(RColorBrewer)

#     DT <- as.data.table(df)

#     # REQUIRED COLUMNS
#     required_cols <- c("BP","LP","SNP","variable_prob","cs",
#                        "locus_chr","locus_start","locus_end")
#     missing <- setdiff(required_cols, names(DT))
#     if (length(missing) > 0)
#         stop("Missing columns: ", paste(missing, collapse=", "))

#     # CS NORMALIZATION
#     DT[, cs := trimws(as.character(cs))]
#     suppressWarnings(DT[, cs := as.integer(cs)])
#     DT[variable_prob > 0 & is.na(cs), cs := -2]   # malformed credible sets
#     DT[, cs_plot := ifelse(variable_prob > 0, cs, -1)]
#     DT[, cs_plot := factor(cs_plot, levels = sort(unique(cs_plot)))]

#     # OUTPUT FOLDER
#     if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

#     # LOCUS ID
#     DT[, locus_id := sprintf("chr%s_%d_%d", locus_chr, locus_start, locus_end)]
#     loci <- unique(DT[,.(locus_id,locus_chr,locus_start,locus_end)])

#     # COLOR PALETTE
#     cs_levels <- levels(DT$cs_plot)
#     pal <- rep("grey80", length(cs_levels))
#     names(pal) <- cs_levels
#     if ("-2" %in% cs_levels) pal["-2"] <- "grey60"
#     cs_real <- cs_levels[as.numeric(as.character(cs_levels)) >= 0]
#     pal[cs_real] <- brewer.pal(8, "Dark2")[seq_along(cs_real)]

#     # SHADED CS REGIONS
#     get_cs_regions <- function(df) {
#         df[cs >= 0, .(
#             xmin = min(BP),
#             xmax = max(BP),
#             cs_plot = factor(unique(cs), levels = cs_levels)
#         ), by = cs]
#     }

#     # ONE LOCUS
#     plot_one_locus <- function(i) {

#         this_locus <- loci[i]
#         locus_df   <- DT[locus_id == this_locus$locus_id]
#         credible_df <- locus_df[variable_prob > 0]
#         if (nrow(credible_df) == 0L) return(NULL)

#         # LABEL TOP 5 PER CS (cs >= 0)
#         label_df <- credible_df[cs >= 0][order(-variable_prob)][, head(.SD,5), by=cs]
#         label_df[, cs_plot := factor(cs, levels = cs_levels)]

#         cs_regions <- get_cs_regions(credible_df)

#         n_snps <- nrow(locus_df)
#         bg_size <- if (n_snps > 5000) 0.005 else if (n_snps > 2000) 0.008 else 0.01

#         # LP PANEL
#         p_lp <- ggplot() +
#             geom_rect(
#                 data = cs_regions,
#                 aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = cs_plot),
#                 alpha = 0.15, inherit.aes = FALSE
#             ) +
#             geom_point(
#                 data = locus_df,
#                 aes(x = BP, y = LP),
#                 alpha = 0.18, size = bg_size, color = "grey80"
#             ) +
#             geom_point(
#                 data = credible_df,
#                 aes(x = BP, y = LP, color = cs_plot),
#                 size = 1.5, alpha = 0.95
#             ) +
#             # ---------------------------
#             # NO LABEL BACKGROUND → geom_text_repel()
#             # ---------------------------
#             geom_text_repel(
#                 data = label_df,
#                 aes(x = BP, y = LP, label = SNP, color = cs_plot),  # Option B: replace "color = cs_plot" with "color = 'black'"
#                 size = 3,
#                 max.overlaps = Inf,
#                 box.padding = 0.4,
#                 point.padding = 0.3,
#                 segment.alpha = 0.7,
#                 segment.color = "black",
#                 label.size = 0
#             ) +
#             scale_color_manual(values = pal) +
#             scale_fill_manual(values = pal) +
#             theme_bw(base_size = 13) +
#             labs(
#                 x = NULL,
#                 y = "LP",
#                 fill  = "Credible Set",
#                 color = "Credible Set",
#                 title = glue("{sample_id} – SuSiE Fine-Mapping\nchr{this_locus$locus_chr}:{this_locus$locus_start}-{this_locus$locus_end}")
#             )

#         # PIP PANEL
#         lift_amt <- 0.03

#         p_pip <- ggplot() +
#             geom_rect(
#                 data = cs_regions,
#                 aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = cs_plot),
#                 alpha = 0.15, inherit.aes = FALSE
#             ) +
#             geom_point(
#                 data = locus_df,
#                 aes(x = BP, y = variable_prob),
#                 size=bg_size, alpha=0.15, color="grey80"
#             ) +
#             geom_point(
#                 data = credible_df,
#                 aes(x = BP, y = variable_prob, color = cs_plot),
#                 size = 1.3, alpha = 0.95
#             ) +
#             geom_segment(
#                 data = label_df,
#                 aes(
#                     x = BP,
#                     y = variable_prob,
#                     xend = BP,
#                     yend = variable_prob + lift_amt,
#                     color = cs_plot
#                 ),
#                 linewidth = 0.35,
#                 alpha = 0.7
#             ) +
#             geom_text_repel(
#                 data = transform(label_df, variable_prob = variable_prob + lift_amt),
#                 aes(x = BP, y = variable_prob, label = SNP, color = cs_plot),  # Option B: set color="black"
#                 size = 3,
#                 max.overlaps = Inf,
#                 segment.alpha = 0.7,
#                 segment.color = "black",
#                 label.size = 0,
#                 point.padding = 0.3,
#                 box.padding = 0.4
#             ) +
#             scale_color_manual(values = pal) +
#             scale_fill_manual(values = pal) +
#             theme_bw(base_size = 13) +
#             theme(legend.position = "none") +
#             labs(x = "Genomic Position (BP)", y = "PIP")

#         combined <- p_lp / p_pip + patchwork::plot_layout(heights = c(2,1))

#         outfile <- sprintf(
#             "%s/%s_chr%s_%d_%d_SUSIE_advanced.png",
#             outdir, sample_id,
#             this_locus$locus_chr,
#             this_locus$locus_start,
#             this_locus$locus_end
#         )
#         ggsave(outfile, combined, width=10, height=9, dpi=300)
#         return(outfile)
#     }

#     # PARALLELIZATION (macOS-safe)
#     cl <- parallel::makeCluster(n_cores, type="PSOCK")
#     parallel::clusterExport(cl,
#         varlist=c("DT","loci","plot_one_locus","pal","cs_levels","sample_id","outdir"),
#         envir=environment()
#     )
#     parallel::clusterEvalQ(cl, {
#         library(data.table); library(ggplot2)
#         library(ggrepel); library(patchwork)
#         library(glue); library(RColorBrewer)
#     })
#     results <- parallel::parLapply(cl, seq_len(nrow(loci)), plot_one_locus)
#     parallel::stopCluster(cl)

#     return(results)
# }
