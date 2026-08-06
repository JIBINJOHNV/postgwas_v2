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

  })

generate_flames_files <- function(
  fitted,
  ld,
  snp_df,
  outfile 
) {
  require(data.table)

  cs_list <- fitted$sets$cs

  if (length(cs_list) == 0) {
    stop("❌ No credible sets found.")
  }

  outdir <- dirname(outfile)
  prefix <- tools::file_path_sans_ext(outfile)

  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
    message("📁 Created folder: ", outdir)
  }

  # Compute log10BF from SuSiE
  log10bf <- fitted$lbf[1] / log(10)

  # -------------------------------------
  # Loop through all credible sets
  # -------------------------------------

  for (cs_name in names(cs_list)) {

    cs <- cs_list[[cs_name]]
    if (is.null(cs) || length(cs) == 0) next

    pip <- fitted$pip[cs]

    # Sort by PIP
    ord <- order(pip, decreasing = TRUE)
    cs  <- cs[ord]
    pip <- pip[ord]

    # SNP names CHR:BP:REF_ALT
    snps <- snp_df[cs, paste0(CHR, ":", BP, ":", ALT, "_", REF)]

    # --------------------------
    # LD Statistics
    # --------------------------
    if (length(cs) == 1) {
      min_ld    <- NA
      mean_ld   <- NA
      median_ld <- NA
    } else {
      ld_sub  <- abs(ld[cs, cs])
      ld_vals <- ld_sub[upper.tri(ld_sub)]
      min_ld    <- min(ld_vals)
      mean_ld   <- mean(ld_vals)
      median_ld <- median(ld_vals)
    }

    # --------------------------
    cs_filename <- sprintf("%s_CS_%s.txt", prefix, cs_name)
    con <- file(cs_filename, "wt")

    # --------------------------
    # Write minimal FLAMES block
    # --------------------------
    writeLines(c(
      paste0("# Post-Pr(# of causal SNPs is ", length(cs), ") = 1"),
      paste0("#log10bf ", round(log10bf, 6), " NA"),
      paste0("#min(|ld|) ", round(min_ld, 6), " NA"),
      paste0("#mean(|ld|) ", round(mean_ld, 6), " NA"),
      paste0("#median(|ld|) ", round(median_ld, 6), " NA"),
      "index cred1 prob1"
    ), con)

    # SNP rows
    for (i in seq_along(cs)) {
      writeLines(
        paste(i, snps[i], round(pip[i], 6)),
        con
      )
    }

    close(con)
  }
}


plot_susie_loci_advanced <- function(
    df,
    sample_id     = "sample",
    outdir        = "susie_plots_advanced",
    n_cores       = 4,
    min_pip_label = 0.1
) {
    require(data.table)
    require(ggplot2)
    require(ggrepel)
    require(parallel)
    require(patchwork)
    require(glue)
    require(RColorBrewer)

    DT <- as.data.table(df)

    # REQUIRED COLUMNS
    required_cols <- c("BP","LP","SNP","variable_prob","cs",
                       "locus_chr","locus_start","locus_end")
    missing <- setdiff(required_cols, names(DT))
    if (length(missing) > 0)
        stop("Missing columns: ", paste(missing, collapse=", "))

    # CS NORMALIZATION
    DT[, cs := trimws(as.character(cs))]
    suppressWarnings(DT[, cs := as.integer(cs)])
    DT[variable_prob > 0 & is.na(cs), cs := -2]   # malformed credible sets
    DT[, cs_plot := ifelse(variable_prob > 0, cs, -1)]
    DT[, cs_plot := factor(cs_plot, levels = sort(unique(cs_plot)))]

    # OUTPUT FOLDER
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

    # LOCUS ID
    DT[, locus_id := sprintf("chr%s_%d_%d", locus_chr, locus_start, locus_end)]
    loci <- unique(DT[,.(locus_id,locus_chr,locus_start,locus_end)])

    # COLOR PALETTE
    cs_levels <- levels(DT$cs_plot)
    pal <- rep("grey80", length(cs_levels))
    names(pal) <- cs_levels
    if ("-2" %in% cs_levels) pal["-2"] <- "grey60"
    cs_real <- cs_levels[as.numeric(as.character(cs_levels)) >= 0]
    pal[cs_real] <- brewer.pal(8, "Dark2")[seq_along(cs_real)]

    # SHADED CS REGIONS
    get_cs_regions <- function(df) {
        df[cs >= 0, .(
            xmin = min(BP),
            xmax = max(BP),
            cs_plot = factor(unique(cs), levels = cs_levels)
        ), by = cs]
    }

    # ONE LOCUS
    plot_one_locus <- function(i) {

        this_locus <- loci[i]
        locus_df   <- DT[locus_id == this_locus$locus_id]
        credible_df <- locus_df[variable_prob > 0]
        if (nrow(credible_df) == 0L) return(NULL)

        # LABEL TOP 5 PER CS (cs >= 0)
        label_df <- credible_df[cs >= 0][order(-variable_prob)][, head(.SD,5), by=cs]
        label_df[, cs_plot := factor(cs, levels = cs_levels)]

        cs_regions <- get_cs_regions(credible_df)

        n_snps <- nrow(locus_df)
        bg_size <- if (n_snps > 5000) 0.005 else if (n_snps > 2000) 0.008 else 0.01

        # LP PANEL
        p_lp <- ggplot() +
            geom_rect(
                data = cs_regions,
                aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = cs_plot),
                alpha = 0.15, inherit.aes = FALSE
            ) +
            geom_point(
                data = locus_df,
                aes(x = BP, y = LP),
                alpha = 0.18, size = bg_size, color = "grey80"
            ) +
            geom_point(
                data = credible_df,
                aes(x = BP, y = LP, color = cs_plot),
                size = 1.5, alpha = 0.95
            ) +
            # ---------------------------
            # NO LABEL BACKGROUND → geom_text_repel()
            # ---------------------------
            geom_text_repel(
                data = label_df,
                aes(x = BP, y = LP, label = SNP, color = cs_plot),  # Option B: replace "color = cs_plot" with "color = 'black'"
                size = 3,
                max.overlaps = Inf,
                box.padding = 0.4,
                point.padding = 0.3,
                segment.alpha = 0.7,
                segment.color = "black",
                label.size = 0
            ) +
            scale_color_manual(values = pal) +
            scale_fill_manual(values = pal) +
            theme_bw(base_size = 13) +
            labs(
                x = NULL,
                y = "LP",
                fill  = "Credible Set",
                color = "Credible Set",
                title = glue("{sample_id} – SuSiE Fine-Mapping\nchr{this_locus$locus_chr}:{this_locus$locus_start}-{this_locus$locus_end}")
            )

        # PIP PANEL
        lift_amt <- 0.03

        p_pip <- ggplot() +
            geom_rect(
                data = cs_regions,
                aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = cs_plot),
                alpha = 0.15, inherit.aes = FALSE
            ) +
            geom_point(
                data = locus_df,
                aes(x = BP, y = variable_prob),
                size=bg_size, alpha=0.15, color="grey80"
            ) +
            geom_point(
                data = credible_df,
                aes(x = BP, y = variable_prob, color = cs_plot),
                size = 1.3, alpha = 0.95
            ) +
            geom_segment(
                data = label_df,
                aes(
                    x = BP,
                    y = variable_prob,
                    xend = BP,
                    yend = variable_prob + lift_amt,
                    color = cs_plot
                ),
                linewidth = 0.35,
                alpha = 0.7
            ) +
            geom_text_repel(
                data = transform(label_df, variable_prob = variable_prob + lift_amt),
                aes(x = BP, y = variable_prob, label = SNP, color = cs_plot),  # Option B: set color="black"
                size = 3,
                max.overlaps = Inf,
                segment.alpha = 0.7,
                segment.color = "black",
                label.size = 0,
                point.padding = 0.3,
                box.padding = 0.4
            ) +
            scale_color_manual(values = pal) +
            scale_fill_manual(values = pal) +
            theme_bw(base_size = 13) +
            theme(legend.position = "none") +
            labs(x = "Genomic Position (BP)", y = "PIP")

        combined <- p_lp / p_pip + patchwork::plot_layout(heights = c(2,1))

        outfile <- sprintf(
            "%s/%s_chr%s_%d_%d_SUSIE_advanced.png",
            outdir, sample_id,
            this_locus$locus_chr,
            this_locus$locus_start,
            this_locus$locus_end
        )
        ggsave(outfile, combined, width=10, height=9, dpi=300)
        return(outfile)
    }

    # PARALLELIZATION (macOS-safe)
    cl <- parallel::makeCluster(n_cores, type="PSOCK")
    parallel::clusterExport(cl,
        varlist=c("DT","loci","plot_one_locus","pal","cs_levels","sample_id","outdir"),
        envir=environment()
    )
    parallel::clusterEvalQ(cl, {
        library(data.table); library(ggplot2)
        library(ggrepel); library(patchwork)
        library(glue); library(RColorBrewer)
    })
    results <- parallel::parLapply(cl, seq_len(nrow(loci)), plot_one_locus)
    parallel::stopCluster(cl)

    return(results)
}


aggregate_and_write_results <- function(
  results_all,
  sample_id,
  analysis_folder,
  verbose = TRUE
) {

  combined_vars <- list()
  combined_cs   <- list()
  qc_rows       <- list()

  for (res in results_all) {
    locus <- res$locus

    if (res$status == "OK") {
      if (!is.null(res$vars)) combined_vars[[length(combined_vars) + 1]] <- res$vars
      if (!is.null(res$cs))   combined_cs[[length(combined_cs) + 1]]   <- res$cs

      qc_rows[[length(qc_rows) + 1]] <- data.frame(
        locus_chr   = locus$chr,
        locus_start = locus$start,
        locus_end   = locus$end,
        stage       = "final",
        converged   = TRUE,
        note        = "OK",
        stringsAsFactors = FALSE
      )

    } else {
      qc_rows[[length(qc_rows) + 1]] <- data.frame(
        locus_chr   = locus$chr,
        locus_start = locus$start,
        locus_end   = locus$end,
        stage       = res$status,
        converged   = FALSE,
        note        = if (!is.null(res$reason)) res$reason else res$status,
        stringsAsFactors = FALSE
      )
    }
  }

  df_vars <- if (length(combined_vars)) data.table::rbindlist(combined_vars, fill = TRUE) else data.table::data.table()
  df_cs   <- if (length(combined_cs))   data.table::rbindlist(combined_cs,   fill = TRUE) else data.table::data.table()

  out_combined <- glue::glue("{analysis_folder}/{sample_id}_SUSIE_combined_results.csv")
  out_credible <- glue::glue("{analysis_folder}/{sample_id}_SUSIE_combined_credibleset.csv")

  qc_file      <- glue::glue("{analysis_folder}/{sample_id}_SuSiE_QC_summary.tsv")

  if (nrow(df_vars) > 0) { data.table::fwrite(df_vars, out_combined) } else {
    msg(verbose,
        glue::glue("⚠ No SuSiE results to write for {sample_id}; combined table is empty."))
  }

  if (length(qc_rows)) { qc_df <- do.call(rbind, qc_rows)
    data.table::fwrite(qc_df, qc_file, sep = "\t")
  } else {
    qc_df <- data.frame()
  }

  if (nrow(df_cs) > 0) { data.table::fwrite(df_cs, out_credible) } else {
    msg(verbose,
        glue::glue("⚠ No SuSiE credible set  to write for {sample_id}; combined table is empty."))
  }
  
  setDT(qc_df)
  failed_df <- qc_df[note == "mhc" | stage == "FAILED"]
  if (nrow(failed_df) > 0) {
      failed_loci <- nrow(failed_df)
      fail_file <- glue::glue("{analysis_folder}/{sample_id}_SuSiE_failed_loci.tsv")
      fwrite(failed_df, fail_file, sep = "\t")
      msg(verbose, glue::glue(
          "⚠ {failed_loci} loci failed or were skipped; see {fail_file} for details."
      ))
  } else {
      msg(verbose, "✅ No failed or MHC-skipped loci found.")
  }


  msg(verbose, glue::glue(
    "✅ SuSiE fine-mapping completed.\n",
    "  • Combined results: {out_combined}\n",
    "  • QC summary:      {qc_file}"
  ))

  return(df_vars)

  # invisible(list(
  #   vars = df_vars,
  #   qc   = if (exists("qc_df")) qc_df else data.frame()
  # ))

}



## ── Helper: safe %||% operator ─────────────────────────────────────
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

## ── Helper: message wrapper ───────────────────────────────────────
msg <- function(verbose, ...) if (isTRUE(verbose)) message(...)



postgwas_ld_matrix <- function(variants,
                            bfile,
                            plink_bin,
                            output_folder,
                            tag,
                            with_alleles = TRUE,
                            logfile = NULL) {

        # -------------------------------------------------------
        # Create output subfolder for LD
        # -------------------------------------------------------
        ld_dir <- file.path(output_folder, "ld_matrix_related")
        dir.create(ld_dir, recursive = TRUE, showWarnings = FALSE)

        # -------------------------------------------------------
        # Create prefix (per-locus output files)
        # -------------------------------------------------------
        safe_tag <- gsub(":", "-", tag)
        prefix <- file.path(ld_dir, paste0("ld_run_", safe_tag))

        # -------------------------------------------------------
        # Logfile setup
        # -------------------------------------------------------
        if (is.null(logfile)) {
        logfile <- file.path(ld_dir, "ld_matrix_local.log")
        }

        file.create(logfile)   # Start clean or overwrite existing


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

      # Detect OS for quoting rules
      shell <- ifelse(Sys.info()[['sysname']] == "Windows", "cmd", "sh")

      # -------------------------------------------------------
      # PLINK call 1 → generate .bim (quiet)
      # -------------------------------------------------------
      fun1 <- paste0(
        shQuote(plink_bin, type = shell),
        " --bfile ", shQuote(bfile, type = shell),
        " --extract ", shQuote(snp_list_file, type = shell),
        " --make-just-bim ",
        " --keep-allele-order ",
        " --out ", shQuote(prefix, type = shell),
        " >> ", shQuote(logfile), " 2>&1"
      )

      system(fun1)

      # Read BIM produced by PLINK
      bim <- read.table(paste0(prefix, ".bim"), stringsAsFactors = FALSE)

      # -------------------------------------------------------
      # PLINK call 2 → LD matrix (quiet)
      # -------------------------------------------------------
      fun2 <- paste0(
        shQuote(plink_bin, type = shell),
        " --bfile ", shQuote(bfile, type = shell),
        " --extract ", shQuote(snp_list_file, type = shell),
        " --r square ",
        " --keep-allele-order ",
        " --out ", shQuote(prefix, type = shell),
        " >> ", shQuote(logfile), " 2>&1"
      )

      system(fun2)

      # -------------------------------------------------------
      # Read LD matrix
      # -------------------------------------------------------
      res <- as.matrix(read.table(paste0(prefix, ".ld"), header = FALSE))

      # Add row/column names
      if (with_alleles) {
        rownames(res) <- colnames(res) <- paste(bim$V2, bim$V5, bim$V6, sep = "_")
      } else {
        rownames(res) <- colnames(res) <- bim$V2
      }

      # -------------------------------------------------------
      # CLEANUP FILES
      # -------------------------------------------------------
      temp_files <- Sys.glob(paste0(prefix, "*"))
      if (length(temp_files) > 0) {
        unlink(temp_files)
      }

      return(res)
    }


    ## ── 1. RAM detection (100% parallel-safe) ───────────────────────
    get_free_memory_gb <- function() {
      os <- Sys.info()[["sysname"]]
      
      if (os == "Darwin") {
        out <- suppressWarnings(system("sysctl hw.memsize", intern = TRUE, ignore.stderr = TRUE))
        if (length(out) == 0) return(Inf)
        total <- as.numeric(sub("hw.memsize: ([0-9]+)", "\\1", out))
        if (is.na(total)) return(Inf)
        return(0.25 * total / 1024^3)  # conservative guess
      }
      
      if (os == "Linux") {
        mem <- tryCatch(readLines("/proc/meminfo"), error = function(e) character(0))
        avail <- grep("^MemAvailable:", mem, value = TRUE)
        if (length(avail) == 0) return(Inf)
        kb <- as.numeric(sub(".* ([0-9]+) kB.*", "\\1", avail))
        return(kb / 1024^2)
      }
      
      Inf  # Windows / unknown → pretend infinite RAM
}

## ── 2. Worker auto-detection ─────────────────────────────────────
auto_detect_workers <- function(min_ram_per_worker_gb = 4, reserve_cores = 1, verbose = TRUE) {
    phys_cores <- parallel::detectCores(logical = FALSE) %||% parallel::detectCores()
    max_by_cores <- max(1, phys_cores - reserve_cores)
    
    free_gb <- tryCatch(get_free_memory_gb(), error = function(e) Inf)
    if (free_gb == Inf || free_gb < 1) {
      if (verbose) message("Free RAM unknown → using ", max_by_cores, " workers (core-based)")
      return(max_by_cores)
    }
    
    max_by_ram <- floor(free_gb / min_ram_per_worker_gb)
    workers <- max(1, min(max_by_cores, max_by_ram))
    
    if (verbose) message("Auto workers: ~", round(free_gb, 1), " GB free → ", workers, " workers")
    workers
}

## ── 3. Memory throttle (Linux only) ─────────────────────────────
wait_for_memory <- function(threshold_used_pct = 85, check_interval_sec = 10, verbose = TRUE) {
    if (Sys.info()[["sysname"]] %in% c("Darwin", "Windows")) return(invisible(TRUE))
    
    repeat {
      mem <- tryCatch(readLines("/proc/meminfo", n = 30), error = function(e) character(0))
      if (length(mem) == 0) { Sys.sleep(check_interval_sec); next }
      
      get_kb <- function(k) as.numeric(sub(".* ([0-9]+) kB.*", "\\1", mem[grepl(paste0("^", k, ":"), mem)])) %||% 0
      
      total <- get_kb("MemTotal")
      available <- get_kb("MemFree") + get_kb("Buffers") + get_kb("Cached") + get_kb("SReclaimable")
      
      if (total == 0) { Sys.sleep(check_interval_sec); next }
      
      used_pct <- 100 * (1 - available / total)
      if (used_pct < threshold_used_pct) break
      
      if (verbose) {
        message("High memory: ", round(used_pct, 1), "% used — waiting ", check_interval_sec, "s...")
      }
      Sys.sleep(check_interval_sec)
    }
    invisible(TRUE)
}

## ── 4. Timeout wrapper ───────────────────────────────────────────
run_with_timeout <- function(expr, timeout, step, tag, verbose = TRUE) {
    if (!is.finite(timeout)) return(tryCatch(expr, error = function(e) NULL))
    
    tryCatch(
      R.utils::withTimeout(expr, timeout = timeout, onTimeout = "error"),
      TimeoutException = function(e) {
        msg(verbose, "Timeout (>", timeout, "s) during ", step, " for ", tag)
        structure(list(timeout = TRUE), class = "timeout")
      },
      error = function(e) {
        msg(verbose, "Error during ", step, " for ", tag, ": ", e$message)
        NULL
      }
    )
}

## ── 5. Fast LD repair ───────────────────────────────────────────
LD_fix_fast <- function(ld.mat, tag = "", verbose = TRUE) {
    ld.mat <- (ld.mat + t(ld.mat)) / 2
    start_t <- Sys.time()
    msg(verbose, "Fast LD repair started for ", tag)
    
    e <- eigen(ld.mat, symmetric = TRUE)
    e$values[e$values < 1e-8] <- 1e-8
    repaired <- e$vectors %*% diag(e$values) %*% t(e$vectors)
    
    if (isTRUE(tryCatch({chol(repaired); TRUE}, error = function(e) FALSE))) {
      msg(verbose, "Fast LD repair OK (", round(difftime(Sys.time(), start_t, units = "secs"), 2), "s)")
      return(repaired)
    }
    
    msg(verbose, "Fast repair failed → falling back to nearPD")
    npd <- Matrix::nearPD(ld.mat, corr = TRUE, keepDiag = TRUE)
    as.matrix(npd$mat)
}



## ── 6. Main locus processor (PASS 1) ─────────────────────────────
process_locus <- function(locus, df, sample_id, ld_ref, plink, analysis_folder,
                          lp_threshold, L, timeout_ld, timeout_susie,
                          skip_mhc, mhc_start, mhc_end, verbose = TRUE) {
    tag <- glue::glue("{locus$chr}:{locus$start}-{locus$end}")
    msg(verbose, "Processing locus ", tag)

    # MHC skip
    if (skip_mhc && locus$chr == "6" &&
        locus$start < mhc_end && locus$end > mhc_start) {
      msg(verbose, "Skipping MHC region")
      return(list(status = "SKIP", reason = "mhc", locus = locus))
    }

    wait_for_memory(verbose = verbose)

    # Subset sumstats
    idx <- which(df$CHR == as.numeric(locus$chr) &
                df$BP >= locus$start & df$BP <= locus$end)
    if (length(idx) == 0) return(list(status = "SKIP", reason = "no_snps", locus = locus))

    selected <- data.table::as.data.table(df[idx, ])
    if (max(selected$LP, na.rm = TRUE) <= lp_threshold)
      return(list(status = "SKIP", reason = "low_signal", locus = locus))

  # LD matrix
  # Run ld_matrix_local
  ld <- postgwas_ld_matrix(
        variants = selected$SNP,
        bfile = ld_ref,
        plink = plink,
        tag=gsub(":", "-", tag),
        with_alleles = FALSE,
        logfile=glue::glue("{analysis_folder}/ld_matrix_related/{sample_id}_{tag}_ld_matrix_local.log"),
        output_folder=analysis_folder )

    selected <- selected[match(rownames(ld), SNP), ]
    set(selected, j = "variable_index", value = seq_len(nrow(selected)))

    if (!identical(selected$SNP, rownames(ld)))
      stop("SNP order mismatch after LD extraction")

    z <- selected$EZ
    n_eff <- median(selected$NEF, na.rm = TRUE)
    if (!is.finite(n_eff) || n_eff <= 0) n_eff <- 1e5

    # SuSiE
    fitted <- tryCatch(
      susieR::susie_rss(z = z, R = ld, n = n_eff, L = L, max_iter = 100,return_correlation = TRUE),
      error = function(e) e
    )

    if (inherits(fitted, "error") || is.null(fitted) || !isTRUE(fitted$converged)) {
      return(list(
        status = "RECOVERY",
        locus = locus, selected = selected, ld.mat = ld, n_eff = n_eff
      ))
    }

      sus <- summary(fitted)
      vars <- if (!is.null(sus$vars)) as.data.table(sus$vars) else data.table()
      cs   <- if (!is.null(sus$cs))   as.data.table(sus$cs)   else data.table()

      if (nrow(vars) > 0 && "variable" %in% names(vars)) {
        setnames(vars, "variable", "variable_index")
        vars[, variable_index := as.integer(variable_index)]
        vars[, `:=`(
          SNP = selected$SNP[variable_index],
          CHR = selected$CHR[variable_index],
          BP  = selected$BP[variable_index],
          REF = selected$REF[variable_index],
          ALT = selected$ALT[variable_index]
        )]
      }

      selected_sus <- merge(selected, vars, by = c("variable_index", "SNP", "CHR", "BP", "REF", "ALT"), all.x = TRUE)
      selected_sus[, `:=`(
        locus_chr   = locus$chr[1L],
        locus_start = locus$start[1L],
        locus_end   = locus$end[1L]
      )]

      cred_df <- selected_sus[!is.na(cs) & cs != -1]

    # Save results
    system(glue::glue("mkdir -p {analysis_folder}/locus_files"), ignore.stdout = TRUE)
    system(glue::glue("mkdir -p {analysis_folder}/rds_files"), ignore.stdout = TRUE)
    cred_file <- glue::glue("{analysis_folder}/locus_files/{sample_id}_{tag}_SUSIE_credible_sets.cred1")
    fwrite(cred_df, cred_file, sep = " ", quote = FALSE)
    saveRDS(fitted, glue::glue("{analysis_folder}/rds_files/{sample_id}_{tag}_SUSIE.rds"))
    #list(status = "OK", fitted = fitted, locus = locus, cred_file = cred_file)
    list(status = "OK", vars = selected_sus, cs = cred_df, fitted = fitted, locus = locus, cred_file = cred_file)
}





## ── 7. Recovery (sequential) ─────────────────────────────────────
## ================================================================
## Module — recover_locus (PASS 2, sequential with timeouts + LD repair)
## ================================================================
recover_locus <- function(
  job,
  sample_id,
  analysis_folder,
  L,
  timeout_ld,
  timeout_susie,
  verbose = TRUE
) {
      locus    <- job$locus
      selected <- job$selected
      ld.mat   <- job$ld.mat
      n_eff    <- job$n_eff
      tag      <- glue::glue("{locus$chr}:{locus$start}-{locus$end}")

      start_rec <- Sys.time()
      msg(verbose, "Recovery started for ", tag)

      
      ## Helper: convert SuSiE fit → standard output
      sus_to_out <- function(fitted, status, reason= " ", suffix = "") {
        sus <- summary(fitted)
        vars <- if (!is.null(sus$vars)) as.data.table(sus$vars) else data.table()
        cs   <- if (!is.null(sus$cs))   as.data.table(sus$cs)   else data.table()

        if (nrow(vars) > 0 && "variable" %in% names(vars)) {
          setnames(vars, "variable", "variable_index")
          vars[, variable_index := as.integer(variable_index)]
          vars[, `:=`(
            SNP = selected$SNP[variable_index],
            CHR = selected$CHR[variable_index],
            BP  = selected$BP[variable_index],
            REF = selected$REF[variable_index],
            ALT = selected$ALT[variable_index]
          )]
        }
        selected_sus <- merge(selected, vars, by = c("variable_index", "SNP", "CHR", "BP", "REF", "ALT"), all.x = TRUE)
        selected_sus[, `:=`(
          locus_chr   = locus$chr[1L],
          locus_start = locus$start[1L],
          locus_end   = locus$end[1L]
        )]

        cred_df <- selected_sus[!is.na(cs) & cs != -1]
        system(glue::glue("mkdir -p {analysis_folder}/locus_files"), ignore.stdout = TRUE)
        system(glue::glue("mkdir -p {analysis_folder}/rds_files"), ignore.stdout = TRUE)
        
        cred_file <- glue::glue("{analysis_folder}/locus_files/{sample_id}_{tag}_SUSIE_credible_sets{suffix}.cred1")
        fwrite(cred_df, cred_file, sep = " ", quote = FALSE)
        saveRDS(fitted, glue::glue("{analysis_folder}/rds_files/{sample_id}_{tag}_SUSIE{suffix}.rds"))
        list(status = status, reason = reason, vars = selected_sus, cs = cred_df, fitted = fitted, locus = locus, cred_file = cred_file)
      }



      ## Step 0 – SuSiE attempt 1 (higher iterations)
      msg(verbose, "   Step0/5: SuSiE attempt 0 (L=", L, ", max_iter=600)")
      fit0 <- run_with_timeout(
        susieR::susie_rss(z = selected$EZ, R = ld.mat, n = n_eff, L = L, max_iter = 600),
        timeout = timeout_susie * 3,
        step    = "SuSiE attempt 1",
        tag     = tag,
        verbose = verbose
      )

      if (!is.null(fit0) && !inherits(fit0, "timeout") && isTRUE(fit0$converged)) {
        elapsed <- round(difftime(Sys.time(), start_rec, units = "secs"), 1)
        msg(verbose, "   Recovery SUCCESS (attempt 0) in ", elapsed, "s for ", tag)
        return(sus_to_out(fit0,status=glue("susie_rss is sucess with L {L} and max_iter 600 ")))
      }



      ## Step 1 – SuSiE attempt 2 (higher iterations + lower L)
      msg(verbose, "   Step1/4: SuSiE attempt 0 (L=", 5, ", max_iter=600)")
      fit0 <- run_with_timeout(
        susieR::susie_rss(z = selected$EZ, R = ld.mat, n = n_eff, L = 5, max_iter = 600),
        timeout = timeout_susie * 3,
        step    = "SuSiE attempt 2",
        tag     = tag,
        verbose = verbose
      )

      if (!is.null(fit0) && !inherits(fit0, "timeout") && isTRUE(fit0$converged)) {
        elapsed <- round(difftime(Sys.time(), start_rec, units = "secs"), 1)
        msg(verbose, "   Recovery SUCCESS (attempt 1) in ", elapsed, "s for ", tag)
        return(sus_to_out(fit0,status=glue("susie_rss is sucess with L 5 and max_iter 600 ")))

      }


      ## Step 2 – LD repair with timeout
      msg(verbose, "   Step 2/4: LD repair (eigen + nearPD fallback)")
      wait_for_memory(verbose = verbose)

      ld_repaired <- run_with_timeout(
        expr    = LD_fix_fast(ld.mat, tag = tag, verbose = verbose),
        timeout = timeout_ld,
        step    = "LD repair",
        tag     = tag,
        verbose = verbose
      )

      if (is.null(ld_repaired) || inherits(ld_repaired, "timeout")) {
        msg(verbose, "   Recovery FAILED at LD repair for ", tag)
        return(list(status = "FAILED", locus = locus, reason = "LD_repair_timeout_or_error"))
      }


      ## Step 3 – SuSiE attempt 1 (higher iterations)
      msg(verbose, "   Step 3/4: SuSiE attempt 1 (L=", L, ", max_iter=300)")
      fit1 <- run_with_timeout(
        susieR::susie_rss(z = selected$EZ, R = ld_repaired, n = n_eff, L = L, max_iter = 300),
        timeout = timeout_susie,
        step    = "SuSiE attempt 1",
        tag     = tag,
        verbose = verbose
      )

      if (!is.null(fit1) && !inherits(fit1, "timeout") && isTRUE(fit1$converged)) {
        elapsed <- round(difftime(Sys.time(), start_rec, units = "secs"), 1)
        msg(verbose, "   Recovery SUCCESS (attempt 1) in ", elapsed, "s for ", tag)
        return(sus_to_out(fit1,status=glue("susie_rss is sucess with repaired ld and  L {L} and max_iter 300 ")))

      }


      ## Step 3 – SuSiE attempt 2 (reduced L)
      L_new <- max(3, floor(L / 2))
      msg(verbose, "   Step 4/4: SuSiE attempt 2 (L=", L_new, ", max_iter=500)")
      fit2 <- run_with_timeout(
        susieR::susie_rss(z = selected$EZ, R = ld_repaired, n = n_eff, L = L_new, max_iter = 600),
        timeout = timeout_susie * 3,   # give it a bit more time
        step    = "SuSiE attempt 2",
        tag     = tag,
        verbose = verbose
      )

      if (!is.null(fit2) && !inherits(fit2, "timeout") && isTRUE(fit2$converged)) {
        elapsed <- round(difftime(Sys.time(), start_rec, units = "secs"), 1)
        msg(verbose, "   Recovery SUCCESS (attempt 2, L=", L_new, ") in ", elapsed, "s for ", tag)
        return(sus_to_out(fit2, suffix = glue::glue("_L{L_new}")))
        return(sus_to_out(fit2,status=glue("susie_rss is sucess with repaired ld and  L {L} and max_iter 600 ")))

      }


      ## Final failure
      msg(verbose, "   Recovery FAILED for ", tag)
      list(status = "FAILED", locus = locus, reason = "all_recovery_attempts_failed")
  }





## ── 8. MAIN FUNCTION (FINAL FIXED VERSION) ───────────────────────
run_susie_finemap_parallel <- function(
  locus_file, sumstat_file, sample_id, ld_ref, plink, SUSIE_Analysis_folder,
  lp_threshold = 7.3, L = 10, 
  workers = "auto",                     # ← changed here
  min_ram_per_worker_gb = 4,
  verbose = TRUE, timeout_ld_seconds = 180, timeout_susie_seconds = 180,
  skip_mhc = TRUE, mhc_start = 25e6, mhc_end = 35e6
) {
    suppressPackageStartupMessages({
      library(data.table); library(glue); library(susieR); library(Matrix)
      library(future); library(future.apply); library(R.utils) ; library(purrr)
    })

    dir.create(SUSIE_Analysis_folder, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(SUSIE_Analysis_folder, "locus_files"), recursive = TRUE, showWarnings = FALSE)

    locus_df <- fread(locus_file)[, .(chr = as.character(CHR), start = as.integer(START), end = as.integer(END))]
    df <- fread(sumstat_file)
    setnames(df, old = c("SNP","REF","ALT","CHR"), new = c("SNP","REF","ALT","CHR"), skip_absent = TRUE)
    df[, CHR := as.character(CHR)]

    options(future.globals.maxSize = +Inf)
    options(future.seed = TRUE)


    ## ── Resolve workers ─────────────────────────────────────────────────
    if (identical(workers, "auto")) {
      workers <- auto_detect_workers(min_ram_per_worker_gb = min_ram_per_worker_gb, verbose = verbose)
    } else if (!is.numeric(workers) || workers < 1) {
      stop("workers must be 'auto' or a positive integer")
    }
    workers <- as.integer(workers)

    ## ── Parallel backend (macOS-safe) ───────────────────────────────────
    os <- Sys.info()[["sysname"]]
    if (os == "Darwin") {
      msg(verbose, "macOS detected → using multicore backend")
      future::plan(multicore, workers = workers, gc = TRUE)
    } else {
      msg(verbose, "Non-macOS → using multisession backend")
      future::plan(multisession, workers = workers, gc = TRUE)
    }

    msg(verbose, "Starting PASS 1 with ", workers, " workers on ", nrow(locus_df), " loci")

    pass1 <- future_lapply(seq_len(nrow(locus_df)), function(i) {
      process_locus(
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
        mhc_start = mhc_start, 
        mhc_end = mhc_end,
        verbose = verbose
      )
    }, future.seed = TRUE)

    # Set to run 2 jobs in parallel
    plan(multisession, workers = max(1, floor(workers / 2)) )   # works on macOS + Ubuntu

    rec_jobs <- purrr::keep(pass1, ~ .x$status == "RECOVERY")
    msg(verbose, length(rec_jobs), " loci need recovery")
    rec_results <- if (length(rec_jobs)) {
      future_lapply(
        rec_jobs,
        recover_locus,
        sample_id = sample_id,
        analysis_folder = SUSIE_Analysis_folder,
        L = L,
        timeout_ld = timeout_ld_seconds,
        timeout_susie = timeout_susie_seconds,
        verbose = verbose,
        future.seed = TRUE
      )
    } else list()

    credible_set_df=aggregate_and_write_results(c(pass1, rec_results), sample_id, SUSIE_Analysis_folder, verbose)
    print(names(credible_set_df))
    plot_susie_loci_advanced(
        df=credible_set_df,
        sample_id   = sample_id,
        outdir      = glue("{SUSIE_Analysis_folder}/susie_plots/"),
        n_cores     = workers,
        min_pip_label = 0.001  )
    future::plan(sequential)
    msg(verbose, "SuSiE fine-mapping finished for ", sample_id)
    invisible(TRUE)
}


# run_susie_finemap(
#   locus_file = "/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/PGC3_SCZ_european/PGC3_SCZ_european_LDpruned_EUR_sig.tsv",
#   sumstat_file = "/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/PGC3_SCZ_european/finemap/PGC3_SCZ_european_finemap.tsv",
#   sample_id = "PGC3_SCZ_european",
#   ld_ref = "/Users/JJOHN41/Documents/software_resources/resourses/postgwas/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001",
#   plink = "/Users/JJOHN41/Documents/software_resources/softwares/plink",
#   SUSIE_Analysis_folder = "/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/PGC3_SCZ_european/finemap/SuSiE_analysis",
#   lp_threshold = 7.3
# )


#   locus_file = "/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/PGC3_SCZ_european/PGC3_SCZ_european_LDpruned_EUR_sig.tsv"
#   sumstat_file = "/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/PGC3_SCZ_european/finemap/PGC3_SCZ_european_finemap.tsv"
#   sample_id = "PGC3_SCZ_european"
#   ld_ref = "/Users/JJOHN41/Documents/software_resources/resourses/postgwas/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001"
#   plink = "/Users/JJOHN41/Documents/software_resources/softwares/plink"
#   SUSIE_Analysis_folder = "/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/PGC3_SCZ_european/finemap/SuSiE_analysis"
#   lp_threshold = 7.3
#  L = 10
#   workers = "auto"                  
#   min_ram_per_worker_gb = 4
#   verbose = TRUE
#   timeout_ld_seconds = 180 
#   timeout_susie_seconds = 180
#   skip_mhc = TRUE
#   mhc_start = 25e6
#   mhc_end = 35e6



# source('/Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/finemapping/susie_test.r', chdir = TRUE)
# source('/Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/finemapping/susie.r', chdir = TRUE)

# run_susie_finemap_parallel(
#   locus_file = "/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/PGC3_SCZ_european/PGC3_SCZ_european_LDpruned_EUR_sig.tsv",
#   sumstat_file = "/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/PGC3_SCZ_european/finemap/PGC3_SCZ_european_finemap.tsv",
#   sample_id = "PGC3_SCZ_european",
#   ld_ref = "/Users/JJOHN41/Documents/software_resources/resourses/postgwas/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001",
#   plink = "/Users/JJOHN41/Documents/software_resources/softwares/plink",
#   SUSIE_Analysis_folder = "/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/PGC3_SCZ_european/finemap/SuSiE_analysis",
#   lp_threshold = 7.3
# )



# locus         = locus_df[i, ]
# df            = df
# sample_id     = sample_id
# ld_ref        = ld_ref
# plink         = plink
# analysis_folder = SUSIE_Analysis_folder
# lp_threshold  = lp_threshold
# L             = 10
# timeout_ld    = 300
# timeout_susie = 300
# skip_mhc      = TRUE

      