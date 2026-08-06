
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
  timeout_seconds = 500
) {

  # ---- guard: empty variants ----
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

  shell <- ifelse(Sys.info()[["sysname"]] == "Windows", "cmd", "sh")

  # -------------------------------------------------------
  # Helper: run PLINK safely
  # -------------------------------------------------------
  run_plink_safe <- function(cmd, step_name) {
    if (requireNamespace("processx", quietly = TRUE)) {
      res <- tryCatch(
        processx::run(
          "sh",
          args = c("-c", cmd),
          timeout = timeout_seconds,
          error_on_status = TRUE
        ),
        error = function(e) {
          warning(step_name, " failed: ", e$message,
                  " | see log: ", logfile)
          NULL
        }
      )
      return(!is.null(res))
    } else {
      status <- system(cmd)
      if (status != 0) {
        warning(step_name, " failed (exit=", status,
                ") | see log: ", logfile)
        return(FALSE)
      }
      TRUE
    }
  }

  # -------------------------------------------------------
  # PLINK call 1 → generate .bim
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

  if (!run_plink_safe(fun1, "PLINK .bim generation")) {
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
  fun2 <- paste0(
    shQuote(plink_bin, type = shell),
    " --bfile ", shQuote(bfile, type = shell),
    " --extract ", shQuote(snp_list_file, type = shell),
    " --r square ",
    " --keep-allele-order ",
    " --out ", shQuote(prefix, type = shell),
    " >> ", shQuote(logfile), " 2>&1"
  )

  if (!run_plink_safe(fun2, "PLINK LD calculation")) {
    return(NULL)
  }

  ld_file <- paste0(prefix, ".ld")
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

  if (with_alleles) {
    rownames(res) <- colnames(res) <- paste(bim$V2, bim$V5, bim$V6, sep = "_")
  } else {
    rownames(res) <- colnames(res) <- bim$V2
  }

  # -------------------------------------------------------
  # Cleanup
  # -------------------------------------------------------
  unlink(Sys.glob(paste0(prefix, "*")), force = TRUE)

  res
}


run_susie_hard_timeout <- function(
  z, R, n, L, max_iter,
  timeout_s,
  tag,
  log_msg,
  susie_verbose = TRUE
) {
  # Absolute safety net — nothing escapes this function
  tryCatch({

    stopifnot(requireNamespace("processx", quietly = TRUE))

    tmp_input  <- tempfile(paste0("susie_input_", tag, "_"), fileext = ".rds")
    tmp_output <- tempfile(paste0("susie_out_",   tag, "_"), fileext = ".rds")
    tmp_script <- tempfile(paste0("susie_run_",   tag, "_"), fileext = ".R")

    tmp_stdout <- tempfile(paste0("susie_stdout_", tag, "_"), fileext = ".log")
    tmp_stderr <- tempfile(paste0("susie_stderr_", tag, "_"), fileext = ".log")

    on.exit({
      unlink(
        c(tmp_input, tmp_output, tmp_script, tmp_stdout, tmp_stderr),
        force = TRUE
      )
    }, add = TRUE)

    saveRDS(list(
      z = z,
      R = R,
      n = n,
      L = L,
      max_iter = max_iter,
      susie_verbose = susie_verbose
    ), tmp_input)

    writeLines(c(
      "suppressPackageStartupMessages(library(susieR))",
      sprintf("inp <- readRDS(%s)", deparse(tmp_input)),
      "fit <- susieR::susie_rss(",
      "  z = inp$z,",
      "  R = inp$R,",
      "  n = inp$n,",
      "  L = inp$L,",
      "  max_iter = inp$max_iter,",
      "  return_correlation = TRUE,",
      "  verbose = inp$susie_verbose",
      ")",
      sprintf("saveRDS(fit, %s)", deparse(tmp_output))
    ), tmp_script)

    log_msg(
      "SuSiE START | p=", length(z),
      " L=", L,
      " max_iter=", max_iter,
      " n_eff=", round(n, 2),
      " timeout=", timeout_s, "s"
    )

    rscript <- Sys.which("Rscript")
    if (!nzchar(rscript)) {
      log_msg("❌ SuSiE FAILED: Rscript not found in PATH")
      return(NULL)
    }

    t0 <- Sys.time()

    px <- processx::process$new(
      command = rscript,
      args = c("--vanilla", tmp_script),
      stdout = tmp_stdout,
      stderr = tmp_stderr,
      cleanup = TRUE
    )

    repeat {
      if (!px$is_alive()) break

      if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > timeout_s) {
        log_msg("❌ SuSiE HARD TIMEOUT → terminating PID ", px$get_pid())

        try(px$terminate(), silent = TRUE)
        Sys.sleep(0.5)

        if (px$is_alive()) {
          log_msg("❌ SuSiE still alive → killing PID ", px$get_pid())
          try(px$kill(), silent = TRUE)
        }

        return(NULL)
      }

      Sys.sleep(0.25)
    }

    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
    log_msg("SuSiE FINISHED in ", elapsed, "s")

    if (!file.exists(tmp_output)) {
      log_msg("❌ SuSiE FAILED: no output produced")

      if (file.size(tmp_stderr) > 0) {
        err <- readLines(tmp_stderr, warn = FALSE)
        log_msg("SuSiE stderr (tail):\n", paste(tail(err, 50), collapse = "\n"))
      }

      return(NULL)
    }

    fit <- readRDS(tmp_output)

    if (!inherits(fit, "susie") || !isTRUE(fit$converged)) {
      log_msg("❌ SuSiE FAILED: did not converge")
      return(NULL)
    }

    fit

  }, error = function(e) {
    # This guarantees the worker NEVER dies
    log_msg("❌ SuSiE INTERNAL ERROR: ", conditionMessage(e))
    NULL
  })
}

annotate_susie <- function(
  fitted,
  selected,
  locus
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
    }

    # ======================================================
    # 3. Merge selected SNPs with SuSiE vars
    # ======================================================
    selected_sus <- merge(
        selected, vars,
        by = c("variable_index", "SNP", "CHR", "BP", "REF", "ALT"),
        all.x = TRUE
    )
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
        for (cs_name in names(fitted$sets$cs)) {
        snp_idxs <- fitted$sets$cs[[cs_name]]
        cs_index[snp_idxs] <- cs_name
        }
    }
    selected_sus[, credible_set := cs_index]

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
        is_pure       = pur$min.abs.corr > 0.5
        )

        selected_sus <- merge(
        selected_sus,
        purity_dt,
        by.x = "credible_set",
        by.y = "cs",
        all.x = TRUE
        )
    }

    # ======================================================
    # 7. Add SuSiE model-level diagnostics
    # ======================================================
    selected_sus[, susie_converged := fitted$converged]
    selected_sus[, susie_niter     := fitted$niter]
    selected_sus[, susie_sigma2    := fitted$sigma2]
    selected_sus[, susie_log10BF1  := fitted$lbf[1] / log(10)]

    # ======================================================
    # 8. Add SNP-level SuSiE signals: PIP, mu, mu2, alpha, lbf
    # ======================================================
    ## 8A — Add PIP
    if (!is.null(fitted$pip)) {
        selected_sus[, global_pip := fitted$pip]
    }
    ## 8B — Add lbf_variable (L x J)
    if (!is.null(fitted$lbf_variable)) {
        L <- nrow(fitted$lbf_variable)
        for (k in 1:L) {
        selected_sus[[paste0("lbf_L", k)]] <- fitted$lbf_variable[k, ]
        }
    }
    ## 8C — Add mu
    if (!is.null(fitted$mu)) {
        L <- nrow(fitted$mu)
        for (k in 1:L) {
        selected_sus[[paste0("mu_L", k)]] <- fitted$mu[k, ]
        }
    }
    ## 8D — Add mu2
    if (!is.null(fitted$mu2)) {
        L <- nrow(fitted$mu2)
        for (k in 1:L) {
        selected_sus[[paste0("mu2_L", k)]] <- fitted$mu2[k, ]
        }
    }
    ## 8E — Add alpha
    if (!is.null(fitted$alpha)) {
        L <- nrow(fitted$alpha)
        for (k in 1:L) {
        selected_sus[[paste0("alpha_L", k)]] <- fitted$alpha[k, ]
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
    min_pip_label  = 0.10,
    bg_size        = 0.05,
    bg_alpha       = 0.15,
    cs_size        = 1.5,
    cs_alpha       = 0.15,
    max_ld_dim     = 5000   # max dimension for LD matrix (downsample if larger)
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
      n_col <- 256L
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
          y    = seq(0, 1, length.out = 100),
          fill = seq(0, 1, length.out = 100)
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
              legend.title     = element_text(size = 8),
              legend.text      = element_text(size = 7),
              legend.key.size  = unit(0.3, "cm")
          )

      ld_leg <- cowplot::get_legend(p_ld_leg)
      ld_leg$widths  <- ld_leg$widths  * 0.55
      ld_leg$heights <- ld_leg$heights * 0.55

      # ------------------------------------------------------------
      # CS LEGEND
      # ------------------------------------------------------------
      p_cs_leg_source <- ggplot(
          data.frame(cs_plot = factor(cs_levels, levels = cs_levels))
      ) +
          geom_point(aes(x = 1, y = cs_plot, color = cs_plot), size = 3) +
          scale_color_manual(values = pal, name = "Credible Set") +
          theme_void() +
          theme(
              legend.position  = "right",
              legend.title     = element_text(size = 8),
              legend.text      = element_text(size = 7),
              legend.key.size  = unit(0.3, "cm")
          )

      cs_leg <- cowplot::get_legend(p_cs_leg_source)
      cs_leg$widths  <- cs_leg$widths  * 0.55
      cs_leg$heights <- cs_leg$heights * 0.55

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
      ][order(-variable_prob)][, head(.SD, 5), by = cs]

      if (nrow(label_df) > 0)
          label_df[, cs_plot := factor(cs, levels = cs_levels)]

      n_snps <- nrow(locus_df)
      bg_size2 <- if (n_snps > 5000) bg_size*0.6 else if (n_snps > 2000) bg_size*0.8 else bg_size

      lift_amt <- 0.03

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
                  size = 3, segment.alpha = 0.7, max.overlaps = Inf
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
                          linewidth = 0.3
                      ),
                      geom_text_repel(
                          data = transform(label_df, variable_prob = variable_prob + lift_amt),
                          aes(x = BP, y = variable_prob, label = SNP, color = cs_plot),
                          size = 3, segment.alpha = 0.7, max.overlaps = Inf
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
      ) + patchwork::plot_layout(heights = c(2, 1, 1))

      combined <- left_panels | legend_column
      combined <- combined + patchwork::plot_layout(widths = c(5, 0.8))

      # ------------------------------------------------------------
      # SAVE PNG + PDF
      # ------------------------------------------------------------
      file_base <- glue("{outdir}/{sample_id}_{locus_id}_SUSIE_LD")

      ggsave(glue("{file_base}.png"), combined, width = 11, height = 12, dpi = 300)
      ggsave(glue("{file_base}.pdf"), combined, width = 11, height = 12)

      invisible(TRUE)
  }








generate_flames_files <- function(
  fitted,
  ld,
  snp_df,
  outfile 
) {
  suppressPackageStartupMessages({ require(data.table) })
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


aggregate_and_write_results <- function(
  results_all,
  sample_id,
  analysis_folder,
  verbose = TRUE
) {

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

  if (nrow(df_vars) > 0) {
    data.table::fwrite(df_vars, paste0(out_combined, ".gz"), compress = "gzip")
  } else {
    msg(verbose, glue::glue("⚠ No SuSiE results to write for {sample_id}; combined table is empty."))
  }

  if (length(qc_rows)) { 
    qc_df <- do.call(rbind, qc_rows)
    data.table::fwrite(qc_df, qc_file, sep = "\t")
  } else {
    qc_df <- data.frame()
  }

  if (nrow(df_cs) > 0) {
    data.table::fwrite(df_cs, out_credible) 
  } else {
    msg(verbose, glue::glue("⚠ No SuSiE credible set to write for {sample_id}; combined CS is empty."))
  }
  
  data.table::setDT(qc_df)
  failed_df <- qc_df[note == "mhc" | stage == "FAILED"]
  if (nrow(failed_df) > 0) {
    failed_loci <- nrow(failed_df)
    fail_file <- glue::glue("{analysis_folder}/{sample_id}_SuSiE_failed_loci.tsv")
    data.table::fwrite(failed_df, fail_file, sep = "\t")
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
