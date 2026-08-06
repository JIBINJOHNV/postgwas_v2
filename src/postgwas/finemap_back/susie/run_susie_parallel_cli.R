#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(glue)
})


# =============================================================
# Load auxiliary susie.r located in the same directory
# =============================================================

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)

  # Find '--file=' argument used by Rscript
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) == 0) {
    stop("Unable to determine script location using commandArgs()", call. = FALSE)
  }

  script_path <- sub("^--file=", "", file_arg)
  normalizePath(dirname(script_path))
}

script_dir <- get_script_dir()
susie_file <- file.path(script_dir, "susie.r")

if (!file.exists(susie_file)) {
  stop(glue::glue("Required file not found: {susie_file}"), call. = FALSE)
}

source(susie_file)



# -------------------------------------------------------------
# Define CLI
# -------------------------------------------------------------
parser <- ArgumentParser(
  description = "Parallel SuSiE fine-mapping with LD matrix construction",
  formatter_class = "argparse.RawTextHelpFormatter"
)

parser$add_argument("--locus_file", required = TRUE,
                    help = "Input locus file (tsv): chr, start, end per row")
parser$add_argument("--sumstat_file", required = TRUE,
                    help = "GWAS summary statistics for fine-mapping")
parser$add_argument("--sample_id", required = TRUE,
                    help = "Sample identifier for output files")

parser$add_argument("--ld_ref", required = TRUE,
                    help = "Reference PLINK prefix for LD (1000G EUR, etc.)")
parser$add_argument("--plink", required = TRUE,
                    help = "Path to plink binary")

parser$add_argument("--SUSIE_Analysis_folder", required = TRUE,
                    help = "Output folder for SuSiE results")

parser$add_argument("--lp_threshold", type = "double", default = 7.3,
                    help = "LP threshold for variant filtering (default = 7.3)")
parser$add_argument("--L", type = "integer", default = 10,
                    help = "Maximum number of causal components (default = 10)")

parser$add_argument("--workers", default = "auto",
                    help = "Number of workers: integer or 'auto' (default = auto)")
parser$add_argument("--min_ram_per_worker_gb", type = "double", default = 4,
                    help = "Min RAM per worker in GB when using auto workers (default = 4)")

parser$add_argument("--timeout_ld_seconds", type = "integer", default = 180,
                    help = "Timeout for LD computation per locus (default 180s)")
parser$add_argument("--timeout_susie_seconds", type = "integer", default = 180,
                    help = "Timeout for SuSiE per locus (default 180s)")

parser$add_argument("--skip_mhc", action = "store_true",
                    help = "Skip MHC region (default TRUE)")

parser$add_argument("--finemap_mhc_chrom",default = '6',
                    help = "Chromosome of MHC region (default 6)")

parser$add_argument("--mhc_start", type = "double", default = 25e6,
                    help = "MHC start position")
parser$add_argument("--mhc_end", type = "double", default = 35e6,
                    help = "MHC end position")

parser$add_argument("--verbose", action = "store_true",
                    help = "Enable verbose logging")

args <- parser$parse_args()

# -------------------------------------------------------------
# Load your function
# -------------------------------------------------------------

# -------------------------------------------------------------
# Run
# -------------------------------------------------------------
run_susie_finemap_parallel(
  locus_file            = args$locus_file,
  sumstat_file          = args$sumstat_file,
  sample_id             = args$sample_id,
  ld_ref                = args$ld_ref,
  plink                 = args$plink,
  SUSIE_Analysis_folder = args$SUSIE_Analysis_folder,
  lp_threshold          = args$lp_threshold,
  L                     = args$L,
  workers               = args$workers,
  min_ram_per_worker_gb = args$min_ram_per_worker_gb,
  verbose               = args$verbose,
  timeout_ld_seconds    = args$timeout_ld_seconds,
  timeout_susie_seconds = args$timeout_susie_seconds,
  skip_mhc              = args$skip_mhc,
  mhc_chrom             = args$finemap_mhc_chrom,
  mhc_start             = args$mhc_start,
  mhc_end               = args$mhc_end
)
