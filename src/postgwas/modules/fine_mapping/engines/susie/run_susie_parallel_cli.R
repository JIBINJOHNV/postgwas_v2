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
parser$add_argument("--sumstat_manifest", required = TRUE,
                    help = "Manifest mapping each locus to its exact summary-statistics file")
parser$add_argument("--sample_id", required = TRUE,
                    help = "Sample identifier for output files")

parser$add_argument("--ld_ref", required = TRUE,
                    help = "Reference PLINK prefix for LD (1000G EUR, etc.)")
parser$add_argument("--plink", required = TRUE,
                    help = "Path to plink binary")

parser$add_argument("--SUSIE_Analysis_folder", required = TRUE,
                    help = "Output folder for SuSiE results")
parser$add_argument(
  "--resolved_configuration_file", required = TRUE,
  help = "Schema-validated resolved fine-mapping configuration JSON"
)

args <- parser$parse_args()

# -------------------------------------------------------------
# Load your function
# -------------------------------------------------------------

# -------------------------------------------------------------
# Run
# -------------------------------------------------------------
run_susie_finemap_parallel(
  locus_file            = args$locus_file,
  sumstat_manifest      = args$sumstat_manifest,
  sample_id             = args$sample_id,
  ld_ref                = args$ld_ref,
  plink                 = args$plink,
  SUSIE_Analysis_folder = args$SUSIE_Analysis_folder,
  configuration_file    = args$resolved_configuration_file
)
