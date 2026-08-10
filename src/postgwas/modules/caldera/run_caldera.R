args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(paste(
    "Expected: REPOSITORY UPSTREAM_SCRIPT POPS_FILE CREDIBLE_SET_FILE",
    "ASSEMBLY OUTPUT_FILE DELIMITER"
  ))
}

repository <- normalizePath(args[[1L]], mustWork = TRUE)
upstream_script <- normalizePath(args[[2L]], mustWork = TRUE)
pops_file <- normalizePath(args[[3L]], mustWork = TRUE)
credible_set_file <- normalizePath(args[[4L]], mustWork = TRUE)
assembly <- suppressWarnings(as.integer(args[[5L]]))
output_file <- args[[6L]]
delimiter <- args[[7L]]

source(upstream_script)
result <- caldera(
  pops_file = pops_file,
  cs_file = credible_set_file,
  assembly = assembly,
  caldera_path = repository
)
data.table::fwrite(result, output_file, sep = delimiter, quote = FALSE, na = "NA")
