#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(
    paste(
      "Usage: export_rdata_columns.R INPUT_RDATA OUTPUT_TSV",
      "OBJECT_NAME COLUMN_1 COLUMN_2"
    ),
    call. = FALSE
  )
}

source_path <- args[[1]]
output_path <- args[[2]]
object_name <- args[[3]]
columns <- args[4:5]

objects <- new.env(parent = emptyenv())
loaded <- load(source_path, envir = objects)
if (!(object_name %in% loaded)) {
  stop(sprintf("RData object %s was not found in %s", object_name, source_path), call. = FALSE)
}

value <- objects[[object_name]]
missing_columns <- setdiff(columns, colnames(value))
if (length(missing_columns) > 0) {
  stop(
    sprintf("RData object is missing columns: %s", paste(missing_columns, collapse = ", ")),
    call. = FALSE
  )
}

write.table(
  value[, columns, drop = FALSE],
  file = output_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = ""
)
