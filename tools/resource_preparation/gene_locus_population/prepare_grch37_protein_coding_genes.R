#!/usr/bin/env Rscript

# Download Ensembl gene tables for GRCh37 and GRCh38 and write the two
# gene-resource layouts used by gene-locus/population-based analyses.

# User-editable settings.
OUTPUT_DIRECTORY <- "/Users/JJOHN41/Documents/software_resources/resourses/postgwas_v2/gene_locus_population"
ANNOTATION_FILENAME <- "protein_coding_autosome_X_gene_annotation.txt"
TSS_FILENAME <- "protein_coding_autosome_X_gene_TSS.txt"
ENSEMBL_DATASET <- "hsapiens_gene_ensembl"
GENOME_BUILDS <- c(GRCh37 = 37L, GRCh38 = 38L)
TARGET_CHROMOSOMES <- c(as.character(1:22), "X")

if (!requireNamespace("biomaRt", quietly = TRUE)) {
  stop("The biomaRt R package is required. Install it with BiocManager::install('biomaRt').")
}

dir.create(OUTPUT_DIRECTORY, recursive = TRUE, showWarnings = FALSE)

for (build_name in names(GENOME_BUILDS)) {
  mart <- biomaRt::useEnsembl(
    biomart = "genes",
    dataset = ENSEMBL_DATASET,
    GRCh = GENOME_BUILDS[[build_name]]
  )

  genes <- biomaRt::getBM(
    attributes = c(
      "ensembl_gene_id",
      "external_gene_name",
      "chromosome_name",
      "start_position",
      "end_position",
      "strand",
      "gene_biotype"
    ),
    mart = mart
  )

  genes <- genes[
    genes$gene_biotype == "protein_coding" &
      genes$chromosome_name %in% TARGET_CHROMOSOMES,
    ,
    drop = FALSE
  ]

  genes <- genes[
    order(
      match(genes$chromosome_name, TARGET_CHROMOSOMES),
      genes$start_position,
      genes$ensembl_gene_id
    ),
    ,
    drop = FALSE
  ]

  build_directory <- file.path(OUTPUT_DIRECTORY, build_name)
  dir.create(build_directory, recursive = TRUE, showWarnings = FALSE)

  annotation <- data.frame(
    ENSGID = genes$ensembl_gene_id,
    CHR = genes$chromosome_name,
    START = genes$start_position,
    END = genes$end_position,
    STRAND = ifelse(genes$strand == 1, "+", "-"),
    NAME = genes$external_gene_name,
    stringsAsFactors = FALSE
  )

  write.table(
    annotation,
    file.path(build_directory, ANNOTATION_FILENAME),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  tss <- data.frame(
    ENSGID = genes$ensembl_gene_id,
    NAME = genes$external_gene_name,
    CHR = genes$chromosome_name,
    START = genes$start_position,
    END = genes$end_position,
    TSS = ifelse(genes$strand == 1, genes$start_position, genes$end_position),
    stringsAsFactors = FALSE
  )

  write.table(
    tss,
    file.path(build_directory, TSS_FILENAME),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )

  message(sprintf("[%s] Wrote %d protein-coding genes to %s", build_name, nrow(genes), build_directory))
}
