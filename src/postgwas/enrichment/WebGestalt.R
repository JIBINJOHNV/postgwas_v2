run_WebGestalt_ORA_NTA <- function(
  my_genes,
  outputDirectory,
  background_genes = NULL,
  organism = "hsapiens",
  ora_sigMethod = "fdr",
  ora_fdrMethod = "BH",
  ora_fdrThr = 0.99,
  nta_fdrThr = 0.05,
  nThreads = 1
) {

  library('WebGestaltR')
  if (!dir.exists(outputDirectory)) {
    dir.create(outputDirectory, recursive = TRUE)
  }

  message("🧬 Number of input genes: ", length(my_genes)

  )

  # ============================================================
  # ORA PART
  # ============================================================
  message("\n========== Running ORA ==========")

  gs <- listGeneSet()
  gs_filtered <- gs[!grepl("network", gs$name, ignore.case = TRUE), ]

  for (db in gs_filtered$name) {
    message("Running ORA for database: ", db)

    try({
      WebGestaltR(
        enrichMethod      = "ORA",
        organism          = organism,
        enrichDatabase    = db,

        interestGene      = my_genes,
        interestGeneType  = "genesymbol",

        # Background handling
        referenceSet      = if (is.null(background_genes)) "genome_protein-coding" else NULL,
        referenceGene     = background_genes,
        referenceGeneType = if (!is.null(background_genes)) "genesymbol" else NULL,

        # ORA parameters
        minNum            = 1,
        maxNum            = 5000,
        sigMethod         = ora_sigMethod,
        fdrMethod         = ora_fdrMethod,
        fdrThr            = ora_fdrThr,
        reportNum         = 20,

        # Output
        isOutput          = TRUE,
        outputDirectory   = outputDirectory,
        projectName       = db,

        # Misc
        nThreads          = nThreads,
        hostName          = "https://www.webgestalt.org/"
      )
    }, silent = TRUE)
  }

  # ============================================================
  # NTA PART
  # ============================================================
  message("\n========== Running NTA ==========")

  int_databases <- c(
    "network_CORUM",
    "network_CORUMA",
    "network_FunMap",
    "network_FunMap_DenseModules",
    "network_FunMap_HierarchalModules",
    "network_Kinase_phosphosite",
    "network_Kinase_target",
    "network_PPI_BIOGRID"
  )

  for (db in int_databases) {
    message("Running NTA for database: ", db)

    try({
      WebGestaltR(
        enrichMethod       = "NTA",
        organism           = organism,
        enrichDatabase     = db,

        interestGene       = my_genes,
        interestGeneType   = "genesymbol",

        # REQUIRED for NTA
        networkConstructionMethod = "Network_Retrieval_Prioritization",

        # Network parameters
        neighborNum        = 10,
        highlightType      = "Seeds",
        highlightSeedNum   = 10,

        # Significance
        sigMethod          = "fdr",
        fdrMethod          = "BH",
        fdrThr             = nta_fdrThr,

        # Clustering
        setCoverNum        = 10,
        useWeightedSetCover = FALSE,
        useAffinityPropagation = FALSE,
        usekMedoid         = TRUE,
        kMedoid_k          = 25,

        # Output
        isOutput           = TRUE,
        outputDirectory    = outputDirectory,
        projectName        = paste0("NTA_", db),

        # Misc
        nThreads           = nThreads,
        hostName           = "https://www.webgestalt.org/"
      )
    }, silent = TRUE)
  }

  message("\n✅ WebGestalt ORA + NTA analysis completed.")
}


# run_WebGestalt_ORA_NTA(
#   my_genes = c("TP53", "BRCA1", "EGFR", "TNF", "APOE", "IL6"),
#   outputDirectory = "/Users/JJOHN41/Documents/developing_software/data/oudir/PGC3_SCZ_european/enrichment_analysis/webgestalt/",
#   nThreads=10
# )
