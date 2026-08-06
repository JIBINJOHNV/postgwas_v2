
# https://dsigdb.tanlab.org/DSigDBv1.0/download.html
# https://pmc.ncbi.nlm.nih.gov/articles/PMC4668778/ 



import rpy2.robjects as ro
from rpy2.robjects.vectors import StrVector
from pathlib import Path


def run_dsigd_ora_webgestalt(
    genes,
    gmt_file,
    output_dir,
    background_genes=None,
    reference_set="genome_protein-coding",
    project_name="DSigDB_ORA",
    fdr_thr=0.05,
):
    """
    Run ORA using WebGestaltR with a custom GMT (e.g. DSigDB_All.gmt).

    Background logic:
      - if background_genes is provided → use referenceGene
      - else → use referenceSet (configurable, NOT hard-coded)
    Args:
        genes (list[str])
        gmt_file (str | Path)
        output_dir (str | Path)
        background_genes (list[str] | None)
        reference_set (str): WebGestalt reference set name
        project_name (str | None)
        fdr_thr (float)
    """
    gmt_file = str(Path(gmt_file).resolve())
    output_dir = str(Path(output_dir).resolve())
    
    # Load WebGestaltR
    ro.r("library(WebGestaltR)")
    # Ensure ZIP command is set correctly for the environment
    ro.r('Sys.setenv(R_ZIPCMD = "/opt/conda/envs/enricher/bin/zip")')

    # ---------------------------------------------------------
    # Convert inputs to R objects
    # ---------------------------------------------------------
    r_genes = StrVector(genes)
    
    # Handle background genes (list -> StrVector OR None -> ro.NULL)
    r_background = StrVector(background_genes) if background_genes else ro.NULL
    
    # FIX: Handle project_name (str -> str OR None -> ro.NULL)
    # The error happened because rpy2 cannot convert Python 'None' to R 'NULL' automatically.
    if project_name is None:
        r_project_name = ro.NULL
    else:
        r_project_name = project_name

    r_reference_set = reference_set

    # ---------------------------------------------------------
    # Define R function
    # ---------------------------------------------------------
    ro.r("""
        run_dsigd_orra <- function(
          interestGene,
          referenceGene,
          referenceSet,
          gmt_file,
          output_dir,
          project_name,
          fdr_thr
        ) {
          if (!is.null(referenceGene)) {
            WebGestaltR(
              enrichMethod       = "ORA",
              organism           = "hsapiens",
              enrichDatabase     = "others",
              enrichDatabaseFile = gmt_file,
              enrichDatabaseType = "genesymbol",
              interestGene       = interestGene,
              interestGeneType   = "genesymbol",
              referenceGene      = referenceGene,
              referenceGeneType  = "genesymbol",
              sigMethod          = "fdr",
              fdrMethod          = "BH",
              fdrThr             = fdr_thr,
              minNum             = 1,
              maxNum             = 5000,
              isOutput           = TRUE,
              outputDirectory    = output_dir,
              projectName        = project_name
            )
          } else {
            WebGestaltR(
              enrichMethod       = "ORA",
              organism           = "hsapiens",
              enrichDatabase     = "others",
              enrichDatabaseFile = gmt_file,
              enrichDatabaseType = "genesymbol",
              interestGene       = interestGene,
              interestGeneType   = "genesymbol",
              referenceSet       = referenceSet,
              referenceGeneType  = "genesymbol",
              sigMethod          = "fdr",
              fdrMethod          = "BH",
              fdrThr             = fdr_thr,
              minNum             = 1,
              maxNum             = 5000,
              isOutput           = TRUE,
              outputDirectory    = output_dir,
              projectName        = project_name
            )
          }
        }
    """)

    # ---------------------------------------------------------
    # Execute
    # ---------------------------------------------------------
    ro.r("run_dsigd_orra")(
        r_genes,
        r_background,
        r_reference_set,
        gmt_file,
        output_dir,
        r_project_name,  # Use the safely converted variable
        fdr_thr,
    )
    
    print(f"✅ DSigDB ORA completed → {output_dir}")