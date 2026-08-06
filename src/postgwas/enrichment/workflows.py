from postgwas.enrichment.biogrid import fetch_biogrid_interactions 
from postgwas.enrichment.david import run_david_enrichment
from postgwas.enrichment.dgidb import query_dgidb
from postgwas.enrichment.enricher import run_enrichr
from postgwas.enrichment.gprofiler2 import run_gprofiler_with_raw_stats
from postgwas.enrichment.omnipath_wrapper import run_omnipath_interactions
from postgwas.enrichment.toppgene import run_topgene_enrichment
from postgwas.enrichment.stringdb import run_stringdb_analysis
from postgwas.enrichment.DSigDB import run_dsigd_ora_webgestalt
import rpy2.robjects as ro
from rpy2.robjects import StrVector
from pathlib import Path
import traceback




# 1. SETUP: You MUST paste your key here
# Register at: https://webservice.thebiogrid.org/

# output_dir="/Users/JJOHN41/Documents/developing_software/data/oudir/PGC3_SCZ_european/enrichment/"
# biogrid_access_key = "86cf1c2c7c8b2972cc5e6b31c45e8584"
# sample_id='PGC3_SCZ_european'
# my_genes = ["BRCA1", "TP53", "EGFR", "MYC"]
# email_id="johnjibinv@gmail.com"

def run_multisource_enrichment_pipeline(
    gene_list,
    output_dir,
    sample_id,
    biogrid_access_key,
    david_email,
    dsigdb_gmt,
    score_threshold=400,
    fdr_thr=0.999,
):
    """
    Run a comprehensive multi-database enrichment & interaction pipeline.
    Each module runs in isolation with its own output folder.
    """

    base_out = Path(output_dir)
    base_out.mkdir(parents=True, exist_ok=True)

    print(f"\n🧬 Running enrichment pipeline for: {sample_id}")
    print(f"🧬 Total input genes: {len(gene_list)}")


    # ---------------------------------------------------------
    # 7. OmniPath
    # ---------------------------------------------------------
    try:
        print("\n🔹 OmniPath")
        out = base_out / "omnipath"
        out.mkdir(exist_ok=True)

        run_omnipath_interactions(
            gene_list=gene_list,
            output_dir=out,
            sample_id=sample_id,
            save=True,
        )

    except Exception:
        print("❌ OmniPath failed")
        traceback.print_exc()
    # ---------------------------------------------------------
    # 1. DSigDB ORA (WebGestaltR)
    # ---------------------------------------------------------
    try:
        print("\n🔹 DSigDB ORA (WebGestaltR)")
        out = base_out / "dsigdb"
        out.mkdir(exist_ok=True)

        run_dsigd_ora_webgestalt(
            genes=gene_list,
            gmt_file=dsigdb_gmt,
            output_dir=out,
            background_genes=None,
            reference_set="genome_protein-coding",
            project_name=None,  # avoid zip bug
            fdr_thr=fdr_thr,
        )
    except Exception:
        print("❌ DSigDB ORA failed")
        traceback.print_exc()

    # ---------------------------------------------------------
    # 2. BioGRID
    # ---------------------------------------------------------
    try:
        print("\n🔹 BioGRID interactions")
        out = base_out / "biogrid"
        out.mkdir(exist_ok=True)

        df_biogrid = fetch_biogrid_interactions(
            gene_list=gene_list,
            access_key=biogrid_access_key,
            tax_id=9606,
            output_dir=out,
            sample_id=sample_id,
        )

        if not df_biogrid.empty:
            df_phys = df_biogrid[
                df_biogrid["Experimental System Type"] == "physical"
            ]
            print(f"  Total: {len(df_biogrid)} | Physical: {len(df_phys)}")
        else:
            print("  No BioGRID interactions found")

    except Exception:
        print("❌ BioGRID failed")
        traceback.print_exc()

    # ---------------------------------------------------------
    # 3. DAVID
    # ---------------------------------------------------------
    try:
        print("\n🔹 DAVID enrichment")
        out = base_out / "david"
        out.mkdir(exist_ok=True)

        df_david = run_david_enrichment(
            gene_list=gene_list,
            email=david_email,
            id_type="OFFICIAL_GENE_SYMBOL",
            output_dir=out,
            sample_id=sample_id,
            save=True,
        )

        if not df_david.empty:
            print(f"  DAVID terms: {len(df_david)}")
        else:
            print("  DAVID returned no results")

    except Exception:
        print("❌ DAVID failed")
        traceback.print_exc()

    # ---------------------------------------------------------
    # 4. DGIdb
    # ---------------------------------------------------------
    try:
        print("\n🔹 DGIdb drug–gene interactions")
        out = base_out / "dgidb"
        out.mkdir(exist_ok=True)

        df_dgidb = query_dgidb(
            gene_list=gene_list,
            output_dir=out,
            sample_id=sample_id,
            save=True,
        )

        if not df_dgidb.empty:
            print(f"  DGIdb interactions: {len(df_dgidb)}")
        else:
            print("  No DGIdb interactions")

    except Exception:
        print("❌ DGIdb failed")
        traceback.print_exc()

    # ---------------------------------------------------------
    # 5. Enrichr
    # ---------------------------------------------------------
    try:
        print("\n🔹 Enrichr")
        out = base_out / "enrichr"
        out.mkdir(exist_ok=True)

        df_enrichr = run_enrichr(
            gene_list=gene_list,
            output_dir=out,
            sample_id=sample_id,
            save=True,
        )

        if not df_enrichr.empty:
            sig = df_enrichr[df_enrichr["AdjPValue"] <= 0.05]
            print(f"  Significant: {sig.shape[0]}")

    except Exception:
        print("❌ Enrichr failed")
        traceback.print_exc()

    # ---------------------------------------------------------
    # 6. gProfiler
    # ---------------------------------------------------------
    try:
        print("\n🔹 gProfiler2")
        out = base_out / "gprofiler"
        out.mkdir(exist_ok=True)

        df_gprof = run_gprofiler_with_raw_stats(
            gene_list=gene_list,
            output_dir=out,
            sample_id=sample_id,
            save=True,
        )

        if not df_gprof.empty:
            sig = df_gprof[df_gprof["p_value"] <= 0.05]
            print(f"  Significant: {sig.shape[0]}")

    except Exception:
        print("❌ gProfiler failed")
        traceback.print_exc()

    # ---------------------------------------------------------
    # 8. ToppGene
    # ---------------------------------------------------------
    try:
        print("\n🔹 ToppGene")
        out = base_out / "toppgene"
        out.mkdir(exist_ok=True)

        run_topgene_enrichment(
            gene_list=gene_list,
            output_dir=out,
            sample_id=sample_id,
            save=True,
        )

    except Exception:
        print("❌ ToppGene failed")
        traceback.print_exc()

    # ---------------------------------------------------------
    # 9. STRINGdb
    # ---------------------------------------------------------
    try:
        print("\n🔹 STRINGdb")
        out = base_out / "stringdb"
        out.mkdir(exist_ok=True)

        run_stringdb_analysis(
            gene_list=gene_list,
            output_dir=out,
            sample_id=sample_id,
            score_threshold=score_threshold,
            plot_network=True,
        )

    except Exception:
        print("❌ STRINGdb failed")
        traceback.print_exc()

    print(f"\n✅ Pipeline completed for {sample_id}")



# run_multisource_enrichment_pipeline(
#     gene_list=["BRCA1", "TP53", "EGFR", "MYC"],
#     output_dir="/Users/JJOHN41/Documents/developing_software/data/oudir/PGC3_SCZ_european/enrichment/",
#     sample_id="PGC3_SCZ_european",
#     biogrid_access_key="YOUR_BIOGRID_KEY",
#     david_email="johnjibinv@gmail.com",
#     dsigdb_gmt="/path/to/DSigDB_All.gmt",
# )
