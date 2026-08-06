import pandas as pd
from pathlib import Path
import zeep
import time



def david_symbols_to_entrez(genes, email):
    """
    Convert gene symbols to ENTREZ IDs using DAVID SOAP API only.
    Args:
        genes (list[str]): Gene symbols
        email (str): DAVID registered email
    Returns:
        list[str]: ENTREZ Gene IDs
    """
    url = "https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService?wsdl"
    client = zeep.Client(url)
    if client.service.authenticate(email) != "true":
        raise RuntimeError("DAVID authentication failed")
    # Normalize
    gene_str = ",".join(genes)
    list_name = f"convert_{int(time.time())}"
    # Upload symbols
    mapped = client.service.addList(
        gene_str,
        "GENE_SYMBOL",
        list_name,
        0
    )
    if mapped < 1:
        print("⚠️ DAVID failed to map any symbols")
        return []
    # Get ENTREZ IDs
    entrez_ids = client.service.getGeneIds()
    return list(map(str, entrez_ids))



def run_david_enrichment(
    gene_list,
    email,
    id_type="OFFICIAL_GENE_SYMBOL",
    output_dir=None,
    sample_id=None,
    save=True,
):
    url = "https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService?wsdl"
    try:
        print("🔗 Connecting to DAVID Web Service…")
        client = zeep.Client(url)
        # ---------------------------------------------------------
        # Authenticate
        # ---------------------------------------------------------
        if client.service.authenticate(email) != "true":
            print(
                f"❌ DAVID authentication failed for {email}. "
                "Register at https://david.ncifcrf.gov/webservice/register.html"
            )
            return pd.DataFrame()
        print("✅ Authentication successful.")
        # ---------------------------------------------------------
        # Normalize ID type (SOAP-safe)
        # ---------------------------------------------------------
        id_type = id_type.upper()
        if id_type == "OFFICIAL_GENE_SYMBOL":
            id_type = "GENE_SYMBOL"
        # ---------------------------------------------------------
        # Force species (important)
        # ---------------------------------------------------------
        try:
            client.service.setSpecies("Homo sapiens")
        except Exception:
            pass  # older DAVID backends ignore this
        # ---------------------------------------------------------
        # Upload gene list (unique name avoids state conflicts)
        # ---------------------------------------------------------
        genes_str = ",".join(map(str, gene_list))
        list_name = f"Python_List_{int(time.time())}"
        mapped = client.service.addList(
            genes_str,
            id_type,
            list_name,
            0,
        )
        if mapped < 1:
            print(
                "⚠️ DAVID mapped 0 genes.\n"
                "👉 Strongly recommended: convert symbols → ENTREZ IDs\n"
                "👉 Then rerun with id_type='ENTREZ_GENE_ID'"
            )
            return pd.DataFrame()
        print(f"📥 Uploaded {mapped} mapped genes.")
        # ---------------------------------------------------------
        # Fetch enrichment
        # ---------------------------------------------------------
        print("📊 Fetching DAVID enrichment report…")
        records = client.service.getChartReport(0.1, 2)
        if not records:
            print("ℹ️ No significant enrichment found.")
            return pd.DataFrame()
        df = pd.DataFrame(
            [
                {
                    "Category": r["categoryName"],
                    "Term": r["termName"],
                    "Count": r["listHits"],
                    "Percent": r["percent"],
                    "PValue": r["ease"],
                    "Bonferroni": r["bonferroni"],
                    "Benjamini": r["benjamini"],
                    "FDR": r["fdr"],
                    "Genes": r["geneIds"],
                }
                for r in records
            ]
        )
        # ---------------------------------------------------------
        # Save output
        # ---------------------------------------------------------
        if save and output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            sample_id = sample_id or "david"
            out = output_dir / f"{sample_id}_david_enrichment.tsv"
            df.to_csv(out, sep="\t", index=False)
            print(f"💾 DAVID enrichment saved to: {out}")
        return df
    except Exception as e:
        print(f"❌ DAVID SOAP error: {e}")
        return pd.DataFrame()



# --- Main Execution ---
if __name__ == "__main__":
    # 1. Your Gene List
    my_genes = ["TP53", "BRCA1", "EGFR", "TNF", "APOE", "IL6"]
    
    # 2. YOUR REGISTERED EMAIL
    # You MUST replace this with your actual registered email
    email_id = "jjohn41@northwell.edu" 
    
    # 3. Run Analysis
    df = run_david_enrichment( gene_list=my_genes,
            email=email_id,
            id_type="GENE_SYMBOL",
            output_dir=None,
            sample_id=None,
            save=True)
    
    if not df.empty:
        # Sort by P-Value
        df = df.sort_values(by="P-Value")
        
        print("\nTop 5 Enriched Terms:")
        print(df[['Category', 'Term', 'P-Value', 'Benjamini']].head(5).to_string(index=False))
        
        # df.to_csv("david_results.csv", index=False)
        
        