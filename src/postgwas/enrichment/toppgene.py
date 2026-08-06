import requests
import pandas as pd
from pathlib import Path



## --- ToppGene Enrichment Module --- ## https://toppgene.cchmc.org/ 
def get_entrez_ids(gene_symbols):
    """ Converts symbols to Entrez IDs via ToppGene Lookup. """
    url = "https://toppgene.cchmc.org/API/lookup"
    try:
        resp = requests.post(url, json={"Symbols": gene_symbols}).json()
        return [g["Entrez"] for g in resp.get("Genes", []) if "Entrez" in g]
    except Exception as e:
        print(f"Lookup failed: {e}")
        return []

def perform_enrichment(entrez_ids, categories):
    """ Performs enrichment and returns the raw JSON response. """
    url = "https://toppgene.cchmc.org/API/enrich"
    # Configure parameters: PValue 1.0 allows non-significant results
    cat_params = [{
        "Type": cat,
        "PValue": 1.0, 
        "MinGenes": 1,
        "MaxGenes": 2000,
        "MaxResults": 50,
        "Correction": "FDR"
    } for cat in categories]
    payload = {"Genes": entrez_ids, "Categories": cat_params}
    try:
        resp = requests.post(url, json=payload)
        resp.raise_for_status()
        return resp.json()
    except Exception as e:
        print(f"Enrichment failed: {e}")
        return None

def results_to_dataframe(json_results):
    """ Converts ToppGene JSON output to a clean Pandas DataFrame. """
    if not json_results or "Annotations" not in json_results:
        return pd.DataFrame()
    annotations = json_results["Annotations"]
    df = pd.DataFrame(annotations)
    # 1. Clean up the 'Genes' column
    # The API returns 'Genes' as a list of dictionaries. We want a string of symbols.
    if 'Genes' in df.columns:
        df['Genes'] = df['Genes'].apply(
            lambda x: ", ".join([g['Symbol'] for g in x]) if isinstance(x, list) else x
        )
    # 2. Select and Reorder useful columns
    # We prioritize these columns if they exist
    target_cols = ['Category', 'ID', 'Name', 'PValue', 'QValueFDRBH', 'QValueFDRBY', 'QValueBonferroni', 'TotalGenes', 
     'GenesInTerm', 'GenesInQuery', 'GenesInTermInQuery', 'Source', 'URL', 'Genes']
    # Filter to only columns that actually exist in the response
    final_cols = [c for c in target_cols if c in df.columns]
    return df[final_cols]


def run_topgene_enrichment(
    gene_list,
    output_dir=None,
    sample_id=None,
    save=True,
):
    """
    Run ToppGene enrichment analysis.
    Args:
        my_genes (list): List of gene symbols.
        output_dir (str | Path | None): Directory to save ToppGene results.
        sample_id (str | None): Sample/dataset identifier for filenames.
        save (bool): Whether to save output to disk.
    Returns:
        pd.DataFrame: ToppGene enrichment results.
    """
    ids = get_entrez_ids(gene_list)
    print(f"\t\t\t\t Tptal genes provided is {len(gene_list)} ")
    print(f"\t\t\t\t Tptal genes converted to  entrez id is {len(ids)} ")
    if not ids:
        print("❌ No Entrez IDs mapped. Skipping ToppGene enrichment.")
        return pd.DataFrame()
    categories = [
        "GeneOntologyBiologicalProcess",
        "GeneOntologyMolecularFunction",
        "GeneOntologyCellularComponent",
        "Pathway",
        "HumanPheno",
        "MousePheno",
        "Disease",
        "Drug",
        "Domain",
        "Interaction",
        "MicroRNA",
        "TFBS",
        "Pubmed",
        "Cytoband",
        "GeneFamily",
        "Coexpression",
        "CoexpressionAtlas",
        "ToppCell",
        "Computational",
    ]
    try:
        print(f"🧬 Running ToppGene enrichment for {len(ids)} genes…")
        raw_json = perform_enrichment(ids, categories)
        df_results = results_to_dataframe(raw_json)
        if df_results.empty:
            print("ℹ️ No ToppGene enrichment results returned.")
            return df_results
        # ---------------------------------------------------------
        # Sort results
        # ---------------------------------------------------------
        if {"Category", "PValue"}.issubset(df_results.columns):
            df_results = df_results.sort_values(
                by=["Category", "PValue"]
            )
        # ---------------------------------------------------------
        # Save output (optional)
        # ---------------------------------------------------------
        if save and output_dir is not None:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            if sample_id is None:
                sample_id = "toppgene"
            out_file = output_dir / f"{sample_id}_toppgene_enrichment.tsv"
            df_results.to_csv(out_file, sep="\t", index=False)
            print(f"💾 ToppGene results saved to: {out_file}")
        return df_results
    except Exception as e:
        print(f"❌ ToppGene enrichment failed: {e}")
        return pd.DataFrame()




# --- Main Execution ---
if __name__ == "__main__":
    # 1. Input Data
    my_genes = ["TP53", "BRCA1", "EGFR", "TNF", "APOE"]
    
    df_results=run_topgene_enrichment(my_genes)
    if not df_results.empty:
        # Sort by Category and then Significance (P-value)
        df_results = df_results.sort_values(by=["Category", "PValue"])
        # Display first few rows
        print("Top 5 Results:")
        print(df_results.head().to_string(index=False))
        # OPTIONAL: Save to CSV
        # df_results.to_csv("toppgene_results.csv", index=False)
        # print("\nSaved results to 'toppgene_results.csv'")
    else:
        print("No results returned.")
            
