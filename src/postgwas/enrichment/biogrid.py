import os
import io
import requests
import pandas as pd
from pathlib import Path



def fetch_biogrid_interactions(
    gene_list,
    access_key,
    tax_id=9606,
    output_dir=None,
    sample_id=None,
    save=True,
):
    """
    Fetch physical interactions from BioGRID for a list of genes.
    Args:
        gene_list (list): List of gene symbols (e.g., ['TP53', 'BRCA1'])
        access_key (str): BioGRID API key
        tax_id (int): 9606 = Human, 10090 = Mouse
        output_dir (str | Path | None): Directory to save output TSV
        sample_id (str | None): Sample or dataset identifier used in filename
        save (bool): Whether to save results to disk

    Returns:
        pd.DataFrame: BioGRID interaction table
    """
    base_url = "https://webservice.thebiogrid.org/interactions/"
    genes_query = "|".join(gene_list)
    params = {
        "accesskey": access_key,
        "geneList": genes_query,
        "taxId": tax_id,
        "searchNames": "true",
        "includeHeader": "true",
        "throughputTag": "any",
        "format": "tab2",
    }
    try:
        print(f"🔍 Querying BioGRID for {len(gene_list)} genes…")
        response = requests.get(base_url, params=params, timeout=60)
        response.raise_for_status()
        df = pd.read_csv(io.StringIO(response.text), sep="\t")
        # ---------------------------------------------------------
        # Save output (optional)
        # ---------------------------------------------------------
        if save and output_dir is not None:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            if sample_id is None:
                sample_id = "biogrid"
            out_file = output_dir / f"{sample_id}_biogrid_interactions.tsv"
            df.to_csv(out_file, sep="\t", index=False)
            print(f"💾 BioGRID interactions saved to: {out_file}")
        return df
    except Exception as e:
        print(f"❌ BioGRID request failed: {e}")
        return pd.DataFrame()


# # --- Main Execution ---
# if __name__ == "__main__":
#     # 1. SETUP: You MUST paste your key here
#     # Register at: https://webservice.thebiogrid.org/
#     MY_ACCESS_KEY = "86cf1c2c7c8b2972cc5e6b31c45e8584"
    
#     my_genes = ["BRCA1", "TP53", "EGFR", "MYC"]

#     if MY_ACCESS_KEY == "PASTE_YOUR_ACCESS_KEY_HERE":
#         print("[!] Error: You must replace MY_ACCESS_KEY with a valid BioGRID key.")
#     else:
#         # 2. Fetch Data
#         df_biogrid = fetch_biogrid_interactions(my_genes, MY_ACCESS_KEY)

#         if not df_biogrid.empty:
#             # 3. Filter for PHYSICAL interactions only
#             # BioGRID contains 'genetic' interactions (e.g. synthetic lethality), 
#             # which you might want to remove if you only care about binding.
#             df_phys = df_biogrid[df_biogrid['Experimental System Type'] == 'physical']
            
#             print(f"\nTotal Interactions: {len(df_biogrid)}")
#             print(f"Physical Interactions: {len(df_phys)}\n")
            
#             # 4. Clean up columns for display
#             # BioGRID returns A LOT of columns. Let's select the useful ones.
#             cols = [
#                 'Official Symbol Interactor A', 
#                 'Official Symbol Interactor B', 
#                 'Experimental System Name', 
#                 'Publication Source',  # Pubmed ID or Author
#                 'Throughput'
#             ]
            
#             # Display Top 10
#             print(df_phys[cols].head(10).to_string(index=False))
            
#             # OPTIONAL: Save to CSV
#             # df_phys.to_csv("biogrid_physical_interactions.csv", index=False)
#         else:
#             print("No interactions found.")