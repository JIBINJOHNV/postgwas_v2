import requests
import pandas as pd
from pathlib import Path


def query_dgidb(
    gene_list,
    output_dir=None,
    sample_id=None,
    save=True,
):
    """
    Queries DGIdb v5 (GraphQL) for drug–gene interactions.
    Args:
        gene_list (list): List of gene symbols.
        output_dir (str | Path | None): Directory to save DGIdb results.
        sample_id (str | None): Sample/dataset identifier for filenames.
        save (bool): Whether to save output to disk.
    Returns:
        pd.DataFrame: Drug–gene interaction table.
    """
    url = "https://dgidb.org/api/graphql"
    query = """
    query getInteractions($genes: [String!]!) {
      genes(names: $genes) {
        nodes {
          name
          interactions {
            interactionScore
            interactionTypes {
              type
              directionality
            }
            drug {
              name
              conceptId
              approved
              drugApprovalRatings {
                rating
              }
            }
            sources {
              sourceDbName
              pmid
            }
          }
        }
      }
    }
    """
    variables = {"genes": gene_list}
    try:
        print(f"💊 Querying DGIdb for {len(gene_list)} genes…")
        response = requests.post(
            url,
            json={"query": query, "variables": variables},
            timeout=60,
        )
        response.raise_for_status()
        data = response.json()
        rows = []
        if data.get("data") and data["data"].get("genes"):
            for gene_node in data["data"]["genes"]["nodes"]:
                gene_symbol = gene_node["name"]
                for interaction in gene_node.get("interactions", []):
                    drug = interaction.get("drug", {})
                    drug_name = drug.get("name")
                    approved = drug.get("approved")
                    score = interaction.get("interactionScore")
                    types = [
                        t["type"]
                        for t in interaction.get("interactionTypes", [])
                        if t.get("type")
                    ]
                    type_str = ", ".join(types) if types else "n/a"
                    sources = [
                        s["sourceDbName"]
                        for s in interaction.get("sources", [])
                        if s.get("sourceDbName")
                    ]
                    source_str = ", ".join(sorted(set(sources)))
                    pmids = [
                        str(s["pmid"])
                        for s in interaction.get("sources", [])
                        if s.get("pmid")
                    ]
                    pmid_str = ", ".join(sorted(set(pmids)))
                    rows.append({
                        "Gene": gene_symbol,
                        "Drug": drug_name,
                        "InteractionType": type_str,
                        "Approved": approved,
                        "Score": score,
                        "Sources": source_str,
                        "PMIDs": pmid_str,
                    })
        df = pd.DataFrame(rows)
        # ---------------------------------------------------------
        # Save output (optional)
        # ---------------------------------------------------------
        if save and output_dir is not None:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            if sample_id is None:
                sample_id = "dgidb"
            out_file = output_dir / f"{sample_id}_dgidb_interactions.tsv"
            df.to_csv(out_file, sep="\t", index=False)
            print(f"💾 DGIdb interactions saved to: {out_file}")
        return df
    except Exception as e:
        print(f"❌ DGIdb request failed: {e}")
        return pd.DataFrame()


# --- Main Execution ---
if __name__ == "__main__":
    # 1. Input Genes
    my_genes = ["EGFR", "BRAF", "TP53", "TNF", "BRCA1"]
    # 2. Run Query
    df = query_dgidb(gene_list=my_genes, output_dir=None,sample_id=None, save=True,)
    if not df.empty:
        # 3. Filter Results
        # Example: Show only FDA Approved drugs that are "inhibitors"
        approved_inhibitors = df[
            (df['Approved'] == True) & 
            (df['Interaction Type'].str.contains("inhibitor", case=False))
        ]
        print(f"\nTotal Interactions Found: {len(df)}")
        print(f"FDA Approved Inhibitors: {len(approved_inhibitors)}\n")
        # Display Top 10
        display_cols = ['Gene', 'Drug', 'Interaction Type', 'Score', 'Sources']
        print(approved_inhibitors[display_cols].head(10).to_string(index=False))
        # 4. Save to CSV
        # df.to_csv("dgidb_results.csv", index=False)
        # print("\nSaved to dgidb_results.csv")
    else:
        print("No interactions found or API error.")