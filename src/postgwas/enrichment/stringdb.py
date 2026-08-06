import requests
import pandas as pd
import networkx as nx
from pathlib import Path
import matplotlib.pyplot as plt


__all__ = [
    "run_stringdb_analysis",
    "get_string_ids",
    "run_string_enrichment",
    "get_string_ppi_network",
    "get_network_image",
    "analyze_hubs",
    "graph_to_hub_dataframe",
]


def graph_to_hub_dataframe(G: nx.Graph) -> pd.DataFrame:
    """Convert a NetworkX graph with node attributes into a DataFrame."""
    rows = []
    for node, attrs in G.nodes(data=True):
        row = {"Gene": node}
        row.update(attrs)  # degree, betweenness, closeness, etc.
        rows.append(row)
    return pd.DataFrame(rows)


def get_string_ids(gene_list, species: int = 9606):
    """Map gene symbols to STRING protein IDs."""
    url = "https://string-db.org/api/json/get_string_ids"
    params = {
        "identifiers": "\r".join(gene_list),
        "species": species,
        "limit": 1,
        "echo_query": 1,
    }
    try:
        response = requests.post(url, data=params, timeout=60)
        response.raise_for_status()
        data = response.json()
        mapped_ids = [item["stringId"] for item in data if "stringId" in item]
        return sorted(set(mapped_ids))
    except Exception as e:
        print(f"❌ STRING ID mapping failed: {e}")
        return []


def run_string_enrichment(string_ids, species: int = 9606) -> pd.DataFrame:
    """Perform functional enrichment on mapped STRING IDs."""
    url = "https://string-db.org/api/json/enrichment"
    params = {"identifiers": "\r".join(string_ids), "species": species}
    try:
        response = requests.post(url, data=params, timeout=60)
        response.raise_for_status()
        return pd.DataFrame(response.json())
    except Exception as e:
        print(f"❌ STRING enrichment failed: {e}")
        return pd.DataFrame()


def get_network_image(
    string_ids,
    species: int = 9606,
    output_dir=None,
    sample_id=None,
    filename=None,
    image_format: str = "svg",   # "svg" (best) or "png"
    overwrite: bool = False,
):
    """
    Download a static STRING network image (SVG recommended).
    Returns Path to saved image, or None on failure.
    """
    fmt = image_format.lower()
    if fmt not in {"svg", "png"}:
        raise ValueError("image_format must be 'svg' or 'png'")

    # STRING endpoints:
    # PNG -> /api/image/network
    # SVG -> /api/svg/network
    endpoint = "svg" if fmt == "svg" else "image"
    url = f"https://string-db.org/api/{endpoint}/network"

    params = {
        "identifiers": "\r".join(string_ids),
        "species": species,
        "add_white_nodes": 0,
        "network_flavor": "evidence",
        "hide_disconnected_nodes": 1,
    }

    try:
        #print(f"🖼️ Downloading STRING network image ({fmt.upper()})…")
        response = requests.post(url, data=params, timeout=60)
        response.raise_for_status()

        output_dir = Path(output_dir) if output_dir else Path(".")
        output_dir.mkdir(parents=True, exist_ok=True)

        if filename:
            out_path = output_dir / filename
        else:
            base = sample_id or "string_network"
            out_path = output_dir / f"{base}_string_network.{fmt}"

        # avoid overwriting unless explicitly allowed
        if out_path.exists() and not overwrite:
            stem, suffix = out_path.stem, out_path.suffix
            i = 1
            while True:
                candidate = output_dir / f"{stem}_v{i}{suffix}"
                if not candidate.exists():
                    out_path = candidate
                    break
                i += 1

        with open(out_path, "wb") as f:
            f.write(response.content)

        #print(f"💾 Network image saved to: {out_path}")
        return out_path

    except Exception as e:
        print(f"\t\t\t\t❌ STRING image download failed: {e}")
        return None


def get_string_ppi_network(gene_list, species: int = 9606, score_threshold: int = 400) -> pd.DataFrame:
    """Fetch PPI interactions from STRING (edge list)."""
    url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": "%0d".join(gene_list),
        "species": species,
        "required_score": score_threshold,
    }
    try:
        #print(f"\t\t\t\tFetching PPI network for {len(gene_list)} genes...")
        response = requests.post(url, data=params, timeout=60)
        response.raise_for_status()
        data = response.json()
        if not data:
            return pd.DataFrame()

        df = pd.DataFrame(data)

        needed = [
            "preferredName_A", "preferredName_B", "score",
            "nscore", "fscore", "pscore", "ascore", "escore", "dscore", "tscore"
        ]
        df_clean = df[needed].copy()

        df_clean = df_clean.rename(columns={
            "preferredName_A": "Source",
            "preferredName_B": "Target",
            "score": "Total_Score",
            "nscore": "neighborhood_score",
            "fscore": "fusion_score",
            "pscore": "cooccurence_score",
            "ascore": "coexpression_score",
            "escore": "experimental_score",
            "dscore": "database_score",
            "tscore": "textmining_score",
        })
        return df_clean

    except Exception as e:
        print(f"❌ Error fetching STRING PPI: {e}")
        return pd.DataFrame()


def analyze_hubs(df_ppi: pd.DataFrame) -> nx.Graph:
    """Build graph and compute node metrics."""
    G = nx.from_pandas_edgelist(df_ppi, source="Source", target="Target", edge_attr=True)
    nx.set_node_attributes(G, dict(G.degree()), "degree")
    nx.set_node_attributes(G, nx.betweenness_centrality(G), "betweenness")
    nx.set_node_attributes(G, nx.closeness_centrality(G), "closeness")
    return G


def run_stringdb_analysis(
    gene_list,
    output_dir=None,
    sample_id=None,
    score_threshold: int = 400,
    plot_network: bool = True,
):
    """
    Run STRING enrichment + PPI + hubs + save outputs.
    Returns dict with dataframes/graph/paths.
    """
    output_dir = Path(output_dir) if output_dir else Path(".")
    output_dir.mkdir(parents=True, exist_ok=True)
    sample_id = sample_id or "string"

    #print(f"\t\t\t\t🧬 Total Input genes: {len(gene_list)}")

    results = {
        "string_ids": [],
        "enrichment": None,
        "ppi": None,
        "graph": None,
        "network_image": None,
        "network_plot": None,
        "hub_table": None,
    }

    # 1) IDs
    mapped_ids = get_string_ids(gene_list)
    if not mapped_ids:
        print("❌ No genes could be mapped to STRING IDs.")
        return results

    #print(f"✅ Successfully mapped {len(mapped_ids)} proteins.")
    results["string_ids"] = mapped_ids

    # 2) Enrichment
    df_enrich = run_string_enrichment(mapped_ids)
    if isinstance(df_enrich, pd.DataFrame) and not df_enrich.empty:
        # fdr is usually present; guard just in case
        if "fdr" in df_enrich.columns:
            df_enrich = df_enrich.sort_values(by="fdr")
        enrich_file = output_dir / f"{sample_id}_stringdb_enrichment.tsv"
        df_enrich.to_csv(enrich_file, sep="\t", index=False)
        #print(f"💾 Enrichment saved to: {enrich_file}")
        results["enrichment"] = df_enrich
    else:
        print("\t\t\t\tℹ️ No enrichment results returned.")

    # 3) STRING network image (best quality: SVG)
    img_path = get_network_image(
        mapped_ids,
        output_dir=output_dir,
        sample_id=sample_id,
        image_format="svg",
    )
    results["network_image"] = img_path

    # 4) PPI
    df_ppi = get_string_ppi_network(gene_list, score_threshold=score_threshold)
    if df_ppi.empty:
        #print("ℹ️ No STRING interactions returned.")
        return results

    ppi_file = output_dir / f"{sample_id}_stringdb_interactions.tsv"
    df_ppi.to_csv(ppi_file, sep="\t", index=False)
    #print(f"💾 PPI network saved to: {ppi_file}")
    results["ppi"] = df_ppi

    # 5) Hubs
    G = analyze_hubs(df_ppi)
    results["graph"] = G

    df_hubs = graph_to_hub_dataframe(G)
    if "degree" in df_hubs.columns:
        df_hubs = df_hubs.sort_values(by="degree", ascending=False)

    hub_file = output_dir / f"{sample_id}_stringdb_ppi_hub_analysis.tsv"
    df_hubs.to_csv(hub_file, sep="\t", index=False)
    #print(f"💾 Hub analysis saved to: {hub_file}")
    results["hub_table"] = hub_file

    # 6) Plot (save as SVG)
    if plot_network and G.number_of_nodes() > 0:
        plt.figure(figsize=(8, 8))
        pos = nx.spring_layout(G, seed=42)
        nx.draw(G, pos, with_labels=True, node_size=1500, font_weight="bold")
        plt.title("STRING Protein–Protein Interaction Network")

        fig_path = output_dir / f"{sample_id}_stringdb_ppi_network.svg"
        plt.savefig(fig_path, format="svg", bbox_inches="tight")
        plt.close()

        #print(f"💾 Network plot saved to: {fig_path}")
        results["network_plot"] = fig_path

    return results















# --- Main Execution ---
if __name__ == "__main__":
    # 1. Input Genes
    my_genes = ["TP53", "BRCA1", "EGFR", "TNF", "AKT1", "MYC"]
    print(f"Input: {my_genes}")

    # 2. Map to STRING IDs
    mapped_ids = get_string_ids(my_genes)
    
    if mapped_ids:
        print(f"Successfully mapped {len(mapped_ids)} proteins.")
        
        # 3. Run Enrichment
        df = run_string_enrichment(mapped_ids)
        
        if not df.empty:
            # Clean and Sort Results
            # Columns usually: category, term, number_of_genes, fdr, description
            df_clean = df[['category', 'term', 'description', 'fdr', 'number_of_genes']]
            df_clean = df_clean.sort_values(by='fdr')
            
            print("\nTop 5 Enriched Terms:")
            print(df_clean.head().to_string(index=False))
            
            # Save CSV
            # df_clean.to_csv("string_enrichment_results.csv", index=False)
        
        # 4. Get Network Image (The cool part of STRING)
        #get_network_image(mapped_ids,output_dir=output_dir,sample_id=sample_id)
    else:
        print("No genes could be mapped to STRING IDs.")

    # 2. Get Interactions
    # Score 400 is "Medium Confidence". Increase to 700 for "High Confidence".
    df_ppi = get_string_ppi_network(my_genes, score_threshold=400)
    
    if not df_ppi.empty:
        print(f"\nFound {len(df_ppi)} interactions.")
        
        # Show data
        print(df_ppi[['Source', 'Target', 'Total_Score']].head().to_string(index=False))
        
        # 3. Find Hubs
        G = analyze_hubs(df_ppi)
        
        # 4. Save to CSV for Cytoscape (Optional)
        # df_ppi.to_csv("string_interactions.csv", index=False)
        # print("\nSaved interactions to 'string_interactions.csv' (Load this into Cytoscape!)")
        
        # 5. Simple Visualization (Optional)
        plt.figure(figsize=(8, 8))
        pos = nx.spring_layout(G, seed=42) # Consistent layout
        nx.draw(G, pos, with_labels=True, node_color='skyblue', node_size=1500, font_weight='bold')
        plt.title("Protein-Protein Interaction Network")
        plt.show()
    else:
        print("No interactions returned.")