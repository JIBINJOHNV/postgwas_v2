import py4cytoscape as p4c
import pandas as pd

def run_genemania_in_cytoscape(gene_list, organism="Homo sapiens"):
    """
    Controls a running instance of Cytoscape to perform GeneMANIA analysis.
    Requires Cytoscape to be OPEN on your desktop.
    """
    # 1. Check connection to Cytoscape
    try:
        p4c.cytoscape_ping()
        print("Successfully connected to Cytoscape!")
    except Exception:
        print("Error: Could not connect. Make sure Cytoscape is OPEN.")
        return None
    # 2. Format genes for GeneMANIA (pipe-separated)
    genes_query = "|".join(gene_list)
    print(f"Querying GeneMANIA for {len(gene_list)} genes...")
    # 3. Run the GeneMANIA Command via CyREST
    # This command triggers the search inside the GUI
    # limit: Number of related genes to find (e.g., 20)
    try:
        cmd = f'genemania search organism="{organism}" genes="{genes_query}" limit=20'
        res = p4c.commands.run(cmd)
        print("GeneMANIA search complete. Network created in Cytoscape.")
        # 4. Extract the Node Table (The Genes & Scores) back to Python
        # This gets the data from the newly created network window
        nodes = p4c.tables.get_table_columns(table='node')
        # Filter columns to keep it clean
        # GeneMANIA usually adds columns like 'Gene Name', 'Score', 'Description'
        if not nodes.empty:
            cols = [c for c in nodes.columns if c in ['name', 'Gene Name', 'Score', 'description']]
            return nodes[cols]
        return pd.DataFrame()
    except Exception as e:
        print(f"GeneMANIA Command Failed: {e}")
        return None

# --- Main Execution ---
if __name__ == "__main__":
    my_genes = ["TP53", "BRCA1", "EGFR", "TNF"]
    
    # Cytoscape must be running for this to work!
    df_results = run_genemania_in_cytoscape(my_genes)
    
    if df_results is not None and not df_results.empty:
        print("\nTop Genes in Network:")
        print(df_results.head(10).to_string(index=False))
        
        # Optional: Export Network Image
        # p4c.export_image("genemania_network.png")