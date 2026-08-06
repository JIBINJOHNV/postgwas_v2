import requests
import pandas as pd

def query_intact(gene_symbol):
    """
    Queries EMBL-EBI IntAct for verified molecular interactions.
    """
    # PSICQUIC is the standard REST standard for interaction data
    url = f"http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/{gene_symbol}"
    try:
        response = requests.get(url)
        # The raw text is in MITAB format (standard interaction format)
        lines = response.text.splitlines()
        parsed_data = []
        for line in lines:
            parts = line.split("\t")
            if len(parts) > 10:
                # MITAB columns are complex, extracting IDs usually at index 0 and 1
                id_a = parts[0].split(":")[-1] # uniprotkb:P04637 -> P04637
                id_b = parts[1].split(":")[-1]
                method = parts[6].split("(")[0].replace('psi-mi:"', '')
                
                parsed_data.append({
                    "Query": gene_symbol,
                    "Interactor_A": id_a,
                    "Interactor_B": id_b,
                    "Detection_Method": method
                })
        return pd.DataFrame(parsed_data)
    except Exception as e:
        print(f"IntAct Error: {e}")
        return pd.DataFrame()

# --- Usage ---
# IntAct works best one gene at a time via this API
df_intact = query_intact("BRCA1")
print(df_intact.head())