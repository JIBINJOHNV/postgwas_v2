import pandas as pd
from pathlib import Path
import omnipath as op

# https://omnipathdb.org/ ; https://omnipath.readthedocs.io/en/latest/api.html 
# 


import pandas as pd
from pathlib import Path
import os
import omnipath as op
import omnipath.interactions 

def run_omnipath_interactions(
    gene_list,
    output_dir=None,
    sample_id=None,
    save=True,
):
    """
    Fetches interactions using the official OmniPath client.
    """
    
    # ---------------------------------------------------------
    # ✅ FIX: Explicitly set the cache directory
    # ---------------------------------------------------------
    # We try to use the safe directory defined in main.py (via XDG_CACHE_HOME),
    # but if that is missing, we create a safe local temp folder.
    try:
        # Get the path (defaulting to a safe /tmp/user_ID folder if env var is missing)
        uid = os.getuid()
        safe_cache = os.environ.get("XDG_CACHE_HOME", f"/tmp/user_{uid}_omnipath_cache")
        
        # Ensure it exists
        Path(safe_cache).mkdir(parents=True, exist_ok=True)
        
        # Apply the setting found in documentation
        op.options.cache = safe_cache
        print(f"⚙️  OmniPath cache set to: {op.options.cache}")
        
    except Exception as e:
        print(f"⚠️  Could not configure OmniPath cache: {e}")

    try:
        print("🔗 Fetching OmniPath data via official client…")

        # Fetch full human interaction network
        df = op.interactions.OmniPath.get(
            genesymbols=True,
            organism="human",
        )

        # ---------------------------------------------------------
        # Filter for genes of interest
        # ---------------------------------------------------------
        if gene_list:
            src_col = (
                "source_genesymbol"
                if "source_genesymbol" in df.columns
                else "source"
            )
            tgt_col = (
                "target_genesymbol"
                if "target_genesymbol" in df.columns
                else "target"
            )

            df = df[
                df[src_col].isin(gene_list) |
                df[tgt_col].isin(gene_list)
            ]

        # ---------------------------------------------------------
        # Save output (optional)
        # ---------------------------------------------------------
        if save and output_dir is not None:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            if sample_id is None:
                sample_id = "omnipath"

            out_file = output_dir / f"{sample_id}_omnipath_interactions.tsv"
            df.to_csv(out_file, sep="\t", index=False)

            print(f"💾 OmniPath interactions saved to: {out_file}")

        return df

    except Exception as e:
        print(f"❌ OmniPath client error: {e}")
        return pd.DataFrame()


# --- Usage ---
if __name__ == "__main__":
    my_genes = ["TP53", "EGFR", "MYC", "BRCA1"]
    
    df_result = run_omnipath_interactions(my_genes)
    
    if not df_result.empty:
        print(f"\nSuccess! Found {len(df_result)} interactions.")
        print(df_result.head())
    else:
        print("No interactions found.")
        