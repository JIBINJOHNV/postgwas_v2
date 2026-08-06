from gprofiler import GProfiler
import pandas as pd
from scipy.stats import hypergeom
from pathlib import Path

# https://biit.cs.ut.ee/gprofiler/page/apis


def run_gprofiler_with_raw_stats(
    gene_list,
    output_dir=None,
    sample_id=None,
    save=True,
):
    """
    Run g:Profiler enrichment and compute raw (unadjusted) p-values.

    Args:
        gene_list (list): List of gene symbols.
        output_dir (str | Path | None): Directory to save g:Profiler results.
        sample_id (str | None): Sample/dataset identifier for filenames.
        save (bool): Whether to save output to disk.

    Returns:
        pd.DataFrame: g:Profiler enrichment results with raw p-values added.
    """

    gp = GProfiler(return_dataframe=True)

    results = gp.profile(
        organism="hsapiens",
        query=gene_list,
        sources=[
            "GO:BP", "GO:MF", "GO:CC",
            "KEGG", "REAC", "WP",
            "TF", "MIRNA", "HPA",
            "CORUM", "HP",
        ],
        user_threshold=1.0,                  # allow non-significant
        significance_threshold_method="fdr",
        all_results=True,
        no_evidences=False,
    )

    if results.empty:
        print("ℹ️ No g:Profiler enrichment results.")
        return results

    # ---------------------------------------------------------
    # Calculate raw (unadjusted) p-values using hypergeometric SF
    # ---------------------------------------------------------
    results["raw_p_value"] = results.apply(
        lambda row: hypergeom.sf(
            row["intersection_size"] - 1,
            row["effective_domain_size"],
            row["term_size"],
            row["query_size"],
        ),
        axis=1,
    )

    # ---------------------------------------------------------
    # Save output (optional)
    # ---------------------------------------------------------
    if save and output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if sample_id is None:
            sample_id = "gprofiler"

        out_file = output_dir / f"{sample_id}_gprofiler_enrichment.tsv"
        results.to_csv(out_file, sep="\t", index=False)

        print(f"💾 g:Profiler results saved to: {out_file}")

    return results




# --- Main Execution ---
if __name__ == "__main__":
    # 1. Input Genes
    my_genes = ["TP53", "BRCA1", "EGFR", "TNF", "APOE", "IL6"]
    print(f"Input Genes: {my_genes}")
    # 2. Run Analysis
    print("Running g:Profiler enrichment...")
    df = run_gprofiler_with_raw_stats(my_genes)
    if not df.empty:
        # 3. Clean and Sort Data
        # Select relevant columns
        cols = ['source', 'native', 'name', 'p_value', 'term_size', 'query_size', 'intersection_size', 'intersections']
    
        # 'intersections' is a list of gene names involved in the term
        # Let's clean it up if it exists in the dataframe
        if 'intersections' in df.columns:
            # Flatten list to string
            df['intersections'] = df['intersections'].apply(lambda x: ", ".join(x) if isinstance(x, list) else str(x))

        # Filter columns that exist
        final_cols = [c for c in cols if c in df.columns]
        df_clean = df[final_cols].sort_values(by=['source', 'p_value'])

        # 4. Display Results
        print("\nTop 5 Results:")
        print(df_clean.head().to_string(index=False))

        # OPTIONAL: Save to CSV
        # df_clean.to_csv("gprofiler_results.csv", index=False)
        # print("\nResults saved to gprofiler_results.csv")
    else:
        print("No results returned.")