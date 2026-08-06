import pandas as pd
import numpy as np
import glob,os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.cluster.hierarchy import fcluster
from sklearn.cluster import AgglomerativeClustering

def robust_read_csv(file, id_col, comment='#'):
    """
    Tries multiple parsing strategies and cleans whitespace to handle MAGMA/GWAS formats.
    """
    strategies = [
        {'sep': r'\s+', 'engine': 'python'}, # Whitespace (MAGMA, PLINK)
        {'sep': None, 'engine': 'python'},   # Auto-detect (TSV/CSV)
        {'sep': ',', 'engine': 'c'}          # Standard CSV
    ]
    for strategy in strategies:
        try:
            # Test with a few rows
            df = pd.read_csv(file, **strategy, comment=comment, nrows=10)
            df.columns = df.columns.str.strip() 
            if id_col in df.columns:
                full_df = pd.read_csv(file, **strategy, comment=comment)
                full_df.columns = full_df.columns.str.strip()
                # Cast ID to string to ensure reliable merging
                full_df[id_col] = full_df[id_col].astype(str).str.strip()
                return full_df
        except Exception:
            continue
    raise ValueError(f"Could not find ID column '{id_col}' in {file}")

def create_gene_matrices(pattern, target_cols, id_col='VARIABLE', p_col=None, 
                         beta_col=None, se_col=None, output_prefix="merged_results"):
    """
    Extracts data, merges into matrices, and saves each matrix to a CSV file.
    """
    # Create output directory if it doesn't exist
    file_list = glob.glob(pattern)
    if not file_list:
        print(f"Warning: No files found for pattern: {pattern}")
        return {}
    calc_z = (beta_col is not None and se_col is not None)
    actual_targets = target_cols + ['Z_SCORE'] if calc_z else target_cols
    merged_dfs = {col: None for col in actual_targets}
    for file in file_list:
        path_parts = file.split(os.sep)
        sample_name = path_parts[0] if path_parts[0] != '.' else path_parts[1]
        try:
            t_df = robust_read_csv(file, id_col, comment='#')
            if calc_z:
                if beta_col in t_df.columns and se_col in t_df.columns:
                    t_df['Z_SCORE'] = pd.to_numeric(t_df[beta_col], errors='coerce') / \
                                      pd.to_numeric(t_df[se_col], errors='coerce')
            for col in actual_targets:
                if col not in t_df.columns:
                    continue
                temp_df = t_df[[id_col, col]].copy()
                temp_df[col] = pd.to_numeric(temp_df[col], errors='coerce')
                temp_df = temp_df.dropna(subset=[col])
                if col == p_col:
                    temp_df = temp_df.sort_values(col, ascending=True)
                else:
                    temp_df['abs_sort'] = temp_df[col].abs()
                    temp_df = temp_df.sort_values('abs_sort', ascending=False)
                temp_df = temp_df.drop_duplicates(subset=[id_col]).drop(columns=['abs_sort'], errors='ignore')
                temp_df = temp_df.rename(columns={col: sample_name})
                if merged_dfs[col] is None:
                    merged_dfs[col] = temp_df
                else:
                    merged_dfs[col] = pd.merge(merged_dfs[col], temp_df, on=id_col, how='outer', suffixes=('', '_dup'))
                    merged_dfs[col] = merged_dfs[col].loc[:, ~merged_dfs[col].columns.str.endswith('_dup')]
        except Exception as e:
            print(f"Error processing {sample_name} ({file}): {e}")
    # Final Processing and CSV Writing
    for col in actual_targets:
        if merged_dfs[col] is not None:
            merged_dfs[col] = merged_dfs[col].set_index(id_col)
            fill_val = 1.0 if col == p_col else 0.0
            merged_dfs[col] = merged_dfs[col].fillna(fill_val)
            # Create a safe filename (e.g., "BETA_STD_merged.csv")
            file_name = f"{output_prefix}_{col.replace(' ', '_')}_merged.csv"
            file_path = os.path.join(file_name)
            merged_dfs[col].to_csv(file_path)
            print(f"Successfully saved: {file_path}")
    return merged_dfs

def get_elbow_k_refined(
    linkage_matrix,
    prefix="analysis",
    visualize=True,
    sensitivity=0.8,
    k_min=2,
    k_max=8,
    min_cluster_size=5,
    max_cluster_frac=0.70
):
    """
    Elbow + size-regularized clustering tuned for N≈60.
    Stable, balanced, biologically interpretable.
    """
    # ------------------------------------------------------------
    # 1. Distances and true K axis
    # ------------------------------------------------------------
    y_raw = linkage_matrix[:, 2][::-1]
    n = len(y_raw)
    # True number of clusters corresponding to each merge
    k_values = np.arange(n, 0, -1)   # n, n-1, ..., 1
    # ------------------------------------------------------------
    # 2. Gentle smoothing (important for N=60)
    # ------------------------------------------------------------
    window = max(3, min(7, n // 12))
    y = (
        pd.Series(y_raw)
        .rolling(window=window, center=True)
        .mean()
        .bfill()
        .ffill()
        .values
    )
    # ------------------------------------------------------------
    # 3. Raw elbow detection in K-space
    # ------------------------------------------------------------
    line_start = np.array([k_values[0], y[0]])
    line_end   = np.array([k_values[-1], y[-1]])
    vec = line_end - line_start
    norm = np.linalg.norm(vec)
    if norm == 0:
        elbow_k = k_min
    else:
        unit = vec / norm
        distances = np.zeros(n)
        for i in range(n):
            p = np.array([k_values[i], y[i]])
            v = p - line_start
            proj = np.dot(v, unit) * unit
            distances[i] = np.linalg.norm(v - proj)
        penalty = np.exp(-np.arange(n) / (n / sensitivity))
        adjusted = distances * penalty
        elbow_k = int(k_values[np.argmax(adjusted)])
    # ------------------------------------------------------------
    # 4. Balanced K selection (core logic)
    # ------------------------------------------------------------
    best_k = None
    best_score = -np.inf
    n_items = linkage_matrix.shape[0] + 1
    for k in range(k_min, k_max + 1):
        labels = fcluster(linkage_matrix, k, criterion="maxclust")
        sizes = np.bincount(labels)[1:]
        # Hard constraints
        if sizes.min() < min_cluster_size:
            continue
        if sizes.max() / n_items > max_cluster_frac:
            continue
        # Balanced score (robust)
        cv = sizes.std() / sizes.mean()
        balance = np.exp(-cv)
        # Prefer close to elbow
        elbow_pref = np.exp(-abs(k - elbow_k) / 2.5)
        score = balance * elbow_pref
        if score > best_score:
            best_score = score
            best_k = k
    # Fallback (very rare now)
    if best_k is None:
        best_k = min(max(elbow_k, 3), 6)
    # ------------------------------------------------------------
    # 5. Visualization
    # ------------------------------------------------------------
    if visualize:
        plt.figure(figsize=(6, 4))
        plt.plot(k_values, y_raw, color='lightgrey', alpha=0.4, label='Raw')
        plt.plot(k_values, y, color='black', lw=2, label='Smoothed')
        plt.axvline(x=elbow_k, color='blue', linestyle='--', label=f'Raw elbow K={elbow_k}')
        plt.axvline(x=best_k, color='red', linestyle='--', lw=2, label=f'Final K={best_k}')
        plt.xlabel("Number of clusters (K)")
        plt.ylabel("Linkage height")
        plt.title(f"Elbow tuned for N≈60: {prefix}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{prefix}_elbow_refined.png", dpi=300)
        plt.close()
    return int(best_k)


def run_cluster_analysis_constrained(df, prefix="analysis", q_threshold=0.80, min_samples=5):
    """
    Optimized to prevent K=1 by using multi-stage fallback constraints.
    """
    # 1. Preparation & Filtering
    matrix = df.copy().fillna(0)
    matrix['abs_max'] = matrix.abs().max(axis=1)
    thresh = matrix['abs_max'].quantile(q_threshold)
    filtered = matrix[matrix['abs_max'] > thresh].drop(columns=['abs_max'])
    if filtered.empty or len(filtered.columns) < 2:
        print(f"[{prefix}] Not enough data to cluster.")
        return None
    # 2. Hierarchical Linkage
    row_linkage = linkage(filtered, method='ward')
    col_linkage = linkage(filtered.T, method='ward')
    # 3. Determine Initial K
    initial_col_k = get_elbow_k_refined(col_linkage, prefix=f"{prefix}_cols", 
                                        sensitivity=1.2, # Higher sensitivity to find splits
                                        min_cluster_size=min_samples)
    final_col_k = 1
    col_assignments = np.ones(len(filtered.columns), dtype=int)
    # --- STRATEGY A: Try strict constraint (min_samples=5) ---
    for k_try in range(max(initial_col_k, 2), 1, -1):
        labels = fcluster(col_linkage, t=k_try, criterion='maxclust')
        if pd.Series(labels).value_counts().min() >= min_samples:
            col_assignments, final_col_k = labels, k_try
            break
    # --- STRATEGY B: Fallback (min_samples=3) ---
    if final_col_k == 1:
        for k_try in range(max(initial_col_k, 2), 1, -1):
            labels = fcluster(col_linkage, t=k_try, criterion='maxclust')
            if pd.Series(labels).value_counts().min() >= 3:
                col_assignments, final_col_k = labels, k_try
                print(f"[{prefix}] Relaxed constraint used: min_samples=3")
                break
    # --- STRATEGY C: Force Split (K=2) ---
    if final_col_k == 1:
        col_assignments = fcluster(col_linkage, t=2, criterion='maxclust')
        final_col_k = 3
        print(f"[{prefix}] Forced split used: K=2 (Constraint could not be met)")
    # 4. Feature (Row) Clustering
    row_k = get_elbow_k_refined(row_linkage, prefix=f"{prefix}_rows", sensitivity=0.8)
    row_assignments = fcluster(row_linkage, t=row_k, criterion='maxclust')
    # 5. Save Results
    df_cols = pd.DataFrame({'Sample': filtered.columns, 'Cluster': col_assignments})
    df_cols.to_csv(f"{prefix}_sample_clusters.csv", index=False)
    df_rows = pd.DataFrame({'Feature_ID': filtered.index, 'Cluster': row_assignments})
    df_rows.to_csv(f"{prefix}_feature_clusters.csv", index=False)
    # 6. Reporting
    print(f"\n" + "="*40 + f"\nSUCCESS: {prefix.upper()}\n" + "-"*40)
    print(df_cols['Cluster'].value_counts().sort_index().to_string())
    print("="*40)
    # 7. Visualization
    sns.clustermap(filtered, cmap='RdBu_r', center=0, 
                   row_linkage=row_linkage, col_linkage=col_linkage,
                   figsize=(12, 10), xticklabels=True)
    plt.savefig(f"{prefix}_clustermap.png", dpi=300, bbox_inches='tight')
    plt.close()
    return df_cols

magma_genes=create_gene_matrices(pattern='*/05_magma/*.genes.out', 
                     target_cols=['ZSTAT', 'P'], id_col='GENE',p_col='P',output_prefix="magma_genes")

magma_tissue=create_gene_matrices(pattern='*/06_magma_covar/*.gsa.out', 
                     target_cols=['BETA_STD', 'P'], id_col='FULL_NAME',p_col='P',output_prefix="magma_tissue")

ld_pruned=create_gene_matrices(pattern='*/04_ld_clump/*LDpruned_EUR.tsv', 
                     target_cols=['P_value'], id_col='LDblock',
                     beta_col='BETA', se_col='SE',p_col='P_value', output_prefix="ld_pruned")

pops_pred=create_gene_matrices(pattern='*/07_pops/*pops.preds', 
                     target_cols=['PoPS_Score'], id_col='ENSGID', output_prefix="pops_pred")

flames=create_gene_matrices(pattern='*/09_flames/FLAMES_scores.raw', 
                     target_cols=['FLAMES_scaled'], id_col='symbol', output_prefix="flames")



run_cluster_analysis_constrained(df=magma_genes['ZSTAT'], prefix="magma_gene_zstat", q_threshold=0.75)
run_cluster_analysis_constrained(df=magma_genes['P'], prefix="magma_gene_P", q_threshold=0.75)

run_cluster_analysis_constrained(df=magma_tissue['BETA_STD'], prefix="magma_tissue_beta", q_threshold=0.75)
run_cluster_analysis_constrained(df=magma_tissue['P'], prefix="magma_tissue_P", q_threshold=0.75)

run_cluster_analysis_constrained(df=ld_pruned['Z_SCORE'], prefix="ld_zscore", q_threshold=0.75)
run_cluster_analysis_constrained(df=pops_pred['PoPS_Score'], prefix="PoPS_Score", q_threshold=0.75)
run_cluster_analysis_constrained(df=flames['FLAMES_scaled'], prefix="FLAMES_scaled", q_threshold=0.0001)




files=glob.glob('*_sample_clusters.csv')

def create_consensus_matrix(file_list, output_name="consensus_clusters.csv"):
    consensus_df = None
    
    for file in file_list:
        if not os.path.exists(file):
            print(f"Skipping {file}: File not found.")
            continue
        # 1. Load the CSV
        df = pd.read_csv(file)
        # 2. Generate the new column name
        # Remove '_sample_clusters.csv' from the filename
        new_col_name = file.replace('_sample_clusters.csv', '')
        # 3. Rename 'Cluster' column
        df = df.rename(columns={'Cluster': new_col_name})
        # 4. Merge into the master dataframe
        if consensus_df is None:
            consensus_df = df
        else:
            # Outer merge ensures that if a sample is missing in one analysis, it still appears
            consensus_df = pd.merge(consensus_df, df, on='Sample', how='outer')
    if consensus_df is not None:
        # Save the result
        consensus_df.to_csv(output_name, index=False)
        print(f"Successfully created consensus matrix: {output_name}")
        return consensus_df
    else:
        print("No files were merged.")
        return None
    
    
# Execute the merger
master_clusters = create_consensus_matrix(files)

master_clusters.to_csv("Phenotype_combined_clusters.csv", index=False)
# Display a preview
print(master_clusters.head())