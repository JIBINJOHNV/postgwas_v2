import os
import sys
import logging
import subprocess
from pathlib import Path
from datetime import datetime
from typing import List

import polars as pl

# =============================================================================
# USER CONFIG
# =============================================================================
input_files='/Users/JJOHN41/Downloads/ukb_ppp_testing/vcf_files/ukb_ppp_testing_input_susmstat.csv'
input_df=pd.read_csv(input_files)

TOL = {
    "BETA": 1e-6,
    "EAF": 1e-4,
    "SE": 1e-6,
    "P": 1e-4,
    "Zscore": 1e-4,
    'INFO': 1e-4,
    "TotalSampleSize": 0.0,
    "NCONTROL": 0.0,
    "NEffective": 0.0,
}

# =============================================================================
# LOGGING & HELPERS
# =============================================================================

def setup_logging(output_dir: str) -> str:
    os.makedirs(output_dir, exist_ok=True)
    log_file = os.path.join(output_dir, f"qc_validation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.handlers = []
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    fh = logging.FileHandler(log_file)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(formatter)
    logger.addHandler(sh)
    return log_file



def ensure_minus_log10_p(df: pl.DataFrame, col_name: str = "P_FINAL") -> pl.DataFrame:
    """
    Checks if a column is raw P-values (max <= 1.0).
    If it is raw, converts to -log10(P). 
    If it is already transformed (contains values > 1.1), leaves it as is.
    """
    # 1. Get the max value of the column to detect format
    # We use a small buffer (1.1) to account for floating point precision
    max_p = df.select(pl.col(col_name).max()).to_series()[0]
    if max_p is None:
        return df # Handle empty or all-null columns
    if max_p <= 1.1:
        print(f"Detected RAW P-values (Max: {max_p:.4f}). Converting {col_name} to -log10...")
        # Convert raw P to -log10, handling potential zeros to avoid infinity
        df = df.with_columns(
            pl.col(col_name).log10().mul(-1).alias(col_name)
        )
    else:
        print(f"Detected ALREADY TRANSFORMED P-values (Max: {max_p:.4f}). No conversion needed.")
    return df


def ensure_beta_format(df: pl.DataFrame, col_name: str = "BETA_RAW") -> pl.DataFrame:
    """
    Detects if a column is Odds Ratio (OR) or Beta.
    If OR (all positive, median ~1), converts to Beta via natural log.
    If already Beta (contains negatives, median ~0), leaves it as is.
    """
    # 1. Calculate min and median to determine the scale
    # We use .select() for speed on 16M rows
    stats = df.select([
        pl.col(col_name).min().alias("min_val"),
        pl.col(col_name).median().alias("median_val")
    ])
    min_val = stats.get_column("min_val")[0]
    median_val = stats.get_column("median_val")[0]
    if min_val is None:
        return df
    # 2. Logic: Odds Ratios are strictly positive and usually center around 1.0
    # Betas (Log-Odds) can be negative and center around 0.0
    if min_val >= 0 and median_val > 0.5:
        print(f"Detected ODDS RATIOS (Min: {min_val:.4f}, Median: {median_val:.4f}).")
        print(f"Converting {col_name} to Log-Odds (Beta)...")
        # log() in Polars is the natural logarithm (ln), which is the standard for Log-Odds
        df = df.with_columns(
            pl.col(col_name).log().alias(col_name)
        )
    else:
        print(f"Detected BETA/LOG-ODDS (Min: {min_val:.4f}, Median: {median_val:.4f}). No conversion needed.")
    return df


def require_columns(df: pl.DataFrame, cols: List[str], df_name: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{df_name} missing required columns: {missing}")

def standardize_chrom_expr(col_name: str) -> pl.Expr:
    return pl.col(col_name).cast(pl.Utf8).str.strip().str.replace(r"^chr", "", literal=False)

def allele_expr(col_name: str) -> pl.Expr:
    return pl.col(col_name).cast(pl.Utf8).str.strip().str.to_uppercase()

def palindromic_expr(ref_col: str, alt_col: str) -> pl.Expr:
    return (
        ((pl.col(ref_col) == "A") & (pl.col(alt_col) == "T")) |
        ((pl.col(ref_col) == "T") & (pl.col(alt_col) == "A")) |
        ((pl.col(ref_col) == "C") & (pl.col(alt_col) == "G")) |
        ((pl.col(ref_col) == "G") & (pl.col(alt_col) == "C"))
    )

def variant_type_expr(ref_col: str, alt_col: str) -> pl.Expr:
    return pl.when((pl.col(ref_col).str.len_chars() > 1) | (pl.col(alt_col).str.len_chars() > 1)).then(pl.lit("INDEL")).otherwise(pl.lit("SNP"))

def safe_corr(df: pl.DataFrame, col1: str, col2: str, method: str = "pearson") -> float:
    tmp = df.select([col1, col2]).drop_nulls()
    if tmp.height < 2 or tmp.select(pl.col(col1).n_unique()).item() < 2 or tmp.select(pl.col(col2).n_unique()).item() < 2:
        return float("nan")
    try:
        return tmp.select(pl.corr(col1, col2, method=method)).item()
    except Exception:
        return float("nan")

def safe_sign_concordance(df: pl.DataFrame, col1: str, col2: str) -> float:
    tmp = df.select([col1, col2]).drop_nulls()
    if tmp.height == 0: return float("nan")
    return tmp.select(((pl.col(col1).sign() == pl.col(col2).sign()).cast(pl.Float64).mean() * 100)).item()

# =============================================================================
# STEPS
# =============================================================================

def load_raw_sumstats(raw_tsv_path: str, tsv_map: dict, optional_map: dict) -> pl.DataFrame:
    logging.info("[STEP 1] Reading raw TSV...")
    df = pl.read_csv(raw_tsv_path, separator="\t", dtypes={"CHROM": pl.Utf8}, infer_schema_length=10000, null_values=["", "NA", "NaN", "nan", "."])
    
    rename_map = {v: k for k, v in tsv_map.items() if v in df.columns}
    for k, v in optional_map.items():
        if v and v in df.columns: rename_map[v] = k
    df = df.rename(rename_map)
    df = df.with_columns([
        (pl.col("CHROM").cast(pl.Utf8) + "_" + pl.col("POS").cast(pl.Utf8) + "_" + pl.col("REF") + "_" + pl.col("ALT")).alias("ID_RAW"),
        (pl.col("CHROM").cast(pl.Utf8) + "_" + pl.col("POS").cast(pl.Utf8) + "_" + pl.col("ALT") + "_" + pl.col("REF")).alias("ID_RAW_REVERSE")
    ])
    exprs = [
        standardize_chrom_expr("CHROM").alias("CHROM"),
        pl.col("POS").cast(pl.Int64, strict=False),
        allele_expr("REF").alias("REF"),
        allele_expr("ALT").alias("ALT"),
        pl.col("BETA").cast(pl.Float64, strict=False).alias("BETA_RAW"),
        pl.col("SE").cast(pl.Float64, strict=False).alias("SE_RAW"),
        pl.col("P").cast(pl.Float64, strict=False).alias("P_RAW"),
        pl.col("EAF").cast(pl.Float64, strict=False).alias("EAF_RAW"),
        pl.col("NCONTROL").cast(pl.Float64, strict=False).alias("NCONTROL_RAW"),
        pl.col("INFO").cast(pl.Float64, strict=False).alias("INFO_RAW"),
        variant_type_expr("REF", "ALT").alias("VARIANT_TYPE_RAW"),
        palindromic_expr("REF", "ALT").alias("PALINDROMIC_RAW")
    ]
    df = df.with_columns(exprs).drop_nulls(subset=["CHROM", "POS", "REF", "ALT"])
    df = df.drop(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'EAF', 'INFO', 'NCONTROL', 'BETA', 'SE', 'P'])
    return df
    

def extract_vcf_to_tsv(vcf_path: str, out_tsv_path: str) -> None:
    logging.info("[STEP 2] Extracting VCF fields...")
    cmd = f"(echo -e 'CHROM\\tPOS\\tID\\tREF\\tALT\\tES\\tSE\\tEZ\\tLP\\tAF\\tSI\\tSS\\tNCO\\tNC\\tNEF'; bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t[%ES]\\t[%SE]\\t[%EZ]\\t[%LP]\\t[%AF]\\t[%SI]\\t[%SS]\\t[%NCO]\\t[%NC]\\t[%NEF]\\n' '{vcf_path}') > '{out_tsv_path}'"
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")

def load_vcf_table(vcf_tsv_path: str) -> pl.DataFrame:
    df = pl.read_csv(vcf_tsv_path, separator="\t", dtypes={"CHROM": pl.Utf8}, infer_schema_length=10000)
    df = df.with_columns([
        (pl.col("CHROM").cast(pl.Utf8) + "_" + pl.col("POS").cast(pl.Utf8) + "_" + pl.col("REF") + "_" + pl.col("ALT")).alias("ID_VCF"),
        (pl.col("CHROM").cast(pl.Utf8) + "_" + pl.col("POS").cast(pl.Utf8) + "_" + pl.col("ALT") + "_" + pl.col("REF")).alias("ID_VCF_REVERSE")
    ]).rename({
        "ES": "BETA_VCF", "SE": "SE_VCF", "EZ": "Zscore_VCF", "LP": "P_VCF", "AF": "EAF_VCF", 
        "SS": "TotalSampleSize_VCF", "NCO": "NCONTROL_VCF",  "NEF": "NEffective_VCF",'SI':'INFO_VCF'
    })
    try: 
       df = df.rename({"NC":"NCASE_VCF"})
    except: 
        pass
    df=df.with_columns([
        variant_type_expr("REF", "ALT").alias("VARIANT_TYPE_VCF"),
        palindromic_expr("REF", "ALT").alias("PALINDROMIC_VCF")
    ])
    df=df.drop(['CHROM', 'POS', 'ID', 'REF', 'ALT'])
    return df

def merge_direct(raw_df: pl.DataFrame, vcf_df: pl.DataFrame) -> pl.DataFrame:
    merged = raw_df.join(vcf_df, left_on="ID_RAW", right_on="ID_VCF", how="inner")
    merged = merged.with_columns([
        pl.lit("DIRECT").alias("MATCH_TYPE"),
        pl.col("ID_RAW").alias("ID_VCF"),
        pl.col("BETA_VCF").alias("BETA_FINAL"),
        pl.col("SE_VCF").alias("SE_FINAL"),
        pl.col("Zscore_VCF").alias("Zscore_FINAL"),
        pl.col("P_VCF").alias("P_FINAL"),
        pl.col("EAF_VCF").alias("EAF_FINAL"),
        pl.col("INFO_VCF").alias("INFO_FINAL"),
        pl.col("TotalSampleSize_VCF").alias("TotalSampleSize_FINAL"),
        pl.col("NCONTROL_VCF").alias("NCONTROL_FINAL"),
        pl.col("NEffective_VCF").alias("NEffective_FINAL"),
        pl.col("VARIANT_TYPE_RAW").alias("VARIANT_TYPE"),
        pl.col("PALINDROMIC_RAW").alias("PALINDROMIC")
    ])
    try: 
        merged = merged.rename({"NCASE_VCF":"NCASE_FINAL"})
    except: 
        pass
    return merged



def merge_flip(raw_df: pl.DataFrame, vcf_df: pl.DataFrame) -> pl.DataFrame:
    merged = raw_df.join(vcf_df, left_on="ID_RAW", right_on="ID_VCF_REVERSE", how="inner")
    return merged.with_columns([
        pl.lit("FLIP").alias("MATCH_TYPE"),
        pl.col("ID_RAW").alias("ID_VCF_REVERSE"),
        (-pl.col("BETA_VCF")).alias("BETA_FINAL"),
        (1.0 - pl.col("EAF_VCF")).alias("EAF_FINAL"),
        (-pl.col("Zscore_VCF")).alias("Zscore_FINAL"),
        pl.col("SE_VCF").alias("SE_FINAL"),
        pl.col("P_VCF").alias("P_FINAL"),
        pl.col("TotalSampleSize_VCF").alias("TotalSampleSize_FINAL"),
        pl.col("NCONTROL_VCF").alias("NCONTROL_FINAL"),
        pl.col("NEffective_VCF").alias("NEffective_FINAL"),
        pl.col("INFO_VCF").alias("INFO_FINAL"),
        pl.col("VARIANT_TYPE_RAW").alias("VARIANT_TYPE"),
        pl.col("PALINDROMIC_RAW").alias("PALINDROMIC")
    ])




def add_qc_columns(df: pl.DataFrame) -> pl.DataFrame:
    if df.height == 0: return df
    # 1. Calculate Differences using your exact column names
    df = df.with_columns([
        (pl.col("BETA_RAW") - pl.col("BETA_FINAL")).abs().alias("BETA_DIFF"),
        (pl.col("EAF_RAW") - pl.col("EAF_FINAL")).abs().alias("EAF_DIFF"),
        (pl.col("SE_RAW") - pl.col("SE_FINAL")).abs().alias("SE_DIFF"),
        (pl.col("P_RAW") - pl.col("P_FINAL")).abs().alias("P_DIFF"),
        ((pl.col("BETA_RAW") / pl.col("SE_RAW")) - pl.col("Zscore_FINAL")).abs().alias("Zscore_DIFF"),
        (pl.col("NCONTROL_RAW") - pl.col("NCONTROL_FINAL")).abs().alias("NCONTROL_DIFF"),
        (pl.col("INFO_RAW") - pl.col("INFO_FINAL")).abs().alias("INFO_DIFF"),
    ])
    # 2. Define the metrics based on your DIFF columns
    metrics = ["BETA", "EAF", "SE", "P", "Zscore", "NCONTROL"]
    # 3. Flag specific reasons based on Tolerance (TOL)
    for col in metrics:
        df = df.with_columns(
            pl.when(pl.col(f"{col}_DIFF") > TOL[col])
            .then(pl.lit(f"{col}_DIFF_GT_CUTOFF"))
            .otherwise(pl.lit(""))
            .alias(f"REASON_{col}")
        )
    # 4. Concatenate all reasons into 'REASON_ALL' with regex cleaning
    return (
        df.with_columns(
            pl.concat_str([pl.col(f"REASON_{c}") for c in metrics], separator="|")
            .str.replace_all(r"\|+", "|")
            .str.replace(r"^\|", "")
            .str.replace(r"\|$", "")
            .alias("REASON_ALL")
        )
        .with_columns(
            pl.when(pl.col("REASON_ALL") == "")
            .then(pl.lit("OK"))
            .otherwise(pl.col("REASON_ALL"))
            .alias("REASON_ALL")
        )
    )


def summarize_gwas_stats_polars(df: pl.DataFrame) -> pl.DataFrame:
    """
    Calculates GWAS statistics for Grouped and Global data.
    Returns separate 'Group' and 'Metric' columns.
    """
    # 1. Setup Columns
    raw_cols = ['BETA_RAW', 'SE_RAW', 'P_RAW', 'EAF_RAW', 'NCONTROL_RAW', 'INFO_RAW']
    final_cols = [
        'BETA_FINAL', 'SE_FINAL', 'Zscore_FINAL', 'P_FINAL',  
        'EAF_FINAL', 'INFO_FINAL', 'TotalSampleSize_FINAL', 'NCONTROL_FINAL', 'NEffective_FINAL'
    ]
    target_cols = [c for c in (raw_cols + final_cols) if c in df.columns]
    threshold_checks = [
        ("P_RAW", 7.3), ("P_FINAL", 7.3),
        ("INFO_RAW", 1.0), ("INFO_FINAL", 1.0)
    ]
    # 2. Aggregation Helper
    def get_aggs():
        exprs = []
        for col in target_cols:
            exprs.extend([
                pl.col(col).mean().alias(f"{col}__mean"),
                pl.col(col).median().alias(f"{col}__median"),
                pl.col(col).max().alias(f"{col}__max"),
                pl.col(col).min().alias(f"{col}__min"),
                (pl.col(col).std() / pl.col(col).count().sqrt()).alias(f"{col}__sem")
            ])
        for col_name, thresh in threshold_checks:
            if col_name in df.columns:
                exprs.extend([
                    (pl.col(col_name) > thresh).sum().alias(f"{col_name}_gt_{thresh}__threshold_count"),
                    ((pl.col(col_name) > thresh).sum() / pl.count() * 100).alias(f"{col_name}_gt_{thresh}__threshold_percentage")
                ])
        return exprs
    # 3. Calculate Stats
    grouped_res = df.group_by("VARIANT_TYPE").agg(get_aggs()).rename({"VARIANT_TYPE": "Group"})
    global_res = df.select(get_aggs()).with_columns(pl.lit("GLOBAL").alias("Group"))
    # 4. Combine and Reshape
    combined_res = pl.concat([grouped_res, global_res], how="diagonal")
    melted = combined_res.melt(id_vars="Group")
    # Split the "variable" column into the actual Metric name and the Stat type
    melted = melted.with_columns([
        pl.col("variable").str.splitn("__", 2).struct.field("field_0").alias("Metric"),
        pl.col("variable").str.splitn("__", 2).struct.field("field_1").alias("Stat")
    ])
    # 5. Pivot Stats into Columns
    # positional: values, index, columns
    final_df = melted.pivot("value", ["Group", "Metric"], "Stat")
    # 6. Final Schema Alignment
    required_cols = ["Group", "Metric", "mean", "median", "max", "min", "sem", "threshold_count", "threshold_percentage"]
    # Ensure all columns exist (important for the threshold rows)
    for c in required_cols:
        if c not in final_df.columns:
            final_df = final_df.with_columns(pl.lit(None).cast(pl.Float64).alias(c))
    return final_df.select(required_cols).sort(["Group", "Metric"])



def compute_set_qc(
    raw_df: pl.DataFrame,
    vcf_df: pl.DataFrame,
    final_df: pl.DataFrame
) -> (pl.DataFrame, pl.DataFrame, pl.DataFrame):
    """
    Computes set-theory counts (Raw, VCF, Common, Only-Raw, Only-VCF) 
    grouped by Variant Type and Global totals.
    """
    # 1. Standardize the Grouping Column for each set
    # We create a 'GROUP' column so we can join different dataframes later
    raw_processed = raw_df.with_columns(pl.col("VARIANT_TYPE_RAW").alias("GROUP"))
    vcf_processed = vcf_df.with_columns(pl.col("VARIANT_TYPE_VCF").alias("GROUP"))
    final_processed = final_df.with_columns(pl.col("VARIANT_TYPE").alias("GROUP"))
    # 2. Identify Unique Keys for Anti-Joins
    matched_raw_keys = final_processed.select(["ID_RAW", "ID_RAW_REVERSE"]).unique()
    matched_vcf_keys = final_processed.select(["ID_VCF", "ID_VCF_REVERSE"]).unique()
    # 3. Calculate 'Only' Sets (Anti-joins)
    # Variants in Raw but not in the final matched set
    raw_only = raw_processed.join(matched_raw_keys, on=["ID_RAW", "ID_RAW_REVERSE"], how="anti")
    # Variants in VCF but not in the final matched set
    vcf_only = vcf_processed.join(matched_vcf_keys, on=["ID_VCF", "ID_VCF_REVERSE"], how="anti")
    # 4. Helper function to aggregate counts by GROUP
    def get_grp_counts(df: pl.DataFrame, col_name: pl.Expr):
        return (
            df.group_by("GROUP")
            .agg(pl.len().alias(col_name))
            .with_columns(pl.col("GROUP").cast(pl.Utf8))
        )
    # 5. Compute Count Tables
    raw_counts = get_grp_counts(raw_processed, "RAW_TOTAL")
    vcf_counts = get_grp_counts(vcf_processed, "VCF_TOTAL")
    common_counts = get_grp_counts(final_processed, "COMMON_FINAL")
    raw_only_counts = get_grp_counts(raw_only, "RAW_ONLY")
    vcf_only_counts = get_grp_counts(vcf_only, "VCF_ONLY")
    # 6. Join all summaries into one table
    # We use a full join to ensure we don't lose any variant types (e.g. if Indels only exist in one set)
    summary_grouped = (
        raw_counts
        .join(vcf_counts, on="GROUP", how="full", coalesce=True)
        .join(common_counts, on="GROUP", how="full", coalesce=True)
        .join(raw_only_counts, on="GROUP", how="full", coalesce=True)
        .join(vcf_only_counts, on="GROUP", how="full", coalesce=True)
        .fill_null(0)
    )
    # 7. Calculate Global Totals
    global_row = summary_grouped.select([
        pl.lit("GLOBAL").alias("GROUP"),
        pl.col("RAW_TOTAL").sum(),
        pl.col("VCF_TOTAL").sum(),
        pl.col("COMMON_FINAL").sum(),
        pl.col("RAW_ONLY").sum(),
        pl.col("VCF_ONLY").sum(),
    ])
    # 8. Assemble Final DataFrame
    summary_final = (
        pl.concat([global_row, summary_grouped])
        .with_columns([
            (pl.col("COMMON_FINAL") / pl.col("RAW_TOTAL") * 100).alias("RETAINED_PCT"),
            (pl.col("RAW_ONLY") / pl.col("RAW_TOTAL") * 100).alias("LOST_PCT")
        ])
        .sort("GROUP")
    )
    return raw_only, vcf_only, summary_final




def generate_qc_stats_polars(df: pl.DataFrame, label: str, total_raw_variants: int,TOL=TOL):
    """
    Calculates QC statistics with fixed schema compatibility for diagonal concatenation.
    Returns:
    1. A Summary DataFrame (Global + Groupby VARIANT_TYPE)
    2. A 'Mismatch' DataFrame of failed variants.
    """
    # 1. Identify Mismatches
    mismatch_conditions = [
        (pl.col(f"{m}_DIFF").abs() > TOL[m]) 
        for m in TOL if f"{m}_DIFF" in df.columns and "INFO" not in m
    ] 
    # Use fold for efficiency
    is_mismatch = pl.fold(acc=pl.lit(False), function=lambda acc, x: acc | x, exprs=mismatch_conditions)
    mismatch_df = df.filter(is_mismatch)
    # 2. Aggregation Helper
    def get_qc_aggs():
        af_raw_col = "EAF_RAW" if "EAF_RAW" in df.columns else "EAF"
        exprs = [
            pl.lit(label).alias("REPORT_LABEL"),
            pl.len().alias("MATCHED_N"),
            is_mismatch.sum().alias("TOTAL_MISMATCHES"),
            (is_mismatch.sum() / pl.len() * 100).alias("MISMATCH_PCT"),
            # Correlations using pl.corr(a, b) for max compatibility
            pl.corr(pl.col("BETA_RAW"), pl.col("BETA_FINAL")).alias("BETA_CORR"),
            pl.corr(pl.col(af_raw_col), pl.col("EAF_FINAL")).alias("EAF_CORR"),
            pl.corr(pl.col("SE_RAW"), pl.col("SE_FINAL")).alias("SE_CORR"),
            (((pl.col("BETA_RAW") * pl.col("BETA_FINAL")) > 0).sum() / pl.len() * 100).alias("SIGN_CONCORDANCE")
        ]
        for m, t in TOL.items():
            if f"{m}_DIFF" in df.columns:
                exprs.append((pl.col(f"{m}_DIFF").abs() > t).sum().alias(f"{m}_MISMATCH_N"))
        return exprs
    # 3. Calculate Grouped and Global separately
    # Explicitly casting Group to String to avoid SchemaError
    grouped_stats = (
        df.group_by("VARIANT_TYPE")
        .agg(get_qc_aggs())
        .with_columns(pl.col("VARIANT_TYPE").cast(pl.Utf8).alias("GROUP"))
        .drop("VARIANT_TYPE")
    )
    global_stats = (
        df.select(get_qc_aggs())
        .with_columns(pl.lit("GLOBAL").cast(pl.Utf8).alias("GROUP"))
    )
    # 4. Final Summary Table
    summary_df = pl.concat([global_stats, grouped_stats], how="diagonal")
    # Reorder to put GROUP at the front
    cols = ["REPORT_LABEL", "GROUP"] + [c for c in summary_df.columns if c not in ["REPORT_LABEL", "GROUP"]]
    summary_df = summary_df.select(cols)
    return summary_df, mismatch_df



def plot_info_histogram_three_lane(df: pl.DataFrame, col_name: str = "INFO_FINAL", out_dir: str = "./", sample_name: str = "sample"):
    # 1. Efficiently filter data using Polars
    # Top Data: Full distribution (dropping nulls)
    data_all = df.select(col_name).drop_nulls()[col_name].to_numpy()
    # Middle Data: Zoom into values > 1.0
    data_gt1 = df.filter(pl.col(col_name) > 1.0).select(col_name).to_series().to_numpy()
    # Bottom Data: Specific outliers > 1.1
    data_gt11 = df.filter(pl.col(col_name) > 1.1).select(col_name).to_series().to_numpy()
    # 2. Initialize Three-Lane Plot
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 15), sharex=False)
    # --- TOP LANE: Full Distribution ---
    ax1.hist(data_all, bins=100, color='#3498db', alpha=0.8)
    ax1.set_title(f'Overall Distribution: {col_name}', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Variant Count', fontsize=11)
    ax1.grid(axis='y', linestyle='--', alpha=0.3)
    ax1.axvline(0.3, color='red', linestyle='--', label='QC Cutoff (0.3)')
    ax1.legend(loc='upper left')
    # --- MIDDLE LANE: Values > 1.0 ---
    if len(data_gt1) > 0:
        ax2.hist(data_gt1, bins=50, color='#e67e22', alpha=0.8)
        ax2.set_title(f'Zoom: {col_name} > 1.0', fontsize=13)
    else:
        ax2.text(0.5, 0.5, "No variants found with INFO > 1.0", ha='center', va='center')
    ax2.set_ylabel('Variant Count', fontsize=11)
    ax2.grid(axis='y', linestyle='--', alpha=0.3)
    # --- BOTTOM LANE: Values > 1.1 ---
    if len(data_gt11) > 0:
        ax3.hist(data_gt11, bins=50, color='#e74c3c', alpha=0.8)
        ax3.set_title(f'Critical Outliers: {col_name} > 1.1', fontsize=13, color='darkred')
    else:
        ax3.text(0.5, 0.5, "No variants found with INFO > 1.1", ha='center', va='center')
    ax3.set_xlabel('INFO Score', fontsize=12)
    ax3.set_ylabel('Variant Count', fontsize=11)
    ax3.grid(axis='y', linestyle='--', alpha=0.3)
    # 3. Final Layout and Save
    plt.tight_layout()
    output_path = os.path.join(out_dir, f"{sample_name}_info_score_3lane.png")
    plt.savefig(output_path, dpi=300)
    plt.close()

# =============================================================================
# MAIN
# =============================================================================

def harmonisation_diagnosis_pipeline(input_df,index):
    input_dict=input_df.iloc[index,].to_dict()
    sample_name = input_dict['gwas_outputname']
    raw_tsv_file = input_dict['sumstat_file']
    #vcf_file = f"{input_dict['output_folder']}/00_harmonised_sumstat/{sample_name}_gwas2vcf_GRCh38_merged.vcf.gz"
    vcf_file = f"{input_dict['output_folder']}/00_harmonised_sumstat/{sample_name}_GRCh38_merged.vcf.gz"
    out_dir = f"{input_dict['output_folder']}/00_harmonised_sumstat/harmonisation_diagnosis/"
    os.makedirs(out_dir, exist_ok=True)
    tsv_map = {
        "CHROM":input_dict['chr_col'],
        "POS": input_dict['pos_col'],
        "REF": input_dict['oa_col'],
        "ALT": input_dict['ea_col'],
        "BETA": input_dict['beta_or_col'],
        "SE": input_dict['se_col'],
        "EAF": input_dict['eaf_col'],
        "P": input_dict['pval_col'],
        "NCONTROL": input_dict['ncontrol_col'],
    }
    optional_tsv_cols = {
        "ID": input_dict['snp_id_col'],
        "INFO": input_dict['imp_info_col'],
        "NCASE":input_dict['ncase_col'],
    }
    log_file = setup_logging(out_dir)
    logging.info(f"Starting QC for {sample_name}")
    raw_df = load_raw_sumstats(raw_tsv_file, tsv_map, optional_tsv_cols)
    raw_df =ensure_minus_log10_p(raw_df, col_name = "P_RAW")
    raw_df = ensure_beta_format(raw_df, col_name = "BETA_RAW")
    plot_info_histogram_three_lane(raw_df, col_name='INFO_RAW',out_dir=out_dir,sample_name=sample_name)
    total_raw_variants = raw_df.height
    extracted_vcf_tsv = os.path.join(out_dir, f"{sample_name}_vcf_extracted.tsv")
    extract_vcf_to_tsv(vcf_file, extracted_vcf_tsv)
    vcf_df = load_vcf_table(extracted_vcf_tsv).unique(subset=["ID_VCF"])
    # Merge Logic
    direct_df = add_qc_columns(merge_direct(raw_df, vcf_df))
    # Remaining for Flip
    flip_df = add_qc_columns(merge_flip(raw_df, vcf_df))
    final_df = pl.concat([direct_df, flip_df], how="diagonal_relaxed")
    os.system(f"rm -f {extracted_vcf_tsv}")
    # Stats & Problems 
    summary_df, mismatch_df = generate_qc_stats_polars(final_df, "FINAL_COMBINED", total_raw_variants)
    raw_only_df, vcf_only_df, set_stats = compute_set_qc(raw_df, vcf_df, final_df)
    final_stats_df = summarize_gwas_stats_polars(final_df).to_pandas()
    variant_count_df=pd.concat([set_stats.to_pandas().T,summary_df.to_pandas().T, ]).reset_index().drop_duplicates()
    variant_count_df.columns=variant_count_df[variant_count_df['index']=='GROUP'].iloc[0,]
    variant_count_df=variant_count_df[~variant_count_df['GROUP'].isin(['REPORT_LABEL','GROUP','MATCHED_N'])]
    variant_count_df['GROUP']=variant_count_df['GROUP'].replace({'COMMON_FINAL': 'RAW_VCF_COMMON','SIGN_CONCORDANCE':'BETA_SIGN_CONCORDANCE'})
    rename_map = {
            'RAW_TOTAL': 'INPUT_SUMSTATS_N',
            'VCF_TOTAL': 'HARMONISED_SUMSTATS_N',
            'RAW_VCF_COMMON': 'INPUT_HARMONISED_COMMON_N',
            'RAW_ONLY': 'Lost_During_Harmonisation_N',
            'VCF_ONLY': 'Gained_During_Harmonisation_N',
            'RETAINED_PCT': 'INPUT_HARMONISED_COMMON_%',
            'LOST_PCT': 'Lost_During_Harmonisation_%',
            'TOTAL_MISMATCHES': 'Total_Mismatches_N',
            'MISMATCH_PCT': 'Total_Mismatches_%'
        }
    variant_count_df['GROUP'] = variant_count_df['GROUP'].replace(rename_map)
    # Output Files 
    raw_only_df.write_csv(os.path.join(out_dir, f"{sample_name}_raw_only_variants.csv"))
    vcf_only_df.write_csv(os.path.join(out_dir, f"{sample_name}_vcf_only_variants.csv"))
    variant_count_df.to_csv(os.path.join(out_dir, f"{sample_name}_variant_counts.csv"), index=False)
    final_stats_df.to_csv(os.path.join(out_dir, f"{sample_name}_qc_summary.csv"), index=False)
    mismatch_df.write_csv(os.path.join(out_dir, f"{sample_name}_mismatch_variants.csv"))




index=0

harmonisation_diagnosis_pipeline(input_df,index)