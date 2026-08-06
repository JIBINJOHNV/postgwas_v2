import os
from typing import Dict, List, Optional
import polars as pl
from pathlib import Path


def summarize_gwas_dataframe(df: pl.DataFrame, sample_column_dict: Dict, chr: str) -> pl.DataFrame:
    """Generate Polars DataFrame summary for all columns."""
    summary_records = []
    for key, col_name in sample_column_dict.items():
        if col_name == "NA" or col_name not in df.columns:
            continue
        dtype = df.schema[col_name]
        n_missing = df.select(pl.col(col_name).is_null().sum().alias("n_missing")).item()
        record = {
            "chromosome": chr,
            "key": key,
            "column_name": col_name,
            "dtype": str(dtype),
            "n_missing": n_missing,
            "n_unique": None,
            "min": None,
            "max": None,
            "mean": None,
            "median": None,
            "std": None
        }
        if dtype in (pl.Float32, pl.Float64, pl.Int32, pl.Int64, pl.UInt32, pl.UInt64):
            stats = df.select([
                pl.col(col_name).min().alias("min"),
                pl.col(col_name).max().alias("max"),
                pl.col(col_name).mean().alias("mean"),
                pl.col(col_name).median().alias("median"),
                pl.col(col_name).std().alias("std"),
                pl.col(col_name).n_unique().alias("n_unique")
            ]).to_dicts()[0]
            record.update(stats)
        elif dtype == pl.Utf8:
            record["n_unique"] = df.select(pl.col(col_name).n_unique()).item()
        summary_records.append(record)
    return pl.DataFrame(summary_records)


def export_gwas_sumstat(
    df: pl.DataFrame,
    sample_column_dict: Dict[str, str],
    output_dir: str,
    gwas_outputname: str,
    chr: str,
    required_columns_in_sumstat: Optional[List[str]] = None,
    optional_columns_in_sumstat: Optional[List[str]] = None,
) -> pl.DataFrame:
    """
    Export GWAS summary statistics to TSV and generate a column mapping dictionary using **column indices**.

    This function:
    1. Maps dataset-specific columns to standard GWAS format.
    2. Validates required columns.
    3. Exports:
       - Main summary stats TSV (for downstream tools like munging)
       - Column mapping TSV (for reproducibility)
       - Summary metadata TSV
    4. Returns a status DataFrame.

    Args:
        df: Input Polars DataFrame with GWAS results.
        sample_column_dict: Dictionary mapping standard names to actual column names in `df`.
            Expected keys: 'chr_col', 'pos_col', 'ea_col', 'oa_col', 'snp_id_col', 'beta_col',
            'se_col', 'imp_z_col', 'pval_col', 'imp_info_col', 'ncase_col', 'ncontrol_col', 'neff_col'
        output_dir: Base output directory.
        gwas_outputname: Prefix for output files.
        chr: Chromosome identifier (e.g., '1', 'X').
        required_columns_in_sumstat: List of required standard columns (defaults to common set).
        optional_columns_in_sumstat: List of optional standard columns.

    Returns:
        Polars DataFrame with export summary (one row), including status and paths.
    """
    # === Default column sets ===
    if required_columns_in_sumstat is None:
        required_columns_in_sumstat = [
            "CHR", "BP", "A1", "A2", "SNP", "BETA", "SE", "Z", "LP",'FRQ', "N_CON", "NEFF"
        ]
    if optional_columns_in_sumstat is None:
        optional_columns_in_sumstat = ["N_CAS", "INFO"]
    
    standard_columns = [
        "CHR", "BP", "A1", "A2", "SNP", "BETA", "SE", "Z", "LP", 'FRQ',"INFO", "N_CAS", "N_CON", "NEFF" ]
    
    # === Column mappings from input dict ===
    mapping_keys = [
        'chr_col', 'pos_col', 'ea_col', 'oa_col', 'snp_id_col', 'beta_col', 'se_col','imp_z_col', 'pval_col','eaf_col', 
        'imp_info_col', 'ncase_col', 'ncontrol_col', 'neff_col' ]
    
    mapped_columns = [ sample_column_dict.get(key, "NA") for key in mapping_keys ]
    
    # === Build mapping DataFrame (standard → source column) ===
    mapping_data = [
        {"mapped_column": std, "standard_column": src}
        for std, src in zip(mapped_columns,standard_columns)
        if std != "NA"
    ]
    
    mapped_df = pl.DataFrame(mapping_data)
    # === Validate required columns ===
    available_std_cols = set(mapped_df["standard_column"].to_list())
    missing_required = set(required_columns_in_sumstat) - available_std_cols
    
    if missing_required:
        print(f"⚠️ Missing required columns: {', '.join(sorted(missing_required))} — export failed.")
        return pl.DataFrame([{
            "chromosome": chr,
            "status": "missing_required_columns",
            "missing_columns": ", ".join(sorted(missing_required)),
            "column_map_path": None,
            "num_rows": 0,
            "num_cols": 0
        }])
    
    # === Setup output paths ===
    output_path = Path(output_dir) / gwas_outputname
    output_path.mkdir(parents=True, exist_ok=True)
    column_map_tsv = output_path / f"{gwas_outputname}_chr{chr}_column_mapping.tsv"
    gwas_output_tsv = output_path / f"{gwas_outputname}_chr{chr}_vcf_input.tsv"
    summary_tsv = output_path / f"{gwas_outputname}_chr{chr}_summary.tsv"
    
    # === Select and reorder columns based on mapping ===
    selected_cols = mapped_df["mapped_column"].to_list()
    df_export = df.select(selected_cols)
    # === Write outputs ===
    # 1. Column mapping (for documentation)
    mapped_df.write_csv(column_map_tsv, separator="\t")
    # 2. Main GWAS summary stats
    df_export.write_csv(gwas_output_tsv, separator="\t")
    # === Generate summary ===
    summary_df = summarize_gwas_dataframe(df_export, sample_column_dict, chr)
    summary_df = summary_df.with_columns([
        pl.lit(str(column_map_tsv)).alias("column_map_path"),
        pl.lit(df_export.height).alias("num_rows"),
        pl.lit(df_export.width).alias("num_cols"),
        pl.lit("success").alias("status")
    ])
    # Write summary
    summary_df.write_csv(summary_tsv, separator="\t")
    # === Logging ===
    print(f"✅ GWAS summary stats exported → {gwas_output_tsv}")
    print(f"📋 Column mapping saved → {column_map_tsv}")
    print(f"📈 Processing summary saved → {summary_tsv}")
    print(f"🎯 Successfully processed chr{chr} ({df_export.height:,} variants).")
    return summary_df

