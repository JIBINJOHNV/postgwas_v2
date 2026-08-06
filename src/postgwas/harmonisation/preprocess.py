import pandas as pd
import os
import polars as pl
from concurrent.futures import ProcessPoolExecutor

def split_chr_pos(df: pl.DataFrame, sample_gwas_dict: dict) -> list:
    """
    Split GWAS dataframe into per-chromosome files and return file paths for downstream parallel processing.

    Parameters
    ----------
    df : pl.DataFrame
        Full GWAS summary statistics dataframe.
    sample_gwas_dict : dict
        Dictionary with 'chr_col' and 'gwas_outputname'.
    output_folder : str
        Path to output folder.

    Returns
    -------
    list
        List of chromosome file paths for downstream analysis.
    """
    output_folder =sample_gwas_dict['output_folder']
    os.makedirs(output_folder, exist_ok=True)
    chr_col = sample_gwas_dict.get("chr_col")
    gwas_name = sample_gwas_dict.get("gwas_outputname", "gwas_output")
    
    if chr_col not in df.columns:
        raise KeyError(f"❌ Chromosome column '{chr_col}' not found in dataframe.")
    # Ensure consistent formatting
    df = df.with_columns(pl.col(chr_col).cast(pl.Utf8).str.strip_chars())
    #chromosomes = df.select(chr_col).unique().to_series().to_list()
    df = df.filter(pl.col(chr_col).is_not_null())
    chromosomes = (
        df.select(chr_col)
        .unique()
        .to_series()
        .to_list()
    )
    
    output_files = []
    
    for chrom in chromosomes:
        #subset = df.filter(pl.col(chr_col) == chrom)
        subset = df.filter(pl.col(chr_col).eq(chrom))
        out_path = os.path.join(output_folder, f"{gwas_name}_chr{chrom}.tsv")
        subset.write_csv(out_path, separator="\t")
        output_files.append(out_path)
        #print(f"✅ Saved chromosome {chrom} ({len(subset)} rows) → {out_path}")
    #print(f"🎉 Split complete: {len(output_files)} chromosome files saved.")
    return output_files
