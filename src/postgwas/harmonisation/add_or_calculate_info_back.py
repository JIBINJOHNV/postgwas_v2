import polars as pl
import io, os
from pathlib import Path
from typing import Tuple, Dict

def add_or_calculate_info(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    info_file: str = "NA",
    info_column: str = "info",
    default_info_file: str = "NA",
    default_info_column: str = "info"
) -> Tuple[pl.DataFrame, Dict, dict]:
    """
    Add or merge imputation INFO score column into GWAS summary dataframe using Polars.
    """
    #print(f"--- Info score calculation started for Chr {chromosome} ---")

    # -------------------------------------------------------
    # Setup chromosome-specific log file
    # -------------------------------------------------------
    gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
    output_dir = sample_column_dict.get("output_folder", ".")
    log_dir = Path(output_dir) / "logs"
    
    # Create log directory immediately
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_add_or_calculate_info.log"

    # Initialize the log file (overwrite if exists)
    with open(log_file, "w") as f:
        f.write(f"Log started for Chromosome {chromosome}\n")

    def log_print(*args):
        """Writes to log file IMMEDIATELY so you can see progress."""
        msg = " ".join(str(a) for a in args)
        # Append to file on disk immediately
        with open(log_file, "a") as f:
            f.write(msg + "\n")

    # -------------------------------------------------------
    # MAIN LOGIC
    # -------------------------------------------------------
    qc_info = {"initial_variants": df.height}

    log_print(f"\n🧩 Starting INFO score harmonization for Chr {chromosome}...")

    # Extract mapping
    chr_col = sample_column_dict["chr_col"]
    pos_col = sample_column_dict["pos_col"]
    ea_col = sample_column_dict["ea_col"]
    oa_col = sample_column_dict["oa_col"]
    imp_info_col = sample_column_dict.get("imp_info_col", "NA")

    # ===============================================================
    # STEP 1: Use internal INFO column if present
    # ===============================================================
    if imp_info_col != "NA" and imp_info_col in df.columns:
        log_print(f"🔍 Found internal INFO column '{imp_info_col}' — validating numeric values.")

        df = df.with_columns(pl.col(imp_info_col).cast(pl.Float64, strict=False))

        invalid_before = df.filter(
            (pl.col(imp_info_col) < 0.0) | (pl.col(imp_info_col) > 1.05)
        ).height

        df = df.with_columns(pl.col(imp_info_col).clip(0.0, 1.0).alias(imp_info_col))

        stats = df.select([
            pl.col(imp_info_col).min().alias("min"),
            pl.col(imp_info_col).max().alias("max"),
            pl.col(imp_info_col).mean().alias("mean"),
            pl.col(imp_info_col).median().alias("median"),
            pl.len().alias("total")
        ]).to_dicts()[0]

        qc_info.update({
            "info_source": "internal",
            "invalid_values_clipped": invalid_before,
            "min_info": stats["min"],
            "max_info": stats["max"],
            "mean_info": stats["mean"],
            "median_info": stats["median"],
            "variants_before": qc_info["initial_variants"],
            "variants_after": stats["total"],
        })

        log_print(
            f"✅ INFO (internal): range={stats['min']:.3f}–{stats['max']:.3f}, "
            f"mean={stats['mean']:.3f}, median={stats['median']:.3f}, n={stats['total']:,}"
        )
        return df, qc_info, sample_column_dict

    # ===============================================================
    # STEP 2: Choose external INFO source
    # ===============================================================
    info_file_to_use = None
    info_col_to_use = None
    info_source = "NA"
    
    potential_chr_file = f"{info_file}_chr{chromosome}.tsv.gz"

    if info_file != "NA" and os.path.exists(info_file):
        info_file_to_use = info_file
        info_col_to_use = info_column
        info_source = "user provided"
    elif info_file != "NA" and os.path.exists(potential_chr_file):
        info_file_to_use = potential_chr_file
        info_col_to_use = info_column
        info_source = "user provided (chrom split)"
        log_print(f"ℹ️ Found chromosome split file: {info_file_to_use}")
    elif default_info_file != "NA" and os.path.exists(default_info_file):
        info_file_to_use = default_info_file
        info_col_to_use = default_info_column
        info_source = "default"
        log_print(f"⚠️ Primary INFO file missing. Using Default: {default_info_file}")
    else:
        err_msg = (
            f"\n🛑 PIPELINE STOPPED: INFO score not available.\n"
            f"• Chromosome: {chromosome}\n"
            f"• Checked user file: {info_file}\n"
            f"• Checked chrom split: {potential_chr_file}\n"
            f"• Checked default: {default_info_file}\n"
        )
        log_print(err_msg)
        raise FileNotFoundError(err_msg)

    log_print(
        f"📂 Using {info_source} INFO file: {os.path.basename(info_file_to_use)} "
        f"(column: '{info_col_to_use}')"
    )

    # ===============================================================
    # Helper: Load external INFO file
    # ===============================================================
    def load_info_file(path: str, colname: str) -> pl.DataFrame:
        log_print("⏳ Loading INFO file from disk... (This may take a moment)")
        if os.path.getsize(path) == 0:
            raise ValueError(f"File is empty: {path}")

        info_df = (
            pl.read_csv(
                path,
                separator="\t",
                null_values=["NA", ".", "-"]
            )
            .select(["CHROM", "POS", "ALT", "REF", colname])
            .with_columns([
                pl.col("POS").cast(pl.Int64),
                pl.col("CHROM")
                .cast(pl.String)
                .str.replace("^chr", "", literal=False),
                pl.col("ALT").str.to_uppercase(),
                pl.col("REF").str.to_uppercase(),
                pl.col(colname)
                .cast(pl.Float64, strict=False)
                .clip(0.0, 1.0),
            ])
            .rename({
                "CHROM": chr_col,
                "POS": pos_col,
                "ALT": ea_col,
                "REF": oa_col,
            })
        )

        log_print(f"✅ INFO file loaded → {info_df.height:,} variants")
        return info_df


    # ===============================================================
    # STEP 3: Merge INFO
    # ===============================================================
    info_df = load_info_file(info_file_to_use, info_col_to_use)
    #print(f"📊 INFO file loaded: {info_df.height:,} variants")
    log_print("⏳ Merging INFO scores with GWAS data...")
    
    # Normalize alleles in input df
    df = df.with_columns([
        pl.col(ea_col).cast(pl.String).str.to_uppercase(),
        pl.col(oa_col).cast(pl.String).str.to_uppercase(),
    ])

    # 3.1 Direct match
    df1 = df.join(info_df, on=[chr_col, pos_col, ea_col, oa_col], how="inner")

    # 3.2 Flipped match (Rename columns, DO NOT calculate 1-info)
    info_df_flip = info_df.rename({ea_col: oa_col, oa_col: ea_col})
    df2 = df.join(info_df_flip, on=[chr_col, pos_col, ea_col, oa_col], how="inner")

    # 3.3 Combine
    df_merged = pl.concat([df1, df2], how="vertical").unique(
        subset=[chr_col, pos_col, ea_col, oa_col]
    )

    #print(f"Total varinat after merging with infor file is : {df_merged.height:,} variants")
    # ===============================================================
    # STEP 4: Statistics
    # ===============================================================
    if df_merged.height == 0:
        log_print("⚠️ WARNING: 0 variants matched between GWAS and INFO file!")
        # Return original DF to avoid crashing, but warn
        return df, qc_info, sample_column_dict

    stats = df_merged.select([
        pl.col(info_col_to_use).min().alias("min"),
        pl.col(info_col_to_use).max().alias("max"),
        pl.col(info_col_to_use).mean().alias("mean"),
        pl.col(info_col_to_use).median().alias("median"),
        pl.len().alias("total")
    ]).to_dicts()[0]

    qc_info.update({
        "info_source": info_source,
        "info_file_used": info_file_to_use,
        "info_column_used": info_col_to_use,
        "min_info": stats["min"],
        "max_info": stats["max"],
        "mean_info": stats["mean"],
        "median_info": stats["median"],
        "variants_before": qc_info["initial_variants"],
        "variants_after": stats["total"],
    })

    sample_column_dict["imp_info_col"] = info_col_to_use

    log_print(
        f"✅ INFO merged ({info_source}): {stats['total']:,} variants "
        f"({stats['min']:.3f}–{stats['max']:.3f}, mean={stats['mean']:.3f}, median={stats['median']:.3f})"
    )
    log_print("🎯 INFO harmonization complete.\n")

    return df_merged, qc_info, sample_column_dict