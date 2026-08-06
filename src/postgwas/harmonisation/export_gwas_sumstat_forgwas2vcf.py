import os
import json
import polars as pl
import os, json, io
from pathlib import Path
from typing import Dict, Tuple,List



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


# import polars as pl
# import os
# import json
# import io
# from pathlib import Path
# from typing import Dict

# def export_gwas_sumstat(
#     df: pl.DataFrame,
#     sample_column_dict: Dict,
#     output_dir: str,
#     gwas_outputname: str,
#     chromosome: str,
#     genome_build: str = "GRCh37"
# ) -> pl.DataFrame:
#     """
#     Export GWAS summary stats (.tsv) and create a .dict JSON file.
#     Merged Version: Strict typing + Original variable/file naming integrity.
#     """

#     # -------------------------------------------------------
#     # Setup log file (per chromosome) - EXACT NAMES FROM OLD
#     # -------------------------------------------------------
#     safe_outputname = gwas_outputname if gwas_outputname not in [None, "", "NA"] else "GWAS"
#     log_root = Path(output_dir) / "logs"
#     log_root.mkdir(parents=True, exist_ok=True)
#     log_file = log_root / f"{safe_outputname}_chr{chromosome}_export.log"

#     log_buffer = io.StringIO()

#     def log_print(*args):
#         msg = " ".join(str(a) for a in args)
#         log_buffer.write(msg + "\n") 

#     log_print(f"\n🧩 Processing chromosome {chromosome} for {gwas_outputname}...")

#     # Column Lists - EXACT NAMES FROM OLD
#     required_columns_in_sumstat = [
#         'chr_col', 'pos_col', 'snp_id_col', 'ea_col', 'oa_col',
#         'eaf_col', 'beta_col', 'se_col', 'imp_z_col',
#         'pval_col', 'ncontrol_col'
#     ]
#     optional_columns_in_sumstat = ['ncase_col', 'imp_info_col']

#     # New Strict Type Mapping
#     type_mapping = {
#         'chr_col': pl.Utf8, 'pos_col': pl.Int64, 'snp_id_col': pl.Utf8,
#         'ea_col': pl.Utf8, 'oa_col': pl.Utf8, 'eaf_col': pl.Float64,
#         'beta_col': pl.Float64, 'se_col': pl.Float64, 'imp_z_col': pl.Float64,
#         'pval_col': pl.Float64, 'ncontrol_col': pl.Int64, 'ncase_col': pl.Int64,
#         'imp_info_col': pl.Float64
#     }

#     # -------------------------------------------------------
#     # Required/Optional Logic - RESTORED FROM OLD
#     # -------------------------------------------------------
#     required_pairs = {
#         k: sample_column_dict[k]
#         for k in required_columns_in_sumstat
#         if k in sample_column_dict and sample_column_dict[k] != 'NA'
#     }

#     missing = set(required_columns_in_sumstat) - set(required_pairs.keys())
#     if missing:
#         # Restored the console prints specifically for the missing column error
#         log_print(f"⚠️ Missing required columns: {', '.join(missing)} — exiting.")
#         print(f"⚠️ Missing required columns: {', '.join(missing)} — exiting.")
#         print(sample_column_dict)
        
#         with open(log_file, "w") as f:
#             f.write(log_buffer.getvalue())

#         return pl.DataFrame([{
#             "chromosome": chromosome,
#             "status": "missing_required_columns",
#             "missing_columns": ", ".join(missing)
#         }])

#     optional_pairs = {
#         k: sample_column_dict[k]
#         for k in optional_columns_in_sumstat
#         if k in sample_column_dict and sample_column_dict[k] != 'NA'
#     }

#     all_pairs = {**required_pairs, **optional_pairs}

#     # -------------------------------------------------------
#     # Select & Cast - INTEGRATED NEW TYPE SAFETY
#     # -------------------------------------------------------
#     selected_cols = list(all_pairs.values())
#     df2 = df.select(selected_cols)

#     # Apply proper data types based on mapping
#     cast_exprs = []
#     for key, col_name in all_pairs.items():
#         if key in type_mapping:
#             cast_exprs.append(pl.col(col_name).cast(type_mapping[key], strict=False))
    
#     df2 = df2.with_columns(cast_exprs)

#     # -------------------------------------------------------
#     # Filtering - RESTORED OLD VERBOSE LOGS & LOGIC
#     # -------------------------------------------------------
#     se_col_name = sample_column_dict['se_col']
#     beta_col_name = sample_column_dict['beta_col']

#     old_count = len(df2)
#     df2 = df2.filter(
#         (pl.col(se_col_name) > 0) &
#         (pl.col(se_col_name).is_not_null()) &
#         (pl.col(se_col_name).is_finite()) &
#         (pl.col(beta_col_name).is_not_null()) &
#         (pl.col(beta_col_name).is_finite())
#     )
#     new_count = len(df2)
#     n_dropped = old_count - new_count

#     if n_dropped > 0:
#         # Restored exact print string for console
#         print(f"\t\t\t\t Dropped {n_dropped} variants due to invalid '{se_col_name} (0,null of Infinity)' or '{beta_col_name} (null or Infinity)' values from chromosome{chromosome}.")

#     # -------------------------------------------------------
#     # Export TSV - EXACT FILENAME RESTORED
#     # -------------------------------------------------------
#     os.makedirs(output_dir, exist_ok=True)
#     tsv_path = f"{output_dir}/{gwas_outputname}_chr{chromosome}_vcf_input.tsv"
#     df2.write_csv(tsv_path, separator="\t")
#     log_print(f"💾 Saved subset to {tsv_path}")

#     # -------------------------------------------------------
#     # JSON Dict - RESTORED EXACT KEY LOGIC
#     # -------------------------------------------------------
#     df_columns = df2.columns
#     params = {}
#     for key, col in all_pairs.items():
#         if key == "snp_id_col":
#             params["snp_col"] = df_columns.index(col)
#         elif key == "beta_col":
#             params["beta_col"] = df_columns.index(col)
#         else:
#             params[key] = df_columns.index(col)

#     params.update({"delimiter": "\t", "header": True, "build": genome_build})

#     dict_path = f"{output_dir}/{gwas_outputname}_chr{chromosome}.dict"
#     with open(dict_path, "w") as outfile:
#         json.dump(params, outfile, indent=4)
#     log_print(f"🧾 Metadata dictionary saved → {dict_path}")

#     # -------------------------------------------------------
#     # Summary & Return - RESTORED ALL METADATA COLUMNS
#     # -------------------------------------------------------
#     # Note: summarize_gwas_dataframe must be available in scope
#     summary_df = summarize_gwas_dataframe(df2, sample_column_dict, chromosome)

#     summary_df = summary_df.with_columns([
#         pl.lit(tsv_path).alias("tsv_path"),
#         pl.lit(dict_path).alias("dict_path"),
#         pl.lit(df2.height).alias("num_rows"),
#         pl.lit(df2.width).alias("num_cols"),
#         pl.lit("success").alias("status")
#     ])

#     # Restored writing the summary file to disk
#     summary_path = f"{output_dir}/{gwas_outputname}_chr{chromosome}_summary.tsv"
#     summary_df.write_csv(summary_path, separator="\t")

#     log_print(f"📈 Summary saved → {summary_path}")
#     log_print(f"🎯 Completed export + summary for chr{chromosome}.\n")

#     # Finalize Logs
#     with open(log_file, "w") as f:
#         f.write(log_buffer.getvalue())

#     return summary_df

def export_gwas_sumstat(
    df: pl.DataFrame,
    sample_column_dict: Dict,
    output_dir: str,
    gwas_outputname: str,
    chromosome: str,
    genome_build: str = "GRCh37"
) -> pl.DataFrame:
    """
    Export GWAS summary stats (.tsv) and create a .dict JSON file
    using column **indices (0, 1, 2, ...)** instead of column names.

    Silent version:
        - NO print() to stdout
        - All log output written to per-chromosome logfile
        - Multiprocessing safe
    """

    # -------------------------------------------------------
    # Setup log file (per chromosome)
    # -------------------------------------------------------
    safe_outputname = gwas_outputname if gwas_outputname not in [None, "", "NA"] else "GWAS"

    log_root = Path(output_dir) / "logs"
    log_root.mkdir(parents=True, exist_ok=True)

    log_file = log_root / f"{safe_outputname}_chr{chromosome}_export.log"

    # Internal buffer — silent
    log_buffer = io.StringIO()

    def log_print(*args):
        msg = " ".join(str(a) for a in args)
        log_buffer.write(msg + "\n")   # *** no console output ***

    # -------------------------------------------------------
    # MAIN EXECUTION
    # -------------------------------------------------------
    log_print(f"\n🧩 Processing chromosome {chromosome} for {gwas_outputname}...")

    required_columns_in_sumstat = [
        'chr_col', 'pos_col', 'snp_id_col', 'ea_col', 'oa_col',
        'eaf_col', 'beta_col', 'se_col', 'imp_z_col',
        'pval_col', 'ncontrol_col'
    ]

    optional_columns_in_sumstat = ['ncase_col', 'imp_info_col']

    # -------------------------------------------------------
    # Required column mapping
    # -------------------------------------------------------
    required_pairs = {
        k: sample_column_dict[k]
        for k in required_columns_in_sumstat
        if k in sample_column_dict and sample_column_dict[k] != 'NA'
    }

    missing = set(required_columns_in_sumstat) - set(required_pairs.keys())
    if missing:
        log_print(f"⚠️ Missing required columns: {', '.join(missing)} — exiting.")
        print(f"⚠️ Missing required columns: {', '.join(missing)} — exiting.")
        print(sample_column_dict)
        # write log
        with open(log_file, "w") as f:
            f.write(log_buffer.getvalue())

        return pl.DataFrame([{
            "chromosome": chromosome,
            "status": "missing_required_columns",
            "missing_columns": ", ".join(missing)
        }])

    # -------------------------------------------------------
    # Optional columns
    # -------------------------------------------------------
    optional_pairs = {
        k: sample_column_dict[k]
        for k in optional_columns_in_sumstat
        if k in sample_column_dict and sample_column_dict[k] != 'NA'
    }

    all_pairs = {**required_pairs, **optional_pairs}

    # -------------------------------------------------------
    # Select columns
    # -------------------------------------------------------
    selected_cols = list(all_pairs.values())
    df2 = df.select(selected_cols)

    # Cast to string for safe writing
    df2 = df2.with_columns([pl.col(c).cast(pl.Utf8) for c in df2.columns])

    # Extract column names from your dictionary
    # 1. Define the desired data types for your logical columns
    type_mapping = {
        'se_col': pl.Float64,
        'beta_col': pl.Float64,
        'eaf_col': pl.Float64,
        'pval_col': pl.Float64,
        'imp_z_col': pl.Float64,
        'imp_info_col': pl.Float64,
        'chr_col': pl.Utf8,
        'pos_col': pl.Int64,
        'snp_id_col': pl.Utf8,
        'ea_col': pl.Utf8,
        'oa_col': pl.Utf8,
        'ncontrol_col': pl.Int64,
        'ncase_col': pl.Int64
    }

    # 2. Build a list of casting expressions safely
    cast_exprs = []
    
    for logical_key, dtype in type_mapping.items():
        # Check if the key exists in the dict
        if logical_key in sample_column_dict:
            actual_col_name = sample_column_dict[logical_key]
            
            # Check if the column is valid AND actually exists in df2
            if actual_col_name != 'NA' and actual_col_name in df2.columns:
                cast_exprs.append(pl.col(actual_col_name).cast(dtype, strict=False))

    # 3. Apply the casts all at once
    if cast_exprs:
        df2 = df2.with_columns(cast_exprs)





    # -------------------------------------------------------
    # Filtering - RESTORED OLD VERBOSE LOGS & LOGIC
    # -------------------------------------------------------
    se_col_name = sample_column_dict['se_col']
    beta_col_name = sample_column_dict['beta_col']

    old_count=len(df2)
    # Apply filtering
    df2 = df2.filter(
        # --- Standard Error Checks ---
        (pl.col(se_col_name) > 0) &                # Removes 0.0 and negatives (CRITICAL for FINEMAP)
        (pl.col(se_col_name).is_not_null()) &      # Removes NaNs/None
        (pl.col(se_col_name).is_finite()) &        # Removes Infinity (Inf)
        # --- Beta Checks ---
        (pl.col(beta_col_name).is_not_null()) &    # Removes NaNs/None
        (pl.col(beta_col_name).is_finite())       # Removes Infinity (Inf)
    )
    new_count=len(df2)

    # (Optional) Check how many were dropped
    n_dropped = old_count - new_count
    if n_dropped > 0:
        print(f"\t\t\t\t Dropped {n_dropped} variants due to invalid '{se_col_name} (0,null of Infinity)' or '{beta_col_name} (null or Infinity)' values from chromosome{chromosome}.")
        
    # -------------------------------------------------------
    # Prepare output folder
    # -------------------------------------------------------
    final_outdir = output_dir
    os.makedirs(final_outdir, exist_ok=True)

    # -------------------------------------------------------
    # Write TSV
    # -------------------------------------------------------
    tsv_path = f"{final_outdir}/{gwas_outputname}_chr{chromosome}_vcf_input.tsv"
    df2.write_csv(tsv_path, separator="\t")

    log_print(f"💾 Saved subset to {tsv_path}")

    # -------------------------------------------------------
    # Build JSON mapping using column indices
    # -------------------------------------------------------
    df_columns = df2.columns
    params = {}

    for key, col in all_pairs.items():
        if key == "snp_id_col":
            params["snp_col"] = df_columns.index(col)
        elif key == "beta_col":
            params["beta_col"] = df_columns.index(col)
        else:
            params[key] = df_columns.index(col)
    # # Fill missing optional fields
    # for opt_key in optional_columns_in_sumstat:
    #     if opt_key not in params:
    #         params[opt_key] = "NA"

    params.update({
        "delimiter": "\t",
        "header": True,
        "build": genome_build
    })

    # Write dict file
    dict_path = f"{final_outdir}/{gwas_outputname}_chr{chromosome}.dict"
    with open(dict_path, "w") as outfile:
        json.dump(params, outfile, indent=4)

    log_print(f"🧾 Metadata dictionary saved → {dict_path}")

    # -------------------------------------------------------
    # Build summary
    # -------------------------------------------------------
    summary_df = summarize_gwas_dataframe(df2, sample_column_dict, chromosome)

    summary_df = summary_df.with_columns([
        pl.lit(tsv_path).alias("tsv_path"),
        pl.lit(dict_path).alias("dict_path"),
        pl.lit(df2.height).alias("num_rows"),
        pl.lit(df2.width).alias("num_cols"),
        pl.lit("success").alias("status")
    ])

    summary_path = f"{final_outdir}/{gwas_outputname}_chr{chromosome}_summary.tsv"
    summary_df.write_csv(summary_path, separator="\t")

    log_print(f"📈 Summary saved → {summary_path}")
    log_print(f"🎯 Completed export + summary for chr{chromosome}.\n")

    # -------------------------------------------------------
    # Write log file
    # -------------------------------------------------------
    with open(log_file, "w") as f:
        f.write(log_buffer.getvalue())

    return summary_df
