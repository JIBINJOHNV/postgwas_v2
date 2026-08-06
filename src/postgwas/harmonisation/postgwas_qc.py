import pandas as pd
import json
from pathlib import Path

def flatten_dict(d, parent_key='', sep='_'):
    """Recursively flattens a dictionary."""
    items = []
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

def qc_json_to_dataframes(json_path):
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    rows = []
    global_metrics = {}

    for key, value in data.items():
        # Identify chromosome keys (numeric strings or X/Y/MT)
        if (key.isdigit() or key in ["X", "Y", "MT"]) and isinstance(value, dict):
            for section_name, section_content in value.items():
                if isinstance(section_content, dict):
                    # Flatten deep structures (like sample_size_stats_N_min)
                    flat_section = flatten_dict(section_content)
                    for metric_name, val in flat_section.items():
                        rows.append({
                            "section": section_name,
                            "metric": metric_name,
                            "chromosome": f"chr{key}",
                            "value": val
                        })
                else:
                    # Handle top-level values within a chromosome (e.g., exit_code)
                    rows.append({
                        "section": "summary",
                        "metric": section_name,
                        "chromosome": f"chr{key}",
                        "value": section_content
                    })
        else:
            # Root level global values
            global_metrics[key] = value

    # Create the main DataFrame from chromosome data
    if not rows:
        return pd.DataFrame()

    df_long = pd.DataFrame(rows)
    
    # Pivot so chromosomes are columns
    df_pivot = df_long.pivot(
        index=["section", "metric"], 
        columns="chromosome", 
        values="value"
    ).reset_index()

    # Create Global rows and append them
    global_rows = []
    for k, v in global_metrics.items():
        global_rows.append({
            "section": "GLOBAL",
            "metric": k,
            "GLOBAL_VALUE": v
        })
    
    if global_rows:
        df_global = pd.DataFrame(global_rows)
        # Combine the pivoted data with the global summary rows
        df_final = pd.concat([df_pivot, df_global], axis=0, ignore_index=True)
        return df_final

    return df_pivot


# def qc_json_to_dataframes(qc_file):
#     qc_file = Path(qc_file)
#     with qc_file.open() as f:
#         qc = json.load(f)
#     # Capture totals
#     total_variant_infile = qc.get("total_variant_infile", None)
#     total_variant_read   = qc.get("total_variant_read", None)
#     chrom_keys = [k for k in qc.keys() if k.isdigit()]
#     long_records = []
#     for chrom in chrom_keys:
#         chrom_data = qc[chrom]
#         for section_name, section_data in chrom_data.items():
#             # --- SPECIAL CASE: sample_size_qc ---
#             if section_name == "sample_size_qc":
#                 # 1. Handle simple metrics
#                 if isinstance(section_data, dict):
#                     for k, v in section_data.items():
#                         if k != "sample_size_stats":
#                             long_records.append({
#                                 "chromosome": chrom,
#                                 "section": section_name,
#                                 "metric": k,
#                                 "value": v
#                             })
#                     # 2. Handle nested stats safely
#                     stats = section_data.get("sample_size_stats", {})
#                     # FIX: Check if 'stats' is actually a dictionary
#                     if isinstance(stats, dict):
#                         for group_name, group_values in stats.items():
#                             if isinstance(group_values, dict):
#                                 for stat_name, stat_value in group_values.items():
#                                     long_records.append({
#                                         "chromosome": chrom,
#                                         "section": section_name,
#                                         "metric": f"{group_name}_{stat_name}",
#                                         "value": stat_value
#                                     })
#                     else:
#                         print(f"⚠️ Warning: 'sample_size_stats' in Chr{chrom} is not a dict. Value: {stats}")
#             # --- NORMAL SECTIONS ---
#             elif isinstance(section_data, dict):
#                 for k, v in section_data.items():
#                     if isinstance(v, dict):
#                         for k2, v2 in v.items():
#                             long_records.append({
#                                 "chromosome": chrom,
#                                 "section": section_name,
#                                 "metric": f"{k}_{k2}",
#                                 "value": v2
#                             })
#                     else:
#                         long_records.append({
#                             "chromosome": chrom,
#                             "section": section_name,
#                             "metric": k,
#                             "value": v
#                         })
#     # If no data found, return empty DF
#     if not long_records:
#         return pd.DataFrame()
#     qc_long = pd.DataFrame(long_records)
#     # Pivot
#     df = qc_long.pivot(
#         index=["section", "metric"],
#         columns="chromosome",
#         values="value"
#     ).reset_index()
#     # Clean columns
#     df.columns = [
#         "_".join(col).strip() if isinstance(col, tuple) else col
#         for col in df.columns
#     ]
#     df.columns = [c.replace("value_", "") for c in df.columns]
#     # Add 'chr' prefix for numeric columns
#     df.columns = [f"chr{c}" if str(c).isdigit() else c for c in df.columns]
#     # Add totals
#     df["total_variant_infile"] = total_variant_infile
#     df["total_variant_read"]   = total_variant_read
#     return df

# def qc_json_to_dataframes(qc_file):
#     """
#     Convert QC summary JSON into a single flattened dataframe:
#         rows    = (section, metric)
#         columns = chr1, chr2, …, chr22 + total_variant_infile, total_variant_read
#     Handles nested:
#         sample_size_qc → sample_size_stats → Nca/Nco/Neff
#     """
#     qc_file = Path(qc_file)
#     # Load JSON
#     with qc_file.open() as f:
#         qc = json.load(f)
#     # ---- NEW: capture top-level totals (if present) ----
#     total_variant_infile = qc.get("total_variant_infile", None)
#     total_variant_read   = qc.get("total_variant_read", None)
#     # Identify chromosome keys
#     chrom_keys = [k for k in qc.keys() if k.isdigit()]
#     long_records = []
#     # -------------------------------------------
#     # Build long-format table
#     # -------------------------------------------
#     for chrom in chrom_keys:
#         chrom_data = qc[chrom]
#         for section_name, section_data in chrom_data.items():
#             # SPECIAL CASE: sample_size_qc
#             if section_name == "sample_size_qc":
#                 # top-level simple metrics
#                 for k, v in section_data.items():
#                     if k != "sample_size_stats":
#                         long_records.append({
#                             "chromosome": chrom,
#                             "section": section_name,
#                             "metric": k,
#                             "value": v
#                         })
#                 # nested stats: Nca / Nco / Neff
#                 stats = section_data.get("sample_size_stats", {})
#                 for group_name, group_values in stats.items():
#                     for stat_name, stat_value in group_values.items():
#                         long_records.append({
#                             "chromosome": chrom,
#                             "section": section_name,
#                             "metric": f"{group_name}_{stat_name}",
#                             "value": stat_value
#                         })
#             # NORMAL SECTIONS
#             elif isinstance(section_data, dict):
#                 for k, v in section_data.items():
#                     if isinstance(v, dict):
#                         # flatten one nested dict level
#                         for k2, v2 in v.items():
#                             long_records.append({
#                                 "chromosome": chrom,
#                                 "section": section_name,
#                                 "metric": f"{k}_{k2}",
#                                 "value": v2
#                             })
#                     else:
#                         long_records.append({
#                             "chromosome": chrom,
#                             "section": section_name,
#                             "metric": k,
#                             "value": v
#                         })
#     qc_long = pd.DataFrame(long_records)
#     # -------------------------------------------
#     # Create wide pivot table
#     # -------------------------------------------
#     df = qc_long.pivot(
#         index=["section", "metric"],
#         columns="chromosome",
#         values="value"
#     ).reset_index()
#     # Flatten MultiIndex columns
#     df.columns = [
#         "_".join(col).strip() if isinstance(col, tuple) else col
#         for col in df.columns
#     ]
#     # Remove "value_" prefix (from pivot)
#     df.columns = [c.replace("value_", "") for c in df.columns]
#     # Add chr prefix for numeric chromosome columns
#     df.columns = [
#         f"chr{c}" if c.isdigit() else c
#         for c in df.columns
#     ]
#     # ---- NEW: add total variant columns (same value for all rows) ----
#     df["total_variant_infile"] = total_variant_infile
#     df["total_variant_read"]   = total_variant_read
#     return df




# def postgwas_qc(qc_results,sample_column_dict) -> pd.DataFrame:
#     """
#     Extracts EVERY useful QC metric from post-GWAS (MungeSumstats-like) report.
#     Returns a beautiful, publication-ready table with genome-wide row.
#     """
#     sample_name = sample_column_dict["gwas_outputname"]

#     output_folder = (
#         Path(sample_column_dict["output_folder"]) 
#         / "qc_summary" )
#     output_folder.mkdir(parents=True, exist_ok=True)
#     records = []
#     for chr_str in [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']:  # in case non-autosomal later
#         if chr_str not in qc_results:
#             continue
#         c = qc_results[chr_str]
#         # Helper to safely extract nested values
#         eaf = c['eaf_qc']
#         info = c['info_qc']
#         beta = c['or_qc']
#         pval = c['pval_qc']
#         z = c['z_qc']
#         record = {
#             "CHR": int(chr_str),
#             "N_Variants": eaf['final_total_variants'],
#             # === Allele Frequency ===
#             "Mean_EAF": round(eaf['final_mean_eaf'], 5),
#             "Median_EAF": None,  # not directly given, but mean ≈ median here
#             "Min_EAF": eaf['final_min_eaf'],
#             # === Imputation Quality ===
#             "Mean_INFO": round(info['mean_info'], 4),
#             "Median_INFO": info.get('median_info', None),  # exists in your data!
#             "Min_INFO": info['min_info'],
#             # === Effect Sizes (Beta) ===
#             "Mean_BETA": round(beta['mean_beta'], 6),
#             "SD_BETA": round(beta['std_beta'], 5),
#             "Min_BETA": round(beta['min_beta'], 5),
#             "Max_BETA": round(beta['max_beta'], 5),
#             "%_Negative_BETA": round(beta['negative_fraction'], 4),
#             # === Z-scores ===
#             "Mean_Z": round(z['z_mean'], 4),
#             "SD_Z": round(z['z_std'], 4),
#             "Min_Z": round(z['z_min'], 2),
#             "Max_Z": round(z['z_max'], 2),
#             # === P-values ===
#             "Median_P": pval['median'],
#             "Max_P": pval['max'],
#             # === Extra flags ===
#             "INFO_source": info['info_source'],
#             "Effect_Type": beta['effect_type'].upper(),
#         }
        
#         # Extract median_info if present (it is!)
#         if 'median_info' in info and info['median_info'] is not None:
#             record["Median_INFO"] = round(info['median_info'], 4)
#         records.append(record)
#     df = pd.DataFrame(records)
#     df = df.sort_values("CHR").reset_index(drop=True)
#     # === Genome-wide summary row ===
#     total_vars = qc_results['total_variant_read']
#     gw = {
#         "CHR": "Genome-wide",
#         "N_Variants": total_vars,
#         "Mean_EAF": round(df["Mean_EAF"].mean(), 5),
#         "Mean_INFO": round(df["Mean_INFO"].mean(), 4),
#         "Median_INFO": round(df["Median_INFO"].mean(), 4) if "Median_INFO" in df else None,
#         "Min_INFO": round(df["Min_INFO"].min(), 4),
#         "Mean_BETA": round(df["Mean_BETA"].mean(), 6),
#         "SD_BETA": round(df["SD_BETA"].mean(), 5),
#         "%_Negative_BETA": round(df["%_Negative_BETA"].mean(), 4),
#         "Mean_Z": round(df["Mean_Z"].mean(), 4),
#         "SD_Z": round(df["SD_Z"].mean(), 4),
#         "Min_Z": round(df["Min_Z"].min(), 2),
#         "Max_Z": round(df["Max_Z"].max(), 2),
#         "Median_P": round(df["Median_P"].median(), 5),
#         "Effect_Type": "BETA",
#     }
#     gw.update({k: None for k in ["Min_EAF", "Min_BETA", "Max_BETA", "Max_P", "INFO_source"]})
#     df_gw = pd.DataFrame([gw])
#     df_full = pd.concat([df, df_gw], ignore_index=True)
#     df_full['total_variants_in_the_input_files'] = qc_results['total_variant_infile']-1 # Sort columns alphabetically
#     df_full['total_variant_read_by_python'] = qc_results['total_variant_read']# Sort columns alphabetically
#     df_full.to_csv(f"{output_folder}/{sample_name}_gwas2vcf_summary.tsv", sep="\t", index=False)
#     return df_full
