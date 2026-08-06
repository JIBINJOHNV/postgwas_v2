import os,sys
import polars as pl
from pathlib import Path
from typing import Tuple, Dict, Optional, Callable
import sys
from datetime import datetime

# -------------------------------------------------------
# 1. Helper: Setup Logging
# -------------------------------------------------------
def setup_logging(output_dir: str, gwas_outputname: str, chromosome: str) -> Tuple[Path, Callable]:
    log_dir = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    
    # Handle chromosome suffix cleanly
    chr_suffix = f"_chr{chromosome}" if chromosome else ""
    log_file = log_dir / f"{gwas_outputname}{chr_suffix}_eaf.log"

    # Initialize log file with Start Time
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "w") as f:
        f.write(f"--- Log started: {gwas_outputname} (Chr {chromosome}) at {start_time} ---\n")

    def log_print(*args, to_screen: bool = False):
        """
        Writes message with timestamp to log file.
        If to_screen=True, prints the raw message to the console.
        """
        # 1. Prepare Timestamp
        now = datetime.now().strftime("%H:%M:%S")
        msg_content = " ".join(str(a) for a in args)
        log_entry = f"[{now}] {msg_content}"

        # 2. Write to File (Safely)
        try:
            with open(log_file, "a") as f:
                f.write(log_entry + "\n")
        except Exception as e:
            # If logging fails (e.g. disk full), print to stderr so we know
            print(f"⚠️ LOGGING FAILED: {e}", file=sys.stderr)
        
        # 3. Print to Screen (Optional)
        if to_screen:
            print(msg_content) # Print clean message (without timestamp) to console

    return log_file, log_print

# -------------------------------------------------------
# 2. Helper: Check EAF Stats (Range & MAF behavior)
# -------------------------------------------------------
def check_eaf_statistics(
    df: pl.DataFrame, 
    colname: str, 
    log_print: Callable,
    maf_cutoff: float = 0.95
) -> Tuple[bool, bool]:
    """
    Returns (is_valid_range, is_suspicious_maf)
    """
    if colname not in df.columns: return False, False
    total = df.height
    if total == 0: return False, False

    # Check 1: Range 0-1
    out_of_range = df.filter((pl.col(colname) < 0.0) | (pl.col(colname) > 1.05)).height
    if (out_of_range / total) > 0.05:
        error_msg = (
            f"❌ CRITICAL ERROR: During '{colname}' check, observed too many invalid variants "
            f"with Allele Frequency > 1.0 or < 0.0.\n"
            f"   Found {out_of_range} bad rows out of {total} ({out_of_range/total:.1%})."
        )
        # 2. Log it (optional, but good practice)
        log_print(error_msg, to_screen=True)
        # 3. STOP the process completely
        sys.exit(1)
    df = df.with_columns(pl.col(colname).clip(0.0, 1.0).alias(colname))

    #Remove rows where the column is null (missing)
    df = df.filter(pl.col(colname).is_not_null())
    #Calculate how many were lost
    n_missing = total - df.height
    #Check if we lost too many (> 5%)
    if total > 0 and (n_missing / total) > 0.05:
        error_msg = (
            f"❌ CRITICAL ERROR: During '{colname}' check, observed too many missing (null) values.\n"
            f"   Found {n_missing} missing rows out of {total} ({n_missing/total:.1%})."
        )
        # Log and Exit
        log_print(error_msg, to_screen=True)
        sys.exit(1)
    
    # Check 2: MAF Behavior
    low_freq = df.filter(pl.col(colname) <= 0.5).height
    low_freq_pct = low_freq / total
    
    is_suspicious = False
    if low_freq_pct > maf_cutoff:
        log_print(f"\t\t\t ⚠️  {colname}: {low_freq_pct*100:.2f}% values ≤ 0.5 (Looks like MAF).",)
        is_suspicious = True
        
    return True, is_suspicious,out_of_range,n_missing,df,low_freq_pct

# -------------------------------------------------------
# 3. Helper: Cross-Check against Default Reference
# -------------------------------------------------------
def check_if_data_is_actually_maf(
    df: pl.DataFrame,
    default_path: str,
    default_af_col: str,
    target_eaf_col: str,
    std_cols: dict,
    chromosome: str,
    log_print: Callable
) -> bool:
    """
    Returns TRUE if the data is confirmed to be MAF (Error).
    Returns FALSE if the data is confirmed to be Rare (Safe).
    """
    log_print(f"🔍 Validating against Default Reference to distinguish MAF vs Rare...")

    # Handle split chromosomes
    load_path = default_path
    if not os.path.exists(load_path) and os.path.exists(f"{default_path}_chr{chromosome}.tsv.gz"):
        load_path = f"{default_path}_chr{chromosome}.tsv.gz"
        
    if not os.path.exists(load_path):
        log_print(f"⚠️ Default Reference not found: {load_path}. Cannot verify.", to_screen=True)
        return True # Fail safe: Assume error if we can't verify suspicious data

    # Load Default Ref (Standard Columns: CHROM, POS, REF, ALT, <AF_Col>)
    try:
        # Try Tab then Space
        try: ref_df = pl.read_csv(load_path, separator="\t", has_header=True)
        except: ref_df = pl.read_csv(load_path, separator=" ", has_header=True)
            
        req_cols = ["CHROM", "POS", "REF", "ALT", default_af_col]
        if not all(c in ref_df.columns for c in req_cols):
            log_print(f"❌ Default file missing columns {req_cols}", to_screen=True)
            return True # Fail safe
            
        # Cast Types
        ref_df = ref_df.with_columns([
            pl.col("CHROM").cast(pl.Utf8).str.replace("0X", "X"),
            pl.col("POS").cast(pl.Int64),
            pl.col("REF").str.to_uppercase(),
            pl.col("ALT").str.to_uppercase(),
            pl.col(default_af_col).cast(pl.Float64)
        ])
    except Exception as e:
        log_print(f"❌ Error loading default ref: {e}", to_screen=True)
        return True

    # Prepare Internal DF
    df = df.with_columns([
        pl.col(std_cols["ea"]).str.to_uppercase(),
        pl.col(std_cols["oa"]).str.to_uppercase()
    ])

    # MERGE: Default(ALT/REF) == Data(EA/OA) -> Use AF
    m1 = df.join(
        ref_df,
        left_on=[std_cols["chr"], std_cols["pos"], std_cols["ea"], std_cols["oa"]],
        right_on=["CHROM", "POS", "ALT", "REF"],
        how="inner"
    ).select([pl.col(target_eaf_col).alias("TEST_EAF"), pl.col(default_af_col).alias("REF_EAF")])

    # MERGE: Default(REF/ALT) == Data(EA/OA) -> Use 1 - AF
    m2 = df.join(
        ref_df,
        left_on=[std_cols["chr"], std_cols["pos"], std_cols["ea"], std_cols["oa"]],
        right_on=["CHROM", "POS", "REF", "ALT"],
        how="inner"
    ).select([pl.col(target_eaf_col).alias("TEST_EAF"), (1 - pl.col(default_af_col)).alias("REF_EAF")])

    matched_df = pl.concat([m1, m2])
    
    if matched_df.height < 5:
        log_print("⚠️ Not enough matches with Default Reference.", to_screen=True)
        return True # Fail safe

    # LOGIC:
    # We look at variants where OUR data says EAF <= 0.5 (Suspiciously Low)
    suspicious_subset = matched_df.filter(pl.col("TEST_EAF") <= 0.5)
    total_sus = suspicious_subset.height
    
    # Check what Reference says:
    # If Reference says > 0.6, it means our data is wrong (It's showing MAF for a Common variant)
    confirmed_common_in_ref = suspicious_subset.filter(pl.col("REF_EAF") > 0.60).height
    
    error_rate = confirmed_common_in_ref / total_sus if total_sus > 0 else 0
    
    log_print(f"   • Overlap variants: {matched_df.height}")
    log_print(f"   • Variants with Input ≤ 0.5: {total_sus}")
    log_print(f"   • Of those, Default Ref says > 0.6: {confirmed_common_in_ref} ({error_rate*100:.1f}%)")

    if error_rate > 0.50:
        log_print("\t\t\t🛑 RESULT: Default Reference DISAGREES. (Input is MAF, Ref is Common).", to_screen=True)
        return True  # Yes, it IS a MAF error
    else:
        log_print("\t\t\t✔ RESULT: Default Reference AGREES. (Input is correct, variants are Rare).", to_screen=True)
        return False # No, it is NOT an error

# -------------------------------------------------------
# 4. Helper: Load Generic External
# -------------------------------------------------------
def load_and_merge_external(
    df: pl.DataFrame, path: str, eaf_col: str, colmap: dict, std_cols: dict, log_print: Callable, chromosome: str
) -> Tuple[pl.DataFrame, str]:
    
    load_path = path
    if not os.path.exists(load_path) and os.path.exists(f"{path}_chr{chromosome}.tsv.gz"):
        load_path = f"{path}_chr{chromosome}.tsv.gz"
        
    if not os.path.exists(load_path): 
        msg = f"❌ CRITICAL ERROR: Internal EAF column missing. The provided External EAF file was also not found: {load_path}"
        log_print(msg, to_screen=True)
        raise sys.exit(1) 
        
    log_print(f"📂 Loading external file: {os.path.basename(load_path)}")
    
    try:
        try: ext_df = pl.read_csv(load_path, separator="\t", has_header=True)
        except: ext_df = pl.read_csv(load_path, separator=" ", has_header=True)
    except:
        import pandas as pd
        ext_df = pl.from_pandas(pd.read_csv(load_path, delim_whitespace=True))
        
    #Identify which required columns are missing
    missing_cols = []
    for k in ["chr", "pos", "a1", "a2"]:
        actual_col_name = colmap[k]
        if actual_col_name not in ext_df.columns:
            missing_cols.append(actual_col_name)
    if missing_cols:
        msg = f"❌ CRITICAL ERROR: External reference file is missing required columns: {missing_cols}"
        log_print(msg, to_screen=True)
        sys.exit(1)

    ext_df = (
        ext_df.select([colmap["chr"], colmap["pos"], colmap["a1"], colmap["a2"], eaf_col])
        .rename({
            colmap["chr"]: std_cols["chr"], colmap["pos"]: std_cols["pos"],
            colmap["a1"]: std_cols["ea"], colmap["a2"]: std_cols["oa"]
        })
        .with_columns([
            pl.col(std_cols["chr"]).cast(pl.Utf8).str.replace("0X", "X"),
            pl.col(std_cols["pos"]).cast(pl.Int64),
            pl.col(std_cols["ea"]).str.to_uppercase(),
            pl.col(std_cols["oa"]).str.to_uppercase(),
            pl.col(eaf_col).cast(pl.Float64)
        ])
    )

    df = df.with_columns([pl.col(std_cols["ea"]).str.to_uppercase(), pl.col(std_cols["oa"]).str.to_uppercase()])
    join_cols = [std_cols["chr"], std_cols["pos"], std_cols["ea"], std_cols["oa"]]
    
    df1 = df.join(ext_df, on=join_cols, how="inner")
    ext_flip = ext_df.rename({std_cols["ea"]: std_cols["oa"], std_cols["oa"]: std_cols["ea"]}).with_columns((1 - pl.col(eaf_col)).alias(eaf_col))
    df2 = df.join(ext_flip, on=join_cols, how="inner")
    
    return pl.concat([df1, df2]).unique(join_cols), eaf_col

# -------------------------------------------------------
# 5. MAIN FUNCTION
# -------------------------------------------------------
def add_or_calculate_eaf(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    eaffile: str = "NA",
    default_eaf_file: str = "NA",
    default_eaf_eafcolumn: str = "EAF",
    maf_eaf_decision_cutoff: float = 0.95,
    external_eaf_colmap: Optional[Dict] = None,
) -> Tuple[pl.DataFrame, Dict, dict]:

    # 1. Setup
    gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
    output_dir = sample_column_dict.get("output_folder", ".")
    _, log_print = setup_logging(output_dir, gwas_outputname, chromosome)
    
    qc_info = {"Total_number_of_variants": df.height}
    log_print("\n🧩 Starting EAF Harmonization...")

    try:
        std_cols = {
            "chr": sample_column_dict["chr_col"], "pos": sample_column_dict["pos_col"],
            "ea": sample_column_dict["ea_col"], "oa": sample_column_dict["oa_col"]
        }
        internal_eaf = sample_column_dict.get("eaf_col", "NA")
        if internal_eaf == "NA": internal_eaf = None
    except KeyError:
        log_print("❌ Missing columns in config", to_screen=True)
        raise SystemExit(1)

    if external_eaf_colmap is None:
        external_eaf_colmap = {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF"}

    final_df = None
    final_eaf_col = None
    
    # ---------------------------------------------------
    # STEP 1: Internal EAF
    # ---------------------------------------------------
    if internal_eaf and internal_eaf in df.columns:
        valid_range, is_suspicious,out_of_range,n_missing,df,low_freq_pct = check_eaf_statistics(df, internal_eaf, log_print, maf_eaf_decision_cutoff)
        stats = df.select([
            pl.col(internal_eaf).min().alias("min"),
            pl.col(internal_eaf).max().alias("max"),
            pl.col(internal_eaf).mean().alias("mean"),
            pl.col(internal_eaf).median().alias("median"),
            pl.len().alias("total")
        ]).to_dicts()[0]
        qc_info.update({
            "Internal_AF": "present",
            "variants_with_InternalAF_before_removing_Missing_AF": qc_info["Total_number_of_variants"],
            "variants_with_InternalAF_with_outofrange_AF_(clipped)": out_of_range,
            "variants_with_InternalAF_with_missing_AF":n_missing,
            "variants_with_InternalAF_After_removing_Missing_AF": stats["total"],
            "variants_with_InternalAF_AfterQC_min_AF": stats["min"],
            "variants_with_InternalAF_AfterQC_max_AF": stats["max"],
            "variants_with_InternalAF_AfterQC_mean_AF": stats["mean"],
            "variants_with_InternalAF_AfterQC_median_AF": stats["median"],
        })
        if not valid_range:
            log_print("🛑 Internal EAF invalid (>1.0). Skipping.", to_screen=True)
        elif not is_suspicious:
            log_print("✔ Internal EAF Valid. Using it.")
            final_df = df; final_eaf_col = internal_eaf
            qc_info.update({"Internal_AF_Identified_AS":"EAF" })
        else:
            log_print(f"\t\t\t ⚠️  QC Warning: Column '{internal_eaf}' resembles MAF rather than EAF ({low_freq_pct*100:.2f}% ≤ 0.5).", to_screen=True)
            log_print(f"\t\t\t ℹ️  Validating against default eaf reference panel: {default_eaf_file}...", to_screen=True)
            # It's suspicious. Check if it is MAF error.
            is_maf_error = check_if_data_is_actually_maf(
                df, default_eaf_file, default_eaf_eafcolumn, internal_eaf, std_cols, chromosome, log_print
            )
            
            if is_maf_error:
                # Returns True -> Confirmed MAF -> EXIT
                log_print(f"\t\t\t 🛑 CRITICAL: Column '{internal_eaf}' verified as MAF (not EAF). Pipeline stopped to prevent errors.", to_screen=True)
                qc_info.update({"Internal_AF_Identified_AS":"MAF" })
                raise SystemExit(1)
            else:
                # Returns False -> Confirmed Rare -> USE
                log_print("\t\t\t ✔ Verification Passed: Internal column confirmed as rare variants.", to_screen=True)
                final_df = df; final_eaf_col = internal_eaf

    # ---------------------------------------------------
    # STEP 2: External EAF (Fallback)
    # ---------------------------------------------------
    if final_df is None:
        log_print("ℹ Using External Reference...")
        target_file = eaffile if (eaffile != "NA" and os.path.exists(eaffile)) else default_eaf_file
        # 2. Detailed Log Message
        if target_file == default_eaf_file:
            log_print(f"\t\t\t-> No INTERNAL EAF column was present in the summary statistics file.")
            log_print(f"\t\t\t-> User did not provide external EAF file (or file not found), so default EAF file is being used: {target_file}")
        else:
            log_print(f"\t\t\t-> No INTERNAL EAF column was present in the summary statistics file.")
            log_print(f"\t\t\t-> User provided EAF file is being used: {target_file}")

        # 2. Configure columns based on selection
        target_col = default_eaf_eafcolumn if target_file == default_eaf_file else sample_column_dict.get("eafcolumn", "EAF")
        target_map = {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF"} if target_file == default_eaf_file else external_eaf_colmap
        try:
            merged_df, new_col_name = load_and_merge_external(
                df, target_file, target_col, target_map, std_cols, log_print, chromosome
            )            
            qc_info.update({
                "Internal_AF": "Absent",
                "User_provided_AF":"present",
                "variants_with_External_AF":merged_df.height,
                })
            valid_range, is_suspicious,out_of_range,n_missing,merged_df,low_freq_pct  = check_eaf_statistics(merged_df, new_col_name, log_print, maf_eaf_decision_cutoff)

            stats = merged_df.select([
                pl.col(new_col_name).min().alias("min"),
                pl.col(new_col_name).max().alias("max"),
                pl.col(new_col_name).mean().alias("mean"),
                pl.col(new_col_name).median().alias("median"),
                pl.len().alias("total")
            ]).to_dicts()[0]
            qc_info.update({
                "variants_with_External_AF_with_outofrange_AF_(clipped)": out_of_range,
                "variants_with_External_AF_with_missing_AF":n_missing,
                "variants_with_External_AF_After_removing_Missing_AF": stats["total"],
                "variants_with_External_AF_AfterQC_min_AF": stats["min"],
                "variants_with_External_AF_AfterQC_max_AF": stats["max"],
                "variants_with_External_AF_AfterQC_mean_AF": stats["mean"],
                "variants_with_External_AF_AfterQC_median_AF": stats["median"],
            })
            
            if not is_suspicious:
                final_df = merged_df; final_eaf_col = new_col_name
            elif target_file == default_eaf_file:
                log_print("ℹ Default Ref is suspicious/low. Trusting it as truth.")
                final_df = merged_df; final_eaf_col = new_col_name
            else:
                # Suspicious User File. Verify.
                is_maf_error = check_if_data_is_actually_maf(
                    merged_df, default_eaf_file, default_eaf_eafcolumn, new_col_name, std_cols, chromosome, log_print
                )
                if is_maf_error:
                    log_print("\t\t\t 🛑 CRITICAL: External File verified as MAF. Pipeline stopped.", to_screen=True)
                    raise SystemExit(1)
                else:
                    log_print("\t\t\t ✔ Verification Passed: External file confirmed as rare.", to_screen=True)
                    final_df = merged_df; final_eaf_col = new_col_name

        except Exception as e:
            log_print(f"🛑 Error: {e}", to_screen=True)
            raise SystemExit(1)

    # Finalize
    final_df = final_df.with_columns([
        pl.when(pl.col(final_eaf_col) <= 0.5).then(pl.col(final_eaf_col)).otherwise(1 - pl.col(final_eaf_col)).alias("zmaf"),
        pl.col(final_eaf_col).clip(0.0, 1.0)
    ])
    sample_column_dict["eaf_col"] = final_eaf_col
    qc_info["final_variants"] = final_df.height
    log_print(f"✅ Done. Variants: {final_df.height}")
    
    return final_df, qc_info, sample_column_dict