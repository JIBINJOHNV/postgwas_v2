import os
import sys
from pathlib import Path
from typing import Tuple, Dict, Optional, Callable, Any
from datetime import datetime
import polars as pl


# -------------------------------------------------------
# 1. Helper: Setup Logging
# -------------------------------------------------------
def setup_logging(output_dir: str, gwas_outputname: str, chromosome: str) -> Tuple[Path, Callable, list]:
    log_dir = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    chr_suffix = f"_chr{chromosome}" if chromosome else ""
    log_file = log_dir / f"{gwas_outputname}{chr_suffix}_eaf.log"
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "w") as f:
        f.write(f"--- Log started: {gwas_outputname} (Chr {chromosome}) at {start_time} ---\n")
    log_buffer = []  # 🔴 NEW
    
    def log_print(*args, to_screen: bool = False):  # 🔴 force default False
        now = datetime.now().strftime("%H:%M:%S")
        msg_content = " ".join(str(a) for a in args)
        log_entry = f"[{now}] {msg_content}"
        # write to file immediately
        try:
            with open(log_file, "a") as f:
                f.write(log_entry + "\n")
        except Exception as e:
            print(f"⚠️ LOGGING FAILED: {e}", file=sys.stderr)
        # store in buffer
        log_buffer.append(msg_content)
        # ONLY print if explicitly forced (we won’t use this anywhere)
        if to_screen:
            print(msg_content)
    return log_file, log_print, log_buffer


# -------------------------------------------------------
# 2. Helper: Compute detailed AF/EAF statistics
# -------------------------------------------------------
def summarize_af_column(
    df: pl.DataFrame,
    colname: str,
) -> Dict[str, Any]:
    """
    Comprehensive descriptive statistics for AF/EAF-like columns.
    Works on the DataFrame exactly as passed in.
    """
    if colname not in df.columns:
        return {
            "column_present": False,
            "total_rows": df.height,
        }
    total_rows = df.height
    if total_rows == 0:
        return {
            "column_present": True,
            "total_rows": 0,
            "non_null_rows": 0,
            "null_rows": 0,
        }
    stats_row = df.select([
        pl.len().alias("total_rows"),
        pl.col(colname).is_not_null().sum().alias("non_null_rows"),
        pl.col(colname).is_null().sum().alias("null_rows"),
        (pl.col(colname) < 0.0).sum().alias("lt_0_rows"),
        (pl.col(colname) > 1.0).sum().alias("gt_1_rows"),
        (pl.col(colname) > 1.05).sum().alias("gt_1_05_rows"),
        (pl.col(colname) < 0.0).sum().alias("invalid_negative_rows"),
        ((pl.col(colname) < 0.0) | (pl.col(colname) > 1.05)).sum().alias("invalid_total_rows"),
        (pl.col(colname) == 0.0).sum().alias("eq_0_rows"),
        (pl.col(colname) == 1.0).sum().alias("eq_1_rows"),
        (pl.col(colname) <= 0.5).sum().alias("le_0_5_rows"),
        (pl.col(colname) < 0.01).sum().alias("lt_0_01_rows"),
        ((pl.col(colname) >= 0.01) & (pl.col(colname) <= 0.99)).sum().alias("between_0_01_and_0_99_rows"),
        ((pl.col(colname) >= 0.4) & (pl.col(colname) <= 0.6)).sum().alias("between_0_4_and_0_6_rows"),
        pl.col(colname).min().alias("min"),
        pl.col(colname).max().alias("max"),
        pl.col(colname).mean().alias("mean"),
        pl.col(colname).median().alias("median"),
        pl.col(colname).quantile(0.25).alias("q25"),
        pl.col(colname).quantile(0.75).alias("q75"),
    ]).to_dicts()[0]
    non_null_rows = stats_row["non_null_rows"] or 0
    total_rows = stats_row["total_rows"] or 0
    def frac(n, d):
        return float(n) / float(d) if d else 0.0
    stats_row["null_fraction"] = frac(stats_row["null_rows"], total_rows)
    stats_row["invalid_total_fraction"] = frac(stats_row["invalid_total_rows"], total_rows)
    stats_row["le_0_5_fraction_among_total"] = frac(stats_row["le_0_5_rows"], total_rows)
    stats_row["le_0_5_fraction_among_non_null"] = frac(stats_row["le_0_5_rows"], non_null_rows)
    stats_row["column_present"] = True
    return stats_row


# -------------------------------------------------------
# 3. Helper: Validate + clean EAF stats
# -------------------------------------------------------
def check_eaf_statistics(
    df: pl.DataFrame,
    colname: str,
    log_print: Callable,
    maf_cutoff: float = 0.95,
    invalid_fraction_cutoff: float = 0.05,
    missing_fraction_cutoff: float = 0.05,
) -> Tuple[bool, bool, int, int, pl.DataFrame, float, Dict[str, Any]]:
    """
    Validate and clean AF/EAF column.
    Returns
    -------
    (
        is_valid_range,
        is_suspicious_maf,
        out_of_range_count,
        missing_count,
        cleaned_df,
        low_freq_pct_after_cleaning,
        stats_dict
    )
    """
    if colname not in df.columns:
        stats = {"column_present": False, "total_rows": df.height}
        return False, False, 0, df.height, df, 0.0, stats
    total_before = df.height
    if total_before == 0:
        stats = summarize_af_column(df, colname)
        return False, False, 0, 0, df, 0.0, stats
    initial_stats = summarize_af_column(df, colname)
    out_of_range = df.filter((pl.col(colname) < 0.0) | (pl.col(colname) > 1.05)).height
    out_of_range_fraction = out_of_range / total_before if total_before else 0.0
    if out_of_range_fraction > invalid_fraction_cutoff:
        error_msg = (
            f"❌ CRITICAL ERROR: During '{colname}' check, observed too many invalid variants "
            f"with Allele Frequency > 1.05 or < 0.0.\n"
            f"   Found {out_of_range} bad rows out of {total_before} ({out_of_range_fraction:.1%})."
        )
        log_print(error_msg, to_screen=True)
    cleaned_df = (
        df.with_columns(
            pl.when(pl.col(colname).is_not_null())
            .then(pl.col(colname).clip(0.0, 1.0))
            .otherwise(None)
            .alias(colname)
        ))
    n_missing = total_before - cleaned_df.height
    missing_fraction = n_missing / total_before if total_before else 0.0
    if total_before > 0 and missing_fraction > missing_fraction_cutoff:
        warn_msg = (
            f"❌ ⚠️ WARNING: During '{colname}' check, observed high missing (null) values.\n"
            f"   Found {n_missing} missing rows out of {total_before} ({missing_fraction:.1%})"
        )
        # log_print handles the file writing, and the terminal will interpret the colors
        log_print(warn_msg, to_screen=True)
    total_after = cleaned_df.height
    low_freq = cleaned_df.filter(pl.col(colname) <= 0.5).height
    low_freq_pct = low_freq / total_after if total_after else 0.0
    is_suspicious = False
    if low_freq_pct > maf_cutoff:
        log_print(
            f"\t\t\t⚠️  {colname}: {low_freq_pct*100:.2f}% values ≤ 0.5 after cleaning (looks like MAF).",
            to_screen=True
        )
        is_suspicious = True
    final_stats = summarize_af_column(cleaned_df, colname)
    combined_stats = {
        "initial": initial_stats,
        "final": final_stats,
        "total_before_cleaning": total_before,
        "total_after_cleaning": total_after,
        "out_of_range_count": out_of_range,
        "out_of_range_fraction": out_of_range_fraction,
        "missing_count_removed": n_missing,
        "missing_fraction_removed": missing_fraction,
        "low_freq_count_after_cleaning": low_freq,
        "low_freq_fraction_after_cleaning": low_freq_pct,
        "maf_like_flag": is_suspicious,
    }
    return True, is_suspicious, out_of_range, n_missing, cleaned_df, low_freq_pct, combined_stats


# -------------------------------------------------------
# 4. Helper: Cross-Check against Default Reference
# -------------------------------------------------------
def check_if_data_is_actually_maf(
    df: pl.DataFrame,
    default_path: str,
    default_af_col: str,
    target_eaf_col: str,
    std_cols: dict,
    chromosome: str,
    log_print: Callable
) -> Tuple[Optional[bool], Dict[str, Any]]:
    """
    Returns
    -------
    (result, stats)

    result:
        True  -> confirmed MAF-like error
        False -> looks like valid EAF / rare-variant-compatible
        None  -> could not determine
    """
    stats: Dict[str, Any] = {
        "reference_path_requested": default_path,
        "reference_path_used": None,
        "reference_found": False,
        "overlap_variants": 0,
        "suspicious_subset_variants": 0,
        "confirmed_common_in_ref": 0,
        "error_rate": None,
        "decision": "undetermined",
    }
    log_print("              🔍 Validating against Default Reference to distinguish MAF vs Rare...", to_screen=True)
    load_path = default_path
    chr_specific_path = f"{default_path}_chr{chromosome}.tsv.gz"
    if not os.path.exists(load_path) and os.path.exists(chr_specific_path):
        load_path = chr_specific_path
    stats["reference_path_used"] = load_path
    if not os.path.exists(load_path):
        log_print(f"              ⚠️ Default Reference not found: {load_path}. Cannot verify suspicious data.", to_screen=True)
        stats["decision"] = "reference_not_found"
        return None, stats
    stats["reference_found"] = True
    try:
        try:
            ref_df = pl.read_csv(load_path, separator="\t", has_header=True)
        except Exception:
            ref_df = pl.read_csv(load_path, separator=" ", has_header=True)
        req_cols = ["CHROM", "POS", "REF", "ALT", default_af_col]
        missing_cols = [c for c in req_cols if c not in ref_df.columns]
        if missing_cols:
            log_print(f"❌ Default file missing columns {missing_cols}", to_screen=True)
            stats["decision"] = "reference_missing_required_columns"
            return None, stats
        ref_df = ref_df.with_columns([
            pl.col("CHROM").cast(pl.Utf8).str.replace("0X", "X"),
            pl.col("POS").cast(pl.Int64),
            pl.col("REF").cast(pl.Utf8).str.to_uppercase(),
            pl.col("ALT").cast(pl.Utf8).str.to_uppercase(),
            pl.col(default_af_col).cast(pl.Float64),
        ])
    except Exception as e:
        log_print(f"❌ Error loading default ref: {e}", to_screen=True)
        stats["decision"] = "reference_load_error"
        stats["reference_load_error"] = str(e)
        return None, stats
    try:
        data_df = df.with_columns([
            pl.col(std_cols["ea"]).cast(pl.Utf8).str.to_uppercase(),
            pl.col(std_cols["oa"]).cast(pl.Utf8).str.to_uppercase(),
            pl.col(std_cols["chr"]).cast(pl.Utf8).str.replace("0X", "X"),
            pl.col(std_cols["pos"]).cast(pl.Int64),
            pl.col(target_eaf_col).cast(pl.Float64),
        ])
    except Exception as e:
        log_print(f"❌ Error preparing internal dataframe for MAF/EAF validation: {e}", to_screen=True)
        stats["decision"] = "input_preparation_error"
        stats["input_preparation_error"] = str(e)
        return None, stats
    m1 = data_df.join(
        ref_df,
        left_on=[std_cols["chr"], std_cols["pos"], std_cols["ea"], std_cols["oa"]],
        right_on=["CHROM", "POS", "ALT", "REF"],
        how="inner",
    ).select([
        pl.col(target_eaf_col).alias("TEST_EAF"),
        pl.col(default_af_col).alias("REF_EAF"),
    ])
    m2 = data_df.join(
        ref_df,
        left_on=[std_cols["chr"], std_cols["pos"], std_cols["ea"], std_cols["oa"]],
        right_on=["CHROM", "POS", "REF", "ALT"],
        how="inner",
    ).select([
        pl.col(target_eaf_col).alias("TEST_EAF"),
        (1 - pl.col(default_af_col)).alias("REF_EAF"),
    ])
    matched_df = pl.concat([m1, m2], how="vertical") if (m1.height > 0 or m2.height > 0) else pl.DataFrame({
        "TEST_EAF": [],
        "REF_EAF": [],
    })
    stats["overlap_variants"] = matched_df.height
    if matched_df.height < 5:
        log_print("              ⚠️ Not enough matches with Default Reference.", to_screen=True)
        stats["decision"] = "insufficient_overlap"
        return None, stats
    suspicious_subset = matched_df.filter(pl.col("TEST_EAF") <= 0.5)
    total_sus = suspicious_subset.height
    confirmed_common_in_ref = suspicious_subset.filter(pl.col("REF_EAF") > 0.60).height
    error_rate = confirmed_common_in_ref / total_sus if total_sus > 0 else 0.0
    stats["suspicious_subset_variants"] = total_sus
    stats["confirmed_common_in_ref"] = confirmed_common_in_ref
    stats["error_rate"] = error_rate 
    log_print(f"              • Overlap variants: {matched_df.height}", to_screen=True)
    log_print(f"              • Variants with Input ≤ 0.5: {total_sus}", to_screen=True)
    log_print(
        f"              • Of those, Default Ref says > 0.6: {confirmed_common_in_ref} ({error_rate*100:.1f}%)",
        to_screen=True
    )
    if error_rate > 0.50:
        log_print("              🛑 RESULT: Default Reference DISAGREES. (Input is MAF, Ref is Common).", to_screen=True)
        stats["decision"] = "confirmed_maf_error"
        return True, stats
    else:
        log_print("              ✔ RESULT: Default Reference AGREES. (Input is correct, variants are Rare).", to_screen=True)
        stats["decision"] = "confirmed_not_maf"
        return False, stats


# -------------------------------------------------------
# 5. Helper: Load Generic External
# ------------------------------------------------------
def load_and_merge_external(
    df: pl.DataFrame,
    path: str,
    eaf_col: str,
    colmap: dict,
    std_cols: dict,
    log_print: Callable,
    chromosome: str
) -> Tuple[pl.DataFrame, str, Dict[str, Any]]:
    stats: Dict[str, Any] = {
        "external_rows_loaded": 0,
        "input_rows_before_merge": df.height,
        "direct_match_rows": 0,
        "flip_match_rows": 0,
        "unmatched_rows": 0,
        "total_missing": 0,
        "missing_fraction": 0.0,
    }
    # --------------------------------------------------
    # Load file
    # -------------------------------------------------- 
    load_path = path
    if not os.path.exists(load_path) and os.path.exists(f"{path}_chr{chromosome}.tsv.gz"):
        load_path = f"{path}_chr{chromosome}.tsv.gz"
    if not os.path.exists(load_path):
        msg = f"❌ CRITICAL ERROR: External EAF file not found: {load_path}"
        log_print(msg, to_screen=True)
        sys.exit(1)
    log_print(f"              📂 Loading external file: {os.path.basename(load_path)}")
    try:
        try:
            ext_df = pl.read_csv(load_path, separator="\t", has_header=True)
        except:
            ext_df = pl.read_csv(load_path, separator=" ", has_header=True)
    except:
        import pandas as pd
        ext_df = pl.from_pandas(pd.read_csv(load_path, delim_whitespace=True))
    stats["external_rows_loaded"] = ext_df.height
    # --------------------------------------------------
    # Validate columns
    # --------------------------------------------------
    missing_cols = []
    for k in ["chr", "pos", "a1", "a2"]:
        if colmap[k] not in ext_df.columns:
            missing_cols.append(colmap[k])
    if missing_cols:
        msg = f"❌ CRITICAL ERROR: Missing columns in external file: {missing_cols}"
        log_print(msg, to_screen=True)
        sys.exit(1)
    # --------------------------------------------------
    # Standardize external
    # --------------------------------------------------
    ext_df = (
        ext_df.select([
            colmap["chr"], colmap["pos"], colmap["a1"], colmap["a2"], eaf_col
        ])
        .rename({
            colmap["chr"]: std_cols["chr"],
            colmap["pos"]: std_cols["pos"],
            colmap["a1"]: std_cols["ea"],
            colmap["a2"]: std_cols["oa"],
        })
        .with_columns([
            pl.col(std_cols["chr"]).cast(pl.Utf8).str.replace("0X", "X"),
            pl.col(std_cols["pos"]).cast(pl.Int64),
            pl.col(std_cols["ea"]).str.to_uppercase(),
            pl.col(std_cols["oa"]).str.to_uppercase(),
            pl.col(eaf_col).cast(pl.Float64),
        ])
    ) # drop duplicates
    subset_cols = [std_cols["chr"], std_cols["pos"], std_cols["ea"], std_cols["oa"]]
    ext_df = ext_df.unique(subset=subset_cols, keep="first")
    # --------------------------------------------------
    # Standardize input df
    # --------------------------------------------------
    df = df.with_columns([
        pl.col(std_cols["chr"]).cast(pl.Utf8).str.replace("0X", "X"),
        pl.col(std_cols["pos"]).cast(pl.Int64),
        pl.col(std_cols["ea"]).str.to_uppercase(),
        pl.col(std_cols["oa"]).str.to_uppercase(),
    ])
    join_cols = [
        std_cols["chr"],
        std_cols["pos"],
        std_cols["ea"],
        std_cols["oa"],
    ]
    # ==================================================
    # STEP 1: DIRECT MATCH
    # ==================================================
    df1 = df.join(ext_df, on=join_cols, how="inner")
    stats["direct_match_rows"] = df1.height
    # Remaining after direct match
    remaining_df = df.join(df1.select(join_cols), on=join_cols, how="anti")
    # ==================================================
    # STEP 2: FLIP MATCH (ONLY ON REMAINING)
    # ==================================================
    ext_flip = (
        ext_df
        .rename({
            std_cols["ea"]: std_cols["oa"],
            std_cols["oa"]: std_cols["ea"]
        })
        .with_columns((1 - pl.col(eaf_col)).alias(eaf_col))
    )
    df2 = remaining_df.join(ext_flip, on=join_cols, how="inner")
    stats["flip_match_rows"] = df2.height
    # Remaining after flip match
    remaining_df2 = remaining_df.join(
        df2.select(join_cols),
        on=join_cols,
        how="anti"
    )
    # ==================================================
    # STEP 3: UNMATCHED → assign NULL EAF
    # ==================================================
    unmatched_df = remaining_df2.with_columns(
        pl.lit(None).cast(pl.Float64).alias(eaf_col)
    )
    stats["unmatched_rows"] = unmatched_df.height
    stats["total_missing"] = unmatched_df.height
    stats["missing_fraction"] = unmatched_df.height / df.height
    # ==================================================
    # FINAL CONCAT (preserve all input variants)
    # ==================================================
    final_df = pl.concat([df1, df2, unmatched_df], how="vertical")
    log_print(f"\t\t\t\tDirect matches   : {stats['direct_match_rows']}")
    log_print(f"\t\t\t\tFlip matches     : {stats['flip_match_rows']}")
    log_print(f"\t\t\t\tUnmatched        : {stats['unmatched_rows']}")
    log_print(f"\t\t\t\tMissing EAF frac : {stats['missing_fraction']:.4f}")
    return final_df, eaf_col, stats


# -------------------------------------------------------
# 6. Helper: Write structured statistics into qc_info
# -------------------------------------------------------
def add_prefixed_stats(
    qc_info: Dict[str, Any],
    prefix: str,
    stats: Dict[str, Any]
) -> None:
    for key, value in stats.items():
        if isinstance(value, dict):
            for sub_key, sub_val in value.items():
                qc_info[f"{prefix}_{key}_{sub_key}"] = sub_val
        else:
            qc_info[f"{prefix}_{key}"] = value


# -------------------------------------------------------
# 7. MAIN FUNCTION
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

    gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
    output_dir = sample_column_dict.get("output_folder", ".")
    merge_stats = {}
    PREFIX = " " * 16

    # 🔴 Create subdirectory for EAF QC outputs
    eaf_qc_dir = os.path.join(output_dir, "eaf_qc")
    os.makedirs(eaf_qc_dir, exist_ok=True)

    _, base_log_print, log_buffer = setup_logging(output_dir, gwas_outputname, chromosome)

    def log_print(*args, to_screen: bool = False):
        base_log_print(*args, to_screen=False)

    qc_info: Dict[str, Any] = {
        "chromosome": chromosome,
        "gwas_outputname": gwas_outputname,
        "Total_number_of_variants": df.height,
        "final_status": "started",
        "decision_source": None,
        "final_eaf_col": None,
    }

    log_print(f"\n📊 [Chr {chromosome}] Starting EAF Harmonization...")

    try:
        std_cols = {
            "chr": sample_column_dict["chr_col"],
            "pos": sample_column_dict["pos_col"],
            "ea": sample_column_dict["ea_col"],
            "oa": sample_column_dict["oa_col"],
        }

        internal_eaf = sample_column_dict.get("eaf_col", "NA")
        if internal_eaf == "NA":
            internal_eaf = None

        if external_eaf_colmap is None:
            external_eaf_colmap = {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF"}

        final_df = None
        final_eaf_col = None

        # ------------------------------
        # STEP 1: Internal EAF
        # ------------------------------
        if internal_eaf and internal_eaf in df.columns:
            log_print("=" * 80)
            log_print(f"{PREFIX}🔎 Internal EAF column detected: {internal_eaf}")
            
            valid_range, is_suspicious, out_of_range, n_missing, cleaned_df, low_freq_pct, stats = \
                check_eaf_statistics(df, internal_eaf, log_print, maf_eaf_decision_cutoff)

            qc_info["Internal_AF"] = "present"
            add_prefixed_stats(qc_info, "internal_af", stats)

            if valid_range and not is_suspicious:
                final_df = cleaned_df
                final_eaf_col = internal_eaf
                qc_info["decision_source"] = "internal_eaf_direct"

            elif valid_range and is_suspicious:
                log_print(f"              ⚠️ Column '{internal_eaf}' looks like MAF ({low_freq_pct*100:.2f}%)")

                is_maf_error, maf_stats = check_if_data_is_actually_maf(
                    cleaned_df,
                    default_eaf_file,
                    default_eaf_eafcolumn,
                    internal_eaf,
                    std_cols,
                    chromosome,
                    log_print,
                )

                add_prefixed_stats(qc_info, "internal_af_reference_check", maf_stats)

                if is_maf_error is True:
                    qc_info["final_status"] = "failed"
                    raise RuntimeError(f"{internal_eaf} is MAF not EAF")

                final_df = cleaned_df
                final_eaf_col = internal_eaf
                qc_info["decision_source"] = "internal_eaf_checked"

        else:
            qc_info["Internal_AF"] = "Absent"
            log_print(f"{PREFIX}ℹ No internal EAF column present")

        # ------------------------------
        # STEP 2: External fallback
        # ------------------------------
        if final_df is None:
            log_print(f"{PREFIX}ℹ Using External Reference")

            user_external_exists = (eaffile != "NA" and os.path.exists(eaffile))
            target_file = eaffile if user_external_exists else default_eaf_file

            target_col = (
                default_eaf_eafcolumn
                if target_file == default_eaf_file
                else sample_column_dict.get("eafcolumn", "EAF")
            )

            target_map = (
                {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF"}
                if target_file == default_eaf_file
                else external_eaf_colmap
            )

            merged_df, new_col_name, merge_stats = load_and_merge_external(
                df, target_file, target_col, target_map, std_cols, log_print, chromosome
            )

            add_prefixed_stats(qc_info, "external_merge", merge_stats)

            _, _, _, _, cleaned_df, _, stats = check_eaf_statistics(
                merged_df, new_col_name, log_print, maf_eaf_decision_cutoff
            )

            add_prefixed_stats(qc_info, "external_af", stats)

            final_df = cleaned_df
            final_eaf_col = new_col_name
            qc_info["decision_source"] = "external_eaf"

        # ------------------------------
        # FINALIZE
        # ------------------------------
        if final_df is None or final_eaf_col is None:
            raise RuntimeError("❌ Failed to determine EAF column")

        # 🔴 QC BEFORE CLIP
        missing_eaf_df = final_df.filter(pl.col(final_eaf_col).is_null())
        out_of_bound_df = final_df.filter(
            (pl.col(final_eaf_col) < 0) | (pl.col(final_eaf_col) > 1)
        )

        missing_eaf_count = missing_eaf_df.height
        out_of_bound_count = out_of_bound_df.height

        # Save QC files inside subdirectory
        missing_file = os.path.join(
            eaf_qc_dir, f"{gwas_outputname}_missing_eaf_chr{chromosome}.tsv"
        )
        outofbound_file = os.path.join(
            eaf_qc_dir, f"{gwas_outputname}_outofbound_eaf_chr{chromosome}.tsv"
        )

        if missing_eaf_count > 0:
            missing_eaf_df.write_csv(missing_file, separator="\t")
            log_print(f"{PREFIX}🧾 Missing EAF variants saved → {missing_file}")

        if out_of_bound_count > 0:
            out_of_bound_df.write_csv(outofbound_file, separator="\t")
            log_print(f"{PREFIX}🧾 Out-of-bound EAF variants saved → {outofbound_file}")

        qc_info["missing_eaf_count"] = missing_eaf_count
        qc_info["out_of_bound_eaf_count"] = out_of_bound_count

        # Original logic unchanged
        final_df = final_df.with_columns([
            pl.when(pl.col(final_eaf_col).is_null())
            .then(None)
            .when(pl.col(final_eaf_col) <= 0.5)
            .then(pl.col(final_eaf_col))
            .otherwise(1 - pl.col(final_eaf_col))
            .alias("zmaf"),

            pl.col(final_eaf_col).clip(0.0, 1.0).alias(final_eaf_col),
        ])

        final_stats = summarize_af_column(final_df, final_eaf_col)
        add_prefixed_stats(qc_info, "final_eaf", final_stats)

        sample_column_dict["eaf_col"] = final_eaf_col
        qc_info["final_eaf_col"] = final_eaf_col
        qc_info["final_variants"] = final_df.height
        qc_info["final_status"] = "success"
        log_print(f"{PREFIX}Initial variants              :{qc_info['Total_number_of_variants']}")
        log_print(f"{PREFIX}Missing EAF variants          :{missing_eaf_count}")
        log_print(f"{PREFIX}Out-of-bound EAF variants     :{out_of_bound_count}")
        log_print(f"{PREFIX}Final variants                :{qc_info['final_variants']}")
        log_print(f"{PREFIX}Variants removed              :{qc_info['Total_number_of_variants'] - qc_info['final_variants']}")
        log_print(f"{PREFIX}✅ EAF Harmonization completed successfully for chromosome chr{chromosome}")
        log_print(" ")

    except Exception as e:
        qc_info["final_status"] = "failed"
        qc_info["error_message"] = str(e)
        raise

    finally:
        print("\n".join(log_buffer))

    return final_df, qc_info, sample_column_dict
# def add_or_calculate_eaf(
#     chromosome: str,
#     df: pl.DataFrame,
#     sample_column_dict: dict,
#     eaffile: str = "NA",
#     default_eaf_file: str = "NA",
#     default_eaf_eafcolumn: str = "EAF",
#     maf_eaf_decision_cutoff: float = 0.95,
#     external_eaf_colmap: Optional[Dict] = None,
# ) -> Tuple[pl.DataFrame, Dict, dict]:
#     gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
#     output_dir = sample_column_dict.get("output_folder", ".")
#     merge_stats = {}
#     # --------------------------------------------------
#     # 🔴 Setup logging (buffered)
#     # --------------------------------------------------
#     _, base_log_print, log_buffer = setup_logging(output_dir, gwas_outputname, chromosome)
#     # 🔴 HARD override → NEVER print during execution
#     def log_print(*args, to_screen: bool = False):
#         base_log_print(*args, to_screen=False)
    
#     qc_info: Dict[str, Any] = {
#         "chromosome": chromosome,
#         "gwas_outputname": gwas_outputname,
#         "Total_number_of_variants": df.height,
#         "final_status": "started",
#         "decision_source": None,
#         "final_eaf_col": None,
#     }
#     log_print(f"\n📊 [Chr {chromosome}] Starting EAF Harmonization...")
#     try:
#         # --------------------------------------------------
#         # Config validation
#         # --------------------------------------------------
#         std_cols = {
#             "chr": sample_column_dict["chr_col"],
#             "pos": sample_column_dict["pos_col"],
#             "ea": sample_column_dict["ea_col"],
#             "oa": sample_column_dict["oa_col"],
#         }
#         internal_eaf = sample_column_dict.get("eaf_col", "NA")
#         if internal_eaf == "NA":
#             internal_eaf = None
#         if external_eaf_colmap is None:
#             external_eaf_colmap = {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF"}
#         final_df: Optional[pl.DataFrame] = None
#         final_eaf_col: Optional[str] = None
#         # --------------------------------------------------
#         # STEP 1: Internal EAF
#         # --------------------------------------------------
#         if internal_eaf and internal_eaf in df.columns: 
#             log_print(f"              🔎 Internal EAF column detected: {internal_eaf}")
#             valid_range, is_suspicious, out_of_range, n_missing, cleaned_df, low_freq_pct, stats = \
#                 check_eaf_statistics(df, internal_eaf, log_print, maf_eaf_decision_cutoff)
#             qc_info["Internal_AF"] = "present"
#             add_prefixed_stats(qc_info, "internal_af", stats)
#             if valid_range and not is_suspicious:
#                 final_df = cleaned_df
#                 final_eaf_col = internal_eaf
#                 qc_info["decision_source"] = "internal_eaf_direct"
#             elif valid_range and is_suspicious:
#                 log_print(f"              ⚠️ Column '{internal_eaf}' looks like MAF ({low_freq_pct*100:.2f}%)")
#                 is_maf_error, maf_stats = check_if_data_is_actually_maf(
#                     cleaned_df,
#                     default_eaf_file,
#                     default_eaf_eafcolumn,
#                     internal_eaf,
#                     std_cols,
#                     chromosome,
#                     log_print,
#                 )
#                 add_prefixed_stats(qc_info, "internal_af_reference_check", maf_stats)
#                 if is_maf_error is True:
#                     qc_info["final_status"] = "failed"
#                     raise RuntimeError(f"{internal_eaf} is MAF not EAF")
#                 final_df = cleaned_df
#                 final_eaf_col = internal_eaf
#                 qc_info["decision_source"] = "internal_eaf_checked"
#         else:
#             qc_info["Internal_AF"] = "Absent"
#             log_print("              ℹ No internal EAF column present")
#         # --------------------------------------------------
#         # STEP 2: External fallback
#         # --------------------------------------------------
#         if final_df is None:
#             log_print("              ℹ Using External Reference")       
#             user_external_exists = (eaffile != "NA" and os.path.exists(eaffile))
#             target_file = eaffile if user_external_exists else default_eaf_file
#             target_col = (
#                 default_eaf_eafcolumn
#                 if target_file == default_eaf_file
#                 else sample_column_dict.get("eafcolumn", "EAF")
#             )
#             target_map = (
#                 {"chr": "CHROM", "pos": "POS", "a1": "ALT", "a2": "REF"}
#                 if target_file == default_eaf_file
#                 else external_eaf_colmap
#             )
#             merged_df, new_col_name, merge_stats = load_and_merge_external(
#                 df, target_file, target_col, target_map, std_cols, log_print, chromosome
#             )
#             #qc_info["merged_rows"] = merge_stats["matched_rows"] 
#             add_prefixed_stats(qc_info, "external_merge", merge_stats)
#             _, _, _, _, cleaned_df, _, stats = check_eaf_statistics(
#                 merged_df, new_col_name, log_print, maf_eaf_decision_cutoff
#             )
#             add_prefixed_stats(qc_info, "external_af", stats)
#             final_df = cleaned_df 
#             final_eaf_col = new_col_name
#             qc_info["decision_source"] = "external_eaf"
#         # --------------------------------------------------
#         # FINALIZE
#         # --------------------------------------------------
#         if final_df is None or final_eaf_col is None:
#             raise RuntimeError("❌ Failed to determine EAF column")
#         final_df = final_df.with_columns([
#             pl.when(pl.col(final_eaf_col) <= 0.5)
#             .then(pl.col(final_eaf_col))
#             .otherwise(1 - pl.col(final_eaf_col))
#             .alias("zmaf"),
#             pl.col(final_eaf_col).clip(0.0, 1.0).alias(final_eaf_col),
#         ])
#         final_stats = summarize_af_column(final_df, final_eaf_col)
#         add_prefixed_stats(qc_info, "final_eaf", final_stats)
#         sample_column_dict["eaf_col"] = final_eaf_col
#         qc_info["final_eaf_col"] = final_eaf_col
#         qc_info["final_variants"] = final_df.height 
#         qc_info["final_status"] = "success"
#         log_print(f"              Initial variants:{qc_info['Total_number_of_variants']}") 
#         log_print(f"              Final variants:{qc_info['final_variants']}")
#         log_print(f"              Variants with no EAF: {merge_stats.get('total_missing', 0)}")
#         log_print(f"              Variants with no EAF fraction: {merge_stats.get('missing_fraction', 0)}")
#         log_print(f"              Variants removed:{qc_info['Total_number_of_variants'] - qc_info['final_variants']}")
#         #log_print(f"              Merged variants:{qc_info['merged_rows']}")
#         log_print(f"              ✅ EAF Harmonization completed successfully for chromosome chr{chromosome}")
#         log_print(" ")
#     except Exception as e:
#         qc_info["final_status"] = "failed"
#         qc_info["error_message"] = str(e)
#         raise
#     finally:
#         # 🔥 ALWAYS PRINT LOGS AT END
#         print("\n".join(log_buffer))
#     return final_df, qc_info, sample_column_dict