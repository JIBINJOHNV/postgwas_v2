import polars as pl
import pandas as pd
import pyarrow,sys
import gzip, lzma, zipfile, csv, subprocess, os,json,shlex,argparse
from typing import Tuple
from io import BytesIO,TextIOWrapper
from pathlib import Path
import yaml
from importlib import resources
from postgwas.harmonisation.chr_pos_process import fix_chr_pos_allele_column
from postgwas.utils.detect_delimter import get_polars_separator
import bz2
import io



# ----------------------------------------------------------------------
# 4. Count data lines – **safe** shell command (list form)
# ----------------------------------------------------------------------
def count_data_lines(path: str, skip_hash: bool = True) -> int:
    import subprocess
    import shlex
    import zipfile
    from pathlib import Path
    # ------------------------------
    # ZIP handling (always skip ##)
    # ------------------------------
    if Path(path).suffix.lower() == ".zip":
        try:
            with zipfile.ZipFile(path) as zf:
                files = [f for f in zf.namelist()
                         if f.lower().endswith(('.csv', '.tsv', '.txt', '.gz'))]
                if not files:
                    return 0
                with zf.open(files[0]) as f:
                    count = 0
                    for line in f:
                        if line.startswith(b"##"):   # ✅ ALWAYS skip
                            continue
                        count += 1
                return count
        except Exception as e:
            print(f"Zip count failed: {e}")
            return 0
    # ------------------------------
    # Compression handling
    # ------------------------------
    if path.endswith(".gz"):
        decomp = "gzip -dc"
    elif path.endswith((".xz", ".lzma")):
        decomp = "xz -dc"
    else:
        decomp = "cat"
    # ------------------------------
    # ALWAYS remove ##
    # ------------------------------
    cmd_str = f"{decomp} {shlex.quote(path)} | grep -v '^##' | wc -l"
    try:
        result = subprocess.run(
            cmd_str,
            shell=True,
            capture_output=True,
            text=True,
            check=True
        )
        return int(result.stdout.strip())
    except Exception as e:
        print(f"Shell count failed ({e}) — using Polars count only")
        return 0




# =========================================================
# 1. REMOVE ## LINES (SAFE + CONSISTENT)
# =========================================================
def remove_double_hash(input_file):
    p = Path(input_file)
    def is_gz(path):
        return str(path).endswith((".gz", ".bgz"))
    opener = gzip.open if is_gz(p) else open
    has_double_hash = False
    with opener(p, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("##"):
                has_double_hash = True
                break
    if not has_double_hash:
        print("ℹ️ No '##' lines found → nothing done")
        return input_file  # ✅ FIX: return original instead of None
    # output path
    if is_gz(p):
        out_path = str(p).replace(".gz", "").replace(".bgz", "") + ".tsv"
    else:
        out_path = str(p) + ".tsv"
    with opener(p, "rt", encoding="utf-8", errors="replace") as fin, \
         open(out_path, "w") as fout:
        for line in fin:
            if not line.startswith("##"):
                fout.write(line)
    print(f"✅ '##' lines removed → {out_path}")
    return out_path


# =========================================================
# 2. DELIMITER DETECTION (ROBUST)
# =========================================================
def detect_delimiter(file_path):
    def is_gz(path):
        return path.endswith((".gz", ".bgz"))
    opener = gzip.open if is_gz(file_path) else open
    lines = []
    with opener(file_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("##") or not line.strip():
                continue
            lines.append(line)
            if len(lines) >= 500:  # enough sample
                break
    if not lines:
        return "\t"  # fallback
    sample = "".join(lines)
    candidates = ["\t", ",", ";", " "]
    best_sep = "\t"
    best_score = -1
    for sep in candidates:
        try:
            df = pl.read_csv(
                io.StringIO(sample),
                separator=sep,
                infer_schema_length=50
            )
            ncols = len(df.columns)
            # realistic GWAS structure
            if 5 <= ncols <= 100:
                if ncols > best_score:
                    best_score = ncols
                    best_sep = sep
        except Exception:
            continue
    return best_sep


# =========================================================
# 3. SMART LOADER (CORE FIXED FUNCTION)
# =========================================================
def smart_load_dataframe(file_path: str, sample_column_dict: dict):
    # ----------------------------
    # 1. Remove ##
    # ----------------------------
    new_file = remove_double_hash(file_path)
    # ----------------------------
    # 2. Detect delimiter
    # ----------------------------
    sep = detect_delimiter(new_file)
    print(f"🔍 Detected delimiter: '{repr(sep)}'")
    # ----------------------------
    # 3. Read file (FIXED)
    # ----------------------------
    try:
        df = pl.read_csv(
            new_file,
            separator=sep,
            infer_schema_length=0,   # ✅ FIX (0 is dangerous)
            ignore_errors=False,
            truncate_ragged_lines=False, # ✅ FIX
            has_header=True,
            null_values=["NA", "NAN", ""]
        )
    except Exception as e:
        print(f"⚠️ Primary read failed: {e}")
        print("⚠️ Retrying with low_memory=True...")
        df = pl.read_csv(
            new_file,
            separator=sep,
            infer_schema_length=0,
            low_memory=True,
            truncate_ragged_lines=False,
            has_header=True
        )
    # ----------------------------
    # 4. Required columns
    # ----------------------------
    keys = [
        "chr_col","pos_col","snp_id_col","ea_col","oa_col","eaf_col",
        "beta_or_col","se_col","imp_z_col","pval_col",
        "ncontrol_col","ncase_col","imp_info_col"
    ]
    required_cols = []
    for k in keys:
        v = sample_column_dict.get(k)
        if v and str(v).strip().upper() not in {"", "NA", "NAN", "NONE"}:
            required_cols.append(v)
    # ----------------------------
    # 5. Missing check
    # ----------------------------
    missing_cols = [c for c in required_cols if c not in df.columns]
    total_cols = len(df.columns)
    # Extended structural checks
    unnamed_columns = [c for c in df.columns if c.lower().startswith("unnamed")]
    empty_columns = [c for c in df.columns if str(c).strip() == ""]
    duplicated_columns = [c for c in df.columns if str(c).startswith("_duplicated")]
    # ----------------------------
    # 6. HARD FIX (SHELL)
    # ----------------------------
    # Trigger ONLY on strong delimiter-failure signals
    if (
        missing_cols
        or total_cols == 1
        or len(unnamed_columns) > 5
        or len(empty_columns) > 0
        or len(duplicated_columns) > 5
    ):
        print("⚠️ Detected file structure issues:")
        print(f"   Missing columns: {missing_cols}")
        print(f"   Total columns : {total_cols}")
        print(f"   Unnamed cols  : {unnamed_columns}")
        print(f"   Empty cols    : {empty_columns}")
        print(f"   Duplicated cols: {duplicated_columns}")
        print(f"   Required columns in the GWAS file: {required_cols}")
        print(f"   Columns in the GWAS file: {df.columns}")
        print("   This is likely due to inconsistent whitespace or delimiters.")
        print("🔧 Applying whitespace normalization...")
        fixed_file = str(Path(file_path).with_suffix("")) + ".fixed.tsv"
        if file_path.endswith((".gz", ".bgz")):
            cmd = f'zcat {shlex.quote(file_path)}| tr "\\t" " " | tr -s " " | tr " " "\\t" > {shlex.quote(fixed_file)}'
        else:
            cmd = f'cat {shlex.quote(file_path)} | tr "\\t" " " | tr -s " " | tr " " "\\t" > {shlex.quote(fixed_file)}'
        subprocess.run(cmd, shell=True, check=True)
        df = pl.read_csv(
            fixed_file,
            separator="\t",
            infer_schema_length=0,
            truncate_ragged_lines=False,
            has_header=True,
            null_values=["NA", "NAN", ""]
        )
        missing_cols = [c for c in required_cols if c not in df.columns]
        if missing_cols:
            raise RuntimeError(f"❌ Required columns still missing: {missing_cols}")
    # ----------------------------
    # 7. FINAL SANITY
    # ----------------------------
    if df.height == 0:
        raise RuntimeError("❌ Empty dataframe")
    print(f"\t\t\t📄 Loaded {df.height:,} rows × {len(df.columns)} columns")
    return df, sample_column_dict



def clean_and_optimize_dataframe(df: pl.DataFrame) -> pl.DataFrame:
    """
    1. Trim whitespace from all string columns
    2. Replace empty strings with null
    3. Preserve numeric precision (no unsafe downcasting)
    """
    df = (
        df.lazy()
        # Step 1: trim whitespace
        .with_columns(
            pl.col(pl.String).str.strip_chars()
        )
        # Step 2: empty string → null
        .with_columns(
            pl.col(pl.String).replace("", None)
        )
        .collect()
    )
    return df


# ======================================================================
# HELPER 1: DataFrame Cleanup (Headers, Types, Chromosome)
# ======================================================================
def clean_dataframe(df: pl.DataFrame, chr_col: str) -> pl.DataFrame:
    """
    Applied immediately after ANY read (Attempt 1 or Attempt 2).
    Trims whitespace, fixes types, and standardizes chromosome column.
    """
    if df.height == 0: return df
    # 1. Clean Headers
    df = df.rename({c: c.strip() for c in df.columns})
    # 2. Clean String Values
    for col in df.columns:
        if df[col].dtype in (pl.Utf8, pl.String):
            df = df.with_columns(pl.col(col).str.strip_chars())
    # 3. Type Conversion
    valid_cols = []
    for col in df.columns:
        if col == chr_col:
            valid_cols.append(pl.col(col)) # Skip casting for CHR
            continue
        numeric_series = df[col].cast(pl.Float64, strict=False)
        if numeric_series.null_count() > df[col].null_count():
            valid_cols.append(pl.col(col)) # Keep original if casting fails
        else:
            valid_cols.append(numeric_series)
    df = df.with_columns(valid_cols)
    # 4. Chromosome Fix (12.0 -> 12)
    if chr_col and chr_col in df.columns:
        df = df.with_columns(
            pl.coalesce(
                pl.col(chr_col).cast(pl.Float64, strict=False).cast(pl.Int64, strict=False).cast(pl.String),
                pl.col(chr_col)
            ).alias(chr_col)
        )
    return df



# =========================================================
# 4. READ SUMSTATS (CRITICAL FIXES)
# =========================================================

import os
import polars as pl


def process_impinfo_column(df: pl.DataFrame, sample_column_dict: dict):
    """
    Process IMPINFO column:
    - If imp_info_col is missing/NA -> do nothing
    - If IMPINFO is single-valued -> do nothing
    - If multi-valued and external infofile+infocolumn are provided -> do nothing
    - Otherwise compute row-wise median and update sample_column_dict['imp_info_col']
    """

    def clean_col(x):
        if x is None:
            return None
        x = str(x).strip()
        if x.upper() in {"", "NA", "NAN", "NONE"}:
            return None
        return x

    imp_col = clean_col(sample_column_dict.get("imp_info_col"))
    infofile = clean_col(sample_column_dict.get("infofile"))
    infocolumn = clean_col(sample_column_dict.get("infocolumn"))

    if not imp_col or imp_col not in df.columns:
        if infofile and infocolumn:
            print("ℹ️ IMPINFO column not present or missing in the file and external infofile not provided")
        return df, sample_column_dict

    # detect whether any row has comma-separated multi-values
    has_multi = df.select(
        pl.col(imp_col)
        .cast(pl.Utf8, strict=False)
        .str.contains(",")
        .any()
    ).item()

    if not has_multi:
        return df, sample_column_dict

    # if external info resource is provided, prefer that
    if infofile and infocolumn:
        print(f"ℹ️ Multi-value IMPINFO detected in {imp_col} column; → But external infofile provided that will be used")
        print(f"   Infofile: {infofile}")
        print(f"   Infocolumn: {infocolumn}")
        return df, sample_column_dict

    print(f"⚙️ Multi-value IMPINFO detected in '{imp_col}' column and no external infofile provided")
    print(f"➡️ Computing median value per row for '{imp_col}' column...")
    new_col = f"{imp_col}_median"

    df = df.with_columns(
        pl.col(imp_col)
        .cast(pl.Utf8, strict=False)
        .str.split(",")
        .list.eval(
            pl.when(
                pl.element().cast(pl.Utf8, strict=False).str.to_uppercase() == "NA"
            )
            .then(None)
            .otherwise(pl.element().cast(pl.Float64, strict=False))
        )
        .list.drop_nulls()
        .list.median()
        .alias(new_col)
    )

    sample_column_dict["imp_info_col"] = new_col

    missing_count = df.select(pl.col(new_col).is_null().sum()).item()
    print(f"✅ Created column '{new_col}' and updated sample_column_dict")
    print(f"⚠️ Median IMPINFO missing in {missing_count:,} rows")
    return df, sample_column_dict



def read_sumstats(sumstat_file: str, output_dir: str, sample_column_dict: dict):
    import os
    import polars as pl

    os.makedirs(output_dir, exist_ok=True)

    shell_cnt = count_data_lines(sumstat_file, skip_hash=False) - 1
    df, sample_column_dict = smart_load_dataframe(sumstat_file, sample_column_dict)
    python_read_count = df.height

    removed_coords = 0
    non_standard_allele_count = 0
    removed_duplicates_count = 0
    removed_missing_count = 0

    print(f"Sumstat reading successful: rows loaded: {df.height:,} Columns: {df.columns}\n")

    # =========================================================
    # 🔴 EARLY MISSINGNESS FILTER (NEW POSITION)
    # =========================================================
    keys = [
        "chr_col", "pos_col", "snp_id_col", "ea_col", "oa_col", "eaf_col",
        "beta_or_col", "se_col", "imp_z_col", "pval_col",
        "ncontrol_col", "ncase_col", "imp_info_col"
    ]

    required_cols = []
    for k in keys:
        v = sample_column_dict.get(k)
        if v and str(v).strip().upper() not in {"", "NA", "NAN", "NONE"}:
            required_cols.append(v)

    required_cols = [c for c in required_cols if c in df.columns]

    if required_cols:
        total_rows = df.height

        missing_stats = df.select([
            pl.col(c).is_null().sum().alias(c) for c in required_cols
        ])

        print("📊 Missing values per required column (EARLY QC):")
        for c in required_cols:
            missing_count = missing_stats[c][0]
            missing_pct = (missing_count / total_rows) * 100
            print(f"   {c}: {missing_count:,} ({missing_pct:.2f}%)")

        before_rows = total_rows
        threshold = len(required_cols) * 0.5

        remove_mask = (
            pl.sum_horizontal(
                [pl.col(c).is_not_null().cast(pl.Int64) for c in required_cols]
            ) <= threshold
        )

        removed_df = df.filter(remove_mask)

        gwas_name = sample_column_dict.get("gwas_outputname", "gwas")
        missing_file = os.path.join(output_dir, f"{gwas_name}_missing_filtered_early.tsv")

        if removed_df.height > 0:
            removed_df.write_csv(missing_file, separator="\t")
            print(f"🧾 Early missing-filtered variants saved → {missing_file}")

        df = df.filter(~remove_mask)

        after_rows = df.height
        removed_missing_count = before_rows - after_rows

        print("🧹 Early missingness filtering:")
        print(f"   Before: {before_rows:,}")
        print(f"   After : {after_rows:,}")
        print(f"   Removed: {removed_missing_count:,}")

    # =========================================================
    # COLUMN FIX
    # =========================================================
    df, sample_column_dict = fix_chr_pos_allele_column(
        chromosome="All_Chrs",
        df=df,
        sample_column_dict=sample_column_dict,
        drop_mt=True
    )

    chr_col = sample_column_dict.get("chr_col")
    pos_col = sample_column_dict.get("pos_col")
    ea_col = sample_column_dict.get("ea_col")
    oa_col = sample_column_dict.get("oa_col")

    df = clean_and_optimize_dataframe(df)
    df = clean_dataframe(df, chr_col)

    print(f"\t\t\t✅ Recovery Successful: {df.height:,} rows loaded.")

    # =========================================================
    # CHR FIX
    # =========================================================
    valid_chr = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]

    if chr_col in df.columns:
        df = df.with_columns(
            pl.col(chr_col)
            .cast(pl.Utf8, strict=False)
            .str.strip_chars()
            .str.replace(r"(?i)^chr", "", literal=False)
            .str.to_uppercase()
            .alias(chr_col)
        )

        before_chr = df.height

        chr_filter = pl.col(chr_col).is_in(valid_chr)
        if pos_col in df.columns:
            chr_filter = chr_filter & pl.col(pos_col).is_not_null()

        df = df.filter(chr_filter)

        removed_coords = before_chr - df.height
        print(f"🧩 Removed {removed_coords:,} invalid coordinate rows")

    # =========================================================
    # POS FIX
    # =========================================================
    if pos_col in df.columns:
        df = df.with_columns(
            pl.col(pos_col).cast(pl.Int64, strict=False).alias(pos_col)
        )

    # =========================================================
    # ALLELE FILTER
    # =========================================================
    if ea_col in df.columns and oa_col in df.columns:
        valid_mask = (
            pl.col(ea_col).cast(pl.Utf8, strict=False).str.to_uppercase().str.contains(r"^[ACGT]+$")
            &
            pl.col(oa_col).cast(pl.Utf8, strict=False).str.to_uppercase().str.contains(r"^[ACGT]+$")
        )

        invalid_df = df.filter(~valid_mask)
        df = df.filter(valid_mask)

        non_standard_allele_count = invalid_df.height

        gwas_name = sample_column_dict.get("gwas_outputname", "gwas")
        invalid_file = os.path.join(output_dir, f"{gwas_name}_invalid_alleles.tsv")

        if invalid_df.height > 0:
            invalid_df.write_csv(invalid_file, separator="\t")
            print(f"🧾 Invalid alleles saved → {invalid_file}")

    print(f"🧩 Removed {non_standard_allele_count:,} invalid allele rows")

    # =========================================================
    # DUPLICATES
    # =========================================================
    subset_cols = [chr_col, pos_col, ea_col, oa_col]
    subset_cols = [c for c in subset_cols if c in df.columns]

    if len(subset_cols) == 4:
        before = df.height

        dup_df = df.filter(pl.struct(subset_cols).is_duplicated())
        df = df.sort(subset_cols).unique(subset=subset_cols, keep="first")

        removed_duplicates_count = before - df.height

        gwas_name = sample_column_dict.get("gwas_outputname", "gwas")
        dup_file = os.path.join(output_dir, f"{gwas_name}_duplicates.tsv")

        if dup_df.height > 0:
            dup_df.write_csv(dup_file, separator="\t")
            print(f"🧾 Duplicates saved → {dup_file}")

    print(f"🧩 Removed {removed_duplicates_count:,} duplicate variants")

    # =========================================================
    # IMPINFO PROCESSING
    # =========================================================
    df, sample_column_dict = process_impinfo_column(df, sample_column_dict)

    if abs(shell_cnt - python_read_count) > 1:
        print("⚠️ Shell vs Polars row mismatch")

    return (
        df,
        shell_cnt,
        python_read_count,
        sample_column_dict,
        removed_coords,
        non_standard_allele_count,
        removed_duplicates_count,
        removed_missing_count,
    )
    
# def read_sumstats(sumstat_file: str, output_dir: str, sample_column_dict: dict):
#     os.makedirs(output_dir, exist_ok=True)
#     shell_cnt = count_data_lines(sumstat_file, skip_hash=False) - 1
#     df, sample_column_dict = smart_load_dataframe(sumstat_file, sample_column_dict)
#     python_read_count = df.height
#     removed_coords = 0
#     non_standard_allele_count = 0
#     removed_duplicates_count = 0
#     print(f"Sumstat reading successful: rows loaded: {df.height} Columns: {df.columns}") 
#     print(" \n")
#     # ------------------------------
#     # COLUMN FIX
#     # ------------------------------
#     df, sample_column_dict = fix_chr_pos_allele_column(
#         chromosome="All_Chrs",
#         df=df,
#         sample_column_dict=sample_column_dict,
#         drop_mt=True
#     )
#     print(f"{df.head(5)}")
#     chr_col = sample_column_dict.get("chr_col")
#     pos_col = sample_column_dict.get("pos_col")
#     ea_col = sample_column_dict.get("ea_col")
#     oa_col = sample_column_dict.get("oa_col")
#     df = clean_and_optimize_dataframe(df)
#     df = clean_dataframe(df, chr_col)
#     print(f"\t\t\t✅ Recovery Successful: {df.height} rows loaded.")
#     # ------------------------------
#     # CHR FIX (SAFE)
#     # ------------------------------
#     valid_chr = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
#     if chr_col in df.columns:
#         df = df.with_columns(
#             pl.col(chr_col)
#             .cast(pl.Utf8, strict=False)
#             .str.strip_chars()
#             .str.replace(r"(?i)^chr", "", literal=False)
#             .str.to_uppercase()
#         )
#         before_chr = df.height
#         df = df.filter(
#             pl.col(chr_col).is_in(valid_chr) &
#             pl.col(pos_col).is_not_null()
#         )
#         removed_coords = before_chr - df.height
#         print(f"🧩 Removed {removed_coords:,} invalid coordinate rows")
#     # ------------------------------
#     # POS FIX
#     # ------------------------------
#     if pos_col in df.columns:
#         df = df.with_columns(
#             pl.col(pos_col).cast(pl.Int64, strict=False)
#         )
#     # ------------------------------
#     # ALLELE FILTER (FIXED LOGIC)
#     # ------------------------------
#     if ea_col in df.columns and oa_col in df.columns:
#         before = df.height
#         valid_mask = (
#             pl.col(ea_col).str.to_uppercase().str.contains(r"^[ACGT]+$") &
#             pl.col(oa_col).str.to_uppercase().str.contains(r"^[ACGT]+$")
#         )
#         invalid_df = df.filter(~valid_mask)
#         df = df.filter(valid_mask)
#         non_standard_allele_count = invalid_df.height
#         gwas_name = sample_column_dict.get("gwas_outputname", "gwas")
#         invalid_file = os.path.join(output_dir, f"{gwas_name}_invalid_alleles.tsv")
#         if invalid_df.height > 0:
#             invalid_df.write_csv(invalid_file, separator="\t")
#             print(f"🧾 Invalid alleles saved → {invalid_file}")
#     print(f"🧩 Removed {non_standard_allele_count:,} invalid allele rows")
#     # ------------------------------
#     # DUPLICATES (FIXED)
#     # ------------------------------
#     subset_cols = [chr_col, pos_col, ea_col, oa_col]
#     if all(c in df.columns for c in subset_cols):
#         before = df.height
#         dup_df = df.filter(pl.struct(subset_cols).is_duplicated())
#         df = df.unique(subset=subset_cols, keep="first")
#         removed_duplicates_count = before - df.height
#         gwas_name = sample_column_dict.get("gwas_outputname", "gwas")
#         dup_file = os.path.join(output_dir, f"{gwas_name}_duplicates.tsv")
#         if dup_df.height > 0:
#             dup_df.write_csv(dup_file, separator="\t")
#             print(f"🧾 Duplicates saved → {dup_file}")
#     print(f"🧩 Removed {removed_duplicates_count:,} duplicate variants")
#     # ------------------------------
#     # FINAL CHECK
#     # ------------------------------
#     if df.height < 5:
#         print(f"⚠️ Too few rows: {df.height}")
#         return (
#             df, shell_cnt, python_read_count,
#             sample_column_dict,
#             removed_coords,
#             non_standard_allele_count,
#             removed_duplicates_count
#         )
#     if abs(shell_cnt - python_read_count) > 1:
#         print("⚠️ Shell vs Polars row mismatch")
#     return (
#         df, shell_cnt, python_read_count,
#         sample_column_dict,
#         removed_coords,
#         non_standard_allele_count,
#         removed_duplicates_count
#     )




# # ----------------------------------------------------------------------
# # 1. Helper: correct opener
# # ----------------------------------------------------------------------
# def _opener(path: str, mode: str = "rt"):
#     p = Path(path)
#     ext = p.suffix.lower()
#     # Use standard library openers based on extension
#     if ext == ".gz":
#         return gzip.open(path, mode)
#     if ext == ".bz2":
#         return bz2.open(path, mode)
#     if ext in {".xz", ".lzma"}:
#         return lzma.open(path, mode)
#     if ext == ".zip":
#         return None 
#     return open(path, mode)


# ----------------------------------------------------------------------
# 2. Delimiter detection – look at the *header* line
# ----------------------------------------------------------------------
# def detect_delimiter(path: str) -> str:
#     opener = _opener(path)
#     if opener is None: # Zip
#         try:
#             with zipfile.ZipFile(path) as zf:
#                 csvs = [f for f in zf.namelist() if f.lower().endswith(('.csv', '.tsv'))]
#                 if not csvs: return "\t"
#                 txt = zf.read(csvs[0]).decode(errors="ignore")
#                 header = next((l for l in txt.splitlines() if not l.startswith("##")), "")
#         except:
#             return "\t"
#     else:
#         try:
#             with opener(path, "rt", errors="ignore") as f:
#                 header = next((l for l in f if not l.startswith("##")), "")
#         except:
#             return "\t"
#     if not header.strip(): return "\t"
#     try:
#         return csv.Sniffer().sniff(header, delimiters="\t,; ").delimiter
#     except csv.Error:
#         return "\t"



# # ----------------------------------------------------------------------
# # 3. Detect ## metadata
# # ----------------------------------------------------------------------
# def has_double_hash(path: str) -> bool:
#     opener = _opener(path)
#     if opener is None: # Zip
#         try:
#             with zipfile.ZipFile(path) as zf:
#                 csvs = [f for f in zf.namelist() if f.lower().endswith(('.csv', '.tsv'))]
#                 if not csvs: return False
#                 txt = zf.read(csvs[0]).decode(errors="ignore")
#                 return any(l.startswith("##") for l in txt.splitlines()[:50])
#         except:
#             return False
#     try:
#         with opener(path, "rt", errors="ignore") as f:
#             for _ in range(50):
#                 line = f.readline()
#                 if not line: break
#                 if line.startswith("##"): return True
#         return False
#     except:
#         return False


# # ----------------------------------------------------------------------
# # 4. Count data lines – **safe** shell command (list form)
# # ----------------------------------------------------------------------
# def count_data_lines(path: str, skip_hash: bool) -> int:
#     if Path(path).suffix.lower() == ".zip":
#         try:
#             with zipfile.ZipFile(path) as zf:
#                 csvs = [f for f in zf.namelist() if f.lower().endswith(('.csv', '.tsv'))]
#                 if not csvs: return 0
#                 txt = zf.read(csvs[0]).decode(errors="ignore")
#                 lines = txt.splitlines()
#                 if skip_hash:
#                     lines = [l for l in lines if not l.startswith("##")]
#                 return len(lines)
#         except Exception as e:
#             print(f"Zip count failed: {e}")
#             return 0
#     # macOS-safe: use "< path" for compressed files
#     if path.endswith(".gz"):
#         decomp = "zcat <"
#     elif path.endswith((".xz", ".lzma")):
#         decomp = "xzcat <"
#     else:
#         decomp = "cat"
#     cmd_str = f"{decomp} {shlex.quote(path)}"
#     if skip_hash:
#         cmd_str += " | grep -v '^##'"
#     cmd_str += " | wc -l"
#     try:
#         result = subprocess.run(cmd_str, shell=True, capture_output=True, text=True, check=True)
#         return int(result.stdout.strip())
#     except Exception as e:
#         print(f"Shell count failed ({e}) — using Polars count only")
#         return 0



# # ----------------------------------------------------------------------
# # 5. Helper: Clean Stream + Count
# # ----------------------------------------------------------------------
# def vcf_to_polars_stream(path: str) -> Tuple[BytesIO, int]:
#     """
#     1. Removes ## metadata lines.
#     2. Counts valid lines on the fly.
#     """
#     buf = BytesIO()
#     row_count = 0
#     opener = _opener(path) 
#     if opener is None and not path.endswith('.zip'): opener = open
#     try:
#         with opener(path, "rt", errors="ignore") as f:
#             for line in f:
#                 if line.startswith("##"):
#                     continue
#                 buf.write(line.encode("utf-8"))
#                 row_count += 1
#     except Exception as e:
#         print(f"   ❌ Error streaming file: {e}")
#         return BytesIO(), 0
#     buf.seek(0)
#     data_count = max(0, row_count - 1)
#     return buf, data_count


# # ======================================================================
# # HELPER 1: DataFrame Cleanup (Headers, Types, Chromosome)
# # ======================================================================
# def clean_dataframe(df: pl.DataFrame, chr_col: str) -> pl.DataFrame:
#     """
#     Applied immediately after ANY read (Attempt 1 or Attempt 2).
#     Trims whitespace, fixes types, and standardizes chromosome column.
#     """
#     if df.height == 0: return df
#     # 1. Clean Headers
#     df = df.rename({c: c.strip() for c in df.columns})
#     # 2. Clean String Values
#     for col in df.columns:
#         if df[col].dtype in (pl.Utf8, pl.String):
#             df = df.with_columns(pl.col(col).str.strip_chars())
#     # 3. Type Conversion
#     valid_cols = []
#     for col in df.columns:
#         if col == chr_col:
#             valid_cols.append(pl.col(col)) # Skip casting for CHR
#             continue
#         numeric_series = df[col].cast(pl.Float64, strict=False)
#         if numeric_series.null_count() > df[col].null_count():
#             valid_cols.append(pl.col(col)) # Keep original if casting fails
#         else:
#             valid_cols.append(numeric_series)
#     df = df.with_columns(valid_cols)
#     # 4. Chromosome Fix (12.0 -> 12)
#     if chr_col and chr_col in df.columns:
#         df = df.with_columns(
#             pl.coalesce(
#                 pl.col(chr_col).cast(pl.Float64, strict=False).cast(pl.Int64, strict=False).cast(pl.String),
#                 pl.col(chr_col)
#             ).alias(chr_col)
#         )
#     return df

# # ======================================================================
# # HELPER 2: Validation (Is the file broken?)
# # ======================================================================
# def is_bad_parsing(df: pl.DataFrame, expected_cols: int) -> bool:
#     if df.width < 2 and expected_cols > 1: return True # Swallowed
#     if df.width > (expected_cols + 5): return True     # Ghost columns
#     if df.height > 0:
#         null_ratio = df.null_count().sum_horizontal().sum() / (df.height * df.width)
#         if null_ratio > 0.5: return True # Null flood
#     return False

# # ======================================================================
# # HELPER 3: Repair (Linux Shell Fix)
# # ======================================================================
# def normalize_whitespace_file(input_path: str, output_dir: str) -> str:
#     base_name = Path(input_path).name
#     while Path(base_name).suffix: base_name = Path(base_name).stem
#     clean_path = os.path.join(output_dir, f"{base_name}_cleaned.tsv")
#     if input_path.lower().endswith(".gz"): cat_cmd = "gzip -dc"
#     elif input_path.lower().endswith(".zip"): cat_cmd = "unzip -p"
#     elif input_path.lower().endswith((".xz", ".lzma")): cat_cmd = "xz -dc"
#     else: cat_cmd = "cat"
#     # Fix: zcat -> remove ## -> squeeze spaces -> tab
#     cmd = (
#         f"{cat_cmd} {shlex.quote(input_path)} | "
#         f"grep -v '^##' | "
#         f"tr -s '[:blank:]' '\\t' > {shlex.quote(clean_path)}"
#     )
#     try:
#         subprocess.run(cmd, shell=True, check=True)
#         if os.path.exists(clean_path) and os.path.getsize(clean_path) > 0:
#             return clean_path
#     except: pass
#     return ""

# # ======================================================================
# # HELPER 4: Header Counter
# # ======================================================================
# def get_expected_cols(path: str) -> int:
#     try:
#         # (Same opener logic as before) ...
#         opener = gzip.open if path.lower().endswith((".gz", ".bgz")) else open
#         if path.lower().endswith(".zip"): return 0
        
#         with opener(path, "rt", errors="ignore") as f:
#             for line in f:
#                 if not line.startswith("##") and line.strip():
#                     # FIX: Replace commas (and semicolons) with spaces first
#                     clean_line = line.replace(",", " ").replace(";", " ").replace("\t", " ")
#                     return len(clean_line.split())
#     except: 
#         return 0
#     return 0


# def clean_and_optimize_dataframe(df: pl.DataFrame) -> pl.DataFrame:
#     """
#     1. Trims whitespace from ALL string columns.
#     2. Replaces empty strings "" with null (optional, but good practice).
#     3. Attempts to convert string columns that look like numbers into actual Integers/Floats.
#     """
#     return (
#         df.lazy()
#         # Step 1: Trim whitespace from all string columns
#         .with_columns(
#             pl.col(pl.String).str.strip_chars()
#         )
#         # Step 2 (Optional): Convert empty strings to null
#         .with_columns(
#             pl.col(pl.String).replace("", None)
#         )
#         # Step 3: Collect immediately to perform type inference
#         .collect()
#     )

# # ======================================================================
# # MAIN FUNCTION
# # ======================================================================
# def read_sumstats(sumstat_file: str, output_dir: str, sample_column_dict: dict) -> Tuple[pl.DataFrame, int, int]:
#     chr_col = sample_column_dict.get("chr_col")
#     chr_pos_col = sample_column_dict.get("chr_pos_col")
#     os.makedirs(output_dir, exist_ok=True)
#     # ------------------------------
#     # Pre-flight Checks
#     # ------------------------------
#     sumstat_file, delim = get_polars_separator(sumstat_file, output_dir)
#     skip_hash = has_double_hash(sumstat_file)
#     expected_cols = get_expected_cols(sumstat_file)
#     df = pl.DataFrame()
#     shell_cnt = 0
#     # ------------------------------
#     # ATTEMPT 1: Optimistic Read
#     # ------------------------------
#     try:
#         suffix = Path(sumstat_file).suffix.lower()
#         if suffix == ".zip":
#             with zipfile.ZipFile(sumstat_file) as zf:
#                 csvs = [f for f in zf.namelist() if f.lower().endswith(('.csv', '.tsv'))]
#                 if not csvs:
#                     raise ValueError("No CSV/TSV in .zip")
#                 txt = zf.read(csvs[0]).decode(errors="ignore")
#                 lines = txt.splitlines()
#                 if skip_hash:
#                     lines = [l for l in lines if not l.startswith("##")]
#                 shell_cnt = max(0, len(lines) - 1)
#                 df = pl.read_csv(
#                     BytesIO("\n".join(lines).encode("utf-8")),
#                     separator=delim,
#                     has_header=True,
#                     quote_char=None,
#                     infer_schema_length=0,
#                     ignore_errors=True,
#                     null_values=["NA", "na", ".", ""],
#                     truncate_ragged_lines=True,
#                 )
#         else:
#             if skip_hash:
#                 # vcf_to_polars_stream should handle decompression internally 
#                 # or be updated to use _opener
#                 source, shell_cnt = vcf_to_polars_stream(sumstat_file)
#                 df = pl.read_csv(
#                     source,
#                     separator=delim,
#                     has_header=True,
#                     quote_char=None,
#                     infer_schema_length=0,
#                     ignore_errors=True,
#                     null_values=["NA", "na", ".", ""],
#                     truncate_ragged_lines=True,
#                 )
#             else:
#                 shell_cnt = count_data_lines(sumstat_file, skip_hash=False)
#                 if shell_cnt >= 1:
#                     shell_cnt -= 1
#                 try:
#                     df = pl.read_csv(
#                         sumstat_file,
#                         separator=delim,
#                         has_header=True,
#                         quote_char=None,
#                         infer_schema_length=0,
#                         ignore_errors=True,
#                         null_values=["NA", "na", ".", ""],
#                         truncate_ragged_lines=True,
#                     )
#                 except (OSError, ValueError):
#                     # fallback for unusual compression / broken files
#                     with _opener(sumstat_file, mode="rb") as f:
#                         df = pl.read_csv(
#                             f,
#                             separator=delim,
#                             has_header=True,
#                             quote_char=None,
#                             infer_schema_length=0,   # 🔴 keep SAME
#                             ignore_errors=True,
#                             null_values=["NA", "na", ".", ""],
#                             truncate_ragged_lines=True,
#                         )
#         df, sample_column_dict = fix_chr_pos_allele_column(
#             chromosome="All_Chrs",
#             df=df,
#             sample_column_dict=sample_column_dict,
#             drop_mt=True
#         )
#         df = clean_and_optimize_dataframe(df)
#     except Exception as e:
#         print(f"\t\t\t ⚠️ Attempt 1 Error: {e}")
#         df = pl.DataFrame()
#     # ------------------------------
#     # ATTEMPT 2: Repair
#     # ------------------------------
#     if is_bad_parsing(df, expected_cols):
#         print(f"\t\t\t⚠️  Standard parsing failed (Rows: {df.height}, Cols: {df.width}). Switching to Repair Mode...")
#         clean_file = normalize_whitespace_file(sumstat_file, output_dir)
#         if clean_file:
#             print(f"\t\t\t🔄 Reading repaired file: {clean_file}")
#             try:
#                 df = pl.read_csv(
#                     clean_file,
#                     separator="\t",
#                     has_header=True,
#                     infer_schema_length=0,
#                     null_values=["NA", "na", ".", ""],
#                     truncate_ragged_lines=True
#                 )
#                 df, sample_column_dict = fix_chr_pos_allele_column(
#                     chromosome="All_Chrs",
#                     df=df,
#                     sample_column_dict=sample_column_dict,
#                     drop_mt=True
#                 )
#                 df = clean_dataframe(df, chr_col)
#                 shell_cnt = df.height
#                 print(f"\t\t\t✅ Recovery Successful: {df.height} rows loaded.")
#             except Exception as e:
#                 print(f"❌ Failed to read repaired file: {e}")
#         else:
#             print("❌ Repair failed (could not create clean file).")
    
#     python_read_count=df.height
#     # ------------------------------
#     # FINAL COLUMN FIXES
#     # ------------------------------ 
#     chr_col = sample_column_dict.get("chr_col")
#     pos_col = sample_column_dict.get("pos_col")
#     ea_col = sample_column_dict.get("ea_col", "effect_allele")
#     oa_col = sample_column_dict.get("oa_col", "other_allele")
#     # 🔹 Fix Chromosome (chr-safe, X/Y safe)
#     if chr_col and chr_col in df.columns:
#         df = df.with_columns(
#             pl.col(chr_col)
#             .cast(pl.Utf8, strict=False)
#             .fill_null("")
#             .str.strip_chars()
#             # 1. Remove the .0 suffix first (turns "11.0" into "11")
#             .str.replace(r"\.0$", "", literal=False)
#             # 2. Remove "chr" prefix
#             .str.replace(r"(?i)^chr", "", literal=False)
#             .str.to_uppercase()
#             .alias(chr_col)
#         )
#         df = df.with_columns(
#             pl.when(pl.col(chr_col).is_in(["X", "23"])).then(pl.lit("X"))
#             .when(pl.col(chr_col).is_in(["Y", "24"])).then(pl.lit("Y"))
#             .when(pl.col(chr_col).is_in(["MT", "M", "25"])).then(pl.lit("MT"))
#             .otherwise(pl.col(chr_col))
#             .alias(chr_col)
#         )
#         # 🔴 Validate chromosomes
#         valid_chr = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
#         # Now "11" will be in valid_chr, but "11.0" would not have been.
#         invalid_chr_df = df.filter(~pl.col(chr_col).is_in(valid_chr))
#         if invalid_chr_df.height > 0:
#             print("❌ Invalid chromosome values detected:")
#             print(invalid_chr_df.select(chr_col).unique())
        
#         df = df.filter(
#             pl.col(chr_col).is_not_null() & 
#             pl.col(pos_col).is_not_null() &
#             pl.col(chr_col).is_in(valid_chr)
#         )
#         removed_coords = python_read_count - df.height
#         print(f"🧩 Removed {removed_coords:,} variants due to NULL coordinates or invalid chromosomes.")
    
#         # 🔴 Optional (recommended): autosomes only
#         # AUTOSOMES_ONLY = True
#         # if AUTOSOMES_ONLY:
#         #     df = df.filter(pl.col(chr_col).is_in([str(i) for i in range(1, 23)]))
#     df = clean_and_optimize_dataframe(df)
#     if pos_col in df.columns:
#         df = df.with_columns(
#             pl.col(pos_col)
#             .cast(pl.Float64, strict=False)
#             .cast(pl.Int64, strict=False)
#         )
#     # --- NEW: ALLELE FILTERING & LOGGING ---
#     # Filter for standard DNA bases (preventing the gwas2vcf 'null allele' crash)
#     df = df.filter(
#         pl.col(ea_col).str.to_uppercase().str.contains(r"^[ACGT]+$") & 
#         pl.col(oa_col).str.to_uppercase().str.contains(r"^[ACGT]+$")
#     )
#     non_standard_allele_count = python_read_count - df.height 
#     print(f"🧩 Filtered out {non_standard_allele_count:,} variants with non-standard alleles (I, D, ., etc.)")
#     # ------------------------------ 
#     # FINAL QC
#     # ------------------------------
#     if df.height < 5:
#         print(f"⚠️  WARNING: File is empty or has too few rows ({df.height}). Exiting.")
#         return df, shell_cnt, python_read_count, sample_column_dict, removed_coords, non_standard_allele_count
#     if abs(shell_cnt - python_read_count) > 1:
#         print(f"⚠️ Warning: Mismatch between line-count ({shell_cnt}) and Polars rows ({python_read_count}).")
#     # ------------------------------ 
#     # SORT BY CHROMOSOME AND POSITION
#     # ------------------------------ 
#     with_dup_count=df.height
#     chr_col = sample_column_dict.get("chr_col")
#     pos_col = sample_column_dict.get("pos_col")
#     ea_col = sample_column_dict.get("ea_col")
#     oa_col = sample_column_dict.get("oa_col")
#     df = df.sort([chr_col, pos_col])    
#     # Define the columns to check for duplicates
#     subset_cols = [chr_col, pos_col, ea_col, oa_col]
#     # Remove duplicates, keeping the first occurrence
#     df = df.unique(subset=subset_cols, keep="first")
#     removed_duplicates_count = with_dup_count - df.height
#     print(f"🧩 Filtered out {removed_duplicates_count:,} variants with duplicate variants[chr,pos,ea,oa]")
    
#     indent = "\t" * 3
#     print(f"{indent} First 5 rows of sumstats {sumstat_file}:")
#     preview_str = str(df.head())
#     preview_str = "\n".join(indent + line for line in preview_str.splitlines())
#     print(preview_str)
#     print(f"{indent}")
#     print(f"{indent} Columns of sumstats {sumstat_file}:")
#     cols = ", ".join(df.columns)
#     print(f"{indent}[COLUMNS] ({len(df.columns)} total):\n{indent}{cols}") 
#     print(" \n")
#     return df, shell_cnt, python_read_count, sample_column_dict, removed_coords, non_standard_allele_count,removed_duplicates_count



def check_file_truncation(file_path):
    # -------------------------
    # 1. HARD FAIL
    # -------------------------
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"❌ FILE NOT FOUND: {file_path}")
    if os.path.getsize(file_path) == 0:
        raise ValueError(f"❌ EMPTY FILE: {file_path}")
    warnings = []
    # -------------------------
    # 2. GZIP integrity ONLY
    # -------------------------
    if file_path.endswith(".gz"):
        if subprocess.run(["gzip", "-t", file_path]).returncode != 0:
            warnings.append("GZIP corrupted or truncated")
    # -------------------------
    # 3. VERY LIGHT tail check
    # -------------------------
    try:
        cmd = f"zcat {shlex.quote(file_path)} | tail -n 1" if file_path.endswith(".gz") else f"tail -n 1 {shlex.quote(file_path)}"
        last_line = subprocess.check_output(cmd, shell=True).decode().strip()
        # ONLY strong signal → abrupt ending
        if last_line.endswith("."):
            warnings.append("Last line ends abruptly → possible truncation")
    except:
        pass  # ignore silently (minimal philosophy)
    # -------------------------
    # 4. PRINT only if needed
    # -------------------------
    if warnings:
        print("\n⚠️ Possible file truncation detected:")
        for w in warnings:
            print(f"   - {w}")
    return warnings



def find_resource_file_path(input_file, resource_folder, grch_version, chromosome):
    """
    Resolve user EAF resource path.

    Logic:
      1. If input_file == "NA" -> return "NA"
      2. If input_file is a full existing path -> return it
      3. Otherwise build:
           {resource_folder}/{grch_version}/external_af/vcf_files/
           {grch_version}_{input_file}_freq_chr{chromosome}.vcf.gz
         If exists -> return it
      4. Else -> return "NA"
    """

    # 1) User explicitly disabled EAF
    if input_file == "NA":
        return "NA"

    # 2) User gave a full path
    if os.path.exists(input_file):
        return input_file
    
    if os.path.exists(f"{input_file}_chr{chromosome}.tsv.gz"):
        return f"{input_file}_chr{chromosome}.tsv.gz"
    
    # 3) Construct the standard resource path (EXACT pattern you gave)
    constructed_path = (
        f"{resource_folder}/{grch_version}/external_af/vcf_files/"
        f"{grch_version}_{input_file}_freq_chr{chromosome}.vcf.gz"
    )

    if os.path.exists(constructed_path):
        return constructed_path

    # 4) Nothing found
    print(f"⚠️ Warning: No valid EAF file found for '{input_file}'.")
    print(f"   Tried full path      : {input_file}")
    print(f"   Tried constructed path: {input_file}_chr{chromosome}.tsv.gz")
    return "NA"



def read_config(csv_path):
    """
    Read a configuration file (CSV or TSV) and return a list of dictionaries.
    Automatically detects delimiter and cleans whitespace.

    Parameters
    ----------
    csv_path : str
        Path to configuration file.

    Returns
    -------
    list[dict]
        List of configuration dictionaries.
    """

    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"❌ The file '{csv_path}' was not found.")

    # -------------------------------
    # Auto-detect separator
    # -------------------------------
    df = pd.read_csv(csv_path, sep=None, engine="python")

    # Replace NA values
    df = df.fillna("NA")

    # -------------------------------
    # Clean column names
    # -------------------------------
    df.columns = df.columns.str.strip()

    # -------------------------------
    # Clean string cells
    # -------------------------------
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    # -------------------------------
    # Convert to list of dictionaries
    # -------------------------------
    cfg_list = df.to_dict(orient="records")

    # Final safety strip
    cfg_list = [
        {k.strip(): (v.strip() if isinstance(v, str) else v) for k, v in cfg.items()}
        for cfg in cfg_list
    ]

    return cfg_list

def load_default_config(filename: str):
    """
    Load a YAML config file.

    New behaviour:
        • If `filename` is a valid filesystem path → load directly.
        • Else → try loading from postgwas.config package resource.
    """

    # 1️⃣ Try direct filesystem path
    if os.path.exists(filename) and os.path.isfile(filename):
        try:
            with open(filename, "r") as f:
                return yaml.safe_load(f)
        except Exception as e:
            raise RuntimeError(
                f"❌ ERROR: Failed to load config from path:\n  {filename}\nReason: {e}"
            )

    # 2️⃣ Fallback to package resource
    try:
        with resources.open_text("postgwas.config", filename) as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        raise FileNotFoundError(
            f"❌ Config file '{filename}' not found.\n"
            f"Tried filesystem path AND postgwas.config resource."
        )
    except Exception as e:
        raise RuntimeError(
            f"❌ ERROR: Failed to load YAML from package resource '{filename}'.\nReason: {e}"
        )
