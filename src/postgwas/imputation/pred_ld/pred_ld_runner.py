#!/usr/bin/env python3

from __future__ import annotations

import functools
import os, re,subprocess,time,shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,ProcessPoolExecutor, as_completed
import pandas as pd
import numpy as np
import polars as pl
from scipy.stats import norm
from postgwas.utils.main import safe_thread_count 
from typing import List, Optional, Tuple
import threading,psutil,sys
from datetime import datetime
    
"""
pred_ld_runner.py — multiprocessing-safe PRED-LD engine for PostGWAS

Designed to run inside Docker AND locally on macOS/Linux.

Features:
---------
✓ macOS-safe multiprocessing (spawn)
✓ Linux-safe multiprocessing (fork)
✓ Worker function defined at top-level (pickle-safe)
✓ Automatic pred_ld.py discovery
✓ Interleaved chromosome scheduling (big + small → reduced peak RAM)
✓ Workspace isolation in HOME (always writable in Docker)
✓ Merges logs and validates outputs
✓ Auto thread reduction using safe_thread_count()
"""

# =====================================================================
# 1. PRED-LD PARALLEL RUNNER
# =====================================================================



# =============================================================================
#  MEMORY HELPERS
# =============================================================================

def _mem_from_psutil() -> Optional[Tuple[float, float]]:
    """Try to read total/free RAM using psutil. Returns (total_gb, free_gb)."""
    try:
        vm = psutil.virtual_memory()
        total_gb = vm.total / (1024 ** 3)
        free_gb = vm.available / (1024 ** 3)
        return total_gb, free_gb
    except Exception:
        return None


def _mem_from_proc() -> Optional[Tuple[float, float]]:
    """Linux fallback: parse /proc/meminfo. Returns (total_gb, free_gb)."""
    try:
        meminfo = Path("/proc/meminfo")
        if not meminfo.exists():
            return None

        info = {}
        with meminfo.open() as f:
            for line in f:
                parts = line.split(":")
                if len(parts) != 2:
                    continue
                key = parts[0].strip()
                val = parts[1].strip().split()[0]
                try:
                    info[key] = float(val)
                except ValueError:
                    continue

        # Values in kB → convert to GB
        total_kb = info.get("MemTotal")
        # Use MemAvailable if present; else fall back to (MemFree + Cached)
        avail_kb = info.get("MemAvailable")
        if avail_kb is None:
            avail_kb = info.get("MemFree", 0.0) + info.get("Cached", 0.0)

        if total_kb is None or avail_kb is None:
            return None

        total_gb = total_kb / (1024 ** 2)
        free_gb = avail_kb / (1024 ** 2)
        return total_gb, free_gb
    except Exception:
        return None


def get_memory_status() -> Tuple[Optional[float], Optional[float]]:
    """
    Unified RAM checker.

    Returns:
        (total_gb, free_gb) or (None, None) if cannot be determined.
    """
    res = _mem_from_psutil()
    if res is not None:
        return res
    res = _mem_from_proc()
    if res is not None:
        return res
    return None, None


# =============================================================================
#  CHROMOSOME HELPERS
# =============================================================================

BIG_CHROMS = {"1", "2", "3", "4"}


def is_big_chr(chr_id: str) -> bool:
    """Return True if chromosome is considered 'big' for scheduling."""
    chr_id = chr_id.upper().replace("CHR", "")
    try:
        n = int(chr_id)
        return str(n) in BIG_CHROMS
    except ValueError:
        return False


def build_interleaved_chr_order(chromosomes: List) -> List[str]:
    """
    Use a fixed, RAM-balanced chromosome execution order:

    1, 22, 21, 20, 2, 19, 18, 17, 3, 16, 15, 14,
    4, 13, 12, 11, 5, 10, 9, 8, 7, 6, X

    - chrX always last (if present)
    - Missing chromosomes are skipped safely
    """
    desired_order = [
        "1", "22", "21", "20","19","18", 
        "2", "17", "16", "15", "14", 
        "3", "13","12", "11", 
        "4", "10", "9",  "8",
        "5", "7", "6", "X"
    ]
    # Normalize input
    chroms = {str(c).upper().replace("CHR", "") for c in chromosomes}
    # Preserve order, include only chromosomes the user actually provided
    final_order = [c for c in desired_order if c in chroms]
    return final_order


# =============================================================================
#  MAIN PRED-LD PARALLEL RUNNER
# ============================================================================
def run_pred_ld_parallel(
    predld_input_dir: str,
    output_folder: str,
    output_prefix: str,
    pred_ld_ref: str,
    chromosomes: List = list(range(1, 23)) + ["X"],
    r2threshold: float = 0.8,
    maf: float = 0.001,
    population: str = "EUR",
    ref: str = "TOP_LD",
    threads: int = 6,
) -> bool:
    """
    Parallel PRED-LD runner with logging to files.
    - Standard output is silenced.
    - Logs are saved to {output_prefix}_chr{ID}_predld.log.
    - Global info saved to {output_prefix}_master.log.
    - Screen output: Only a final summary of failed chromosomes.
    """

    # -------------------------------------------------------------------------
    # Helper: Write to log file with timestamp
    # -------------------------------------------------------------------------
    def write_log(file_path: Path, message: str):
        """Appends a timestamped message to the specified log file."""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        try:
            with file_path.open("a") as f:
                f.write(f"[{timestamp}] {message}\n")
        except Exception:
            pass

    # -------------------------------------------------------------------------
    # 0. Locate pred_ld.py dynamically
    # -------------------------------------------------------------------------
    module_root = Path(__file__).resolve().parent
    pred_ld_script = module_root / "pred_ld.py"
    if not pred_ld_script.exists():
        raise FileNotFoundError(f"ERROR: pred_ld.py not found at {pred_ld_script}")

    # -------------------------------------------------------------------------
    # 1. Normalize paths & Setup Master Log
    # -------------------------------------------------------------------------
    predld_input_dir = Path(predld_input_dir).expanduser().resolve()
    output_folder = Path(output_folder).expanduser().resolve()
    pred_ld_ref = Path(pred_ld_ref).expanduser().resolve()

    output_folder.mkdir(parents=True, exist_ok=True)
    
    master_log_file = output_folder / f"{output_prefix}_master.log"
    if master_log_file.exists():
        master_log_file.unlink()

    # -------------------------------------------------------------------------
    # 2. Workspace Setup
    # -------------------------------------------------------------------------
    workspace = output_folder / ".predld_work"
    if workspace.exists():
        shutil.rmtree(workspace, ignore_errors=True)
    workspace.mkdir(parents=True, exist_ok=True)

    local_ref = workspace / "ref"
    shutil.copytree(pred_ld_ref, local_ref, dirs_exist_ok=True)

    # -------------------------------------------------------------------------
    # 3. Threads / workers & chromosome order
    # -------------------------------------------------------------------------
    threads = safe_thread_count(threads, gb_per_thread=30)
    #max_workers = min(2, max(1, threads))
    max_workers = threads
    chr_order = build_interleaved_chr_order(chromosomes)

    write_log(master_log_file, f"Running PRED-LD with {max_workers} parallel workers")
    write_log(master_log_file, f"Input dir : {predld_input_dir}")
    write_log(master_log_file, f"Output dir: {output_folder}")
    write_log(master_log_file, f"Workspace : {workspace}")
    write_log(master_log_file, f"Ref LD    : {pred_ld_ref}")
    write_log(master_log_file, f"Chromosomes: {', '.join(chr_order)}")

    # -------------------------------------------------------------------------
    # 4. Shared state for scheduling
    # -------------------------------------------------------------------------
    lock = threading.Lock()
    running_chroms: set[str] = set()
    running_big: bool = False
    succeeded: list[str] = []
    failed: list[Tuple[str, str]] = []

    def wait_until_allowed(chr_id: str, log_file: Path) -> None:
        nonlocal running_big
        while True:
            with lock:
                total_gb, free_gb = get_memory_status()
                ram_known = (total_gb is not None and free_gb is not None)

                if ram_known:
                    threshold = min(20.0, total_gb * 0.80)
                    free_ok = free_gb >= threshold
                else:
                    threshold = None
                    free_ok = True

                this_big = is_big_chr(chr_id)
                big_ok = (not this_big) or (not running_big)

                if free_ok and big_ok:
                    running_chroms.add(chr_id)
                    if this_big:
                        running_big = True
                    
                    if ram_known:
                        write_log(log_file, f"🚀 Started chr{chr_id}: free RAM {free_gb:.1f} GB (threshold {threshold:.1f} GB)")
                    else:
                        write_log(log_file, f"🚀 Started chr{chr_id}: RAM unknown.")
                    print(f"            🚀 Started imputation for chr{chr_id} using PRED-LD")
                    return

                if ram_known and not free_ok:
                    write_log(log_file, f"⏳ Waiting for RAM: free={free_gb:.1f} GB, needed={threshold:.1f} GB")
                if this_big and not big_ok:
                    write_log(log_file, "⏳ Waiting: another big chromosome is running.")

            time.sleep(10)

    def release_chr(chr_id: str) -> None:
        nonlocal running_big
        with lock:
            running_chroms.discard(chr_id)
            if is_big_chr(chr_id):
                if not any(is_big_chr(c) for c in running_chroms):
                    running_big = False

    # -------------------------------------------------------------------------
    # 5. Worker (SILENT MODE)
    # -------------------------------------------------------------------------
    def worker(chr_id: str) -> None:
        chr_input = predld_input_dir / f"{output_prefix}_chr{chr_id}_predld_input.tsv"
        log_file = output_folder / f"{output_prefix}_chr{chr_id}_predld.log"

        if log_file.exists():
            log_file.unlink()
        write_log(log_file, f"--- Processing Chromosome {chr_id} ---")

        if not chr_input.exists():
            msg = f"PRED-LD input missing for chr{chr_id}: {chr_input}"
            write_log(log_file, f"⚠️ SKIPPED: {msg}")
            # REMOVED: print(f"⚠️ chr{chr_id} skipped...")
            failed.append((chr_id, msg))
            return

        wait_until_allowed(chr_id, log_file)

        try:
            cmd = [
                "python", str(pred_ld_script),
                "--file-path", str(chr_input),
                "--pop", population,
                "--ref", ref,
                "--ref_dir", str(local_ref),
                "--r2threshold", str(r2threshold),
                "--maf", str(maf),
                "--out_dir", str(output_folder),
            ]

            with log_file.open("a") as log_f:
                log_f.write("\n--- SUBPROCESS OUTPUT START ---\n")
                log_f.flush()
                proc = subprocess.run(cmd, stdout=log_f, stderr=log_f)
                log_f.write("\n--- SUBPROCESS OUTPUT END ---\n")

            if proc.returncode != 0:
                reason = f"     PRED-LD failed for chr{chr_id} (exit={proc.returncode}).        See log: {log_file}"
                write_log(log_file, f"❌ {reason}")
                # REMOVED: print(f"❌ {reason}")
                failed.append((chr_id, reason))
                return

            expected_out = output_folder / f"imputation_results_chr{chr_id}.txt"
            if not expected_out.exists():
                reason = f"PRED-LD output missing for chr{chr_id}. See log: {log_file}"
                write_log(log_file, f"❌ {reason}")
                # REMOVED: print(f"❌ {reason}")
                failed.append((chr_id, reason))
                return

            write_log(log_file, f"✅ chr{chr_id} finished successfully.")
            succeeded.append(chr_id)

        except Exception as e:
            reason = f"Exception processing chr{chr_id}: {str(e)}"
            write_log(log_file, f"❌ {reason}")
            # REMOVED: print(f"❌ {reason}")
            failed.append((chr_id, reason))
        finally:
            release_chr(chr_id)

    # -------------------------------------------------------------------------
    # 6. Launch workers
    # -------------------------------------------------------------------------
    t0 = time.time()
    threads_list: List[threading.Thread] = []

    for chr_id in chr_order:
        t = threading.Thread(target=worker, args=(chr_id,), daemon=True)
        threads_list.append(t)

    active: List[threading.Thread] = []
    for t in threads_list:
        while True:
            active = [thr for thr in active if thr.is_alive()]
            if len(active) < max_workers:
                t.start()
                active.append(t)
                break
            time.sleep(2)

    for t in threads_list:
        t.join()

    # -------------------------------------------------------------------------
    # 7. Merge Logs
    # -------------------------------------------------------------------------
    merged_log = output_folder / f"{output_prefix}_predld_combined.log"
    write_log(master_log_file, f"Merging chromosome logs into {merged_log}")
    
    with merged_log.open("w") as fout:
        if master_log_file.exists():
            with master_log_file.open() as f:
                fout.write(f.read())
                fout.write("\n" + "="*80 + "\n")
        
        for lf in sorted(output_folder.glob(f"{output_prefix}_chr*_predld.log")):
            with lf.open() as f:
                fout.write(f"\nLOG: {lf.name}\n")
                fout.write(f.read())

    # -------------------------------------------------------------------------
    # 8. Summary (Screen Output Logic)
    # -------------------------------------------------------------------------
    elapsed_min = (time.time() - t0) / 60.0
    write_log(master_log_file, f"\n🎯 PRED-LD completed in {elapsed_min:.2f} minutes.")
    
    # --- CHANGED: Only print one summary line for failures ---
    if failed:
        write_log(master_log_file, "⚠️ Failed / skipped chromosomes:")
        # Sort keys numerically if possible, else string sort
        try:
            # Attempt to sort numerically (1, 2, 10 instead of 1, 10, 2)
            failed_ids = sorted([str(c) for c, r in failed], key=lambda x: int(x) if x.isdigit() else 999)
        except:
            failed_ids = sorted([str(c) for c, r in failed])
            
        # Log details to file
        for chr_id, reason in failed:
            write_log(master_log_file, f"   chr{chr_id}: {reason}")
        
        # Print CONCISE summary to screen
        print(f"\n          ❌ PRED-LD Failed for chromosomes: {', '.join(failed_ids)}")
        print(f"            See detailed logs in: {merged_log}\n")
    else:
        print("\n           ✅Summary stiatistics Imputation using PRED-LD software finished.")
        print(" ")
        print(" ")

    return bool(succeeded)


def process_pred_ld_results_all_parallel(
    folder_path: str,
    output_path: Optional[str] = None,
    gwas2vcf_resource_folder: Optional[str] = None,
    output_prefix: Optional[str] = None,
    corr_method: str = "pearson",
    threads: int = 6,
):
    """
    Production-ready PRED-LD postprocessing pipeline.
    - Parallel Polars implementation
    - Safe dtype coercion for all PRED-LD fields
    - Accurate imputed-SNP sample-size propagation
    - Duplicate SNP diagnostics + correlation
    - Final TSV + correlation summary
    - gwas2vcf config generation
    - Cleanup of intermediate files
    """
    # =====================================================================
    # Safe helpers
    # =====================================================================
    def coerce_dtypes(df: pl.DataFrame, dtype_map: dict) -> pl.DataFrame:
        """Coerce columns to specific dtypes (safe casting)."""
        for col, dtype in dtype_map.items():
            if col in df.columns:
                df = df.with_columns(pl.col(col).cast(dtype, strict=False))
        return df

    def safe_drop(df: pl.DataFrame, cols) -> pl.DataFrame:
        """
        Drop only columns that exist.
        cols can be a string or list of strings.
        """
        if isinstance(cols, str):
            cols = [cols]
        existing = [c for c in cols if c in df.columns]
        if not existing:
            return df
        return df.drop(existing)

    # =====================================================================
    # Expected column types
    # =====================================================================
    data_types = {
        "chr": pl.Utf8, "snp": pl.Utf8,
        "A1": pl.Utf8, "A2": pl.Utf8,
        "pos": pl.Float64,
        "beta": pl.Float64, "SE": pl.Float64, "z": pl.Float64,
        "imputed": pl.Float64,
        "R2": pl.Float64, "NC": pl.Float64, "SS": pl.Float64,
        "AF": pl.Float64, "LP": pl.Float64, "SI": pl.Float64,
    }

    info_types = {
        "pos1": pl.Float64, "pos2": pl.Float64,
        "R2": pl.Float64, "Dprime": pl.Float64,
        "ALT_AF1": pl.Float64, "ALT_AF2": pl.Float64,
        "+/-corr": pl.Utf8,
        "rsID1": pl.Utf8, "rsID2": pl.Utf8,
        "REF1": pl.Utf8, "ALT1": pl.Utf8,
        "REF2": pl.Utf8, "ALT2": pl.Utf8
    }

    # =====================================================================
    # Discover chromosome files
    # =====================================================================
    folder_path = Path(folder_path)
    all_files = os.listdir(folder_path)

    imp_files = [f for f in all_files if f.startswith("imputation_results_chr")]
    info_files = [f for f in all_files if f.startswith("LD_info_TOP_LD_chr")]

    if not imp_files or not info_files:
        raise FileNotFoundError("❌ Missing imputation results or LD info files.")

    def _chr_sort_key(x):
        return (0, int(x)) if x.isdigit() else (1, {"X": 23, "Y": 24, "MT": 25}.get(x.upper(), 99))

    chromosomes = sorted(
        set(re.findall(r"chr(\w+)", " ".join(imp_files))),
        key=_chr_sort_key
    )

    #print(f"🧬 Detected chromosomes: {', '.join(chromosomes)}")

    # =====================================================================
    # Per-chromosome worker
    # =====================================================================
    results, summaries = [], []

    def process_chr(chr_id):
        try:
            imp_file = folder_path / f"imputation_results_chr{chr_id}.txt"
            info_file = folder_path / f"LD_info_TOP_LD_chr{chr_id}.txt"

            if not imp_file.exists() or not info_file.exists():
                print(f"            ⚠️ chr{chr_id}: missing files, skipping.")
                return None, None

            # -------- Load files --------
            data_df = coerce_dtypes(
                pl.read_csv(imp_file, separator="\t", infer_schema_length=5000),
                data_types
            )
            info_df = coerce_dtypes(
                pl.read_csv(info_file, separator="\t", infer_schema_length=5000),
                info_types
            )

            # -------------------------------------------------------------
            # Duplicate SNP diagnostics
            # -------------------------------------------------------------
            duplicated = data_df.filter(pl.col("snp").is_duplicated())

            imputed_dups = duplicated.filter(pl.col("imputed") == 1)
            imputed_dups = safe_drop(imputed_dups, ["NC", "SS", "AF", "LP", "SI"])

            nonimp_dups = duplicated.filter(pl.col("imputed") != 1)
            nonimp_dups = safe_drop(nonimp_dups, ["NC", "SS", "AF", "LP", "SI"])

            beta_corr = z_corr = np.nan
            dup_beta = (
                imputed_dups
                .select(["snp", "beta", "z"])
                .rename({"beta": "beta_imp", "z": "z_imp"})
                .join(
                    nonimp_dups
                    .select(["snp", "beta", "z"])
                    .rename({"beta": "beta_nonimp", "z": "z_nonimp"}),
                    on="snp",
                    how="inner"
                )
            )

            if dup_beta.height >= 2:
                beta_corr = np.corrcoef(
                    dup_beta["beta_imp"].to_numpy(),
                    dup_beta["beta_nonimp"].to_numpy()
                )[0, 1]

                z_corr = np.corrcoef(
                    dup_beta["z_imp"].to_numpy(),
                    dup_beta["z_nonimp"].to_numpy()
                )[0, 1]
            else:
                beta_corr = z_corr = np.nan
            # -------------------------------------------------------------
            # NON-IMPUTED summary statistics
            # -------------------------------------------------------------
            nonimp = data_df.filter(pl.col("imputed") == 0)
            nonimp = safe_drop(nonimp, "R2")

            if "LP" in nonimp.columns:
                nonimp = nonimp.with_columns(
                    (10 ** (-pl.col("LP"))).alias("p_value")
                )
                nonimp = safe_drop(nonimp, "LP")

            # Leave NC, SS, AF, SI for later propagation

            # -------------------------------------------------------------
            # IMPUTED SNPs (initial cleaning)
            # -------------------------------------------------------------
            imputed = (
                data_df.filter(pl.col("imputed") == 1)
                .filter(~pl.col("snp").is_in(nonimp["snp"]))
            )
            imputed = safe_drop(imputed, ["R2", "NC", "SS", "AF", "LP", "SI"])

            # =================================================================
            # ⭐ SAMPLE-SIZE PROPAGATION (NC, SS, AF, SI for imputed SNPs)
            # =================================================================
            donor_stats = (
                nonimp
                .select(["snp", "NC", "SS", "AF", "SI"])
                .drop_nulls(subset=["NC", "SS"])
            )

            linked = info_df.join(
                donor_stats,
                left_on="rsID1",
                right_on="snp",
                how="inner"
            )

            imputed_stats = (
                linked.group_by("rsID2")
                .agg([
                    pl.col("NC").mean(),
                    pl.col("SS").mean(),
                    pl.col("AF").mean(),
                    pl.col("SI").mean(),
                ])
                .rename({"rsID2": "snp"})
            )

            imputed = imputed.join(imputed_stats, on="snp", how="left")

            # Compute p-values for imputed SNPs
            if "z" in imputed.columns:
                imputed = imputed.with_columns(
                    (2 * pl.Series(norm.sf(np.abs(imputed["z"].to_numpy()))))
                    .alias("p_value")
                )

            # -------------------------------------------------------------
            # Combine non-imputed + imputed
            # -------------------------------------------------------------
            nonimp = nonimp.with_columns(pl.lit("no").alias("imputed"))
            imputed = imputed.with_columns(pl.lit("yes").alias("imputed"))

            final_df = pl.concat([nonimp, imputed], how="diagonal")

            # NC/SS rounding and ncase
            for col in ["NC", "SS"]:
                if col in final_df.columns:
                    final_df = final_df.with_columns(
                        pl.col(col).round(0).cast(pl.Int64, strict=False)
                    )

            if "SS" in final_df.columns and "NC" in final_df.columns:
                final_df = final_df.with_columns(
                    pl.when(pl.col("SS") > pl.col("NC"))
                    .then((pl.col("SS") - pl.col("NC")).cast(pl.Int64, strict=False))
                    .otherwise(0)
                    .alias("ncase_col")
                )

            summary = {
                "chromosome": chr_id,
                "imputed_markers": imputed.height,
                "n_dups_imputed": imputed_dups.height,
                "n_dups_nonimputed": nonimp_dups.height,
                "beta_corr": float(beta_corr) if beta_corr == beta_corr else None,
                "z_corr": float(z_corr) if z_corr == z_corr else None,
            }

            return final_df, summary

        except Exception as e:
            print(f"            ❌ Error processing chr{chr_id}: {e}")
            return None, None

    # =====================================================================
    # Parallel execution
    # =====================================================================
    with ThreadPoolExecutor(max_workers=threads) as exe:
        futures = {exe.submit(process_chr, c): c for c in chromosomes}

        for fut in as_completed(futures):
            df, summ = fut.result()
            if df is not None:
                results.append(df)
            if summ is not None:
                summaries.append(summ)

    if not results:
        raise RuntimeError("❌ No chromosome processed successfully!")

    combined_df = pl.concat(results, how="diagonal")
    corr_df = pl.DataFrame(summaries)

    # =====================================================================
    # Output setup
    # =====================================================================
    output_folder = output_path or folder_path
    os.makedirs(output_folder, exist_ok=True)
    output_prefix = output_prefix or "PRED_LD"

    imput_harmonisation_folder = Path(output_folder).parent / "imputed_harmonised"
    imput_harmonisation_folder = str(imput_harmonisation_folder)

    combined_path = f"{output_folder}/{output_prefix}_PREDLD_allchr.tsv"
    corr_path = f"{output_folder}/{output_prefix}_PREDLD_correlations.tsv"
    config_path = f"{output_folder}/{output_prefix}_gwas2vcf_config.csv"

    combined_df.write_csv(combined_path, separator="\t")
    corr_df.write_csv(corr_path, separator="\t")
    os.system(f"gzip {combined_path}")

    # =====================================================================
    # Build gwas2vcf config
    # =====================================================================
    columns = [
        "sumstat_file", "gwas_outputname",
        "chr_col", "pos_col", "snp_id_col",
        "ea_col", "oa_col", "eaf_col",
        "beta_or_col", "se_col", "imp_z_col",
        "pval_col", "ncontrol_col", "ncase_col",
        "ncontrol", "ncase", "imp_info_col",
        "delimiter", "infofile", "infocolumn",
        "eaffile", "eafcolumn", "liftover",
        "chr_pos_col", "resourse_folder", "output_folder"
    ]

    cfg = pd.DataFrame([[None] * len(columns)], columns=columns)

    cfg.loc[0, [
        "sumstat_file", "gwas_outputname",
        "output_folder", "resourse_folder", "delimiter"
    ]] = [
        f"{combined_path}.gz",
        f"{output_prefix}_imputed",
        imput_harmonisation_folder,
        gwas2vcf_resource_folder,
        "\t"
    ]

    mapping = {
        "chr_col": "chr",
        "pos_col": "pos",
        "snp_id_col": "snp",
        "ea_col": "A1",
        "oa_col": "A2",
        "eaf_col": "AF",
        "beta_or_col": "beta",
        "se_col": "SE",
        "pval_col": "p_value",
        "ncontrol_col": "NC",
        "ncase_col": "ncase_col",
        "imp_info_col": "SI",
    }

    for cfg_col, df_col in mapping.items():
        cfg.loc[0, cfg_col] = df_col if df_col in combined_df.columns else "NA"

    cfg.to_csv(config_path, index=False)

    # =====================================================================
    # Cleanup original files
    # =====================================================================
    subprocess.run(
        f"(pigz -p {threads} -c {folder_path}/imputation_results_chr*.txt 2>/dev/null "
        f"|| gzip -1 -c {folder_path}/imputation_results_chr*.txt) "
        f"> {output_folder}/{output_prefix}_imputation_results.txt.gz",
        shell=True,
        check=False,
    )

    subprocess.run(
        f"(pigz -p {threads} -c {folder_path}/LD_info_TOP_LD_chr*.txt 2>/dev/null "
        f"|| gzip -1 -c {folder_path}/LD_info_TOP_LD_chr*.txt) "
        f"> {output_folder}/{output_prefix}_LD_info_TOP_LD.txt.gz",
        shell=True,
        check=False,
    )
    subprocess.run(f"rm -f {folder_path}/imputation_results_chr*.txt", shell=True, check=False)
    subprocess.run(f"rm -f {folder_path}/LD_info_TOP_LD_chr*.txt", shell=True, check=False)
    subprocess.run(f"rm -f {output_folder}/{output_prefix}_chr*_predld.log", shell=True, check=False)
    return combined_df, corr_df,config_path
