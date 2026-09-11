#!/usr/bin/env python3

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
import gzip
import os
from pathlib import Path
import re
import shutil
import subprocess
import threading
import time
from typing import List, Optional, Sequence, Tuple

import pandas as pd
import numpy as np
import polars as pl
import psutil
from scipy.stats import norm, spearmanr

from postgwas.core.execution.runtime import safe_thread_count
from postgwas.core.reference_resources import require_file_inventory
from postgwas.core.validation_reporting import register_file_availability_bundle

"""Resource-bounded chromosome workers and post-processing for PRED-LD.

The configured schedule interleaves chromosomes to reduce peak RAM. Workers
read the validated reference in place, preserve complete per-chromosome logs,
and require every configured chromosome output before consolidation.
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

def _normalise_chromosome(chr_id: str) -> str:
    return str(chr_id).upper().removeprefix("CHR")


def is_big_chr(chr_id: str, large_chromosomes: Sequence[str]) -> bool:
    """Return True if chromosome is considered 'big' for scheduling."""
    return _normalise_chromosome(chr_id) in {
        _normalise_chromosome(chromosome)
        for chromosome in large_chromosomes
    }


def build_interleaved_chr_order(
    chromosomes: Sequence[str],
    preferred_order: Sequence[str],
) -> List[str]:
    """Apply the configured RAM-balanced order without dropping chromosomes."""
    chroms = list(dict.fromkeys(
        _normalise_chromosome(chromosome) for chromosome in chromosomes
    ))
    available = set(chroms)
    ordered = [
        _normalise_chromosome(chromosome)
        for chromosome in preferred_order
        if _normalise_chromosome(chromosome) in available
    ]
    ordered_set = set(ordered)
    return ordered + [
        chromosome for chromosome in chroms if chromosome not in ordered_set
    ]


# =============================================================================
#  MAIN PRED-LD PARALLEL RUNNER
# ============================================================================
def pred_ld_script_path() -> Path:
    """Return the bundled PRED-LD entry point or fail before analysis output."""
    script = Path(__file__).resolve().parent / "pred_ld.py"
    if not script.is_file() or script.stat().st_size <= 0:
        raise FileNotFoundError("PRED-LD entry point is missing or empty: %s" % script)
    return script


def validate_pred_ld_reference(
    reference_directory: str | Path,
    *,
    population: str,
    chromosomes: Sequence[str],
    mode: str,
    subdirectory_template: str,
    file_template: str,
    file_kinds: Sequence[str],
) -> tuple[Path, tuple[Path, ...]]:
    """Validate the exact files consumed by the bundled supported reference mode."""
    root = Path(reference_directory).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError("PRED-LD reference directory does not exist: %s" % root)
    if mode != "TOP_LD":
        raise ValueError(
            "PostGWAS post-processing currently supports only PRED-LD reference "
            "mode TOP_LD; resolved mode was %r." % mode
        )
    reference_root = root / subdirectory_template.format(
        mode=mode,
        population=population,
    )
    files = tuple(
        reference_root / file_template.format(
            population=population,
            chromosome=chromosome,
            kind=kind,
        )
        for chromosome in chromosomes
        for kind in file_kinds
    )
    require_file_inventory(
        files, "PRED-LD TOP_LD reference",
        missing_message="PRED-LD TOP_LD reference files are missing or empty",
        error_type=FileNotFoundError,
    )
    register_file_availability_bundle(
        files,
        "PRED-LD chromosome reference bundle",
        (
            ("count", "chromosomes", list(chromosomes)),
            (
                "success",
                "pred_ld_reference_files",
                "%d / %d" % (len(files), len(files)),
            ),
            ("count", "pred_ld_reference_file_kinds", list(file_kinds)),
            ("analysis", "pred_ld_reference_mode", mode),
            ("genetic", "population", population),
            ("info", "pred_ld_reference_directory", reference_root, True),
        ),
    )
    return root, files


def run_pred_ld_parallel(
    predld_input_dir: str,
    output_folder: str,
    output_prefix: str,
    pred_ld_ref: str,
    chromosomes: Sequence[str],
    r2threshold: float,
    maf: float,
    population: str,
    ref: str,
    threads: int,
    memory_gb: float,
    pred_ld_script: str | Path,
    python_executable: str,
    memory_gb_per_worker: float,
    free_memory_threshold_gb: float,
    free_memory_threshold_fraction: float,
    memory_poll_seconds: float,
    worker_poll_seconds: float,
    large_chromosomes: Sequence[str],
    preferred_chromosome_order: Sequence[str],
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

    max_workers = safe_thread_count(
        threads,
        gb_per_thread=memory_gb_per_worker,
        available_ram_gb=memory_gb,
        enforce_memory_budget=True,
    )
    pred_ld_script = Path(pred_ld_script).expanduser().resolve()
    if not pred_ld_script.is_file() or pred_ld_script.stat().st_size <= 0:
        raise FileNotFoundError("PRED-LD entry point is missing or empty: %s" % pred_ld_script)

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
    # 2. Threads / workers & chromosome order
    # -------------------------------------------------------------------------
    chr_order = build_interleaved_chr_order(
        chromosomes, preferred_chromosome_order,
    )

    write_log(master_log_file, f"Running PRED-LD with {max_workers} parallel workers")
    write_log(
        master_log_file,
        f"Resolved memory budget: {memory_gb} GB; "
        f"per-worker reservation: {memory_gb_per_worker} GB",
    )
    write_log(master_log_file, f"Input dir : {predld_input_dir}")
    write_log(master_log_file, f"Output dir: {output_folder}")
    write_log(master_log_file, f"Ref LD    : {pred_ld_ref}")
    write_log(master_log_file, f"Chromosomes: {', '.join(chr_order)}")

    # -------------------------------------------------------------------------
    # 3. Shared state for scheduling
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
                    threshold = min(
                        free_memory_threshold_gb,
                        total_gb * free_memory_threshold_fraction,
                    )
                    free_ok = free_gb >= threshold
                else:
                    threshold = None
                    free_ok = True

                this_big = is_big_chr(chr_id, large_chromosomes)
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

            time.sleep(memory_poll_seconds)

    def release_chr(chr_id: str) -> None:
        nonlocal running_big
        with lock:
            running_chroms.discard(chr_id)
            if is_big_chr(chr_id, large_chromosomes):
                if not any(
                    is_big_chr(chromosome, large_chromosomes)
                    for chromosome in running_chroms
                ):
                    running_big = False

    # -------------------------------------------------------------------------
    # 4. Worker (SILENT MODE)
    # -------------------------------------------------------------------------
    def worker(chr_id: str) -> None:
        chr_input = predld_input_dir / f"{output_prefix}_chr{chr_id}_pred_ld_input.tsv"
        log_file = output_folder / f"{output_prefix}_chr{chr_id}_predld.log"

        if log_file.exists():
            log_file.unlink()
        write_log(log_file, f"--- Processing Chromosome {chr_id} ---")

        if not chr_input.is_file() or chr_input.stat().st_size <= 0:
            msg = f"PRED-LD input missing or empty for chr{chr_id}: {chr_input}"
            write_log(log_file, f"⚠️ SKIPPED: {msg}")
            # REMOVED: print(f"⚠️ chr{chr_id} skipped...")
            failed.append((chr_id, msg))
            return

        wait_until_allowed(chr_id, log_file)

        try:
            cmd = [
                str(python_executable), str(pred_ld_script),
                "--file-path", str(chr_input),
                "--pop", population,
                "--ref", ref,
                "--ref_dir", str(pred_ld_ref),
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

            expected_outputs = (
                output_folder / f"imputation_results_chr{chr_id}.txt",
                output_folder / f"LD_info_TOP_LD_chr{chr_id}.txt",
            )
            missing_outputs = [
                path for path in expected_outputs
                if not path.is_file() or path.stat().st_size <= 0
            ]
            if missing_outputs:
                reason = (
                    "PRED-LD output missing or empty for chr%s: %s. See log: %s"
                    % (
                        chr_id,
                        ", ".join(str(path) for path in missing_outputs),
                        log_file,
                    )
                )
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
    # 5. Launch workers
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
            time.sleep(worker_poll_seconds)

    for t in threads_list:
        t.join()

    elapsed_min = (time.time() - t0) / 60.0
    merged_log = output_folder / f"{output_prefix}_predld_combined.log"
    failed_ids = sorted(
        (str(chromosome) for chromosome, _reason in failed),
        key=lambda value: (
            (0, int(value)) if value.isdigit() else (1, value)
        ),
    )
    write_log(
        master_log_file,
        "PRED-LD chromosome execution finished in %.2f minutes: "
        "successful=%d failed=%d"
        % (elapsed_min, len(succeeded), len(failed)),
    )
    if failed:
        write_log(master_log_file, "Failed or missing chromosomes:")
        for chr_id, reason in sorted(
            failed,
            key=lambda item: (
                (0, int(str(item[0])))
                if str(item[0]).isdigit()
                else (1, str(item[0]))
            ),
        ):
            write_log(master_log_file, f"   chr{chr_id}: {reason}")
    write_log(master_log_file, f"Merging chromosome logs into {merged_log}")

    # -------------------------------------------------------------------------
    # 6. Merge Logs
    # -------------------------------------------------------------------------
    with merged_log.open("w") as fout:
        if master_log_file.exists():
            with master_log_file.open() as f:
                fout.write(f.read())
                fout.write("\n" + "="*80 + "\n")

        for lf in sorted(output_folder.glob(f"{output_prefix}_chr*_predld.log")):
            with lf.open() as f:
                fout.write(f"\nLOG: {lf.name}\n")
                fout.write(f.read())

    if failed:
        print(f"\n          ❌ PRED-LD Failed for chromosomes: {', '.join(failed_ids)}")
        print(f"            See detailed logs in: {merged_log}\n")
        raise RuntimeError(
            "PRED-LD did not complete every configured chromosome; failed or "
            "missing chromosomes: %s. See %s"
            % (", ".join(failed_ids), merged_log)
        )
    else:
        print(
            "\n          ✅ PRED-LD completed all %d configured chromosomes.\n"
            % len(succeeded)
        )

    if not succeeded:
        raise RuntimeError(
            "PRED-LD produced no successful chromosome result. See %s"
            % merged_log
        )
    return True


def process_pred_ld_results_all_parallel(
    folder_path: str,
    output_path: str,
    output_prefix: str,
    harmonised_dataset_id: str,
    corr_method: str,
    threads: int,
):
    """
    Production-ready PRED-LD postprocessing pipeline.
    - Parallel Polars implementation
    - Safe dtype coercion for all PRED-LD fields
    - Accurate imputed-SNP sample-size propagation
    - Duplicate SNP diagnostics + correlation
    - Final TSV + correlation summary
    - Harmonisation sample-sheet generation
    - Cleanup of intermediate files
    """
    if corr_method not in {"pearson", "spearman"}:
        raise ValueError(
            "Unsupported PRED-LD correlation method: %r" % corr_method
        )

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

    def paired_correlation(first, second) -> float:
        """Apply the resolved QC correlation method to aligned observations."""
        if corr_method == "pearson":
            return float(np.corrcoef(first, second)[0, 1])
        if corr_method == "spearman":
            return float(spearmanr(first, second).statistic)
        raise AssertionError("validated PRED-LD correlation method was lost")

    def write_gzip(source: Path, destination: Path) -> None:
        """Compress one file without invoking a shell or an external compressor."""
        with source.open("rb") as input_handle, gzip.open(
            destination, "wb", compresslevel=1,
        ) as output_handle:
            shutil.copyfileobj(input_handle, output_handle)

    def archive_text_files(sources: Sequence[Path], destination: Path) -> None:
        """Concatenate chromosome text files into one deterministic gzip stream."""
        with gzip.open(destination, "wb", compresslevel=1) as output_handle:
            for source in sources:
                with source.open("rb") as input_handle:
                    shutil.copyfileobj(input_handle, output_handle)

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
    postprocessing_failures: list[tuple[str, str]] = []

    def process_chr(chr_id):
        try:
            imp_file = folder_path / f"imputation_results_chr{chr_id}.txt"
            info_file = folder_path / f"LD_info_TOP_LD_chr{chr_id}.txt"

            if any(
                not path.is_file() or path.stat().st_size <= 0
                for path in (imp_file, info_file)
            ):
                return None, None, (
                    "required chromosome result pair is incomplete: %s, %s"
                    % (imp_file, info_file)
                )

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
                beta_corr = paired_correlation(
                    dup_beta["beta_imp"].to_numpy(),
                    dup_beta["beta_nonimp"].to_numpy()
                )

                z_corr = paired_correlation(
                    dup_beta["z_imp"].to_numpy(),
                    dup_beta["z_nonimp"].to_numpy()
                )
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

            imputed = imputed.join(
                imputed_stats, on="snp", how="left", coalesce=True,
            )

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

            return final_df, summary, None

        except Exception as e:
            return None, None, "%s: %s" % (type(e).__name__, e)

    # =====================================================================
    # Parallel execution
    # =====================================================================
    with ThreadPoolExecutor(max_workers=threads) as exe:
        futures = {exe.submit(process_chr, c): c for c in chromosomes}

        for fut in as_completed(futures):
            chromosome = futures[fut]
            df, summ, failure = fut.result()
            if df is not None:
                results.append(df)
            if summ is not None:
                summaries.append(summ)
            if failure is not None:
                postprocessing_failures.append((chromosome, failure))

    if postprocessing_failures:
        postprocessing_failures.sort(
            key=lambda item: _chr_sort_key(str(item[0])),
        )
        combined_log = folder_path / f"{output_prefix}_predld_combined.log"
        with combined_log.open("a", encoding="utf-8") as handle:
            handle.write("\nPRED-LD post-processing failed:\n")
            for chromosome, reason in postprocessing_failures:
                handle.write("  chr%s: %s\n" % (chromosome, reason))
        raise RuntimeError(
            "PRED-LD post-processing failed for chromosome(s) %s. See %s"
            % (
                ", ".join(
                    chromosome for chromosome, _reason in postprocessing_failures
                ),
                combined_log,
            )
        )

    if not results:
        raise RuntimeError("❌ No chromosome processed successfully!")

    combined_df = pl.concat(results, how="diagonal")
    corr_df = pl.DataFrame(summaries)

    # =====================================================================
    # Output setup
    # =====================================================================
    output_folder = Path(output_path).expanduser().resolve()
    output_folder.mkdir(parents=True, exist_ok=True)

    combined_path = output_folder / f"{output_prefix}_PREDLD_allchr.tsv"
    corr_path = output_folder / f"{output_prefix}_PREDLD_correlations.tsv"
    sample_sheet_path = output_folder / f"{output_prefix}_harmonisation_sample_sheet.csv"

    combined_df.write_csv(combined_path, separator="\t")
    corr_df.write_csv(corr_path, separator="\t")
    combined_gzip = Path(str(combined_path) + ".gz")
    write_gzip(combined_path, combined_gzip)
    combined_path.unlink()

    # =====================================================================
    # Build the canonical harmonisation sample sheet
    # =====================================================================
    required = {"chr", "pos", "A1", "A2", "AF", "beta", "SE", "p_value", "SI"}
    missing = sorted(required - set(combined_df.columns))
    if missing:
        raise RuntimeError(
            "PRED-LD output cannot be harmonised because required columns are missing: "
            + ", ".join(missing)
        )
    row = {
        "config_version": 2,
        "dataset_id": harmonised_dataset_id,
        "input_file": str(combined_gzip),
        "chromosome_column": "chr",
        "position_column": "pos",
        "variant_id_column": "snp" if "snp" in combined_df.columns else None,
        "effect_allele_column": "A1",
        "other_allele_column": "A2",
        "effect_allele_frequency_column": "AF",
        "effect_column": "beta",
        "standard_error_column": "SE",
        "p_value_column": "p_value",
        "imputation_info_column": "SI",
        "control_count_column": "NC" if "NC" in combined_df.columns else None,
        "case_count_column": "ncase_col" if "ncase_col" in combined_df.columns else None,
    }
    pd.DataFrame([row]).to_csv(sample_sheet_path, index=False)

    # =====================================================================
    # Cleanup original files
    # =====================================================================
    imputation_sources = sorted(folder_path.glob("imputation_results_chr*.txt"))
    information_sources = sorted(folder_path.glob("LD_info_TOP_LD_chr*.txt"))
    archive_text_files(
        imputation_sources,
        output_folder / f"{output_prefix}_imputation_results.txt.gz",
    )
    archive_text_files(
        information_sources,
        output_folder / f"{output_prefix}_LD_info_TOP_LD.txt.gz",
    )
    for source in (*imputation_sources, *information_sources):
        source.unlink()
    for log_file in output_folder.glob(f"{output_prefix}_chr*_predld.log"):
        log_file.unlink()
    return combined_df, corr_df, str(sample_sheet_path)


__all__ = [
    "pred_ld_script_path",
    "process_pred_ld_results_all_parallel",
    "run_pred_ld_parallel",
    "validate_pred_ld_reference",
]
