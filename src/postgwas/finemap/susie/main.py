import subprocess
import shlex
from pathlib import Path
import pandas as pd
import csv

from postgwas.finemap.defaults import (
    DEFAULT_GENOME_BUILD,
    DEFAULT_MHC_CHROMOSOME,
    DEFAULT_MHC_END,
    DEFAULT_MHC_START,
    DEFAULT_SUSIE_LD_TIMEOUT_SECONDS,
    DEFAULT_SUSIE_LP_THRESHOLD,
    DEFAULT_SUSIE_MAX_CAUSAL_COMPONENTS,
    DEFAULT_SUSIE_MIN_RAM_PER_WORKER_GB,
    DEFAULT_SUSIE_TIMEOUT_SECONDS,
)

def validate_locus_file(path, locus_type="range"):
    # --- Step 1: read first few lines to guess delimiter ---
    try:
        with open(path, "r") as f:
            sample = f.read(2048)
    except Exception as e:
        raise ValueError(f"Cannot open locus file '{path}': {e}") from e

    # Try Sniffer first
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=[",", "\t", ";", " "])
        sep = dialect.delimiter
    except Exception:
        if "\t" in sample:
            sep = "\t"
        elif "," in sample:
            sep = ","
        elif ";" in sample:
            sep = ";"
        else:
            sep = r"\s+"

    # --- Step 2: read using inferred delimiter ---
    try:
        df = pd.read_csv(path, sep=sep, engine="python")
    except Exception as e:
        raise ValueError(
            f"Cannot parse locus file '{path}' with inferred separator '{sep}': {e}"
        ) from e

    # --- Step 3: normalize column names ---
    df.columns = [c.strip().upper() for c in df.columns]

    # --- Step 4: validate required columns ---
    locus_type = str(locus_type).strip().lower()
    required = {"CHR", "POS"} if locus_type == "point" else {"CHR", "START", "END"}
    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            "Locus file is missing required columns.\n"
            f"   Missing: {', '.join(missing)}\n"
            f"   Required: {', '.join(sorted(required))}\n"
            f"   Found: {', '.join(df.columns)}\n"
            f"   Inferred delimiter: '{sep}'"
        )

    #print(f"✅ Locus file loaded successfully using delimiter '{sep}'.")
    if df.empty:
        raise ValueError(f"Locus file '{path}' contains no loci")

    return df

def run_susie(
    locus_file,
    sumstat_file,
    sample_id,
    ld_ref,
    plink,
    output_folder,
    lp_threshold=DEFAULT_SUSIE_LP_THRESHOLD,
    L=DEFAULT_SUSIE_MAX_CAUSAL_COMPONENTS,
    workers="auto",
    min_ram_per_worker_gb=DEFAULT_SUSIE_MIN_RAM_PER_WORKER_GB,
    timeout_ld_seconds=DEFAULT_SUSIE_LD_TIMEOUT_SECONDS,
    timeout_susie_seconds=DEFAULT_SUSIE_TIMEOUT_SECONDS,
    skip_mhc=False,
    finemap_mhc_chrom=DEFAULT_MHC_CHROMOSOME,
    mhc_start=DEFAULT_MHC_START,
    mhc_end=DEFAULT_MHC_END,
    verbose=False,
    genome_build=DEFAULT_GENOME_BUILD,
):

    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    log_file = output_folder / f"{sample_id}_run_susie.log"

    script_dir = Path(__file__).resolve().parent
    rscript_file = script_dir / "run_susie_parallel_cli.R"

    if not rscript_file.exists():
        raise FileNotFoundError(f"❌ ERROR: Cannot find R script at {rscript_file}")

    cmd = [
        "Rscript", str(rscript_file),
        "--locus_file", str(locus_file),
        "--sumstat_file", str(sumstat_file),
        "--sample_id", str(sample_id),
        "--ld_ref", str(ld_ref),
        "--plink", str(plink),
        "--SUSIE_Analysis_folder", str(output_folder),
        "--lp_threshold", str(lp_threshold),
        "--L", str(L),
        "--workers", str(workers),
        "--min_ram_per_worker_gb", str(min_ram_per_worker_gb),
        "--timeout_ld_seconds", str(timeout_ld_seconds),
        "--timeout_susie_seconds", str(timeout_susie_seconds),
        "--finemap_mhc_chrom", str(finemap_mhc_chrom),
        "--mhc_start", str(mhc_start),   # ✅ FIX
        "--mhc_end", str(mhc_end),       # ✅ FIX
        "--genome_build", str(genome_build),
    ]

    if skip_mhc:
        cmd.append("--skip_mhc")

    if verbose:
        cmd.append("--verbose")

    with open(log_file, "w") as log:
        log.write("CMD: " + shlex.join(cmd) + "\n")
        log.write("--- R OUTPUT ---\n")
        log.flush()
        try:
            subprocess.run(
                cmd,
                stdout=log,
                stderr=log,
                check=True
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"SuSiE failed (exit code {e.returncode}); see {log_file}"
            ) from e

    flames_input = output_folder / "flames_input"

    return {
        "status": "success",
        "log_file": str(log_file),
        "output_dir": str(output_folder),
        "flames_input": str(flames_input),
    }
