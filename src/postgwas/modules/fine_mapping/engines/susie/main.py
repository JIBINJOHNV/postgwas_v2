import csv
import shlex
import subprocess
from pathlib import Path

import pandas as pd

from postgwas.core.r_runtime import resolve_r_runtime


SUSIE_REQUIRED_R_PACKAGES = (
    "argparse",
    "glue",
    "data.table",
    "susieR",
    "Matrix",
    "R.utils",
    "ggplot2",
    "ggrepel",
    "processx",
    "jsonlite",
)


def resolve_susie_r_runtime(rscript, timeout_seconds):
    """Resolve one R runtime and validate its complete SuSiE package stack."""
    return resolve_r_runtime(
        rscript, timeout_seconds, SUSIE_REQUIRED_R_PACKAGES, label="SuSiE",
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
    sumstat_manifest,
    sample_id,
    ld_ref,
    plink,
    output_folder,
    resolved_configuration_file,
    rscript,
    r_environment,
):

    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    log_file = output_folder / f"{sample_id}_run_susie.log"

    script_dir = Path(__file__).resolve().parent
    rscript_file = script_dir / "run_susie_parallel_cli.R"
    resolved_configuration_file = Path(resolved_configuration_file).resolve()

    if not rscript_file.exists():
        raise FileNotFoundError(f"❌ ERROR: Cannot find R script at {rscript_file}")
    if not resolved_configuration_file.is_file():
        raise FileNotFoundError(
            "Resolved SuSiE configuration JSON is missing: "
            f"{resolved_configuration_file}"
        )

    cmd = [
        str(rscript), str(rscript_file),
        "--locus_file", str(locus_file),
        "--sumstat_manifest", str(sumstat_manifest),
        "--sample_id", str(sample_id),
        "--ld_ref", str(ld_ref),
        "--plink", str(plink),
        "--SUSIE_Analysis_folder", str(output_folder),
        "--resolved_configuration_file", str(resolved_configuration_file),
    ]

    with open(log_file, "w") as log:
        log.write("CMD: " + shlex.join(cmd) + "\n")
        log.write("--- R OUTPUT ---\n")
        log.flush()
        try:
            subprocess.run(
                cmd,
                stdout=log,
                stderr=log,
                check=True,
                env=r_environment,
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
