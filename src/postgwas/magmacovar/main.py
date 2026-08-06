import subprocess
import pandas as pd
from pathlib import Path
import datetime
import platform


# ============================================================
# Helper: Compare first-column gene lists (no column-name checks)
# ============================================================
def run_magma_covariates(
    magma_bin: str,
    gene_results_file: str,
    covariates_file: str,
    output_prefix: str,
    model: str = None,
    direction: str = None,
    log_file: str = None,
):
    """
    Run MAGMA gene-property (covariate) analysis with correct model/direction logic.

    Rules:
      - If model is None --> no --model flag
      - If direction is None --> no direction= token
      - direction must always follow --model <value>
      - MAGMA default direction is two-sided (do not pass explicitly)
    """

    gene_results_file = Path(gene_results_file)
    covariates_file = Path(covariates_file)
    output_prefix = Path(output_prefix)

    out_dir = output_prefix.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    if log_file is None:
        log_file = output_prefix.with_suffix(".log")
    else:
        log_file = Path(log_file)

    def write_log(msg: str):
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open(log_file, "a") as lf:
            lf.write(f"[{timestamp}] {msg}\n")

    write_log("===== MAGMA Covariates Analysis Started =====")
    write_log(f"MAGMA binary: {magma_bin}")
    write_log(f"Output prefix: {output_prefix}")

    # ---------------------------------------------------------
    # Check required files
    # ---------------------------------------------------------
    if not gene_results_file.exists():
        raise FileNotFoundError(f"Gene results file not found: {gene_results_file}")
    if not covariates_file.exists():
        raise FileNotFoundError(f"Covariates file not found: {covariates_file}")

    # ---------------------------------------------------------
    # Build MAGMA command
    # ---------------------------------------------------------
    cmd = [
        magma_bin,
        "--gene-results", str(gene_results_file),
        "--gene-covar", str(covariates_file),
        "--out", str(output_prefix)
    ]

    # -----------------------
    # MODEL + DIRECTION LOGIC
    # -----------------------

    # Case A — model is provided
    if model:
        write_log(f"Using model: {model}")
        cmd += ["--model", model]

        # direction only valid if model provided
        if direction:
            write_log(f"Using direction: {direction}")
            cmd.append(f"direction-covar={direction}")

    # Case B — model None, direction provided (invalid in MAGMA)
    elif direction:
        write_log(f"⚠ WARNING: direction='{direction}' ignored because no model was provided.")
        # Do NOT add direction without model

    # Case C — both omitted → use default
    else:
        write_log("Using default MAGMA model (no --model passed).")

    # ---------------------------------------------------------
    # Log final command
    # ---------------------------------------------------------
    final_cmd_str = " ".join(cmd)
    write_log("Final MAGMA command:")
    write_log(final_cmd_str)
    #print("🔧 MAGMA CMD:", final_cmd_str)

    # ---------------------------------------------------------
    # Execute MAGMA
    # ---------------------------------------------------------
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )

        if result.stdout:
            write_log("MAGMA STDOUT:\n" + result.stdout)
        if result.stderr:
            write_log("MAGMA STDERR:\n" + result.stderr)

        write_log("MAGMA execution completed successfully.")

    except subprocess.CalledProcessError as e:
        write_log("MAGMA FAILED")
        write_log("STDOUT:\n" + (e.stdout or ""))
        write_log("STDERR:\n" + (e.stderr or ""))
        raise RuntimeError("MAGMA execution failed. Check log for details.") from e

    write_log("===== MAGMA Covariates Analysis Completed =====")

