import os
import sys
import subprocess
from pathlib import Path
import datetime


# =========================================================
# COLORS (Preserved Exact Class)
# =========================================================
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'


# =========================================================
# UTILITY PRINT HELPERS (Preserved)
# =========================================================
def print_step(msg):
    print(f"\n      {Colors.HEADER}=================================================={Colors.ENDC}")
    print(f"        {Colors.BOLD}STEP: {msg}{Colors.ENDC}")
    print(f"        {Colors.HEADER}=================================================={Colors.ENDC}")


def print_error(step_name, specific_msg):
    print(f"\n      {Colors.FAIL}❌ CRITICAL ERROR IN: {step_name}{Colors.ENDC}")
    print(f"        {Colors.WARNING}Details: {specific_msg}{Colors.ENDC}")
    print(f"        Check the log file for full traceback.\n")


# =========================================================
# PATH NORMALIZATION
# =========================================================
def _d(p: Path) -> str:
    return str(p)


# =========================================================
# SUBPROCESS RUNNER (Replaces Docker Runner)
# =========================================================
def run_subprocess(command_list, log_file, step_name):
    """
    Execute shell command inside the container and stream output ONLY into the log file.
    Replaces the old 'run_docker' function but keeps exact logging logic.
    """
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "a") as lf:
        lf.write(f"[{ts}] STARTING: {step_name}\n")
        lf.write(f"[{ts}] COMMAND: {' '.join(command_list)}\n")

    print(f"        {Colors.CYAN}Running {step_name}...{Colors.ENDC}")

    try:
        # We use subprocess.run with stdout redirected to the log file handle.
        # This achieves the same "silence on console, verbose in log" behavior.
        with open(log_file, "a") as lf:
            proc = subprocess.run(
                command_list,
                stdout=lf,              # Redirect stdout to log file
                stderr=subprocess.STDOUT, # Redirect stderr to stdout (so it also goes to log)
                text=True,
                check=False             # We check return code manually below to raise custom error
            )

        if proc.returncode != 0:
            raise subprocess.CalledProcessError(proc.returncode, command_list)

        print(f"            {Colors.GREEN}✔ {step_name} Completed Successfully{Colors.ENDC}")

    except subprocess.CalledProcessError as e:
        print_error(step_name, f"LDSC process exited with code {e.returncode}")
        print(f"   See log for details: {log_file}")
        raise e


# =========================================================
# HELPER: LOG PARSER (Preserved)
# =========================================================
def extract_ldsc_metrics(log_file):
    """
    Parses an LDSC log file to extract h2, Intercept, and Ratio.
    Returns a formatted string or 'Not found'.
    """
    metrics = {
        "h2": "N/A",
        "intercept": "N/A",
        "ratio": "N/A"
    }
    
    path = Path(log_file)
    if not path.exists():
        return metrics

    try:
        with path.open("r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("Total Liability scale h2:"):
                    metrics["h2"] = line.split(":", 1)[1].strip()
                elif line.startswith("Total Observed scale h2:"):
                    metrics["h2"] = line.split(":", 1)[1].strip()
                elif line.startswith("Intercept:"):
                    metrics["intercept"] = line.split(":", 1)[1].strip()
                elif line.startswith("Ratio:"):
                    metrics["ratio"] = line.split(":", 1)[1].strip()
    except Exception:
        pass
        
    return metrics


# =========================================================
# MAIN EXECUTION WRAPPER
# =========================================================
def run_ldsc(
    sumstats_tsv: str,
    out_prefix: str,
    hm3_snplist: str,
    ldscore_dir: str,
    info_min: float = 0.7,
    maf_min: float = 0.01,
    samp_prev: float = None,
    pop_prev: float = None,
):
    """
    Full LDSC workflow (Embedded Mode):
      1) munge_sumstats.py  → prefix.sumstats.gz
      2) ldsc.py --h2 (Liability scale) -> OPTIONAL
      3) ldsc.py --h2 (Observed scale) -> ALWAYS RUN

    All steps run inside the 'ldsc' conda environment within this container.
    """

    # -----------------------------------------------------
    # Resolve & validate files
    # -----------------------------------------------------
    try:
        sumstats_tsv = Path(sumstats_tsv).expanduser().resolve(strict=True)
        hm3_snplist = Path(hm3_snplist).expanduser().resolve(strict=True)
        ldscore_dir = Path(ldscore_dir).expanduser().resolve(strict=True)
        out_prefix = Path(out_prefix).expanduser().resolve()
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        print_error("PRE-FLIGHT CHECKS", f"File not found or invalid: {e}")
        sys.exit(1)

    # -----------------------------------------------------
    # Log file
    # -----------------------------------------------------
    log_file = out_prefix.with_suffix(".ldsc.log")
    with open(log_file, "w") as lf:
        lf.write("=== LDSC LOG START (Embedded Mode) ===\n")

    # -----------------------------------------------------
    # Mount logic -> REMOVED
    # (No longer needed because we are already inside the container)
    # -----------------------------------------------------
    
    # Define absolute path to LDSC scripts inside your Docker image
    # Assumes you cloned LDSC to /opt/ldsc in your Dockerfile
    LDSC_ROOT = "/opt/ldsc"

    # ----------------------------------------------------------------------------
    # STEP 1 — MUNGE SUMSTATS
    # ----------------------------------------------------------------------------
    # print_step("1. Munge Summary Statistics") # Kept commented out as per original

    munged_prefix = out_prefix
    munged_file = f"{munged_prefix}.sumstats.gz"

    # COMMAND CHANGED: Uses 'micromamba run' instead of 'docker run'
    cmd_munge = [
        "micromamba", "run", "-n", "ldsc",
        "python", f"{LDSC_ROOT}/munge_sumstats.py",
        "--sumstats", _d(sumstats_tsv),
        "--out", _d(munged_prefix),
        "--merge-alleles", _d(hm3_snplist),
        "--info-min", str(info_min),
        "--maf-min", str(maf_min)
    ]

    run_subprocess(cmd_munge, log_file, "MUNGE_SUMSTATS")

    if not Path(munged_file).exists():
        print_error("MUNGE_SUMSTATS", f"Output missing: {munged_file}")
        sys.exit(1)

    # ----------------------------------------------------------------------------
    # STEP 2 — LDSC (Liability-scale) [OPTIONAL]
    # ----------------------------------------------------------------------------
    out_liab_log = "Skipped"
    liab_metrics = {"h2": "N/A", "intercept": "N/A", "ratio": "N/A"}

    if samp_prev is not None and pop_prev is not None:
        # print_step("        2. LDSC Heritability (Liability Scale)")

        out_liab = f"{munged_prefix}_Liability_scale_h2"

        # COMMAND CHANGED: Uses 'micromamba run' instead of 'docker run'
        cmd_liability = [
            "micromamba", "run", "-n", "ldsc",
            "python", f"{LDSC_ROOT}/ldsc.py",
            "--h2", _d(munged_file),
            "--ref-ld-chr", _d(ldscore_dir) + "/",
            "--w-ld-chr", _d(ldscore_dir) + "/",
            "--out", str(out_liab),
            "--samp-prev", str(samp_prev),
            "--pop-prev", str(pop_prev),
        ]

        run_subprocess(cmd_liability, log_file, "LDSC_Liability_scale_HERITABILITY")
        out_liab_log = f"{out_liab}.log"
        liab_metrics = extract_ldsc_metrics(out_liab_log)
    else:
        print(f"\n      {Colors.WARNING}⚠️ Skipping Liability Scale Analysis (Missing samp_prev/pop_prev){Colors.ENDC}")

    # ----------------------------------------------------------------------------
    # STEP 3 — LDSC (Observed scale) [ALWAYS RUN]
    # ----------------------------------------------------------------------------
    # print_step("3. LDSC Heritability (Observed Scale)")

    out_obs = f"{munged_prefix}_h2"

    # COMMAND CHANGED: Uses 'micromamba run' instead of 'docker run'
    cmd_obs = [
        "micromamba", "run", "-n", "ldsc",
        "python", f"{LDSC_ROOT}/ldsc.py",
        "--h2", _d(munged_file),
        "--ref-ld-chr", _d(ldscore_dir) + "/",
        "--w-ld-chr", _d(ldscore_dir) + "/",
        "--out", str(out_obs),
    ]

    run_subprocess(cmd_obs, log_file, "LDSC_Observed_scale_HERITABILITY")
    out_obs_log = f"{out_obs}.log"
    obs_metrics = extract_ldsc_metrics(out_obs_log)
    
    # ----------------------------------------------------------------------------
    # COMPLETED (Preserved Exact Formatting)
    # ----------------------------------------------------------------------------
    # Use \t\t to force exactly two tabs of indentation
    print(f"\n\t\t{Colors.GREEN}{Colors.BOLD}🎉 LDSC PIPELINE COMPLETED!{Colors.ENDC}")
    
    # Print Table-like Summary
    print(f"\n\t\t{Colors.BOLD}🔹 HERITABILITY RESULTS:{Colors.ENDC}")
    
    # Indentation inside f-string preserved
    print(f"\t\t{'-' * 60}") 
    
    print(f"\t\t{'METRIC':<25} | {'OBSERVED SCALE':<15} | {'LIABILITY SCALE':<15}")
    print(f"\t\t{'-' * 60}")
    print(f"\t\t{'h2 (Heritability)':<25} | {obs_metrics['h2']:<15} | {liab_metrics['h2']:<15}")
    print(f"\t\t{'Intercept':<25} | {obs_metrics['intercept']:<15} | {liab_metrics['intercept']:<15}")
    print(f"\t\t{'Ratio':<25} | {obs_metrics['ratio']:<15} | {liab_metrics['ratio']:<15}")
    print(f"\t\t{'-' * 60}")

    print(f"\n\t\t{Colors.BOLD}📂 OUTPUT FILES:{Colors.ENDC}")
    # Indent list items preserved
    print(f"\t\t   • Munged:   {munged_file}")
    print(f"\t\t   • Logs:     {log_file}")
    print(f"\t\t   • H2 (Obs): {out_obs_log}")
    if samp_prev is not None:
        print(f"\t\t   • H2 (Liab):{out_liab_log}")
    print()
    return {
        "munged_sumstats": munged_file,
        "h2_liability": out_liab_log,
        "h2_observed": out_obs_log
    }