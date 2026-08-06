import os
import sys
import subprocess
from pathlib import Path
import datetime


# =========================================================
# COLORS
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
# UTILITY PRINT HELPERS
# =========================================================
def print_step(msg):
    print(f"\n      {Colors.HEADER}=================================================={Colors.ENDC}")
    print(f"        {Colors.BOLD}STEP: {msg}{Colors.ENDC}")
    print(f"        {Colors.HEADER}=================================================={Colors.ENDC}")


def print_error(step_name, specific_msg):
    print(f"\n      {Colors.FAIL}❌ CRITICAL ERROR IN: {step_name}{Colors.ENDC}")
    print(f"        {Colors.WARNING}Details: {specific_msg}{Colors.ENDC}")
    print(f"        Check the log file for full LDSC traceback.\n")


# =========================================================
# PATH NORMALIZATION
# =========================================================
def infer_mount_root(*paths):
    """
    Determine a common parent directory to mount into the LDSC container.
    """
    resolved = [str(Path(p).expanduser().resolve()) for p in paths]
    return Path(os.path.commonpath(resolved))


def _d(p: Path) -> str:
    return str(p)


# =========================================================
# DOCKER RUNNER
# =========================================================
# =========================================================
# DOCKER RUNNER
# =========================================================
def run_docker(command_list, log_file, step_name):
    """
    Execute docker command and stream output ONLY into the log file.
    """
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "a") as lf:
        lf.write(f"[{ts}] STARTING: {step_name}\n")
        lf.write(f"[{ts}] COMMAND: {' '.join(command_list)}\n")

    print(f"        {Colors.CYAN}Running {step_name}...{Colors.ENDC}")

    try:
        with subprocess.Popen(
            command_list,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, # Redirect stderr to stdout so we capture errors too
            text=True,
            bufsize=1
        ) as proc, open(log_file, "a") as lf:

            for line in proc.stdout:
                # print(line, end="")  <--- REMOVED THIS LINE to silence console output
                lf.write(line)

            exit_code = proc.wait()

            if exit_code != 0:
                raise subprocess.CalledProcessError(exit_code, command_list)

        print(f"            {Colors.GREEN}✔ {step_name} Completed Successfully{Colors.ENDC}")

    except subprocess.CalledProcessError as e:
        print_error(step_name, f"LDSC Docker process exited with code {e.returncode}")
        # Optionally print the last few lines of the log if it failed
        print(f"   See log for details: {log_file}")
        raise e


# =========================================================
# HELPER: LOG PARSER
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
    docker_image: str = "jibinjv/ldsc:1.0.1",
    platform: str = "linux/amd64",
    info_min: float = 0.7,
    maf_min: float = 0.01,
    samp_prev: float = None,  # Changed default to None
    pop_prev: float = None,   # Changed default to None
):
    """
    Full LDSC workflow:
      1) munge_sumstats.py  → prefix.sumstats.gz
      2) ldsc.py --h2 (Liability scale) -> OPTIONAL (Run only if samp_prev & pop_prev provided)
      3) ldsc.py --h2 (Observed scale) -> ALWAYS RUN

    All steps run inside an external LDSC docker image.
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
        lf.write("=== LDSC LOG START ===\n")

    # -----------------------------------------------------
    # Mount logic
    # -----------------------------------------------------
    mount_root = infer_mount_root(sumstats_tsv, out_prefix, hm3_snplist, ldscore_dir)
    print(f"        [{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Mount root inferred: {mount_root}")

    # Helper for docker paths
    mount_arg = f"{mount_root}:{mount_root}"

    # ----------------------------------------------------------------------------
    # STEP 1 — MUNGE SUMSTATS
    # ----------------------------------------------------------------------------
    #print_step("1. Munge Summary Statistics")

    munged_prefix = out_prefix
    munged_file = f"{munged_prefix}.sumstats.gz"

    cmd_munge = [
        "docker", "run",
        f"--platform={platform}",
        "--rm",
        "-v", mount_arg,
        docker_image,
        "python",
        "munge_sumstats.py",
        "--sumstats", _d(sumstats_tsv),
        "--out", _d(munged_prefix),
        "--merge-alleles", _d(hm3_snplist),
        "--info-min", str(info_min),
        "--maf-min", str(maf_min)
    ]

    run_docker(cmd_munge, log_file, "MUNGE_SUMSTATS")

    if not Path(munged_file).exists():
        print_error("MUNGE_SUMSTATS", f"Output missing: {munged_file}")
        sys.exit(1)

    # ----------------------------------------------------------------------------
    # STEP 2 — LDSC (Liability-scale) [OPTIONAL]
    # ----------------------------------------------------------------------------
    out_liab_log = "Skipped"
    liab_metrics = {"h2": "N/A", "intercept": "N/A", "ratio": "N/A"}

    if samp_prev is not None and pop_prev is not None:
        #print_step("        2. LDSC Heritability (Liability Scale)")

        out_liab = f"{munged_prefix}_Liability_scale_h2"

        cmd_liability = [
            "docker", "run",
            f"--platform={platform}",
            "--rm",
            "-v", mount_arg,
            docker_image,
            "python",
            "ldsc.py",
            "--h2", munged_file,
            "--ref-ld-chr", _d(ldscore_dir) + "/",
            "--w-ld-chr", _d(ldscore_dir) + "/",
            "--out", out_liab,
            "--samp-prev", str(samp_prev),
            "--pop-prev", str(pop_prev),
        ]

        run_docker(cmd_liability, log_file, "LDSC_Liability_scale_HERITABILITY")
        out_liab_log = f"{out_liab}.log"
        liab_metrics = extract_ldsc_metrics(out_liab_log)
    else:
        print(f"\n      {Colors.WARNING}⚠️ Skipping Liability Scale Analysis (Missing samp_prev/pop_prev){Colors.ENDC}")

    # ----------------------------------------------------------------------------
    # STEP 3 — LDSC (Observed scale) [ALWAYS RUN]
    # ----------------------------------------------------------------------------
    #print_step("3. LDSC Heritability (Observed Scale)")

    out_obs = f"{munged_prefix}_h2"

    cmd_obs = [
        "docker", "run",
        f"--platform={platform}",
        "--rm",
        "-v", mount_arg,
        docker_image,
        "python",
        "ldsc.py",
        "--h2", munged_file,
        "--ref-ld-chr", _d(ldscore_dir) + "/",
        "--w-ld-chr", _d(ldscore_dir) + "/",
        "--out", out_obs,
    ]

    run_docker(cmd_obs, log_file, "LDSC_Observed_scale_HERITABILITY")
    out_obs_log = f"{out_obs}.log"
    obs_metrics = extract_ldsc_metrics(out_obs_log)
    # ----------------------------------------------------------------------------
    # COMPLETED
    # ----------------------------------------------------------------------------
    # Use \t\t to force exactly two tabs of indentation
    print(f"\n\t\t{Colors.GREEN}{Colors.BOLD}🎉 LDSC PIPELINE COMPLETED!{Colors.ENDC}")
    
    # Print Table-like Summary
    print(f"\n\t\t{Colors.BOLD}🔹 HERITABILITY RESULTS:{Colors.ENDC}")
    
    # Fix: Put indentation INSIDE the f-string so the dashes are indented
    print(f"\t\t{'-' * 60}") 
    
    print(f"\t\t{'METRIC':<25} | {'OBSERVED SCALE':<15} | {'LIABILITY SCALE':<15}")
    print(f"\t\t{'-' * 60}")
    print(f"\t\t{'h2 (Heritability)':<25} | {obs_metrics['h2']:<15} | {liab_metrics['h2']:<15}")
    print(f"\t\t{'Intercept':<25} | {obs_metrics['intercept']:<15} | {liab_metrics['intercept']:<15}")
    print(f"\t\t{'Ratio':<25} | {obs_metrics['ratio']:<15} | {liab_metrics['ratio']:<15}")
    print(f"\t\t{'-' * 60}")

    print(f"\n\t\t{Colors.BOLD}📂 OUTPUT FILES:{Colors.ENDC}")
    # Indent list items slightly more (2 tabs + 3 spaces for alignment)
    print(f"\t\t   • Munged:   {munged_file}")
    print(f"\t\t   • Logs:     {log_file}")
    print(f"\t\t   • H2 (Obs): {out_obs_log}")
    if samp_prev is not None:
        print(f"\t\t   • H2 (Liab):{out_liab_log}")
    print()
    return {
        "munged_sumstats": munged_file,
        "h2_liability": out_liab_log,
        "h2_observed": out_obs_log    }