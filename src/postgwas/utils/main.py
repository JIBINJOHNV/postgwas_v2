import os
import multiprocessing
from datetime import datetime
from pathlib import Path
import argparse
import glob
import sys
import subprocess
import shutil
from typing import Dict


## auto_detect_workers
## create_default_log(prefix="log", sample_id=None, outdir=None)
## validate_path( must_exist=True, must_be_file=False, must_be_dir=False, create_if_missing=False, must_not_be_empty=False,
    ## allowed_suffixes=None, dir_must_have_files=False, dir_require_nonempty_files=None, dir_nonempty_suffixes=None,
    # dir_forbid_empty_files=False)
    
## detect_total_memory_gb()

## validate_plink_prefix
## safe_thread_count(requested_threads, gb_per_thread = 20)
## run_cmd(cmd)




def validate_alphanumeric(value: str) -> str:
    import argparse, re
    if value is None:
        raise argparse.ArgumentTypeError("Value cannot be empty.")
    v = value.strip()
    if not re.match(r'^[A-Za-z0-9_-]+$', v):
        raise argparse.ArgumentTypeError(
            f"Invalid value '{value}': must contain only letters, digits, and underscores."
        )
    return v



def require_executable(name: str):
    """
    Ensure an external executable exists in PATH.
    """
    if shutil.which(name) is None:
        raise EnvironmentError(
            f"❌ Required executable not found in PATH: '{name}'\n"
            f"👉 Please install it or activate the correct environment."
        )


# ------------------------------------------------------------
# Reusable command runner with clean failure messages
# ------------------------------------------------------------
def run_cmd(cmd, check=True, shell=None, capture_output=True, text=True):
    """
    Simple, generalised subprocess runner.
    - If shell is not provided:
        * list command  -> shell=False
        * string command -> shell=True
    """
    # Auto-detect shell if not specified
    if shell is None:
        shell = isinstance(cmd, str)

    try:
        result = subprocess.run(
            cmd,
            check=check,
            shell=shell,
            capture_output=capture_output,
            text=text
        )
        return result

    except FileNotFoundError:
        raise FileNotFoundError(
            f"❌ Executable not found: {cmd}\n"
            "Ensure it is installed and available in your PATH."
        )

    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"❌ Command failed with exit code {e.returncode}:\n"
            f"   {cmd}\n"
            "--------------- stderr ---------------\n"
            f"{e.stderr}"
        ) from e



def infer_mount_root(*paths):
    """
    Infer the deepest common parent directory across all input paths.
    """
    resolved = [str(Path(p).expanduser().resolve()) for p in paths]
    common = os.path.commonpath(resolved)
    return Path(common)


def auto_detect_workers():
    """
    Universal CPU auto-detection.
    - Works on macOS, Linux, Windows
    - Respects SLURM limits
    - Respects Docker / cgroup v1 & v2 CPU quotas
    - Leaves 1–2 CPUs free for system stability
    """

    # ----------------------------------------------------
    # 1) SLURM limits (authoritative on HPC)
    # ----------------------------------------------------
    slurm_vars = [
        "SLURM_CPUS_PER_TASK",
        "SLURM_CPUS_ON_NODE",
        "SLURM_JOB_CPUS_PER_NODE",
    ]

    for var in slurm_vars:
        val = os.environ.get(var)
        if val:
            # Handle formats like "16", "16(x2)"
            try:
                cpus = int(val.split("(")[0])
                return max(1, cpus - 1)
            except ValueError:
                pass

    # ----------------------------------------------------
    # 2) cgroup v2 (modern Docker / Kubernetes)
    # ----------------------------------------------------
    try:
        cpu_max_path = "/sys/fs/cgroup/cpu.max"
        if os.path.exists(cpu_max_path):
            quota, period = open(cpu_max_path).read().strip().split()
            if quota != "max":
                cpus = int(quota) // int(period)
                if cpus > 0:
                    return max(1, cpus - 2)
    except Exception:
        pass

    # ----------------------------------------------------
    # 3) cgroup v1 (older Docker)
    # ----------------------------------------------------
    try:
        quota_path = "/sys/fs/cgroup/cpu/cpu.cfs_quota_us"
        period_path = "/sys/fs/cgroup/cpu/cpu.cfs_period_us"

        if os.path.exists(quota_path) and os.path.exists(period_path):
            quota = int(open(quota_path).read())
            period = int(open(period_path).read())

            if quota > 0 and period > 0:
                cpus = quota // period
                if cpus > 0:
                    return max(1, cpus - 2)
    except Exception:
        pass

    # ----------------------------------------------------
    # 4) OS fallback (macOS / Linux / Windows)
    # ----------------------------------------------------
    try:
        cpus = multiprocessing.cpu_count()
        return max(1, cpus - 1)
    except Exception:
        return 1


def auto_detect_ram_gb():
    """
    Universal RAM auto-detection.
    - Respects SLURM memory limits
    - Respects Docker / cgroup v1 & v2 memory limits
    - Works on Linux, macOS, Windows
    - Returns total usable RAM in GB (float)
    """

    # ----------------------------------------------------
    # 1) SLURM memory limits (authoritative on HPC)
    # ----------------------------------------------------
    slurm_mem = os.environ.get("SLURM_MEM_PER_NODE")
    if slurm_mem and slurm_mem.isdigit():
        # SLURM reports MB
        return float(slurm_mem) / 1024

    slurm_mem_cpu = os.environ.get("SLURM_MEM_PER_CPU")
    if slurm_mem_cpu and slurm_mem_cpu.isdigit():
        slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
        if slurm_cpus and slurm_cpus.isdigit():
            return (float(slurm_mem_cpu) * int(slurm_cpus)) / 1024

    # ----------------------------------------------------
    # 2) cgroup v2 memory limit (modern Docker / Kubernetes)
    # ----------------------------------------------------
    try:
        mem_max_path = "/sys/fs/cgroup/memory.max"
        if os.path.exists(mem_max_path):
            val = open(mem_max_path).read().strip()
            if val != "max":
                return int(val) / 1024**3
    except Exception:
        pass

    # ----------------------------------------------------
    # 3) cgroup v1 memory limit (older Docker)
    # ----------------------------------------------------
    try:
        mem_limit_path = "/sys/fs/cgroup/memory/memory.limit_in_bytes"
        if os.path.exists(mem_limit_path):
            val = int(open(mem_limit_path).read())
            # Docker sometimes reports a huge number when unlimited
            if val < 1 << 60:
                return val / 1024**3
    except Exception:
        pass

    # ----------------------------------------------------
    # 4) Linux host RAM (/proc/meminfo)
    # ----------------------------------------------------
    try:
        if os.path.exists("/proc/meminfo"):
            with open("/proc/meminfo") as f:
                for line in f:
                    if line.startswith("MemTotal"):
                        kb = int(line.split()[1])
                        return kb / 1024 / 1024
    except Exception:
        pass

    # ----------------------------------------------------
    # 5) macOS
    # ----------------------------------------------------
    try:
        ram_bytes = int(
            subprocess.check_output(["sysctl", "-n", "hw.memsize"]).strip()
        )
        return ram_bytes / 1024**3
    except Exception:
        pass

    # ----------------------------------------------------
    # 6) Absolute fallback
    # ----------------------------------------------------
    return 0.0



def safe_thread_count(requested_threads, gb_per_thread = 20):
    """
    Ensure the number of threads does not exceed available memory.
    Example: If each thread needs 20GB, auto-adjust threads based on RAM.
    Parameters
    ----------
    requested_threads : int
        Number of threads requested by user.
    gb_per_thread : int
        Minimum GB of RAM required per thread.
    Returns
    -------
    int
        Safe number of threads based on system memory.
    """
    try:
        import psutil
        total_ram_gb = psutil.virtual_memory().total / (1024**3)
    except ImportError:
        print("⚠️ psutil not installed → falling back to requested threads.")
        return requested_threads
    # Max threads based on memory
    max_threads = int(total_ram_gb // gb_per_thread)
    if max_threads < 1:
        print(f"❌ Not enough RAM: {total_ram_gb:.1f} GB available, "
              f"but require ≥ {gb_per_thread} GB per thread.")
        return 1
    if requested_threads > max_threads:
        print(f"            ⚠️ Reducing threads from {requested_threads} → {max_threads} "
              f"            (RAM available: {total_ram_gb:.1f} GB; {gb_per_thread} GB/thread)")
        return max_threads
    
    print(f"    ✔ Using {requested_threads} threads "
          f"    (RAM available: {total_ram_gb:.1f} GB; {gb_per_thread} GB/thread)")
    return requested_threads




def create_default_log(prefix="log", sample_id=None, outdir=None):
    """
    Create a unified default log filename for any module.

    Parameters
    ----------
    prefix : str
        Module name prefix (e.g., 'flames', 'susie', 'magma').
    sample_id : str or None
        Optional sample or phenotype name (e.g., 'PGC3_SCZ').
    outdir : str or Path or None
        Output directory where log should be placed. Default: current working directory.

    Returns
    -------
    str : Full path to log file.
    """
    if outdir is None:
        outdir = Path.cwd()
    else:
        outdir = Path(outdir).expanduser().resolve()

    outdir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    if sample_id:
        filename = f"{prefix}_{sample_id}_{timestamp}.log"
    else:
        filename = f"{prefix}_{timestamp}.log"

    return str(outdir / filename)


def validate_path(
    must_exist=False,
    must_be_file=False,
    must_be_dir=False,
    create_if_missing=False,
    must_not_be_empty=False,
    allowed_suffixes=None,

    # ---- NEW DIRECTORY RULES ----
    dir_must_have_files=False,
    dir_require_nonempty_files=None,        
    dir_nonempty_suffixes=None,            
    dir_forbid_empty_files=False,          
):
    """
    Generalized path validator with expanded directory/file checks.
    """

    def _validator(path_str):
        p = Path(path_str)

        # ------------------------------------------------------------
        # CREATE DIRECTORY IF MISSING
        # ------------------------------------------------------------
        if create_if_missing:
            if p.exists() and not p.is_dir():
                raise argparse.ArgumentTypeError(f"❌ Expected a directory but got file: {p}")
            if not p.exists():
                try:
                    p.mkdir(parents=True, exist_ok=True)
                except Exception as e:
                    raise argparse.ArgumentTypeError(f"❌ Cannot create directory {p}: {e}")
            return p

        # ------------------------------------------------------------
        # BASIC EXISTENCE + TYPE CHECKS
        # ------------------------------------------------------------
        if must_exist and not p.exists():
            raise argparse.ArgumentTypeError(f"❌ Path does not exist: {p}")

        if must_be_file and not p.is_file():
            raise argparse.ArgumentTypeError(f"❌ Expected a file but got: {p}")

        if must_be_dir and not p.is_dir():
            raise argparse.ArgumentTypeError(f"❌ Expected a directory but got: {p}")

        # ------------------------------------------------------------
        # FILE-SPECIFIC CHECKS
        # ------------------------------------------------------------
        if p.is_file():
            if must_not_be_empty and p.stat().st_size == 0:
                raise argparse.ArgumentTypeError(f"❌ File is empty: {p}")

            if allowed_suffixes:
                suffix = "".join(p.suffixes)
                allowed = {s if s.startswith(".") else f".{s}" for s in allowed_suffixes}
                if suffix not in allowed:
                    raise argparse.ArgumentTypeError(
                        f"❌ Invalid suffix '{suffix}' for {p}. Allowed: {sorted(allowed)}"
                    )
            return p

        # ------------------------------------------------------------
        # DIRECTORY-SPECIFIC CHECKS
        # ------------------------------------------------------------
        if p.is_dir():
            files = [f for f in p.iterdir() if f.is_file()]

            # ---- Require directory not to be empty ----
            if dir_must_have_files and len(files) == 0:
                raise argparse.ArgumentTypeError(f"❌ Directory contains no files: {p}")

            # ---- STRICT MODE: no empty files allowed ----
            if dir_forbid_empty_files:
                empty_files = [f for f in files if f.stat().st_size == 0]
                if empty_files:
                    raise argparse.ArgumentTypeError(
                        "❌ Directory contains empty files:\n" +
                        "\n".join(f"   - {f}" for f in empty_files)
                    )

            # ---- NEW RULE: At least N non-empty files ----
            if dir_require_nonempty_files is not None:
                # Decide which files count
                if dir_nonempty_suffixes:
                    allowed = {s if s.startswith(".") else f".{s}"
                               for s in dir_nonempty_suffixes}
                    relevant = [
                        f for f in files
                        if "".join(f.suffixes) in allowed
                    ]
                else:
                    relevant = files

                nonempty = [f for f in relevant if f.stat().st_size > 0]

                if len(nonempty) < dir_require_nonempty_files:
                    raise argparse.ArgumentTypeError(
                        f"❌ Directory requires at least "
                        f"{dir_require_nonempty_files} non-empty files "
                        f"(found {len(nonempty)}).\nDirectory: {p}"
                    )

            return p

        return p

    return _validator


# ------------------------------
# Reuse argparse validator manually
# ------------------------------
def apply_validator(validator, path_str):
    """
    Apply an argparse-style validator manually.
    Useful when validating paths from config/YAML (not from CLI).
    """
    try:
        return validator(path_str)
    except argparse.ArgumentTypeError as e:
        raise ValueError(str(e))




def validate_plink_prefix(prefix: str):
    prefix = Path(prefix)

    bed = prefix.with_suffix(".bed")
    bim = prefix.with_suffix(".bim")
    fam = prefix.with_suffix(".fam")

    missing = [p for p in (bed, bim, fam) if not p.exists()]
    if missing:
        missing_list = "\n".join([f"  - {str(m)}" for m in missing])
        raise FileNotFoundError(
            f"PLINK reference panel prefix is incomplete.\n"
            f"Missing files:\n{missing_list}"
        )
    return prefix






def validate_prefix_files(prefix: str, suffixes):
    """
    General-purpose validator:
    Accepts prefix + list of suffixes or wildcard patterns.
    Example suffixes:
        [".bed", ".bim", ".fam"]
        ["*.vcf.gz"]
        ["*_scores.txt.gz"]
    """
    prefix = Path(prefix)

    if isinstance(suffixes, str):
        suffixes = [suffixes]

    missing = []

    for suf in suffixes:

        # 1️⃣ Wildcard patterns: "*.vcf.gz", "*txt"
        if "*" in suf:
            pattern = str(prefix) + suf
            files = [Path(f) for f in glob.glob(pattern)]

        # 2️⃣ Direct suffixes: ".bed", ".bim", ".fam"
        else:
            # DO NOT use .with_suffix() → breaks prefixes with dots
            f = Path(str(prefix) + suf)
            files = [f] if f.exists() else []

        # Validate existence + non-empty
        valid = [f for f in files if f.exists() and f.stat().st_size > 0]

        if not valid:
            missing.append(str(prefix) + suf)

    if missing:
        raise FileNotFoundError("Missing required files:\n" + "\n".join(missing))

    return str(prefix)


def detect_total_memory_gb():
    """
    Safely detect total physical memory (RAM) in GB.
    Works on macOS, Linux, and Windows — never fails.
    """

    # Try psutil (best, cross-platform)
    try:
        import psutil
        return psutil.virtual_memory().total // (1024**3)
    except Exception:
        pass

    # Linux fallback: /proc/meminfo
    if os.path.exists("/proc/meminfo"):
        try:
            with open("/proc/meminfo") as f:
                for line in f:
                    if line.startswith("MemTotal:"):
                        kb = int(line.split()[1])
                        return (kb * 1024) // (1024**3)
        except Exception:
            pass

    # macOS fallback: sysctl
    if sys.platform == "darwin":
        try:
            out = subprocess.check_output(["sysctl", "-n", "hw.memsize"])
            return int(out.strip()) // (1024**3)
        except Exception:
            pass

    # Windows fallback (psutil usually works here)
    if sys.platform == "win32":
        try:
            import psutil
            return psutil.virtual_memory().total // (1024**3)
        except Exception:
            pass

    # Last safe fallback
    return 4  # assume at least 4GB





def decide_magma_batches_from_annot(
    annot_file: str,
    min_genes_per_batch: int = 1000,
    ram_per_cpu_gb: int = 16,
):
    """
    Decide MAGMA batch count using:
      - gene count from MAGMA .genes.annot
      - auto-detected CPUs (SLURM / Docker / local)
      - auto-detected RAM (SLURM / Docker / local)
    """
    def count_genes_from_magma_annot(annot_file):
        """
        Count genes from a MAGMA .genes.annot file.
        """
        count = 0
        with open(annot_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                count += 1
        return count

    # ----------------------------------------------------
    # 1) Count genes from MAGMA annotation file
    # ----------------------------------------------------
    total_genes = count_genes_from_magma_annot(annot_file)

    if total_genes <= 0:
        return 1

    # ----------------------------------------------------
    # 2) Auto-detect resources
    # ----------------------------------------------------
    cpus = auto_detect_workers()
    ram_gb = auto_detect_ram_gb()

    # ----------------------------------------------------
    # 3) Limit workers by RAM (≥ 8 GB per CPU)
    # ----------------------------------------------------
    max_workers_by_ram = max(1, int(ram_gb // ram_per_cpu_gb))

    # Effective usable workers
    workers = max(1, min(cpus, max_workers_by_ram))

    # ----------------------------------------------------
    # 4) Gene-based upper bound on batches
    #    (minimum 1000 genes per batch)
    # ----------------------------------------------------
    max_batches_by_genes = max(1, total_genes // min_genes_per_batch)

    # ----------------------------------------------------
    # 5) Final safe batch count
    # ----------------------------------------------------
    total_batches = max(1, min(workers, max_batches_by_genes))

    return total_batches
