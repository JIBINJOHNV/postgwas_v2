import subprocess
from pathlib import Path

from postgwas.finemap.defaults import (
    DEFAULT_BGEN_BITS,
    DEFAULT_CREDIBLE_SET_COVERAGE,
    DEFAULT_EXTERNAL_TOOL_THREADS,
    DEFAULT_PLINK_MEMORY_MB,
)

def _run_logged(command, log_file, cwd=None):
    """Run a command and retain stdout/stderr for a reproducible failure report."""
    result = subprocess.run(command, cwd=cwd, capture_output=True, text=True)
    with Path(log_file).open("a") as handle:
        handle.write("CMD: " + " ".join(map(str, command)) + "\n")
        handle.write("--- STDOUT ---\n" + (result.stdout or ""))
        handle.write("\n--- STDERR ---\n" + (result.stderr or ""))
        handle.write(f"\nExit Code: {result.returncode}\n")
    if result.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {result.returncode}; see {log_file}")
    return result


def run_plink_extraction(
    bfile_prefix: str,
    snp_file: Path,
    out_prefix: Path,
    plink_binary: str = "plink2",
    plink_memory_mb: int = DEFAULT_PLINK_MEMORY_MB,
):
    """Run PLINK 2 extraction and verify all BED outputs.

    ``plink_memory_mb`` is explicit and overrideable so the executed resource
    limit is captured by the caller's run configuration.
    """
    log_file = out_prefix.with_suffix(".plink.log")
    if log_file.exists():
        log_file.unlink()
    _run_logged([
        str(plink_binary), "--bfile", str(bfile_prefix),
        "--extract", str(snp_file),
        "--make-bed", "--out", str(out_prefix),
        "--silent", "--memory", str(int(plink_memory_mb))
    ], log_file)
    for suffix in [".bed", ".bim", ".fam"]:
        output = out_prefix.with_suffix(suffix)
        if not output.is_file() or output.stat().st_size == 0:
            raise RuntimeError(f"PLINK 2 did not create a usable {output.name}; see {log_file}")


def run_plink_to_bgen(
    bfile_prefix: Path,
    out_prefix: Path,
    plink_binary: str = "plink2",
    bgen_bits: int = DEFAULT_BGEN_BITS,
):
    """Convert the reconciled PLINK BED files to BGEN v1.2."""
    log_file = out_prefix.with_suffix(".plink.log")
    if log_file.exists():
        log_file.unlink()
    _run_logged([
        str(plink_binary), "--bfile", str(bfile_prefix),
        "--export", "bgen-1.2", f"bits={int(bgen_bits)}",
        "--out", str(out_prefix),
        "--silent"
    ], log_file)
    bgen_file = out_prefix.with_suffix(".bgen")
    if not bgen_file.is_file() or bgen_file.stat().st_size == 0:
        raise RuntimeError(f"PLINK 2 did not create a usable BGEN file; see {log_file}")


def run_bgen_indexing(bgen_file: Path):
    """Indexes BGEN file using bgenix."""
    log_file = bgen_file.with_suffix(".bgenix.log")
    if log_file.exists():
        log_file.unlink()
    _run_logged([
        "bgenix", "-g", str(bgen_file), "-index", "-clobber"
    ], log_file)
    bgi_file = Path(f"{bgen_file}.bgi")
    if not bgi_file.is_file() or bgi_file.stat().st_size == 0:
        raise RuntimeError(f"bgenix did not create a usable index; see {log_file}")




def run_ldstore(master_file: Path, n_threads: int = DEFAULT_EXTERNAL_TOOL_THREADS):
    """
    Runs LDstore to compute the LD matrix.
    Saves commands and logs to a file for debugging.
    """
    master_path = Path(master_file).resolve()
    # Create a specific log file for LDstore operations
    log_file = master_path.with_suffix(".ldstore.log")
    
    # Helper to log and run
    def _run_step(cmd, step_name):
        with open(log_file, "a") as f:
            f.write(f"\n{'='*20}\nRunning {step_name}\n{'='*20}\n")
            f.write(f"CMD: {' '.join(cmd)}\n\n")
        
        # Run process
        res = subprocess.run(cmd, capture_output=True, text=True)
        
        # Save Output to log
        with open(log_file, "a") as f:
            f.write("--- STDOUT ---\n")
            f.write(res.stdout)
            f.write("\n--- STDERR ---\n")
            f.write(res.stderr)
            f.write(f"\nExit Code: {res.returncode}\n")
            
        return res

    # --- Step 1: Compute correlations (write-bcor) ---
    cmd_bcor = [
        "ldstore", 
        "--in-files", str(master_path),
        "--read-only-bgen", "--write-bcor",
        "--n-threads", str(n_threads)
    ]
    
    # Clear previous log if it exists
    if log_file.exists():
        log_file.unlink()
        
    res_bcor = _run_step(cmd_bcor, "Step 1: Write BCOR")

    # Strict Error Checking
    if res_bcor.returncode != 0:
        raise RuntimeError(f"LDstore (write-bcor) failed. See log: {log_file}")
        
    if "Error" in res_bcor.stderr or "Error" in res_bcor.stdout:
         raise RuntimeError(f"LDstore (write-bcor) error detected in output. See log: {log_file}")

    # Verify BCOR file creation
    bcor_path = master_path.with_name(master_path.name.replace(".ldstore.master", ".bcor"))
    if not bcor_path.exists() or bcor_path.stat().st_size == 0:
         raise RuntimeError(f"LDstore succeeded but created empty/missing BCOR file: {bcor_path}")


    # --- Step 2: Convert to Text Matrix ---
    cmd_text = [
        "ldstore", 
        "--in-files", str(master_path),
        "--bcor-to-text"
    ]
    
    res_text = _run_step(cmd_text, "Step 2: BCOR to Text")

    if res_text.returncode != 0:
        raise RuntimeError(f"LDstore (bcor-to-text) failed. See log: {log_file}")
        
    if "Error" in res_text.stderr:
        raise RuntimeError(f"LDstore (bcor-to-text) error detected. See log: {log_file}")


def run_finemap_binary(
    master_file: Path,
    config: dict,
    n_threads: int = DEFAULT_EXTERNAL_TOOL_THREADS,
):
    """
    Runs the FINEMAP v1.4 binary.
    Saves command and logs to a file for debugging.
    """
    if float(config["prob_cred_set"]) != DEFAULT_CREDIBLE_SET_COVERAGE:
        raise ValueError("FINEMAP credible-set coverage must be exactly 0.95")
    master_path = Path(master_file).resolve()
    # Log file will be something like: locus_id.finemap.log
    log_file = master_path.with_suffix(".finemap.log")
    
    # Use parent directory as working dir to avoid path length issues
    work_dir = master_path.parent
    master_filename = master_path.name

    cmd = [
        "finemap", 
        "--sss", 
        "--in-files", master_filename,
        "--n-threads", str(n_threads),
        "--n-causal-snps", str(config["n_causal_snps"]),
        "--prob-cred-set", str(config["prob_cred_set"]),
        "--n-iter", str(config["n_iter"]),
        "--n-conv-sss", str(config["n_conv_sss"]),
        "--prob-conv-sss-tol", str(config["prob_conv_sss_tol"]),
    ]
    
    # Add optional flags if they exist in config
    if config.get("prior_k0"): 
        cmd.extend(["--prior-k0", str(config["prior_k0"])])
    if config.get("prior_k"): 
        cmd.append("--prior-k")
    if config.get("prior_snps"): 
        cmd.append("--prior-snps")
    if config.get("std_effects"): 
        cmd.append("--std-effects")

    # 1. Log the Command
    with open(log_file, "w") as f:
        f.write(f"{'='*20}\nRunning FINEMAP\n{'='*20}\n")
        f.write(f"Work Dir: {work_dir}\n")
        f.write(f"CMD: {' '.join(cmd)}\n\n")

    # 2. Run Process
    # We run inside 'work_dir' so we can just pass the filename, avoiding long path errors
    res = subprocess.run(cmd, cwd=str(work_dir), capture_output=True, text=True)
    
    # 3. Log the Output
    with open(log_file, "a") as f:
        f.write("--- STDOUT ---\n")
        f.write(res.stdout)
        f.write("\n--- STDERR ---\n")
        f.write(res.stderr)
        f.write(f"\nExit Code: {res.returncode}\n")

    # 4. Error Checking
    if res.returncode != 0:
        # Check for non-fatal "No causal configuration" message
        if "No causal configuration found" in res.stdout:
            return False, "No causal SNPs found"
            
        # If it's a real crash, raise error pointing to the log
        raise RuntimeError(f"FINEMAP execution failed. See log: {log_file}")
        
    return True, "Success"
