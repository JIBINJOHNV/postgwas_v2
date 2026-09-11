import subprocess
from pathlib import Path
from time import monotonic

from postgwas.core.processes import (
    SupervisedProcessTimeout,
    run_supervised_process,
)

from postgwas.modules.fine_mapping.defaults import (
    DEFAULT_BGEN_BITS,
    DEFAULT_EXTERNAL_TOOL_THREADS,
    DEFAULT_PLINK_MEMORY_MB,
)


class FinemapExternalToolTimeout(TimeoutError):
    """A FINEMAP-engine external stage exceeded its resolved timeout."""

    def __init__(
        self,
        reason,
        stage,
        timeout_error,
        log_file,
        configured_timeout_seconds=None,
    ):
        self.reason = str(reason)
        self.stage = str(stage)
        self.timeout_seconds = float(
            timeout_error.timeout_seconds
            if configured_timeout_seconds is None
            else configured_timeout_seconds
        )
        self.elapsed_seconds = timeout_error.elapsed_seconds
        self.forced_termination = timeout_error.forced_termination
        self.log_file = Path(log_file)
        super().__init__(
            "%s exceeded %.12g seconds (elapsed %.3f seconds; "
            "forced_termination=%s); see %s"
            % (
                self.stage,
                self.timeout_seconds,
                self.elapsed_seconds,
                self.forced_termination,
                self.log_file,
            )
        )


def _run_logged(
    command,
    log_file,
    cwd=None,
    *,
    timeout_seconds=None,
    termination_grace_seconds=None,
    timeout_reason=None,
    stage=None,
    configured_timeout_seconds=None,
    check=True,
):
    """Run a command and retain stdout/stderr for a reproducible failure report."""
    with Path(log_file).open("a") as handle:
        handle.write("CMD: " + " ".join(map(str, command)) + "\n")
        if timeout_seconds is not None:
            handle.write(
                "Timeout Seconds: %.12g\nTermination Grace Seconds: %.12g\n"
                % (float(timeout_seconds), float(termination_grace_seconds))
            )
            if configured_timeout_seconds is not None:
                handle.write(
                    "Configured Stage Timeout Seconds: %.12g\n"
                    % float(configured_timeout_seconds)
                )
    try:
        if timeout_seconds is None:
            result = subprocess.run(
                command, cwd=cwd, capture_output=True, text=True
            )
        else:
            result = run_supervised_process(
                command,
                cwd=cwd,
                timeout_seconds=float(timeout_seconds),
                termination_grace_seconds=float(termination_grace_seconds),
            )
    except SupervisedProcessTimeout as exc:
        with Path(log_file).open("a") as handle:
            handle.write("--- STDOUT ---\n" + exc.stdout)
            handle.write("\n--- STDERR ---\n" + exc.stderr)
            handle.write(
                "\nStatus: TIMEOUT\nElapsed Seconds: %.3f\n"
                "Forced Termination: %s\n"
                % (exc.elapsed_seconds, exc.forced_termination)
            )
        raise FinemapExternalToolTimeout(
            timeout_reason,
            stage,
            exc,
            log_file,
            configured_timeout_seconds=configured_timeout_seconds,
        ) from exc
    with Path(log_file).open("a") as handle:
        handle.write("--- STDOUT ---\n" + (result.stdout or ""))
        handle.write("\n--- STDERR ---\n" + (result.stderr or ""))
        handle.write(f"\nExit Code: {result.returncode}\n")
    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {result.returncode}; see {log_file}")
    return result


def run_plink_extraction(
    bfile_prefix: str,
    snp_file: Path,
    out_prefix: Path,
    *,
    plink_binary: str,
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
    *,
    plink_binary: str,
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


def run_bgen_indexing(bgen_file: Path, *, bgenix_binary: str):
    """Indexes BGEN file using bgenix."""
    log_file = bgen_file.with_suffix(".bgenix.log")
    if log_file.exists():
        log_file.unlink()
    _run_logged([
        str(bgenix_binary), "-g", str(bgen_file), "-index", "-clobber"
    ], log_file)
    bgi_file = Path(f"{bgen_file}.bgi")
    if not bgi_file.is_file() or bgi_file.stat().st_size == 0:
        raise RuntimeError(f"bgenix did not create a usable index; see {log_file}")




def run_ldstore(
    master_file: Path,
    threads: int = DEFAULT_EXTERNAL_TOOL_THREADS,
    *,
    ldstore_binary: str,
    timeout_seconds: float,
    termination_grace_seconds: float,
):
    """
    Runs LDstore to compute the LD matrix.
    Saves commands and logs to a file for debugging.
    """
    master_path = Path(master_file).resolve()
    # Create a specific log file for LDstore operations
    log_file = master_path.with_suffix(".ldstore.log")

    stage_started = monotonic()

    # Helper to keep both LDstore subprocesses inside one locus-stage budget.
    def _run_step(cmd, step_name):
        with open(log_file, "a") as f:
            f.write(f"\n{'='*20}\nRunning {step_name}\n{'='*20}\n")
        elapsed_seconds = monotonic() - stage_started
        remaining_seconds = float(timeout_seconds) - elapsed_seconds
        if remaining_seconds <= 0:
            with open(log_file, "a") as f:
                f.write(
                    "Status: TIMEOUT\nElapsed Seconds: %.3f\n"
                    "Forced Termination: False\n" % elapsed_seconds
                )
            timeout_error = SupervisedProcessTimeout(
                command=cmd,
                timeout_seconds=timeout_seconds,
                elapsed_seconds=elapsed_seconds,
                stdout="",
                stderr="",
                forced_termination=False,
            )
            raise FinemapExternalToolTimeout(
                "ldstore_timeout", "LDstore", timeout_error, log_file
            )

        return _run_logged(
            cmd,
            log_file,
            timeout_seconds=remaining_seconds,
            termination_grace_seconds=termination_grace_seconds,
            timeout_reason="ldstore_timeout",
            stage="LDstore",
            configured_timeout_seconds=timeout_seconds,
            check=False,
        )

    # --- Step 1: Compute correlations (write-bcor) ---
    cmd_bcor = [
        str(ldstore_binary),
        "--in-files", str(master_path),
        "--read-only-bgen", "--write-bcor",
        "--n-threads", str(threads)
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
        str(ldstore_binary),
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
    threads: int = DEFAULT_EXTERNAL_TOOL_THREADS,
    *,
    finemap_binary: str,
    timeout_seconds: float,
    termination_grace_seconds: float,
):
    """
    Runs the FINEMAP v1.4 binary.
    Saves command and logs to a file for debugging.
    """
    master_path = Path(master_file).resolve()
    # Log file will be something like: locus_id.finemap.log
    log_file = master_path.with_suffix(".finemap.log")

    # Use parent directory as working dir to avoid path length issues
    work_dir = master_path.parent
    master_filename = master_path.name

    algorithm = str(config["algorithm"])
    if algorithm not in {"sss", "cond"}:
        raise ValueError("FINEMAP algorithm must be 'sss' or 'cond'")

    cmd = [
        str(finemap_binary),
        f"--{algorithm}",
        "--in-files", master_filename,
        "--n-threads", str(threads),
        "--n-causal-snps", str(config["n_causal_snps"]),
        "--prob-cred-set", str(config["prob_cred_set"]),
        "--prior-std", str(config["prior_std"]),
    ]
    # The probability endpoint 1 means no SNP filtering. FINEMAP 1.4.2 uses
    # that endpoint when omitted but rejects an explicit --pvalue-snps 1.0.
    if float(config["pvalue_snps"]) != 1.0:
        cmd.extend(["--pvalue-snps", str(config["pvalue_snps"])])
    if algorithm == "sss":
        cmd.extend([
            "--n-iter", str(config["n_iter"]),
            "--n-conv-sss", str(config["n_conv_sss"]),
            "--prob-conv-sss-tol", str(config["prob_conv_sss_tol"]),
            "--n-configs-top", str(config["n_configs_top"]),
            "--corr-config", str(config["corr_config"]),
        ])
    else:
        cmd.extend(["--cond-pvalue", str(config["cond_pvalue"])])

    # Add optional flags if they exist in config
    if config.get("prior_k0"):
        cmd.extend(["--prior-k0", str(config["prior_k0"])])
    if config.get("prior_k"):
        cmd.append("--prior-k")
    if config.get("prior_snps"):
        cmd.append("--prior-snps")
    if config.get("std_effects"):
        cmd.append("--std-effects")
    if config.get("force_n_samples"):
        cmd.append("--force-n-samples")

    # 1. Log the Command
    with open(log_file, "w") as f:
        f.write(f"{'='*20}\nRunning FINEMAP\n{'='*20}\n")
        f.write(f"Work Dir: {work_dir}\n")
        f.write("Resolved SNP p-value threshold: %s%s\n" % (
            config["pvalue_snps"],
            " (native unfiltered endpoint; flag omitted)"
            if float(config["pvalue_snps"]) == 1.0 else "",
        ))

    # 2. Run inside work_dir so the master filename remains short.
    res = _run_logged(
        cmd,
        log_file,
        cwd=str(work_dir),
        timeout_seconds=timeout_seconds,
        termination_grace_seconds=termination_grace_seconds,
        timeout_reason="finemap_timeout",
        stage="FINEMAP",
        check=False,
    )

    # FINEMAP 1.4.2 may report its native error diagnostic with exit status 0.
    native_errors = [
        line.strip()
        for output in (res.stdout, res.stderr)
        for line in (output or "").splitlines()
        if line.lstrip().startswith("Error :")
    ]
    if native_errors:
        raise RuntimeError(
            "FINEMAP reported an error: %s. See log: %s"
            % ("; ".join(native_errors), log_file)
        )

    # 4. Error Checking
    if res.returncode != 0:
        # Check for non-fatal "No causal configuration" message
        if "No causal configuration found" in res.stdout:
            return False, "No causal SNPs found"

        # If it's a real crash, raise error pointing to the log
        raise RuntimeError(f"FINEMAP execution failed. See log: {log_file}")

    return True, "Success"
