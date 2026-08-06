
import os
import subprocess
from pathlib import Path
from typing import Tuple


import subprocess
import io
from pathlib import Path
from typing import Tuple


def gwastovcf(
    gwas_outputname: str,
    chromosome: str,
    grch_version: str,
    output_folder: str,
    fasta: str,
    dbsnp: str,
    aliasfile: str = "NA",
    main_script_path: str = "/app/main.py",
) -> Tuple[str, int]:
    """
    Silent + logged version of GWAS→VCF conversion.
    All stdout removed; all text captured to per-chromosome log file.

    Returns:
        (command_string, exit_code)
    """

    # -------------------------------------------------------
    # Setup paths
    # -------------------------------------------------------
    outdir = Path(output_folder)
    outdir.mkdir(parents=True, exist_ok=True)

    input_tsv = outdir / f"{gwas_outputname}_chr{chromosome}_vcf_input.tsv"
    output_vcf = outdir / f"{gwas_outputname}_chr{chromosome}_{grch_version}.vcf"
    json_dict = outdir / f"{gwas_outputname}_chr{chromosome}.dict"
    error_log = outdir / f"{gwas_outputname}_gwastovcf_chr{chromosome}_errors.txt"

    # -------------------------------------------------------
    # Setup per-chromosome log file
    # -------------------------------------------------------
    log_dir = Path(output_folder) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_gwastovcf.log"

    log_buffer = io.StringIO()

    def log_print(*args):
        log_buffer.write(" ".join(str(a) for a in args) + "\n")

    # -------------------------------------------------------
    # Input validation
    # -------------------------------------------------------
    for file_path, label in [
        (fasta, "FASTA"),
        (dbsnp, "dbSNP"),
        (main_script_path, "Main script")
    ]:
        if not Path(file_path).exists():
            msg = f"❌ Missing {label} file → {file_path}"
            log_print(msg)
            with open(log_file, "w") as f:
                f.write(log_buffer.getvalue())
            raise FileNotFoundError(msg)

    # -------------------------------------------------------
    # Build command
    # -------------------------------------------------------
    cmd = [
        "python", main_script_path,
        "--data", str(input_tsv),
        "--ref", str(fasta),
        "--out", str(output_vcf),
        "--id", str(gwas_outputname),
        "--json", str(json_dict)
    ]

    if aliasfile != "NA" and Path(aliasfile).exists():
        cmd += ["--alias", str(aliasfile)]

    command_str = " ".join(cmd)

    log_print(f"🧬 Running gwas2vcf for chr{chromosome}")
    log_print(f"Command: {command_str}")
    log_print(f"Output VCF: {output_vcf}")

    # -------------------------------------------------------
    # Run command
    # -------------------------------------------------------
    try:
        with open(error_log, "w") as errfile:
            result = subprocess.run(cmd, stderr=errfile, text=True)
        exit_code = result.returncode
    except Exception as e:
        log_print(f"❌ Execution failed: {e}")
        exit_code = 1

    # -------------------------------------------------------
    # Write execution status
    # -------------------------------------------------------
    if exit_code == 0:
        log_print(f"✅ chr{chromosome} GWAS-to-VCF conversion completed successfully.")
    else:
        log_print(f"⚠️ chr{chromosome} conversion failed. See error log: {error_log}")

    # -------------------------------------------------------
    # Write log file
    # -------------------------------------------------------
    with open(log_file, "w") as f:
        f.write(log_buffer.getvalue())

    return command_str, exit_code
