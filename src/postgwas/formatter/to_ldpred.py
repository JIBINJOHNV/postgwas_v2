import subprocess
import os
import sys
import shutil
import traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

# ============================================================
# 1. WORKER FUNCTION
# ============================================================
import shutil
import subprocess
import os
from datetime import datetime

def run_ldpred_chromosome(args):
    """
    Worker function to process a single chromosome.
    Saves the executed command to a log file for debugging.
    """
    # Unpack arguments
    chrom, sumstat_vcf, sample_name, output_folder, bcftools_path = args

    chr_output = output_folder / f"{sample_name}_chr{chrom}_predld_input.tsv"
    log_file = output_folder / f"{sample_name}_chr{chrom}_bcftools.log"
    
    # Check for 'sed'
    if shutil.which("sed") is None:
        return f"❌ Chr {chrom} CRITICAL: 'sed' command not found in PATH."

    # Construct command using /bin/sh syntax
    cmd = f"""
    (
      echo "snp\\tchr\\tpos\\tA1\\tA2\\tbeta\\tSE\\tNC\\tSS\\tAF\\tLP\\tSI" && \
      bcftools view \
          -r {chrom} \
          --min-alleles 2 --max-alleles 2 "{sumstat_vcf}" | \
      bcftools query -f '%CHROM:%POS:%REF:%ALT\\t%CHROM\\t%POS\\t%ALT\\t%REF\\t[%ES]\\t[%SE]\\t[%NC]\\t[%SS]\\t[%AF]\\t[%LP]\\t[%SI]\\n' | \
      sed 's|:|_|g'
    ) > "{chr_output}"
    """
    
    # --- NEW: Save the command to a log file ---
    try:
        with log_file.open("w") as f:
            f.write(f"# Timestamp: {datetime.now()}\n")
            f.write(f"# Output File: {chr_output}\n")
            f.write("# Executed Command:\n")
            f.write(cmd)
    except Exception as e:
        # Don't crash the worker just because logging failed, but return a warning
        return f"⚠️ Chr {chrom} Log Error: {str(e)}"

    # Execute
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            env=os.environ.copy() 
        )
        
        if result.returncode != 0:
            # Append stderr to the log file for easier debugging
            with log_file.open("a") as f:
                f.write("\n\n# EXECUTION FAILED:\n")
                f.write(result.stderr)
            return f"❌ Chr {chrom} failed (Exit {result.returncode}). See {log_file.name}"
            
        return f"✅ Chr {chrom} done"

    except Exception as e:
        return f"❌ Chr {chrom} PYTHON ERROR: {str(e)}"

# ============================================================
# 2. MANAGER FUNCTION
# ============================================================

def vcf_to_ldpred(
    sumstat_vcf: str,
    output_folder: str,
    sample_name: str,
    bcftools_path: str,
    chrom_list=range(1, 23),
    nthreads: int = 8,
):
    """
    Orchestrates parallel execution.
    """
    output_folder = Path(output_folder) / "ldpred"
    output_folder.mkdir(parents=True, exist_ok=True)
    
    # --- CRITICAL FIX FOR [Errno 2] FileNotFoundError ---
    # If the script was launched from a directory that was later deleted,
    # multiprocessing will crash when trying to get os.getcwd().
    # We force the CWD to be the valid output folder we just created.
    try:
        os.getcwd()
    except (FileNotFoundError, OSError):
        print(f"   ⚠️  [Warning] Current working directory appears invalid/deleted.")
        print(f"       Switching CWD to: {output_folder}")
        os.chdir(output_folder)
    # ----------------------------------------------------

    # Prepare arguments
    task_args = [
        (chrom, sumstat_vcf, sample_name, output_folder, bcftools_path)
        for chrom in chrom_list
    ]
    
    messages = []
    with ProcessPoolExecutor(max_workers=nthreads) as executor:
        futures = {
            executor.submit(run_ldpred_chromosome, args): args[0]
            for args in task_args
        }
        
        for future in as_completed(futures):
            chrom = futures[future]
            try:
                msg = future.result()
                messages.append(msg)
                if "❌" in msg:
                    print(f"      {msg}")
            except Exception as e:
                print(f"💥 Worker Exception in Chr {chrom}:")
                traceback.print_exc() 
                raise e 
    
    return {
        "ldpred_folder": str(output_folder),
        "status": "completed",
        "messages": messages
    }