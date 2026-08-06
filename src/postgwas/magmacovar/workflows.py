import argparse
from pathlib import Path
import sys
import shutil

from postgwas.utils.main import (
    validate_path,
    detect_total_memory_gb,
    validate_alphanumeric,
    validate_prefix_files,
)

# MAGMA covariate engine (step 6)
from postgwas.magmacovar.main import run_magma_covariates

import argparse
import sys
import shutil
from pathlib import Path
# assuming run_magma_covariates is imported




# ==========================================================
# DIRECT MODE ENGINE
# ==========================================================
def run_magma_covar_direct(args: argparse.Namespace, ctx=None):
    """
    Run ONLY the MAGMA covariate (gene-property) analysis step.
    """
    print("\n           === Running MAGMA Covariate Analysis ===")
    
    # ---------------------------------------------------------
    # SMART MAGMA DETECTION
    # ---------------------------------------------------------
    final_magma_path = None

    # 1. Check if user provided a specific path
    if args.magma:
        if Path(args.magma).exists():
            final_magma_path = args.magma
        else:
            print(f"   ⚠️  [WARNING] The provided path does not exist: {args.magma}")
            print("       Attempting to auto-detect 'magma' in system PATH instead...")

    # 2. If we don't have a valid path yet (either not provided OR provided path was wrong), search PATH
    if not final_magma_path:
        detected_magma = shutil.which("magma")
        if detected_magma:
            print(f"   ✅  Found 'magma' in system PATH: {detected_magma}")
            final_magma_path = detected_magma

    # 3. Final Validation: If still missing, CRASH.
    if not final_magma_path:
        print("\n❌ [ERROR] MAGMA executable not found!")
        print("   The pipeline failed to find MAGMA in the provided argument or the system $PATH.")
        print("   Please either:")
        print("   1. Add 'magma' to your system environment, OR")
        print("   2. Provide a valid full path using --magma /path/to/magma")
        sys.exit(1)

    # Update args with the confirmed valid path
    args.magma = final_magma_path
    

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Prefix = <outdir>/<sample_id>
    output_prefix = str(outdir / args.sample_id)

    if args.covar_model:
        print("         ⚠ WARNING: Using a custom covariate model (--covar_model).")
        print("         Ensure this matches your scientific intent.\n")

    # ---------------------------------------------------------
    # Call MAGMA covariates engine
    # ---------------------------------------------------------
    run_magma_covariates(
        magma_bin=str(args.magma),
        gene_results_file=str(args.magama_gene_assoc_raw),
        covariates_file=str(args.covariates),
        output_prefix=output_prefix,
        model=args.covar_model,
        direction=args.covar_direction,
        log_file=getattr(args, "log_file", None),
    )
    
    output_file = f'{output_prefix}.gsa.out'
    
    print("           🎉 MAGMA Covariate Analysis Completed.")
    print(" ")
    
    if ctx is not None:
        ctx["magma_covar"] = output_file
        
    return output_file