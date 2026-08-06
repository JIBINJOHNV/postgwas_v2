from pathlib import Path
import shutil
import sys
from postgwas.formatter.cli import get_formatter_parser
from postgwas.gene_assoc.magma_main import magma_analysis_pipeline
from postgwas.utils.main import run_cmd,require_executable



import sys
import shutil
from pathlib import Path
# assuming require_executable, magma_analysis_pipeline are imported from elsewhere

# =====================================================================
#   DIRECT MODE ENGINE
# =====================================================================
def run_magma_direct(args, ctx=None):
    print("\n               === Running MAGMA Gene & Geneset Analysis ===")
    
    final_magma_path = None

    # 1. Check if user provided a specific path
    if args.magma:
        if Path(args.magma).exists():
            final_magma_path = args.magma
        else:
            print(f"    ⚠️[WARNING] The provided path does not exist: {args.magma}")
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
    
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    outputs = magma_analysis_pipeline(
        output_dir=args.outdir,
        sample_id=args.sample_id,
        ld_ref=args.ld_ref,
        gene_loc_file=args.gene_loc_file,
        snp_loc_file=args.snp_loc_file,
        pval_file=args.pval_file,
        geneset_file=args.geneset_file,
        log_file=f"{args.outdir}/{args.sample_id}_magma_gene_assoc.log",
        threads=args.nthreads,
        window_upstream=args.window_upstream,
        window_downstream=args.window_downstream,
        gene_model=args.gene_model,
        n_sample_col=args.n_sample_col,
        seed=args.seed,
        magma=args.magma, # Now guaranteed to be valid
    )
    
    if ctx is not None:
        ctx["magma_gene"] = outputs
        
    print("\n           🎉 MAGMA analysis Completed.")
    print(" ")
    return outputs



