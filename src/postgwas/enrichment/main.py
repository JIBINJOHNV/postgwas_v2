import sys
import os
import subprocess
from pathlib import Path

# =========================================================
# 🛑 CRITICAL SETUP: MUST RUN BEFORE ANY OTHER IMPORTS
# =========================================================
uid = os.getuid()
user_scratch_dir = f"/tmp/user_{uid}_postgwas_scratch"

custom_tmp = f"{user_scratch_dir}/tmp"
xdg_cache = f"{user_scratch_dir}/cache"
mpl_cache = f"{user_scratch_dir}/matplotlib"
omnipath_cache = f"{xdg_cache}/omnipathdb"

for d in (custom_tmp, xdg_cache, mpl_cache, omnipath_cache):
    Path(d).mkdir(parents=True, exist_ok=True)

# --- Global temp handling ---
os.environ["TMPDIR"] = custom_tmp
os.environ["TEMP"] = custom_tmp
os.environ["TMP"] = custom_tmp

# --- Cache handling ---
os.environ["XDG_CACHE_HOME"] = xdg_cache
os.environ["MPLCONFIGDIR"] = mpl_cache

# --- OmniPath (THIS IS THE MISSING PIECE) ---
os.environ["OMNIPATH_CACHE_DIR"] = omnipath_cache


# =========================================================
# NOW IT IS SAFE TO IMPORT YOUR PACKAGE
# =========================================================
try:
    from postgwas.enrichment.cli import get_geneset_parser
except ImportError as e:
    current_env = os.environ.get("CONDA_DEFAULT_ENV", "unknown")
    if current_env == "enricher":
        print(f"\n❌ FATAL ERROR inside 'enricher' environment.")
        print(f"   Could not import 'get_geneset_parser'.\n   Python Error: {e}")
        sys.exit(1)
        
    print(f"⚠️  CLI parser import failed in '{current_env}'. Switching to 'enricher'...")
    # Pass the 'enricher' environment, but we don't need to manually pass env vars
    # here because the child process is launched below in the dispatcher.
    sys.exit(subprocess.call(["micromamba", "run", "-n", "enricher", "python"] + sys.argv))

def main():
    # 1. Initialize Parser & Show Help
    parser = get_geneset_parser(add_help=True)
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    try:
        args = parser.parse_args()
    except SystemExit:
        sys.exit(0)

    # =========================================================
    # PHASE 2: ENVIRONMENT DISPATCHER
    # =========================================================
    current_env = os.environ.get("CONDA_DEFAULT_ENV", "unknown")
    target_env = "enricher"
    
    script_path = Path(__file__).resolve()
    project_root = str(script_path.parent.parent.parent) 

    if current_env != target_env:
        print(f"🚀 Arguments valid. Switching to '{target_env}' environment for analysis...")

        # We must explicitly pass our custom env vars to the child process
        cmd = [
            "micromamba", "run", "-n", target_env,
            "env", f"PYTHONPATH={project_root}",
            "env", f"TMPDIR={custom_tmp}",         # <--- Fixes general temp files
            "env", f"TEMP={custom_tmp}",           # <--- Windows/Legacy compatibility
            "env", f"TMP={custom_tmp}",            # <--- Windows/Legacy compatibility
            "env", f"MPLCONFIGDIR={mpl_cache}",    # <--- Fixes Matplotlib
            "env", f"XDG_CACHE_HOME={xdg_cache}",  # <--- Fixes OmniPath
            "python", str(script_path)
        ] + sys.argv[1:]

        try:
            subprocess.check_call(cmd)
            return 
        except subprocess.CalledProcessError as e:
            sys.exit(e.returncode)

    # =========================================================
    # PHASE 3: HEAVY WORKER LOGIC (Runs ONLY inside 'enricher')
    # =========================================================
    print(f"⚙️  Starting Analysis in '{current_env}'...")

    try:
        from postgwas.enrichment.workflows import run_multisource_enrichment_pipeline
        from postgwas.enrichment.utils import load_gene_list
    except ImportError as e:
        print(f"❌ Critical Import Error in worker: {e}")
        sys.exit(1)

    # Run Analysis
    try:
        gene_list = load_gene_list(args.gene_inputfile)
        print(f"🧬 Loaded {len(gene_list)} genes.")
        
        run_multisource_enrichment_pipeline(
            gene_list=gene_list,
            output_dir=args.outdir,
            sample_id=args.sample_id,
            biogrid_access_key=args.biogrid_key,
            david_email=args.david_email,
            dsigdb_gmt=args.dsigdb_gmt,
            score_threshold=getattr(args, 'string_score', 400),
            fdr_thr=getattr(args, 'fdr_thr', 0.05),
        )
    except Exception as e:
        print(f"❌ Pipeline failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()