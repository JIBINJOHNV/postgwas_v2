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
from postgwas.modules.enrichment.cli import get_geneset_parser
from postgwas.core.ui import print_screen_message


def _configured_enrichment_python() -> Path:
    value = os.environ.get("POSTGWAS_ENRICHMENT_PYTHON")
    if not value:
        installed_runtime = (
            Path(sys.prefix)
            / "share"
            / "postgwas"
            / "environments"
            / "enrichment"
            / "bin"
            / "python"
        )
        if installed_runtime.is_file():
            value = str(installed_runtime)
        else:
            raise RuntimeError(
                "The pathway-enrichment runtime is not configured. Run "
                "tools/setup/install_postgwas.sh --all-tools and activate the "
                "created environment."
            )
    python = Path(value).expanduser().resolve()
    if not python.is_file() or not os.access(python, os.X_OK):
        raise RuntimeError(
            "The configured pathway-enrichment Python is not executable: %s"
            % python
        )
    return python


def _enrichment_runtime_environment(target_python: Path) -> dict[str, str]:
    runtime_prefix = target_python.parent.parent
    runtime_bin = target_python.parent
    runtime_r_home = runtime_prefix / "lib" / "R"
    runtime_r = runtime_bin / "R"
    if not runtime_r.is_file() or not os.access(runtime_r, os.X_OK):
        raise RuntimeError(
            "The pathway-enrichment runtime has no executable R installation: %s"
            % runtime_r
        )
    if not runtime_r_home.is_dir():
        raise RuntimeError(
            "The pathway-enrichment runtime has no R home directory: %s"
            % runtime_r_home
        )

    environment = os.environ.copy()
    environment.update({
        "POSTGWAS_ENRICHMENT_PYTHON": str(target_python),
        "PATH": os.pathsep.join(
            (str(runtime_bin), environment.get("PATH", ""))
        ),
        "R_HOME": str(runtime_r_home),
        "R_LIBS_USER": "/dev/null",
        "PYTHONNOUSERSITE": "1",
        "PYTHONPATH": "",
        "TMPDIR": custom_tmp,
        "TEMP": custom_tmp,
        "TMP": custom_tmp,
        "MPLCONFIGDIR": mpl_cache,
        "XDG_CACHE_HOME": xdg_cache,
        "OMNIPATH_CACHE_DIR": omnipath_cache,
    })
    return environment


def main() -> int:
    # 1. Initialize Parser & Show Help
    parser = get_geneset_parser(add_help=True)

    try:
        args = parser.parse_args()
    except SystemExit as exc:
        return int(exc.code or 0)

    # =========================================================
    # PHASE 2: ENVIRONMENT DISPATCHER
    # =========================================================
    try:
        target_python = _configured_enrichment_python()
        runtime_environment = _enrichment_runtime_environment(target_python)
    except RuntimeError as exc:
        print_screen_message("error", str(exc), stderr=True)
        return 1

    if Path(sys.executable).resolve() != target_python:
        print_screen_message("run", "Starting the isolated pathway-enrichment runtime...")
        completed = subprocess.run(
            [str(target_python), str(Path(__file__).resolve()), *sys.argv[1:]],
            env=runtime_environment,
            check=False,
        )
        return completed.returncode

    os.environ.update(runtime_environment)

    # =========================================================
    # PHASE 3: HEAVY WORKER LOGIC (Runs ONLY inside 'enricher')
    # =========================================================
    print_screen_message("run", "Starting pathway-enrichment analysis...")

    try:
        from postgwas.modules.enrichment.service import (
            run_multisource_enrichment_pipeline,
        )
        from postgwas.modules.enrichment.utils import load_gene_list
    except ImportError as e:
        print_screen_message("error", f"Critical Import Error in worker: {e}", stderr=True)
        sys.exit(1)

    # Run Analysis
    try:
        gene_list = load_gene_list(args.gene_input_file)

        run_multisource_enrichment_pipeline(
            gene_list=gene_list,
            output_dir=args.output_directory,
            sample_id=args.dataset_id,
            biogrid_access_key=args.biogrid_key,
            david_email=args.david_email,
            dsigdb_gmt=args.dsigdb_gmt,
            reference_set=args.reference_set,
            score_threshold=getattr(args, 'string_score', 400),
            fdr_thr=getattr(args, 'fdr_thr', 0.05),
        )
    except Exception as e:
        print_screen_message("error", f"Pipeline failed: {e}", stderr=True)
        return 1

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
