import subprocess
import argparse
from pathlib import Path
import sys
import pandas as pd
import time


def run_flames_direct(args):
    """
    Run FLAMES using precomputed:
      - SuSiE credible set files (directory)
      - MAGMA gene-level output
      - MAGMA covariate output
      - PoPS score file
    """
    # ==============================================================
    # Helper: Create FLAMES indexfile
    # ==============================================================
    def create_flames_indexfile(finemap_dir: Path, outdir: Path, indexfile_path: str):
        finemap_dir = Path(finemap_dir)

        if not finemap_dir.exists() or not finemap_dir.is_dir():
            raise NotADirectoryError(
                f"ERROR: --finemap_cred_dir must be a folder, not a file.\n"
                f"Provided: {finemap_dir}"
            )

        outdir.mkdir(parents=True, exist_ok=True)

        # Collect SuSiE credible set files
        finemap_files = sorted([
            f for f in finemap_dir.iterdir()
            if f.is_file() and not f.name.startswith(".")
        ])

        if len(finemap_files) == 0:
            raise RuntimeError(
                f"No SuSiE credible-set files found in folder: {finemap_dir}"
            )

        rows = []
        for f in finemap_files:
            prefix = f.stem
            annotfile = outdir / f"{prefix}_locus_1.txt"
            rows.append({
                "Filename": str(f.resolve()),
                "Annotfiles": str(annotfile.resolve())
            })

        df = pd.DataFrame(rows)
        df.to_csv(indexfile_path, sep="\t", index=False)

        return indexfile_path

    # ==============================================================
    # PREPARE DIRECT INPUTS
    # ==============================================================

    # Folder containing ALL SuSiE credible-set files
    finemap_dir = Path(args.finemap_cred_dir)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Create indexfile
    indexfile_path = outdir / f"{args.sample_id}_flames_index.tsv"

    indexfile = create_flames_indexfile(
        finemap_dir=finemap_dir,
        outdir=outdir,
        indexfile_path=str(indexfile_path)
    )

    # FLAMES inputs
    # Convert string "False" (from Docker/Argparse) back to boolean False
    cmd_vep = args.cmd_vep
    if str(cmd_vep).lower() == "false":
        cmd_vep = False

    vep_cache = args.vep_cache
    if str(vep_cache).lower() == "false":
        vep_cache = False
        
    annotation_dir = args.flames_annot_dir
    pops          = args.pops_score_file
    magma_z       = args.magma_genes_out
    magma_tissue  = args.magma_tissue_covar_results
    cmd_vep = args.cmd_vep
    vep_cache = args.vep_cache
    tabix = args.tabix
    CADD_file = args.CADD_file
    
    # Defaults for FLAMES model
    filename  = "FLAMES_scores"
    modelpath = Path(__file__).parent / "model"
    snp_col   = "cred1"
    prob_col  = "prob1"
    distance  = 750000
    weight    = 0.725

    flames_py = Path(__file__).parent / "FLAMES.py"

    # ==============================================================
    # STEP 1 — FLAMES annotate
    # ==============================================================
    cmd_annot = f"""
    python "{flames_py}" annotate \
        --annotation_dir "{annotation_dir}" \
        --pops "{pops}" \
        --magma_z "{magma_z}" \
        --magma_tissue "{magma_tissue}" \
        --indexfile "{indexfile}" \
        --prob_col "{prob_col}" \
        --SNP_col "{snp_col}" \
        --cmd_vep "{cmd_vep}" \
        --vep_cache "{vep_cache}" \
        --tabix "{tabix}" \
        --CADD_file  "{CADD_file}"
        
    """
    print("             🔥 Running FLAMES annotation…")

    MAX_RETRIES = 10          # total attempts
    RETRY_SLEEP = 150        # seconds (2 minutes)

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            subprocess.run(cmd_annot, shell=True, check=True)
            break  # ✅ success → exit loop

        except subprocess.CalledProcessError:
            print(
                f"⚠️  FLAMES annotation failed "
                f"(attempt {attempt}/{MAX_RETRIES}).\n"
                f"    Retrying in {RETRY_SLEEP // 60} minutes… "
                f"(likely temporary VEP / network issue)\n\n"
                "👉 Suggested fix:\n"
                "   Re-run FLAMES using local tools instead of APIs:\n"
                "     --cmd-vep --vep-cache\n"
                "     --CADD_file <path> --tabix <path>\n\n"
                "   This avoids VEP / CADD server outages and is recommended\n"
                "   for large or long-running analyses.\n"
            )

            if attempt < MAX_RETRIES:
                time.sleep(RETRY_SLEEP)
            else:
                print(
                    "❌ All retries exhausted. Exiting pipeline.\n"
                )
                raise

    # ==============================================================
    # STEP 2 — FLAMES scoring
    # ==============================================================
    cmd_score = f"""
    python "{flames_py}" FLAMES \
        --indexfile "{indexfile}" \
        --outdir "{outdir}" \
        --filename "{filename}" \
        --distance {distance} \
        --weight {weight} \
        --modelpath "{modelpath}"
    """
    print("             🔥 Running FLAMES scoring…")
    subprocess.run(cmd_score, shell=True, check=True)
    print("             🎉 FLAMES Analysis completed successfully!")
    print(" ")

