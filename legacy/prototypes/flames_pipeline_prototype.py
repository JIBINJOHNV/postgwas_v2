#!/usr/bin/env python3

from pathlib import Path
import os
import subprocess
from postgwas.modules.formatting.exporters.to_magma import vcf_to_magma
from postgwas.gene_tests.magma_main import magma_analysis_pipeline
from postgwas.modules.magmacovar.main import run_magma_covariates
from postgwas.modules.fine_mapping.cli import main as finemap_main
from postgwas.modules.pops.pops import pops_main, get_pops_args
from postgwas.modules.fine_mappingping.main import validate_locus_file, run_susie

def run_flames(
    annotation_dir,
    pops,
    magma_z,
    magma_tissue,
    indexfile,
    outdir,
    filename="FLAMES_scores",
    modelpath=None,
    snp_col="cred1",
    prob_col="prob1",
    distance=750000,
    weight=0.725,
):

    # Auto-detect FLAMES.py path
    flames_py = Path(__file__).parent / "FLAMES.py"

    # -----------------------------
    # 🔥 STEP 1 — FLAMES annotate
    # -----------------------------
    cmd_annot = f"""
    python "{flames_py}" annotate \
        --annotation_dir "{annotation_dir}" \
        --pops "{pops}" \
        --magma_z "{magma_z}" \
        --magma_tissue "{magma_tissue}" \
        --indexfile "{indexfile}" \
        --prob_col "{prob_col}" \
        --SNP_col "{snp_col}"
    """

    print("             🔥 Running FLAMES annotation…")
    subprocess.run(cmd_annot, shell=True, check=True)

    # -----------------------------
    # 🔥 STEP 2 — FLAMES scoring
    # -----------------------------
    cmd_score = f"""
    python "{flames_py}" FLAMES \
        --indexfile "{indexfile}" \
        --outdir "{outdir}" \
        --filename "{filename}" \
        --distance {distance} \
        --weight {weight} \
        --modelpath "{modelpath if modelpath else Path(__file__).parent / 'resources' / 'model'}"
    """
    print("                 FLAMES annotation… completed")
    print("             🔥 Running FLAMES scoring…")
    subprocess.run(cmd_score, shell=True, check=True)

    print("             🎉 FLAMES annotation + scoring completed.")


if __name__ == "__main__":
    raise SystemExit(
        "This legacy FLAMES prototype is not an executable entry point. "
        "Use `postgwas flames` instead."
    )
