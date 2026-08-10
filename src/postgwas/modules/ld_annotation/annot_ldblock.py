#!/usr/bin/env python3
"""
annot_ldblock.py — Annotate GWAS VCFs with population-specific LD blocks.
"""

import subprocess,os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from postgwas.core.execution.runtime import run_cmd,require_executable





# ------------------------------------------------------------
# Main Function
# ------------------------------------------------------------
def annotate_ldblocks(
    vcf_path: str,
    output_directory: str,
    genome_build: str,
    ld_dir: str,
    populations=("EUR", "AFR", "EAS"),
    max_memory_gb='6G',
    threads=2,
    dataset_id="dataset",
) -> dict:
    """
    Sequentially annotate a VCF with multiple population LD-block INFO tags.
    Overwrites the input file safely. Cross-platform (macOS / Ubuntu).
    """
    # -------------------------------------------------
    # Dependency checks (FAIL FAST)
    # -------------------------------------------------
    require_executable("bcftools")
    require_executable("tabix")

    current_vcf = Path(vcf_path)
    if not current_vcf.exists():
        raise FileNotFoundError(f"❌ Input VCF not found: {current_vcf}")


    step_output_directory = Path(output_directory)
    step_output_directory.mkdir(parents=True, exist_ok=True)

    annotated_vcf = step_output_directory / f"{dataset_id}_ldblock.vcf.gz"
    os.system(f"cp {current_vcf} {annotated_vcf}")
    tmp_vcf = annotated_vcf.with_suffix(".tmp.vcf.gz")

    for pop in populations:
        # -------------------------------
        # Validate LD block file
        # -------------------------------
        bed = Path(ld_dir) / f"{genome_build}_{pop}_ldetect.bed.gz"
        if not bed.exists():
            raise FileNotFoundError(
                f"❌ LD block BED file missing:\n"
                f"   {bed}\n"
                f"Expected: {genome_build}_{pop}_ldetect.bed.gz\n"
                f"Check directory: {ld_dir}"
            )

        # -------------------------------
        # Build INFO header
        # -------------------------------
        info_id = f"{pop}_LDblock"
        info_desc = f"LD Block ID from Berisa & Pickrell 2016 ({pop}, {genome_build})"
        header_line = f'##INFO=<ID={info_id},Number=1,Type=String,Description="{info_desc}">'

        # -------------------------------
        # bcftools annotate command
        # -------------------------------
        cmd = [
            "bcftools", "annotate",
            "--threads", str(threads),
            "-a", str(bed),
            "-c", f"CHROM,FROM,TO,{info_id}",
            "-H", header_line,
            "-O", "z",
            "-o", str(tmp_vcf),
            str(annotated_vcf),
        ]
        run_cmd(cmd, shell=False)
        # Replace original VCF with annotated version
        if annotated_vcf.exists():
            annotated_vcf.unlink()
        tmp_vcf.rename(annotated_vcf)
        # Re-index
        run_cmd(["tabix", "-f", "-p", "vcf", str(annotated_vcf)], shell=False)
    return {
        "annotated_vcf": annotated_vcf
    }
