

from pathlib import Path
import os

def clean_intermediate_files(output_dir: str, gwas_outputname):
    """
    Concatenate chromosome-wise VCFs for:
      1. GRCh37 successful liftovers
      2. GRCh38 successful liftovers
      3. *_notlifted.vcf.gz files (failed liftovers)
    Then delete all chromosome-level VCF files (and .tbi indices).

    Parameters
    ----------
    output_dir : str
        Directory containing per-chromosome VCF files.
    gwas_outputname : str
        Prefix used in filenames.
    mode : str, optional
        'concurrent' (default): Run GRCh37, GRCh38, and notlifted merges in parallel.
        'sequential': Run one after another.
    Returns
    -------
    dict
        Dictionary of merged output paths by category.
    """
    outdir = Path(output_dir)
    
    qc_dir = outdir / "qc_summary"
    qc_dir.mkdir(parents=True, exist_ok=True)
    # --- Merge log-like files (non-VCF) ---
    try:
        os.system(f"cat {outdir}/*_errors.txt > {qc_dir}/{gwas_outputname}_gwas2vcf_errors.txt")
    except Exception as e:
        print(f"⚠️ Skipped merging error logs: {e}")
    try:
        os.system(f"cat {outdir}/*.dict > {qc_dir}/{gwas_outputname}_gwas2vcf.dict")
    except Exception as e:
        print(f"⚠️ Skipped merging dict files: {e}")
    try:
        os.system(f"cat {outdir}/*_summary.tsv > {qc_dir}/{gwas_outputname}_gwas2vcf_summary.tsv")
    except Exception as e:
        print(f"⚠️ Skipped merging summary files: {e}")
    
    try:
        threads = os.cpu_count() or 4
        # 1. Create the specific output directory 
        target_dir = f"{qc_dir}/{gwas_outputname}_gwas2vcf_input"
        os.makedirs(target_dir, exist_ok=True)
        # 2. Move files into the new directory
        os.system(f"mv {outdir}/*_vcf_input.tsv {target_dir}/")
        # 3. Compress files individually inside the new folder
        os.system(f"pigz -p {threads} {target_dir}/*_vcf_input.tsv")
        # os.system(
        #     f"pigz -p {threads} -c {outdir}/*_vcf_input.tsv > "
        #     f"{qc_dir}/{gwas_outputname}_gwas2vcf_input.tsv.gz"
        # )
        # os.system(f"rm {outdir}/*_vcf_input.tsv")
    except Exception as e:
        print(f"⚠️ Skipped merging gwas2vcf files: {e}")
    
    try:
        os.system(f"rm {output_dir}/{gwas_outputname}*_chr*.tsv")
    except Exception as e:
        print(f"⚠️ Skipped merging gwas2vcf files: {e}")
    
    # --- Clean up intermediate per-chromosome logs ---
    try:
        os.system(f"rm -f {outdir}/*gwastovcf_chr*_errors.txt {outdir}/*_chr*.dict {outdir}/*_chr*_summary.tsv")
    except Exception as e:
        print(f"⚠️ Cleanup failed: {e}")
    # Placeholder: return paths for later use
    return {
        "errors": f"{qc_dir}/{gwas_outputname}_gwas2vcf_errors.txt",
        "dict": f"{qc_dir}/{gwas_outputname}_gwas2vcf.dict",
        "summary": f"{qc_dir}/{gwas_outputname}_gwas2vcf_summary.tsv",
    }