import os
import polars as pl
from typing import Dict, Tuple

def gwastovcf(
    gwas_outputname: str,
    chr: str,
    genome_build: str,
    output_folder: str,
    fasta: str,
    dbsnp: str,
    aliasfile: str = "NA",
    main_script_path: str = "/app/main.py"
) -> Tuple[str, int]:
    """
    Build and execute the GWAS-to-VCF conversion command.

    Parameters
    ----------
    gwas_outputname : str
        GWAS dataset base name.
    chr : str
        Chromosome identifier.
    genome_build : str
        Genome build (e.g., GRCh37 or GRCh38).
    output_folder : str
        Path to the output directory.
    fasta : str
        Reference FASTA file path.
    dbsnp : str
        dbSNP VCF reference path.
    aliasfile : str, default="NA"
        Optional alias file path.
    main_script_path : str, default="/app/main.py"
        Path to the main GWAS-to-VCF conversion script.

    Returns
    -------
    Tuple[str, int]
        - The constructed shell command string.
        - The exit code (0 = success, nonzero = failure).
    """
    # --- Prepare paths ---
    output_subdir = os.path.join(output_folder, gwas_outputname)
    os.makedirs(output_subdir, exist_ok=True)
    
    input_file_name = f"{output_subdir}/{gwas_outputname}_chr{chr}_vcf_input.tsv"
    output_file_name = f"{output_subdir}/{gwas_outputname}_chr{chr}_build_{genome_build}_vcf"
    dict_file_name = f"{output_subdir}/{gwas_outputname}_chr{chr}.dict"
    error_log = f"{output_subdir}/{gwas_outputname}_gwastovcf_chr{chr}_errors.txt"
    
    # --- Validate required reference files ---
    for path, label in [(fasta, "FASTA"), (dbsnp, "dbSNP")]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"❌ Required {label} file not found → {path}")
    
    # --- Build command ---
    command_parts = [
        "python", main_script_path,
        "--data", input_file_name,
        "--ref", fasta,
        "--dbsnp", dbsnp,
        "--out", output_file_name,
        "--id", gwas_outputname,
        "--json", dict_file_name
    ]
    
    if aliasfile != "NA" and os.path.exists(aliasfile):
        command_parts += ["--alias", aliasfile]
    
    command = " ".join(command_parts) + f" 2>{error_log}"
    
    # --- Display and execute ---
    print(f"\n🚀 Running GWAS-to-VCF conversion for chr{chr}...\n{command}\n")
    
    exit_code = os.system(command)
    if exit_code == 0:
        print(f"✅ Conversion completed successfully for chr{chr}.")
    else:
        print(f"⚠️ Conversion failed for chr{chr}. Check error log: {error_log}")