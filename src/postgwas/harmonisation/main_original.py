import os,re
import yaml
import sys,threading
import json
import subprocess
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple
from concurrent.futures import ProcessPoolExecutor
from postgwas.utils.main import validate_path,validate_prefix_files,apply_validator

import polars as pl
import pandas as pd
import textwrap
from pathlib import Path
import argparse

# ----------------------------------------------------------------------
# Thread-safe print
# ----------------------------------------------------------------------
print_lock = threading.Lock()


from postgwas.harmonisation.io import read_sumstats, read_config, find_resource_file_path
from postgwas.harmonisation.validator import validate_gwas_config
from postgwas.harmonisation.preprocess import split_chr_pos
from postgwas.harmonisation.chr_pos_process import fix_chr_pos_column
from postgwas.harmonisation.genome_build import genome_build
from postgwas.harmonisation.qc_utils import check_eaf_or_maf, compute_effective_n
from postgwas.harmonisation.parallel import parallel_split_and_run
from postgwas.harmonisation.reference import load_reference_build
from postgwas.harmonisation.eaf_maf import add_or_calculate_eaf
from postgwas.harmonisation.case_control import harmonize_sample_sizes
from postgwas.harmonisation.beta_or_or import is_beta_or_or
from postgwas.harmonisation.beta_and_se_from_z import calculate_beta_and_se_from_z
from postgwas.harmonisation.calculate_z_from_beta_se import calculate_z_from_beta_se
from postgwas.harmonisation.detect_pvalue_type_process import (
    detect_and_convert_pval,
    detect_pval_type,
    convert_pval_to_mlogp,
)
from postgwas.harmonisation.calculate_se_from_beta_pvalue import calculate_se_from_beta_pvalue
from postgwas.harmonisation.add_or_calculate_info import add_or_calculate_info
from postgwas.harmonisation.add_or_fix_snp_column import add_or_fix_snp_column
from postgwas.harmonisation.sumstat_to_vcf import (
    run_bcftools_munge,
    run_bcftools_annot,
    concat_vcfs_by_build,
)
from postgwas.harmonisation.export_gwas_sumstat_forgwas2vcf import export_gwas_sumstat
from postgwas.harmonisation.gwastovcf_gwas2vcf import gwastovcf
from postgwas.harmonisation.clean_intermediate import clean_intermediate_files
from postgwas.harmonisation.postgwas_qc import qc_json_to_dataframes
from postgwas.harmonisation.io import read_sumstats, read_config,load_default_config


from postgwas.sumstat_filter.sumstat_filter import filter_gwas_vcf_bcftools
from postgwas.qc_summary.main import run_qc_summary
from postgwas.harmonisation.utilities import combine_logs_per_chromosome
# ===============================================================
# Logging setup helper
# ===============================================================

def _get_chr_logger(sample_id: str, chromosome: str, log_dir: Path) -> logging.Logger:
    """
    Create / fetch a dedicated logger for a chromosome.
    Writes to {log_dir}/{sample_id}_chr{chromosome}.log
    """
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{sample_id}_chr{chromosome}.log"
    logger_name = f"postgwas.chr{chromosome}.{sample_id}"
    logger = logging.getLogger(logger_name)
    # Avoid attaching multiple handlers if called multiple times
    if not logger.handlers:
        logger.setLevel(logging.INFO)
        handler = logging.FileHandler(log_file, mode="w")
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.propagate = False  # no duplicate logs to root
    return logger


# ===============================================================
# 1️⃣ Resource Map Builder
# ===============================================================

def build_resource_map(
    chromosome: str,
    grch_version: str,
    resource_folder: str,
    user_eaf_file: str,
    default_eaf_file: str,
    default_comparison_af_file: str,
    user_info_file: Optional[str] = "NA",
    default_info_file: Optional[str] = "NA",
    dbsnp: str = "dbSNP",
    user_eaf_column: Optional[str] = "NA",
    default_eaf_column: Optional[str] = "NA",
    default_comparison_af_column: Optional[str] = "NA",
    user_info_column: Optional[str] = "NA",
    default_info_column: Optional[str] = "NA",
) -> Dict[str, Any]:
    """
    Construct all resource file paths for one chromosome.
    User-provided EAF/INFO files override defaults if provided.
    """
    # --- Allele frequency (EAF) resources ---
    user_eaf_path = find_resource_file_path(
        input_file=user_eaf_file,
        resource_folder=resource_folder,
        grch_version=grch_version,
        chromosome=chromosome,
    )
    default_eaf_path = (
        f"{resource_folder}/{grch_version}/default_af/tab_files/"
        f"{grch_version}_{default_eaf_file}_freq_chr{chromosome}.tsv.gz"
    )
    default_comparison_af_path = (
        f"{resource_folder}/{grch_version}/external_af/vcf_files/"
        f"{grch_version}_{default_comparison_af_file}_freq_chr{chromosome}.vcf.gz"
    )
    # --- INFO score resources ---
    user_info_path = find_resource_file_path(
        input_file=user_info_file,
        resource_folder=resource_folder,
        grch_version=grch_version,
        chromosome=chromosome,
    )
    default_info_path = (
        f"{resource_folder}/{grch_version}/default_infoscore/tab_files/"
        f"{grch_version}_{default_info_file}_infoscore_chr{chromosome}.tsv.gz"
    )
    # --- Core genome references ---
    genome_fasta_path = (
        f"{resource_folder}/{grch_version}/fasta_files/{grch_version}_chr{chromosome}.fa"
    )
    dbsnp_path = (
        f"{resource_folder}/{grch_version}/dbSNP/vcf_files/"
        f"{grch_version}_{dbsnp}_chr{chromosome}.vcf.gz"
    )
    annot_path = (
        f"{resource_folder}/{grch_version}/gff_files/{grch_version}_ensembl.gff3.gz"
    )
    
    # --- Chain and target FASTA for liftover ---
    if grch_version == "GRCh37":
        target_fasta = f"{resource_folder}/GRCh38/fasta_files/GRCh38_chr{chromosome}.fa"
        chain_file = f"{resource_folder}/chain_files/GRCh37_to_GRCh38.chain"
    else:
        target_fasta = f"{resource_folder}/GRCh37/fasta_files/GRCh37_chr{chromosome}.fa"
        chain_file = f"{resource_folder}/chain_files/GRCh38_to_GRCh37.chain"

    file_validator = validate_path( must_exist=True, must_be_file=True, must_not_be_empty=True )
    errors = []
    for path in [default_eaf_path, default_info_path, default_comparison_af_path, dbsnp_path ]:
        try:
            file_validator(path)
        except argparse.ArgumentTypeError as e:
            errors.append(f"❌ File invalid: {path}\n   → {e}")
    if errors:
        raise ValueError(
            "\n❌ One or more default files are invalid:\n" +
            "\n".join(errors)
        )
    
    return {
        # EAF
        "user_eaf_file": user_eaf_path,
        "default_eaf_file": default_eaf_path,
        "default_comparison_af_file": default_comparison_af_path,
        "user_eaf_column": user_eaf_column,
        "default_eaf_column": default_eaf_column,
        "default_comparison_af_column": default_comparison_af_column,
        # INFO
        "user_info_file": user_info_path,
        "default_info_file": default_info_path,
        "user_info_column": user_info_column,
        "default_info_column": default_info_column,
        # Genome, dbSNP, chain
        "genome_fasta_file": genome_fasta_path,
        "dbsnp_file": dbsnp_path,
        "annot_path": annot_path,
        "target_fasta": target_fasta,
        "chain_file": chain_file,
    }


# ===============================================================
# 2️⃣ Per-Chromosome Harmonization Pipeline
# ===============================================================

def process_one_chromosome(
    chromosome: str,
    chr_file: str,
    sample_column_dict: Dict[str, Any],
    resource_folder: str,
    grch_version: str,
    user_eaf_file: str,
    default_eaf_file: str,
    default_comparison_af_file: str,
    user_info_file: Optional[str] = "NA",
    default_info_file: Optional[str] = "NA",
    user_eaf_column: Optional[str] = "NA",
    default_eaf_column: Optional[str] = "NA",
    default_comparison_af_column: Optional[str] = "NA",
    user_info_column: Optional[str] = "NA",
    default_info_column: Optional[str] = "NA",
    dbsnp: str = "dbSNP",
    output_dir: str = ".",
    threads: int = 5,
    gwastovcf_main_script_path: str = "/app/main.py",
) -> Tuple[str, Dict[str, Any]]:
    """
    Run full harmonisation + GWAS→VCF + bcftools annotation for a single chromosome.
    Returns:
        (chromosome, qc_dict)
    """
    sample_id = sample_column_dict["gwas_outputname"]
    output_dir_path = Path(output_dir)
    log_dir = output_dir_path / "logs"
    logger = _get_chr_logger(sample_id=sample_id, chromosome=chromosome, log_dir=log_dir)
    qc_dict: Dict[str, Any] = {}
    try:
        logger.info(f"🚀 Starting harmonization pipeline for chromosome {chromosome}")
        logger.info(f"Reading per-chromosome file: {chr_file}")
        # ------------------------------
        # 1. Load per-chromosome file
        # ------------------------------
        df = pl.read_csv(chr_file, separator="\t")
        logger.info(f"Loaded {df.height:,} rows × {df.width} columns")
        df = fix_chr_pos_column(chromosome=chromosome, df=df, sample_column_dict=sample_column_dict, drop_mt=True)
        logger.info("Chromosome/POS columns fixed (MT dropped if present).")
        # ------------------------------
        # 2. Build resource map
        # ------------------------------
        res = build_resource_map(
            chromosome=chromosome,
            grch_version=grch_version,
            resource_folder=resource_folder,
            user_eaf_file=user_eaf_file,
            default_eaf_file=default_eaf_file,
            default_comparison_af_file=default_comparison_af_file,
            user_info_file=user_info_file,
            default_info_file=default_info_file,
            user_eaf_column=user_eaf_column,
            default_eaf_column=default_eaf_column,
            default_comparison_af_column=default_comparison_af_column,
            user_info_column=user_info_column,
            default_info_column=default_info_column,
            dbsnp=dbsnp,
        )
        logger.info("Resource map constructed.")
        # ------------------------------
        # 3. Harmonization steps
        # ------------------------------
        df, eaf_qc, sample_column_dict = add_or_calculate_eaf(
            chromosome=chromosome,
            df=df,
            sample_column_dict=sample_column_dict,
            eaffile=res["user_eaf_file"],
            default_eaf_file=res["default_eaf_file"],
            default_eaf_eafcolumn=res["default_eaf_column"],
        )
        logger.info("EAF harmonization completed.")
        df, size_qc, sample_column_dict = harmonize_sample_sizes(chromosome=chromosome,df=df, sample_column_dict=sample_column_dict)
        logger.info("Sample size harmonization completed.")
        df, beta_qc, sample_column_dict = calculate_beta_and_se_from_z(chromosome=chromosome,df=df, sample_column_dict=sample_column_dict)
        logger.info("BETA/SE from Z completed.")
        df, or_qc, sample_column_dict = is_beta_or_or(chromosome=chromosome,df=df, sample_column_dict=sample_column_dict)
        logger.info("Effect-type (BETA/OR) harmonization completed.")
        df, pval_qc, sample_column_dict = detect_and_convert_pval(chromosome=chromosome,df=df, sample_column_dict=sample_column_dict)
        logger.info("P-value detection and conversion completed.")
        df, se_qc, sample_column_dict = calculate_se_from_beta_pvalue(chromosome=chromosome,df=df, sample_column_dict=sample_column_dict)
        logger.info("SE from BETA/P-value calculation completed.")
        df, z_qc, sample_column_dict = calculate_z_from_beta_se(chromosome=chromosome,df=df, sample_column_dict=sample_column_dict)
        logger.info("Z from BETA/SE calculation completed.")
        df, info_qc, sample_column_dict = add_or_calculate_info(
            chromosome=chromosome,
            df=df,
            sample_column_dict=sample_column_dict,
            default_info_file=res["default_info_file"],
            info_file=res["user_info_file"],
            default_info_column=res["default_info_column"],
            info_column=res["user_info_column"],
        )
        logger.info("INFO score harmonization completed.")
        df, sample_column_dict = add_or_fix_snp_column(chromosome=chromosome, df=df, sample_column_dict=sample_column_dict)
        logger.info("SNP column added/fixed.")
        # ------------------------------
        # 4. Export cleaned GWAS sumstats
        # ------------------------------
        export_gwas_summary_df = export_gwas_sumstat(
            df=df,
            sample_column_dict=sample_column_dict,
            output_dir=str(output_dir_path),
            gwas_outputname=sample_column_dict["gwas_outputname"],
            chromosome=chromosome,
            genome_build=grch_version,
        )
        logger.info("Exported cleaned GWAS summary statistics.")
        # ------------------------------
        # 5. Convert GWAS to VCF (gwastovcf)
        # ------------------------------
        gwastovcf_command_str, gwastovcf_exit_code = gwastovcf(
            gwas_outputname=sample_column_dict["gwas_outputname"],
            chromosome=chromosome,
            grch_version=grch_version,
            output_folder=str(output_dir_path),
            fasta=res["genome_fasta_file"],
            dbsnp=res["dbsnp_file"],
            aliasfile="NA",
            main_script_path=gwastovcf_main_script_path,
        )
        logger.info(
            f"GWAS-to-VCF conversion finished with exit code {gwastovcf_exit_code}. "
            f"Command: {gwastovcf_command_str}"
        )
        # ------------------------------
        # 6. Annotate VCF (bcftools)
        # ------------------------------
        logger.info("Starting bcftools annotation + liftover.")
        run_bcftools_annot(
            output_dir=str(output_dir_path),
            gwas_outputname=sample_column_dict["gwas_outputname"],
            chromosome=chromosome,
            external_eaf_file=res["default_comparison_af_file"],
            default_dbsnp_file=res["dbsnp_file"],
            genome_fasta_file=res["genome_fasta_file"],
            target_genome_fasta_file=res["target_fasta"],
            gff_file=res["annot_path"],
            chain_file=res["chain_file"],
            grch_version=grch_version,
            threads=threads,
        )
        logger.info("bcftools annotation + liftover completed.")
        # ------------------------------
        # 7. QC summary
        # ------------------------------
        qc_dict = {
            "eaf_qc": eaf_qc,
            "sample_size_qc": size_qc,
            "beta_qc": beta_qc,
            "or_qc": or_qc,
            "pval_qc": pval_qc,
            "se_qc": se_qc,
            "z_qc": z_qc,
            "info_qc": info_qc,
            "gwastovcf_exit_code": gwastovcf_exit_code,
        }
        logger.info(f"✅ Completed chromosome {chromosome}")
    except Exception as e:
        logger.exception(f"❌ ERROR in chromosome {chromosome}: {e}")
        qc_dict = {"error": str(e)}
    return chromosome, qc_dict




# ----------------------------------------------------------------------
# Thread-safe print
# ----------------------------------------------------------------------
print_lock = threading.Lock()

def safe_print(*args, **kwargs):
    with print_lock:
        print(*args, **kwargs)

# ----------------------------------------------------------------------
# Helper: Chromosome extractor
# ----------------------------------------------------------------------
def extract_chromosome_from_filename(filepath: str) -> str:
    """
    Extract '1'..'22','X','Y','MT' from filenames like:
    PGC3_SCZ_european_chr5_vcf_input.tsv
    """
    m = re.search(r"chr(\d+|X|Y|MT)", Path(filepath).stem, re.I)
    if not m:
        raise ValueError(f"Cannot extract chromosome from filename: {filepath}")
    return m.group(1).upper()


# ===============================================================
#  MULTIPROCESSING VERSION — Python 3.8 SAFE
# ===============================================================

def gwas_to_vcf_parallel(
    sumstat_file,
    sample_column_dict,
    output_dir,
    resource_folder,
    default_eaf_file="1000G",
    default_info_file="1000G",
    default_comparison_af_file="1000G",
    dbsnp="dbsnp155",
    required_columns_in_sumstat=None,
    optional_columns_in_sumstat=None,
    max_workers=5,
    user_eaf_file="NA",
    user_eaf_column="EUR",
    user_info_file="NA",
    user_info_column="INFO",
    default_eaf_column="EUR",
    default_comparison_af_column="EUR",
    default_info_column="INFO",
    gwastovcf_main_script_path="/app/main.py",
    grch37_file=None,
    grch38_file=None,
):
    # ------------------------------
    # Create output path
    # ------------------------------
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    # ------------------------------
    # Read file + validate config
    # ------------------------------
    df, file_cvariant_count, polars_rows = read_sumstats(
        sumstat_file=sumstat_file,
        output_dir=str(output_dir_path),
    )
    # safe_print(
    #     f"Reading: {sumstat_file}\n"
    #     f"Polars rows: {polars_rows:,} | Shell-count variants: {file_cvariant_count:,}"
    # )
    validate_gwas_config(sample_column_dict, df)
    # -----------------------------
    # Genome build inference
    # ------------------------------
    genome_build_info = genome_build(df, grch37_file, grch38_file, sample_column_dict)
    grch_version = genome_build_info["inferred_build"]
    safe_print(f"🧬 Inferred genome build: {grch_version}")
    # ------------------------------
    # Split by chromosome
    # ------------------------------
    chr_files = split_chr_pos(df, sample_gwas_dict=sample_column_dict)
    final_chr_files = []
    for f in chr_files:
        try:
            _ = extract_chromosome_from_filename(f)
            final_chr_files.append(f)
        except Exception:
            continue
    #safe_print(f"🧩 Found {len(final_chr_files)} chromosome files to process.\n")
    # ------------------------------
    # QC container
    # ------------------------------
    per_chr_qc = {}
    per_chr_qc["total_variant_infile"] = file_cvariant_count
    per_chr_qc["total_variant_read"] = polars_rows
    # (You’re not using required/optional_columns_in_sumstat here, that's fine.)
    # ------------------------------
    # Run each chromosome in multiprocessing
    # ------------------------------
    try:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_chr = {}
            for chr_file in final_chr_files:
                chrom = extract_chromosome_from_filename(chr_file)
                future = executor.submit(
                    process_one_chromosome,
                    chrom,
                    chr_file,
                    sample_column_dict.copy(),
                    resource_folder,
                    grch_version,
                    user_eaf_file,
                    default_eaf_file,
                    default_comparison_af_file,
                    user_info_file,
                    default_info_file,
                    user_eaf_column,
                    default_eaf_column,
                    default_comparison_af_column,
                    user_info_column,
                    default_info_column,
                    dbsnp,
                    output_dir,
                    5,
                    gwastovcf_main_script_path,
                )
                future_to_chr[future] = chrom
            # --------------------------
            # Collect results safely
            # --------------------------
            for future in as_completed(future_to_chr):
                chrom = future_to_chr[future]
                try:
                    chrom_res, qc = future.result()
                    per_chr_qc[chrom_res] = qc
                    if "error" in qc:
                        safe_print(f"❌ Chromosome {chrom_res} completed with ERROR.")
                    else:
                        pass
                        #safe_print(f"✅ Completed chromosome {chrom_res}")
                except Exception as exc:
                    #safe_print(f"❌ [ERROR] Chromosome {chrom} crashed: {exc}")
                    per_chr_qc[chrom] = {"error": str(exc)}
    finally:
        pass
        #safe_print("\n📦 Running concatenation and cleanup steps...\n")
    return per_chr_qc




def run_harmonisation_pipeline(
    sample_column_dict,
    default_cfg,
    nthreads: int
):
    """
    Full harmonisation + gwas2vcf + merge + cleanup pipeline.
    Clean and modular version to be called from cli.py.
    """
    # # -----------------------------------------------------
    # # Load DEFAULT config
    # # -----------------------------------------------------
    try:
        default_cfg = load_default_config(default_cfg)
        if not default_cfg:
            raise ValueError("load_config returned empty output.")
    except Exception as e:
        print(f"❌ ERROR: Failed to load default config file '{default_cfg}'.")
        print(f"   Reason: {e}")
        sys.exit(1)
    # ---------------------------------------------------------
    #  Extract all defaults safely
    # ---------------------------------------------------------
    default_eaf            = default_cfg["default_eaf"]
    default_info           = default_cfg["default_info"]
    default_comparison_af  = default_cfg["default_comparison_af"]

    default_eaf_column     = default_cfg["default_eaf_column"]
    default_info_column    = default_cfg["default_info_column"]
    comparison_cols        = default_cfg["default_comparison_af_column"]
    dbsnp                  = default_cfg["default_dbsnp"]
    required_cols          = default_cfg["required_columns_in_sumstat"]
    optional_cols          = default_cfg["optional_columns_in_sumstat"]
    gwastovcf_script       = default_cfg["gwastovcf_main_script_path"]
    
    # Build clean path using pathlib
    output_folder = (Path(sample_column_dict['output_folder'])/sample_column_dict['gwas_outputname']/"1_harmonisation")
    # Save normalized string version
    sample_column_dict['output_folder'] = str(output_folder)

    # Try to create it safely
    try:
        output_folder.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        print("\n❌ ERROR: Failed to create output directory.")
        print(f"   Path: {output_folder}")
        print(f"   Reason: {e}\n")
        sys.exit(1)


    resource_folder = Path(sample_column_dict["resourse_folder"])
    # ---- Validate existence ----
    if not resource_folder.exists():
        raise ValueError(
            f"❌ ERROR: Resource folder does not exist:\n   {resource_folder}"
        )
    # ---- Validate directory ----
    if not resource_folder.is_dir():
        raise ValueError(
            f"❌ ERROR: Resource folder is not a directory:\n   {resource_folder}"
        )
    # ---- Validate non-empty ----
    if not any(resource_folder.iterdir()):
        raise ValueError(
            f"❌ ERROR: Resource folder is empty:\n   {resource_folder}"
        )

    file_validator = validate_path( must_exist=True, must_be_file=True, must_not_be_empty=True )
    errors = []
    for path in [sample_column_dict["sumstat_file"],
                 f"{sample_column_dict['resourse_folder']}/GRCh37_38_check_files/GRCh37_check_file.tsv",
                 f"{sample_column_dict['resourse_folder']}/GRCh37_38_check_files/GRCh38_check_file.tsv" ]:
        try:
            file_validator(path)
        except argparse.ArgumentTypeError as e:
            errors.append(f"❌ File invalid: {path}\n   → {e}")
    if errors:
        raise ValueError(
            "\n❌ One or more default files are invalid:\n" +
            "\n".join(errors)
        )
    
    # ---------------------------------------------------------
    # STEP 1 — FULL harmonisation + per-chromosome GWAS2VCF
    # ---------------------------------------------------------
    qc_results = gwas_to_vcf_parallel(
            sumstat_file = sample_column_dict["sumstat_file"],
            sample_column_dict = sample_column_dict,
            output_dir = sample_column_dict['output_folder'],
            resource_folder = sample_column_dict['resourse_folder'],

            default_eaf_file = default_eaf,
            default_info_file = default_info,
            default_comparison_af_file = default_comparison_af, 
            dbsnp = dbsnp,
            required_columns_in_sumstat = required_cols,
            optional_columns_in_sumstat = optional_cols,
            max_workers = nthreads,

            user_eaf_file = sample_column_dict['eaffile'],
            user_eaf_column = sample_column_dict['eafcolumn'],
            user_info_file = sample_column_dict['infofile'],
            user_info_column = sample_column_dict['infocolumn'],

            default_eaf_column = default_eaf_column,
            default_comparison_af_column = comparison_cols,
            default_info_column = default_info_column,

            grch37_file = f"{sample_column_dict['resourse_folder']}/GRCh37_38_check_files/GRCh37_check_file.tsv",
            grch38_file = f"{sample_column_dict['resourse_folder']}/GRCh37_38_check_files/GRCh38_check_file.tsv",

            gwastovcf_main_script_path = gwastovcf_script 
    )
    
    def save_qc_results(qc_results, out_file):
        out_path = Path(out_file)
        try:
            # Try to save as JSON (best for structured data)
            with out_path.open("w", encoding="utf-8") as f:
                json.dump(qc_results, f, indent=4)
            print(f"✅ QC results saved to JSON: {out_path}")
            return
        except TypeError:
            # JSON failed (object not serializable)
            pass
        # Fallback: save as plain text
        with out_path.open("w", encoding="utf-8") as f:
            f.write(str(qc_results))
    
    save_qc_results(qc_results=qc_results,out_file=f"{sample_column_dict['output_folder']}/{sample_column_dict['gwas_outputname']}_QC_sumamry.txt")
        
    
    try:
        postgwas_qc_df=qc_json_to_dataframes(qc_file=f"{sample_column_dict['output_folder']}/{sample_column_dict['gwas_outputname']}_QC_sumamry.txt")
    except:
        print("qc_json_to_dataframes function failed")



    # ---------------------------------------------------------
    # STEP 2 — Merge per-chromosome VCFs
    # ---------------------------------------------------------
    concat_vcfs_by_build(
        output_dir = sample_column_dict["output_folder"],
        gwas_outputname = sample_column_dict["gwas_outputname"],
        mode = "concurrent"
    )

    # ---------------------------------------------------------
    # STEP 3 — Cleanup temporary files
    # ---------------------------------------------------------
    clean_intermediate_files(
        output_dir = sample_column_dict["output_folder"],
        gwas_outputname = sample_column_dict["gwas_outputname"]
    )
    #print("\n🎉 Started variant summarisation.")
    outdir = Path(sample_column_dict["output_folder"])
    #print(outdir)
    #print(outdir / f"{sample_column_dict['gwas_outputname']}_GRCh37_merged.vcf.gz")


    raw_vcf_qc_df= run_qc_summary(
        vcf_path=outdir / f"{sample_column_dict['gwas_outputname']}_GRCh37_merged.vcf.gz",
        qc_outdir=outdir,
        sample_id=f"{sample_column_dict['gwas_outputname']}_GRCh37_raw",
        external_af_name = "EUR",
        allelefreq_diff_cutoff= 0.2,
        n_threads= nthreads,
        bcftools_bin= "bcftools"
    )
    raw_vcf_qc_df.columns=['index','raw_variant_count']
    qc_passed_vcf = filter_gwas_vcf_bcftools(
        vcf_path=outdir / f"{sample_column_dict['gwas_outputname']}_GRCh37_merged.vcf.gz",
        output_folder=str(outdir),
        output_prefix=f"{sample_column_dict['gwas_outputname']}_GRCh37_raw",
        pval_cutoff=None,
        maf_cutoff=None,
        allelefreq_diff_cutoff=0.2,
        info_cutoff=0.7,
        external_af_name=comparison_cols,
        include_indels=True,
        include_palindromic=True,
        palindromic_af_lower=0.4,
        palindromic_af_upper=0.6,
        remove_mhc=False,
        mhc_chrom="6",
        mhc_start=25000000,
        mhc_end=34000000,
        threads=nthreads,
        max_mem="12G"
    )
    print(qc_passed_vcf)
    qc_passed_vcf_df= run_qc_summary(
        vcf_path=outdir / f"{sample_column_dict['gwas_outputname']}_GRCh37_raw_filtered.vcf.gz",
        qc_outdir=outdir,
        sample_id=f"{sample_column_dict['gwas_outputname']}_GRCh37_filtered",
        external_af_name = "EUR",
        allelefreq_diff_cutoff= 0.2,
        n_threads= nthreads,
        bcftools_bin= "bcftools"
    )
    qc_passed_vcf_df.columns=['index','qc_passed_variant_count']
    qc_outdir = outdir/"qc_summary"
    qc_outdir.mkdir(parents=True, exist_ok=True)
    qc_file = qc_outdir / f"{sample_column_dict['gwas_outputname']}_final_qc_summary.csv"

    postgwas_qc_df.to_csv(f"{qc_outdir}/{sample_column_dict['gwas_outputname']}_inputfile_QC_summary.csv",index=None)
    final_qc_df =pd.merge(raw_vcf_qc_df,qc_passed_vcf_df,on="index",how="outer")
    #final_qc_df = pd.concat([final_qc_df,new_df], axis=0)
    # Ensure summary columns exist
    final_qc_df["total_variant_infile"] = pd.NA
    final_qc_df["total_variant_read"] = pd.NA
    # Create an empty row with correct length
    values = [pd.NA] * final_qc_df.shape[1]
    # Fill the "index" column
    values[final_qc_df.columns.get_loc("index")] = "variant_count"
    # Fill summary values
    values[final_qc_df.columns.get_loc("total_variant_infile")] = postgwas_qc_df.loc[0, "total_variant_infile"]
    values[final_qc_df.columns.get_loc("total_variant_read")] = postgwas_qc_df.loc[0, "total_variant_read"]
    # Append row
    final_qc_df.loc[len(final_qc_df)] = values
    final_qc_df.to_csv(qc_file, sep=",",index=None)
    os.system(f"rm {outdir}/{sample_column_dict['gwas_outputname']}_GRCh37_filtered*")
    os.system(f"rm {outdir}/{sample_column_dict['gwas_outputname']}_GRCh37_raw_filtered*")
    os.system(f"rm {outdir}/{sample_column_dict['gwas_outputname']}_GRCh37_raw_qc_summary.tsv")
    os.system(f"rm {outdir}/{sample_column_dict['gwas_outputname']}_gwas2vcf_summary.tsv")
    os.system(f"rm {outdir}/{sample_column_dict['gwas_outputname']}_*.stats")
    combine_logs_per_chromosome(f"{outdir}/logs/")
    return final_qc_df
