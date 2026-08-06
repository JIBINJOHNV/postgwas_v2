import os
import re
import sys
import json
import yaml
import argparse
import logging
import threading
import subprocess
import textwrap
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context

import polars as pl
import pandas as pd

# ----------------------------------------------------------------------
# Thread-safe print
# ----------------------------------------------------------------------
print_lock = threading.Lock()


def safe_print(*args, **kwargs):
    """Thread-safe / process-safe print helper."""
    with print_lock:
        print(*args, **kwargs, flush=True)


# ----------------------------------------------------------------------
# PostGWAS imports
# ----------------------------------------------------------------------
from postgwas.utils.main import validate_path, validate_prefix_files, apply_validator

from postgwas.harmonisation.io import (
    read_sumstats,
    read_config,
    find_resource_file_path,
    load_default_config,
    check_file_truncation
)
from postgwas.harmonisation.validator import validate_gwas_config
from postgwas.harmonisation.preprocess import split_chr_pos
from postgwas.harmonisation.chr_pos_process import fix_chr_pos_allele_column
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
from postgwas.harmonisation.calculate_se_from_beta_pvalue import (
    calculate_se_from_beta_pvalue,
)
from postgwas.harmonisation.add_or_calculate_info import add_or_calculate_info
from postgwas.harmonisation.add_or_fix_snp_column import add_or_fix_snp_column
from postgwas.harmonisation.sumstat_to_vcf import (
    run_bcftools_munge,
    run_bcftools_annot,
    concat_vcfs_by_build,
)
from postgwas.harmonisation.export_gwas_sumstat_forgwas2vcf import (
    export_gwas_sumstat,
)
from postgwas.harmonisation.gwastovcf_gwas2vcf import gwastovcf
from postgwas.harmonisation.clean_intermediate import clean_intermediate_files
from postgwas.harmonisation.postgwas_qc import qc_json_to_dataframes
from postgwas.harmonisation.utilities import combine_logs_per_chromosome

from postgwas.sumstat_filter.sumstat_filter import filter_gwas_vcf_bcftools
from postgwas.qc_summary.main import run_qc_summary
from postgwas.harmonisation.print_qc_summary import print_qc_summary


# ===============================================================
# Logging helpers
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
        formatter = logging.Formatter(
            "%(asctime)s - %(levelname)s - %(message)s"
        )
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
        target_fasta = (
            f"{resource_folder}/GRCh38/fasta_files/GRCh38_chr{chromosome}.fa"
        )
        chain_file = f"{resource_folder}/chain_files/GRCh37_to_GRCh38.chain"
    else:
        target_fasta = (
            f"{resource_folder}/GRCh37/fasta_files/GRCh37_chr{chromosome}.fa"
        )
        chain_file = f"{resource_folder}/chain_files/GRCh38_to_GRCh37.chain"
    # --- Validate core default files exist & are non-empty ---
    file_validator = validate_path(
        must_exist=True,
        must_be_file=True,
        must_not_be_empty=True,
    )
    errors = []
    for path in [
        default_eaf_path,
        default_info_path,
        default_comparison_af_path,
        dbsnp_path,
    ]:
        try:
            file_validator(path)
        except argparse.ArgumentTypeError as e:
            errors.append(f"❌ File invalid: {path}\n   → {e}")
    if errors:
        raise ValueError(
            "\n❌ One or more default files are invalid:\n" + "\n".join(errors)
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
    maf_eaf_decision_cutoff=0.95,
    extrnal_eaf_colmap={"chr": "CHROM","pos": "POS","a1": "ALT","a2": "REF"},

) -> Tuple[str, Dict[str, Any]]:
    """
    Run full harmonisation + GWAS→VCF + bcftools annotation for a single chromosome.
    Returns:
        (chromosome, qc_dict)
    """
    sample_id = sample_column_dict["gwas_outputname"]
    output_dir_path = Path(output_dir)
    invalide_info_path=output_dir_path / "invalid_info"
    invalide_info_path.mkdir(parents=True, exist_ok=True)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    log_dir = output_dir_path / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = _get_chr_logger(sample_id=sample_id, chromosome=chromosome, log_dir=log_dir)

    qc_dict: Dict[str, Any] = {}
    try:
        logger.info(f"🚀 Starting harmonization pipeline for chromosome {chromosome}")
        # ------------------------------
        # 1. Load per-chromosome file
        # ------------------------------
        #df = pl.read_csv(chr_file, separator="\t")
        df = pl.read_csv(
            chr_file, 
            separator="\t",
            null_values=["", "NA", "NaN", ".", "null"],
            infer_schema_length=10000
        )
        logger.info(f"Loaded {df.height:,} rows × {df.width} columns")
        df,sample_column_dict = fix_chr_pos_allele_column(
            chromosome=chromosome,
            df=df,
            sample_column_dict=sample_column_dict,
            drop_mt=True,
        )
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
            dbsnp=dbsnp
            )
        logger.info(f"Resource map constructed for {chromosome}")
        # ------------------------------
        # 3. Harmonisation steps
        # ------------------------------
        df, eaf_qc, sample_column_dict = add_or_calculate_eaf(
            chromosome=chromosome,
            df=df,
            sample_column_dict=sample_column_dict,
            eaffile=res["user_eaf_file"],
            default_eaf_file=res["default_eaf_file"],
            default_eaf_eafcolumn=res["default_eaf_column"],
            maf_eaf_decision_cutoff=maf_eaf_decision_cutoff,
            external_eaf_colmap=extrnal_eaf_colmap
        )
        logger.info("EAF harmonisation completed.")
        df, size_qc, sample_column_dict = harmonize_sample_sizes(
            chromosome=chromosome,
            df=df,
            sample_column_dict=sample_column_dict,
        )
        logger.info("Sample size harmonisation completed.")
        df, beta_qc, sample_column_dict = calculate_beta_and_se_from_z(
            chromosome=chromosome,
            df=df,
            sample_column_dict=sample_column_dict,
        )
        logger.info("BETA/SE from Z completed.")
        df, or_qc, sample_column_dict = is_beta_or_or(
            chromosome=chromosome,
            df=df,
            sample_column_dict=sample_column_dict,
        )
        logger.info("Effect-type (BETA/OR) harmonisation completed.")
        df, pval_qc, sample_column_dict = detect_and_convert_pval(
            chromosome=chromosome,
            df=df,
            sample_column_dict=sample_column_dict,
        )
        logger.info("P-value detection and conversion completed.")
        df, se_qc, sample_column_dict = calculate_se_from_beta_pvalue(
            chromosome=chromosome,
            df=df,
            sample_column_dict=sample_column_dict,
        )
        logger.info("SE from BETA/P-value calculation completed.")
        df, z_qc, sample_column_dict = calculate_z_from_beta_se(
            chromosome=chromosome,
            df=df,
            sample_column_dict=sample_column_dict,
        )
        logger.info("Z from BETA/SE calculation completed.")
        df, info_qc, sample_column_dict = add_or_calculate_info(
            chromosome=chromosome,
            df=df,
            sample_column_dict=sample_column_dict,
            info_file=res["user_info_file"],              
            info_column=res["user_info_column"],
            default_info_file=res["default_info_file"],
            default_info_column=res["default_info_column"],
            output_folder=str(invalide_info_path),
            output_prefix=sample_column_dict["gwas_outputname"],
            log_dir=str(log_dir),
        )
        logger.info("INFO score harmonisation completed.")
        df, sample_column_dict = add_or_fix_snp_column(
            chromosome=chromosome,
            df=df,
            sample_column_dict=sample_column_dict,
        )
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
        # 7. QC summary (per-chromosome)
        # ------------------------------
        qc_dict = {
            "eaf_qc": eaf_qc,
            "sample_size_qc": size_qc,
            "beta_and_se_from_zscore_qc": beta_qc, 
            "beta_or_oddsratio_qc": or_qc,
            "pval_detection_and_conversion_qc": pval_qc,
            "se_from_beta_pvalue_qc": se_qc,
            "z_from_beta_se_qc": z_qc,
            "info_qc": info_qc,
            "gwastovcf_exit_code": gwastovcf_exit_code,
        }
        logger.info(f"      ✅ Completed chromosome {chromosome}")
    except Exception as e:
        logger.exception(f"❌ ERROR in chromosome {chromosome}: {e}")
        qc_dict = {"error": str(e)}
    return chromosome, qc_dict


# ----------------------------------------------------------------------
# Helper: Chromosome extractor
# ----------------------------------------------------------------------
def extract_chromosome_from_filename(filepath: str) -> str:
    filename = Path(filepath).name 
    # Exact match for "_chr" (case-sensitive)
    m = re.search(
        r"_chr(\d+|X|Y|MT)\.tsv", 
        filename
    )
    if not m:
        raise ValueError(
            f"Cannot extract chromosome from filename: {filepath}"
        ) 
    extracted_chr = m.group(1).upper()
    return extracted_chr

# ===============================================================
#  MULTIPROCESSING VERSION — Python 3.8 / Docker SAFE
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
    maf_eaf_decision_cutoff=0.95,
    extrnal_eaf_colmap={"chr": "CHROM","pos": "POS","a1": "ALT","a2": "REF"}
):
    """
    Full per-chromosome parallel harmonisation + GWAS2VCF.

    External behaviour is retained; internal implementation is made more robust
    for Linux/Docker (Python 3.8) using spawn-safe ProcessPoolExecutor and
    better logging/progress reporting.
    """
    safe_print(f"          🔹 Starting harmonisation + GWAS-to-VCF pipeline with max_workers={max_workers} ")
    safe_print(f"            for {sumstat_file}")
    # ------------------------------
    # Create output path
    # ------------------------------
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    # ------------------------------
    # Read file + validate config
    # ------------------------------
    df, file_cvariant_count, polars_rows,sample_column_dict,removed_coords, non_standard_allele_count,removed_duplicates_count,removed_missing_count = read_sumstats(
        sumstat_file=sumstat_file,
        output_dir=str(output_dir_path),
        sample_column_dict=sample_column_dict  
    )
    safe_print("")
    safe_print(f"           Total variants in input file                                          : {file_cvariant_count:,}")
    safe_print(f"           Total variants read by the python                                     : {polars_rows:,}")
    safe_print(f"           Total variants removed due to missing values in required columns      : {removed_missing_count:,}")
    safe_print(f"           Total variants removed due to NULL coordinates or invalid chromosomes : {removed_coords:,}")
    safe_print(f"           Total variants removed due to non-standard alleles                    : {non_standard_allele_count:,}")
    safe_print(f"           Total variants removed due to duplicate variants[chr,pos,ea,oa]       : {removed_duplicates_count:,}")
    safe_print(f"           Total variants remaining for harmonization                            : {df.height:,}")
    safe_print(" ")
    
    # ------------------------------
    # QC container
    # ------------------------------
    per_chr_qc: Dict[str, Any] = {}
    # 1. Call the function ONCE and unpack the results
    is_valid, error_list = validate_gwas_config(sample_column_dict, df)
    # 2. Check the boolean flag
    if not is_valid:
        print(f"\nExiting due to {len(error_list)} validation error(s):")
        for err in error_list:
            print(f" - {err}")
        sys.exit(1)
    # Genome build inference
    # ------------------------------
    genome_build_info = genome_build(
        df,
        grch37_file,
        grch38_file,
        sample_column_dict,
    )
    grch_version = genome_build_info["inferred_build"]
    safe_print(" ")
    safe_print(
        f"\033[96m           🧬 Inferred genome build:\033[0m "
        f"\033[93m{grch_version}\033[0m"
    )
    safe_print(
    f"                  🧬 Build: {genome_build_info['inferred_build']} | "
    f"Total variants: {genome_build_info['total_variants']:,}" )
    safe_print(
        f"                  🧬 GRCh37 matches: {genome_build_info['grch37_matches']:,} "
        f"({genome_build_info['percent_grch37']}%)  |  "
        f"GRCh38 matches: {genome_build_info['grch38_matches']:,} "
        f"({genome_build_info['percent_grch38']}%)"
    )
    safe_print(" ")
    if grch_version == "Ambiguous":
        raise RuntimeError(
            "❌ Genome build inference failed (Ambiguous).\n"
            "   Please specify the genome build explicitly using:\n"
            "     --genome-build GRCh37  or  --genome-build GRCh38\n"
            "   Aborting pipeline to avoid invalid gwas2vcf resource usage."
        )
    # ------------------------------
    # Split by chromosome
    # ------------------------------
    chr_files = split_chr_pos(df, sample_gwas_dict=sample_column_dict)
    final_chr_files: List[str] = []
    for f in chr_files:
        try:
            _ = extract_chromosome_from_filename(f)
            final_chr_files.append(f)
        except Exception:
            # skip files without chr pattern
            continue
    n_chr = len(final_chr_files) 
    max_workers = min(n_chr, max_workers) 
    safe_print(
        f"          🔹 [INFO] Processing {n_chr} chromosomes ({final_chr_files}"
    )
    safe_print(
        f"          🔹in parallel using {max_workers} workers" )
    safe_print(" ")
    if n_chr == 0:
        raise RuntimeError(
            "No per-chromosome files found after splitting. "
            "Check that 'CHR' and 'POS' columns are correctly mapped."
        )
    # ------------------------------
    # Run each chromosome in multiprocessing
    # ------------------------------
    # Use spawn context for better safety on Linux / Docker
    mp_ctx = get_context("spawn")
    errors_seen: List[str] = []
    try:
        with ProcessPoolExecutor(
            max_workers=max_workers,
            mp_context=mp_ctx,
        ) as executor:
            future_to_chr: Dict[Any, str] = {}
            for chr_file in final_chr_files:
                chrom = extract_chromosome_from_filename(chr_file)
                # clone dict so each process works on its own copy
                scd = sample_column_dict.copy()
                future = executor.submit(
                    process_one_chromosome,
                    chrom,
                    chr_file,
                    scd,
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
                    5,  # threads per chromosome for bcftools etc.
                    gwastovcf_main_script_path,
                    maf_eaf_decision_cutoff=maf_eaf_decision_cutoff,
                    extrnal_eaf_colmap=extrnal_eaf_colmap,
                )
                future_to_chr[future] = chrom
            # --------------------------
            # Collect results safely
            # --------------------------
            completed = 0
            total = len(future_to_chr)
            for future in as_completed(future_to_chr):
                chrom = future_to_chr[future]
                completed += 1
                try:
                    chrom_res, qc = future.result()
                    per_chr_qc[chrom_res] = qc
                    if "error" in qc:
                        msg = f"❌ Chromosome {chrom_res} completed with ERROR."
                        safe_print(f"{msg} See logs/{sample_column_dict['gwas_outputname']}_chr{chrom_res}.log")
                        errors_seen.append(msg + " " + qc.get("error", ""))
                    else:
                        safe_print(f"                ✅ [{completed}/{total}] Completed chromosome {chrom_res}")
                except Exception as exc:
                    msg = f"❌ [ERROR] Chromosome {chrom} crashed: {exc}"
                    safe_print(msg)
                    per_chr_qc[chrom] = {"error": str(exc)}
                    errors_seen.append(msg) 
    finally:
        per_chr_qc["total_variant_infile"] = file_cvariant_count
        per_chr_qc["total_variant_read"] = polars_rows
        per_chr_qc["total_variant_removed_null_coords"] = removed_coords
        per_chr_qc["total_variant_removed_non_standard_alleles"] = non_standard_allele_count
        per_chr_qc["total_variant_remaining_for_harmonisation"] = df.height
        per_chr_qc["total_variant_removed_missing_values"] = removed_missing_count
        if errors_seen:
            safe_print("\n          ⚠️ Some chromosomes finished with errors:")
            for e in errors_seen:
                safe_print("   -", e)
        else:
            safe_print("\n          🎉All chromosomes completed without errors.")
    return per_chr_qc,grch_version


# ===============================================================
#  High-level pipeline wrapper
# ===============================================================

def _save_qc_results(qc_results: Dict[str, Any], out_file: Path) -> None:
    """
    Save QC dictionary to disk.
    Tries JSON first, falls back to plain-text representation.
    """
    out_file = Path(out_file)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    try:
        with out_file.open("w", encoding="utf-8") as f:
            json.dump(qc_results, f, indent=4)
        #safe_print(f"✅ QC results saved to JSON: {out_file}")
        return
    except TypeError:
        # JSON failed (object not serializable)
        pass
    # Fallback: save as plain text
    with out_file.open("w", encoding="utf-8") as f:
        f.write(str(qc_results))
    #safe_print(f"✅ QC results saved as plain text: {out_file}")


def run_harmonisation_pipeline(
    sample_column_dict: Dict[str, Any],
    default_cfg: str,
    nthreads: int,
):
    """
    Full harmonisation + gwas2vcf + merge + cleanup pipeline.
    Clean and modular version to be called from cli.py.
    External behaviour (arguments, outputs, filenames) is retained.
    Implementation is made more robust for Docker / Linux.
    """
    BOLD = "\033[1m"
    RESET = "\033[0m"
    
    # -----------------------------------------------------
    # Load DEFAULT config
    # -----------------------------------------------------
    try:
        default_cfg_obj = load_default_config(default_cfg)
        if not default_cfg_obj:
            raise ValueError("load_default_config returned empty output.")
    except Exception as e:
        safe_print(f"❌ ERROR: Failed to load default config file '{default_cfg}'.")
        safe_print(f"   Reason: {e}")
        sys.exit(1)
    # Extract all defaults safely
    default_eaf = default_cfg_obj["default_eaf"]
    default_info = default_cfg_obj["default_info"]
    default_comparison_af = default_cfg_obj["default_comparison_af_file"]
    default_eaf_column = default_cfg_obj["default_eaf_column"]
    default_info_column = default_cfg_obj["default_info_column"]
    comparison_cols = default_cfg_obj["default_comparison_af_column"]
    dbsnp = default_cfg_obj["default_dbsnp"]
    required_cols = default_cfg_obj["required_columns_in_sumstat"]
    optional_cols = default_cfg_obj["optional_columns_in_sumstat"]
    #gwastovcf_script = default_cfg_obj["gwastovcf_main_script_path"]
    maf_eaf_decision_cutoff =default_cfg_obj['maf_eaf_decision_cutoff']
    extrnal_eaf_colmap=default_cfg_obj["extrnal_eaf_colmap"]
    # Build clean path using pathlib
    output_folder = Path(sample_column_dict["output_folder"])
    gwastovcf_script = str(Path(__file__).parent / "gwas2vcf" / "main.py")
    print("")
    safe_print(f"            To compare GWAS Alternative allele frequency with External allele freqency use : {comparison_cols} column")
    safe_print(f"            present in the dataset {default_comparison_af}")
    safe_print("")
    # Save normalized output_folder  back to dict
    output_folder = output_folder / "00_harmonised_sumstat"
    sample_column_dict["output_folder"] = str(output_folder)
    # Clean whitespace from keys AND values
    sample_column_dict = {
        k.strip(): (v.strip() if isinstance(v, str) else v) 
        for k, v in sample_column_dict.items()
    }
    # Try to create it safely
    try:
        output_folder.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        safe_print("\n❌ ERROR: Failed to create output directory.")
        safe_print(f"   Path: {output_folder}")
        safe_print(f"   Reason: {e}\n")
        sys.exit(1)
    # ---------------------------------------------------------
    # Validate resource folder
    # ---------------------------------------------------------
    resource_folder = Path(sample_column_dict["resourse_folder"])
    if not resource_folder.exists():
        raise ValueError(
            f"❌ ERROR: Resource folder does not exist:\n   {resource_folder}"
        )
    if not resource_folder.is_dir():
        raise ValueError(
            f"❌ ERROR: Resource folder is not a directory:\n   {resource_folder}"
        )
    if not any(resource_folder.iterdir()):
        raise ValueError(
            f"❌ ERROR: Resource folder is empty:\n   {resource_folder}"
        )
    # Validate sumstats + GRCh37/38 check files
    file_validator = validate_path(
        must_exist=True,
        must_be_file=True,
        must_not_be_empty=True,
    )
    errors = []
    grch37_check = (
        f"{sample_column_dict['resourse_folder']}/GRCh37_38_check_files/GRCh37_check_file.tsv"
    )
    grch38_check = (
        f"{sample_column_dict['resourse_folder']}/GRCh37_38_check_files/GRCh38_check_file.tsv"
    )

    for path in [
        sample_column_dict["sumstat_file"],
        grch37_check,
        grch38_check,
    ]:
        try:
            file_validator(path)
        except argparse.ArgumentTypeError as e:
            errors.append(f"❌ File invalid: {path}\n   → {e}")
    if errors:
        raise ValueError(
            "\n❌ One or more required files are invalid:\n" + "\n".join(errors)
        )
    # --------------------------------------------------------
    # STEP 0.5 — Check GWAS file truncation OR EXISTS
    # ---------------------------------------------------------
    check_file_truncation(
        file_path=sample_column_dict["sumstat_file"]
    )
    # ---------------------------------------------------------
    # STEP 1 — FULL harmonisation + per-chromosome GWAS2VCF
    # ---------------------------------------------------------
    qc_results,grch_version = gwas_to_vcf_parallel(
        sumstat_file=sample_column_dict["sumstat_file"],
        sample_column_dict=sample_column_dict,
        output_dir=sample_column_dict["output_folder"],
        resource_folder=sample_column_dict["resourse_folder"],
        default_eaf_file=default_eaf,
        default_info_file=default_info,
        default_comparison_af_file=default_comparison_af,
        dbsnp=dbsnp,
        required_columns_in_sumstat=required_cols,
        optional_columns_in_sumstat=optional_cols,
        max_workers=nthreads,
        user_eaf_file=sample_column_dict["eaffile"],
        user_eaf_column=sample_column_dict["eafcolumn"],
        user_info_file=sample_column_dict["infofile"],
        user_info_column=sample_column_dict["infocolumn"],
        default_eaf_column=default_eaf_column,
        default_comparison_af_column=comparison_cols,
        default_info_column=default_info_column,
        grch37_file=grch37_check,
        grch38_file=grch38_check,
        gwastovcf_main_script_path=gwastovcf_script,
        maf_eaf_decision_cutoff=maf_eaf_decision_cutoff,
        extrnal_eaf_colmap=extrnal_eaf_colmap
    )
    qc_summary_path = (
        Path(sample_column_dict["output_folder"])
        / f"{sample_column_dict['gwas_outputname']}_QC_summary.txt"
    )
    _save_qc_results(qc_results=qc_results, out_file=qc_summary_path)
    try:
        postgwas_qc_df2 = qc_json_to_dataframes(
            json_path=str(qc_summary_path)
        )
    except Exception as e:
        safe_print("⚠️ qc_json_to_dataframes function failed; continuing without detailed input QC.")
        safe_print(f"⚠️ {e}")
        # fabricate minimal info to not break downstream summary
        postgwas_qc_df2 = pd.DataFrame(
            [
                {
                    "total_variant_infile": qc_results.get(
                        "total_variant_infile", pd.NA
                    ),
                    "total_variant_read": qc_results.get(
                        "total_variant_read", pd.NA
                    ),
                    "total_variant_removed_null_coords": qc_results.get(
                        "total_variant_removed_null_coords", pd.NA
                    ),
                    "total_variant_removed_non_standard_alleles": qc_results.get(
                        "total_variant_removed_non_standard_alleles", pd.NA
                    ),
                    "total_variant_remaining_for_harmonisation": qc_results.get(
                        "total_variant_remaining_for_harmonisation", pd.NA
                    ),
                }
            ]
        )
    

    # ---------------------------------------------------------
    # STEP 2 — Merge per-chromosome VCFs
    # ---------------------------------------------------------
    #safe_print("\n📦 Concatenating per-chromosome VCFs...")
    concated_vcf_files = concat_vcfs_by_build(
        output_dir=sample_column_dict["output_folder"],
        gwas_outputname=sample_column_dict["gwas_outputname"],
        grch_version=grch_version
    )

    # ---------------------------------------------------------
    # STEP 3 — Cleanup temporary files
    # ---------------------------------------------------------
    #safe_print("🧹 Cleaning intermediate files...")
    clean_intermediate_files(
        output_dir=sample_column_dict["output_folder"],
        gwas_outputname=sample_column_dict["gwas_outputname"],
    )

    outdir = Path(sample_column_dict["output_folder"])

    # ---------------------------------------------------------
    # STEP 4 — Raw VCF QC
    # ---------------------------------------------------------
    raw_vcf_path = (
        outdir
        / f"{sample_column_dict['gwas_outputname']}_GRCh37_merged.vcf.gz"
    )
    #safe_print(f"\n📊 Running QC summary on raw VCF: {raw_vcf_path}")

    raw_vcf_qc_df = run_qc_summary(
        vcf_path=raw_vcf_path,
        qc_outdir=outdir,
        sample_id=f"{sample_column_dict['gwas_outputname']}_GRCh37_raw",
        external_af_name=comparison_cols,
        allelefreq_diff_cutoff=0.2,
        n_threads=nthreads,
        bcftools_bin="bcftools",
    )
    raw_vcf_qc_df.columns = ["index", "raw_variant_count"]

    # ---------------------------------------------------------
    # STEP 5 — Filtering VCF
    # ---------------------------------------------------------
    #safe_print("🧪 Filtering VCF based on INFO / AF / MHC settings...") 
    gwas_to_vcf_qc_df = pd.read_csv( Path(outdir) / "qc_summary" / f"{sample_column_dict['gwas_outputname']}_gwas2vcf_summary.tsv",sep="\t")
    gwas_to_vcf_qc_df.columns=[x.strip() for x in gwas_to_vcf_qc_df.columns]
    gwas_to_vcf_qc_df = gwas_to_vcf_qc_df[  gwas_to_vcf_qc_df['chromosome'] != "chromosome" ]
    snp_id_df = gwas_to_vcf_qc_df[  gwas_to_vcf_qc_df['key'] == "snp_id_col" ]
    total_variants_in_gwastovcf = int(snp_id_df['num_rows'].astype(int).sum())
    gwas_to_vcf_qc_df.to_csv(Path(outdir) / "qc_summary" / f"{sample_column_dict['gwas_outputname']}_gwas2vcf_summary.tsv",sep="\t",index=None)
    postgwas_qc_df2.to_csv(Path(sample_column_dict["output_folder"]) / "qc_summary" / f"{sample_column_dict['gwas_outputname']}_gwas2vcf_summary2.tsv",sep="\t",index=None)
    try:
        missing_eaf_count = (
            postgwas_qc_df2.loc[postgwas_qc_df2['metric'] == 'missing_eaf_count']
            .iloc[:, 2:]
            .apply(pd.to_numeric, errors='coerce')
            .sum(axis=1)
            .iloc[0]  
        )
        invalid_beta_se_count = (
            postgwas_qc_df2.loc[postgwas_qc_df2['metric'] == 'variants_removed_invalid_beta_se']
            .iloc[:, 2:]
            .apply(pd.to_numeric, errors='coerce')
            .sum(axis=1)
            .iloc[0]  
        ) 
    except:
        missing_eaf_count = None
        invalid_beta_se_count = None

    indent = "\t" * 3
    print(
        f"\n{indent}{BOLD}🔎 NOTE:{RESET}\n"
        f"{indent}{BOLD}Variant filtering is performed only to assess the total number of high-quality variants.{RESET}\n"
    )    
    variant_qc_summary={
                    "total_variant_infile": qc_results.get( "total_variant_infile", pd.NA),
                    "total_variant_read": qc_results.get( "total_variant_read", pd.NA ),
                    "total_variant_in_vcf_input":total_variants_in_gwastovcf,
                    "total_variant_removed_null_coords": qc_results.get(
                        "total_variant_removed_null_coords", pd.NA
                    ),
                    "total_variant_removed_non_standard_alleles": qc_results.get(
                        "total_variant_removed_non_standard_alleles", pd.NA
                    ),
                    "total_variant_remaining_for_harmonisation": qc_results.get(
                        "total_variant_remaining_for_harmonisation", pd.NA
                    ),
                    "total_variant_with_missing_eaf": int(missing_eaf_count),
                    "total_variant_with_invalid_beta_se": int(invalid_beta_se_count)
                    }
    
    qc_passed_vcf = filter_gwas_vcf_bcftools(
        vcf_path=raw_vcf_path,
        output_folder=str(outdir),
        output_prefix=f"{sample_column_dict['gwas_outputname']}", 
        pval_cutoff=None,
        maf_cutoff=0.01,
        allelefreq_diff_cutoff=0.2,
        info_cutoff=0.7,
        info_max=1.05,
        info_missing="remove",
        external_af_name=comparison_cols,
        include_indels=False,
        exclude_palindromic=True,
        palindromic_af_lower=0.4,
        palindromic_af_upper=0.6,
        remove_mhc=True,
        mhc_chrom="6",
        mhc_start=25000000,
        mhc_end=34000000,
        threads=nthreads,
        max_mem="12G",
        extra_options=variant_qc_summary
    )
    #safe_print(f"✅ Filtered VCF written to: {qc_passed_vcf}")

    filtered_vcf_path = qc_passed_vcf['filtered_vcf']


    qc_passed_vcf_df = run_qc_summary(
        vcf_path=filtered_vcf_path,
        qc_outdir=outdir,
        sample_id=f"{sample_column_dict['gwas_outputname']}_GRCh37_filtered",
        external_af_name=comparison_cols,
        allelefreq_diff_cutoff=0.2,
        n_threads=nthreads,
        bcftools_bin="bcftools",
    )
    qc_passed_vcf_df.columns = ["index", "qc_passed_variant_count"]

    # ---------------------------------------------------------
    # STEP 6 — Final QC summary
    # ---------------------------------------------------------
    qc_outdir = outdir / "qc_summary"
    qc_outdir.mkdir(parents=True, exist_ok=True)
    qc_file = (
        qc_outdir
        / f"{sample_column_dict['gwas_outputname']}_final_qc_summary.csv"
    )


    final_qc_df = pd.merge(
        raw_vcf_qc_df, qc_passed_vcf_df, on="index", how="outer"
    )  

    final_qc_df.loc[final_qc_df["index"] == "num_records", "total_variant_infile"] = variant_qc_summary["total_variant_infile"]
    final_qc_df.loc[final_qc_df["index"] == "num_records", "total_variant_read"] = variant_qc_summary["total_variant_read"]
    final_qc_df.loc[final_qc_df["index"] == "num_records", "total_variant_input_to_gwas2vcf"] = variant_qc_summary["total_variant_in_vcf_input"]
    final_qc_df.to_csv(qc_file, sep=",", index=None)

    # ---------------------------------------------------------
    # STEP 7 — Cleanup VCF QC intermediates
    # ---------------------------------------------------------
    #safe_print("🧹 Cleaning QC intermediate files...")
    os.system(
        f"rm -f {outdir}/{sample_column_dict['gwas_outputname']}_GRCh37_filtered*"
    )
    os.system(
        f"rm -f {outdir}/{sample_column_dict['gwas_outputname']}_GRCh37_raw_filtered*"
    )
    os.system(
        f"rm -f {outdir}/{sample_column_dict['gwas_outputname']}_GRCh37_raw_qc_summary.tsv"
    )
    os.system(
        f"rm -f {outdir}/{sample_column_dict['gwas_outputname']}_gwas2vcf_summary.tsv"
    )
    os.system(
        f"rm -f {outdir}/{sample_column_dict['gwas_outputname']}_*.stats"
    )

    # Combine logs into a single summary log 
    PREFIX = "\t\t"
    combine_logs_per_chromosome(f"{outdir}/logs/")

    print_qc_summary(raw_vcf_qc_df, columns="raw_variant_count") 
    print(f"{PREFIX}🧬 {'Total variant with missing EAF':45s}: {int(missing_eaf_count):,}")
    print_qc_summary(qc_passed_vcf_df, columns="qc_passed_variant_count")

    safe_print(f"\n          🎉 Harmonisation + GWAS2VCF pipeline completed successfully for sample: {sample_column_dict['gwas_outputname']}.\n")
    print(" ")
    return {
        "GRCh37":f"{outdir}/{sample_column_dict['gwas_outputname']}_GRCh37_merged.vcf.gz",
        "GRCh38":f"{outdir}/{sample_column_dict['gwas_outputname']}_GRCh38_merged.vcf.gz",
        "gwas2vcf":f"{outdir}/{sample_column_dict['gwas_outputname']}_gwas2vcf_{grch_version}_merged.vcf.gz"
    }