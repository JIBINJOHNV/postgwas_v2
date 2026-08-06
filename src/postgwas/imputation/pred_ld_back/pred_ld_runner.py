
import argparse
import os

from postgwas.imputation.pred_ld import pred_ld_functions


def main():
    parser = argparse.ArgumentParser(
        prog="postgwas impute pred-ld",
        description="Run PRED-LD summary statistics imputation and postprocessing (correlation, merging)."
    )

    # ─── Input arguments ────────────────────────────────────────────────────────────
    parser.add_argument("--sumstats", required=True,
                        help="Input GWAS summary statistics file in VCF format (full path).")
    parser.add_argument("--ref-ld", required=True,
                        help="Reference LD panel directory for PRED-LD.")
    parser.add_argument("--out-folder", required=True,
                        help="Output folder for imputation results.")
    parser.add_argument("--output-prefix", required=True,
                        help="Prefix for output files (e.g., PGC3_SCZ_GRCh38).")
    parser.add_argument("--gwas2vcf-resource", required=False, default=None,
                        help="Optional resource folder for downstream GWAS2VCF conversion.")

    # ─── Optional params ────────────────────────────────────────────────────────────
    parser.add_argument("--r2threshold", type=float, default=0.8,
                        help="LD r2 threshold for tagging (default: 0.8).")
    parser.add_argument("--maf", type=float, default=0.001,
                        help="Minor allele frequency threshold (default: 0.001).")
    parser.add_argument("--population", type=str, default="EUR",
                        help="Population code for LD reference (default: EUR).")
    parser.add_argument("--ref", type=str, default="TOP_LD",
                        help="Reference panel name (default: TOP_LD).")
    parser.add_argument("--threads", type=int, default=6,
                        help="Number of threads for parallel execution (default: 6).")
    parser.add_argument("--corr-method", type=str, default="pearson",
                        choices=["pearson", "spearman"],
                        help="Correlation method to use for per-chromosome QC (default: pearson).")
    parser.add_argument("--verbose", action="store_true",
                        help="Enable detailed logging.")

    args = parser.parse_args()
    
    pred_ld_functions.process_data(file_path, r2threshold, population, maf_input, ref_file)


    ref= "TOP_LD",
    process_data(file_path, 
                 r2threshold, 
                 population, 
                 maf_input, 
                 ref_file="TOP_LD")


sumstats='/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/PGC3_SCZ_european/6_format_inputs/ldpred/PGC3_SCZ_european_chr21_predld_input.tsv'