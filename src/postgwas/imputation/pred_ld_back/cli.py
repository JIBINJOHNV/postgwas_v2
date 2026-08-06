import argparse
import os
from .pred_ld_runner_back import run_pred_ld_parallel
from .pred_ld_result_processor import process_pred_ld_results_all_parallel

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

    # ─── Step 1: Run PRED-LD imputation ─────────────────────────────────────────────
    print("🚀 Starting PRED-LD imputation pipeline...")
    run_pred_ld_parallel(
        sumstat_vcf=args.sumstats,
        output_folder=args.out_folder,
        output_prefix=args.output_prefix,
        pred_ld_ref=args.ref_ld,
        r2threshold=args.r2threshold,
        maf=args.maf,
        population=args.population,
        ref=args.ref,
        threads=args.threads
    )

    # ─── Step 2: Process and summarize results ──────────────────────────────────────
    print("\n📊 Running post-processing (correlations + merging)...")
    combined_df, corr_df = process_pred_ld_results_all_parallel(
        folder_path=args.out_folder,
        output_prefix=args.output_prefix,
        corr_method=args.corr_method,
        threads=args.threads
    )

    # ─── Optional: save GWAS2VCF spec ───────────────────────────────────────────────
    if args.gwas2vcf_resource:
        print("\n🧩 Preparing GWAS2VCF specification file...")
        spec_path = os.path.join(args.out_folder, f"{args.output_prefix}_gwas2vcf_spec.tsv")
        # Minimal spec table (if helper exists)
        combined_df.to_csv(spec_path, sep="\t", index=False)
        print(f"💾 Saved GWAS2VCF spec file → {spec_path}")

    print("\n✅ All steps completed successfully!")

if __name__ == "__main__":
    main()





# sumstat_vcf="/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/PGC3_SCZ_european/PGC3_SCZ_european_GRCh38_merged_filtered.vcf.gz"
# output_folder="/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/PGC3_SCZ_european/"
# output_prefix="PGC3_SCZ_european_GRCh38"
# pred_ld_reff="/Users/JJOHN41/Documents/software_resources/resourses/postgwas/imputation/pred-ld/ref/"
# gwas2vcf_resource_folder="/Users/JJOHN41/Documents/software_resources/resourses/gwas2vcf/"

# run_pred_ld_parallel(
#     sumstat_vcf=sumstat_vcf,
#     output_folder=output_folder,
#     output_prefix=output_prefix,
#     pred_ld_ref=pred_ld_reff,
#     r2threshold= 0.8,
#     maf= 0.001,
#     population= "EUR",
#     ref= "TOP_LD",
#     threads= 3  # number of parallel jobs
# )

# process_pred_ld_results_all_parallel(
#     folder_path=output_folder,
#     output_path=output_folder,
#     output_prefix=output_prefix,
#     corr_method= "pearson",
#     gwas2vcf_resource_folder=gwas2vcf_resource_folder
# )