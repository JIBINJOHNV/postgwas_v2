import os
import pandas as pd
from pathlib import Path

from postgwas.qc_summary.qc_summary import (
    run_bcftools_stats,
    parse_bcftools_stats,
    bcftools_essential_summary
)

from postgwas.sumstat_filter.sumstat_filter import filter_gwas_vcf_bcftools




import pandas as pd
from pathlib import Path  # <--- REQUIRED IMPORT

def run_qc_summary(
    vcf_path: str,
    qc_outdir: str,
    sample_id: str,
    external_af_name: str = "EUR",
    allelefreq_diff_cutoff: float = 0.2,
    n_threads: int = 5,
    bcftools_bin: str = "bcftools",
):
    """
    Pipeline wrapper that runs:
        1. run_bcftools_stats
        2. parse_bcftools_stats
        3. bcftools_essential_summary
        4. Construct final one-column summary table

    Returns:
        pandas.DataFrame — final QC table (index = metrics)
    """

    # ------------------------------------------------------------------
    # 1️⃣ Run bcftools QC stats
    # ------------------------------------------------------------------
    raw_data_qc1 = run_bcftools_stats(
        vcf_path=vcf_path,
        external_af_name=external_af_name,
        allelefreq_diff_cutoff=allelefreq_diff_cutoff,
        threads=n_threads,
        bcftools_bin=bcftools_bin
    )
    print(" ") 
    print("\t\t\tpostgwas qc summary module started running")
    print(" ")
    
    # ------------------------------------------------------------------
    # 2️⃣ Parse bcftools stats into sectioned DataFrames
    # ------------------------------------------------------------------
    stats_file = f"{vcf_path}.stats"
    data_dfs = parse_bcftools_stats(stats_file)

    # ------------------------------------------------------------------
    # 3️⃣ Extract essential summary (your original function)
    # ------------------------------------------------------------------
    data_df = bcftools_essential_summary(data_dfs)

    # ------------------------------------------------------------------
    # 4️⃣ Reformat into expected table (transpose → one column)
    # ------------------------------------------------------------------
    data_df = data_df.T
    data_df.columns = ["raw_variants"]

    # ------------------------------------------------------------------
    # 5️⃣ Append two additional metrics
    # ------------------------------------------------------------------
    data_df.loc["variants_with_missing_extrnal_af"] = raw_data_qc1["missing_extrnal_af"]
    data_df.loc["variants_with_EAF_diff_gt_cutoff"] = raw_data_qc1["af_diff_failed"]
    data_df = data_df.reset_index()

    # ------------------------------------------------------------------
    # 6️⃣ Save table (FIXED PATH HANDLING)
    # ------------------------------------------------------------------
    qc_outdir = Path(qc_outdir)  # <--- FIX: Convert string to Path object
    qc_outdir.mkdir(parents=True, exist_ok=True)
    
    qc_file = qc_outdir / f"{sample_id}_qc_summary.tsv"
    data_df.to_csv(qc_file, sep="\t", index=None)

    return data_df











# def run_full_vcf_qc_and_filter(
#     vcf_path: str,
#     output_folder: str,
#     output_prefix: str,
#     pval_cutoff=None,
#     maf_cutoff=None,
#     allelefreq_diff_cutoff=0.2,
#     info_cutoff=0.7,
#     external_af_name="EUR",
#     include_indels=True,
#     include_palindromic=True,
#     palindromic_af_lower=0.4,
#     palindromic_af_upper=0.6,
#     remove_mhc=False,
#     mhc_chrom="6",
#     mhc_start=25000000,
#     mhc_end=34000000,
#     threads=5,
#     max_mem="5G"
# ):
#     """
#     Run:
#       - Raw bcftools QC
#       - Filtering pipeline
#       - Post-filter bcftools QC
#       - Produce final summary DataFrame

#     Returns:
#         final_qc_df (pandas.DataFrame)
#         filtered_vcf_path (str)
#     """
#     output_folder = Path(output_folder)
#     output_folder.mkdir(parents=True, exist_ok=True)
#     # ------------------------------
#     # 1️⃣ RAW QC
#     # ------------------------------
#     raw_stats = run_bcftools_stats(
#         vcf_path=vcf_path,
#         external_af_name=external_af_name,
#         allelefreq_diff_cutoff=allelefreq_diff_cutoff,
#         threads=threads
#     )
#     raw_dfs = parse_bcftools_stats(stats_file=f"{vcf_path}.stats")
#     raw_df = bcftools_essential_summary(raw_dfs).T
#     raw_df.columns = ["raw_variants"]
#     raw_df.loc["variants_with_missing_extrnal_af"] = raw_stats["missing_extrnal_af"]
#     raw_df.loc["variants_with_EAF_diff_gt_0.2"] = raw_stats["af_diff_failed"]
#     # ------------------------------
#     # 2️⃣ FILTERING
#     # ------------------------------
#     filtered_vcf = filter_gwas_vcf_bcftools(
#         vcf_path=vcf_path,
#         output_folder=str(output_folder),
#         output_prefix=output_prefix,
#         pval_cutoff=pval_cutoff,
#         maf_cutoff=maf_cutoff,
#         allelefreq_diff_cutoff=allelefreq_diff_cutoff,
#         info_cutoff=info_cutoff,
#         external_af_name=external_af_name,
#         include_indels=include_indels,
#         include_palindromic=include_palindromic,
#         palindromic_af_lower=palindromic_af_lower,
#         palindromic_af_upper=palindromic_af_upper,
#         remove_mhc=remove_mhc,
#         mhc_chrom=mhc_chrom,
#         mhc_start=mhc_start,
#         mhc_end=mhc_end,
#         threads=threads,
#         max_mem=max_mem
#     )
#     # ------------------------------
#     # 3️⃣ POST-FILTER QC
#     # ------------------------------
#     filt_stats = run_bcftools_stats(
#         vcf_path=filtered_vcf,
#         external_af_name=external_af_name,
#         allelefreq_diff_cutoff=allelefreq_diff_cutoff,
#         threads=threads
#     )
#     filt_dfs = parse_bcftools_stats(stats_file=f"{filtered_vcf}.stats")
#     filt_df = bcftools_essential_summary(filt_dfs).T
#     filt_df.columns = ["qc_passed_variants"]
#     filt_df.loc["variants_with_missing_extrnal_af"] = filt_stats["missing_extrnal_af"]
#     filt_df.loc["variants_with_EAF_diff_gt_0.2"] = filt_stats["af_diff_failed"]
#     # ------------------------------
#     # 4️⃣ BUILD FINAL SUMMARY TABLE
#     # ------------------------------
#     final_qc_df = pd.concat([raw_df, filt_df], axis=1)
#     # Save table
#     qc_outdir = output_folder / "qc_summary"
#     qc_outdir.mkdir(exist_ok=True)
#     qc_file = qc_outdir / f"{output_prefix}_final_qc_summary.tsv"
#     final_qc_df.to_csv(qc_file, sep="\t")
#     p = Path(filtered_vcf)
#     # delete the VCF
#     p.unlink(missing_ok=True)
#     # delete index if exists
#     idx1 = p.with_suffix(p.suffix + ".tbi")
#     idx1.unlink(missing_ok=True)
#     return final_qc_df
