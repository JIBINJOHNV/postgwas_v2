import pandas as pd
import re
import subprocess



def run_bcftools_stats(
    vcf_path: str,
    external_af_name: str = "EUR",
    allelefreq_diff_cutoff: float = 0.2,
    threads: int = 5,
    bcftools_bin: str = "bcftools"
):
    """
    Run three bcftools operations (NO changes to user logic):

      1. INFO/<external_af_name> == "."
      2. abs(FORMAT/AF - INFO/<external_af_name>) <= cutoff
      3. bcftools stats

    bcftools_bin:
        Path to bcftools binary. Default: 'bcftools'.
    """
    results = {}

    # ----------------------------------------------------------------------
    # 1️⃣ Count sites where INFO/<external_af_name> == "."
    # ----------------------------------------------------------------------
    cmd_missing = (
        f'{bcftools_bin} view --threads {threads} --no-header '
        f'-i \'INFO/{external_af_name}=="."\' "{vcf_path}" | wc -l'
    )
    missing_count = int(subprocess.check_output(cmd_missing, shell=True).decode().strip())
    results["missing_extrnal_af"] = missing_count

    # ----------------------------------------------------------------------
    # 2️⃣ Count sites where abs(FORMAT/AF - INFO/<external_af_name>) <= cutoff
    # ----------------------------------------------------------------------
    cmd_diff = (
        f'{bcftools_bin} view --threads {threads} --no-header '
        f'-e \'abs(FORMAT/AF - INFO/{external_af_name}) <= {allelefreq_diff_cutoff}\' '
        f'"{vcf_path}" | wc -l'
    )
    diff_count = int(subprocess.check_output(cmd_diff, shell=True).decode().strip())
    results["af_diff_failed"] = diff_count

    # ----------------------------------------------------------------------
    # 3️⃣ Run bcftools stats
    # ----------------------------------------------------------------------
    stats_file = f"{vcf_path}.stats"
    cmd_stats = f'{bcftools_bin} stats "{vcf_path}" > "{stats_file}"'
    subprocess.run(cmd_stats, shell=True, check=True)
    results["stats_file"] = stats_file

    return results




def parse_bcftools_stats(stats_file: str):
    """
    Parse bcftools stats output into multiple pandas DataFrames.
    (NO logic changes requested by user)
    """
    sections = {
        "SN": [],
        "TSTV": [],
        "SiS": [],
        "AF": [],
        "QUAL": [],
        "IDD": [],
        "ST": [],
        "DP": []
    }

    with open(stats_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            tag = parts[0]
            if tag in sections:
                sections[tag].append(parts[1:])

    dfs = {}

    # SN
    if sections["SN"]:
        df_sn = pd.DataFrame(sections["SN"], columns=["id", "key", "value"])
        df_sn["value"] = pd.to_numeric(df_sn["value"], errors="ignore")
        dfs["SN"] = df_sn

    # All other sections (USE EXACT HEADERS FROM YOUR ORIGINAL CODE)
    header_map = {
        "TSTV": ["id", "ts", "tv", "ts/tv", "ts_1st", "tv_1st", "ts/tv_1st"],
        "SiS": ["id", "allele_count", "num_snps", "num_transitions",
                "num_transversions", "num_indels",
                "repeat_consistent", "repeat_inconsistent", "not_applicable"],
        "AF": ["id", "allele_freq", "num_snps", "num_ts", "num_tv",
               "num_indels", "repeat_consistent", "repeat_inconsistent", "not_applicable"],
        "QUAL": ["id", "quality", "num_snps", "ts_1st", "tv_1st", "num_indels"],
        "IDD": ["id", "length", "num_sites", "num_genotypes", "mean_vaf"],
        "ST": ["id", "type", "count"],
        "DP": ["id", "depth_bin", "num_genotypes", "frac_genotypes",
               "num_sites", "frac_sites"]
    }

    for sec, header in header_map.items():
        rows = sections.get(sec, [])
        if rows:
            df = pd.DataFrame(rows, columns=header)
            for c in header:
                df[c] = pd.to_numeric(df[c], errors="ignore")
            dfs[sec] = df

    return dfs



def bcftools_essential_summary(dfs: dict) -> pd.DataFrame:
    """
    Combine key bcftools-stats metrics into a single flat DataFrame.
    (NO logic changes — exact same behaviour)
    """
    summary = {}

    # SN
    if "SN" in dfs:
        sn = dfs["SN"].set_index("key")["value"]
        summary.update({
            "num_samples": sn.get("number of samples:", None),
            "num_records": sn.get("number of records:", None),
            "num_snps": sn.get("number of SNPs:", None),
            "num_indels": sn.get("number of indels:", None),
            "num_mnps": sn.get("number of MNPs:", None),
            "num_others": sn.get("number of others:", None),
            "num_multiallelic": sn.get("number of multiallelic sites:", None),
            "num_multiallelic_snp": sn.get("number of multiallelic SNP sites:", None),
        })

    # TSTV
    if "TSTV" in dfs:
        tstv = dfs["TSTV"].iloc[0]
        summary.update({
            "ts": tstv["ts"],
            "tv": tstv["tv"],
            "ts_tv_ratio": tstv["ts/tv"]
        })

    # SiS
    if "SiS" in dfs:
        sis = dfs["SiS"].iloc[0]
        summary.update({
            "singleton_count": sis["allele_count"],
            "singleton_snps": sis["num_snps"],
            "singleton_ts": sis["num_transitions"],
            "singleton_tv": sis["num_transversions"],
            "singleton_indels": sis["num_indels"]
        })

    # AF
    if "AF" in dfs:
        af = dfs["AF"].iloc[0]
        summary.update({
            "af_0_snps": af["num_snps"],
            "af_0_ts": af["num_ts"],
            "af_0_tv": af["num_tv"],
            "af_0_indels": af["num_indels"]
        })

    return pd.DataFrame([summary])




































# def run_bcftools_stats(vcf_path: str,
#                        external_af_name: str = "EUR",
#                        allelefreq_diff_cutoff: float = 0.2,
#                        threads: int = 5):
#     """
#     Run three bcftools-based operations:

#       1. Count sites where INFO/<external_af_name> == "."
#       2. Count sites where abs(FORMAT/AF - INFO/<external_af_name>) <= cutoff
#       3. Run bcftools stats and save <vcf>.stats

#     Returns a dict with all results.
#     """
#     results = {}
#     # --------------------------------------------------------
#     # 1️⃣ Count sites where INFO/EUR == "."
#     # --------------------------------------------------------
#     cmd_missing = (
#         f"""bcftools view --threads {threads} --no-header -i 'INFO/{external_af_name}=="."' {vcf_path} | wc -l """
#     )
    
#     missing_count = int(subprocess.check_output(cmd_missing, shell=True).decode().strip())
#     results["missing_extrnal_af"] = missing_count
#     # --------------------------------------------------------
#     # 2️⃣ Count sites where abs(FORMAT/AF - INFO/EUR) <= cutoff
#     # --------------------------------------------------------
#     cmd_diff = (f"""bcftools view --threads {threads} --no-header -e 'abs(FORMAT/AF - INFO/{external_af_name})\
#         <= {allelefreq_diff_cutoff}' {vcf_path} | wc -l """ )
#     diff_count = int(subprocess.check_output(cmd_diff, shell=True).decode().strip())
#     results["af_diff_failed"] = diff_count
#     # --------------------------------------------------------
#     # 3️⃣ Run bcftools stats → save .stats file
#     # --------------------------------------------------------
#     stats_file = f"{vcf_path}.stats"
#     cmd_stats = f"bcftools stats {vcf_path} > {stats_file}"
#     subprocess.run(cmd_stats, shell=True, check=True)
#     results["stats_file"] = stats_file
#     return results




# def parse_bcftools_stats(stats_file: str):
#     """
#     Parse bcftools stats output into multiple pandas DataFrames.
#     Returns:
#         dict of DataFrames:
#             {"SN": df_summary,
#                 "TSTV": df_tstv,
#                 "SiS": df_sis,
#                 "AF": df_af,
#                 "QUAL": df_qual,
#                 "IDD": df_idd,
#                 "ST": df_st,
#                 "DP": df_dp
#             }
#     """
#     sections = {
#         "SN": [],
#         "TSTV": [],
#         "SiS": [],
#         "AF": [],
#         "QUAL": [],
#         "IDD": [],
#         "ST": [],
#         "DP": []
#     }
#     with open(stats_file) as f:
#         for line in f:
#             if line.startswith("#"):
#                 continue  # ignore comments
#             parts = line.strip().split("\t")
#             tag = parts[0]
#             if tag in sections:
#                 sections[tag].append(parts[1:])  # store columns after the tag
#     # Convert each section to DataFrame
#     dfs = {}
#     # SN section must be converted to key-value table
#     if sections["SN"]:
#         df_sn = pd.DataFrame(sections["SN"], columns=["id", "key", "value"])
#         df_sn["value"] = pd.to_numeric(df_sn["value"], errors="ignore")
#         dfs["SN"] = df_sn
#     # Generic sections requiring simple table structure
#     for sec, header in {
#         "TSTV": ["id", "ts", "tv", "ts/tv", "ts_1st", "tv_1st", "ts/tv_1st"],
#         "SiS": ["id", "allele_count", "num_snps", "num_transitions", "num_transversions",
#                 "num_indels", "repeat_consistent", "repeat_inconsistent", "not_applicable"],
#         "AF": ["id", "allele_freq", "num_snps", "num_ts", "num_tv", "num_indels",
#                "repeat_consistent", "repeat_inconsistent", "not_applicable"],
#         "QUAL": ["id", "quality", "num_snps", "ts_1st", "tv_1st", "num_indels"],
#         "IDD": ["id", "length", "num_sites", "num_genotypes", "mean_vaf"],
#         "ST": ["id", "type", "count"],
#         "DP": ["id", "depth_bin", "num_genotypes", "frac_genotypes",
#                "num_sites", "frac_sites"]
#     }.items():
#         if sections[sec]:
#             df = pd.DataFrame(sections[sec], columns=header)
#             # Convert numeric columns where appropriate
#             for c in header:
#                 df[c] = pd.to_numeric(df[c], errors="ignore")
#             dfs[sec] = df
#     return dfs



# def bcftools_essential_summary(dfs: dict) -> pd.DataFrame:
#     """
#     Combine key bcftools-stats metrics into a single flat DataFrame.
#     Input: dict of dfs returned by parse_bcftools_stats()
#     Output: one-row essential summary DataFrame
#     """
#     summary = {}
#     # ----------------------------
#     # SN — "Summary Numbers"
#     # ----------------------------
#     if "SN" in dfs:
#         sn = dfs["SN"].set_index("key")["value"]
#         summary.update({
#             "num_samples": sn.get("number of samples:", None),
#             "num_records": sn.get("number of records:", None),
#             "num_snps": sn.get("number of SNPs:", None),
#             "num_indels": sn.get("number of indels:", None),
#             "num_mnps": sn.get("number of MNPs:", None),
#             "num_others": sn.get("number of others:", None),
#             "num_multiallelic": sn.get("number of multiallelic sites:", None),
#             "num_multiallelic_snp": sn.get("number of multiallelic SNP sites:", None),
#         })
#     # ----------------------------
#     # TSTV
#     # ----------------------------
#     if "TSTV" in dfs:
#         tstv = dfs["TSTV"].iloc[0]
#         summary.update({
#             "ts": tstv["ts"],
#             "tv": tstv["tv"],
#             "ts_tv_ratio": tstv["ts/tv"]
#         })
#     # ----------------------------
#     # SiS — Singleton stats
#     # ----------------------------
#     if "SiS" in dfs:
#         sis = dfs["SiS"].iloc[0]
#         summary.update({
#             "singleton_count": sis["allele_count"],
#             "singleton_snps": sis["num_snps"],
#             "singleton_ts": sis["num_transitions"],
#             "singleton_tv": sis["num_transversions"],
#             "singleton_indels": sis["num_indels"]
#         })
#     # ----------------------------
#     # AF — allele frequency table (first entry)
#     # ----------------------------
#     if "AF" in dfs:
#         af = dfs["AF"].iloc[0]
#         summary.update({
#             "af_0_snps": af["num_snps"],
#             "af_0_ts": af["num_ts"],
#             "af_0_tv": af["num_tv"],
#             "af_0_indels": af["num_indels"]
#         })
#     # ----------------------------
#     # QUAL — QC by quality
#     # ----------------------------
#     if "QUAL" in dfs:
#         q = dfs["QUAL"].iloc[0]
#         summary.update({
#             "qual_missing_snps": q["num_snps"],
#             "qual_indels": q["num_indels"]
#         })
#     # ----------------------------
#     # ST — substitution types (sum)
#     # ----------------------------
#     if "ST" in dfs:
#         summary["total_substitutions"] = dfs["ST"]["count"].sum()
#     return pd.DataFrame([summary])






# sumstat_vcf="/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/PGC3_SCZ_european/PGC3_SCZ_european_GRCh38_merged.vcf.gz"
# n_threads=5
# allelefreq_diff_cutoff=0.2
# external_af_name="EUR"
# output_folder="/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/PGC3_SCZ_european/"



# raw_data_qc1=run_bcftools_stats(vcf_path=sumstat_vcf,
#                        external_af_name=external_af_name,
#                        allelefreq_diff_cutoff=allelefreq_diff_cutoff,
#                        threads=n_threads)


# raw_data_dfs=parse_bcftools_stats(stats_file=f"{sumstat_vcf}.stats")
# raw_data_df=bcftools_essential_summary(raw_data_dfs)
# raw_data_df=raw_data_df.T
# raw_data_df.columns=["raw_variants"]
# raw_data_df.loc["variants_with_missing_extrnal_af"] = raw_data_qc1['missing_extrnal_af']
# raw_data_df.loc["variants_with_EAF_diff_gt_0.2"] = raw_data_qc1['af_diff_failed']


# qc_vcf=filter_gwas_vcf_bcftools(
#     vcf_path=sumstat_vcf,
#     output_folder=output_folder,
#     output_prefix=sumstat_vcf[:-5],
#     pval_cutoff = None,
#     maf_cutoff = None,
#     allelefreq_diff_cutoff = 0.2,
#     info_cutoff=0.7,
#     external_af_name = "EUR",
#     include_indels = True,
#     include_palindromic = True,
#     palindromic_af_lower = 0.4,
#     palindromic_af_upper = 0.6,
#     remove_mhc = False,
#     mhc_chrom = "6",
#     mhc_start = 25000000,
#     mhc_end = 34000000,
#     threads = 5,
#     max_mem="5G")


# filt_data_qc1=run_bcftools_stats(vcf_path=qc_vcf,
#                        external_af_name=external_af_name,
#                        allelefreq_diff_cutoff=allelefreq_diff_cutoff,
#                        threads=n_threads)

# filt_data_dfs=parse_bcftools_stats(stats_file=f"{qc_vcf}.stats")
# filt_data_df=bcftools_essential_summary(filt_data_dfs)
# filt_data_df=filt_data_df.T
# filt_data_df.columns=["qc_passed_variants"]
# filt_data_df.loc["variants_with_missing_extrnal_af"] = filt_data_qc1['missing_extrnal_af']
# filt_data_df.loc["variants_with_EAF_diff_gt_0.2"] = filt_data_qc1['af_diff_failed']
# final_qc_df=pd.concat([raw_data_df, filt_data_df], axis=1)