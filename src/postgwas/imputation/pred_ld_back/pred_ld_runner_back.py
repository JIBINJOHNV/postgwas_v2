
import os, re,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
import numpy as np
import polars as pl
from scipy.stats import norm






def run_pred_ld_parallel(
    sumstat_vcf: str,
    output_folder: str,
    output_prefix: str,
    pred_ld_ref: str,
    chromosomes: list = list(range(1, 23)) + ["X"],
    r2threshold: float = 0.8,
    maf: float = 0.001,
    population: str = "EUR",
    ref: str = "TOP_LD",
    threads: int = 6,
    pred_ld="/Users/JJOHN41/Documents/developing_software/postgwas/src/PRED-LD/pred_ld.py"
):
    """
    Parallelized PRED-LD preprocessing and imputation pipeline.
    Uses bcftools to extract variants per chromosome and runs pred_ld.py concurrently.
    """
    t0 = time.time()
    output_folder = os.path.abspath(output_folder.strip())
    sumstat_vcf = os.path.abspath(sumstat_vcf.strip())
    pred_ld_ref = os.path.abspath(pred_ld_ref.strip())
    os.makedirs(output_folder, exist_ok=True)
    
    def process_chrom(chrom):
        print(f"🧬 Starting chromosome {chrom}...")
        chr_output = Path(output_folder) / f"{output_prefix}_chr{chrom}_predld_input.tsv"
        bcftools_cmd = f"""
        (
          echo -e "snp\\tchr\\tpos\\tA1\\tA2\\tbeta\\tSE\\tNC\\tSS\\tAF\\tLP\\tSI" && \
          bcftools view \
              -r {chrom} \
              --min-alleles 2 --max-alleles 2 "{sumstat_vcf}" | \
          bcftools query -f '%CHROM:%POS:%REF:%ALT\\t%CHROM\\t%POS\\t%ALT\\t%REF\\t[%ES]\\t[%SE]\\t[%NC]\\t[%SS]\\t[%AF]\\t[%LP]\\t[%SI]\\n' | \
          sed 's|:|_|g'
        ) > "{chr_output}"
        """
        # --- Step 1: Run bcftools
        bcftools_log = Path(output_folder) / f"{output_prefix}_chr{chrom}_bcftools.log"
        with open(bcftools_log, "w") as log:
            subprocess.run(bcftools_cmd, shell=True, check=True, executable="/bin/bash", stdout=log, stderr=log)
        # --- Step 2: Run pred_ld.py
        pred_ld_cmd = [
            "python", pred_ld,
            "--file-path", str(chr_output),
            "--pop", population,
            "--ref", ref,
            "--ref_dir", pred_ld_ref,
            "--r2threshold", str(r2threshold),
            "--maf", str(maf),
            "--out_dir", output_folder
        ]
        predld_log = Path(output_folder) / f"{output_prefix}_chr{chrom}_predld.log"
        with open(predld_log, "w") as log:
            subprocess.run(pred_ld_cmd, check=True, stdout=log, stderr=log)
        
        expected_out = Path(output_folder) / f"imputation_results_chr{chrom}.txt"
        if not expected_out.exists():
            raise FileNotFoundError(f"⚠️ PRED-LD output not found for chr{chrom}")
        
        print(f"✅ Finished chromosome {chrom}")
        return chrom
    
    # --- Run in parallel (I/O bound so threads are fine)
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(process_chrom, c): c for c in chromosomes}
        for future in as_completed(futures):
            chrom = futures[future]
            try:
                future.result()
            except Exception as e:
                print(f"❌ Chromosome {chrom} failed: {e}")
    
    output_folder = Path(output_folder)
    prefix = output_prefix
    # Combine all bcftools logs
    subprocess.run(f"cat {output_folder}/{prefix}_chr*_bcftools.log > {output_folder}/{prefix}_bcftools.log", shell=True, check=False)
    subprocess.run(f"rm -f {output_folder}/{prefix}_chr*_bcftools.log", shell=True, check=False)
    # Combine all predld logs
    subprocess.run(f"cat {output_folder}/{prefix}_chr*_predld.log > {output_folder}/{prefix}_predld.log", shell=True, check=False)
    subprocess.run(f"rm -f {output_folder}/{prefix}_chr*_predld.log", shell=True, check=False)
    subprocess.run(f"rm -f {output_folder}/*_predld_input.tsv", shell=True, check=False)
    print(f"\n🎯 All chromosome jobs completed in {(time.time() - t0)/60:.2f} minutes.")




def process_pred_ld_results_all_parallel(
    folder_path: str,
    output_path: str | None = None,
    gwas2vcf_resource_folder: str | None = None,
    output_prefix: str | None = None,
    corr_method: str = "pearson",
    threads: int = 6,
):
    """Parallel Polars version of process_pred_ld_results_all with safe dtype coercion and CSV export."""
    # ---------- Helper: Safe type coercion ----------
    def coerce_dtypes(df: pl.DataFrame, dtype_map: dict) -> pl.DataFrame:
        """Coerce columns to specific dtypes (safe casting)."""
        for col, dtype in dtype_map.items():
            if col in df.columns:
                df = df.with_columns(pl.col(col).cast(dtype, strict=False))
        return df
    # ---------- Expected column types ----------
    data_types = {
        "chr": pl.Utf8, "snp": pl.Utf8, "A1": pl.Utf8, "A2": pl.Utf8,
        "pos": pl.Float64,
        "beta": pl.Float64, "SE": pl.Float64, "z": pl.Float64,
        "imputed": pl.Float64,
        "R2": pl.Float64, "NC": pl.Float64, "SS": pl.Float64,
        "AF": pl.Float64, "LP": pl.Float64, "SI": pl.Float64,
    }
    info_types = {
        "pos1": pl.Float64, "pos2": pl.Float64,
        "R2": pl.Float64, "Dprime": pl.Float64,
        "ALT_AF1": pl.Float64, "ALT_AF2": pl.Float64,
        "+/-corr": pl.Utf8, "rsID1": pl.Utf8, "rsID2": pl.Utf8,
        "REF1": pl.Utf8, "ALT1": pl.Utf8, "REF2": pl.Utf8, "ALT2": pl.Utf8,
    }
    # ---------- Locate files ----------
    all_files = os.listdir(folder_path)
    imp_files = [f for f in all_files if f.startswith("imputation_results_chr") and f.endswith(".txt")]
    info_files = [f for f in all_files if f.startswith("LD_info_TOP_LD_chr") and f.endswith(".txt")]
    if not imp_files or not info_files:
        raise FileNotFoundError("❌ Could not find imputation or LD info files in the folder.")
    
    def _chr_sort_key(x):
        return (0, int(x)) if x.isdigit() else (1, {"X": 23, "Y": 24, "MT": 25, "M": 25}.get(x.upper(), 99))
    
    chromosomes = sorted(set(re.findall(r"chr(\w+)", " ".join(imp_files))), key=_chr_sort_key)
    print(f"🧬 Found chromosomes: {', '.join(chromosomes)}")
    results, summaries = [], []
    # ---------- Per-chromosome worker ----------
    def process_chr(chr_id):
        try:
            imp_file = os.path.join(folder_path, f"imputation_results_chr{chr_id}.txt")
            info_file = os.path.join(folder_path, f"LD_info_TOP_LD_chr{chr_id}.txt")
            if not (os.path.exists(imp_file) and os.path.exists(info_file)):
                print(f"⚠️ Missing files for chr{chr_id}, skipping…")
                return None, None
            
            info_df = coerce_dtypes(pl.read_csv(info_file, separator="\t", infer_schema_length=5000), info_types)
            data_df = coerce_dtypes(pl.read_csv(imp_file, separator="\t", infer_schema_length=5000), data_types)
            # Duplicated SNPs
            duplicated_df = data_df.filter(pl.col("snp").is_duplicated())
            imputed_dups = (
                duplicated_df.filter(pl.col("imputed") == 1)
                .drop(["NC", "SS", "AF", "LP", "SI"], strict=False)
                .sort(["chr", "pos", "snp"])
            )
            not_imputed_dups = (
                duplicated_df.filter(pl.col("imputed") != 1)
                .drop(["NC", "SS", "AF", "LP", "SI"], strict=False)
                .sort(["chr", "pos", "snp"])
            )
            beta_corr, z_corr = np.nan, np.nan
            if imputed_dups.height > 0 and not_imputed_dups.height > 0:
                beta_corr = np.corrcoef(imputed_dups["beta"].to_numpy(), not_imputed_dups["beta"].to_numpy())[0, 1]
                z_corr = np.corrcoef(imputed_dups["z"].to_numpy(), not_imputed_dups["z"].to_numpy())[0, 1]
            # Non-imputed data
            not_imputed_df = (
                data_df.filter(pl.col("imputed") == 0)
                .drop(["imputed", "R2"], strict=False)
                .with_columns((10 ** (-pl.col("LP").cast(pl.Float64, strict=False))).alias("p_value"))
                .drop("LP", strict=False)
            )
            # Imputed data
            imputed_df = (
                data_df.filter(pl.col("imputed") == 1)
                .filter(~pl.col("snp").is_in(not_imputed_df["snp"]))
                .drop(["imputed", "R2", "NC", "SS", "AF", "LP", "SI"], strict=False)
            )
            # Merge LD info with not-imputed stats
            merged = info_df.join(
                not_imputed_df.select(["snp", "NC", "SS", "AF", "SI"]),
                left_on="rsID1", right_on="snp", how="inner"
            )
            imputed_stats = (
                merged.group_by("rsID2")
                .agg([
                    pl.col("NC").mean(),
                    pl.col("SS").mean(),
                    pl.col("AF").mean(),
                    pl.col("SI").mean(),
                ])
                .rename({"rsID2": "snp"})
            )
            imputed_df = imputed_df.join(imputed_stats, on="snp", how="left")
            imputed_df = imputed_df.with_columns(
                (2 * pl.Series(norm.sf(np.abs(imputed_df["z"].to_numpy())))).alias("p_value")
            )
            # Combine
            not_imputed_df = not_imputed_df.with_columns(pl.lit("no").alias("imputed"))
            imputed_df = imputed_df.with_columns(pl.lit("yes").alias("imputed"))
            final_df = pl.concat([not_imputed_df, imputed_df], how="diagonal")
            # Round NC/SS and derive ncase_col
            for col in ["NC", "SS"]:
                if col in final_df.columns:
                    final_df = final_df.with_columns(pl.col(col).round(0).cast(pl.Int64, strict=False))
            
            final_df = final_df.with_columns(
                pl.when(pl.col("SS") > pl.col("NC"))
                .then((pl.col("SS") - pl.col("NC")).cast(pl.Int64, strict=False))
                .otherwise(0)
                .alias("ncase_col")
            )
            summary = {
                "chromosome": chr_id,
                "n_imputed_dups": imputed_dups.height,
                "n_nonimputed_dups": not_imputed_dups.height,
                "imputed_markers": imputed_df.height,
                "beta_corr": beta_corr,
                "z_corr": z_corr,
            }
            print(f"✅ chr{chr_id}: β_corr={beta_corr:.3f} | z_corr={z_corr:.3f}")
            return final_df, summary
        
        except Exception as e:
            print(f"❌ Error in chr{chr_id}: {e}")
            return None, None
    # ---------- Parallel execution ----------
    with ThreadPoolExecutor(max_workers=threads) as executor:
        future_to_chr = {executor.submit(process_chr, c): c for c in chromosomes}
        for future in as_completed(future_to_chr):
            df, summ = future.result()
            if df is not None:
                results.append(df)
            if summ is not None:
                summaries.append(summ)
    if not results:
        raise ValueError("❌ No chromosome processed successfully!")
    # ---------- Combine results ----------
    combined_df = pl.concat(results, how="diagonal")
    corr_df = pl.DataFrame(summaries)
    # ---------- Output paths ----------
    output_folder = output_path or folder_path
    os.makedirs(output_folder, exist_ok=True)
    output_prefix = output_prefix or "PRED_LD"
    combined_path = f"{output_folder}/{output_prefix}_PREDLD_allchr.tsv"
    corr_path = f"{output_folder}/{output_prefix}_PREDLD_correlations.tsv"
    config_path = f"{output_folder}/{output_prefix}_gwas2vcf_config.csv"
    # ---------- Save outputs ----------
    combined_df.write_csv(combined_path, separator="\t")
    corr_df.write_csv(corr_path, separator="\t")
    os.system(f"gzip {combined_path}")
    print(f"\n💾 Saved combined results → {combined_path}")
    print(f"💾 Saved correlation summary → {corr_path}")
    # ---------- Create gwas2vcf config ----------
    columns = [
        "sumstat_file", "gwas_outputname", "chr_col", "pos_col", "snp_id_col",
        "ea_col", "oa_col", "eaf_col", "beta_or_col", "se_col", "imp_z_col",
        "pval_col", "ncontrol_col", "ncase_col", "ncontrol", "ncase",
        "imp_info_col", "delimiter", "infofile", "infocolumn", "eaffile",
        "eafcolumn", "liftover", "chr_pos_col", "resourse_folder", "output_folder"
    ]
    gwas2vcf_spec_df = pd.DataFrame([[None] * len(columns)], columns=columns)
    # Basic fields
    gwas2vcf_spec_df.loc[0, [
        "sumstat_file", "gwas_outputname", "output_folder", "resourse_folder", "delimiter"
    ]] = [f"{combined_path}.gzip", f"{output_prefix}_imputed", output_folder, gwas2vcf_resource_folder, "\t"]
    # Expected mappings → check existence
    mapping = {
        "chr_col": "chr", "pos_col": "pos", "snp_id_col": "snp",
        "ea_col": "A1", "oa_col": "A2", "eaf_col": "AF",
        "beta_or_col": "beta", "se_col": "SE", "pval_col": "p_value",
        "ncontrol_col": "NC", "ncase_col": "ncase_col", "imp_info_col": "SI"
    }
    for spec_col, df_col in mapping.items():
        gwas2vcf_spec_df.loc[0, spec_col] = df_col if df_col in combined_df.columns else "NA"
    gwas2vcf_spec_df.to_csv(config_path, index=False)
    print(f"💾 Saved GWAS2VCF config → {config_path}")
    print("🎯 All chromosomes processed successfully (parallel).")
    # Combine all per-chromosome imputation result files
    subprocess.run(f"cat {folder_path}/imputation_results_chr*.txt | gzip -c > {folder_path}/{output_prefix}_imputation_results.txt.gz", shell=True, check=False)
    subprocess.run(f"cat {folder_path}/LD_info_TOP_LD_chr*.txt | gzip -c > {folder_path}/{output_prefix}_LD_info_TOP_LD.txt.gz", shell=True, check=False)
    # Remove the individual chromosome files
    subprocess.run(f"rm -f {folder_path}/imputation_results_chr*.txt", shell=True, check=False)
    subprocess.run(f"rm -f {folder_path}/LD_info_TOP_LD_chr*.txt", shell=True, check=False)
    return combined_df, corr_df