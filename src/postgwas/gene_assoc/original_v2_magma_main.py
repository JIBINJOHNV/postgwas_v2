import subprocess,os
import logging
from pathlib import Path
import numpy as np
import pandas as pd
import polars as pl
from statsmodels.stats.multitest import multipletests
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import Pool, cpu_count
from functools import partial
from typing import Optional

from postgwas.utils.main import decide_magma_batches_from_annot

# ==========================================
# 0. LOGGING HELPER
# ==========================================

def setup_logger(log_file: str, name: str = "magma_logger"):
    """
    Configures a logger to write to a file ONLY (no console output).
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    
    # Clear existing handlers to avoid duplicate logs if function is re-run
    if logger.hasHandlers():
        logger.handlers.clear()
        
    # File Handler
    file_handler = logging.FileHandler(log_file)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    
    logger.addHandler(file_handler)
    logger.propagate = False # Prevent bubbling up to root logger (which might print to screen)
    
    return logger

def run_subprocess_with_logging(cmd: list, logger: logging.Logger):
    """
    Helper to run subprocess commands and redirect output/errors to the log file.
    """
    try:
        # capture_output=True prevents printing to screen
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        # Log standard output if verbose logging is desired, or just success
        # logger.info(f"CMD Output: {result.stdout}") 
    except subprocess.CalledProcessError as e:
        logger.error(f"❌ Command failed: {' '.join(cmd)}")
        logger.error(f"Error details: {e.stderr}")
        raise e



# ==========================================
# 1. MAGMA EXECUTION
# ==========================================
def run_magma_analysis(
    magma_analysis_folder: str,
    sample_id: str,
    ld_ref: str,
    gene_loc_file: str,
    snp_loc_file: str,
    pval_file: str,
    log_file: str,

    window_upstream: int = 35,
    window_downstream: int = 10,
    gene_model: str = "snp-wise=mean",
    n_sample_col: str = "N_COL",
    num_cores: int = 6,
    seed: int = 10,
    magma: str = "magma",
    geneset_file: str = None
):
    """
    Run full MAGMA analysis pipeline:
      1. SNP → gene annotation
      2. Gene analysis (batched)
      3. Merge results
      4. OPTIONAL: gene-set enrichment

    ALWAYS RETURNS:
      (gene_annot_out, merged_prefix)
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed
    import subprocess
    from pathlib import Path
    
    logger = setup_logger(log_file)
    magma_analysis_folder = Path(magma_analysis_folder)
    magma_analysis_folder.mkdir(parents=True, exist_ok=True)
    # Output prefixes
    gene_annot_out = magma_analysis_folder / f"{sample_id}_magma"
    merged_prefix  = magma_analysis_folder / f"{sample_id}_magma_{window_upstream}up_{window_downstream}down"
    merged_raw     = merged_prefix.with_suffix(".genes.raw")
    # ---------------------------------------------------------
    # Step 1 — Annotation
    # ---------------------------------------------------------
    logger.info("🔹 [Step 1/4] Running MAGMA annotation...")
    print("                 🔹 [Step 1/4] Running MAGMA annotation...")
    annotate_cmd = [
        magma,
        "--annotate",
        f"window={window_upstream},{window_downstream}",
        "--snp-loc", str(snp_loc_file),
        "--gene-loc", str(gene_loc_file),
        "--out", str(gene_annot_out),
    ]
    run_subprocess_with_logging(annotate_cmd, logger)
    
    annot_path = f"{gene_annot_out}.genes.annot"

    with open(annot_path, "r") as f:
        n_lines = sum(1 for _ in f)

    if n_lines > 1000:
        num_batches = decide_magma_batches_from_annot(annot_file=annot_path)
        #num_batches = min(num_batches, num_cores)
    else:
        num_batches = 1

    # ---------------------------------------------------------
    # Step 2 — Gene analysis
    # ---------------------------------------------------------
    logger.info(
        f"🔹 [Step 2/4] Running MAGMA gene analysis "
        f"({'single run' if num_batches == 1 else f'{num_batches} batches'})..."
    )
    print(
        f"                 🔹 [Step 2/4] Running MAGMA gene analysis "
        f"({'single run' if num_batches == 1 else f'{num_batches} batches'})..."
    )

    def run_batch(i=None):
        batch_cmd = [
            magma,
            "--bfile", str(ld_ref),
            "--gene-annot", f"{gene_annot_out}.genes.annot",
            "--pval", str(pval_file), f"ncol={n_sample_col}",
            "--gene-model", gene_model,
            "--out", str(gene_annot_out),
        ]

        # ✅ ONLY add --batch when MAGMA allows it
        if num_batches > 1:
            batch_cmd.extend(["--batch", str(i), str(num_batches),])

        subprocess.run(batch_cmd, check=True, capture_output=True, text=True)

    # ---------------------------------------------------------
    # Execute
    # ---------------------------------------------------------
    if num_batches == 1:
        run_batch()
        os.system(f"mv {gene_annot_out}.genes.raw {merged_prefix}.genes.raw")
        os.system(f"mv {gene_annot_out}.genes.out {merged_prefix}.genes.out")
    else:
        with ThreadPoolExecutor(max_workers=num_cores) as executor:
            futures = [
                executor.submit(run_batch, i)
                for i in range(1, num_batches + 1)
            ]
            for f in as_completed(futures):
                f.result()

    # ---------------------------------------------------------
    # Step 3 — Merge MAGMA batches (ONLY if needed)
    # ---------------------------------------------------------
    if num_batches > 1:
        logger.info("🔹 [Step 3/4] Merging MAGMA batch results...")
        print("                 🔹 [Step 3/4] Merging MAGMA batch results...")

        merge_cmd = [
            magma,
            "--merge", str(gene_annot_out),
            "--out", str(merged_prefix),
        ]
        run_subprocess_with_logging(merge_cmd, logger)

        logger.info("🗑️ Cleaning up intermediate batch files...")
        for f in magma_analysis_folder.glob(f"{sample_id}_magma.batch*"):
            try:
                f.unlink()
            except Exception as e:
                logger.warning(f"Could not delete {f}: {e}")

    else:
        logger.info("🔹 [Step 3/4] Single batch run — merge step skipped.")
        print("                 🔹 [Step 3/4] Single batch run — merge step skipped.")


    # ---------------------------------------------------------
    # Step 4 — Gene-set enrichment (optional)
    # ---------------------------------------------------------
    if geneset_file and Path(geneset_file).exists():
        logger.info("🔹 [Step 4/4] Running MAGMA gene-set analysis...")
        print("                 🔹 [Step 4/4] Running MAGMA gene-set analysis...")
        geneset_cmd = [
            magma,
            "--gene-results", str(merged_raw),
            "--set-annot", str(geneset_file),
            "--out", str(magma_analysis_folder / sample_id),
            "--seed", str(seed),
        ]
        run_subprocess_with_logging(geneset_cmd, logger)
        logger.info("✅ MAGMA full analysis completed successfully.")
    else:
        logger.info("ℹ️ No geneset file provided — skipping gene-set analysis.")
        logger.info("✅ MAGMA gene-level analysis completed successfully.")
    # ALWAYS return key files for PoPS
    return {
        "gene_annot": str(gene_annot_out),
        "merged_prefix": str(merged_prefix),
    }


# ==========================================
# 2. STATISTICAL CORRECTION
# ==========================================

def correct_p_values(
    magma_analysis_folder: str, 
    sample_id: str, 
    log_file: str,
    output_file: Optional[str] = None,
) -> pl.DataFrame:
    logger = setup_logger(log_file)
    folder = Path(magma_analysis_folder)
    if output_file is None:
        output_file = folder / f"{sample_id}.gsa_corrected.tsv"
    
    gsa_out = folder / f"{sample_id}.gsa.out"
    logger.info(f"📉 Reading MAGMA output from: {gsa_out}")
    # Read Data
    try:
        df = pd.read_fwf( gsa_out, comment="#", infer_nrows=500,)
        df = pl.from_pandas(df)
    except Exception as e:
        logger.error(f"❌ Failed to read MAGMA output file: {e}")
        raise e
    # Cleanup columns
    df = df.select([col for col in df.columns if col != ""])
    required_cols = ["VARIABLE", "P"]
    if not all(col in df.columns for col in required_cols):
        msg = f"⚠️ Missing columns. Found: {df.columns}"
        logger.error(msg)
        raise ValueError(msg)
    
    if "FULL_NAME" not in df.columns:
        df = df.rename({"VARIABLE": "FULL_NAME"})
    n_tests = df.height
    p_values = df["P"].to_numpy()
    logger.info(f"📉 applying multiple testing correction for {n_tests} tests...")
    # Global Corrections
    def get_fdr(p): return multipletests(p, method="fdr_bh")[1]
    def get_holm(p): return multipletests(p, method="holm")[1]
    
    df = df.with_columns([
        pl.col("P").mul(n_tests).clip(upper_bound=1.0).alias("P_bonferroni_corr"),
        (1 - (1 - pl.col("P")).pow(n_tests)).alias("P_sidak_corr"),
        pl.Series(get_holm(p_values)).alias("P_holm_corr"),
        pl.Series(get_fdr(p_values)).alias("P_fdr_bh_corr")
    ])
    # Category Specific Corrections
    categories = {
        "go_kegg_reactome": r"^(GOBP_|GOCC_|GOMF_|KEGG_|REACTOME_)",
        "go": r"^(GOBP_|GOCC_|GOMF_)",
        "reactome": r"^REACTOME_",
        "kegg": r"^KEGG_"
    }
    for prefix, pattern in categories.items():
        mask = pl.col("FULL_NAME").str.contains(pattern)
        subset = df.filter(mask)
        if subset.height > 0:
            sub_p = subset["P"].to_numpy()
            sub_bonf = np.minimum(sub_p * len(sub_p), 1.0)
            sub_fdr = multipletests(sub_p, method="fdr_bh")[1]
            
            temp_df = subset.select("FULL_NAME").with_columns([
                pl.Series(sub_bonf).alias(f"{prefix}_P_bonferroni_corr"),
                pl.Series(sub_fdr).alias(f"{prefix}_P_fdr_bh_corr")
            ])
            df = df.join(temp_df, on="FULL_NAME", how="left")
        else:
            df = df.with_columns([
                pl.lit(None).alias(f"{prefix}_P_bonferroni_corr"),
                pl.lit(None).alias(f"{prefix}_P_fdr_bh_corr")
            ])
    df = df.sort("P")
    df.write_csv(output_file, separator="\t")
    logger.info(f"✅ Correction completed. Results saved to: {output_file}")
    return df


# ==========================================
# 3. GENE SET PROCESSING
# ==========================================

def parse_gmt_file(gmt_path: Path) -> pl.DataFrame:
    data = []
    with open(gmt_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                data.append({"V1": parts[0], "V2": parts[1], "V3": ",".join(parts[2:])})
    return pl.DataFrame(data)

def add_gene_details(gene_string: str, gsa_genes_dict: dict):
    if not gene_string:
        return {
            "CommonGenes": None,
            "CommonGenes_Pvalues": None,
            "TotalGenes": 0,
            "CommonGeneCount": 0,
        }
    gene_list = [g.strip() for g in gene_string.split(",") if g.strip()]
    total_genes = len(gene_list)
    common_genes_data = [
        (g, gsa_genes_dict[g]) for g in gene_list if g in gsa_genes_dict
    ]
    if not common_genes_data:
        return {
            "CommonGenes": None,
            "CommonGenes_Pvalues": None,
            "TotalGenes": total_genes,
            "CommonGeneCount": 0,
        }
    # sort alphabetically for consistent output
    common_genes_data.sort(key=lambda x: x[0])
    genes = [x[0] for x in common_genes_data]
    pvals = [str(x[1]) for x in common_genes_data]
    return {
        "CommonGenes": ",".join(genes),
        "CommonGenes_Pvalues": ",".join(pvals),
        "TotalGenes": total_genes,
        "CommonGeneCount": len(genes),
    }


def process_gene_sets(
    magma_analysis_folder: str, 
    sample_id: str, 
    geneset_file: str, 
    corrected_df: pl.DataFrame,
    log_file: str # <--- NEW ARGUMENT
):
    logger = setup_logger(log_file)
    folder = Path(magma_analysis_folder)
    genes_out_file = folder / f"{sample_id}_magma_35up_10down.genes.out"
    logger.info("📊 Processing gene sets and overlaps...")
    gsa_genes_df = pd.read_fwf( genes_out_file, comment="#", infer_nrows=500,)
    gsa_genes_df = pl.from_pandas(gsa_genes_df).select([col for col in ["GENE", "P"] if col])
    gsa_genes_dict = dict(zip(gsa_genes_df["GENE"], gsa_genes_df["P"]))
    geneset_file = Path(geneset_file)
    logger.info("📂 Parsing GMT file...")
    msigdb_df = parse_gmt_file(geneset_file)
    n_cores = max(1, cpu_count() - 1)
    logger.info(f"🧠 Using {n_cores} cores for parallel overlap computation...")
    func = partial(add_gene_details, gsa_genes_dict=gsa_genes_dict)
    gene_strings = msigdb_df["V3"].to_list()
    with Pool(processes=n_cores) as pool:
        results = pool.map(func, gene_strings)
    result_df = pl.DataFrame(results)
    msigdb_combined = pl.concat([msigdb_df, result_df], how="horizontal")
    logger.info("🔄 Merging statistical results...")
    merged_df = corrected_df.join(
        msigdb_combined,
        left_on="FULL_NAME",
        right_on="V1",
        how="left"
    ).sort("P")
    merged_df = merged_df.with_columns(pl.lit(sample_id).alias("sample_id"))
    output_file = folder / f"{sample_id}.gsa_corrected_annotated.tsv"
    merged_df.write_csv(output_file, separator="\t")
    logger.info(f"✅ Final pipeline output saved to: {output_file}")
    return merged_df


# ==========================================
# 4. MAIN PIPELINE WRAPPER
# ==========================================
def magma_analysis_pipeline(
    output_dir: str,
    sample_id: str,
    ld_ref: str,
    gene_loc_file: str,
    snp_loc_file: str,
    pval_file: str,
    geneset_file: str,
    log_file: str,

    # MAGMA parameters
    threads: int = 6,
    window_upstream: int = 35,
    window_downstream: int = 10,
    gene_model: str = "snp-wise=mean",
    n_sample_col: str = "N_COL",
    seed: int = 10,
    magma: str = "magma"
):
    """
    Full MAGMA analysis pipeline:
        1. Annotate → Gene analysis → Merge → Gene-set analysis
        2. Correct P-values
        3. Add gene overlaps, metadata
    """

    # ---------------------------------------------------------
    # Step 1 — Run MAGMA core analysis
    # ---------------------------------------------------------
    magma_gene_file=run_magma_analysis(
        magma_analysis_folder=output_dir,
        sample_id=sample_id,
        ld_ref=ld_ref,
        gene_loc_file=gene_loc_file,
        snp_loc_file=snp_loc_file,
        pval_file=pval_file,
        geneset_file=geneset_file,
        log_file=log_file,
        n_sample_col=n_sample_col,
        num_cores=threads,
        window_upstream=window_upstream,
        window_downstream=window_downstream,
        gene_model=gene_model,
        seed=seed,
        magma=magma
    )

    if geneset_file and Path(geneset_file).exists():
        # ---------------------------------------------------------
        # Step 2 — Correct gene-set p-values
        # ---------------------------------------------------------
        corrected_df = correct_p_values(
            magma_analysis_folder=output_dir,
            sample_id=sample_id,
            log_file=log_file,
        )

        # ---------------------------------------------------------
        # Step 3 — Process gene-sets + add gene overlaps
        # ---------------------------------------------------------
        gene_set_summary = process_gene_sets(
            magma_analysis_folder=output_dir,
            sample_id=sample_id,
            geneset_file=geneset_file,
            corrected_df=corrected_df,
            log_file=log_file,
        )
        return {
            "magma_pathway": str(gene_set_summary),
            "magma_genes_raw": f"{magma_gene_file['merged_prefix']}.genes.raw",
            "magma_genes_out": f"{magma_gene_file['merged_prefix']}.genes.out",
        }
    else:
        print("         ❌ Pathway analysis not performed because the gene set file was not provided or does not exist.")
        return {
            "magma_genes_raw": f"{magma_gene_file['merged_prefix']}.genes.raw",
            "magma_genes_out": f"{magma_gene_file['merged_prefix']}.genes.out",
            "magma_genes_prefix": f"{magma_gene_file['merged_prefix']}",
        }



# return {
#     "gene_annot": str(gene_annot_out),
#     "merged_prefix": str(merged_prefix),
# }



    # magma_Analysis_folder <- glue("{output_folder}/magma_Analysis/{sample_id}/")
    # dir.create(magma_Analysis_folder, showWarnings = FALSE, recursive = TRUE)

    # df=process_vcf_parallel(output_folder, sample_id, num_threads = 5)
    # df=process_variant_data(df, dbsnp_df = NULL, variant_id_type = variant_id_type, mhc_region = mhc_region)

    # prepare_magma_files(df, magma_Analysis_folder, sample_id) 
    # run_magma_analysis(magma_Analysis_folder, sample_id, ld_ref,gene_loc_file, geneset_file, num_batches = 6, num_cores = 6)
    # corrected_df <- correct_p_values(magma_Analysis_folder, sample_id)

    # process_gene_sets(magma_Analysis_folder, sample_id, geneset_file, corrected_df) 





# def magma_input(sumstat_vcf: str, output_folder: str, sample_name: str):
#     """
#     Convert a single summary-statistics VCF file into a MAGMA-compatible
#     tab-delimited text file.
#     Parameters
#     ----------
#     sumstat_vcf : str
#         Path to the input VCF (e.g., "PGC3_SCZ_european_chr1.vcf.gz")
#     output_folder : str
#         Folder where the converted text file will be written
#     sample_name : str
#         Prefix for output naming (e.g., "PGC3_SCZ_european")
#     Output
#     ------
#     <output_folder>/<sample_name>_<chromosome>_magma_input.txt
#     """
#     vcf_path = Path(sumstat_vcf)
#     output_dir = Path(output_folder)
#     output_dir.mkdir(parents=True, exist_ok=True)

#     # Output file naming
#     output_file = output_dir / f"{sample_name}_magma_input.tsv"

#     # Build bcftools → query → sed pipeline
#     cmd = f"""(
#         printf "SNP\\tCHR\\tBP\\tREF\\tALT\\tBETA\\tSE\\tN_COL\\tAF\\tLP\\n"
#         bcftools view --min-alleles 2 --max-alleles 2 "{vcf_path}" | \
#         bcftools query -f '%CHROM:%POS:%REF:%ALT\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t[%ES]\\t[%SE]\\t[%NEF]\\t[%AF]\\t[%LP]\\n' | \
#         sed 's|:|_|g'
#     ) > "{output_file}"
#     """
#     try:
#         subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
#         print(f"✅ Successfully created: {output_file}")
#     except subprocess.CalledProcessError as e:
#         print(f"❌ Conversion failed for {vcf_path}: {e}")
    
#     return str(output_file)


# def prepare_magma_files(input_file: str, magma_analysis_folder: str, sample_id: str):
#     """
#     Prepare MAGMA input files from a summary-statistic text file (output of magma_input()).
#     Parameters
#     ----------
#     input_file : str
#         Path to the input tab-delimited file (e.g., output of magma_input()).
#     magma_analysis_folder : str
#         Directory where MAGMA input files will be written.
#     sample_id : str
#         Identifier prefix for output filenames (e.g., "PGC3_SCZ_european_1").
#     Outputs
#     --------
#     - <sample_id>_magma_input.tsv
#     - <sample_id>_magma_snp_loc.tsv
#     - <sample_id>_magma_P_val.tsv
#     """
#     input_path = Path(input_file)
#     output_dir = Path(magma_analysis_folder)
#     output_dir.mkdir(parents=True, exist_ok=True)
#     # ---------- Read file ----------
#     df = pl.read_csv(input_path, separator="\t")
#     # ---------- Check for required columns ----------
#     required_cols = ["SNP", "CHR", "BP", "REF", "ALT", "LP", "BETA", "SE", "AF", "N_COL"]
#     missing = [c for c in required_cols if c not in df.columns]
#     if missing:
#         raise ValueError(f"⚠️ Missing columns in {input_file}: {', '.join(missing)}")
#     # ---------- Convert -log10(P) → P ----------
#     df = df.with_columns((10 ** (-pl.col("LP"))).alias("P"))
#     # ---------- Select and reorder ----------
#     df2 = df.select([
#         pl.col("SNP"),
#         pl.col("CHR"),
#         pl.col("BP"),
#         pl.col("REF"),
#         pl.col("ALT"),
#         pl.col("P"),
#         pl.col("BETA"),
#         pl.col("SE"),
#         pl.col("AF"),
#         pl.col("N_COL")
#     ])
#     # ---------- Define output files ----------
#     main_file = output_dir / f"{sample_id}_magma_input.tsv"
#     snp_loc_file = output_dir / f"{sample_id}_magma_snp_loc.tsv"
#     pval_file = output_dir / f"{sample_id}_magma_P_val.tsv"
#     # ---------- Write files ----------
#     df2.write_csv(main_file, separator="\t")
#     df2.select(["SNP", "CHR", "BP"]).write_csv(snp_loc_file, separator="\t")
#     df2.select(["SNP", "CHR", "BP", "P", "N_COL"]).write_csv(pval_file, separator="\t")
#     print(f"✅ MAGMA input files successfully generated in: {output_dir}")
#     print(f"📄 Main file: {main_file.name}")
#     print(f"📄 SNP location file: {snp_loc_file.name}")
#     print(f"📄 P-value file: {pval_file.name}")
#     return {
#         "main_file": str(main_file),
#         "snp_loc_file": str(snp_loc_file),
#         "pval_file": str(pval_file),
#     }
    