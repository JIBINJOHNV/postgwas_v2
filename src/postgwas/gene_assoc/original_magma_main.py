import subprocess
from pathlib import Path
import pandas as pd
import numpy as np
import polars as pl
from statsmodels.stats.multitest import multipletests
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import Pool, cpu_count
from functools import partial

    
    
def run_magma_analysis(
    magma_analysis_folder: str,
    sample_id: str,
    ld_ref: str,
    gene_loc_file: str,
    geneset_file: str,
    num_batches: int = 6,
    num_cores: int = 6,
):
    """
    Run complete MAGMA analysis (annotation → gene analysis → merge → gene-set test).

    Parameters
    ----------
    magma_analysis_folder : str
        Path where MAGMA input and output files are stored
    sample_id : str
        Prefix used in filenames (e.g. 'PGC3_SCZ_european_1')
    ld_ref : str
        Path to LD reference binary prefix (--bfile), e.g. 1000G_EUR/1000G_EUR
    gene_loc_file : str
        Path to gene location file (GRCh37 or GRCh38 MAGMA format)
    geneset_file : str
        Path to gene set annotation file (.gmt or MAGMA .geneset)
    num_batches : int
        Number of genome-wide batches to run in parallel
    num_cores : int
        Number of parallel CPU threads
    """
    magma_analysis_folder = Path(magma_analysis_folder)
    # Input files generated from previous steps
    snp_loc_file   = magma_analysis_folder / f"{sample_id}_magma_snp_loc.tsv"
    pval_file      = magma_analysis_folder / f"{sample_id}_magma_P_val.tsv"
    gene_annot_out = magma_analysis_folder / f"{sample_id}_magma"
    # Expected outputs
    merged_prefix  = magma_analysis_folder / f"{sample_id}_magma_35up_10down"
    merged_raw     = merged_prefix.with_suffix(".genes.raw")
    print("🔹 Running MAGMA annotation...")
    annotate_cmd = [
        "magma",
        f"--annotate", "window=35,10",
        f"--snp-loc", str(snp_loc_file),
        f"--gene-loc", str(gene_loc_file),
        f"--out", str(gene_annot_out)
    ]
    subprocess.run(annotate_cmd, check=True)
    print("🔹 Running MAGMA SNP-wise mean gene analysis (parallel)...")
    def run_batch(i: int):
        batch_cmd = [
            "magma",
            f"--bfile", str(ld_ref),
            f"--gene-annot", f"{gene_annot_out}.genes.annot",
            f"--pval", str(pval_file), "ncol=N_COL",
            "--gene-model", "snp-wise=mean",
            "--batch", str(i), str(num_batches),
            "--out", str(gene_annot_out)
        ]
        subprocess.run(batch_cmd, check=True)
    # Run batches in parallel
    with ThreadPoolExecutor(max_workers=num_cores) as executor:
        futures = [executor.submit(run_batch, i) for i in range(1, num_batches + 1)]
        for _ in as_completed(futures):
            pass
    print("🔹 Merging MAGMA batch results...")
    merge_cmd = [
        "magma",
        "--merge", str(gene_annot_out),
        "--out", str(merged_prefix)
    ]
    subprocess.run(merge_cmd, check=True)
    print("🗑️ Cleaning up intermediate batch files...")
    for f in magma_analysis_folder.glob(f"{sample_id}_magma.batch*"):
        f.unlink()
    print("🔹 Running gene set analysis...")
    geneset_cmd = [
        "magma",
        "--gene-results", str(merged_raw),
        "--set-annot", str(geneset_file),
        "--out", str(magma_analysis_folder / sample_id),
        "--seed", "10"
    ]
    subprocess.run(geneset_cmd, check=True)
    print("✅ MAGMA analysis completed successfully!\n")




def correct_p_values(magma_analysis_folder: str, sample_id: str, output_file: str | None = None):
    """
    Perform multiple testing correction on MAGMA gene set analysis output (.gsa.out)
    and generate a cleaned, tab-delimited corrected results file.

    Parameters
    ----------
    magma_analysis_folder : str
        Folder containing MAGMA .gsa.out output file
    sample_id : str
        Identifier prefix used for the analysis
    output_file : str, optional
        Output .tsv file path (default: <magma_analysis_folder>/<sample_id>.gsa_corrected.tsv)
    Returns
    -------
    pandas.DataFrame
        DataFrame with corrected p-values (Bonferroni, Sidak, Holm, FDR, etc.)
    """
    folder = Path(magma_analysis_folder)
    folder.mkdir(parents=True, exist_ok=True)
    if output_file is None:
        output_file = folder / f"{sample_id}.gsa_corrected.tsv"
    else:
        output_file = Path(output_file)
    gsa_out = folder / f"{sample_id}.gsa.out"
    gsa_tsv = folder / f"{sample_id}.gsa.tsv"
    # ---------- Step 1: Clean up .gsa.out file ----------
    cmd = f"cat {gsa_out} | tr -s ' ' | tr ' ' '\\t' | grep -v '#' > {gsa_tsv}"
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    # ---------- Step 2: Load cleaned data ----------
    df = pd.read_csv(gsa_tsv, sep="\t", comment="#", engine="python")
    if "P" not in df.columns:
        raise ValueError("⚠️ Missing column: 'P' not found in file.")
    if "VARIABLE" not in df.columns:
        raise ValueError("⚠️ Missing column: 'VARIABLE' not found in file.")
    n_tests = len(df)
    P = df["P"].astype(float).values
    # ---------- Step 3: Apply global multiple-testing corrections ----------
    df["P_bonferroni_corr"] = np.minimum(P * n_tests, 1.0)
    df["P_sidak_corr"] = 1 - (1 - P) ** n_tests
    df["P_holm_corr"] = multipletests(P, method="holm")[1]
    df["P_fdr_bh_corr"] = multipletests(P, method="fdr_bh")[1]
    # ---------- Step 4: Category-specific corrections ----------
    def apply_subset_correction(df, pattern, prefix):
        subset = df[df["VARIABLE"].str.match(pattern, na=False)].copy()
        if len(subset) == 0:
            return pd.DataFrame(columns=["VARIABLE", f"{prefix}_P_bonferroni_corr", f"{prefix}_P_fdr_bh_corr"])
        pvals = subset["P"].astype(float).values
        subset[f"{prefix}_P_bonferroni_corr"] = np.minimum(pvals * len(pvals), 1.0)
        subset[f"{prefix}_P_fdr_bh_corr"] = multipletests(pvals, method="fdr_bh")[1]
        return subset[["VARIABLE", f"{prefix}_P_bonferroni_corr", f"{prefix}_P_fdr_bh_corr"]]
    go_kegg_reactome = apply_subset_correction(df, r"^(GOBP_|GOCC_|GOMF_|KEGG_|REACTOME_)", "go_kegg_reactome")
    go_df = apply_subset_correction(df, r"^(GOBP_|GOCC_|GOMF_)", "go")
    reactome_df = apply_subset_correction(df, r"^REACTOME_", "reactome")
    kegg_df = apply_subset_correction(df, r"^KEGG_", "kegg")
    # ---------- Step 5: Merge all corrected data ----------
    merged_df = df.copy()
    for sub_df in [go_kegg_reactome, go_df, reactome_df, kegg_df]:
        merged_df = merged_df.merge(sub_df, on="VARIABLE", how="outer")
    merged_df = merged_df.sort_values(by="P", ascending=True)
    # ---------- Step 6: Save results ----------
    merged_df.to_csv(output_file, sep="\t", index=False)
    print(f"✅ Multiple testing correction completed successfully!")
    print(f"📄 Results saved to: {output_file}")
    return merged_df




def add_gene_details(gene_string: str, gsa_genes_set: set, gsa_genes_df: pl.DataFrame):
    """
    Compute overlap between an input gene list and MAGMA result genes, returning
    common genes, their p-values, and summary counts.
    Parameters
    ----------
    gene_string : str
        Comma-separated gene symbols (e.g., "DISC1,NRG1,GRIN2B").
    gsa_genes_set : set
        Set of gene names found in MAGMA results (for fast intersection).
    gsa_genes_df : pl.DataFrame
        MAGMA gene results dataframe with at least two columns: "GENE" and "P".
    Returns
    -------
    dict
        {
          "CommonGenes": comma-separated gene names (or None),
          "CommonGenes_Pvalues": comma-separated p-values (or None),
          "TotalGenes": total count of genes in input list,
          "CommonGeneCount": count of overlapping genes
        }
    """
    # Handle empty input
    if not isinstance(gene_string, str) or gene_string.strip() == "":
        return {
            "CommonGenes": None,
            "CommonGenes_Pvalues": None,
            "TotalGenes": 0,
            "CommonGeneCount": 0,
        }
    # Split comma-separated genes and clean whitespace
    gene_list = [g.strip() for g in gene_string.split(",") if g.strip()]
    total_genes = len(gene_list)
    # Find intersection with available genes
    common_genes = [g for g in gene_list if g in gsa_genes_set]
    common_gene_count = len(common_genes)
    if common_gene_count == 0:
        return {
            "CommonGenes": None,
            "CommonGenes_Pvalues": None,
            "TotalGenes": total_genes,
            "CommonGeneCount": 0,
        }
    # Filter gsa_genes_df to keep only the intersecting genes
    temp_df = (
        gsa_genes_df.filter(pl.col("GENE").is_in(common_genes))
        .sort("GENE")
        .select(["GENE", "P"])
    )
    # Convert columns to Python lists
    genes = temp_df["GENE"].to_list()
    pvalues = [str(p) for p in temp_df["P"].to_list()]
    # Join into comma-separated strings
    return {
        "CommonGenes": ",".join(genes),
        "CommonGenes_Pvalues": ",".join(pvalues),
        "TotalGenes": total_genes,
        "CommonGeneCount": common_gene_count,
    }




def process_gene_sets(magma_analysis_folder: str, sample_id: str, geneset_file: str, corrected_df: pl.DataFrame):
    """
    Process MAGMA gene sets and merge with corrected enrichment results.

    Parameters
    ----------
    magma_analysis_folder : str
        Folder containing MAGMA outputs (.genes.out, .gsa_corrected.tsv)
    sample_id : str
        Sample identifier (prefix)
    geneset_file : str
        Path to MSigDB .gmt or .tsv file
    corrected_df : pl.DataFrame
        DataFrame containing corrected MAGMA gene-set results

    Returns
    -------
    pl.DataFrame
        Final merged dataframe with gene overlaps and corrected enrichment info
    """
    folder = Path(magma_analysis_folder)
    folder.mkdir(parents=True, exist_ok=True)
    # ---------- 1️⃣ Clean MAGMA gene-level output ----------
    genes_out_file = folder / f"{sample_id}_magma_35up_10down.genes.out"
    temp_file = folder / "temp_file"
    cmd = f"cat {genes_out_file} | tr -s ' ' | grep -v '^#' > {temp_file} && mv {temp_file} {genes_out_file}"
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    # ---------- 2️⃣ Read MAGMA gene-level results ----------
    gsa_genes_df = pl.read_csv(genes_out_file, separator=" ", has_header=True)
    if "GENE" not in gsa_genes_df.columns or "P" not in gsa_genes_df.columns:
        raise ValueError("⚠️ Required columns (GENE, P) not found in MAGMA .genes.out file.")
    gsa_genes_set = set(gsa_genes_df["GENE"].to_list())
    # ---------- 3️⃣ Read geneset file ----------
    geneset_file = Path(geneset_file)
    if geneset_file.suffix == ".gmt":
        genesetfile_tsv = geneset_file.with_suffix(".tsv")
        # Convert .gmt → .tsv (simple copy rename for compatibility)
        subprocess.run(f"cp {geneset_file} {genesetfile_tsv}", shell=True)
    else:
        genesetfile_tsv = geneset_file
    msigdb_df = pl.read_csv(genesetfile_tsv, separator="\t", has_header=False)
    msigdb_df.columns = ["V1", "V2", "V3"]  # V1=name, V2=desc, V3=genes (comma-separated)
    # ---------- 4️⃣ Parallel processing with add_gene_details ----------
    n_cores = max(1, cpu_count() - 1)
    print(f"🧠 Using {n_cores} cores for parallel gene overlap computation...")
    def process_row(row):
        return add_gene_details(row["V3"], gsa_genes_set, gsa_genes_df)
    # Convert to list of dicts for multiprocessing
    rows = msigdb_df.iter_rows(named=True)
    with Pool(processes=n_cores) as pool:
        results = pool.map(partial(process_row), rows)
    # Combine results back into Polars DataFrame
    result_df = pl.DataFrame(results)
    msigdb_combined = pl.concat([msigdb_df, result_df], how="horizontal")
    # ---------- 5️⃣ Merge with corrected_df ----------
    if "FULL_NAME" not in corrected_df.columns:
        raise ValueError("⚠️ 'FULL_NAME' column not found in corrected_df for merging.")
    merged_df = corrected_df.join(
        msigdb_combined,
        left_on="FULL_NAME",
        right_on="V1",
        how="left"
    ).sort("P")
    # ---------- 6️⃣ Add sample ID and save ----------
    merged_df = merged_df.with_columns(pl.lit(sample_id).alias("sample_id"))
    output_file = folder / f"{sample_id}.gsa_corrected.tsv"
    merged_df.write_csv(output_file, separator="\t")
    print(f"✅ Process completed successfully! Output saved to: {output_file}")
    return merged_df




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
    