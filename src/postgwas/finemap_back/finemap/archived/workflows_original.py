import polars as pl
import subprocess
import shutil
import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp



# ==============================================================
# 1. LOAD + PREPARE SUMSTATS
# ==============================================================
def load_and_prepare_finemap_sumstats(
    sig_loci_path: str,
    sumstats_path: str,
    finemap_ld_ref_bim: str,
):
    print(f"   -> Reading sig loci: {sig_loci_path}")
    sig_loci_df = pl.read_csv(sig_loci_path, separator="\t")
    
    print(f"   -> Reading sumstats: {sumstats_path}")
    sumstats_df = pl.read_csv(sumstats_path, separator="\t")
    
    print(f"   -> Reading LD reference BIM: {finemap_ld_ref_bim}")
    ld_ref_df = pl.read_csv(
        finemap_ld_ref_bim,
        separator="\t",
        has_header=False,
        columns=[1],
        new_columns=["variant_id"],
    )

    # AF → MAF transformation (if not already MAF)
    # Ensure we don't have nulls before comparison
    sumstats_df = sumstats_df.with_columns(
        pl.when(pl.col("maf") <= 0.5)
        .then(pl.col("maf"))
        .otherwise(1.0 - pl.col("maf"))
        .alias("maf")
    )

    # Deterministic order
    sumstats_df = sumstats_df.sort(
        by=["chromosome", "position", "allele1", "allele2"]
    )

    # Filter to LD reference intersection
    print("   -> Intersecting with Reference Panel...")
    sumstats_filt = sumstats_df.join(
        ld_ref_df,
        left_on="rsid",
        right_on="variant_id",
        how="inner", # 'inner' is usually safer than 'semi' to ensure data integrity
    )
    
    return sig_loci_df, sumstats_filt




# ==============================================================
# 2. PREPARE FINEMAP INPUTS FOR ONE LOCUS
# ==============================================================
def prepare_finemap_ld_for_locus(
    locus_sumstats: pl.DataFrame,
    locus_id: str,
    outdir: str,
    plink_ref_prefix: str,
    n_threads: int = 4,
):
    if locus_sumstats.is_empty():
        raise ValueError(f"{locus_id}: empty locus")
    if "NEF" not in locus_sumstats.columns:
        raise ValueError(f"{locus_id}: NEF missing")
    
    # Check for bgenix
    if shutil.which("bgenix") is None:
        raise EnvironmentError("❌ 'bgenix' tool not found! LDstore requires BGEN files to be indexed. Please install 'bgenix' (often part of the BGEN suite).")

    # 1. Use Absolute Paths
    outdir_path = Path(outdir).resolve()
    outdir_path.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------
    # File paths
    # --------------------------------------------------
    sumstats_file = outdir_path / f"{locus_id}.z"
    snp_file      = outdir_path / f"{locus_id}.snp"
    plink_out     = outdir_path / f"{locus_id}_plink"
    
    # PLINK2 outputs
    bgen_prefix   = outdir_path / f"{locus_id}_genotypes"
    bgen_file     = outdir_path / f"{locus_id}_genotypes.bgen"
    # Note: bgenix usually creates <file>.bgen.bgi
    bgi_file      = outdir_path / f"{locus_id}_genotypes.bgen.bgi" 
    
    # LDstore files
    bcor_file     = outdir_path / f"{locus_id}.bcor"
    ld_matrix     = outdir_path / f"{locus_id}.ld"
    ldstore_master = outdir_path / f"{locus_id}.ldstore.master"

    # FINEMAP files
    config_file   = outdir_path / f"{locus_id}.config"
    cred_file     = outdir_path / f"{locus_id}.cred"
    log_file      = outdir_path / f"{locus_id}.log"
    master_file   = outdir_path / f"{locus_id}.master"

    # --------------------------------------------------
    # 2. Write sumstats
    # --------------------------------------------------
    cols = ["rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se"]
    locus_sumstats.select(cols).write_csv(sumstats_file, separator=" ")

    n_samples = int(
        locus_sumstats
        .select(pl.col("NEF").drop_nulls().median())
        .item()
    )

    # --------------------------------------------------
    # 3. SNP list & PLINK -> BGEN
    # --------------------------------------------------
    locus_sumstats.select("rsid").write_csv(snp_file, has_header=False)

    # Extract SNPs to BED
    subprocess.run(
        ["plink2", "--bfile", plink_ref_prefix,
         "--extract", str(snp_file),
         "--make-bed", "--out", str(plink_out), 
         "--silent"],
        check=True, stdout=subprocess.DEVNULL
    )

    # Export to BGEN (8-bit)
    subprocess.run(
        ["plink2", "--bfile", str(plink_out),
         "--export", "bgen-1.2", "bits=8",
         "--out", str(bgen_prefix), 
         "--silent"],
        check=True, stdout=subprocess.DEVNULL
    )

    # --------------------------------------------------
    # 4. Generate BGEN Index (Required for LDstore)
    # --------------------------------------------------
    # bgenix -g <file.bgen> -index -clobber
    subprocess.run(
        ["bgenix", "-g", str(bgen_file), "-index", "-clobber"],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )

    # --------------------------------------------------
    # 5. LDstore Master File
    # --------------------------------------------------
    n_ld_samples =   subprocess.run(
        ["wc", "-l", f"{plink_ref_prefix}.fam"],
        capture_output=True, text=True, check=True,)
    
    n_ld_samples = int(n_ld_samples.stdout.strip().split()[0])

    # Assuming same sample size for LD
    # UPDATED HEADER: z;bgen;bgi;bcor;ld;sample
    with ldstore_master.open("w") as f:
        f.write("z;bgen;bgi;bcor;ld;n_samples\n")
        f.write(
            f"{sumstats_file};"
            f"{bgen_file};"
            f"{bgi_file};"
            f"{bcor_file};"
            f"{ld_matrix};"
            f"{n_ld_samples}\n"
        )

    # --------------------------------------------------
    # 6. Run LDstore Steps
    # --------------------------------------------------
    # Step A: BGEN -> BCOR (Multi-threaded)
    subprocess.run(
        [
            "ldstore",
            "--in-files", str(ldstore_master),
            "--read-only-bgen",
            "--write-bcor",
            "--n-threads", str(n_threads),
        ],
        check=True, stdout=subprocess.DEVNULL
    )

    # Step B: BCOR -> LD Matrix (Single-threaded only)
    subprocess.run(
        [
            "ldstore",
            "--in-files", str(ldstore_master),
            "--bcor-to-text",
        ],
        check=True, stdout=subprocess.DEVNULL
    )

    # --------------------------------------------------
    # 7. FINEMAP Master File
    # --------------------------------------------------
    for f in (config_file, cred_file, log_file):
        if f.exists(): f.unlink()
        f.touch()

    with master_file.open("w") as f:
        f.write("z;ld;snp;config;cred;log;n_samples\n")
        f.write(
            f"{sumstats_file};{ld_matrix};{snp_file};"
            f"{config_file};{cred_file};{log_file};{n_samples}\n"
        )

    # --------------------------------------------------
    # Step 8: Run FINEMAP
    # --------------------------------------------------
    print(f"   -> Running FINEMAP for {locus_id}...")
        # inside the master_file to determine where to write the log.
    run_res = subprocess.run(
        [ "finemap", 
            "--sss", 
            "--in-files", str(master_file), 
            "--n-threads", str(n_threads)
        ],
        capture_output=True,
        text=True
    )
    
    # "--log-file", str(locus_id + "_finemap_run.log")

    # Check for success
    if run_res.returncode != 0:
        raise RuntimeError(f"FINEMAP failed:\n{run_res.stderr}\n{run_res.stdout}")
    
    return f"✅ {locus_id} done"




# ==============================================================
# 3. WORKER FUNCTION (pickle-safe)
# ==============================================================

def run_one_locus_finemap(args):
    locus, sumstats_filt, outdir, plink_ref_prefix, n_threads = args
    
    # Access Polars Struct or Dict depending on how it was iterated
    chrom = locus["CHR"] # Adjust column name if your sig_loci file differs
    start = locus["START"]
    end   = locus["END"]
    
    locus_id = f"chr{chrom}_{start}_{end}"
    
    # Filter Sumstats for this window
    locus_sumstats = sumstats_filt.filter(
        (pl.col("chromosome") == chrom) &
        (pl.col("position") >= start) &
        (pl.col("position") <= end)
    )
    
    if locus_sumstats.height < 10:
        return f"⚠️ {locus_id}: skipped (only {locus_sumstats.height} variants found)"
        
    try:
        return prepare_finemap_ld_for_locus(
            locus_sumstats,
            locus_id,
            outdir,
            plink_ref_prefix,
            n_threads,
        )
    except Exception as e:
        import traceback
        return f"❌ {locus_id} failed: {e}\n{traceback.format_exc()}"


# ==============================================================
# 4. MULTI-LOCUS PARALLEL DRIVER
# ==============================================================

def run_finemap_all_loci(
    sig_loci_df,
    sumstats_filt,
    outdir,
    plink_ref_prefix,
    max_workers=4,
    n_threads_ldstore=4,
):
    loci = list(sig_loci_df.iter_rows(named=True))
    
    # Check for required tools
    for tool in ["plink", "plink2", "ldstore"]:
        if shutil.which(tool) is None:
            raise EnvironmentError(f"Tool '{tool}' not found in PATH.")

    tasks = [
        (locus, sumstats_filt, outdir, plink_ref_prefix, n_threads_ldstore)
        for locus in loci
    ]
    
    print(f"🚀 Starting processing for {len(tasks)} loci with {max_workers} workers...")
    
    with ProcessPoolExecutor(
        max_workers=max_workers,
        mp_context=mp.get_context("spawn"),
    ) as executor:
        futures = [executor.submit(run_one_locus_finemap, t) for t in tasks]
        
        for i, fut in enumerate(as_completed(futures)):
            result = fut.result()
            print(f"[{i+1}/{len(tasks)}] {result}")


def run_finemap_pipeline(
    sig_loci_path: str,
    sumstats_path: str,
    finemap_ld_ref_bim: str,
    outdir: str,
    max_workers: int = 4,
    n_threads_ldstore: int = 4,
):
    print("📥 Loading and preparing summary statistics...")
    sig_loci_df, sumstats_filt = load_and_prepare_finemap_sumstats(
        sig_loci_path=sig_loci_path,
        sumstats_path=sumstats_path,
        finemap_ld_ref_bim=finemap_ld_ref_bim,
    )
    
    print(f"🧬 Significant loci loaded: {sig_loci_df.height}")
    print(f"📊 Variants after LD filter: {sumstats_filt.height}")
    
    run_finemap_all_loci(
        sig_loci_df=sig_loci_df,
        sumstats_filt=sumstats_filt,
        outdir=outdir,
        plink_ref_prefix=finemap_ld_ref_bim.replace(".bim", ""),
        max_workers=max_workers,
        n_threads_ldstore=n_threads_ldstore,
    )
    print("✅ FINEMAP pipeline completed.")


# ==============================================================
# 5. ENTRY POINT
# ==============================================================

if __name__ == "__main__":
    # Ensure this runs in a `if __name__` block for multiprocessing
    # run_finemap_pipeline(
    #     sig_loci_path="/Users/JJOHN41/Documents/developing_software/data/oudir/07_ld_clump/PGC3_SCZ_european_LDpruned_EUR_sig.tsv",
    #     sumstats_path="/Users/JJOHN41/Documents/developing_software/data/oudir/PGC3_SCZ_european/06_formatter/PGC3_SCZ_european_finemap_finemap.tsv",
    #     finemap_ld_ref_bim="/Users/JJOHN41/Documents/software_resources/resourses/postgwas/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001.bim",
    #     outdir="/Users/JJOHN41/Documents/developing_software/data/oudir/08_finemap/loci",
    #     max_workers=4,
    #     n_threads_ldstore=4,
    # )