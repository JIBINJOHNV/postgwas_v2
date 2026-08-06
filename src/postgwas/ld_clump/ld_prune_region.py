
import numpy as np
import polars as pl
from pathlib import Path
import subprocess
import os
import datetime
import logging




def ld_clump_by_regions(
    sumstat_vcf: str,
    output_folder: str,
    sample_name: str,
    population: str = "EUR",
    bcftools: str = "bcftools",
    nthreads=4
):
    """
    Extract summary statistics from a VCF and perform LD-block-based pruning.
    
    Changes:
    - FIX: Uses colon ':' separator in bcftools to prevent 'REF_' tag error.
    - FIX: Replaces ':' with '_' in Polars for safe ID generation.
    - Logs command and errors to <sample_name>_ld_clump.log.
    """
    population = population.upper()
    if population not in {"EUR", "AFR", "EAS"}:
        raise ValueError("Population must be one of: 'EUR', 'AFR', or 'EAS'")
    
    vcf_path = Path(sumstat_vcf)
    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Define filenames
    raw_out = output_dir / f"{sample_name}_vcf.tsv.gz"
    pruned_out = output_dir / f"{sample_name}_LDpruned_{population}.tsv"
    pruned_sig_out = output_dir / f"{sample_name}_LDpruned_{population}_sig.tsv"
    log_file = output_dir / f"{sample_name}_ld_clump.log"

    # --------------------- 1️⃣ Extract INFO fields (with Logging) ---------------------
    
    # FIX: Use ':' separators here. bcftools handles %REF: correctly. 
    # We will convert ':' to '_' in the Polars step below.
    cmd = f"""
    {{
        printf "SNP\\tCHR\\tBP\\tREF\\tALT\\tBETA\\tSE\\tNEF\\tAF\\tAFR\\tEAS\\tEUR\\tSAS\\tLP\\tNC\\tNCO\\tEUR_LDblock\\tAFR_LDblock\\tEAS_LDblock\\tBCSQ\\n"
        {bcftools} view --threads {nthreads} --min-alleles 2 --max-alleles 2 "{vcf_path}" | \
        {bcftools} query -f '%CHROM:%POS:%REF:%ALT\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t[%ES]\\t[%SE]\\t[%NEF]\\t[%AF]\\t[%AFR]\\t[%EAS]\\t[%EUR]\\t[%SAS]\\t[%LP]\\t[%NC]\\t[%NCO]\\t%INFO/EUR_LDblock\\t%INFO/AFR_LDblock\\t%INFO/EAS_LDblock\\t[%BCSQ]\\n'
    }} | gzip -c > "{raw_out}"
    """

    with open(log_file, "w") as log:
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log.write(f"[{timestamp}] STARTING LD CLUMPING\n")
        log.write(f"[{timestamp}] COMMAND:\n{cmd}\n\n")
        log.write(f"[{timestamp}] BCFTOOLS OUTPUT (STDERR):\n")
        log.flush()

        try:
            subprocess.run(
                cmd, 
                shell=True, 
                check=True, 
                executable="/bin/bash", 
                stderr=log, 
                stdout=log 
            )
            log.write(f"\n[{timestamp}] ✅ Extraction successful -> {raw_out}\n")
            
        except subprocess.CalledProcessError as e:
            err_msg = f"❌ bcftools extraction failed with exit code {e.returncode}. Check log: {log_file}"
            print(err_msg)
            log.write(f"\n[{timestamp}] {err_msg}\n")
            return None

    # --------------------- 2️⃣ Load with Polars ---------------------
    try:
        df = pl.read_csv(raw_out, separator="\t", null_values=[".", "NA"], infer_schema_length=20000)
    except Exception as e:
        with open(log_file, "a") as log:
            log.write(f"\n❌ Error reading TSV file: {e}\n")
        print(f"❌ Error reading extracted TSV: {e}")
        return None

    # --------------------- 3️⃣ Prepare Columns ---------------------
    ld_field = f"{population}_LDblock"
    if ld_field not in df.columns:
        msg = f"Expected LD block column '{ld_field}' not found in file."
        with open(log_file, "a") as log: log.write(f"\n❌ {msg}\n")
        raise ValueError(msg)

    # Filter empty LD blocks first
    df = df.filter(pl.col(ld_field).is_not_null())

    if df.height == 0:
        msg = "\t\t\t\t  ⚠️ No variants found with LD block annotations. Check population match."
        print(msg)
        with open(log_file, "a") as log: log.write(f"\n{msg}\n")
        return None

    df = (
        df.with_columns([
            # FIX: Convert ':' to '_' to finalize ID format (e.g., 1:100:A:T -> 1_100_A_T)
            pl.col("SNP").str.replace_all(":", "_"),
            pl.col(ld_field).alias("LDblock"),
            pl.col("LP").cast(pl.Float64, strict=False),
        ])
        .filter(pl.col("LP").is_not_null())
        .with_columns([
            (10 ** (-pl.col("LP"))).alias("P_value"),
            pl.col("LDblock").str.split("_").alias("LD_parts")
        ])
        .with_columns([
            pl.col("LD_parts").list.get(1).cast(pl.Int64, strict=False).alias("START"),
            pl.col("LD_parts").list.get(2).cast(pl.Int64, strict=False).alias("END"),
        ])
        .drop("LD_parts")
    )

    # --------------------- 4️⃣ LD-block-based pruning ---------------------
    pruned = (
        df.sort("LP", descending=True)
          .group_by("LDblock", maintain_order=True)
          .head(1)
    )

    # --------------------- 5️⃣ Genome-wide significant subset ---------------------
    pruned_sig = pruned.filter(pl.col("LP") >= 7.3)

    # --------------------- 6️⃣ Save outputs ---------------------
    pruned.write_csv(pruned_out, separator="\t")
    pruned_sig.write_csv(pruned_sig_out, separator="\t")

    with open(log_file, "a") as log:
        log.write(f"✅ LD-block pruned file: {pruned_out} (Rows: {pruned.height})\n")
        log.write(f"✅ Significant SNPs file: {pruned_sig_out} (Rows: {pruned_sig.height})\n")

    return {
        "ldpruned_sig_file": str(pruned_sig_out),
        "ldpruned_file": str(pruned_out),
        "log_file": str(log_file)
    }
