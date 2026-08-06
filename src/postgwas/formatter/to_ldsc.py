import subprocess
from pathlib import Path
import polars as pl
import textwrap
import time


## reason to use effective sample size and sample prevalence 0.5 ; https://groups.google.com/g/ldsc_users/c/yJT-_qSh_44/m/MmKKJYsBAwAJ


## reason to use effective sample size and sample prevalence 0.5

# =========================================================
# LOGGER
# =========================================================
def log_message(log_file: Path, message: str):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    msg = f"[{ts}] {message}"
    print(msg)
    with open(log_file, "a") as lf:
        lf.write(msg + "\n")


# =========================================================
# STEP 1: DROP / FIX N_CAS / N_CON
# ========================================================
def drop_if_mostly_empty(input_file: str, log_file: str, threshold=0.8):
    input_path = Path(input_file)
    log_path = Path(log_file)
    log_message(log_path, f"[STEP] Checking N_CAS / N_CON in: {input_path}")
    # Read the data
    df = pl.read_csv(input_path, separator="\t")
    total_rows = df.height
    cols = ["N_CAS", "N_CON"]
    present_cols = []
    sample_prev = None 
    # ---------------------------
    # 🛠️ CONVERT & DROP IMMEDIATELY
    # ---------------------------
    for c in cols:
        if c not in df.columns:
            log_message(log_path, f"[SKIP] {c} not present")
            continue
        # 1. Cast to Float: Dots/strings become Nulls
        df = df.with_columns(pl.col(c).cast(pl.Float64, strict=False))
        # 2. Check empty fraction
        empty_count = df.select(pl.col(c).is_null().sum()).item()
        empty_fraction = empty_count / total_rows
        # 3. DROP from df immediately if mostly garbage
        if empty_fraction >= threshold:
            df = df.drop(c)
            log_message(log_path, f"[DROP] {c} is {empty_fraction:.2%} invalid -> removed from dataframe")
        else:
            present_cols.append(c)
            log_message(log_path, f"[KEEP] {c}: {empty_count}/{total_rows} ({empty_fraction:.2%}) empty")
    # ---------------------------
    # 🧬 SAMPLE PREVALENCE (Only if BOTH survived)
    # ---------------------------
    if len(present_cols) == 2:
        try:
            # We use median because UKB sample sizes can vary slightly per SNP
            n_cases = df.select(pl.col("N_CAS").median()).item()
            n_controls = df.select(pl.col("N_CON").median()).item()
            if n_cases and n_controls and (n_cases + n_controls) > 0:
                sample_prev = n_cases / (n_cases + n_controls)
                log_message(log_path, f"[SPREV] Computed prevalence = {sample_prev:.4f}")
        except Exception as e:
            log_message(log_path, f"[SPREV] Error: {e}")
    # ---------------------------
    # ⚖️ SIMPLIFIED RENAMING LOGIC
    # ---------------------------
    # If only one column survived (like N_CON), rename it to total N
    if len(present_cols) == 1:
        c = present_cols[0]
        log_message(log_path, f"[FIX] Only {c} is valid. Renaming to 'N'")
        df = df.rename({c: "N"})
    elif len(present_cols) == 0:
        log_message(log_path, "[WARNING] No valid N columns remained after filtering.")
    # ---------------------------
    # SAVE
    # ---------------------------
    df.write_csv(input_path, separator="\t")
    log_message(log_path, "[DONE] N column check complete")
    return str(input_path), sample_prev

# =========================================================
# STEP 2: VCF → LDSC
# =========================================================
def vcf_to_ldssc_munge_input(sumstat_vcf: str, output_folder: str, sample_name: str):

    vcf_path = Path(sumstat_vcf)
    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"{sample_name}_ldsc_input.tsv"
    log_file = output_dir / f"{sample_name}_vcf_to_ldsc.log"

    cmd = textwrap.dedent(f"""
        {{
            printf "SNP\\tA2\\tA1\\tZ\\tP\\tN_CAS\\tFRQ\\tINFO\\tN_CON\\n"
            bcftools view --min-alleles 2 --max-alleles 2 "{vcf_path}" | \
            bcftools query -f '%ID\\t%REF\\t%ALT\\t[%EZ]\\t[%LP]\\t[%NC]\\t[%AF]\\t[%SI]\\t[%NCO]\\n' | \
            sed 's|:|_|g' | \
            awk -F '\\t' 'BEGIN {{ OFS="\\t" }} {{
                raw_p = exp(-$5 * log(10))
                $5 = sprintf("%.6g", raw_p)
                print
            }}'
        }} > "{output_file}"
    """)

    # write command
    with open(log_file, "w") as lf:
        lf.write("### COMMAND\n\n")
        lf.write(cmd + "\n\n")

    log_message(log_file, "[STEP] Running VCF → LDSC")

    try:
        with open(log_file, "a") as lf:
            subprocess.run(
                cmd,
                shell=True,
                executable="/bin/bash",
                check=True,
                stdout=lf,
                stderr=lf,
            )
    except subprocess.CalledProcessError:
        log_message(log_file, "❌ Conversion failed")
        raise RuntimeError(f"vcf_to_ldsc failed → {log_file}")

    log_message(log_file, "[DONE] VCF → LDSC complete")

    return {
        "ldsc_file": str(output_file),
        "log_file": str(log_file),
    }


# =========================================================
# PIPELINE WRAPPER (UNCHANGED STRUCTURE)
# =========================================================
def vcf_to_ldsc(sumstat_vcf, output_folder, sample_name):

    result = vcf_to_ldssc_munge_input(sumstat_vcf, output_folder, sample_name)

    result["ldsc_file"], sample_prev=drop_if_mostly_empty( 
        input_file=result["ldsc_file"],
        log_file=result["log_file"],
        threshold=0.8
    )

    return {
        "ldsc_file": str(result["ldsc_file"]),
        "log_file": str(result["log_file"]),
        "sample_prev": sample_prev,
    }