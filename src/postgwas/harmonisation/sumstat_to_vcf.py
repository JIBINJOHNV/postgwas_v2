from postgwas.harmonisation import gwastovcf_gwas2vcf
from pathlib import Path
import shutil
import subprocess
import io
import re
from concurrent.futures import ThreadPoolExecutor


# ============================================================
# 1. GWAS → VCF (bcftools munge)
# ============================================================
def run_bcftools_munge(
    output_dir: str,
    gwas_outputname: str,
    chr: str,
    genome_fasta_file: str,
    grch_version: str,
) -> str:

    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    input_file = outdir / f"{gwas_outputname}_chr{chr}_vcf_input.tsv"
    column_file = outdir / f"{gwas_outputname}_chr{chr}_column_mapping.tsv"
    output_vcf = outdir / f"{gwas_outputname}_chr{chr}_{grch_version}.vcf.gz"

    if not input_file.exists():
        raise FileNotFoundError(f"[chr{chr}] Missing input: {input_file}")

    cmd = [
        "bash", "-c",
        f"""
        bcftools +munge --no-version \
            --columns-file "{column_file}" \
            --fasta-ref "{genome_fasta_file}" \
            --sample-name "{gwas_outputname}" \
            "{input_file}" |
        bcftools norm -m-any -d exact |
        bcftools sort -Oz -o "{output_vcf}" --write-index=tbi
        """
    ]

    subprocess.run(cmd, check=True)
    return str(output_vcf)


# ============================================================
# HELPER: Run bcftools step
# ============================================================
def _run_bcftools_step(cmd, log_file, step_name, chromosome):
    with open(log_file, "a") as lf:
        lf.write(f"\n[chr{chromosome}] === {step_name} ===\n")
        result = subprocess.run(cmd, stdout=lf, stderr=lf, text=True)

        if result.returncode != 0:
            raise RuntimeError(
                f"[chr{chromosome}] ❌ STEP FAILED: {step_name}\nCheck log: {log_file}"
            )


# ============================================================
# HELPER: Variant count (safe)
# ============================================================
def get_vcf_variant_count(vcf_path: str):
    result = subprocess.run(
        ["bcftools", "index", "-n", vcf_path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode == 0:
        return int(result.stdout.strip())

    return None  # safe fallback


# ============================================================
# HELPER: Logging stats
# ============================================================
def log_variant_stats(log, step, current, previous, original, chromosome):

    prefix = f"[chr{chromosome}]"

    if current is None:
        log(f"{prefix} ⚠️ {step}: variant count unavailable")
        return

    msg = f"{prefix} 📊 {step}: {current:,}"

    if previous:
        drop_prev = previous - current
        pct_prev = (drop_prev / previous) * 100
        msg += f" | Δprev: -{drop_prev:,} ({pct_prev:.2f}%)"

    if original:
        drop_orig = original - current
        pct_orig = (drop_orig / original) * 100
        msg += f" | Δorig: -{drop_orig:,} ({pct_orig:.2f}%)"

    log(msg)

    if previous and current < previous * 0.8:
        log(f"{prefix} 🚨 WARNING: >20% variant drop at {step}")


# ============================================================
# 2. ANNOTATION + LIFTOVER PIPELINE
# ============================================================
def run_bcftools_annot(
    output_dir: str,
    gwas_outputname: str,
    chromosome: str,
    external_eaf_file: str,
    default_dbsnp_file: str,
    genome_fasta_file: str,
    target_genome_fasta_file: str,
    gff_file: str,
    chain_file: str,
    grch_version: str,
    threads: int = 5,
) -> str:

    from pathlib import Path
    import shutil
    import tempfile

    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    log_dir = outdir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_bcftools_annot.log"

    def log(msg):
        with open(log_file, "a") as lf:
            lf.write(msg + "\n")

    input_vcf = outdir / f"{gwas_outputname}_chr{chromosome}_{grch_version}.vcf.gz"

    if not input_vcf.exists():
        raise FileNotFoundError(f"[chr{chromosome}] Missing VCF: {input_vcf}")

    original_vcf = outdir / f"{gwas_outputname}_chr{chromosome}_{grch_version}_ORIGINAL.vcf.gz"

    # =======================================================
    # SAFE COPY (VCF + INDEX)
    # =======================================================
    if not original_vcf.exists():
        shutil.copy2(input_vcf, original_vcf)
        if Path(f"{input_vcf}.tbi").exists():
            shutil.copy2(f"{input_vcf}.tbi", f"{original_vcf}.tbi")

    norm_vcf = outdir / f"{gwas_outputname}_chr{chromosome}_{grch_version}_NORM.vcf.gz"
    id_vcf   = outdir / f"{gwas_outputname}_chr{chromosome}_{grch_version}_ID.vcf.gz"
    af_vcf   = outdir / f"{gwas_outputname}_chr{chromosome}_{grch_version}_EAF.vcf.gz"
    csq_vcf  = outdir / f"{gwas_outputname}_chr{chromosome}_{grch_version}_CSQ.vcf.gz"

    target_build = "GRCh38" if grch_version == "GRCh37" else "GRCh37"
    target_vcf = outdir / f"{gwas_outputname}_chr{chromosome}_{target_build}_CSQ.vcf.gz"
    reject_vcf = outdir / f"{gwas_outputname}_chr{chromosome}_{target_build}_CSQ_notlifted.vcf.gz"

    counts = {}
    sort_tmp_dir = Path(output_dir) / f"temp_chr{chromosome}"
    sort_tmp_dir.mkdir(parents=True, exist_ok=True)

    # =======================================================
    # ORIGINAL COUNT
    # =======================================================
    counts["ORIGINAL"] = get_vcf_variant_count(str(original_vcf))
    log_variant_stats(log, "ORIGINAL", counts["ORIGINAL"], None, counts["ORIGINAL"], chromosome)

    # =======================================================
    # STEP 0: NORMALIZATION 
    # =======================================================
    _run_bcftools_step([
        "bash", "-c",
        f"""
        set -euo pipefail
        bcftools norm --threads {threads} -Ou -m-any -d exact --fasta-ref "{genome_fasta_file}" "{input_vcf}" |
        bcftools annotate --threads {threads} -x ID -Oz -o {norm_vcf} --write-index=tbi
        """
    ], log_file, "STEP0_NORM", chromosome)

    counts["NORM"] = get_vcf_variant_count(str(norm_vcf))
    log_variant_stats(log, "NORM", counts["NORM"], counts["ORIGINAL"], counts["ORIGINAL"], chromosome)

    # =======================================================
    # STEP 1: ID ANNOTATION 
    # =======================================================
    _run_bcftools_step([ 
        "bash", "-c",
        f"""
        set -euo pipefail
        bcftools annotate --threads {threads} \
            --annotations "{default_dbsnp_file}" \
            --columns ID --pair-logic exact \
            "{norm_vcf}" |
        bcftools annotate --threads {threads} \
            --set-id +'%CHROM\\_%POS\\_%REF\\_%ALT' |
        bcftools view --threads {threads} \
            -Oz -o "{id_vcf}" \
            --write-index=tbi
        """
    ], log_file, "STEP1_ID", chromosome)

    # _run_bcftools_step([
    #     "bash", "-c",
    #     f"""
    #     set -euo pipefail
    #     micromamba run -n vep java -jar /opt/snpEff/SnpSift.jar annotate \
    #         -id -tabix "{default_dbsnp_file}" "{norm_vcf}" |
    #     bcftools annotate --threads {threads} \
    #         --set-id +'%CHROM\\_%POS\\_%REF\\_%ALT' |
    #     bcftools view --threads {threads} \
    #         -Oz -o "{id_vcf}" \
    #         --write-index=tbi
    #     """
    # ], log_file, "STEP1_ID", chromosome)

    ##     bcftools norm --threads {threads} -Ou -m-any -d exact --fasta-ref "{genome_fasta_file}" |

    counts["ID"] = get_vcf_variant_count(str(id_vcf))
    log_variant_stats(log, "ID", counts["ID"], counts["NORM"], counts["ORIGINAL"], chromosome)

    # =======================================================
    # STEP 2: EAF ANNOTATION
    # =======================================================
    _run_bcftools_step([
        "bcftools", "annotate",
        "--threads", str(threads),
        "--annotations", str(external_eaf_file),
        "--pair-logic", "exact",
        "--columns", "CHROM,POS,REF,ALT,INFO/AFR,INFO/EAS,INFO/EUR,INFO/SAS",
        "-Oz", "-o", str(af_vcf),
        "--write-index=tbi",
        str(id_vcf)
    ], log_file, "STEP2_AF", chromosome)

    counts["EAF"] = get_vcf_variant_count(str(af_vcf))
    log_variant_stats(log, "EAF", counts["EAF"], counts["ID"], counts["ORIGINAL"], chromosome)

    # =======================================================
    # STEP 3: CSQ
    # =======================================================
    _run_bcftools_step([
        "bcftools", "csq",
        "--fasta-ref", str(genome_fasta_file),
        "--gff-annot", str(gff_file),
        "--threads", str(threads),
        "--unify-chr-names", "-,chr,-",
        "-Oz", "-o", str(csq_vcf),
        "--write-index=tbi",
        str(af_vcf)
    ], log_file, "STEP3_CSQ", chromosome)

    counts["CSQ"] = get_vcf_variant_count(str(csq_vcf))
    log_variant_stats(log, "CSQ", counts["CSQ"], counts["EAF"], counts["ORIGINAL"], chromosome)

    # =======================================================
    # STEP 4: LIFTOVER
    # =======================================================
    _run_bcftools_step([
        "bash", "-c",
        f"""
        set -euo pipefail

        bcftools +liftover "{csq_vcf}" --no-version -Ou -- \
            --src-fasta-ref "{genome_fasta_file}" \
            --fasta-ref "{target_genome_fasta_file}" \
            --chain "{chain_file}" \
            --reject "{reject_vcf}" \
            --reject-type z |

        bcftools view --threads {threads} -e "INFO/SWAP==1 || INFO/SWAP==-1" |

        bcftools norm --threads {threads} -Ou -m-any -d exact --fasta-ref "{target_genome_fasta_file}" |

        bcftools sort -m {threads * 512}M --temp-dir "{sort_tmp_dir}" -Oz -o "{target_vcf}" --write-index=tbi

        # SAFE INDEX REJECT
        if [ -f "{reject_vcf}" ] && [ -s "{reject_vcf}" ]; then
            tabix -f -p vcf "{reject_vcf}" || true
        fi
        """
    ], log_file, "STEP4_LIFTOVER", chromosome)

    # =======================================================
    # POST-LIFTOVER QC
    # =======================================================
    counts["LIFTED"] = get_vcf_variant_count(str(target_vcf))
    counts["REJECTED"] = (
        get_vcf_variant_count(str(reject_vcf))
        if Path(reject_vcf).exists()
        else 0
    )

    log_variant_stats(log, "LIFTED", counts["LIFTED"], counts["CSQ"], counts["ORIGINAL"], chromosome)

    if counts["REJECTED"]:
        log(f"[chr{chromosome}] 📉 REJECTED (not lifted): {counts['REJECTED']:,}")

    # 🚨 LIFTOVER FAILURE WARNING
    if counts["CSQ"] and counts["REJECTED"] is not None:
        failure_rate = counts["REJECTED"] / counts["CSQ"]
        log(f"[chr{chromosome}] 📊 Liftover failure rate: {failure_rate*100:.2f}%")

        if failure_rate > 0.30:
            log(f"[chr{chromosome}] 🚨 WARNING: >30% variants failed liftover")
        elif failure_rate > 0.10:
            log(f"[chr{chromosome}] ⚠️ Notice: >10% variants failed liftover")
    shutil.rmtree(sort_tmp_dir, ignore_errors=True)
    return str(target_vcf)

# ============================================================
# 3. CONCAT FUNCTION
# ============================================================
def concat_vcfs_by_build(
    output_dir: str,
    gwas_outputname: str,
    grch_version: str,
):
    from pathlib import Path
    import subprocess, re, os
    from concurrent.futures import ThreadPoolExecutor

    outdir = Path(output_dir)
    target_build = "GRCh38" if grch_version == "GRCh37" else "GRCh37"

    # -------------------------------------------------------
    # AUTO CPU DISTRIBUTION
    # -------------------------------------------------------
    total_cpus = os.cpu_count() or 4
    max_workers = min(3, total_cpus)
    threads_per_job = max(1, total_cpus // max_workers)

    # print("\t\t\t\t Started CONCATENATING CHROMOSOMAL VCF FILES")
    # print(f"\t\t\t\t🧠 CPU detected: {total_cpus}")
    # print(f"\t\t\t\t⚙️ Parallel jobs: {max_workers}")
    # print(f"\t\t\t\t🔧 Threads per job: {threads_per_job}")

    # -------------------------------------------------------
    # CHR SORT
    # -------------------------------------------------------
    def chr_sort_key(path):
        m = re.search(r"chr(\w+)", path.name)
        if not m:
            return (1000, path.name)

        v = m.group(1).upper()
        chr_map = {"X": 23, "Y": 24, "MT": 25, "M": 25}

        if v.isdigit():
            n = int(v)
            if 1 <= n <= 22:
                return (n, "")

        if v in chr_map:
            return (chr_map[v], "")

        return (1000, v)

    # -------------------------------------------------------
    # 🔥 VALIDATION (AUTO FIX + SKIP BAD FILES)
    # -------------------------------------------------------
    def validate_and_fix_vcf(vcf):
        vcf = Path(vcf)

        if not vcf.exists():
            print(f"❌ Missing VCF: {vcf}")
            return None

        # skip tiny files (header-only / broken)
        if vcf.stat().st_size < 100:
            print(f"⚠️ Skipping empty VCF: {vcf}")
            return None

        tbi = Path(str(vcf) + ".tbi")

        # try to auto-index if missing
        if not tbi.exists():
            print(f"⚠️ Missing index → attempting to index: {vcf}")
            try:
                subprocess.run(
                    ["tabix", "-f", "-p", "vcf", str(vcf)],
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
            except Exception:
                print(f"❌ Failed to index → skipping: {vcf}")
                return None

        # final check
        if not tbi.exists():
            print(f"❌ Still no index → skipping: {vcf}")
            return None

        return vcf

    # -------------------------------------------------------
    # CSQ FILES
    # -------------------------------------------------------
    csq_files = list(outdir.glob(f"{gwas_outputname}_chr*_GRCh*_CSQ.vcf.gz"))

    csq_37 = sorted(
        [validate_and_fix_vcf(f) for f in csq_files if "GRCh37" in f.name],
        key=chr_sort_key,
    )
    csq_37 = [f for f in csq_37 if f]

    csq_38 = sorted(
        [validate_and_fix_vcf(f) for f in csq_files if "GRCh38" in f.name],
        key=chr_sort_key,
    )
    csq_38 = [f for f in csq_38 if f]

    # -------------------------------------------------------
    # NOT LIFTED FILES
    # -------------------------------------------------------
    notlifted_raw = list(
        outdir.glob(f"{gwas_outputname}_chr*_{target_build}_CSQ_notlifted.vcf.gz")
    )

    notlifted_files = sorted(
        [validate_and_fix_vcf(f) for f in notlifted_raw],
        key=chr_sort_key,
    )
    notlifted_files = [f for f in notlifted_files if f]

    #print(f"🔍 Valid notlifted files: {len(notlifted_files)}")

    # -------------------------------------------------------
    # RAW FILES
    # -------------------------------------------------------
    raw_files = [
        f for f in outdir.glob(f"{gwas_outputname}_chr*_{grch_version}.vcf.gz")
        if f.name.endswith(f"_{grch_version}.vcf.gz")
        and all(tag not in f.name for tag in ["_NORM", "_ID", "_EAF", "_CSQ"])
    ]

    raw_files = sorted(
        [validate_and_fix_vcf(f) for f in raw_files],
        key=chr_sort_key,
    )
    raw_files = [f for f in raw_files if f]

    # -------------------------------------------------------
    # MERGE FUNCTION (SAFE)
    # -------------------------------------------------------
    def merge(tag, vcf_list, out):
        if not vcf_list:
            print(f"⚠️ [{tag}] No valid files — skipping")
            return None

        #print(f"\t\t\t\t🔧 [{tag}] Merging {len(vcf_list)} files")

        try:
            subprocess.run([
                "bcftools", "concat",
                "-a",
                "--output-type", "z",
                "--output", str(out),
                "--threads", str(threads_per_job),
                "--write-index=tbi",
                *map(str, vcf_list)
            ], check=True)

            #print(f"✅ [{tag}] → {out}")
            return str(out)

        except subprocess.CalledProcessError:
            print(f"❌ [{tag}] concat failed — skipping")
            return None

    # -------------------------------------------------------
    # PARALLEL EXECUTION
    # -------------------------------------------------------
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {}

        if csq_37:
            futures["grch37"] = executor.submit(
                merge, "GRCh37", csq_37,
                outdir / f"{gwas_outputname}_GRCh37_merged.vcf.gz"
            )

        if csq_38:
            futures["grch38"] = executor.submit(
                merge, "GRCh38", csq_38,
                outdir / f"{gwas_outputname}_GRCh38_merged.vcf.gz"
            )

        if raw_files:
            futures["gwas2vcf"] = executor.submit(
                merge, f"gwas2vcf_{grch_version}", raw_files,
                outdir / f"{gwas_outputname}_gwas2vcf_{grch_version}_merged.vcf.gz"
            )

        if notlifted_files:
            futures["notlifted"] = executor.submit(
                merge, f"notlifted_{target_build}", notlifted_files,
                outdir / f"{gwas_outputname}_notlifted_{target_build}_merged.vcf.gz"
            )

    # -------------------------------------------------------
    # WAIT
    # -------------------------------------------------------
    results = {
        "grch37": futures["grch37"].result() if "grch37" in futures else None,
        "grch38": futures["grch38"].result() if "grch38" in futures else None,
        "gwas2vcf": futures["gwas2vcf"].result() if "gwas2vcf" in futures else None,
        "notlifted": futures["notlifted"].result() if "notlifted" in futures else None,
    }

    # -------------------------------------------------------
    # CLEANUP
    # -------------------------------------------------------
    os.system(f"rm {outdir}/*_chr*_GRCh*_ORIGINAL.vcf.gz*")
    os.system(f"rm {outdir}/*_chr*_GRCh*_NORM.vcf.gz*")
    os.system(f"rm {outdir}/*_chr*_GRCh*_ID.vcf.gz*")
    os.system(f"rm {outdir}/*_chr*_GRCh*_EAF.vcf.gz*")
    os.system(f"rm {outdir}/*_chr*_GRCh*_CSQ.vcf.gz*")
    os.system(f"rm {outdir}/*_chr*_GRCh*.vcf.gz*")

    return results