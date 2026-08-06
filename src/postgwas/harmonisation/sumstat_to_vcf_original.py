import os
import re
import io
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

# ============================================================
# 1. GWAS → VCF via bcftools +munge
# ============================================================

def run_bcftools_munge(
    output_dir: str,
    gwas_outputname: str,
    chr: str,
    genome_fasta_file: str,
    grch_version: str,
) -> str:
    """
    Run bcftools +munge → bcftools norm → bcftools sort to convert GWAS summary stats to VCF.

    Parameters
    ----------
    output_dir : str
        Directory where input/output files are located.
    gwas_outputname : str
        Base GWAS dataset name.
    chr : str
        Chromosome number or label.
    genome_fasta_file : str
        Path to genome FASTA reference file.
    grch_version : str
        Genome build (GRCh37/GRCh38).

    Returns
    -------
    str
        Path to the final VCF.gz file.
    """
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    input_file = outdir / f"{gwas_outputname}_chr{chr}_vcf_input.tsv"
    column_file = outdir / f"{gwas_outputname}_chr{chr}_column_mapping.tsv"
    output_vcf = outdir / f"{gwas_outputname}_chr{chr}_{grch_version}.vcf.gz"

    cmd = (
        f"bcftools +munge --no-version "
        f"--columns-file {column_file} "
        f"--fasta-ref {genome_fasta_file} "
        f"--sample-name {gwas_outputname} "
        f"{input_file} | "
        f"bcftools norm -m-any -d exact | "
        f"bcftools sort -Oz -o {output_vcf} --write-index=tbi"
    )

    print(f"🚀 Running BCFtools munge/norm/sort for chr{chr}:\n{cmd}\n")
    try:
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
        print(f"✅ Completed: {output_vcf}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running bcftools munge/norm/sort:\n{e}")
        return "Error"

    return str(output_vcf)


# ============================================================
# 2. SAFE bcftools annotation + liftover (Option C)
#    - No pipes
#    - Each step separate
#    - Fully logged
# ============================================================

def _run_bcftools_step(args, log_file: Path, step_name: str):
    """
    Helper to run a single bcftools command, capturing stdout/stderr to a log file.
    Raises RuntimeError if return code != 0.
    """
    result = subprocess.run(
        args,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    with log_file.open("a") as lf:
        lf.write(f"\n[STEP {step_name}] CMD: {' '.join(map(str, args))}\n")
        if result.stdout:
            lf.write("[STDOUT]\n")
            lf.write(result.stdout)
            lf.write("\n")
        if result.stderr:
            lf.write("[STDERR]\n")
            lf.write(result.stderr)
            lf.write("\n")
        lf.write(f"[EXIT CODE] {result.returncode}\n")

    if result.returncode != 0:
        raise RuntimeError(
            f"bcftools step '{step_name}' failed with exit code {result.returncode}. "
            f"See log: {log_file}"
        )

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
    """
    bcftools annotation + liftover, preserving the original PostGWAS logic.
    Steps 4–6 use a high-speed piped implementation.
    """

    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------------
    # LOG FILE
    # -------------------------------------------------------
    log_dir = outdir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_bcftools_annot.log"

    def log(msg: str):
        with open(log_file, "a") as lf:
            lf.write(msg + "\n")

    log(f"▶ Starting bcftools annotation for chromosome {chromosome}")

    # -------------------------------------------------------
    # ORIGINAL PATHS
    # -------------------------------------------------------
    input_vcf   = outdir / f"{gwas_outputname}_chr{chromosome}_{grch_version}.vcf.gz"
    output1_vcf = outdir / f"{gwas_outputname}_chr{chromosome}_{grch_version}_ID.vcf.gz"
    output2_vcf = outdir / f"{gwas_outputname}_chr{chromosome}_{grch_version}_EAF.vcf.gz"
    output3_vcf = outdir / f"{gwas_outputname}_chr{chromosome}_{grch_version}_CSQ.vcf.gz"

    if grch_version == "GRCh37":
        target_build = "GRCh38"
    else:
        target_build = "GRCh37"

    target_vcf = outdir / f"{gwas_outputname}_chr{chromosome}_{target_build}.vcf.gz"
    reject_vcf = outdir / f"{gwas_outputname}_chr{chromosome}_{target_build}_notlifted.vcf.gz"
    
    
    tmp_norm = f"{output1_vcf}.pre_norm.vcf.gz"
    tmp_annot = f"{output1_vcf}.annot.vcf.gz"

    # -------------------------------------------------------
    # Pre-normalisation
    # -------------------------------------------------------
    cmd0 = [
        "bash", "-c",
        f"bcftools norm -m-any -d exact \"{input_vcf}\" "
        f"| bgzip -c > \"{tmp_norm}\" "
        f"&& tabix -f -p vcf \"{tmp_norm}\" "
    ]
    _run_bcftools_step(cmd0, log_file, "annotate_norm_split")

    # ❌ OLD (buggy):
    # os.system(f"mv {tmp_norm} > {input_vcf}")
    # ✅ NEW (correct move, also move index)
    os.replace(tmp_norm, input_vcf)
    tbi_tmp = f"{tmp_norm}.tbi"
    tbi_final = f"{input_vcf}.tbi"
    if os.path.exists(tbi_tmp):
        os.replace(tbi_tmp, tbi_final)

    # -------------------------------------------------------
    # Step 1: dbSNP annotation
    # -------------------------------------------------------
    log("🧬 Step 1: dbSNP annotation + multi-allelic split...")

    cmd1 = [
        "bash", "-c",
        f"bcftools annotate "
        f"--threads {threads} "
        f"--annotations \"{default_dbsnp_file}\" "
        f"--columns CHROM,POS,REF,ALT,ID "
        f"\"{input_vcf}\" "
        f"| bcftools norm -m-any -d exact "
        f"| bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' "
        f"| bgzip -c > \"{output1_vcf}\" "
        f"&& tabix -f -p vcf \"{output1_vcf}\""
    ]

    _run_bcftools_step(cmd1, log_file, "annotate_norm_split")

    # -------------------------------------------------------
    # Step 2: AF annotation
    # -------------------------------------------------------
    log("🌍 Step 2/4: AF annotation...")
    cmd2 = [
        "bcftools", "annotate",
        "--threads", str(threads),
        "--annotations", external_eaf_file,
        "--columns", "CHROM,POS,REF,ALT,INFO/AFR,INFO/EAS,INFO/EUR,INFO/SAS",
        "-Oz",
        "-o", str(output2_vcf),
        "--write-index=tbi",
        str(output1_vcf),
    ]
    _run_bcftools_step(cmd2, log_file, "af_annotate")

    # -------------------------------------------------------
    # Step 3: bcftools csq (OVERWRITE input_vcf)
    # -------------------------------------------------------
    log("🧫 Step 3/4: Functional consequence annotation (bcftools csq)...")
    cmd3 = [
        "bcftools", "csq",
        "--fasta-ref", genome_fasta_file,
        "--gff-annot", gff_file,
        "--threads", str(threads),
        "--unify-chr-names", "-,chr,-",
        "-Oz",
        "-o", str(input_vcf),
        "--write-index=tbi",
        str(output2_vcf),
    ]
    _run_bcftools_step(cmd3, log_file, "csq")

    # Optional CSQ copy
    try:
        import shutil
        shutil.copy2(input_vcf, output3_vcf)
        shutil.copy2(str(input_vcf) + ".tbi", str(output3_vcf) + ".tbi")
    except Exception:
        pass

    # -------------------------------------------------------
    # Steps 4–6 combined PIPELINE (FASTEST IMPLEMENTATION)
    # -------------------------------------------------------
    log(f"🧭 Step 4–6: liftover → filter SWAP → sort (pipelined)")

    # ✅ NEW: correct +liftover syntax, no --reject-type z
    pipeline_cmd = f"""
        bcftools +liftover \"{input_vcf}\" --no-version -Ou -- \
            --src-fasta-ref \"{genome_fasta_file}\" \
            --fasta-ref \"{target_genome_fasta_file}\" \
            --chain \"{chain_file}\" \
            --reject \"{reject_vcf}\" \
            --reject-type z \
        | bcftools view -e 'INFO/SWAP==1 || INFO/SWAP==-1' \
        | bcftools norm -m-any -d exact \
        | bcftools sort -Oz \
            -o \"{target_vcf}\" \
            --write-index=tbi
    """

    _run_bcftools_step(["bash", "-c", pipeline_cmd], log_file, "liftover_view_sort")

    log(f"✅ DONE → {target_vcf}")
    return f"Completed bcftools annotation for chr{chromosome}"

# ============================================================
# 3. Concatenate per-chromosome VCFs by build
# ============================================================

def concat_vcfs_by_build(
    output_dir: str,
    gwas_outputname: str,
    mode: str = "concurrent",
):
    """
    Concatenate chromosome-wise VCFs for:
      • GRCh37
      • GRCh38
      • GRCh37_notlifted
      • GRCh38_notlifted

    Writes a single log:
      logs/{gwas_outputname}_concat_vcfs_by_build.log

    Returns
    -------
    dict
        {
          "grch37": Path or None,
          "grch38": Path or None,
          "notlifted_37": Path or None,
          "notlifted_38": Path or None
        }
    """
    outdir = Path(output_dir).resolve()
    log_dir = outdir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{gwas_outputname}_concat_vcfs_by_build.log"

    log_buffer = io.StringIO()

    def log_write(*args):
        msg = " ".join(str(a) for a in args)
        log_buffer.write(msg + "\n")

    def chr_sort_key(path: Path):
        m = re.search(r"chr(\d+|X|Y|MT)", path.stem)
        if not m:
            return 999
        v = m.group(1)
        if v == "X":
            return 23
        if v == "Y":
            return 24
        if v == "MT":
            return 25
        return int(v)

    all_vcfs = list(outdir.glob(f"{gwas_outputname}_chr*_*GRCh*.vcf.gz"))
    if not all_vcfs:
        log_write("❌ No chromosome VCF files found in output directory.")
        with open(log_file, "w") as f:
            f.write(log_buffer.getvalue())
        raise FileNotFoundError("No chromosome VCFs found")

    valid_vcfs = [
        f for f in all_vcfs
        if "_notlifted" not in f.name
        and (f.name.endswith("_GRCh37.vcf.gz") or f.name.endswith("_GRCh38.vcf.gz"))
    ]
    notlifted_vcfs = [f for f in all_vcfs if "_notlifted" in f.name]

    grch37_vcfs = sorted([p for p in valid_vcfs if "GRCh37" in p.name], key=chr_sort_key)
    grch38_vcfs = sorted([p for p in valid_vcfs if "GRCh38" in p.name], key=chr_sort_key)
    notlifted_37 = sorted([p for p in notlifted_vcfs if "GRCh37" in p.name], key=chr_sort_key)
    notlifted_38 = sorted([p for p in notlifted_vcfs if "GRCh38" in p.name], key=chr_sort_key)

    merged_paths = {
        "grch37":       outdir / f"{gwas_outputname}_GRCh37_merged.vcf.gz",
        "grch38":       outdir / f"{gwas_outputname}_GRCh38_merged.vcf.gz",
        "notlifted_37": outdir / f"{gwas_outputname}_GRCh37_notlifted_merged.vcf.gz",
        "notlifted_38": outdir / f"{gwas_outputname}_GRCh38_notlifted_merged.vcf.gz",
    }

    def merge_vcfs(tag: str, vcf_list, out_path: Path):
        if not vcf_list:
            log_write(f"⚠️ No VCF files found for {tag}. Skipping.")
            return None

        log_write(f"🔧 [{tag}] Running bcftools concat on {len(vcf_list)} files…")

        result = subprocess.run(
            [
                "bcftools", "concat",
                "--output-type", "z",
                "--output", str(out_path),
                "--threads", "4",
                "--write-index=tbi",
                *[str(f) for f in vcf_list],
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        if result.stdout.strip():
            log_write(f"[bcftools stdout:{tag}] {result.stdout.strip()}")
        if result.stderr.strip():
            log_write(f"[bcftools stderr:{tag}] {result.stderr.strip()}")

        if result.returncode != 0:
            log_write(f"❌ bcftools concat failed for {tag} (exit {result.returncode})")
            return None

        log_write(f"✅ [{tag}] Merge completed → {out_path}")

        # Clean chromosome-level files
        for f in vcf_list:
            try:
                f.unlink(missing_ok=True)
                tbi = f.with_suffix(f.suffix + ".tbi")
                if tbi.exists():
                    tbi.unlink()
            except Exception:
                pass

        log_write(f"🧹 [{tag}] Cleaned per-chromosome VCF + index files.")
        return str(out_path)

    tasks = []
    results = []

    if mode == "concurrent":
        with ThreadPoolExecutor(max_workers=4) as ex:
            if grch37_vcfs:
                tasks.append(ex.submit(merge_vcfs, "GRCh37", grch37_vcfs, merged_paths["grch37"]))
            if grch38_vcfs:
                tasks.append(ex.submit(merge_vcfs, "GRCh38", grch38_vcfs, merged_paths["grch38"]))
            if notlifted_37:
                tasks.append(ex.submit(merge_vcfs, "GRCh37_notlifted", notlifted_37, merged_paths["notlifted_37"]))
            if notlifted_38:
                tasks.append(ex.submit(merge_vcfs, "GRCh38_notlifted", notlifted_38, merged_paths["notlifted_38"]))
            results = [t.result() for t in tasks]
    else:
        if grch37_vcfs:
            results.append(merge_vcfs("GRCh37", grch37_vcfs, merged_paths["grch37"]))
        if grch38_vcfs:
            results.append(merge_vcfs("GRCh38", grch38_vcfs, merged_paths["grch38"]))
        if notlifted_37:
            results.append(merge_vcfs("GRCh37_notlifted", notlifted_37, merged_paths["notlifted_37"]))
        if notlifted_38:
            results.append(merge_vcfs("GRCh38_notlifted", notlifted_38, merged_paths["notlifted_38"]))

    # Extra cleanup for any leftover chr*.vcf.gz or indexes
    for f in outdir.glob(f"{gwas_outputname}_chr*.vcf.gz*"):
        try:
            f.unlink()
        except Exception:
            pass

    log_write("🎯 All merges completed successfully.")

    with open(log_file, "w") as f:
        f.write(log_buffer.getvalue())

    return merged_paths
