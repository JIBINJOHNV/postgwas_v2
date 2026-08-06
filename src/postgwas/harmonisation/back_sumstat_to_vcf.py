from __future__ import annotations

import io
import os
import re
import shutil
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Optional


# ============================================================
# Helpers
# ============================================================

def _assert_exists(path_like, label: str = "file") -> Path:
    p = Path(path_like)
    if not p.exists():
        raise FileNotFoundError(f"Missing required {label}: {p}")
    return p


def _assert_nonempty(path_like, label: str = "file") -> Path:
    p = Path(path_like)
    if not p.exists():
        raise FileNotFoundError(f"Expected output {label} not found: {p}")
    if p.stat().st_size == 0:
        raise RuntimeError(f"Expected output {label} is empty: {p}")
    return p


def _safe_replace_gz_and_index(tmp_gz: Path, final_gz: Path) -> None:
    """
    Atomically replace final .vcf.gz with tmp .vcf.gz and also move .tbi if present.
    """
    _assert_nonempty(tmp_gz, "temporary gz output")
    os.replace(tmp_gz, final_gz)

    tmp_tbi = Path(str(tmp_gz) + ".tbi")
    final_tbi = Path(str(final_gz) + ".tbi")
    if tmp_tbi.exists():
        os.replace(tmp_tbi, final_tbi)


def _count_variants(vcf_gz: Path) -> int:
    """
    Count non-header records in a VCF/BCF using bcftools view -H.
    """
    result = subprocess.run(
        ["bash", "-c", f'bcftools view -H "{vcf_gz}" | wc -l'],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        return -1
    try:
        return int(result.stdout.strip())
    except Exception:
        return -1


def _run_bcftools_step(args, log_file: Path, step_name: str) -> None:
    """
    Run a single command and capture stdout/stderr to a log file.
    Raises RuntimeError on failure.
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
            if not result.stdout.endswith("\n"):
                lf.write("\n")
        if result.stderr:
            lf.write("[STDERR]\n")
            lf.write(result.stderr)
            if not result.stderr.endswith("\n"):
                lf.write("\n")
        lf.write(f"[EXIT CODE] {result.returncode}\n")

    if result.returncode != 0:
        raise RuntimeError(
            f"bcftools step '{step_name}' failed with exit code {result.returncode}. "
            f"See log: {log_file}"
        )


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

    Returns
    -------
    str
        Path to final VCF.gz
    """
    outdir = Path(output_dir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    log_dir = outdir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{gwas_outputname}_chr{chr}_bcftools_munge.log"

    input_file = outdir / f"{gwas_outputname}_chr{chr}_vcf_input.tsv"
    column_file = outdir / f"{gwas_outputname}_chr{chr}_column_mapping.tsv"
    output_vcf = outdir / f"{gwas_outputname}_chr{chr}_{grch_version}.vcf.gz"

    _assert_exists(input_file, "munge input TSV")
    _assert_exists(column_file, "munge column mapping TSV")
    _assert_exists(genome_fasta_file, "reference FASTA")

    cmd = [
        "bash", "-c",
        (
            f'bcftools +munge --no-version '
            f'--columns-file "{column_file}" '
            f'--fasta-ref "{genome_fasta_file}" '
            f'--sample-name "{gwas_outputname}" '
            f'"{input_file}" '
            f'| bcftools norm -m-any -d exact '
            f'| bcftools sort -Oz -o "{output_vcf}" --write-index=tbi'
        )
    ]

    with log_file.open("a") as lf:
        lf.write(f"▶ Starting bcftools munge for chr{chr}\n")

    _run_bcftools_step(cmd, log_file, f"munge_chr{chr}")

    _assert_nonempty(output_vcf, "munge output VCF.gz")
    with log_file.open("a") as lf:
        lf.write(f"✅ Completed: {output_vcf}\n")

    return str(output_vcf)


# ============================================================
# 2. SAFE bcftools annotation + liftover
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
    """
    bcftools annotation + csq + liftover.
    Preserves your original design:
      - pre-normalize
      - dbSNP annotate
      - append unique allele-aware ID
      - EAF annotation
      - CSQ becomes final same-name source VCF via safe temp replace
      - liftover with reject file in compressed VCF format
    """
    outdir = Path(output_dir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    log_dir = outdir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_bcftools_annot.log"

    def log(msg: str) -> None:
        with log_file.open("a") as lf:
            lf.write(msg + "\n")

    # -------------------------------------------------------
    # Paths
    # -------------------------------------------------------
    input_vcf = outdir / f"{gwas_outputname}_chr{chromosome}_{grch_version}.vcf.gz"
    output1_vcf = outdir / f"{gwas_outputname}_chr{chromosome}_{grch_version}_ID.vcf.gz"
    output2_vcf = outdir / f"{gwas_outputname}_chr{chromosome}_{grch_version}_EAF.vcf.gz"
    output3_vcf = outdir / f"{gwas_outputname}_chr{chromosome}_{grch_version}_CSQ.vcf.gz"

    if grch_version == "GRCh37":
        target_build = "GRCh38"
    else:
        target_build = "GRCh37"

    target_vcf = outdir / f"{gwas_outputname}_chr{chromosome}_{target_build}.vcf.gz"
    reject_vcf = outdir / f"{gwas_outputname}_chr{chromosome}_{target_build}_notlifted.vcf.gz"

    tmp_norm = Path(str(output1_vcf) + ".pre_norm.vcf.gz")
    csq_tmp_vcf = Path(str(input_vcf) + ".csq.tmp.vcf.gz")

    # -------------------------------------------------------
    # Input validation
    # -------------------------------------------------------
    _assert_exists(input_vcf, "input VCF.gz")
    _assert_exists(str(input_vcf) + ".tbi", "input VCF index")
    _assert_exists(default_dbsnp_file, "dbSNP annotation VCF/BCF")
    _assert_exists(external_eaf_file, "external EAF annotation file")
    _assert_exists(genome_fasta_file, "source genome FASTA")
    _assert_exists(target_genome_fasta_file, "target genome FASTA")
    _assert_exists(gff_file, "GFF annotation file")
    _assert_exists(chain_file, "liftover chain file")

    log(f"▶ Starting bcftools annotation for chromosome {chromosome}")

    # -------------------------------------------------------
    # Step 0: Pre-normalisation
    # -------------------------------------------------------
    log("🧹 Step 0: Pre-normalisation...")
    cmd0 = [
        "bash", "-c",
        (
            f'bcftools norm -m-any -d exact "{input_vcf}" '
            f'| bgzip -c > "{tmp_norm}" '
            f'&& tabix -f -p vcf "{tmp_norm}"'
        )
    ]
    _run_bcftools_step(cmd0, log_file, "pre_norm_split")

    _safe_replace_gz_and_index(tmp_norm, input_vcf)
    log(f"✅ Pre-normalised input replaced in place: {input_vcf}")

    # -------------------------------------------------------
    # Step 1: dbSNP annotation + unique appended ID
    # -------------------------------------------------------
    log("🧬 Step 1: dbSNP annotation + multi-allelic split + appended unique ID...")
    cmd1 = [
        "bash", "-c",
        (
            f'bcftools annotate '
            f'--threads {threads} '
            f'--annotations "{default_dbsnp_file}" '
            f'--columns CHROM,POS,REF,ALT,ID '
            f'"{input_vcf}" '
            f'| bcftools norm -m-any -d exact '
            f"| bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' "
            f'| bgzip -c > "{output1_vcf}" '
            f'&& tabix -f -p vcf "{output1_vcf}"'
        )
    ]
    _run_bcftools_step(cmd1, log_file, "annotate_norm_split")
    _assert_nonempty(output1_vcf, "ID-annotated VCF.gz")

    # -------------------------------------------------------
    # Step 2: AF annotation
    # -------------------------------------------------------
    log("🌍 Step 2: AF annotation...")
    cmd2 = [
        "bcftools", "annotate",
        "--threads", str(threads),
        "--annotations", str(external_eaf_file),
        "--columns", "CHROM,POS,REF,ALT,INFO/AFR,INFO/EAS,INFO/EUR,INFO/SAS",
        "-Oz",
        "-o", str(output2_vcf),
        "--write-index=tbi",
        str(output1_vcf),
    ]
    _run_bcftools_step(cmd2, log_file, "af_annotate")
    _assert_nonempty(output2_vcf, "EAF-annotated VCF.gz")

    # -------------------------------------------------------
    # Step 3: CSQ annotation
    # Safe temp → atomic replace to keep final name = input_vcf
    # -------------------------------------------------------
    log("🧫 Step 3: Functional consequence annotation (bcftools csq)...")
    cmd3 = [
        "bcftools", "csq",
        "--fasta-ref", str(genome_fasta_file),
        "--gff-annot", str(gff_file),
        "--threads", str(threads),
        "--unify-chr-names", "-,chr,-",
        "-Oz",
        "-o", str(csq_tmp_vcf),
        "--write-index=tbi",
        str(output2_vcf),
    ]
    _run_bcftools_step(cmd3, log_file, "csq")
    _safe_replace_gz_and_index(csq_tmp_vcf, input_vcf)
    log(f"✅ CSQ output replaced in place: {input_vcf}")

    # Optional CSQ copy
    try:
        shutil.copy2(input_vcf, output3_vcf)
        shutil.copy2(str(input_vcf) + ".tbi", str(output3_vcf) + ".tbi")
        log(f"✅ CSQ copy written: {output3_vcf}")
    except Exception as e:
        log(f"⚠️ Optional CSQ copy failed: {e}")

    # -------------------------------------------------------
    # Step 4–6: liftover → filter SWAP → norm → sort
    # Note: --reject-type z is valid and preserved
    # -------------------------------------------------------
    log("🧭 Step 4–6: liftover → filter SWAP → sort...")
    pipeline_cmd = f"""
        bcftools +liftover "{input_vcf}" --no-version -Ou -- \\
            --src-fasta-ref "{genome_fasta_file}" \\
            --fasta-ref "{target_genome_fasta_file}" \\
            --chain "{chain_file}" \\
            --reject "{reject_vcf}" \\
            --reject-type z \\
        | bcftools view -e 'INFO/SWAP==1 || INFO/SWAP==-1' \\
        | bcftools norm -m-any -d exact \\
        | bcftools sort -Oz -o "{target_vcf}" --write-index=tbi
    """
    _run_bcftools_step(["bash", "-c", pipeline_cmd], log_file, "liftover_view_sort")
    _assert_nonempty(target_vcf, "lifted target VCF.gz")

    # -------------------------------------------------------
    # QC summaries
    # -------------------------------------------------------
    lifted_n = _count_variants(target_vcf)
    rejected_n = _count_variants(reject_vcf) if reject_vcf.exists() else 0
    source_n = _count_variants(input_vcf)

    log(f"[QC] source_variants={source_n}")
    log(f"[QC] lifted_variants={lifted_n}")
    log(f"[QC] rejected_variants={rejected_n}")
    if source_n > 0 and lifted_n >= 0:
        pct = (lifted_n / source_n) * 100
        log(f"[QC] lifted_percent={pct:.2f}")

    log(f"✅ DONE → {target_vcf}")
    return f"Completed bcftools annotation for chr{chromosome}"


# ============================================================
# 3. Concatenate per-chromosome VCFs by build
# ============================================================

def concat_vcfs_by_build(
    output_dir: str,
    gwas_outputname: str,
    mode: str = "concurrent",
    threads: int = 4,
) -> Dict[str, Optional[Path]]:
    """
    Concatenate chromosome-wise VCFs for:
      • GRCh37
      • GRCh38
      • GRCh37_notlifted
      • GRCh38_notlifted

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

    merged_paths: Dict[str, Optional[Path]] = {
        "grch37": outdir / f"{gwas_outputname}_GRCh37_merged.vcf.gz",
        "grch38": outdir / f"{gwas_outputname}_GRCh38_merged.vcf.gz",
        "notlifted_37": outdir / f"{gwas_outputname}_GRCh37_notlifted_merged.vcf.gz",
        "notlifted_38": outdir / f"{gwas_outputname}_GRCh38_notlifted_merged.vcf.gz",
    }

    def merge_vcfs(tag: str, vcf_list: List[Path], out_path: Path) -> Optional[str]:
        if not vcf_list:
            log_write(f"⚠️ No VCF files found for {tag}. Skipping.")
            return None

        log_write(f"🔧 [{tag}] Running bcftools concat on {len(vcf_list)} files...")

        result = subprocess.run(
            [
                "bcftools", "concat",
                "--output-type", "z",
                "--output", str(out_path),
                "--threads", str(threads),
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

        if not out_path.exists() or out_path.stat().st_size == 0:
            log_write(f"❌ [{tag}] output missing or empty: {out_path}")
            return None

        log_write(f"✅ [{tag}] Merge completed → {out_path}")

        # Clean chromosome-level files for this group only
        for f in vcf_list:
            try:
                f.unlink(missing_ok=True)
                tbi = Path(str(f) + ".tbi")
                tbi.unlink(missing_ok=True)
            except Exception as e:
                log_write(f"⚠️ [{tag}] cleanup failed for {f}: {e}")

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

    # Optional extra cleanup for leftover per-chromosome files only
    for f in outdir.glob(f"{gwas_outputname}_chr*.vcf.gz"):
        try:
            f.unlink()
        except Exception:
            pass
        tbi = Path(str(f) + ".tbi")
        try:
            tbi.unlink()
        except Exception:
            pass

    log_write("🎯 All merges completed.")
    log_write(f"results={results}")

    with open(log_file, "w") as f:
        f.write(log_buffer.getvalue())

    # Return None for any merged outputs that were expected but not created
    for key, path in list(merged_paths.items()):
        if path is not None and not path.exists():
            merged_paths[key] = None

    return merged_paths