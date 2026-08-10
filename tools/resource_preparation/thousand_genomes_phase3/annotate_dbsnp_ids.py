#!/usr/bin/env python3
"""Add exact dbSNP IDs or CHROM_POS_REF_ALT fallback IDs to 1000G VCFs."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


# Edit resource and compute settings here.
RESOURCE_DIRECTORY = Path(
    "/Users/JJOHN41/Documents/software_resources/resourses/postgwas_v2/"
    "1000_genomes_phase3/population_vcfs"
)
DBSNP_DIRECTORY = Path(
    "/Users/JJOHN41/Documents/software_resources/resourses/postgwas_v2/"
    "dbsnp_build157_grch37/derived/per_chromosome_biallelic"
)
REFERENCE_DIRECTORY = Path(
    "/Users/JJOHN41/Documents/software_resources/resourses/postgwas/"
    "gwas2vcf/GRCh37/fasta_files"
)
POPULATIONS = ("AFR", "AMR", "EAS", "EUR", "SAS")
CHROMOSOMES = [str(value) for value in range(1, 23)] + ["X"]
INPUT_TEMPLATE = "{population}.GRCh37.phase3.autosomes_chrX.vcf.gz"
DBSNP_TEMPLATE = "dbsnp_build157.GRCh37.chr{chromosome}.biallelic.vcf.gz"
OUTPUT_TEMPLATE = "{population}.GRCh37.phase3.autosomes_chrX.with_dbsnp_rsid.vcf.gz"
BCFTOOLS_EXECUTABLE = "/Users/JJOHN41/miniconda3/bin/bcftools"
WORKERS = 4


def reference_for(chromosome: str) -> Path:
    for candidate in (
        REFERENCE_DIRECTORY / f"GRCh37_chr{chromosome}.fa",
        REFERENCE_DIRECTORY / f"GRCh37_chr{chromosome}",
    ):
        if candidate.is_file():
            for index in (Path(f"{candidate}.fai"), candidate.with_suffix(".fai")):
                if index.is_file():
                    return candidate
    raise FileNotFoundError(f"GRCh37 FASTA and index not found for chr{chromosome}")


def run_chromosome(population: str, chromosome: str, temporary: Path) -> Path:
    source = RESOURCE_DIRECTORY / INPUT_TEMPLATE.format(population=population)
    annotation = DBSNP_DIRECTORY / DBSNP_TEMPLATE.format(chromosome=chromosome)
    reference = reference_for(chromosome)
    output = temporary / f"{population}.chr{chromosome}.vcf.gz"
    view = subprocess.Popen(
        [BCFTOOLS_EXECUTABLE, "view", "--regions", chromosome, "--output-type", "u", "--no-version", str(source)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    if view.stdout is None:
        raise RuntimeError(f"view pipe unavailable: {population} chr{chromosome}")
    norm = subprocess.Popen(
        [BCFTOOLS_EXECUTABLE, "norm", "--fasta-ref", str(reference), "--multiallelics", "-any", "-d", "exact", "--output-type", "u", "--no-version", "-"],
        stdin=view.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    view.stdout.close()
    if norm.stdout is None:
        view.kill(); norm.kill()
        raise RuntimeError(f"norm pipe unavailable: {population} chr{chromosome}")
    fallback = subprocess.Popen(
        [BCFTOOLS_EXECUTABLE, "annotate", "--set-id", "%CHROM_%POS_%REF_%FIRST_ALT", "--output-type", "u", "--no-version", "-"],
        stdin=norm.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    norm.stdout.close()
    if fallback.stdout is None:
        view.kill(); norm.kill(); fallback.kill()
        raise RuntimeError(f"fallback-ID pipe unavailable: {population} chr{chromosome}")
    annotate = subprocess.Popen(
        [BCFTOOLS_EXECUTABLE, "annotate", "--annotations", str(annotation), "--columns", "ID", "--pair-logic", "exact", "--output-type", "z", "--output", str(output), "--no-version", "-"],
        stdin=fallback.stdout, stderr=subprocess.PIPE,
    )
    fallback.stdout.close()
    annotate_return = annotate.wait()
    fallback_return = fallback.wait()
    norm_return = norm.wait()
    view_return = view.wait()
    if any(code != 0 for code in (view_return, norm_return, fallback_return, annotate_return)):
        errors = []
        for process, label in ((view, "view"), (norm, "norm"), (fallback, "set-id"), (annotate, "annotate")):
            if process.stderr is not None:
                errors.append(f"{label}: {process.stderr.read().decode(errors='replace')}")
        raise RuntimeError(f"failed {population} chr{chromosome}: {'; '.join(errors)}")
    subprocess.run([BCFTOOLS_EXECUTABLE, "index", "--force", "--tbi", str(output)], check=True)
    return output


def main() -> int:
    if shutil.which(BCFTOOLS_EXECUTABLE) is None:
        raise RuntimeError("bcftools is unavailable")
    for population in POPULATIONS:
        source = RESOURCE_DIRECTORY / INPUT_TEMPLATE.format(population=population)
        if not source.is_file():
            raise FileNotFoundError(source)
        for chromosome in CHROMOSOMES:
            annotation = DBSNP_DIRECTORY / DBSNP_TEMPLATE.format(chromosome=chromosome)
            if not annotation.is_file() or not Path(f"{annotation}.tbi").is_file():
                raise FileNotFoundError(annotation)
    output_dir = RESOURCE_DIRECTORY / "with_dbsnp_rsid"
    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="rsid_", dir=output_dir) as temp_name:
        temporary = Path(temp_name)
        with ThreadPoolExecutor(max_workers=WORKERS) as executor:
            futures = {
                executor.submit(run_chromosome, population, chromosome, temporary): (population, chromosome)
                for population in POPULATIONS for chromosome in CHROMOSOMES
            }
            for future in as_completed(futures):
                population, chromosome = futures[future]
                print(f"Completed {population} chr{chromosome}", flush=True)
                future.result()
        for population in POPULATIONS:
            parts = [temporary / f"{population}.chr{chromosome}.vcf.gz" for chromosome in CHROMOSOMES]
            output = output_dir / OUTPUT_TEMPLATE.format(population=population)
            partial = Path(f"{output}.partial")
            subprocess.run([BCFTOOLS_EXECUTABLE, "concat", "--output-type", "z", "--output", str(partial), *map(str, parts)], check=True)
            subprocess.run([BCFTOOLS_EXECUTABLE, "index", "--force", "--tbi", str(partial)], check=True)
            partial.replace(output)
            Path(f"{partial}.tbi").replace(Path(f"{output}.tbi"))
            print(f"Published {output}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr, flush=True)
        raise SystemExit(1) from error
