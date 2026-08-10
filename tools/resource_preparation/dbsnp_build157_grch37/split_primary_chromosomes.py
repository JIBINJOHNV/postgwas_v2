#!/usr/bin/env python3
"""Create chromosome-wise biallelic dbSNP GRCh37 reference VCFs."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


# Edit resource and compute settings here.
SOURCE_VCF = Path(
    "/Users/JJOHN41/Documents/software_resources/resourses/postgwas_v2/"
    "dbsnp_build157_grch37/downloads/GCF_000001405.25.gz"
)
CHROMOSOME_MAP = Path(
    "/Users/JJOHN41/Documents/software_resources/resourses/postgwas_v2/"
    "dbsnp_build157_grch37/metadata/GRCh37p13_refseq_to_standard_chromosomes.tsv"
)
OUTPUT_DIRECTORY = Path(
    "/Users/JJOHN41/Documents/software_resources/resourses/postgwas_v2/"
    "dbsnp_build157_grch37/derived/per_chromosome_biallelic"
)
REFERENCE_DIRECTORY = Path(
    "/Users/JJOHN41/Documents/software_resources/resourses/postgwas/"
    "gwas2vcf/GRCh37/fasta_files"
)
BCFTOOLS_EXECUTABLE = "/Users/JJOHN41/miniconda3/bin/bcftools"
TABIX_EXECUTABLE = "/Users/JJOHN41/miniconda3/bin/tabix"
WORKERS = 4
BCFTOOLS_THREADS_PER_WORKER = 1
OUTPUT_TEMPLATE = "dbsnp_build157.GRCh37.chr{chromosome}.biallelic.vcf.gz"
CHROMOSOMES = [str(value) for value in range(1, 23)] + ["X"]


def read_mapping() -> list[tuple[str, str]]:
    rows = []
    for line in CHROMOSOME_MAP.read_text(encoding="utf-8").splitlines():
        accession, chromosome = line.split("\t")
        rows.append((accession, chromosome))
    return rows


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


def convert_one(accession: str, chromosome: str) -> Path:
    output = OUTPUT_DIRECTORY / OUTPUT_TEMPLATE.format(chromosome=chromosome)
    reference = reference_for(chromosome)
    partial = Path(f"{output}.partial")
    partial.unlink(missing_ok=True)
    Path(f"{partial}.tbi").unlink(missing_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{chromosome}.", dir=OUTPUT_DIRECTORY) as temp_dir:
        temporary = Path(temp_dir) / output.name
        annotate_map = CHROMOSOME_MAP
        view = subprocess.Popen(
            [BCFTOOLS_EXECUTABLE, "view", "--regions", accession, "--output-type", "u", "--no-version", str(SOURCE_VCF)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if view.stdout is None:
            raise RuntimeError(f"bcftools view pipe unavailable for chromosome {chromosome}")
        annotate = subprocess.Popen(
            [BCFTOOLS_EXECUTABLE, "annotate", "--rename-chrs", str(annotate_map), "--remove", "INFO", "--output-type", "u", "--no-version", "-"],
            stdin=view.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        view.stdout.close()
        if annotate.stdout is None:
            view.kill()
            annotate.kill()
            raise RuntimeError(f"bcftools annotate pipe unavailable for chromosome {chromosome}")
        norm = subprocess.Popen(
            [BCFTOOLS_EXECUTABLE, "norm", "--fasta-ref", str(reference), "--multiallelics", "-any", "-d", "exact", "--threads", str(BCFTOOLS_THREADS_PER_WORKER), "--output-type", "z", "--output", str(temporary), "--no-version", "-"],
            stdin=annotate.stdout,
            stderr=subprocess.PIPE,
        )
        annotate.stdout.close()
        norm_return = norm.wait()
        annotate_return = annotate.wait()
        view_return = view.wait()
        if view_return or annotate_return or norm_return:
            errors = []
            for process, label in ((view, "view"), (annotate, "annotate"), (norm, "norm")):
                if process.stderr is not None:
                    errors.append(f"{label}: {process.stderr.read().decode(errors='replace')}")
            raise RuntimeError(f"chromosome {chromosome} conversion failed; {'; '.join(errors)}")
        subprocess.run([BCFTOOLS_EXECUTABLE, "index", "--force", "--tbi", str(temporary)], check=True)
        partial.parent.mkdir(parents=True, exist_ok=True)
        temporary.replace(partial)
        Path(f"{temporary}.tbi").replace(Path(f"{partial}.tbi"))
    partial.replace(output)
    Path(f"{partial}.tbi").replace(Path(f"{output}.tbi"))
    return output


def main() -> int:
    if not SOURCE_VCF.is_file() or not CHROMOSOME_MAP.is_file():
        raise FileNotFoundError("source VCF or chromosome map is missing")
    if shutil.which(BCFTOOLS_EXECUTABLE) is None or shutil.which(TABIX_EXECUTABLE) is None:
        raise RuntimeError("bcftools and tabix are required")
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    mapping = [
        (accession, chromosome)
        for accession, chromosome in read_mapping()
        if chromosome in CHROMOSOMES
    ]
    with ThreadPoolExecutor(max_workers=WORKERS) as executor:
        futures = {executor.submit(convert_one, accession, chromosome): chromosome for accession, chromosome in mapping}
        for future in as_completed(futures):
            chromosome = futures[future]
            print(f"Created {future.result()}", flush=True)
            print(f"Validated chromosome {chromosome}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr, flush=True)
        raise SystemExit(1) from error
