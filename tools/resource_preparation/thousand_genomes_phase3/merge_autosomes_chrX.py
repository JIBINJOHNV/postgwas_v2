#!/usr/bin/env python3
"""Create one valid population VCF containing autosomes and chromosome X."""

from __future__ import annotations

from pathlib import Path
import os
import shutil
import subprocess
import sys


# Edit resource and output settings here.
RESOURCE_DIRECTORY = Path(
    "/Users/JJOHN41/Documents/software_resources/resourses/postgwas_v2/"
    "1000_genomes_phase3/population_vcfs"
)
POPULATIONS = ("AFR", "AMR", "EAS", "EUR", "SAS")
BUILD_LABEL = "GRCh37.phase3"
OUTPUT_SUFFIX = "autosomes_chrX"
BCFTOOLS_EXECUTABLE = "bcftools"


def sample_names(path: Path) -> list[str]:
    result = subprocess.run(
        [BCFTOOLS_EXECUTABLE, "query", "-l", str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.splitlines()


def merge_population(population: str) -> Path:
    autosomes = RESOURCE_DIRECTORY / f"{population}.{BUILD_LABEL}.autosomes.vcf.gz"
    chromosome_x = RESOURCE_DIRECTORY / f"{population}.{BUILD_LABEL}.chrX.vcf.gz"
    output = RESOURCE_DIRECTORY / f"{population}.{BUILD_LABEL}.{OUTPUT_SUFFIX}.vcf.gz"
    partial = Path(f"{output}.partial")
    for path in (autosomes, chromosome_x):
        if not path.is_file() or not Path(f"{path}.csi").is_file():
            raise FileNotFoundError(f"missing indexed input: {path}")
    if sample_names(autosomes) != sample_names(chromosome_x):
        raise RuntimeError(f"sample order differs between autosomes and chrX: {population}")

    subprocess.run(
        [BCFTOOLS_EXECUTABLE, "concat", "-Oz", "-o", str(partial), str(autosomes), str(chromosome_x)],
        check=True,
    )
    subprocess.run([BCFTOOLS_EXECUTABLE, "index", "-f", str(partial)], check=True)
    partial.replace(output)
    Path(f"{partial}.csi").replace(Path(f"{output}.csi"))
    return output


def main() -> int:
    if shutil.which(BCFTOOLS_EXECUTABLE) is None:
        raise RuntimeError("bcftools is not available on PATH")
    for population in POPULATIONS:
        output = merge_population(population)
        print(f"Created {output}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr, flush=True)
        raise SystemExit(1) from error
