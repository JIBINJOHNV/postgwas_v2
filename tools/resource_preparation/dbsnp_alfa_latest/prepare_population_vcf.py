#!/usr/bin/env python3
"""Convert NCBI ALFA frequencies into a compact PostGWAS VCF."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone

import pysam


# Edit resource and output settings here.
INPUT_VCF = Path(
    "/Users/JJOHN41/Documents/software_resources/resourses/postgwas_v2/"
    "dbsnp_alfa_latest/freq.vcf.gz"
)
OUTPUT_DIRECTORY = INPUT_VCF.parent / "derived"
OUTPUT_FILENAME_TEMPLATE = "GRCh38_ALFA_freq_chr{chromosome}.vcf.gz"
BGZIP_EXECUTABLE = "bgzip"
TABIX_EXECUTABLE = "tabix"
OUTPUT_QUAL = "100"
OUTPUT_FILTER = "PASS"
AF_DECIMAL_PLACES = 10

# ALFA BioSample IDs mapped to the requested PostGWAS super-populations.
# AFR and AMR are pooled from their component ALFA ancestry groups using AC/AN.
POPULATION_SAMPLE_IDS = {
    "EAS": ("SAMN10492697",),  # East Asian
    "AMR": ("SAMN10492699", "SAMN10492700"),  # Latin American 1 and 2
    "AFR": ("SAMN10492696", "SAMN10492698", "SAMN10492703"),
    "EUR": ("SAMN10492695",),  # European
    "SAS": ("SAMN10492702",),  # South Asian
}
OUTPUT_POPULATIONS = ("EAS", "AMR", "AFR", "EUR", "SAS")

# GRCh38 primary chromosomes in the ALFA VCF.
ACCESSION_TO_CHROMOSOME = {
    **{f"NC_{number:06d}.{version}": str(number) for number, version in (
        (1, 11), (2, 12), (3, 12), (4, 12), (5, 10), (6, 11), (7, 14),
        (8, 11), (9, 12), (10, 11), (11, 11), (12, 13), (13, 11), (14, 9),
        (15, 10), (16, 11), (17, 12), (18, 11), (19, 9), (20, 10),
        (21, 9), (22, 11),
    )},
    "NC_000023.11": "X",
    "NC_000024.10": "Y",
    "NC_012920.1": "MT",
}


def timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def format_frequency(value: float) -> str:
    if value == 0:
        return "0"
    text = f"{value:.{AF_DECIMAL_PLACES}f}".rstrip("0").rstrip(".")
    return text or "0"


def parse_count(value: str) -> int:
    if value in {"", "."}:
        return 0
    return int(value)


def write_header(handle) -> None:
    header = [
        "##fileformat=VCFv4.2",
        "##source=NCBI_ALFA_latest_release_population_frequency",
        "##reference=GRCh38",
        '##INFO=<ID=EAS,Number=A,Type=Float,Description="ALFA East Asian alternate allele frequency">',
        '##INFO=<ID=AMR,Number=A,Type=Float,Description="ALFA pooled Latin American alternate allele frequency">',
        '##INFO=<ID=AFR,Number=A,Type=Float,Description="ALFA pooled African alternate allele frequency">',
        '##INFO=<ID=EUR,Number=A,Type=Float,Description="ALFA European alternate allele frequency">',
        '##INFO=<ID=SAS,Number=A,Type=Float,Description="ALFA South Asian alternate allele frequency">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    handle.write(("\n".join(header) + "\n").encode())


def convert() -> dict[str, int]:
    """Stream ALFA records into chromosome-specific temporary BGZF files."""
    if not INPUT_VCF.is_file():
        raise FileNotFoundError(f"input VCF does not exist: {INPUT_VCF}")
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    rows = defaultdict(int)
    skipped = defaultdict(int)
    with tempfile.TemporaryDirectory(prefix="alfa_population_", dir=OUTPUT_DIRECTORY) as temp_dir:
        writers = {}
        try:
            with pysam.BGZFile(str(INPUT_VCF), "r") as source:
                sample_ids = None
                format_indexes = None
                for raw_line in source:
                    line = raw_line.decode().rstrip("\n")
                    if line.startswith("#CHROM"):
                        columns = line.split("\t")
                        sample_ids = columns[9:]
                        sample_indexes = {name: index for index, name in enumerate(sample_ids)}
                        format_indexes = True
                        required = {sample for ids in POPULATION_SAMPLE_IDS.values() for sample in ids}
                        missing = sorted(required.difference(sample_ids))
                        if missing:
                            raise RuntimeError(f"ALFA VCF is missing population columns: {missing}")
                        continue
                    if line.startswith("#"):
                        continue
                    if sample_ids is None or format_indexes is None:
                        raise RuntimeError("VCF data appeared before the #CHROM header")
                    fields = line.split("\t")
                    record_format_indexes = {
                        name: index for index, name in enumerate(fields[8].split(":"))
                    }
                    if "AN" not in record_format_indexes or "AC" not in record_format_indexes:
                        skipped["missing_an_or_ac"] += 1
                        continue
                    chromosome = ACCESSION_TO_CHROMOSOME.get(fields[0])
                    if chromosome is None:
                        skipped["non_primary_contig"] += 1
                        continue
                    alt_alleles = fields[4].split(",")
                    sample_values = fields[9:]
                    population_counts = {}
                    for population, population_samples in POPULATION_SAMPLE_IDS.items():
                        total_an = 0
                        alt_counts = [0] * len(alt_alleles)
                        for sample_id in population_samples:
                            sample = sample_values[sample_indexes[sample_id]].split(":")
                            total_an += parse_count(sample[record_format_indexes["AN"]])
                            ac = sample[record_format_indexes["AC"]].split(",")
                            for index, value in enumerate(ac[:len(alt_alleles)]):
                                alt_counts[index] += parse_count(value)
                        population_counts[population] = (
                            total_an,
                            alt_counts,
                        )
                    if chromosome not in writers:
                        path = Path(temp_dir) / f"chr{chromosome}.vcf.gz"
                        writers[chromosome] = pysam.BGZFile(str(path), "w")
                        write_header(writers[chromosome])
                    for alt_index, alt in enumerate(alt_alleles):
                        info_values = []
                        for population in OUTPUT_POPULATIONS:
                            total_an, alt_counts = population_counts[population]
                            frequency = alt_counts[alt_index] / total_an if total_an else 0.0
                            info_values.append(f"{population}={format_frequency(frequency)}")
                        output = [chromosome, fields[1], fields[2] or ".", fields[3], alt,
                                  OUTPUT_QUAL, OUTPUT_FILTER, ";".join(info_values)]
                        writers[chromosome].write(("\t".join(output) + "\n").encode())
                        rows[chromosome] += 1
        finally:
            for writer in writers.values():
                writer.close()

        for chromosome in [*(str(i) for i in range(1, 23)), "X", "Y", "MT"]:
            path = Path(temp_dir) / f"chr{chromosome}.vcf.gz"
            if not path.exists():
                continue
            output_vcf = OUTPUT_DIRECTORY / OUTPUT_FILENAME_TEMPLATE.format(chromosome=chromosome)
            partial_output = Path(f"{output_vcf}.partial")
            path.replace(partial_output)
            subprocess.run([TABIX_EXECUTABLE, "-f", "-p", "vcf", str(partial_output)], check=True)
            partial_output.replace(output_vcf)
            Path(f"{partial_output}.tbi").replace(Path(f"{output_vcf}.tbi"))
            print(f"[{timestamp()}] Published {output_vcf}", flush=True)
    return {**rows, "skipped_non_primary_contig": skipped["non_primary_contig"]}


def main() -> int:
    if shutil.which(BGZIP_EXECUTABLE) is None or shutil.which(TABIX_EXECUTABLE) is None:
        raise RuntimeError("bgzip and tabix must be available on PATH")
    counts = convert()
    print(
        f"[{timestamp()}] ALFA population VCFs complete in {OUTPUT_DIRECTORY}",
        flush=True,
    )
    print(f"[{timestamp()}] Row counts: {dict(counts)}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"[{timestamp()}] ERROR: {error}", file=sys.stderr, flush=True)
        raise SystemExit(1) from error
