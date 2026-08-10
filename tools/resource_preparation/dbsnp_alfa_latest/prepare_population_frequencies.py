#!/usr/bin/env python3
"""Create a canonical-GRCh37 VCF of ALFA superpopulation ALT frequencies."""

from __future__ import annotations

from collections.abc import Iterable
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


# Edit resource, population, and compute settings here.
INPUT_VCF = Path(
    "/Users/JJOHN41/Documents/software_resources/resourses/postgwas_v2/"
    "dbsnp_alfa_latest/freq.vcf.gz"
)
OUTPUT_DIRECTORY = Path(
    "/Users/JJOHN41/Documents/software_resources/resourses/postgwas_v2/"
    "dbsnp_alfa_latest/derived"
)
OUTPUT_FILENAME = "dbsnp_alfa_grch37_superpopulation_af.vcf.gz"
LOG_FILENAME = "population_frequency_preparation.jsonl"
BCFTOOLS_EXECUTABLE = "/Users/JJOHN41/miniconda3/envs/postgwas/bin/bcftools"
BGZIP_EXECUTABLE = "/Users/JJOHN41/miniconda3/envs/postgwas/bin/bgzip"
TABIX_EXECUTABLE = "/Users/JJOHN41/miniconda3/envs/postgwas/bin/tabix"
COMPRESSION_THREADS = 4
PARTIAL_SUFFIX = ".partial"
AF_SIGNIFICANT_DIGITS = 10

# Canonical GRCh37 primary-assembly RefSeq accessions in output chromosome order.
GRCH37_CONTIGS = {
    "1": "NC_000001.10",
    "2": "NC_000002.11",
    "3": "NC_000003.11",
    "4": "NC_000004.11",
    "5": "NC_000005.9",
    "6": "NC_000006.11",
    "7": "NC_000007.13",
    "8": "NC_000008.10",
    "9": "NC_000009.11",
    "10": "NC_000010.10",
    "11": "NC_000011.9",
    "12": "NC_000012.11",
    "13": "NC_000013.10",
    "14": "NC_000014.8",
    "15": "NC_000015.9",
    "16": "NC_000016.9",
    "17": "NC_000017.10",
    "18": "NC_000018.9",
    "19": "NC_000019.9",
    "20": "NC_000020.10",
    "21": "NC_000021.8",
    "22": "NC_000022.10",
    "X": "NC_000023.10",
}

# NCBI BioSample columns in the ALFA VCF. AMR pools the two ALFA Latin-American
# strata by allele counts; the other four outputs use the direct ALFA aggregates.
POPULATION_SAMPLE_GROUPS = {
    "EAS": ("SAMN10492697",),
    "AMR": ("SAMN10492699", "SAMN10492700"),
    "AFR": ("SAMN10492703",),
    "EUR": ("SAMN10492695",),
    "SAS": ("SAMN10492702",),
}


def timestamp() -> str:
    """Return an ISO-8601 UTC timestamp."""
    return datetime.now(timezone.utc).isoformat()


def log_event(handle, event: str, **details: object) -> None:
    """Write one durable JSON Lines progress event."""
    record = {"timestamp_utc": timestamp(), "event": event, **details}
    handle.write(json.dumps(record, sort_keys=True) + "\n")
    handle.flush()


def population_af(
    sample_values: Iterable[str], format_keys: list[str], alt_index: int
) -> str:
    """Pool AC and AN across samples and return one ALT-allele frequency."""
    try:
        an_index = format_keys.index("AN")
        ac_index = format_keys.index("AC")
    except ValueError as error:
        raise RuntimeError("ALFA record FORMAT must contain AN and AC") from error

    total_an = 0
    total_ac = 0
    for value in sample_values:
        fields = value.split(":")
        if max(an_index, ac_index) >= len(fields):
            raise RuntimeError(f"malformed ALFA sample value: {value}")
        an_text = fields[an_index]
        ac_values = fields[ac_index].split(",")
        if an_text in {"", "."} or alt_index >= len(ac_values):
            continue
        ac_text = ac_values[alt_index]
        if ac_text in {"", "."}:
            continue
        an = int(an_text)
        ac = int(ac_text)
        if an < 0 or ac < 0 or ac > an:
            raise RuntimeError(f"invalid ALFA allele counts AN={an}, AC={ac}")
        total_an += an
        total_ac += ac

    if total_an == 0:
        return "."
    return format(total_ac / total_an, f".{AF_SIGNIFICANT_DIGITS}g")


def convert_record(
    line: str,
    sample_indexes: dict[str, tuple[int, ...]],
    accession_to_chromosome: dict[str, str],
) -> list[str]:
    """Split one ALFA record into biallelic records with population INFO fields."""
    columns = line.rstrip("\n").split("\t")
    if len(columns) < 10:
        raise RuntimeError("malformed ALFA VCF record with fewer than 10 columns")
    chromosome = accession_to_chromosome.get(columns[0])
    if chromosome is None:
        raise RuntimeError(f"unexpected source contig selected: {columns[0]}")

    alts = columns[4].split(",")
    format_keys = columns[8].split(":")
    records = []
    for alt_index, alt in enumerate(alts):
        if alt in {"", "."}:
            raise RuntimeError(f"invalid ALT allele at {columns[0]}:{columns[1]}")
        info = []
        for population, indexes in sample_indexes.items():
            values = (columns[index] for index in indexes)
            info.append(f"{population}={population_af(values, format_keys, alt_index)}")
        output = [
            chromosome,
            columns[1],
            columns[2],
            columns[3],
            alt,
            columns[5],
            columns[6],
            ";".join(info),
        ]
        records.append("\t".join(output) + "\n")
    return records


def output_header() -> str:
    """Build the valid sites-only output VCF header."""
    lines = [
        "##fileformat=VCFv4.2",
        "##source=NCBI_dbSNP_ALFA_population_frequency_latest_release",
        "##reference=GRCh37.p13",
        "##ALFA_AMR_definition=pooled_LAC_and_LEN_allele_counts",
    ]
    lines.extend(f"##contig=<ID={chromosome}>" for chromosome in GRCH37_CONTIGS)
    for population in POPULATION_SAMPLE_GROUPS:
        lines.append(
            f'##INFO=<ID={population},Number=1,Type=Float,'
            f'Description="ALFA ALT allele frequency for {population}; missing when AN=0">'
        )
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    return "\n".join(lines) + "\n"


def executable(path: str) -> str:
    """Resolve one configured executable or fail with an actionable error."""
    resolved = shutil.which(path)
    if resolved is None:
        raise RuntimeError(f"required executable was not found: {path}")
    return resolved


def main() -> int:
    """Extract GRCh37 ALFA records, calculate AF, compress, index, and publish."""
    if not INPUT_VCF.is_file():
        raise RuntimeError(f"input VCF was not found: {INPUT_VCF}")
    input_index = Path(f"{INPUT_VCF}.tbi")
    if not input_index.is_file():
        raise RuntimeError(f"input Tabix index was not found: {input_index}")

    bcftools = executable(BCFTOOLS_EXECUTABLE)
    bgzip = executable(BGZIP_EXECUTABLE)
    tabix = executable(TABIX_EXECUTABLE)
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    destination = OUTPUT_DIRECTORY / OUTPUT_FILENAME
    partial = Path(f"{destination}{PARTIAL_SUFFIX}")
    partial_index = Path(f"{partial}.tbi")
    for path in (partial, partial_index):
        if path.exists():
            path.unlink()

    selected_samples = tuple(
        sample
        for samples in POPULATION_SAMPLE_GROUPS.values()
        for sample in samples
    )
    regions = ",".join(GRCH37_CONTIGS.values())
    log_path = OUTPUT_DIRECTORY / LOG_FILENAME
    with log_path.open("a", encoding="utf-8") as log_handle:
        log_event(
            log_handle,
            "run_started",
            input=str(INPUT_VCF),
            output=str(destination),
            populations=POPULATION_SAMPLE_GROUPS,
        )
        reader = subprocess.Popen(
            [
                bcftools,
                "view",
                "--regions",
                regions,
                "--samples",
                ",".join(selected_samples),
                "--output-type",
                "v",
                str(INPUT_VCF),
            ],
            stdout=subprocess.PIPE,
            text=True,
        )
        writer = subprocess.Popen(
            [bgzip, "--threads", str(COMPRESSION_THREADS), "--stdout"],
            stdin=subprocess.PIPE,
            stdout=partial.open("wb"),
            text=True,
        )
        if reader.stdout is None or writer.stdin is None:
            raise RuntimeError("failed to open the VCF processing streams")

        sample_indexes: dict[str, tuple[int, ...]] | None = None
        input_records = 0
        output_records = 0
        writer.stdin.write(output_header())
        try:
            for line in reader.stdout:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    names = line.rstrip("\n").split("\t")
                    missing = [sample for sample in selected_samples if sample not in names]
                    if missing:
                        raise RuntimeError(
                            "required ALFA samples are absent: " + ", ".join(missing)
                        )
                    sample_indexes = {
                        population: tuple(names.index(sample) for sample in samples)
                        for population, samples in POPULATION_SAMPLE_GROUPS.items()
                    }
                    continue
                if sample_indexes is None:
                    raise RuntimeError("source VCF column header was not found")
                input_records += 1
                for record in convert_record(
                    line, sample_indexes, {v: k for k, v in GRCH37_CONTIGS.items()}
                ):
                    writer.stdin.write(record)
                    output_records += 1
            reader.stdout.close()
            reader_code = reader.wait()
            writer.stdin.close()
            writer_code = writer.wait()
        except BaseException:
            reader.kill()
            writer.kill()
            reader.wait()
            writer.wait()
            raise
        if reader_code != 0:
            raise RuntimeError(f"bcftools view failed with exit code {reader_code}")
        if writer_code != 0:
            raise RuntimeError(f"bgzip failed with exit code {writer_code}")

        subprocess.run(
            [tabix, "--preset", "vcf", str(partial)],
            check=True,
        )
        subprocess.run([bcftools, "view", "--header-only", str(partial)], check=True)
        subprocess.run([tabix, "--list-chroms", str(partial)], check=True)
        os.replace(partial, destination)
        os.replace(partial_index, Path(f"{destination}.tbi"))
        log_event(
            log_handle,
            "run_completed",
            status="success",
            input_records=input_records,
            output_records=output_records,
        )

    print(
        f"[{timestamp()}] Completed ALFA population VCF: {destination}\n"
        f"[{timestamp()}] Source records: {input_records:,}; "
        f"biallelic output records: {output_records:,}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"[{timestamp()}] ERROR: {error}", file=sys.stderr, flush=True)
        raise SystemExit(1) from error
