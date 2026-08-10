#!/usr/bin/env python3
"""Create PostGWAS FinnGen allele-frequency tab files from the combined AF."""

from __future__ import annotations

from datetime import datetime, timezone
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import BinaryIO

POLARS_THREADS = 2
os.environ.setdefault("POLARS_MAX_THREADS", str(POLARS_THREADS))
import polars as pl


# Edit resource and processing settings here.
INPUT_FILE = Path(
    "/Users/JJOHN41/Documents/software_resources/resourses/postgwas_v2/"
    "finngen_r13_annotations/finngen_R13_annotated_variants_v0.gz"
)
OUTPUT_DIRECTORY = INPUT_FILE.parent / "frequency" / "tab_files"
OUTPUT_FILENAME_TEMPLATE = "GRCh38_finngen_freq_chr{chrom}.tsv.gz"
CHROMOSOMES = [str(value) for value in range(1, 23)] + ["X"]
CHROMOSOME_RENAMES = {"23": "X"}
BATCH_ROWS = 50_000
NULL_VALUES = ["NA", "."]
BGZIP_EXECUTABLE = "/Users/JJOHN41/miniconda3/bin/bgzip"
TABIX_EXECUTABLE = "/Users/JJOHN41/miniconda3/bin/tabix"
BGZIP_THREADS = 1
PARTIAL_MARKER = ".partial"
OUTPUT_HEADER = "CHROM\tPOS\tREF\tALT\tEUR\n"


class FrequencyError(RuntimeError):
    pass


def timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def executable(path: str) -> str:
    found = shutil.which(path)
    if found is None:
        raise FrequencyError(f"required executable is unavailable: {path}")
    return found


class ChromosomeWriter:
    def __init__(self, output_directory: Path, chromosome: str, bgzip: str) -> None:
        final = output_directory / OUTPUT_FILENAME_TEMPLATE.format(chrom=chromosome)
        partial = final.with_name(f"{final.name.removesuffix('.gz')}{PARTIAL_MARKER}.gz")
        partial.unlink(missing_ok=True)
        Path(f"{partial}.tbi").unlink(missing_ok=True)
        self.final = final
        self.partial = partial
        self.handle: BinaryIO = partial.open("wb")
        self.process = subprocess.Popen(
            [bgzip, "--threads", str(BGZIP_THREADS), "--stdout"],
            stdin=subprocess.PIPE,
            stdout=self.handle,
        )
        if self.process.stdin is None:
            self.process.kill()
            self.handle.close()
            raise FrequencyError("bgzip stdin pipe is unavailable")
        self.process.stdin.write(OUTPUT_HEADER.encode())

    def write(self, frame: pl.DataFrame) -> None:
        if self.process.stdin is None or self.process.stdin.closed:
            raise FrequencyError("bgzip stdin pipe closed unexpectedly")
        self.process.stdin.write(frame.write_csv(separator="\t", include_header=False).encode())

    def finish(self, tabix: str, bgzip: str, chromosome: str) -> None:
        if self.process.stdin is not None:
            self.process.stdin.close()
        if self.process.wait() != 0:
            self.handle.close()
            raise FrequencyError(f"bgzip failed for chromosome {chromosome}")
        self.handle.close()
        subprocess.run(
            [tabix, "--sequence", "1", "--begin", "2", "--end", "2", "--skip-lines", "1", "--force", str(self.partial)],
            check=True,
        )
        subprocess.run([bgzip, "--test", str(self.partial)], check=True)
        observed = subprocess.run([tabix, "--list-chroms", str(self.partial)], check=True, text=True, capture_output=True).stdout.splitlines()
        if observed != [chromosome]:
            raise FrequencyError(f"unexpected indexed chromosome(s): {observed}")
        self.partial.replace(self.final)
        Path(f"{self.partial}.tbi").replace(Path(f"{self.final}.tbi"))

    def abort(self) -> None:
        if self.process.stdin is not None and not self.process.stdin.closed:
            self.process.stdin.close()
        if self.process.poll() is None:
            self.process.terminate()
        self.process.wait()
        self.handle.close()
        self.partial.unlink(missing_ok=True)
        Path(f"{self.partial}.tbi").unlink(missing_ok=True)


def transform(batch: pl.DataFrame) -> pl.DataFrame:
    frame = batch.with_columns(
        pl.col("chr").replace(CHROMOSOME_RENAMES).cast(pl.String).alias("CHROM"),
        pl.col("pos").cast(pl.Int64, strict=True).alias("POS"),
        pl.col("ref").cast(pl.String).alias("REF"),
        pl.col("alt").cast(pl.String).alias("ALT"),
        pl.col("AF").cast(pl.Float64, strict=True).alias("EUR"),
    ).select("CHROM", "POS", "REF", "ALT", "EUR")
    invalid = frame.filter(
        ~pl.col("CHROM").is_in(CHROMOSOMES)
        | (pl.col("POS") < 1)
        | pl.col("REF").is_null()
        | pl.col("ALT").is_null()
        | (pl.col("EUR") < 0)
        | (pl.col("EUR") > 1)
    )
    if invalid.height:
        raise FrequencyError(f"batch contains {invalid.height} invalid frequency rows")
    return frame


def main() -> int:
    if not INPUT_FILE.is_file():
        raise FrequencyError(f"input file is missing: {INPUT_FILE}")
    bgzip = executable(BGZIP_EXECUTABLE)
    tabix = executable(TABIX_EXECUTABLE)
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    lazy = pl.scan_csv(
        INPUT_FILE,
        separator="\t",
        schema_overrides={"chr": pl.String, "pos": pl.Int64, "ref": pl.String, "alt": pl.String, "AF": pl.Float64},
        null_values=NULL_VALUES,
        infer_schema_length=10_000,
        low_memory=True,
    ).select("chr", "pos", "ref", "alt", "AF")
    rank = {chromosome: index for index, chromosome in enumerate(CHROMOSOMES)}
    writer: ChromosomeWriter | None = None
    previous_chromosome: str | None = None
    previous_position: int | None = None
    counts = {chromosome: 0 for chromosome in CHROMOSOMES}
    try:
        for raw_batch in lazy.collect_batches(chunk_size=BATCH_ROWS, maintain_order=True, engine="streaming"):
            batch = transform(raw_batch)
            for value in batch["CHROM"].unique(maintain_order=True):
                chromosome = str(value)
                part = batch.filter(pl.col("CHROM") == chromosome)
                if previous_chromosome is not None and rank[chromosome] < rank[previous_chromosome]:
                    raise FrequencyError("input chromosome order is not monotonic")
                if chromosome == previous_chromosome and previous_position is not None and int(part["POS"][0]) < previous_position:
                    raise FrequencyError(f"input positions are not sorted on chromosome {chromosome}")
                if chromosome != previous_chromosome:
                    if writer is not None:
                        writer.finish(tabix, bgzip, previous_chromosome)
                    writer = ChromosomeWriter(OUTPUT_DIRECTORY, chromosome, bgzip)
                    previous_position = None
                if not part["POS"].is_sorted():
                    raise FrequencyError(f"input positions are not sorted within chromosome {chromosome}")
                writer.write(part)
                counts[chromosome] += part.height
                previous_chromosome = chromosome
                previous_position = int(part["POS"][-1])
        if writer is not None and previous_chromosome is not None:
            writer.finish(tabix, bgzip, previous_chromosome)
    except BaseException:
        if writer is not None:
            writer.abort()
        raise
    if not any(counts.values()):
        raise FrequencyError("input contains no variants")
    print(f"[{timestamp()}] Created FinnGen frequency files in {OUTPUT_DIRECTORY}", flush=True)
    print(f"[{timestamp()}] Rows by chromosome: {counts}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FrequencyError, OSError, pl.exceptions.PolarsError, subprocess.CalledProcessError) as error:
        print(f"[{timestamp()}] ERROR: {error}", file=sys.stderr, flush=True)
        raise SystemExit(1) from error
