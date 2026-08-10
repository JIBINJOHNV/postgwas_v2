#!/usr/bin/env python3
"""Create original, mean, and median FinnGen R13 INFO-score tracks."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import BinaryIO

# Polars reads this setting when imported.
POLARS_THREADS = 2
os.environ.setdefault("POLARS_MAX_THREADS", str(POLARS_THREADS))
import polars as pl


# Edit resource and processing settings here.
INPUT_FILE = Path(
    "/Users/JJOHN41/Documents/software_resources/resourses/postgwas_v2/"
    "finngen_r13_annotations/finngen_R13_annotated_variants_v0.gz"
)
OUTPUT_DIRECTORY = INPUT_FILE.parent / "info_scores"
OUTPUT_FILENAME_TEMPLATES = {
    "original": "GRCh38_finngenR13_Original_infoscore_chr{chrom}.tsv.gz",
    "mean": "GRCh38_finngenR13_Mean_infoscore_chr{chrom}.tsv.gz",
    "median": "GRCh38_finngenR13_Median_infoscore_chr{chrom}.tsv.gz",
}
CHROMOSOMES = [str(value) for value in range(1, 23)] + ["X"]
CHROMOSOME_RENAMES = {"23": "X"}
COORDINATE_COLUMNS = ["chr", "pos", "ref", "alt"]
ORIGINAL_INFO_COLUMN = "INFO"
BATCH_INFO_PREFIX = "INFO_"
OUTPUT_HEADER = "CHROM\tPOS\tREF\tALT\tINFO\n"
NULL_VALUES = ["NA", "."]
BATCH_ROWS = 50_000
BGZIP_THREADS = 1
BGZIP_EXECUTABLE = "/Users/JJOHN41/miniconda3/bin/bgzip"
TABIX_EXECUTABLE = "/Users/JJOHN41/miniconda3/bin/tabix"
PARTIAL_MARKER = ".partial"
LOG_FILENAME = "preparation.jsonl"


class FinnGenInfoError(RuntimeError):
    """Raised when an INFO-score track cannot be created safely."""


def timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def log_event(log_path: Path, event: str, **details: object) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    record = {"timestamp_utc": timestamp(), "event": event, **details}
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")


def read_header(source: Path) -> list[str]:
    with gzip.open(source, "rt", encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"), None)
    if not header:
        raise FinnGenInfoError(f"input has no tab-delimited header: {source}")
    if len(header) != len(set(header)):
        raise FinnGenInfoError("input header contains duplicate column names")
    return header


def select_columns(header: list[str]) -> tuple[list[str], list[str]]:
    required = COORDINATE_COLUMNS + [ORIGINAL_INFO_COLUMN]
    missing = [column for column in required if column not in header]
    if missing:
        raise FinnGenInfoError(f"input is missing required columns: {missing}")
    batch_info = [column for column in header if column.startswith(BATCH_INFO_PREFIX)]
    if not batch_info:
        raise FinnGenInfoError(
            f"input contains no columns beginning with {BATCH_INFO_PREFIX!r}"
        )
    return required + batch_info, batch_info


def executable(configured: str) -> str:
    resolved = shutil.which(configured)
    if resolved is None:
        raise FinnGenInfoError(f"required executable is unavailable: {configured}")
    return resolved


def run(command: list[str]) -> None:
    subprocess.run(command, check=True)


def output_paths(output_directory: Path, chrom: str) -> dict[str, Path]:
    return {
        statistic: output_directory / template.format(chrom=chrom)
        for statistic, template in OUTPUT_FILENAME_TEMPLATES.items()
    }


def validate_track(path: Path, chrom: str, bgzip: str, tabix: str) -> None:
    index = Path(f"{path}.tbi")
    if not path.is_file() or path.stat().st_size == 0 or not index.is_file():
        raise FinnGenInfoError(f"track or Tabix index is missing: {path}")
    run([bgzip, "--test", str(path)])
    observed = subprocess.run(
        [tabix, "--list-chroms", str(path)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    ).stdout.splitlines()
    if observed != [chrom]:
        raise FinnGenInfoError(
            f"indexed chromosomes differ for {path}: observed {observed}, "
            f"expected {[chrom]}"
        )


class ChromosomeWriters:
    """Stream three atomic BGZF tracks for one chromosome."""

    def __init__(self, output_directory: Path, chrom: str, bgzip: str) -> None:
        self.chrom = chrom
        self.final_paths = output_paths(output_directory, chrom)
        self.partial_paths: dict[str, Path] = {}
        self.processes: dict[str, subprocess.Popen[bytes]] = {}
        self.output_handles: dict[str, BinaryIO] = {}
        for statistic, final in self.final_paths.items():
            partial = final.with_name(
                f"{final.name.removesuffix('.gz')}{PARTIAL_MARKER}.gz"
            )
            partial.unlink(missing_ok=True)
            Path(f"{partial}.tbi").unlink(missing_ok=True)
            output_handle = partial.open("wb")
            process = subprocess.Popen(
                [bgzip, "--threads", str(BGZIP_THREADS), "--stdout"],
                stdin=subprocess.PIPE,
                stdout=output_handle,
            )
            if process.stdin is None:
                process.kill()
                output_handle.close()
                raise FinnGenInfoError("bgzip stdin pipe is unavailable")
            process.stdin.write(OUTPUT_HEADER.encode("utf-8"))
            self.partial_paths[statistic] = partial
            self.processes[statistic] = process
            self.output_handles[statistic] = output_handle

    def write(self, frames: dict[str, pl.DataFrame]) -> None:
        for statistic, frame in frames.items():
            process = self.processes[statistic]
            if process.stdin is None:
                raise FinnGenInfoError("bgzip stdin pipe closed unexpectedly")
            process.stdin.write(
                frame.write_csv(separator="\t", include_header=False).encode("utf-8")
            )

    def abort(self) -> None:
        for process in self.processes.values():
            if process.stdin is not None and not process.stdin.closed:
                process.stdin.close()
            if process.poll() is None:
                process.terminate()
            process.wait()
        for output_handle in self.output_handles.values():
            output_handle.close()
        for partial in self.partial_paths.values():
            partial.unlink(missing_ok=True)
            Path(f"{partial}.tbi").unlink(missing_ok=True)

    def finish(self, tabix: str, bgzip: str) -> None:
        for statistic, process in self.processes.items():
            if process.stdin is not None:
                process.stdin.close()
            return_code = process.wait()
            self.output_handles[statistic].close()
            if return_code != 0:
                self.abort()
                raise FinnGenInfoError(
                    f"bgzip failed for chromosome {self.chrom}/{statistic}"
                )
        try:
            for statistic, partial in self.partial_paths.items():
                run([
                    tabix,
                    "--sequence",
                    "1",
                    "--begin",
                    "2",
                    "--end",
                    "2",
                    "--skip-lines",
                    "1",
                    "--force",
                    str(partial),
                ])
                validate_track(partial, self.chrom, bgzip, tabix)
            for statistic, partial in self.partial_paths.items():
                final = self.final_paths[statistic]
                os.replace(partial, final)
                os.replace(Path(f"{partial}.tbi"), Path(f"{final}.tbi"))
        except BaseException:
            self.abort()
            raise


def transform_batch(batch: pl.DataFrame, batch_info: list[str]) -> pl.DataFrame:
    transformed = batch.with_columns(
        pl.col("chr")
        .replace(CHROMOSOME_RENAMES)
        .cast(pl.String)
        .alias("CHROM"),
        pl.col("pos").cast(pl.Int64, strict=True).alias("POS"),
        pl.col("ref").cast(pl.String).alias("REF"),
        pl.col("alt").cast(pl.String).alias("ALT"),
        pl.col(ORIGINAL_INFO_COLUMN).cast(pl.Float64, strict=True).alias(
            "INFO_ORIGINAL"
        ),
        pl.mean_horizontal(batch_info).alias("INFO_MEAN"),
        pl.concat_list(batch_info).list.median().alias("INFO_MEDIAN"),
    ).select(
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "INFO_ORIGINAL",
        "INFO_MEAN",
        "INFO_MEDIAN",
    )
    invalid = transformed.filter(
        ~pl.col("CHROM").is_in(CHROMOSOMES)
        | pl.col("POS").is_null()
        | (pl.col("POS") < 1)
        | pl.col("REF").is_null()
        | (pl.col("REF").str.len_chars() < 1)
        | pl.col("ALT").is_null()
        | (pl.col("ALT").str.len_chars() < 1)
    )
    if invalid.height:
        raise FinnGenInfoError(
            f"batch contains {invalid.height} invalid coordinate or allele rows"
        )
    for column in ["INFO_ORIGINAL", "INFO_MEAN", "INFO_MEDIAN"]:
        out_of_range = transformed.select(
            ((pl.col(column) < 0) | (pl.col(column) > 1)).sum()
        ).item()
        if out_of_range:
            raise FinnGenInfoError(
                f"batch contains {out_of_range} {column} values outside [0, 1]"
            )
    return transformed


def statistic_frames(frame: pl.DataFrame) -> dict[str, pl.DataFrame]:
    coordinates = ["CHROM", "POS", "REF", "ALT"]
    return {
        "original": frame.select(
            *coordinates, pl.col("INFO_ORIGINAL").alias("INFO")
        ),
        "mean": frame.select(*coordinates, pl.col("INFO_MEAN").alias("INFO")),
        "median": frame.select(
            *coordinates, pl.col("INFO_MEDIAN").alias("INFO")
        ),
    }


def prepare_info_scores(
    source: Path = INPUT_FILE,
    output_directory: Path = OUTPUT_DIRECTORY,
) -> dict[str, int]:
    if not source.is_file() or source.stat().st_size == 0:
        raise FinnGenInfoError(f"input file is missing or empty: {source}")
    bgzip = executable(BGZIP_EXECUTABLE)
    tabix = executable(TABIX_EXECUTABLE)
    output_directory.mkdir(parents=True, exist_ok=True)
    log_path = output_directory / LOG_FILENAME
    header = read_header(source)
    required_columns, batch_info = select_columns(header)
    log_event(
        log_path,
        "run_started",
        source=str(source),
        source_bytes=source.stat().st_size,
        batch_info_columns=len(batch_info),
        batch_rows=BATCH_ROWS,
        polars_threads=POLARS_THREADS,
        bgzip_threads=BGZIP_THREADS,
    )

    schema = {
        "chr": pl.String,
        "pos": pl.Int64,
        "ref": pl.String,
        "alt": pl.String,
        ORIGINAL_INFO_COLUMN: pl.Float64,
        **{column: pl.Float64 for column in batch_info},
    }
    lazy = pl.scan_csv(
        source,
        separator="\t",
        schema_overrides=schema,
        null_values=NULL_VALUES,
        infer_schema_length=10_000,
        low_memory=True,
    ).select(required_columns)
    rank = {chrom: index for index, chrom in enumerate(CHROMOSOMES)}
    previous_chromosome: str | None = None
    previous_position: int | None = None
    writers: ChromosomeWriters | None = None
    counts = {chrom: 0 for chrom in CHROMOSOMES}
    try:
        for raw_batch in lazy.collect_batches(
            chunk_size=BATCH_ROWS,
            maintain_order=True,
            engine="streaming",
        ):
            batch = transform_batch(raw_batch, batch_info)
            for chrom_value in batch["CHROM"].unique(maintain_order=True):
                chrom = str(chrom_value)
                part = batch.filter(pl.col("CHROM") == chrom)
                if previous_chromosome is not None:
                    if rank[chrom] < rank[previous_chromosome]:
                        raise FinnGenInfoError(
                            f"chromosome order decreased from {previous_chromosome} "
                            f"to {chrom}"
                        )
                    if chrom == previous_chromosome and previous_position is not None:
                        first_position = int(part["POS"][0])
                        if first_position < previous_position:
                            raise FinnGenInfoError(
                                f"positions are unsorted on chromosome {chrom}"
                            )
                if chrom != previous_chromosome:
                    if writers is not None:
                        writers.finish(tabix, bgzip)
                        log_event(
                            log_path,
                            "chromosome_completed",
                            chromosome=previous_chromosome,
                            rows=counts[previous_chromosome],
                        )
                    writers = ChromosomeWriters(output_directory, chrom, bgzip)
                    previous_position = None
                positions = part["POS"]
                if not positions.is_sorted():
                    raise FinnGenInfoError(
                        f"positions are unsorted within a batch on chromosome {chrom}"
                    )
                if writers is None:
                    raise FinnGenInfoError("output writers were not initialized")
                writers.write(statistic_frames(part))
                counts[chrom] += part.height
                previous_chromosome = chrom
                previous_position = int(positions[-1])
        if writers is not None and previous_chromosome is not None:
            writers.finish(tabix, bgzip)
            log_event(
                log_path,
                "chromosome_completed",
                chromosome=previous_chromosome,
                rows=counts[previous_chromosome],
            )
    except BaseException as error:
        if writers is not None:
            writers.abort()
        log_event(
            log_path,
            "run_completed",
            status="failure",
            error_type=type(error).__name__,
            error=str(error),
        )
        raise
    if not any(counts.values()):
        raise FinnGenInfoError("input contains no data rows")
    log_event(log_path, "run_completed", status="success", rows_by_chromosome=counts)
    return counts


def main() -> int:
    counts = prepare_info_scores()
    print(
        f"[{timestamp()}] Completed FinnGen R13 INFO tracks: "
        f"{sum(counts.values()):,} variants",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FinnGenInfoError, OSError, pl.exceptions.PolarsError, subprocess.CalledProcessError) as error:
        print(f"[{timestamp()}] ERROR: {error}", file=sys.stderr, flush=True)
        raise SystemExit(1) from error
