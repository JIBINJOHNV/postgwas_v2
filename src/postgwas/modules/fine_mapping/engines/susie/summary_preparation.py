"""Prepare exact, boundary-preserving summary-statistics inputs for SuSiE.

The source table is streamed once into chromosome partitions. Each partition is
then loaded once and range-indexed to write one exact file per effective locus.
Overlapping primary loci intentionally receive overlapping variants: merging is
reserved for the configured post-primary overlap-resolution round.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd


# These columns are the scientific input contract consumed by the SuSiE model.
REQUIRED_SUMSTAT_COLUMNS = ("CHR", "BP", "SNP", "REF", "ALT", "EZ", "NEF", "LP")
LOCUS_KEY_COLUMNS = ("CHR", "START", "END")
CHROMOSOME_MANIFEST_COLUMNS = (
    "chromosome",
    "filename",
    "n_variants",
    "source_summary_statistics",
)
LOCUS_MANIFEST_COLUMNS = (
    "CHR",
    "START",
    "END",
    "GenomicLocus",
    "Filename",
    "n_variants",
    "chromosome_partition",
    "preparation_policy",
    "source_summary_statistics",
)


def normalize_chromosome(values: pd.Series) -> pd.Series:
    """Match the chromosome normalization used by the R SuSiE implementation."""
    normalized = values.astype("string").str.replace(
        r"^chr", "", regex=True, case=False
    ).str.upper().str.replace(r"\.0$", "", regex=True)
    return normalized.replace(
        {"X": "23", "Y": "24", "XY": "25", "M": "26", "MT": "26"}
    )


def _read_loci(locus_files: list[str | Path]) -> pd.DataFrame:
    frames = [pd.read_csv(path, sep="\t") for path in locus_files]
    loci = pd.concat(frames, ignore_index=True)
    missing = set(LOCUS_KEY_COLUMNS) - set(loci.columns)
    if missing:
        raise ValueError(
            "Prepared SuSiE locus chunks are missing columns: "
            + ", ".join(sorted(missing))
        )
    loci["CHR"] = normalize_chromosome(loci["CHR"])
    for column in ("START", "END"):
        loci[column] = pd.to_numeric(loci[column], errors="coerce")
    invalid = (
        loci["CHR"].isna()
        | loci["CHR"].eq("")
        | loci["START"].isna()
        | loci["END"].isna()
        | loci["START"].lt(1)
        | loci["END"].lt(loci["START"])
    )
    if invalid.any():
        raise ValueError(
            f"Prepared SuSiE locus chunks contain {int(invalid.sum())} invalid interval(s)"
        )
    loci[["START", "END"]] = loci[["START", "END"]].astype("int64")
    if "GenomicLocus" not in loci:
        loci["GenomicLocus"] = (
            "chr" + loci["CHR"] + ":" + loci["START"].astype(str)
            + "-" + loci["END"].astype(str)
        )
    invalid_names = (
        loci["GenomicLocus"].isna()
        | loci["GenomicLocus"].astype(str).str.strip().eq("")
    )
    if invalid_names.any():
        raise ValueError("Prepared SuSiE loci must have non-empty GenomicLocus values")
    if loci[list(LOCUS_KEY_COLUMNS)].duplicated().any():
        raise ValueError("Prepared SuSiE locus chunks contain duplicate intervals")
    return loci


def _validate_partition(
    table: pd.DataFrame, chromosome: str, source: Path
) -> pd.DataFrame:
    missing = set(REQUIRED_SUMSTAT_COLUMNS) - set(table.columns)
    if missing:
        raise ValueError(
            f"Summary statistics {source} are missing required columns: "
            + ", ".join(sorted(missing))
        )
    table = table.copy()
    table["CHR"] = normalize_chromosome(table["CHR"])
    for column in ("BP", "EZ", "NEF", "LP"):
        table[column] = pd.to_numeric(table[column], errors="coerce")
    table["SNP"] = table["SNP"].astype("string")
    table["REF"] = table["REF"].astype("string").str.upper()
    table["ALT"] = table["ALT"].astype("string").str.upper()
    invalid = (
        table["CHR"].isna()
        | table["CHR"].ne(chromosome)
        | table["BP"].isna()
        | table["BP"].lt(1)
        | table["SNP"].isna()
        | table["SNP"].eq("")
        | table["REF"].isna()
        | table["REF"].eq("")
        | table["ALT"].isna()
        | table["ALT"].eq("")
        | ~np.isfinite(table["EZ"])
        | ~np.isfinite(table["NEF"])
        | table["NEF"].le(0)
        | ~np.isfinite(table["LP"])
    )
    if invalid.any():
        raise ValueError(
            f"Summary statistics contain {int(invalid.sum())} invalid row(s) "
            f"on analysed chromosome {chromosome}"
        )
    if table["SNP"].duplicated().any():
        raise ValueError(
            f"Summary statistics contain duplicate SNP identifiers on chromosome {chromosome}"
        )
    if table[["CHR", "BP", "REF", "ALT"]].duplicated().any():
        raise ValueError(
            "Summary statistics contain duplicate chromosome/position/allele "
            f"variants on chromosome {chromosome}"
        )
    table["BP"] = table["BP"].astype("int64")
    return table.sort_values(["BP", "SNP"], kind="mergesort").reset_index(drop=True)


def _load_existing_cache(
    manifest_path: Path, chromosomes: set[str], source: Path
) -> pd.DataFrame:
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
    missing = set(CHROMOSOME_MANIFEST_COLUMNS) - set(manifest.columns)
    if missing:
        raise ValueError(
            "Chromosome summary-statistics manifest is missing columns: "
            + ", ".join(sorted(missing))
        )
    selected = manifest[manifest["chromosome"].isin(chromosomes)].copy()
    if manifest["chromosome"].duplicated().any():
        raise ValueError("Chromosome summary-statistics manifest has duplicate rows")
    recorded_sources = {
        str(Path(value).resolve()) for value in selected["source_summary_statistics"]
    }
    if recorded_sources != {str(source.resolve())}:
        raise ValueError(
            "Chromosome cache was prepared from different summary statistics"
        )
    absent = chromosomes - set(selected["chromosome"])
    if absent:
        raise ValueError(
            "Chromosome cache does not cover joint loci on chromosome(s): "
            + ", ".join(sorted(absent))
        )
    missing_files = [
        value for value in selected["filename"] if not Path(value).is_file()
    ]
    if missing_files:
        raise FileNotFoundError(
            "Chromosome cache manifest references missing file(s): "
            + ", ".join(missing_files)
        )
    return selected


def _stream_chromosome_cache(
    source: Path,
    cache_dir: Path,
    chromosomes: set[str],
    settings: dict,
) -> tuple[pd.DataFrame, int]:
    if cache_dir.exists():
        shutil.rmtree(cache_dir)
    cache_dir.mkdir(parents=True)
    paths = {
        chromosome: cache_dir / settings["chromosome_filename_template"].format(
            chromosome=chromosome
        )
        for chromosome in chromosomes
    }
    counts = {chromosome: 0 for chromosome in chromosomes}
    source_rows = 0
    saw_columns = False
    for chunk in pd.read_csv(
        source,
        sep=settings["separator"],
        chunksize=settings["chunk_rows"],
        compression=(
            None if settings["compression"] == "none" else settings["compression"]
        ),
    ):
        source_rows += len(chunk)
        missing = set(REQUIRED_SUMSTAT_COLUMNS) - set(chunk.columns)
        if missing:
            raise ValueError(
                f"Summary statistics {source} are missing required columns: "
                + ", ".join(sorted(missing))
            )
        saw_columns = True
        chunk = chunk.copy()
        chunk["CHR"] = normalize_chromosome(chunk["CHR"])
        for chromosome, subset in chunk[chunk["CHR"].isin(chromosomes)].groupby(
            "CHR", sort=False
        ):
            chromosome = str(chromosome)
            subset.to_csv(
                paths[chromosome],
                sep=settings["separator"],
                index=False,
                mode="a",
                header=not paths[chromosome].exists(),
            )
            counts[chromosome] += len(subset)
    if not saw_columns:
        raise ValueError(f"Summary-statistics file contains no variants: {source}")
    if source_rows == 0:
        raise ValueError(f"Summary-statistics file contains no variants: {source}")
    for chromosome, path in paths.items():
        if not path.exists():
            pd.DataFrame(columns=REQUIRED_SUMSTAT_COLUMNS).to_csv(
                path, sep=settings["separator"], index=False
            )
    manifest = pd.DataFrame(
        [
            {
                "chromosome": chromosome,
                "filename": str(paths[chromosome].resolve()),
                "n_variants": counts[chromosome],
                "source_summary_statistics": str(source.resolve()),
            }
            for chromosome in sorted(chromosomes)
        ],
        columns=CHROMOSOME_MANIFEST_COLUMNS,
    )
    return manifest, source_rows


def prepare_locus_summary_statistics(
    *,
    source_file: str | Path,
    locus_files: list[str | Path],
    chromosome_cache_directory: str | Path,
    locus_input_directory: str | Path,
    manifest_directory: str | Path,
    quality_control_directory: str | Path,
    settings: dict,
    existing_chromosome_manifest: str | Path | None = None,
) -> dict:
    """Write exact per-locus inputs and return their validated manifest paths."""
    source = Path(source_file).resolve()
    if not source.is_file():
        raise FileNotFoundError(f"SuSiE summary-statistics file is missing: {source}")
    loci = _read_loci(locus_files)
    cache_dir = Path(chromosome_cache_directory).resolve()
    locus_dir = Path(locus_input_directory).resolve()
    manifest_dir = Path(manifest_directory).resolve()
    quality_control_dir = Path(quality_control_directory).resolve()
    manifest_dir.mkdir(parents=True, exist_ok=True)
    quality_control_dir.mkdir(parents=True, exist_ok=True)
    if locus_dir.exists():
        shutil.rmtree(locus_dir)
    locus_dir.mkdir(parents=True)
    chromosomes = set(loci["CHR"].astype(str))
    cache_reused = existing_chromosome_manifest is not None
    source_rows = pd.NA
    if cache_reused:
        chromosome_manifest_path = Path(existing_chromosome_manifest).resolve()
        chromosome_manifest = _load_existing_cache(
            chromosome_manifest_path, chromosomes, source
        )
    else:
        chromosome_manifest, source_rows = _stream_chromosome_cache(
            source, cache_dir, chromosomes, settings
        )
        chromosome_manifest_path = (
            manifest_dir / settings["chromosome_manifest_filename"]
        )
        chromosome_manifest.to_csv(chromosome_manifest_path, sep="\t", index=False)

    partition_by_chromosome = {
        str(row.chromosome): row
        for row in chromosome_manifest.itertuples(index=False)
    }
    locus_rows = []
    total_memberships = 0
    for chromosome, chromosome_loci in loci.groupby("CHR", sort=False):
        partition_record = partition_by_chromosome[str(chromosome)]
        partition_path = Path(partition_record.filename)
        partition = pd.read_csv(partition_path, sep=settings["separator"])
        partition = _validate_partition(partition, str(chromosome), source)
        if len(partition) != int(partition_record.n_variants):
            raise ValueError(
                "Chromosome summary-statistics row count differs from manifest "
                f"for chromosome {chromosome}"
            )
        positions = partition["BP"].to_numpy()
        for locus in chromosome_loci.itertuples(index=False):
            locus_index = len(locus_rows) + 1
            locus_path = locus_dir / settings["locus_filename_template"].format(
                index=locus_index
            )
            left = int(np.searchsorted(positions, int(locus.START), side="left"))
            right = int(np.searchsorted(positions, int(locus.END), side="right"))
            selected = partition.iloc[left:right]
            selected.to_csv(locus_path, sep=settings["separator"], index=False)
            total_memberships += len(selected)
            locus_rows.append(
                {
                    "CHR": str(chromosome),
                    "START": int(locus.START),
                    "END": int(locus.END),
                    "GenomicLocus": str(locus.GenomicLocus),
                    "Filename": str(locus_path.resolve()),
                    "n_variants": len(selected),
                    "chromosome_partition": str(partition_path.resolve()),
                    "preparation_policy": settings["policy"],
                    "source_summary_statistics": str(source),
                }
            )
    locus_manifest_path = manifest_dir / settings["locus_manifest_filename"]
    pd.DataFrame(locus_rows, columns=LOCUS_MANIFEST_COLUMNS).to_csv(
        locus_manifest_path, sep="\t", index=False
    )
    summary_path = quality_control_dir / settings["preparation_summary_filename"]
    cached_variants = int(
        pd.to_numeric(chromosome_manifest["n_variants"], errors="raise").sum()
    )
    pd.DataFrame(
        [
            {
                "policy": settings["policy"],
                "source_read_policy": settings["source_read_policy"],
                "validation_scope": settings["validation_scope"],
                "source_summary_statistics": str(source),
                "source_rows_scanned": source_rows,
                "chromosome_cache_variants": cached_variants,
                "source_rows_outside_analysed_chromosomes": (
                    pd.NA if pd.isna(source_rows) else int(source_rows) - cached_variants
                ),
                "analysed_chromosomes": len(chromosomes),
                "prepared_loci": len(loci),
                "locus_variant_memberships": total_memberships,
                "chromosome_cache_reused": cache_reused,
                "chromosome_manifest": str(chromosome_manifest_path),
                "locus_manifest": str(locus_manifest_path),
            }
        ]
    ).to_csv(summary_path, sep="\t", index=False)
    return {
        "chromosome_manifest": str(chromosome_manifest_path),
        "locus_manifest": str(locus_manifest_path),
        "preparation_summary": str(summary_path),
        "n_loci": len(loci),
        "n_variant_memberships": total_memberships,
        "cache_reused": cache_reused,
    }
