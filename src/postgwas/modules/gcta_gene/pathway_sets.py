#!/usr/bin/env python3
"""Map GMT gene pathways to exact PLINK BIM IDs for GCTA fastBAT sets."""

from __future__ import annotations

import argparse
from bisect import bisect_left, bisect_right
from concurrent.futures import ProcessPoolExecutor
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import errno
import heapq
import io
import math
from pathlib import Path
import pickle
import shlex
import shutil
import sys
import tempfile
import time
from collections.abc import Callable, Iterator

import pyarrow as pa
import pyarrow.parquet as pq
import yaml

from postgwas.core.resource_preparation import (
    ResourcePreparationError,
    create_staging_directory,
    required_file,
    sha256,
    unique_nonempty_values,
    validate_output_filename,
    write_text,
)
from postgwas.core.ui import MeasuredProgress, StageProgress


OFFICIAL_FASTBAT_DOCUMENTATION = (
    "https://yanglab.westlake.edu.cn/software/gcta/#fastBAT"
)


@dataclass(frozen=True)
class Pathway:
    index: int
    name: str
    description: str
    genes: tuple[str, ...]
    duplicate_gene_entries: int


@dataclass(frozen=True)
class GeneInterval:
    chromosome: str
    start: int
    end: int
    gene: str


@dataclass(frozen=True)
class ChromosomeMappingTask:
    chromosome: str
    bim_cache: Path
    intervals: tuple[GeneInterval, ...]
    result_cache: Path


def _normalize_chromosome(value: str, policy: str) -> str:
    if policy == "strip_chr_prefix" and value.lower().startswith("chr"):
        return value[3:]
    return value


def _iter_gmt(path: Path, duplicate_gene_policy: str) -> Iterator[Pathway]:
    """Validate and yield one GMT pathway at a time in source order."""
    names: set[str] = set()
    pathway_index = 0
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                raise ResourcePreparationError(
                    f"GMT line {line_number} must contain pathway, description, "
                    "and at least one gene"
                )
            name = fields[0].strip()
            description = fields[1].strip()
            genes = [value.strip() for value in fields[2:]]
            if not name or any(character.isspace() for character in name):
                raise ResourcePreparationError(
                    f"GMT line {line_number} has an empty or whitespace-containing "
                    "pathway ID"
                )
            if name == "END":
                raise ResourcePreparationError("GMT pathway ID END is reserved by GCTA")
            if name in names:
                raise ResourcePreparationError(f"GMT repeats pathway ID {name!r}")
            if any(not gene for gene in genes):
                raise ResourcePreparationError(
                    f"GMT pathway {name!r} contains an empty gene identifier"
                )
            unique_genes = tuple(dict.fromkeys(genes))
            duplicate_count = len(genes) - len(unique_genes)
            if duplicate_count and duplicate_gene_policy == "error":
                raise ResourcePreparationError(
                    f"GMT pathway {name!r} contains {duplicate_count} duplicate "
                    "gene entries"
                )
            names.add(name)
            yield Pathway(
                pathway_index, name, description, unique_genes, duplicate_count,
            )
            pathway_index += 1
    if pathway_index == 0:
        raise ResourcePreparationError(f"GMT file contains no pathways: {path}")


def _scan_gmt(
    path: Path, duplicate_gene_policy: str,
) -> tuple[int, set[str]]:
    """Count pathways and collect the unique gene universe."""
    pathway_count = 0
    requested_genes: set[str] = set()
    for pathway in _iter_gmt(path, duplicate_gene_policy):
        requested_genes.update(pathway.genes)
        pathway_count += 1
    return pathway_count, requested_genes


def _read_gene_intervals(
    path: Path,
    requested_genes: set[str],
    allowed_chromosomes: set[str],
    chromosome_policy: str,
    window_bp: int,
) -> tuple[dict[str, GeneInterval], dict[str, list[GeneInterval]], set[str]]:
    genes: set[str] = set()
    requested: dict[str, GeneInterval] = {}
    requested_chromosomes: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.split()
            if len(fields) != 4:
                raise ResourcePreparationError(
                    f"gene-list line {line_number} has {len(fields)} fields; expected 4"
                )
            raw_chromosome, start_text, end_text, gene = fields
            chromosome = _normalize_chromosome(raw_chromosome, chromosome_policy)
            try:
                start, end = int(start_text), int(end_text)
            except ValueError as exc:
                raise ResourcePreparationError(
                    f"gene-list line {line_number} has non-integer coordinates"
                ) from exc
            if chromosome not in allowed_chromosomes:
                raise ResourcePreparationError(
                    f"gene-list line {line_number} has disallowed chromosome "
                    f"{raw_chromosome!r} after configured normalization"
                )
            if start < 1 or end < start or not gene:
                raise ResourcePreparationError(
                    f"gene-list line {line_number} has invalid coordinates or gene ID"
                )
            if gene in genes:
                raise ResourcePreparationError(
                    f"gene list repeats gene identifier {gene!r}"
                )
            genes.add(gene)
            if gene in requested_genes:
                requested[gene] = GeneInterval(
                    chromosome,
                    max(1, start - window_bp),
                    end + window_bp,
                    gene,
                )
                requested_chromosomes.add(chromosome)
    if not genes:
        raise ResourcePreparationError(f"gene list contains no genes: {path}")
    intervals: dict[str, list[GeneInterval]] = {}
    for interval in requested.values():
        intervals.setdefault(interval.chromosome, []).append(interval)
    for values in intervals.values():
        values.sort(key=lambda value: (value.start, value.end, value.gene))
    return requested, intervals, requested_chromosomes


def _split_bim_for_gene_mapping(
    bim: Path,
    intervals_by_chromosome: dict[str, list[GeneInterval]],
    allowed_chromosomes: set[str],
    chromosome_policy: str,
    workspace: Path,
    *,
    expected_variant_total: int | None,
    bim_ids_prevalidated: bool,
    analyzable_variant_ids: set[str] | None,
    progress_callback: Callable[[int], None] | None = None,
) -> tuple[dict[str, Path], list[str], dict]:
    """Validate the BIM once and cache only chromosomes containing requested genes."""
    current_chromosome: str | None = None
    completed_chromosomes: set[str] = set()
    last_position = 0
    variant_count = 0
    observed_chromosomes: set[str] = set()
    cache_paths: dict[str, Path] = {}
    cache_handles: dict[str, io.TextIOWrapper] = {}
    cached_variant_ids: list[str] = []
    seen_variant_ids: set[str] | None = set() if not bim_ids_prevalidated else None
    try:
        with bim.open("r", encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, 1):
                if not raw.strip():
                    raise ResourcePreparationError(
                        f"PLINK BIM contains a blank line at {line_number}"
                    )
                fields = raw.split()
                if len(fields) != 6:
                    raise ResourcePreparationError(
                        f"PLINK BIM line {line_number} has {len(fields)} fields; "
                        "expected 6"
                    )
                raw_chromosome, variant_id, _, position_text, _, _ = fields
                chromosome = _normalize_chromosome(
                    raw_chromosome, chromosome_policy,
                )
                if chromosome not in allowed_chromosomes:
                    raise ResourcePreparationError(
                        f"PLINK BIM line {line_number} has disallowed chromosome "
                        f"{raw_chromosome!r} after configured normalization"
                    )
                if not variant_id or variant_id in {".", "END"}:
                    raise ResourcePreparationError(
                        f"PLINK BIM line {line_number} has an empty, '.', or "
                        "reserved END variant ID"
                    )
                if seen_variant_ids is not None:
                    if variant_id in seen_variant_ids:
                        raise ResourcePreparationError(
                            f"PLINK BIM repeats variant ID {variant_id!r} at line "
                            f"{line_number}"
                        )
                    seen_variant_ids.add(variant_id)
                try:
                    position = int(position_text)
                except ValueError as exc:
                    raise ResourcePreparationError(
                        f"PLINK BIM line {line_number} has a non-integer position"
                    ) from exc
                if position < 1:
                    raise ResourcePreparationError(
                        f"PLINK BIM line {line_number} has a non-positive position"
                    )
                if chromosome != current_chromosome:
                    if chromosome in completed_chromosomes:
                        raise ResourcePreparationError(
                            "PLINK BIM chromosome blocks are not contiguous; sort "
                            "the BIM before mapping (repeated chromosome "
                            f"{raw_chromosome!r})"
                        )
                    if current_chromosome is not None:
                        completed_chromosomes.add(current_chromosome)
                    current_chromosome = chromosome
                    last_position = 0
                if position < last_position:
                    raise ResourcePreparationError(
                        "PLINK BIM positions decrease within chromosome "
                        f"{raw_chromosome!r} at line {line_number}; sort the BIM "
                        "before mapping"
                    )
                last_position = position
                observed_chromosomes.add(chromosome)
                variant_count += 1
                if (
                    chromosome in intervals_by_chromosome
                    and (
                        analyzable_variant_ids is None
                        or variant_id in analyzable_variant_ids
                    )
                ):
                    if chromosome not in cache_handles:
                        cache = tempfile.NamedTemporaryFile(
                            mode="w",
                            encoding="utf-8",
                            dir=workspace,
                            delete=False,
                        )
                        cache_handles[chromosome] = cache
                        cache_paths[chromosome] = Path(cache.name)
                    variant_index = len(cached_variant_ids)
                    cached_variant_ids.append(variant_id)
                    cache_handles[chromosome].write(
                        f"{variant_index}\t{position}\n"
                    )
                if progress_callback is not None and (
                    expected_variant_total is None
                    or variant_count < expected_variant_total
                ):
                    progress_callback(variant_count)
    finally:
        for cache_handle in cache_handles.values():
            cache_handle.close()
    if variant_count == 0:
        raise ResourcePreparationError(f"PLINK BIM contains no variants: {bim}")
    return cache_paths, cached_variant_ids, {
        "bim_variants": variant_count,
        "analyzable_bim_variants_cached": len(cached_variant_ids),
        "bim_chromosomes": observed_chromosomes,
    }


def _map_one_chromosome(task: ChromosomeMappingTask) -> dict:
    """Create one compact gene-to-BIM cache using binary searches."""
    variant_indexes: list[int] = []
    positions: list[int] = []
    with task.bim_cache.open("r", encoding="utf-8") as handle:
        for raw in handle:
            index_text, position_text = raw.rstrip("\n").split("\t")
            variant_indexes.append(int(index_text))
            positions.append(int(position_text))

    gene_variants: dict[str, tuple[int, ...]] = {}
    genes_with_variants: list[str] = []
    for interval in task.intervals:
        left = bisect_left(positions, interval.start)
        right = bisect_right(positions, interval.end)
        variants = tuple(variant_indexes[left:right])
        gene_variants[interval.gene] = variants
        if variants:
            genes_with_variants.append(interval.gene)

    merged_intervals: list[tuple[int, int]] = []
    for interval in task.intervals:
        if not merged_intervals or interval.start > merged_intervals[-1][1]:
            merged_intervals.append((interval.start, interval.end))
        else:
            start, end = merged_intervals[-1]
            merged_intervals[-1] = (start, max(end, interval.end))
    mapped_variant_count = sum(
        bisect_right(positions, end) - bisect_left(positions, start)
        for start, end in merged_intervals
    )
    with task.result_cache.open("wb") as handle:
        pickle.dump(gene_variants, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return {
        "chromosome": task.chromosome,
        "mapped_variants": mapped_variant_count,
        "genes_with_variants": genes_with_variants,
    }


def _mapping_worker_plan(
    tasks: list[ChromosomeMappingTask],
    requested_workers: int,
    memory_gb: float,
    worker_memory_multiplier: float,
) -> tuple[int, dict]:
    if requested_workers < 1:
        raise ResourcePreparationError("mapping worker count must be positive")
    if memory_gb <= 0:
        raise ResourcePreparationError("mapping memory must be positive")
    if worker_memory_multiplier <= 1:
        raise ResourcePreparationError(
            "mapping worker-memory multiplier must be greater than one"
        )
    largest_cache = max(task.bim_cache.stat().st_size for task in tasks)
    estimated_worker_bytes = max(
        1, int(largest_cache * worker_memory_multiplier),
    )
    available_bytes = max(1, int(memory_gb * 1024**3))
    memory_workers = available_bytes // estimated_worker_bytes
    if memory_workers < 1:
        raise ResourcePreparationError(
            "configured mapping memory is insufficient for one chromosome "
            "worker: estimated %d bytes, available %d bytes"
            % (estimated_worker_bytes, available_bytes)
        )
    workers = min(len(tasks), requested_workers, memory_workers)
    return workers, {
        "requested_mapping_workers": requested_workers,
        "mapping_workers": workers,
        "largest_chromosome_cache_bytes": largest_cache,
        "estimated_worker_memory_bytes": estimated_worker_bytes,
        "mapping_memory_bytes": available_bytes,
    }


def _map_bim_variants(
    bim: Path,
    intervals_by_chromosome: dict[str, list[GeneInterval]],
    allowed_chromosomes: set[str],
    chromosome_policy: str,
    workspace: Path,
    *,
    expected_variant_total: int | None,
    requested_workers: int,
    memory_gb: float,
    worker_memory_multiplier: float,
    bim_ids_prevalidated: bool,
    analyzable_variant_ids: set[str] | None = None,
    progress_callback: Callable[[int], None] | None = None,
) -> tuple[dict[str, tuple[int, ...]], list[str], dict]:
    """Map BIM variants to genes before expanding reusable pathway memberships."""
    cache_paths, variant_ids, metrics = _split_bim_for_gene_mapping(
        bim,
        intervals_by_chromosome,
        allowed_chromosomes,
        chromosome_policy,
        workspace,
        expected_variant_total=expected_variant_total,
        bim_ids_prevalidated=bim_ids_prevalidated,
        analyzable_variant_ids=analyzable_variant_ids,
        progress_callback=progress_callback,
    )
    tasks: list[ChromosomeMappingTask] = []
    for chromosome, intervals in intervals_by_chromosome.items():
        bim_cache = cache_paths.get(chromosome)
        if bim_cache is None:
            continue
        result_handle = tempfile.NamedTemporaryFile(dir=workspace, delete=False)
        result_cache = Path(result_handle.name)
        result_handle.close()
        tasks.append(ChromosomeMappingTask(
            chromosome=chromosome,
            bim_cache=bim_cache,
            intervals=tuple(intervals),
            result_cache=result_cache,
        ))
    if not tasks:
        metrics.update({
            "bim_variants_mapped_to_requested_genes": 0,
            "genes_with_at_least_one_bim_variant": set(),
            "requested_mapping_workers": requested_workers,
            "mapping_workers": 0,
            "largest_chromosome_cache_bytes": 0,
            "estimated_worker_memory_bytes": 0,
            "mapping_memory_bytes": max(1, int(memory_gb * 1024**3)),
        })
        return {}, variant_ids, metrics

    workers, worker_metrics = _mapping_worker_plan(
        tasks, requested_workers, memory_gb, worker_memory_multiplier,
    )
    try:
        if workers == 1:
            chromosome_metrics = [_map_one_chromosome(task) for task in tasks]
        else:
            with ProcessPoolExecutor(max_workers=workers) as executor:
                chromosome_metrics = list(
                    executor.map(_map_one_chromosome, tasks)
                )
    except ResourcePreparationError:
        raise
    except Exception as exc:
        raise ResourcePreparationError(
            "BIM chromosome mapping failed with %d worker(s): %s"
            % (workers, exc)
        ) from exc

    gene_variants: dict[str, tuple[int, ...]] = {}
    for task in tasks:
        with task.result_cache.open("rb") as handle:
            chromosome_gene_variants = pickle.load(handle)
        duplicate_genes = set(gene_variants) & set(chromosome_gene_variants)
        if duplicate_genes:
            raise ResourcePreparationError(
                "internal gene cache repeats identifiers: "
                + ", ".join(sorted(duplicate_genes)[:10])
            )
        gene_variants.update(chromosome_gene_variants)
    metrics.update(worker_metrics)
    metrics.update({
        "bim_variants_mapped_to_requested_genes": sum(
            value["mapped_variants"] for value in chromosome_metrics
        ),
        "genes_with_at_least_one_bim_variant": {
            gene
            for value in chromosome_metrics
            for gene in value["genes_with_variants"]
        },
    })
    return gene_variants, variant_ids, metrics


def _iter_pathway_gene_variants(
    genes: list[str],
    gene_variants: dict[str, tuple[int, ...]],
    variant_ids: list[str],
):
    """Yield gene memberships in stable BIM order, then gene-ID order."""

    def tagged_variants(gene: str):
        for variant_index in gene_variants[gene]:
            yield variant_index, gene, variant_ids[variant_index]

    streams = (
        tagged_variants(gene)
        for gene in genes
    )
    return heapq.merge(*streams)


def _pathway_gene_groups(
    pathway: Pathway,
    gene_coordinates: dict[str, GeneInterval],
    gene_variants: dict[str, tuple[int, ...]],
) -> tuple[list[str], list[str], list[str]]:
    mapped = [gene for gene in pathway.genes if gene in gene_coordinates]
    missing = [gene for gene in pathway.genes if gene not in gene_coordinates]
    hit = [gene for gene in mapped if gene_variants.get(gene)]
    return mapped, missing, hit


def _unique_pathway_variant_indexes(
    genes: list[str], gene_variants: dict[str, tuple[int, ...]],
) -> Iterator[int]:
    """Yield deduplicated variant indexes in stable BIM order."""
    last_index: int | None = None
    for variant_index in heapq.merge(*(gene_variants[gene] for gene in genes)):
        if variant_index == last_index:
            continue
        yield variant_index
        last_index = variant_index


def _tsv_row_bytes(*values: object) -> int:
    return len(("\t".join(str(value) for value in values) + "\n").encode("utf-8"))


class _ParquetRowWriter:
    """Write string audit rows in configured bounded record batches."""

    def __init__(
        self,
        path: Path,
        columns: tuple[str, ...],
        *,
        compression: str,
        batch_rows: int,
    ):
        self.path = path
        self.columns = columns
        self.batch_rows = batch_rows
        self.schema = pa.schema([(column, pa.string()) for column in columns])
        self.buffers = {column: [] for column in columns}
        codec = None if compression == "none" else compression
        self.writer = pq.ParquetWriter(
            path,
            self.schema,
            compression=codec,
            use_dictionary=True,
        )

    def append(self, *values: object) -> None:
        if len(values) != len(self.columns):
            raise ResourcePreparationError(
                "internal pathway-audit row has the wrong number of columns"
            )
        for column, value in zip(self.columns, values):
            self.buffers[column].append(str(value))
        if len(self.buffers[self.columns[0]]) >= self.batch_rows:
            self.flush()

    def flush(self) -> None:
        if not self.buffers[self.columns[0]]:
            return
        table = pa.Table.from_pydict(self.buffers, schema=self.schema)
        self.writer.write_table(table, row_group_size=self.batch_rows)
        self.buffers = {column: [] for column in self.columns}

    def close(self) -> None:
        self.flush()
        self.writer.close()


def _table_text(fieldnames: list[str], rows: list[dict]) -> str:
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def _readme_text(manifest: dict) -> str:
    resource = manifest["resource"]
    validation = manifest["validation"]
    outputs = manifest["outputs"]
    return f"""# {resource['name']}

This directory contains a GCTA `--fastBAT-set-list` generated from a GMT
gene-set file, a build-matched gene-coordinate list, and a PLINK BIM LD
reference. Variant identifiers are copied verbatim from BIM column 2; they are
never inferred or converted to rsIDs.

Genome build: `{resource['genome_build']}`. The GMT mapping coordinates, BIM
reference, harmonised GWAS, and analysis configuration must use this same
build. A {resource['gene_window_kb']} kb window was applied to both sides of
each gene during gene-to-variant mapping.

Pathways written: {validation['pathways_written']}  
Pathways omitted: {validation['pathways_omitted']}  
GMT genes not found in the coordinate list: {validation['unmapped_gene_entries']}  
Unique BIM variants written across retained pathways: {validation['unique_variants_written']}

`{outputs['pathway_mapping']}` reports gene and variant counts for every input
pathway, and `{outputs['unmapped_genes']}` records each unresolved GMT gene.
The configured audit level is `{manifest['policies']['audit']['level']}`.
Normalized audit output stores pathway-to-gene and gene-to-BIM-ID relationships
once in Parquet instead of repeating them for every pathway-gene-variant
combination. The optional expanded audit is generated only when explicitly
configured. Full provenance and checksums are in `{outputs['manifest']}` and
`{outputs['checksums']}`.

This resource is for fastBAT custom-set analysis. GCTA does not document an
mBAT-combo option equivalent to `--fastBAT-set-list`.

Format source:

- {OFFICIAL_FASTBAT_DOCUMENTATION}
"""


def prepare_resource(
    args: argparse.Namespace,
    *,
    stage_progress: StageProgress | None = None,
    bim_variant_total: int | None = None,
    progress_refresh_seconds: float | None = None,
    mapping_workers: int,
    mapping_memory_gb: float,
    worker_memory_multiplier: float,
    bim_ids_prevalidated: bool = False,
    analyzable_variant_ids: set[str] | None = None,
    analysis_variant_source: Path | None = None,
    maximum_set_variants: int,
    oversized_set_policy: str,
    audit_level: str,
    audit_format: str,
    audit_compression: str,
    audit_batch_rows: int,
    minimum_free_disk_gb: float,
    disk_estimation_safety_factor: float,
) -> dict:
    progress = stage_progress or StageProgress(
        "GMT-to-fastBAT preparation", enabled=False,
    )
    with progress.step(1, 8, "Validate inputs and scan GMT pathways"):
        gmt = required_file(args.gmt, "GMT pathway file")
        gene_list = required_file(args.gene_list, "GCTA gene-coordinate file")
        bim = required_file(args.bim, "PLINK BIM file")
        analysis_source = (
            required_file(analysis_variant_source, "analyzable GWAS variant source")
            if analysis_variant_source is not None else None
        )
        if analyzable_variant_ids is not None:
            if not analyzable_variant_ids:
                raise ResourcePreparationError(
                    "analyzable GWAS/BIM variant universe is empty"
                )
            if analysis_source is None:
                raise ResourcePreparationError(
                    "an analyzable variant universe requires its GWAS source for "
                    "provenance and resume validation"
                )
        apply_set_size_policy = analyzable_variant_ids is not None
        output_directory = Path(args.output_directory).expanduser().resolve()
        metadata_names = {
            args.pathway_mapping_name,
            args.pathway_gene_mapping_name,
            args.gene_variant_mapping_name,
            args.expanded_mapping_name,
            args.unmapped_genes_name,
            args.manifest_name,
            args.readme_name,
            args.checksums_name,
        }
        validate_output_filename(args.output_name, metadata_names)
        for name in metadata_names:
            validate_output_filename(name, {args.output_name})
        if len(metadata_names) != 8:
            raise ResourcePreparationError("resource output filenames must be unique")
        if args.gene_window_kb < 0:
            raise ResourcePreparationError("--gene-window-kb must be nonnegative")
        if bim_variant_total is not None and int(bim_variant_total) < 1:
            raise ResourcePreparationError(
                "validated BIM progress total must be a positive integer"
            )
        if int(mapping_workers) < 1:
            raise ResourcePreparationError("mapping worker count must be positive")
        if float(mapping_memory_gb) <= 0:
            raise ResourcePreparationError("mapping memory must be positive")
        if float(worker_memory_multiplier) <= 1:
            raise ResourcePreparationError(
                "mapping worker-memory multiplier must be greater than one"
            )
        if int(maximum_set_variants) < 1:
            raise ResourcePreparationError(
                "maximum fastBAT set variants must be positive"
            )
        if oversized_set_policy not in {"error", "omit"}:
            raise ResourcePreparationError(
                "oversized fastBAT set policy must be 'error' or 'omit'"
            )
        if audit_level not in {"summary", "normalized", "expanded"}:
            raise ResourcePreparationError(
                "pathway audit level must be summary, normalized, or expanded"
            )
        if audit_format != "parquet":
            raise ResourcePreparationError("pathway audit format must be parquet")
        if audit_compression not in {"zstd", "snappy", "gzip", "none"}:
            raise ResourcePreparationError(
                "pathway audit compression must be zstd, snappy, gzip, or none"
            )
        if int(audit_batch_rows) < 1:
            raise ResourcePreparationError("pathway audit batch size must be positive")
        if float(minimum_free_disk_gb) < 0:
            raise ResourcePreparationError(
                "minimum free disk space must be nonnegative"
            )
        if float(disk_estimation_safety_factor) < 1:
            raise ResourcePreparationError(
                "disk estimation safety factor must be at least one"
            )
        if progress.enabled and (
            progress_refresh_seconds is None
            or float(progress_refresh_seconds) <= 0
        ):
            raise ResourcePreparationError(
                "enabled BIM mapping progress requires a positive refresh interval"
            )
        allowed_chromosomes = unique_nonempty_values(
            args.allowed_chromosomes, "--allowed-chromosomes",
        )
        pathway_count, requested_genes = _scan_gmt(
            gmt, args.duplicate_gene_policy,
        )

    with progress.step(2, 8, "Resolve genes and apply configured boundaries"):
        window_bp = args.gene_window_kb * 1000
        gene_coordinates, intervals, gene_chromosomes = _read_gene_intervals(
            gene_list,
            requested_genes,
            allowed_chromosomes,
            args.chromosome_label_policy,
            window_bp,
        )
        missing_genes = requested_genes - set(gene_coordinates)
        if missing_genes and args.unmapped_gene_policy == "error":
            examples = ", ".join(sorted(missing_genes)[:10])
            raise ResourcePreparationError(
                f"{len(missing_genes)} GMT genes are absent from the coordinate list; "
                f"examples: {examples}"
            )

    staging: Path | None = None
    active_output_path: Path | None = None
    completed_output_pathways = 0
    retained_pathways = 0
    disk_preflight: dict = {}
    try:
        with progress.step(3, 8, "Map analyzable BIM variants to genes"):
            staging = create_staging_directory(output_directory)
            mapping_workspace = Path(tempfile.mkdtemp(dir=staging))
            mapping_progress = MeasuredProgress(
                "BIM-to-pathway mapping progress",
                enabled=progress.enabled,
            )
            mapping_title = "Map analyzable BIM variants to genes"
            mapping_progress.start(mapping_title, total=bim_variant_total)
            last_refresh = time.monotonic()

            def update_mapping_progress(completed: int) -> None:
                nonlocal last_refresh
                now = time.monotonic()
                if now - last_refresh < float(progress_refresh_seconds):
                    return
                displayed = int(completed)
                if bim_variant_total is not None:
                    displayed = min(displayed, int(bim_variant_total))
                mapping_progress.update(
                    displayed,
                    total=bim_variant_total,
                    title=mapping_title,
                )
                last_refresh = now

            try:
                gene_variants, variant_ids, bim_metrics = _map_bim_variants(
                    bim,
                    intervals,
                    allowed_chromosomes,
                    args.chromosome_label_policy,
                    mapping_workspace,
                    expected_variant_total=bim_variant_total,
                    requested_workers=int(mapping_workers),
                    memory_gb=float(mapping_memory_gb),
                    worker_memory_multiplier=float(worker_memory_multiplier),
                    bim_ids_prevalidated=bim_ids_prevalidated,
                    analyzable_variant_ids=analyzable_variant_ids,
                    progress_callback=(
                        update_mapping_progress
                        if mapping_progress.enabled else None
                    ),
                )
                mapped_total = bim_metrics["bim_variants"]
                if (
                    bim_variant_total is not None
                    and mapped_total != int(bim_variant_total)
                ):
                    raise ResourcePreparationError(
                        "PLINK BIM variant count changed after reference validation: "
                        f"expected {int(bim_variant_total)}, observed {mapped_total}"
                    )
                shared_chromosomes = (
                    gene_chromosomes & bim_metrics["bim_chromosomes"]
                )
                if not shared_chromosomes:
                    raise ResourcePreparationError(
                        "gene coordinates and PLINK BIM do not share chromosome labels"
                    )
                if mapping_progress.enabled:
                    update_mapping_progress(mapped_total)
            except BaseException:
                mapping_progress.fail(title=mapping_title)
                raise
            else:
                shutil.rmtree(mapping_workspace)
                mapping_progress.complete(
                    mapped_total,
                    total=mapped_total,
                    title="Validate BIM-to-pathway mapping",
                )

        with progress.step(4, 8, "Validate pathway-to-variant memberships"):
            mapping_rows: list[dict] = []
            unmapped_entry_count = 0
            empty: list[str] = []
            for pathway in _iter_gmt(gmt, args.duplicate_gene_policy):
                mapped_genes, missing, hit_genes = _pathway_gene_groups(
                    pathway, gene_coordinates, gene_variants,
                )
                status = "retained" if hit_genes else "omitted_empty"
                if not hit_genes:
                    empty.append(pathway.name)
                mapping_rows.append({
                    "pathway": pathway.name,
                    "description": pathway.description,
                    "input_genes": len(pathway.genes),
                    "duplicate_gene_entries_removed": pathway.duplicate_gene_entries,
                    "genes_with_coordinates": len(mapped_genes),
                    "genes_with_analyzable_variants": len(hit_genes),
                    "unmapped_genes": len(missing),
                    "unique_analyzable_variants": 0,
                    "status": status,
                })
                unmapped_entry_count += len(missing)
            if empty and args.empty_pathway_policy == "error":
                examples = ", ".join(empty[:10])
                raise ResourcePreparationError(
                    f"{len(empty)} pathways have no analyzable variants; examples: "
                    f"{examples}"
                )
            if len(empty) == pathway_count:
                raise ResourcePreparationError(
                    "no pathways retain analyzable BIM variants"
                )

        with progress.step(5, 8, "Preflight pathway sizes and disk requirements"):
            retained_genes: set[str] = set()
            unique_variant_flags = bytearray(len(variant_ids))
            oversized: list[tuple[str, int]] = []
            total_memberships = 0
            total_gene_variant_memberships = 0
            overlapping_memberships_collapsed = 0
            set_file_bytes = 0
            pathway_gene_rows = 0
            pathway_gene_raw_bytes = _tsv_row_bytes("pathway", "gene")
            expanded_rows = 0
            expanded_raw_bytes = _tsv_row_bytes(
                "pathway", "gene", "variant_id",
            )
            for pathway in _iter_gmt(gmt, args.duplicate_gene_policy):
                row = mapping_rows[pathway.index]
                if row["status"] != "retained":
                    continue
                _, _, hit_genes = _pathway_gene_groups(
                    pathway, gene_coordinates, gene_variants,
                )
                unique_indexes = list(_unique_pathway_variant_indexes(
                    hit_genes, gene_variants,
                ))
                unique_count = len(unique_indexes)
                row["unique_analyzable_variants"] = unique_count
                if (
                    apply_set_size_policy
                    and unique_count > int(maximum_set_variants)
                ):
                    row["status"] = "omitted_oversized"
                    oversized.append((pathway.name, unique_count))
                    continue
                retained_pathways += 1
                retained_genes.update(hit_genes)
                total_memberships += unique_count
                gene_memberships = sum(
                    len(gene_variants[gene]) for gene in hit_genes
                )
                total_gene_variant_memberships += gene_memberships
                overlapping_memberships_collapsed += (
                    gene_memberships - unique_count
                )
                set_file_bytes += len((pathway.name + "\n").encode("utf-8"))
                set_file_bytes += sum(
                    len((variant_ids[index] + "\n").encode("utf-8"))
                    for index in unique_indexes
                )
                set_file_bytes += len("END\n\n".encode("utf-8"))
                for index in unique_indexes:
                    unique_variant_flags[index] = 1
                if audit_level in {"normalized", "expanded"}:
                    for gene in hit_genes:
                        pathway_gene_rows += 1
                        pathway_gene_raw_bytes += _tsv_row_bytes(
                            pathway.name, gene,
                        )
                if audit_level == "expanded":
                    for _, gene, variant_id in _iter_pathway_gene_variants(
                        hit_genes, gene_variants, variant_ids,
                    ):
                        expanded_rows += 1
                        expanded_raw_bytes += _tsv_row_bytes(
                            pathway.name, gene, variant_id,
                        )
            if oversized and oversized_set_policy == "error":
                raise ResourcePreparationError(
                    "%d pathways exceed GCTA's configured %s-variant set limit; "
                    "examples: %s"
                    % (
                        len(oversized),
                        format(int(maximum_set_variants), ","),
                        ", ".join("%s (%d)" % item for item in oversized[:10]),
                    )
                )
            if retained_pathways == 0:
                raise ResourcePreparationError(
                    "no pathways remain after analyzable-variant and set-size validation"
                )

            gene_variant_rows = 0
            gene_variant_raw_bytes = _tsv_row_bytes("gene", "variant_id")
            if audit_level in {"normalized", "expanded"}:
                for gene in retained_genes:
                    for variant_index in gene_variants[gene]:
                        gene_variant_rows += 1
                        gene_variant_raw_bytes += _tsv_row_bytes(
                            gene, variant_ids[variant_index],
                        )
            mapping_columns = [
                "pathway", "description", "input_genes",
                "duplicate_gene_entries_removed", "genes_with_coordinates",
                "genes_with_analyzable_variants", "unmapped_genes",
                "unique_analyzable_variants", "status",
            ]
            mapping_text = _table_text(mapping_columns, mapping_rows)
            unmapped_raw_bytes = _tsv_row_bytes("pathway", "gene")
            for pathway in _iter_gmt(gmt, args.duplicate_gene_policy):
                _, missing, _ = _pathway_gene_groups(
                    pathway, gene_coordinates, gene_variants,
                )
                unmapped_raw_bytes += sum(
                    _tsv_row_bytes(pathway.name, gene) for gene in missing
                )
            audit_estimated_bytes = 0
            if audit_level in {"normalized", "expanded"}:
                audit_estimated_bytes += (
                    pathway_gene_raw_bytes + gene_variant_raw_bytes
                )
            if audit_level == "expanded":
                audit_estimated_bytes += expanded_raw_bytes
            estimated_output_bytes = (
                set_file_bytes
                + len(mapping_text.encode("utf-8"))
                + unmapped_raw_bytes
                + audit_estimated_bytes
            )
            free_bytes = shutil.disk_usage(staging).free
            minimum_free_bytes = math.ceil(
                float(minimum_free_disk_gb) * 1024**3
            )
            required_free_bytes = (
                math.ceil(
                    estimated_output_bytes
                    * float(disk_estimation_safety_factor)
                )
                + minimum_free_bytes
            )
            disk_preflight = {
                "estimated_output_bytes": estimated_output_bytes,
                "estimated_audit_uncompressed_bytes": audit_estimated_bytes,
                "available_bytes": free_bytes,
                "minimum_free_bytes_after_write": minimum_free_bytes,
                "estimation_safety_factor": float(
                    disk_estimation_safety_factor
                ),
                "required_available_bytes": required_free_bytes,
                "retained_pathways": retained_pathways,
                "passed": free_bytes >= required_free_bytes,
            }
            if not disk_preflight["passed"]:
                raise ResourcePreparationError(
                    "insufficient disk space for GMT-to-fastBAT outputs: "
                    "estimated outputs=%s bytes, required available space=%s "
                    "bytes including the configured safety margin, available=%s "
                    "bytes, audit level=%s"
                    % (
                        format(estimated_output_bytes, ","),
                        format(required_free_bytes, ","),
                        format(free_bytes, ","),
                        audit_level,
                    )
                )

        with progress.step(6, 8, "Write fastBAT sets and normalized audit tables"):
            output_file = staging / args.output_name
            mapping_file = staging / args.pathway_mapping_name
            pathway_gene_file = (
                staging / args.pathway_gene_mapping_name
                if audit_level in {"normalized", "expanded"} else None
            )
            gene_variant_file = (
                staging / args.gene_variant_mapping_name
                if audit_level in {"normalized", "expanded"} else None
            )
            expanded_file = (
                staging / args.expanded_mapping_name
                if audit_level == "expanded" else None
            )
            unmapped_file = staging / args.unmapped_genes_name
            relationship_writers: list[_ParquetRowWriter] = []
            pathway_gene_writer = None
            gene_variant_writer = None
            expanded_writer = None
            if pathway_gene_file is not None:
                active_output_path = pathway_gene_file
                pathway_gene_writer = _ParquetRowWriter(
                    pathway_gene_file,
                    ("pathway", "gene"),
                    compression=audit_compression,
                    batch_rows=int(audit_batch_rows),
                )
                relationship_writers.append(pathway_gene_writer)
            if gene_variant_file is not None:
                active_output_path = gene_variant_file
                gene_variant_writer = _ParquetRowWriter(
                    gene_variant_file,
                    ("gene", "variant_id"),
                    compression=audit_compression,
                    batch_rows=int(audit_batch_rows),
                )
                relationship_writers.append(gene_variant_writer)
            if expanded_file is not None:
                active_output_path = expanded_file
                expanded_writer = _ParquetRowWriter(
                    expanded_file,
                    ("pathway", "gene", "variant_id"),
                    compression=audit_compression,
                    batch_rows=int(audit_batch_rows),
                )
                relationship_writers.append(expanded_writer)
            writing_progress = MeasuredProgress(
                "fastBAT pathway writing progress",
                enabled=progress.enabled,
            )
            writing_title = "Write validated pathway sets"
            writing_progress.start(writing_title, total=retained_pathways)
            last_write_refresh = time.monotonic()
            try:
                active_output_path = output_file
                with output_file.open("w", encoding="utf-8") as set_handle:
                    for pathway in _iter_gmt(gmt, args.duplicate_gene_policy):
                        if mapping_rows[pathway.index]["status"] != "retained":
                            continue
                        _, _, hit_genes = _pathway_gene_groups(
                            pathway, gene_coordinates, gene_variants,
                        )
                        unique_indexes = list(_unique_pathway_variant_indexes(
                            hit_genes, gene_variants,
                        ))
                        expected_count = mapping_rows[pathway.index][
                            "unique_analyzable_variants"
                        ]
                        if len(unique_indexes) != expected_count:
                            raise ResourcePreparationError(
                                "pathway variant membership changed after preflight: "
                                + pathway.name
                            )
                        active_output_path = output_file
                        set_handle.write(pathway.name + "\n")
                        set_handle.writelines(
                            variant_ids[index] + "\n" for index in unique_indexes
                        )
                        set_handle.write("END\n\n")
                        if pathway_gene_writer is not None:
                            active_output_path = pathway_gene_file
                            for gene in hit_genes:
                                pathway_gene_writer.append(pathway.name, gene)
                        if expanded_writer is not None:
                            active_output_path = expanded_file
                            for _, gene, variant_id in _iter_pathway_gene_variants(
                                hit_genes, gene_variants, variant_ids,
                            ):
                                expanded_writer.append(
                                    pathway.name, gene, variant_id,
                                )
                        completed_output_pathways += 1
                        now = time.monotonic()
                        if (
                            writing_progress.enabled
                            and now - last_write_refresh
                            >= float(progress_refresh_seconds)
                        ):
                            writing_progress.update(
                                completed_output_pathways,
                                total=retained_pathways,
                                title=writing_title,
                            )
                            last_write_refresh = now
                if gene_variant_writer is not None:
                    active_output_path = gene_variant_file
                    for gene in sorted(retained_genes):
                        for variant_index in gene_variants[gene]:
                            gene_variant_writer.append(
                                gene, variant_ids[variant_index],
                            )
            except BaseException:
                for writer in relationship_writers:
                    try:
                        writer.close()
                    except BaseException:
                        pass
                writing_progress.fail(title=writing_title)
                raise
            else:
                for writer in relationship_writers:
                    writer.close()
                writing_progress.complete(
                    retained_pathways,
                    total=retained_pathways,
                    title="Validate written pathway sets",
                )
            active_output_path = mapping_file
            write_text(mapping_file, mapping_text)
            active_output_path = unmapped_file
            with unmapped_file.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["pathway", "gene"])
                for pathway in _iter_gmt(gmt, args.duplicate_gene_policy):
                    _, missing, _ = _pathway_gene_groups(
                        pathway, gene_coordinates, gene_variants,
                    )
                    writer.writerows((pathway.name, gene) for gene in missing)
            unique_variants_written = sum(unique_variant_flags)
            bim_metrics["overlapping_gene_variant_memberships_collapsed"] = (
                overlapping_memberships_collapsed
            )

        with progress.step(7, 8, "Write provenance and checksums"):
            output_manifest = {
                "pathway_mapping": mapping_file.name,
                "pathway_gene_mapping": (
                    pathway_gene_file.name if pathway_gene_file is not None else None
                ),
                "gene_variant_mapping": (
                    gene_variant_file.name if gene_variant_file is not None else None
                ),
                "expanded_mapping": (
                    expanded_file.name if expanded_file is not None else None
                ),
                "unmapped_genes": unmapped_file.name,
                "manifest": args.manifest_name,
                "readme": args.readme_name,
                "checksums": args.checksums_name,
            }
            manifest = {
                "schema_version": "gcta_fastbat_pathway_resource.v2",
                "resource": {
                    "name": args.resource_name,
                    "genome_build": args.genome_build,
                    "method": "fastbat_set",
                    "gcta_option": "--fastBAT-set-list",
                    "file": output_file.name,
                    "format": "GCTA set blocks: set ID, exact BIM IDs, END",
                    "gene_window_kb": args.gene_window_kb,
                    "sha256": sha256(output_file),
                },
                "sources": {
                    "gmt": {"path": str(gmt), "sha256": sha256(gmt)},
                    "gene_coordinates": {
                        "path": str(gene_list),
                        "sha256": sha256(gene_list),
                    },
                    "plink_bim": {
                        "path": str(bim),
                        "sha256": sha256(bim),
                    },
                    "analysis_variants": (
                        {
                            "path": str(analysis_source),
                            "sha256": sha256(analysis_source),
                            "analyzable_unique_variants": len(
                                analyzable_variant_ids
                            ),
                        }
                        if analysis_source is not None else None
                    ),
                },
                "policies": {
                    "allowed_chromosomes": sorted(allowed_chromosomes),
                    "chromosome_label_policy": args.chromosome_label_policy,
                    "duplicate_gene_policy": args.duplicate_gene_policy,
                    "unmapped_gene_policy": args.unmapped_gene_policy,
                    "empty_pathway_policy": args.empty_pathway_policy,
                    "oversized_set_policy": oversized_set_policy,
                    "maximum_set_variants": int(maximum_set_variants),
                    "maximum_set_variants_applied": apply_set_size_policy,
                    "variant_universe": (
                        "gwas_bim_intersection"
                        if analyzable_variant_ids is not None
                        else "plink_bim_reference"
                    ),
                    "audit": {
                        "level": audit_level,
                        "format": audit_format,
                        "compression": audit_compression,
                        "batch_rows": int(audit_batch_rows),
                    },
                    "disk": {
                        "minimum_free_gb": float(minimum_free_disk_gb),
                        "estimation_safety_factor": float(
                            disk_estimation_safety_factor
                        ),
                    },
                    "variant_identifier_policy": (
                        "copy PLINK BIM column 2 verbatim; never infer rsIDs"
                    ),
                },
                "validation": {
                    "input_pathways": pathway_count,
                    "pathways_written": retained_pathways,
                    "pathways_omitted": len(empty) + len(oversized),
                    "pathways_omitted_empty": len(empty),
                    "pathways_omitted_oversized": len(oversized),
                    "input_unique_genes": len(requested_genes),
                    "genes_with_coordinates": len(gene_coordinates),
                    "unmapped_unique_genes": len(missing_genes),
                    "unmapped_gene_entries": unmapped_entry_count,
                    "bim_variants": bim_metrics["bim_variants"],
                    "analyzable_bim_variants_cached": bim_metrics[
                        "analyzable_bim_variants_cached"
                    ],
                    "analyzable_bim_variants_mapped_to_requested_genes": (
                        bim_metrics["bim_variants_mapped_to_requested_genes"]
                    ),
                    "overlapping_gene_variant_memberships_collapsed": (
                        bim_metrics[
                            "overlapping_gene_variant_memberships_collapsed"
                        ]
                    ),
                    "total_pathway_variant_memberships": total_memberships,
                    "total_pathway_gene_variant_memberships": (
                        total_gene_variant_memberships
                    ),
                    "pathway_gene_audit_rows": pathway_gene_rows,
                    "gene_variant_audit_rows": gene_variant_rows,
                    "expanded_audit_rows": expanded_rows,
                    "unique_variants_written": unique_variants_written,
                    "shared_chromosomes": sorted(shared_chromosomes),
                },
                "outputs": output_manifest,
                "disk_preflight": disk_preflight,
                "generation": {
                    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                    "mapping_algorithm": (
                        "chromosome_bisect_compact_variant_index_cache"
                    ),
                    "requested_mapping_workers": bim_metrics[
                        "requested_mapping_workers"
                    ],
                    "mapping_workers": bim_metrics["mapping_workers"],
                    "largest_chromosome_cache_bytes": bim_metrics[
                        "largest_chromosome_cache_bytes"
                    ],
                    "estimated_worker_memory_bytes": bim_metrics[
                        "estimated_worker_memory_bytes"
                    ],
                    "mapping_memory_bytes": bim_metrics["mapping_memory_bytes"],
                    "command": (
                        getattr(args, "generation_command", None)
                        or shlex.join(sys.argv)
                    ),
                    "utility": str(Path(__file__).resolve()),
                },
                "scientific_sources": [OFFICIAL_FASTBAT_DOCUMENTATION],
            }
            manifest_file = staging / args.manifest_name
            readme_file = staging / args.readme_name
            checksum_file = staging / args.checksums_name
            active_output_path = manifest_file
            write_text(
                manifest_file,
                yaml.safe_dump(manifest, sort_keys=False, allow_unicode=True),
            )
            active_output_path = readme_file
            write_text(readme_file, _readme_text(manifest))
            checksummed = [
                output_file, mapping_file, unmapped_file, manifest_file, readme_file,
            ]
            checksummed.extend(
                path for path in (
                    pathway_gene_file, gene_variant_file, expanded_file,
                )
                if path is not None
            )
            active_output_path = checksum_file
            write_text(
                checksum_file,
                "".join(f"{sha256(path)}  {path.name}\n" for path in checksummed),
            )

        with progress.step(8, 8, "Publish validated pathway resource"):
            staging.replace(output_directory)
    except BaseException as exc:
        disk_full = (
            isinstance(exc, OSError) and exc.errno == errno.ENOSPC
        ) or "no space left on device" in str(exc).lower()
        failure_details = None
        if disk_full and staging is not None:
            try:
                available_bytes = shutil.disk_usage(staging).free
            except OSError:
                available_bytes = None
            bytes_written = sum(
                path.stat().st_size
                for path in staging.rglob("*")
                if path.is_file()
            )
            failure_details = (
                "disk space was exhausted during GMT-to-fastBAT preparation: "
                "failed output=%s, staging directory=%s, bytes written=%s, "
                "estimated output bytes=%s, available bytes=%s, audit level=%s, "
                "completed pathways=%s/%s; the incomplete staging resource was "
                "not published and the next run will restart from the safe "
                "preparation boundary"
                % (
                    active_output_path or "unknown",
                    staging,
                    format(bytes_written, ","),
                    format(disk_preflight.get("estimated_output_bytes", 0), ","),
                    (
                        format(available_bytes, ",")
                        if available_bytes is not None else "unknown"
                    ),
                    audit_level,
                    completed_output_pathways,
                    disk_preflight.get("retained_pathways", retained_pathways),
                )
            )
        if staging is not None:
            shutil.rmtree(staging, ignore_errors=True)
        if failure_details is not None:
            raise ResourcePreparationError(failure_details) from exc
        raise
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert GMT gene pathways to GCTA fastBAT SNP-set blocks using "
            "build-matched gene coordinates and exact PLINK BIM variant IDs."
        )
    )
    parser.add_argument("--gmt", required=True, type=Path, metavar="PATH")
    parser.add_argument("--gene-list", required=True, type=Path, metavar="PATH")
    parser.add_argument("--bim", required=True, type=Path, metavar="PATH")
    parser.add_argument("--output-directory", required=True, type=Path, metavar="PATH")
    parser.add_argument("--output-name", required=True, metavar="NAME")
    parser.add_argument("--resource-name", required=True, metavar="NAME")
    parser.add_argument("--genome-build", required=True, metavar="BUILD")
    parser.add_argument("--gene-window-kb", required=True, type=int, metavar="KB")
    parser.add_argument(
        "--allowed-chromosomes", required=True, nargs="+", metavar="CHR",
    )
    parser.add_argument(
        "--chromosome-label-policy",
        required=True,
        choices=("exact", "strip_chr_prefix"),
    )
    parser.add_argument(
        "--duplicate-gene-policy", required=True, choices=("error", "deduplicate"),
    )
    parser.add_argument(
        "--unmapped-gene-policy", required=True, choices=("error", "report"),
    )
    parser.add_argument(
        "--empty-pathway-policy", required=True, choices=("error", "omit"),
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=argparse.SUPPRESS,
        metavar="N",
        help="Maximum chromosome-mapping workers; omitted value comes from YAML.",
    )
    parser.add_argument(
        "--memory-gb",
        type=float,
        default=argparse.SUPPRESS,
        metavar="GB",
        help="Memory budget for chromosome-mapping workers; omitted value comes from YAML.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    from postgwas.config import load_configuration

    configuration = load_configuration()
    conversion = configuration.modules.gcta_gene.set_annotation.conversion
    names = conversion.output_names
    args.pathway_mapping_name = names.pathway_mapping
    args.pathway_gene_mapping_name = names.pathway_gene_mapping
    args.gene_variant_mapping_name = names.gene_variant_mapping
    args.expanded_mapping_name = names.expanded_mapping
    args.unmapped_genes_name = names.unmapped_genes
    args.manifest_name = names.manifest
    args.readme_name = names.readme
    args.checksums_name = names.checksums
    for value in (args.output_name, args.resource_name, args.genome_build):
        if not str(value).strip():
            parser.error("output name, resource name, and genome build must be non-empty")
    try:
        manifest = prepare_resource(
            args,
            mapping_workers=getattr(
                args, "threads", configuration.execution.threads,
            ),
            mapping_memory_gb=getattr(
                args, "memory_gb", configuration.execution.memory_gb,
            ),
            worker_memory_multiplier=(
                conversion.parallelism.worker_memory_multiplier
            ),
            maximum_set_variants=(
                configuration.modules.gcta_gene.set_annotation.maximum_set_variants
            ),
            oversized_set_policy=(
                configuration.modules.gcta_gene.set_annotation.oversized_set_policy
            ),
            audit_level=conversion.audit.level,
            audit_format=conversion.audit.format,
            audit_compression=conversion.audit.compression,
            audit_batch_rows=conversion.audit.batch_rows,
            minimum_free_disk_gb=conversion.disk.minimum_free_gb,
            disk_estimation_safety_factor=(
                conversion.disk.estimation_safety_factor
            ),
        )
    except (OSError, ResourcePreparationError, yaml.YAMLError) as exc:
        parser.error(str(exc))
    print(
        "Created %s (%d pathways, %d unique BIM variants, %s)"
        % (
            Path(args.output_directory).expanduser().resolve() / args.output_name,
            manifest["validation"]["pathways_written"],
            manifest["validation"]["unique_variants_written"],
            manifest["resource"]["genome_build"],
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
