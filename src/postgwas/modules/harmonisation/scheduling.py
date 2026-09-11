"""Chromosome workload estimates; no study values are read or transformed.

Parquet uncompressed byte counts describe encoded data, not peak RAM. YAML
multipliers cover decoded strings, joins, immutable snapshots and allocation
overhead. Non-Parquet references use a deliberately explicit disk expansion
estimate. These estimates must be calibrated against representative runs.
"""

from __future__ import annotations

import math
from pathlib import Path

import psutil
import pyarrow.parquet as pq

from postgwas.core.execution.scheduling import ResourceQueue
from postgwas.core.paths import configured_output_path


def plan_chromosome_schedule(
    policies, chromosome_files, partition_rows, resource_maps, *,
    requested_workers, output_dir, output_layout, sample_id,
):
    """Resolve a reproducible estimate and CPU allocation before any submission."""
    if not chromosome_files:
        raise ValueError("Adaptive chromosome scheduling requires at least one validated partition.")
    cpu = policies.get("execution.total_cpu_budget")
    memory = policies.get("execution.memory_budget_gb")
    if cpu is None or memory is None:
        raise ValueError(
            "Adaptive chromosome scheduling requires resolved CPU and memory "
            "budgets. Provide --threads and --memory-gb, or set "
            "execution.total_cpu_budget and execution.memory_budget_gb for "
            "a direct library call."
        )
    cpu, memory = int(cpu), float(memory)
    # A GiB is exactly 2**30 bytes. Match existing PostGWAS memory units.
    gib = float(1024**3)
    parent = psutil.Process().memory_info().rss / gib
    headroom = max(parent, memory * policies.get("execution.memory_headroom_fraction"))
    available = memory - headroom
    floor = float(policies.get("execution.worker_memory_floor_gb"))
    margin = float(policies.get("execution.memory_safety_factor"))
    minimum_threads = int(policies.get("execution.min_threads_per_chromosome"))
    maximum_threads = int(policies.get("execution.threads_per_chromosome"))
    if minimum_threads > maximum_threads:
        raise ValueError("execution.min_threads_per_chromosome must not exceed execution.threads_per_chromosome.")
    numeric = [memory, floor, margin, available]
    if cpu < 1 or requested_workers < 1 or any(not math.isfinite(v) or v <= 0 for v in numeric):
        raise ValueError("Adaptive chromosome scheduling has no usable CPU/memory budget after parent headroom.")
    # The lightest jobs bound possible concurrency. Actual admission below also
    # checks the sum, so a large chromosome can reduce simultaneous jobs.
    workers = min(len(chromosome_files), requested_workers, max(1, cpu // minimum_threads))
    threads = min(maximum_threads, max(1, cpu // workers))
    estimates = {}
    evidence = {}
    metadata_cache = {}

    def file_bytes(path):
        path = Path(path).resolve()
        if path not in metadata_cache:
            size = path.stat().st_size
            if path.suffix == ".parquet":
                metadata = pq.read_metadata(path)
                size = sum(metadata.row_group(i).total_byte_size for i in range(metadata.num_row_groups))
            metadata_cache[path] = size
        return metadata_cache[path]

    for chromosome, file in chromosome_files.items():
        rows = int(partition_rows[chromosome])
        if rows < 1:
            raise ValueError("Chromosome %s has no validated partition rows." % chromosome)
        paths = {Path(file).resolve()}
        if policies.get("rejects.enabled"):
            paths.add(configured_output_path(
                output_dir, output_layout["chromosome_source_snapshot"],
                dataset_id=sample_id, chromosome=chromosome,
            ).resolve())
        study_bytes = math.fsum(file_bytes(path) for path in paths)
        references = set()
        for key in policies.get("execution.memory_reference_keys"):
            value = resource_maps[chromosome].get(key)
            if value:
                reference = Path(value).resolve()
                # Preflight permits an absent default panel when neither
                # strand nor MAF validation uses it; the map still names it.
                if key == "default_eaf_file" and not reference.exists():
                    continue
                references.add(reference)
        reference_bytes = (
            math.fsum(file_bytes(path) for path in references - paths)
            * policies.get("execution.reference_memory_multiplier")
        )
        working_bytes = max(
            rows * policies.get("execution.memory_bytes_per_variant"),
            study_bytes * policies.get("execution.study_memory_multiplier"),
        )
        # Include configured sort buffers as well as Python/Polars working data.
        sort_gb = threads * policies.get("vcf.sort_memory_mb_per_thread") / 1024.0
        estimate = margin * (floor + (working_bytes + reference_bytes) / gib + sort_gb)
        estimates[chromosome] = estimate
        evidence[chromosome] = {
            "variants": rows, "study_encoded_bytes": study_bytes,
            "reference_expanded_bytes_estimate": reference_bytes,
            "estimated_memory_gb": estimate,
        }
    # Check *every* chromosome now, including the largest, before starting any.
    ResourceQueue(estimates, available, workers)
    fitting = 0
    used = 0.0
    for cost in sorted(estimates.values()):
        if used + cost > available:
            break
        fitting += 1
        used += cost
    workers = min(workers, fitting)
    return {
        "mode": "adaptive", "workers": workers,
        "threads_per_chromosome": threads,
        "total_cpu_budget": cpu, "memory_budget_gb": memory,
        "parent_rss_gb": parent, "memory_headroom_gb": headroom,
        "worker_memory_budget_gb": available,
        "estimated_memory_gb_by_chromosome": estimates,
        "workload_evidence": evidence,
    }
