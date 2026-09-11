"""Resource admission and scientific-value invariance across worker counts."""

from concurrent.futures import Future
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import polars as pl
import pytest

from postgwas.core.execution.scheduling import ResourceQueue
from postgwas.modules.harmonisation import service
from postgwas.modules.harmonisation.policies import PolicyError, default_policies
from postgwas.modules.harmonisation.scheduling import plan_chromosome_schedule
from postgwas.modules.harmonisation.z_score import derive_z_score_from_effect_and_standard_error


def test_queue_mixes_sizes_refills_and_never_overcommits():
    queue = ResourceQueue({"large": 20, "medium": 18, "small": 7, "tiny": 5}, 30, 3)
    assert queue.admit() == "large"
    assert queue.admit() == "tiny"
    assert queue.admit() is None
    queue.complete("large")
    assert queue.admit() == "medium"
    assert queue.admit() == "small"
    assert queue.reserved == 30
    assert queue.admit() is None


@pytest.mark.parametrize("cost,budget", [(31, 30), (0, 30), (-1, 30), (float("nan"), 30), (1, 0), (1, float("inf"))])
def test_queue_rejects_invalid_or_unrunnable_work_before_submission(cost, budget):
    with pytest.raises(ValueError):
        ResourceQueue({"chr1": cost}, budget, 2)


@pytest.fixture
def workloads(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "postgwas.modules.harmonisation.scheduling.psutil.Process",
        lambda: SimpleNamespace(memory_info=lambda: SimpleNamespace(rss=1024**3)),
    )
    files, counts, resources = {}, {}, {}
    for i, count in enumerate([1000, 700, 500, 300, 200, 100, 50], 1):
        chrom = str(i)
        path = tmp_path / (chrom + ".parquet")
        pl.DataFrame({"BETA": [0.0, 0.1, -0.2, None] * count, "SE": [0.1, 0.2, 0.3, 0.4] * count}).write_parquet(path)
        files[chrom], counts[chrom], resources[chrom] = str(path), count * 4, {}
    return files, counts, resources


def plan(workloads, tmp_path, cpu=11, memory=57, **overrides):
    pol = default_policies().with_overrides({
        "execution.total_cpu_budget": cpu,
        "execution.memory_budget_gb": memory,
        "rejects.enabled": False,
        **overrides,
    })
    return plan_chromosome_schedule(
        pol, *workloads, requested_workers=cpu, output_dir=tmp_path,
        output_layout={}, sample_id="study",
    )


@pytest.mark.parametrize("cpu,memory,workers,threads", [(11, 57, 5, 2), (4, 20, 2, 2), (1, 10, 1, 1), (32, 128, 7, 4)])
def test_adapts_to_resolved_system_budgets(workloads, tmp_path, cpu, memory, workers, threads):
    result = plan(workloads, tmp_path, cpu, memory)
    assert (result["workers"], result["threads_per_chromosome"]) == (workers, threads)
    assert workers * threads <= cpu
    assert result["worker_memory_budget_gb"] < memory
    costs = result["estimated_memory_gb_by_chromosome"]
    assert costs["1"] > costs["7"]


def test_reference_size_and_duplicate_roles_are_accounted_once(workloads, tmp_path):
    base = plan(workloads, tmp_path)
    reference = tmp_path / "reference.parquet"
    pl.DataFrame({"AF": [0.1, 0.2] * 1000}).write_parquet(reference)
    workloads[2]["1"].update(user_eaf_file=str(reference), user_info_file=str(reference))
    result = plan(workloads, tmp_path)
    assert result["estimated_memory_gb_by_chromosome"]["1"] > base["estimated_memory_gb_by_chromosome"]["1"]
    del workloads[2]["1"]["user_info_file"]
    assert plan(workloads, tmp_path) == result


def test_oversized_chromosome_and_missing_reference_stop_planning(workloads, tmp_path):
    with pytest.raises(ValueError, match="exceeding.*Increase --memory-gb"):
        plan(workloads, tmp_path, memory=3)
    workloads[2]["1"]["user_info_file"] = str(tmp_path / "absent.tsv")
    with pytest.raises(FileNotFoundError):
        plan(workloads, tmp_path)


@pytest.mark.parametrize("key,value", [
    ("scheduling_mode", "unknown"), ("min_threads_per_chromosome", 6),
    ("memory_headroom_fraction", 1), ("memory_safety_factor", 0.5),
    ("reference_memory_multiplier", float("nan")), ("worker_memory_floor_gb", float("inf")),
])
def test_invalid_scheduling_policy_is_rejected(key, value):
    with pytest.raises(PolicyError):
        default_policies().with_overrides({"execution." + key: value})


def test_missing_library_budgets_have_actionable_error(workloads, tmp_path):
    with pytest.raises(ValueError, match="Provide --threads and --memory-gb"):
        plan_chromosome_schedule(default_policies(), *workloads, requested_workers=2,
                                 output_dir=tmp_path, output_layout={}, sample_id="study")


def _scientific_worker(**kwargs):
    chrom = kwargs["chromosome"]
    frame = pl.read_parquet(kwargs["chr_file"])
    result, qc, mapping = derive_z_score_from_effect_and_standard_error(
        chrom, frame, {"beta_col": "BETA", "se_col": "SE", "imp_z_col": None},
        policies=default_policies(),
    )
    return chrom, {"values": result.to_dict(as_series=False)}, "", "ok"


def test_round_refills_slots_and_preserves_values(workloads, tmp_path, monkeypatch):
    submitted = []
    monkeypatch.setattr(service, "process_one_chromosome", _scientific_worker)
    monkeypatch.setattr(service, "print_chromosome_summary", lambda **kwargs: None)

    class Executor:
        def submit(self, fn, **kwargs):
            submitted.append(kwargs["chromosome"])
            future = Future()
            future.set_result(fn(**kwargs))
            return future

    @contextmanager
    def executor(*args):
        yield Executor()

    monkeypatch.setattr(service, "_bounded_chromosome_executor", executor)
    monkeypatch.setattr(service, "as_completed", lambda futures: iter(futures))
    files, counts, resources = workloads
    arguments = dict(
        pending=list(files), chr_file_by_chrom=files, resource_maps=resources,
        round_no=1, attempts={}, worker_kwargs={"threads": 2, "sample_column_dict": {}},
        sample_id="study", output_dir=tmp_path, output_layout={},
        screen_order="completion", logger=None,
    )
    serial = service._run_one_round(**arguments, max_workers=1)
    submitted.clear()
    adaptive = plan(workloads, tmp_path)
    parallel = service._run_one_round(**arguments, max_workers=5, schedule=adaptive)
    assert serial == parallel
    assert submitted[:2] == ["1", "7"]
    assert len(submitted) == len(set(submitted)) == len(files)


@pytest.mark.parametrize("during_submit", [False, True])
def test_broken_pool_accounts_for_unstarted_chromosomes(tmp_path, monkeypatch, during_submit):
    from concurrent.futures.process import BrokenProcessPool

    class Executor:
        def submit(self, fn, **kwargs):
            if during_submit:
                raise BrokenProcessPool("worker died")
            future = Future()
            future.set_exception(BrokenProcessPool("worker died"))
            return future

    @contextmanager
    def executor(*args):
        yield Executor()

    monkeypatch.setattr(service, "_bounded_chromosome_executor", executor)
    monkeypatch.setattr(service, "print_chromosome_summary", lambda **kwargs: None)
    monkeypatch.setattr(service, "_log_tail", lambda *args: "worker died")
    result = service._run_one_round(
        pending=["1", "2", "3"], chr_file_by_chrom={c: c for c in ["1", "2", "3"]},
        resource_maps={c: {} for c in ["1", "2", "3"]}, round_no=1, attempts={},
        worker_kwargs={"threads": 1, "sample_column_dict": {}}, max_workers=1,
        sample_id="study", output_dir=tmp_path,
        output_layout={"chromosome_log": "{chromosome}.log"},
        screen_order="completion", logger=None,
    )
    assert set(result) == {"1", "2", "3"}
    assert all(value["status"] == "failed" for value in result.values())


def test_source_snapshot_metadata_is_included(workloads, tmp_path):
    files, counts, resources = workloads
    for chromosome, file in files.items():
        frame = pl.read_parquet(file).with_columns(pl.lit("original").alias("RAW"))
        frame.write_parquet(tmp_path / (chromosome + "_source.parquet"))
    pol = default_policies().with_overrides({
        "execution.total_cpu_budget": 11, "execution.memory_budget_gb": 57,
    })
    with_sources = plan_chromosome_schedule(
        pol, files, counts, resources, requested_workers=11, output_dir=tmp_path,
        output_layout={"chromosome_source_snapshot": "{chromosome}_source.parquet"}, sample_id="study",
    )
    without = plan(workloads, tmp_path)
    assert with_sources["workload_evidence"]["1"]["study_encoded_bytes"] > without["workload_evidence"]["1"]["study_encoded_bytes"]


def test_real_spawn_produces_identical_scientific_values(workloads):
    from multiprocessing import get_context

    files, _, _ = workloads
    kwargs = [{"chromosome": chrom, "chr_file": file} for chrom, file in list(files.items())[:2]]
    expected = [_scientific_worker(**item) for item in kwargs]
    with service._bounded_chromosome_executor(2, 2, get_context("spawn")) as executor:
        observed = [future.result(timeout=60) for future in [executor.submit(_scientific_worker, **item) for item in kwargs]]
    assert observed == expected
