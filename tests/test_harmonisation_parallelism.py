"""CPU, RAM, and Polars-pool contracts for chromosome workers."""

import os
from contextlib import contextmanager
from pathlib import Path
import subprocess
import sys
from multiprocessing import get_context
from unittest.mock import patch

import pytest

from postgwas.modules.harmonisation import service as service_module
from postgwas.modules.harmonisation.policies import default_policies
from postgwas.modules.harmonisation.service import (
    _POLARS_THREAD_ENVIRONMENT_VARIABLE,
    _bounded_chromosome_executor,
    _derive_parallelism,
    _run_one_round,
)


class _RecordingExecutor:
    created = []

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self.created.append(self)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback_value):
        return False


def _worker_polars_size():
    import polars as pl

    return (
        os.environ.get(_POLARS_THREAD_ENVIRONMENT_VARIABLE),
        pl.thread_pool_size(),
    )


def _parallelism_policies(*, cpu, memory):
    return default_policies().with_overrides({
        "execution.total_cpu_budget": cpu,
        "execution.memory_budget_gb": memory,
    })


def test_cpu_budget_preserves_per_chromosome_threads_and_reduces_workers():
    workers, threads = _derive_parallelism(
        _parallelism_policies(cpu=16, memory=100),
        requested_workers=16,
        n_chromosomes=22,
    )

    assert (workers, threads) == (3, 5)


def test_memory_budget_can_reduce_cpu_derived_worker_count():
    workers, threads = _derive_parallelism(
        _parallelism_policies(cpu=16, memory=48),
        requested_workers=16,
        n_chromosomes=22,
    )

    assert (workers, threads) == (2, 5)


def test_cpu_budget_smaller_than_one_worker_reduces_that_worker():
    workers, threads = _derive_parallelism(
        _parallelism_policies(cpu=4, memory=48),
        requested_workers=16,
        n_chromosomes=22,
    )

    assert (workers, threads) == (1, 4)


def test_executor_sets_pre_spawn_limit_and_restores_parent_environment(
    monkeypatch,
):
    _RecordingExecutor.created.clear()
    monkeypatch.setenv(_POLARS_THREAD_ENVIRONMENT_VARIABLE, "11")

    with patch(
        "postgwas.core.polars_runtime.ProcessPoolExecutor",
        _RecordingExecutor,
    ):
        with _bounded_chromosome_executor(3, 5, "spawn-context") as executor:
            assert os.environ[_POLARS_THREAD_ENVIRONMENT_VARIABLE] == "5"
            assert executor.kwargs["max_workers"] == 3
            assert executor.kwargs["mp_context"] == "spawn-context"
            assert executor.kwargs["initargs"] == (5,)

    assert os.environ[_POLARS_THREAD_ENVIRONMENT_VARIABLE] == "11"


@pytest.mark.parametrize(
    ("workers", "threads", "message"),
    [
        (0, 1, "max_workers must be at least 1"),
        (1, 0, "threads_per_worker must be at least 1"),
    ],
)
def test_bounded_polars_executor_rejects_invalid_resource_limits(
    workers, threads, message,
):
    with pytest.raises(ValueError, match=message):
        with _bounded_chromosome_executor(workers, threads, get_context("spawn")):
            pass


def test_fresh_interpreter_uses_inherited_polars_thread_limit():
    project = Path(__file__).resolve().parents[1]
    environment = dict(os.environ)
    environment[_POLARS_THREAD_ENVIRONMENT_VARIABLE] = "2"
    environment["PYTHONPATH"] = str(project / "src")
    completed = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "from postgwas.modules.harmonisation.service import "
                "_initialise_chromosome_worker; "
                "_initialise_chromosome_worker(2); "
                "import polars as pl; print(pl.thread_pool_size())"
            ),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
        cwd=project,
    )

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "2"


def test_real_spawn_pool_bounds_polars_and_restores_unset_parent_value(
    monkeypatch,
):
    monkeypatch.delenv(_POLARS_THREAD_ENVIRONMENT_VARIABLE, raising=False)

    try:
        with _bounded_chromosome_executor(
            1, 2, get_context("spawn")
        ) as executor:
            configured, actual = executor.submit(_worker_polars_size).result(
                timeout=30
            )
    except PermissionError as exc:
        pytest.skip("platform forbids multiprocessing semaphores: %s" % exc)

    assert (configured, actual) == ("2", 2)
    assert _POLARS_THREAD_ENVIRONMENT_VARIABLE not in os.environ


def test_parent_does_not_convert_worker_memory_error_into_retryable_result(
    tmp_path, monkeypatch,
):
    class Future:
        def __init__(self, error=None):
            self.error = error
            self.cancelled = False

        def result(self):
            if self.error is not None:
                raise self.error
            return "2", {}, "", "ok"

        def cancel(self):
            self.cancelled = True
            return True

    memory_future = Future(MemoryError("allocation failed"))
    pending_future = Future()

    class Executor:
        def __init__(self):
            self.futures = [memory_future, pending_future]

        def submit(self, *_args, **_kwargs):
            return self.futures.pop(0)

    @contextmanager
    def fake_executor(*_args, **_kwargs):
        yield Executor()

    class Logger:
        def __init__(self):
            self.errors = []

        def error(self, message):
            self.errors.append(message)

    logger = Logger()
    monkeypatch.setattr(
        service_module, "_bounded_chromosome_executor", fake_executor,
    )
    monkeypatch.setattr(
        service_module, "as_completed", lambda futures: list(futures),
    )
    monkeypatch.setattr(
        service_module, "print_chromosome_summary", lambda **_kwargs: None,
    )
    monkeypatch.setattr(
        service_module, "_log_tail", lambda *_args, **_kwargs: "worker log",
    )

    with pytest.raises(MemoryError, match="without retrying"):
        _run_one_round(
            pending=["1", "2"],
            chr_file_by_chrom={"1": "chr1.parquet", "2": "chr2.parquet"},
            resource_maps={"1": {}, "2": {}},
            round_no=1,
            attempts={"1": 0, "2": 0},
            worker_kwargs={"threads": 1, "sample_column_dict": {}},
            max_workers=2,
            sample_id="study",
            output_dir=tmp_path,
            output_layout={
                "chromosome_log": "logs/{dataset_id}_chr{chromosome}.log",
            },
            screen_order="completion",
            logger=logger,
        )

    assert pending_future.cancelled is True
    assert any("without retrying" in message for message in logger.errors)
