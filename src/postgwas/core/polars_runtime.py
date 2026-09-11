"""Process-scoped Polars thread-pool control shared by PostGWAS modules.

Polars fixes its global thread pool when it is first imported.  A configured
thread limit therefore has to be present in the environment inherited by a new
process; changing the environment in an already-imported process is too late.
"""

from __future__ import annotations

import os
import traceback
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from multiprocessing import get_context
from threading import Lock
from typing import Any, Callable, Iterator, Mapping


POLARS_THREAD_ENVIRONMENT_VARIABLE = "POLARS_MAX_THREADS"
_POLARS_ENVIRONMENT_LOCK = Lock()


def current_polars_thread_runtime(expected_threads: int) -> dict[str, Any]:
    """Validate and describe the effective Polars pool in the current process."""
    expected = int(expected_threads)
    if expected < 1:
        raise ValueError("expected_threads must be at least 1")
    configured = os.environ.get(POLARS_THREAD_ENVIRONMENT_VARIABLE)
    if configured != str(expected):
        raise RuntimeError(
            "Polars worker did not inherit %s=%d before process start "
            "(received %r)."
            % (POLARS_THREAD_ENVIRONMENT_VARIABLE, expected, configured)
        )

    import polars as pl

    actual = int(pl.thread_pool_size())
    if actual != expected:
        raise RuntimeError(
            "Polars worker pool has %d threads, expected %d from %s."
            % (actual, expected, POLARS_THREAD_ENVIRONMENT_VARIABLE)
        )
    return {
        "requested_threads": expected,
        "polars_thread_pool_size": actual,
        "thread_budget_enforced": True,
        "thread_environment_variable": POLARS_THREAD_ENVIRONMENT_VARIABLE,
    }


def initialise_polars_worker(expected_threads: int) -> None:
    """Fail worker startup unless the pre-import Polars limit took effect."""
    current_polars_thread_runtime(expected_threads)


@contextmanager
def _bounded_polars_environment(threads: int) -> Iterator[None]:
    """Scope the pre-spawn environment without racing another local caller."""
    if threads < 1:
        raise ValueError("threads_per_worker must be at least 1")
    with _POLARS_ENVIRONMENT_LOCK:
        previous = os.environ.get(POLARS_THREAD_ENVIRONMENT_VARIABLE)
        os.environ[POLARS_THREAD_ENVIRONMENT_VARIABLE] = str(threads)
        try:
            yield
        finally:
            if previous is None:
                os.environ.pop(POLARS_THREAD_ENVIRONMENT_VARIABLE, None)
            else:
                os.environ[POLARS_THREAD_ENVIRONMENT_VARIABLE] = previous


@contextmanager
def bounded_polars_executor(
    max_workers: int,
    threads_per_worker: int,
    mp_context=None,
) -> Iterator[ProcessPoolExecutor]:
    """Yield a process pool whose workers have bounded Polars thread pools."""
    workers = int(max_workers)
    threads = int(threads_per_worker)
    if workers < 1:
        raise ValueError("max_workers must be at least 1")
    if threads < 1:
        raise ValueError("threads_per_worker must be at least 1")
    context = mp_context or get_context("spawn")
    with _bounded_polars_environment(threads):
        with ProcessPoolExecutor(
            max_workers=workers,
            mp_context=context,
            initializer=initialise_polars_worker,
            initargs=(threads,),
        ) as executor:
            yield executor


def _run_polars_call(
    connection,
    expected_threads: int,
    function: Callable[..., Any],
    args: tuple[Any, ...],
    call_kwargs: dict[str, Any],
) -> None:
    """Child entry point for one pipe-returned, semaphore-free operation."""
    try:
        initialise_polars_worker(expected_threads)
        connection.send(("result", function(*args, **call_kwargs), None))
    except BaseException as exc:
        remote_traceback = traceback.format_exc()
        try:
            connection.send(("exception", exc, remote_traceback))
        except BaseException:
            connection.send((
                "error_text",
                "%s: %s" % (type(exc).__name__, exc),
                remote_traceback,
            ))
    finally:
        connection.close()


def run_in_bounded_polars_process(
    function: Callable[..., Any],
    *args: Any,
    threads: int,
    call_kwargs: Mapping[str, Any] | None = None,
) -> Any:
    """Run one picklable Polars operation in a fresh bounded worker process."""
    worker_threads = int(threads)
    if worker_threads < 1:
        raise ValueError("threads must be at least 1")
    context = get_context("spawn")
    receiver, sender = context.Pipe(duplex=False)
    process = context.Process(
        target=_run_polars_call,
        args=(
            sender,
            worker_threads,
            function,
            tuple(args),
            dict(call_kwargs or {}),
        ),
    )
    status = None
    payload = None
    remote_traceback = None
    with _bounded_polars_environment(worker_threads):
        try:
            process.start()
        except BaseException:
            sender.close()
            receiver.close()
            raise
        sender.close()
        try:
            while not receiver.poll(0.1):
                if not process.is_alive():
                    break
            if receiver.poll():
                status, payload, remote_traceback = receiver.recv()
        except BaseException:
            if process.is_alive():
                process.terminate()
            process.join()
            raise
        finally:
            receiver.close()
        process.join()

    if status == "result" and process.exitcode == 0:
        return payload
    if status == "exception":
        raise payload
    if status == "error_text":
        raise RuntimeError(
            "Bounded Polars worker failed: %s\n%s"
            % (payload, remote_traceback or "")
        )
    raise RuntimeError(
        "Bounded Polars worker exited with code %s without returning a result."
        % process.exitcode
    )
