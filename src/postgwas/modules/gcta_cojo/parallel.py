"""Memory-bounded chromosome execution for genome-wide COJO stepwise selection."""

from __future__ import annotations

import math
import re
import shutil
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

from postgwas.config.models.modules.gcta_cojo import (
    GctaCojoParallelOutputContract,
)
from postgwas.core.paths import configured_output_path, require_nonempty_file
from postgwas.core.ui import (
    MeasuredProgress,
    print_screen_block,
    screen_field,
    screen_line,
)
from postgwas.modules.gcta_cojo.adapters import (
    build_cojo_command,
    run_cojo_command,
)
from postgwas.modules.gcta_cojo.errors import GctaCojoError
from postgwas.modules.gcta_cojo.results import (
    normalize_cojo_results,
    write_empty_cojo_results,
)


_GIB = float(1024**3)


@dataclass(frozen=True)
class ChromosomeWorkerPlan:
    """Resolved worker limits after CPU and measured-memory arbitration."""

    workers: int
    cpu_limit: int
    memory_limit: int
    configured_limit: int
    threads_per_worker: int
    memory_budget_gb: float
    parent_memory_gb: float
    available_worker_memory_gb: float
    estimated_memory_gb_per_worker: float


@dataclass(frozen=True)
class _AttemptPaths:
    prefix: Path
    primary: Path
    ld_matrix: Path
    conditional: Path
    condition_snps: Path
    log: Path


@dataclass(frozen=True)
class _AttemptOutcome:
    chromosome: int | None
    attempt: int
    command: tuple[str, ...]
    paths: _AttemptPaths
    peak_rss_bytes: int
    memory_samples: int
    process_tree_samples: int
    root_only_samples: int
    no_signals: bool
    error: str | None
    error_type: str | None
    phase: str


@dataclass(frozen=True)
class ParallelSlctResult:
    """Validated combined outputs and execution metrics for the caller."""

    metrics: dict
    no_signals: bool
    dry_run: bool


def _parent_rss_gb() -> float:
    try:
        import psutil

        return float(psutil.Process().memory_info().rss) / _GIB
    except (ImportError, OSError, RuntimeError, ValueError):
        return 0.0


def plan_chromosome_workers(
    chromosome_count: int,
    *,
    total_threads: int,
    threads_per_worker: int,
    memory_budget_gb: float,
    parent_memory_gb: float,
    estimated_memory_gb_per_worker: float,
    configured_max_workers,
) -> ChromosomeWorkerPlan:
    """Choose the largest concurrency that fits every resolved resource cap."""
    count = int(chromosome_count)
    total = int(total_threads)
    worker_threads = int(threads_per_worker)
    budget = float(memory_budget_gb)
    parent = max(0.0, float(parent_memory_gb))
    per_worker = float(estimated_memory_gb_per_worker)
    if count < 1:
        raise GctaCojoError("Chromosome worker planning requires at least one job.")
    if total < 1 or worker_threads < 1:
        raise GctaCojoError("Chromosome worker CPU limits must be positive.")
    if not math.isfinite(budget) or budget <= 0:
        raise GctaCojoError("Resolved execution.memory_gb must be positive.")
    if not math.isfinite(per_worker) or per_worker <= 0:
        raise GctaCojoError("Estimated GCTA memory per worker must be positive.")

    cpu_limit = total // worker_threads
    if cpu_limit < 1:
        raise GctaCojoError(
            "GCTA chromosome execution requests %d thread(s) per worker but "
            "execution.threads provides only %d. Reduce "
            "modules.gcta_cojo.chromosome_execution.threads_per_worker or "
            "increase execution.threads." % (worker_threads, total)
        )
    available = budget - parent
    memory_limit = int(available // per_worker) if available > 0 else 0
    if memory_limit < 1:
        raise GctaCojoError(
            "The resolved %.3f GB memory budget cannot safely start another "
            "GCTA chromosome worker: the PostGWAS parent uses %.3f GB and the "
            "measured, safety-adjusted worker requirement is %.3f GB. Increase "
            "--memory-gb or lower the configured chromosome-worker memory floor."
            % (budget, parent, per_worker)
        )
    configured_limit = (
        count
        if configured_max_workers == "auto"
        else int(configured_max_workers)
    )
    workers = min(count, cpu_limit, memory_limit, configured_limit)
    return ChromosomeWorkerPlan(
        workers=workers,
        cpu_limit=cpu_limit,
        memory_limit=memory_limit,
        configured_limit=configured_limit,
        threads_per_worker=worker_threads,
        memory_budget_gb=budget,
        parent_memory_gb=parent,
        available_worker_memory_gb=available,
        estimated_memory_gb_per_worker=per_worker,
    )


def numeric_chromosome_workloads(
    counts: Mapping[str, int],
) -> dict[int, int] | None:
    """Return numeric chromosome workloads, or None when safe fan-out is impossible."""
    workloads: dict[int, int] = {}
    for raw, raw_count in counts.items():
        try:
            chromosome = int(str(raw).strip())
        except ValueError:
            return None
        if chromosome < 1 or str(chromosome) != str(raw).strip():
            return None
        count = int(raw_count)
        if count < 1 or chromosome in workloads:
            return None
        workloads[chromosome] = count
    return workloads or None


def chromosome_parallel_slct_requested(module) -> bool:
    """Return whether resolved policy requests genome-wide slct fan-out."""
    return bool(
        module.chromosome_execution.enabled
        and module.mode == "slct"
        and module.analysis.chromosome is None
    )


def should_parallelize_slct(module, chromosome_counts: Mapping[str, int]) -> bool:
    """Gate fan-out without changing any non-slct or explicit-chromosome run."""
    workloads = numeric_chromosome_workloads(chromosome_counts)
    return bool(
        chromosome_parallel_slct_requested(module)
        and workloads is not None
        and len(workloads) > 1
    )


def _remap_prefix_path(
    combined_prefix: Path,
    combined_path: Path,
    worker_prefix: Path,
) -> Path:
    if not combined_path.name.startswith(combined_prefix.name):
        raise GctaCojoError(
            "Configured GCTA output %s does not share prefix %s."
            % (combined_path, combined_prefix)
        )
    suffix = combined_path.name[len(combined_prefix.name):]
    return worker_prefix.with_name(worker_prefix.name + suffix)


def _attempt_paths(
    prefix: Path,
    *,
    combined_prefix: Path,
    combined_primary: Path,
    combined_ld: Path,
    combined_conditional: Path,
    combined_condition_snps: Path,
    combined_log: Path,
) -> _AttemptPaths:
    return _AttemptPaths(
        prefix=prefix,
        primary=_remap_prefix_path(combined_prefix, combined_primary, prefix),
        ld_matrix=_remap_prefix_path(combined_prefix, combined_ld, prefix),
        conditional=_remap_prefix_path(
            combined_prefix, combined_conditional, prefix,
        ),
        condition_snps=_remap_prefix_path(
            combined_prefix, combined_condition_snps, prefix,
        ),
        log=_remap_prefix_path(combined_prefix, combined_log, prefix),
    )


def _log_reports_no_signals(path: Path, module) -> bool:
    if not path.is_file():
        return False
    pattern = re.compile(module.log_parsing.no_signals_pattern)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        return any(pattern.search(line) for line in handle)


def _run_slct_attempt(
    chromosome: int,
    attempt: int,
    command: Sequence[str],
    paths: _AttemptPaths,
    configuration,
    module,
) -> _AttemptOutcome:
    paths.prefix.parent.mkdir(parents=True, exist_ok=True)
    resources: dict[str, int] = {}
    error = None
    error_type = None
    no_signals = False
    try:
        run_cojo_command(
            command,
            configuration,
            None,
            dry_run=False,
            purpose="GCTA-COJO slct chromosome %d" % chromosome,
            resource_metrics=resources,
            resource_poll_seconds=(
                module.chromosome_execution.memory_poll_interval_seconds
            ),
        )
        require_nonempty_file(
            paths.log,
            "GCTA chromosome %d log" % chromosome,
            error_type=GctaCojoError,
        )
        no_signals = (
            not paths.primary.is_file()
            and _log_reports_no_signals(paths.log, module)
        )
        if not no_signals:
            for path, label in (
                (paths.primary, "primary result"),
                (paths.ld_matrix, "LD matrix"),
            ):
                require_nonempty_file(
                    path,
                    "GCTA chromosome %d %s" % (chromosome, label),
                    error_type=GctaCojoError,
                )
    except Exception as exc:
        error = str(exc)
        error_type = type(exc).__name__
    return _AttemptOutcome(
        chromosome=chromosome,
        attempt=attempt,
        command=tuple(str(value) for value in command),
        paths=paths,
        peak_rss_bytes=int(resources.get("peak_rss_bytes", 0)),
        memory_samples=int(resources.get("memory_samples", 0)),
        process_tree_samples=int(resources.get("process_tree_samples", 0)),
        root_only_samples=int(resources.get("root_only_samples", 0)),
        no_signals=no_signals,
        error=error,
        error_type=error_type,
        phase="selection",
    )


def _run_conditional_attempt(
    attempt: int,
    command: Sequence[str],
    paths: _AttemptPaths,
    configuration,
    module,
    selected_identifiers: set[str],
) -> _AttemptOutcome:
    paths.prefix.parent.mkdir(parents=True, exist_ok=True)
    resources: dict[str, int] = {}
    error = None
    error_type = None
    try:
        run_cojo_command(
            command,
            configuration,
            None,
            dry_run=False,
            purpose="GCTA-COJO genome-wide conditional output reconstruction",
            resource_metrics=resources,
            resource_poll_seconds=(
                module.chromosome_execution.memory_poll_interval_seconds
            ),
        )
        for path, label in (
            (paths.log, "conditional reconstruction log"),
            (paths.conditional, "conditional reconstruction result"),
            (paths.condition_snps, "conditional selected-SNP result"),
        ):
            require_nonempty_file(path, label, error_type=GctaCojoError)
        observed = _read_identifier_column(
            paths.condition_snps,
            module.results.schemas["cond"].identifier_column,
            module.input_validation.table_delimiter_pattern,
        )
        if observed != selected_identifiers:
            raise GctaCojoError(
                "The genome-wide conditional reconstruction retained %d of %d "
                "selected SNPs and returned %d unexpected SNPs."
                % (
                    len(observed & selected_identifiers),
                    len(selected_identifiers),
                    len(observed - selected_identifiers),
                )
            )
    except Exception as exc:
        error = str(exc)
        error_type = type(exc).__name__
    return _AttemptOutcome(
        chromosome=None,
        attempt=attempt,
        command=tuple(str(value) for value in command),
        paths=paths,
        peak_rss_bytes=int(resources.get("peak_rss_bytes", 0)),
        memory_samples=int(resources.get("memory_samples", 0)),
        process_tree_samples=int(resources.get("process_tree_samples", 0)),
        root_only_samples=int(resources.get("root_only_samples", 0)),
        no_signals=False,
        error=error,
        error_type=error_type,
        phase="conditional_reconstruction",
    )


def _read_identifier_column(
    path: Path,
    identifier_column: str,
    delimiter_pattern: str,
) -> set[str]:
    pattern = re.compile(delimiter_pattern)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        header = handle.readline()
        columns = pattern.split(header.strip()) if header else []
        if identifier_column not in columns:
            raise GctaCojoError(
                "GCTA output %s lacks identifier column %s."
                % (path, identifier_column)
            )
        index = columns.index(identifier_column)
        identifiers = set()
        for line_number, raw in enumerate(handle, 2):
            if not raw.strip():
                continue
            fields = pattern.split(raw.strip())
            if len(fields) != len(columns):
                raise GctaCojoError(
                    "GCTA output %s line %d has %d fields; expected %d."
                    % (path, line_number, len(fields), len(columns))
                )
            identifier = fields[index]
            if identifier in identifiers:
                raise GctaCojoError(
                    "GCTA output %s contains duplicate SNP ID %s."
                    % (path, identifier)
                )
            identifiers.add(identifier)
    return identifiers


def _run_selection_round(
    chromosomes: Sequence[int],
    *,
    workers: int,
    attempt_numbers: dict[int, int],
    build_attempt,
    configuration,
    module,
    logger,
    on_outcome: Callable[[_AttemptOutcome], None] | None = None,
) -> list[_AttemptOutcome]:
    prepared = []
    for chromosome in chromosomes:
        attempt_numbers[chromosome] = attempt_numbers.get(chromosome, 0) + 1
        attempt = attempt_numbers[chromosome]
        command, paths = build_attempt(chromosome, attempt)
        logger.record(
            "INPUT",
            "gcta_cojo_chromosome_command",
            chromosome=chromosome,
            attempt=attempt,
            command=command,
        )
        prepared.append((chromosome, attempt, command, paths))

    outcomes = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_attempt = {
            executor.submit(
                _run_slct_attempt,
                chromosome,
                attempt,
                command,
                paths,
                configuration,
                module,
            ): (chromosome, attempt)
            for chromosome, attempt, command, paths in prepared
        }
        for future in as_completed(future_to_attempt):
            outcome = future.result()
            outcomes.append(outcome)
            if outcome.error is None:
                logger.record(
                    "OUTPUT",
                    "gcta_cojo_chromosome_completed",
                    chromosome=outcome.chromosome,
                    attempt=outcome.attempt,
                    no_signals=outcome.no_signals,
                    peak_rss_bytes=outcome.peak_rss_bytes,
                    memory_samples=outcome.memory_samples,
                    process_tree_samples=outcome.process_tree_samples,
                    root_only_samples=outcome.root_only_samples,
                )
            else:
                logger.warning(
                    "GCTA-COJO chromosome %d attempt %d failed: %s: %s"
                    % (
                        outcome.chromosome,
                        outcome.attempt,
                        outcome.error_type,
                        outcome.error,
                    )
                )
            if on_outcome is not None:
                on_outcome(outcome)
    return sorted(outcomes, key=lambda item: int(item.chromosome))


def _run_selection_with_retries(
    chromosomes: Sequence[int],
    *,
    initial_workers: int,
    attempt_numbers: dict[int, int],
    build_attempt,
    configuration,
    module,
    logger,
    on_outcome: Callable[[_AttemptOutcome], None] | None = None,
) -> tuple[dict[int, _AttemptOutcome], list[_AttemptOutcome]]:
    pending = list(chromosomes)
    successes: dict[int, _AttemptOutcome] = {}
    history: list[_AttemptOutcome] = []
    maximum_rounds = module.chromosome_execution.failed_chromosome_retries + 1
    for round_number in range(1, maximum_rounds + 1):
        workers = (
            initial_workers
            if round_number == 1
            else min(module.chromosome_execution.retry_max_workers, len(pending))
        )
        outcomes = _run_selection_round(
            pending,
            workers=workers,
            attempt_numbers=attempt_numbers,
            build_attempt=build_attempt,
            configuration=configuration,
            module=module,
            logger=logger,
            on_outcome=on_outcome,
        )
        history.extend(outcomes)
        pending = []
        for outcome in outcomes:
            if outcome.error is None:
                successes[int(outcome.chromosome)] = outcome
            else:
                pending.append(int(outcome.chromosome))
        if not pending:
            return successes, history
        if round_number < maximum_rounds:
            logger.record(
                "ACTION",
                "gcta_cojo_chromosome_retry",
                retry_round=round_number,
                chromosomes=pending,
                workers=min(
                    module.chromosome_execution.retry_max_workers,
                    len(pending),
                ),
            )
    details = []
    for chromosome in pending:
        latest = next(
            outcome for outcome in reversed(history)
            if outcome.chromosome == chromosome
        )
        details.append(
            "chromosome %d: %s: %s"
            % (chromosome, latest.error_type, latest.error)
        )
    raise GctaCojoError(
        "GCTA-COJO failed after automatically retrying only the failed "
        "chromosome(s): %s" % "; ".join(details)
    )


def _read_table_for_merge(
    path: Path,
    *,
    expected_chromosome: int,
    module,
) -> tuple[list[str], list[list[str]], list[str]]:
    pattern = re.compile(module.input_validation.table_delimiter_pattern)
    schema = module.results.schemas["slct"]
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        header_line = handle.readline()
        header = pattern.split(header_line.strip()) if header_line else []
        if schema.chromosome_column not in header or schema.identifier_column not in header:
            raise GctaCojoError(
                "GCTA chromosome result %s lacks configured chromosome or SNP columns."
                % path
            )
        chromosome_index = header.index(schema.chromosome_column)
        identifier_index = header.index(schema.identifier_column)
        rows = []
        identifiers = []
        for line_number, raw in enumerate(handle, 2):
            if not raw.strip():
                continue
            fields = pattern.split(raw.strip())
            if len(fields) != len(header):
                raise GctaCojoError(
                    "GCTA chromosome result %s line %d has %d fields; expected %d."
                    % (path, line_number, len(fields), len(header))
                )
            try:
                observed_chromosome = int(fields[chromosome_index])
            except ValueError as exc:
                raise GctaCojoError(
                    "GCTA chromosome result %s line %d has non-numeric chromosome %s."
                    % (path, line_number, fields[chromosome_index])
                ) from exc
            if observed_chromosome != expected_chromosome:
                raise GctaCojoError(
                    "GCTA chromosome %d result contains a chromosome %d row: %s"
                    % (expected_chromosome, observed_chromosome, path)
                )
            rows.append(fields)
            identifiers.append(fields[identifier_index])
    if not rows:
        raise GctaCojoError(
            "GCTA chromosome %d primary result contains no selected SNPs."
            % expected_chromosome
        )
    if len(identifiers) != len(set(identifiers)):
        raise GctaCojoError(
            "GCTA chromosome %d primary result contains duplicate SNP IDs."
            % expected_chromosome
        )
    return header, rows, identifiers


def _read_ld_matrix(
    path: Path,
    expected_identifiers: Sequence[str],
    module,
) -> list[list[str]]:
    pattern = re.compile(module.input_validation.table_delimiter_pattern)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        header_line = handle.readline()
        header = pattern.split(header_line.strip()) if header_line else []
        if not header or header[0] != module.results.schemas["slct"].identifier_column:
            raise GctaCojoError("GCTA LD matrix has an invalid header: %s" % path)
        column_ids = header[1:]
        rows = []
        row_ids = []
        for line_number, raw in enumerate(handle, 2):
            if not raw.strip():
                continue
            fields = pattern.split(raw.strip())
            if len(fields) != len(header):
                raise GctaCojoError(
                    "GCTA LD matrix %s line %d has %d fields; expected %d."
                    % (path, line_number, len(fields), len(header))
                )
            row_ids.append(fields[0])
            numeric = []
            for value in fields[1:]:
                try:
                    parsed = float(value)
                except ValueError as exc:
                    raise GctaCojoError(
                        "GCTA LD matrix %s contains a non-numeric value."
                        % path
                    ) from exc
                if not math.isfinite(parsed):
                    raise GctaCojoError(
                        "GCTA LD matrix %s contains a non-finite value." % path
                    )
                numeric.append(value)
            rows.append(numeric)
    expected = list(expected_identifiers)
    if column_ids != expected or row_ids != expected:
        raise GctaCojoError(
            "GCTA LD matrix SNP order does not match its chromosome result: %s"
            % path
        )
    return rows


def _atomic_write_lines(path: Path, lines: Sequence[str], module) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + module.results.atomic_output_suffix)
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            for line in lines:
                handle.write(line)
                if not line.endswith("\n"):
                    handle.write("\n")
        temporary.replace(path)
    except OSError as exc:
        raise GctaCojoError("Cannot write combined GCTA output %s: %s" % (path, exc)) from exc
    finally:
        temporary.unlink(missing_ok=True)


def _merge_primary_and_ld(
    outcomes: Mapping[int, _AttemptOutcome],
    primary_destination: Path,
    ld_destination: Path,
    module,
) -> list[str]:
    delimiter = module.results.normalized_delimiter
    expected_header = None
    primary_lines = []
    blocks = []
    all_identifiers = []
    for chromosome in sorted(outcomes):
        outcome = outcomes[chromosome]
        if outcome.no_signals:
            continue
        header, rows, identifiers = _read_table_for_merge(
            outcome.paths.primary,
            expected_chromosome=chromosome,
            module=module,
        )
        if expected_header is None:
            expected_header = header
            primary_lines.append(delimiter.join(header))
        elif header != expected_header:
            raise GctaCojoError(
                "GCTA chromosome %d result header differs from earlier chromosomes."
                % chromosome
            )
        primary_lines.extend(delimiter.join(row) for row in rows)
        block = _read_ld_matrix(outcome.paths.ld_matrix, identifiers, module)
        blocks.append((identifiers, block))
        all_identifiers.extend(identifiers)
    if not all_identifiers:
        raise GctaCojoError("No chromosome produced a selected SNP to merge.")
    if len(all_identifiers) != len(set(all_identifiers)):
        raise GctaCojoError(
            "Parallel GCTA chromosome results contain duplicate SNP IDs."
        )
    _atomic_write_lines(primary_destination, primary_lines, module)

    ld_lines = [delimiter.join([
        module.results.schemas["slct"].identifier_column,
        *all_identifiers,
    ])]
    offset = 0
    total = len(all_identifiers)
    for identifiers, block in blocks:
        block_size = len(identifiers)
        for row_index, identifier in enumerate(identifiers):
            values = ["0"] * total
            values[offset:offset + block_size] = block[row_index]
            ld_lines.append(delimiter.join([identifier, *values]))
        offset += block_size
    _atomic_write_lines(ld_destination, ld_lines, module)
    return all_identifiers


def _write_selected_snps(path: Path, identifiers: Sequence[str], module) -> None:
    _atomic_write_lines(path, list(identifiers), module)


def _write_combined_log(
    destination: Path,
    history: Sequence[_AttemptOutcome],
    module,
) -> None:
    lines = []
    ordered = sorted(
        history,
        key=lambda item: (
            0 if item.phase == "selection" else 1,
            -1 if item.chromosome is None else item.chromosome,
            item.attempt,
        ),
    )
    for outcome in ordered:
        target = (
            "genome-wide conditional reconstruction"
            if outcome.chromosome is None
            else "chromosome %d" % outcome.chromosome
        )
        status = "COMPLETED" if outcome.error is None else "FAILED"
        lines.extend([
            "===== %s attempt %d: %s =====" % (target, outcome.attempt, status),
            "Command: %s" % " ".join(outcome.command),
        ])
        if outcome.paths.log.is_file():
            lines.append(
                outcome.paths.log.read_text(
                    encoding="utf-8", errors="replace",
                ).rstrip("\n")
            )
        elif outcome.error is not None:
            lines.append("%s: %s" % (outcome.error_type, outcome.error))
        lines.append("")
    _atomic_write_lines(destination, lines, module)


def _conditional_module(module, selected_snp_list: Path):
    document = module.model_dump(mode="python")
    document["mode"] = "cond"
    document["inputs"]["condition_snps"] = str(selected_snp_list)
    document["analysis"]["chromosome"] = None
    return type(module).model_validate(document)


def _print_execution_details(configuration, heading: str, fields) -> None:
    """Emit one aligned live/durable block when subordinate progress is enabled."""
    if not configuration.logging.show_progress:
        return
    label_width = configuration.logging.terminal_label_width
    lines = [screen_line("analysis", heading, indent=6)]
    lines.extend(
        screen_field(
            kind,
            label,
            value,
            indent=10,
            label_width=label_width,
        )
        for kind, label, value in fields
    )
    print_screen_block("\n".join(lines))


def execute_parallel_slct(
    *,
    configuration,
    module,
    logger,
    executable: str,
    summary_file: Path,
    reference_prefix: Path,
    staging: Path,
    combined_prefix: Path,
    combined_primary: Path,
    combined_ld: Path,
    combined_conditional: Path,
    combined_condition_snps: Path,
    combined_log: Path,
    normalized_destination: Path,
    effective_exclude: Path | None,
    chromosome_counts: Mapping[str, int],
    output_contract: GctaCojoParallelOutputContract,
    dry_run: bool,
) -> ParallelSlctResult:
    """Run, retry, and combine chromosome-wise stepwise-selection results."""
    workloads = numeric_chromosome_workloads(chromosome_counts)
    if workloads is None or len(workloads) < 2:
        raise GctaCojoError(
            "Parallel slct requires at least two numeric chromosome labels in "
            "the GCTA-usable .ma/BIM overlap."
        )
    chromosomes = sorted(workloads)
    execution = module.chromosome_execution
    worker_threads = execution.threads_per_worker
    attempt_numbers: dict[int, int] = {}

    def build_attempt(chromosome: int, attempt: int):
        prefix = configured_output_path(
            staging,
            module.output_layout.chromosome_output_prefix,
            error_type=GctaCojoError,
            dataset_id=configuration.run.dataset_id,
            mode=module.mode,
            chromosome=chromosome,
            attempt=attempt,
        )
        paths = _attempt_paths(
            prefix,
            combined_prefix=combined_prefix,
            combined_primary=combined_primary,
            combined_ld=combined_ld,
            combined_conditional=combined_conditional,
            combined_condition_snps=combined_condition_snps,
            combined_log=combined_log,
        )
        command = build_cojo_command(
            executable,
            summary_file,
            reference_prefix,
            prefix,
            worker_threads,
            module,
            exclude_snps=effective_exclude,
            chromosome=chromosome,
        )
        return command, paths

    if dry_run:
        commands = [build_attempt(chromosome, 1)[0] for chromosome in chromosomes]
        if output_contract == "complete":
            selected_snp_list = configured_output_path(
                staging,
                module.output_layout.selected_snp_list,
                error_type=GctaCojoError,
                dataset_id=configuration.run.dataset_id,
                mode=module.mode,
            )
            conditional_template_module = _conditional_module(
                module, selected_snp_list,
            )
            conditional_prefix = configured_output_path(
                staging,
                module.output_layout.conditional_output_prefix,
                error_type=GctaCojoError,
                dataset_id=configuration.run.dataset_id,
                mode=module.mode,
                attempt=1,
            )
            commands.append(build_cojo_command(
                executable,
                summary_file,
                reference_prefix,
                conditional_prefix,
                worker_threads,
                conditional_template_module,
                exclude_snps=effective_exclude,
            ))
        plan = plan_chromosome_workers(
            len(chromosomes),
            total_threads=configuration.execution.threads,
            threads_per_worker=worker_threads,
            memory_budget_gb=configuration.execution.memory_gb,
            parent_memory_gb=_parent_rss_gb(),
            estimated_memory_gb_per_worker=(
                execution.minimum_memory_gb_per_worker
            ),
            configured_max_workers=execution.max_workers,
        )
        return ParallelSlctResult(
            metrics={
                "execution_strategy": "chromosome_parallel",
                "chromosomes": chromosomes,
                "parallel_workers": plan.workers,
                "threads_per_worker": worker_threads,
                "memory_estimate_source": "configured_floor_dry_run",
                "estimated_memory_gb_per_worker": (
                    execution.minimum_memory_gb_per_worker
                ),
                "parallel_output_contract": output_contract,
                "conditional_reconstruction_status": (
                    "planned" if output_contract == "complete"
                    else "not_requested"
                ),
                "conditional_reconstruction_attempts": 0,
                "command": commands,
            },
            no_signals=False,
            dry_run=True,
        )

    pilot = min(chromosomes, key=lambda value: (-workloads[value], value))
    progress = MeasuredProgress(
        "GCTA-COJO chromosome progress",
        enabled=configuration.logging.show_progress,
    )
    phase_title = "Determine safe chromosome parallelism"
    completed_chromosomes: set[int] = set()

    def report_outcome(outcome: _AttemptOutcome) -> None:
        if outcome.error is None:
            completed_chromosomes.add(int(outcome.chromosome))
            progress.update(
                len(completed_chromosomes),
                total=len(chromosomes),
                title=phase_title,
            )
        else:
            progress.set_phase(
                "chr%d attempt %d failed · retry pending"
                % (outcome.chromosome, outcome.attempt)
            )

    progress.start(phase_title, total=len(chromosomes))
    try:
        pilot_successes, history = _run_selection_with_retries(
            [pilot],
            initial_workers=1,
            attempt_numbers=attempt_numbers,
            build_attempt=build_attempt,
            configuration=configuration,
            module=module,
            logger=logger,
            on_outcome=report_outcome,
        )
        pilot_outcome = pilot_successes[pilot]
        observed_pilot_gb = float(pilot_outcome.peak_rss_bytes) / _GIB
        memory_was_measured = bool(
            pilot_outcome.memory_samples > 0 and pilot_outcome.peak_rss_bytes > 0
        )
        estimated_worker_gb = max(
            execution.minimum_memory_gb_per_worker,
            observed_pilot_gb * execution.memory_safety_factor,
        )
        if not memory_was_measured:
            memory_estimate_source = "configured_floor_unmeasured"
        elif pilot_outcome.process_tree_samples > 0:
            memory_estimate_source = "measured_process_tree_pilot"
        else:
            memory_estimate_source = "measured_root_process_pilot"
        effective_max_workers = (
            execution.max_workers if memory_was_measured else 1
        )
        if not memory_was_measured:
            logger.warning(
                "PostGWAS could not measure the pilot GCTA process-tree memory. "
                "Chromosome execution will continue with one worker at a time "
                "instead of risking an unsupported parallel-memory estimate."
            )
        remaining = [value for value in chromosomes if value != pilot]
        plan = plan_chromosome_workers(
            len(remaining),
            total_threads=configuration.execution.threads,
            threads_per_worker=worker_threads,
            memory_budget_gb=configuration.execution.memory_gb,
            parent_memory_gb=_parent_rss_gb(),
            estimated_memory_gb_per_worker=estimated_worker_gb,
            configured_max_workers=effective_max_workers,
        )
        logger.record(
            "DECISION",
            "gcta_cojo_chromosome_workers",
            pilot_chromosome=pilot,
            pilot_usable_variants=workloads[pilot],
            pilot_peak_memory_gb=observed_pilot_gb,
            estimated_memory_gb_per_worker=estimated_worker_gb,
            memory_safety_factor=execution.memory_safety_factor,
            minimum_memory_gb_per_worker=execution.minimum_memory_gb_per_worker,
            memory_estimate_source=memory_estimate_source,
            requested_max_workers=execution.max_workers,
            workers=plan.workers,
            cpu_limit=plan.cpu_limit,
            memory_limit=plan.memory_limit,
            configured_limit=plan.configured_limit,
            threads_per_worker=plan.threads_per_worker,
            memory_budget_gb=plan.memory_budget_gb,
            parent_memory_gb=plan.parent_memory_gb,
        )
        phase_title = "Select chromosomes · up to %d at once" % plan.workers
        progress.set_phase(phase_title)
        _print_execution_details(
            configuration,
            "GCTA-COJO chromosome execution",
            [
                (
                    "count", "Chromosomes scheduled",
                    f"{len(chromosomes):,}",
                ),
                (
                    "count", "Chromosomes running simultaneously",
                    "up to %d" % plan.workers,
                ),
            ],
        )
        remaining_successes, remaining_history = _run_selection_with_retries(
            remaining,
            initial_workers=plan.workers,
            attempt_numbers=attempt_numbers,
            build_attempt=build_attempt,
            configuration=configuration,
            module=module,
            logger=logger,
            on_outcome=report_outcome,
        )
        history.extend(remaining_history)
    except BaseException:
        progress.fail(title="GCTA-COJO chromosome selection failed")
        raise
    else:
        progress.complete(
            len(chromosomes),
            total=len(chromosomes),
            title=phase_title,
        )
    successes = {**pilot_successes, **remaining_successes}
    no_signal_chromosomes = sorted(
        chromosome for chromosome, outcome in successes.items()
        if outcome.no_signals
    )
    all_no_signals = len(no_signal_chromosomes) == len(chromosomes)
    conditional_history: list[_AttemptOutcome] = []
    conditional_reconstruction_status = "not_applicable_no_signals"
    if all_no_signals:
        result_metrics, _ = write_empty_cojo_results(
            normalized_destination, module,
        )
    else:
        selected_identifiers = _merge_primary_and_ld(
            successes,
            combined_primary,
            combined_ld,
            module,
        )
        if output_contract == "complete":
            selected_snp_list = configured_output_path(
                staging,
                module.output_layout.selected_snp_list,
                error_type=GctaCojoError,
                dataset_id=configuration.run.dataset_id,
                mode=module.mode,
            )
            _write_selected_snps(selected_snp_list, selected_identifiers, module)
            selected_set = set(selected_identifiers)
            conditional_module = _conditional_module(module, selected_snp_list)
            conditional_outcome = None
            maximum_attempts = execution.failed_chromosome_retries + 1
            final_title = (
                "Rebuild conditional statistics · 1 process · %d thread/process"
                % worker_threads
            )
            final_progress = MeasuredProgress(
                "GCTA-COJO final model progress",
                enabled=configuration.logging.show_progress,
            )
            final_progress.start(final_title, total=1)
            _print_execution_details(
                configuration,
                "GCTA-COJO final model",
                [
                    (
                        "analysis", "Current phase",
                        "genome-wide conditional reconstruction",
                    ),
                    ("count", "GCTA processes", 1),
                    ("count", "Threads for this process", worker_threads),
                    (
                        "count",
                        "Selected SNPs being conditioned on",
                        len(selected_set),
                    ),
                ],
            )
            try:
                for attempt in range(1, maximum_attempts + 1):
                    conditional_prefix = configured_output_path(
                        staging,
                        module.output_layout.conditional_output_prefix,
                        error_type=GctaCojoError,
                        dataset_id=configuration.run.dataset_id,
                        mode=module.mode,
                        attempt=attempt,
                    )
                    conditional_paths = _attempt_paths(
                        conditional_prefix,
                        combined_prefix=combined_prefix,
                        combined_primary=combined_primary,
                        combined_ld=combined_ld,
                        combined_conditional=combined_conditional,
                        combined_condition_snps=combined_condition_snps,
                        combined_log=combined_log,
                    )
                    conditional_command = build_cojo_command(
                        executable,
                        summary_file,
                        reference_prefix,
                        conditional_prefix,
                        worker_threads,
                        conditional_module,
                        exclude_snps=effective_exclude,
                    )
                    logger.record(
                        "INPUT",
                        "gcta_cojo_conditional_reconstruction_command",
                        attempt=attempt,
                        command=conditional_command,
                    )
                    conditional_outcome = _run_conditional_attempt(
                        attempt,
                        conditional_command,
                        conditional_paths,
                        configuration,
                        module,
                        selected_set,
                    )
                    conditional_history.append(conditional_outcome)
                    if conditional_outcome.error is None:
                        break
                    logger.warning(
                        "GCTA-COJO conditional reconstruction attempt %d failed: "
                        "%s: %s"
                        % (
                            attempt,
                            conditional_outcome.error_type,
                            conditional_outcome.error,
                        )
                    )
                    if attempt < maximum_attempts:
                        logger.record(
                            "ACTION",
                            "gcta_cojo_conditional_reconstruction_retry",
                            next_attempt=attempt + 1,
                        )
                if (
                    conditional_outcome is None
                    or conditional_outcome.error is not None
                ):
                    raise GctaCojoError(
                        "GCTA-COJO genome-wide conditional output reconstruction "
                        "failed after %d attempt(s): %s"
                        % (
                            maximum_attempts,
                            (
                                None
                                if conditional_outcome is None
                                else conditional_outcome.error
                            ),
                        )
                    )
            except BaseException:
                final_progress.fail(title=final_title)
                raise
            else:
                final_progress.complete(1, total=1, title=final_title)
            combined_conditional.parent.mkdir(parents=True, exist_ok=True)
            temporary_conditional = combined_conditional.with_name(
                combined_conditional.name + module.results.atomic_output_suffix
            )
            try:
                shutil.copyfile(
                    conditional_outcome.paths.conditional,
                    temporary_conditional,
                )
                temporary_conditional.replace(combined_conditional)
            finally:
                temporary_conditional.unlink(missing_ok=True)
            conditional_reconstruction_status = "completed"
        else:
            logger.record(
                "SKIP",
                "gcta_cojo_conditional_reconstruction",
                reason="selection_only_output_contract",
                selected_signals=len(selected_identifiers),
                retained_outputs="combined_jma_and_ldr",
                omitted_output="genome_wide_cma",
                scientific_effect=(
                    "none_on_selected_signals_joint_statistics_or_physical_loci"
                ),
            )
            conditional_reconstruction_status = "not_requested"
        result_metrics, _ = normalize_cojo_results(
            combined_primary,
            normalized_destination,
            module,
        )

    history.extend(conditional_history)
    _write_combined_log(combined_log, history, module)
    peak_by_chromosome = {
        str(chromosome): float(outcome.peak_rss_bytes) / _GIB
        for chromosome, outcome in sorted(successes.items())
    }
    result_metrics.update({
        "execution_strategy": "chromosome_parallel",
        "chromosomes": chromosomes,
        "chromosome_usable_variants": {
            str(chromosome): workloads[chromosome]
            for chromosome in chromosomes
        },
        "pilot_chromosome": pilot,
        "pilot_peak_memory_gb": observed_pilot_gb,
        "pilot_process_tree_samples": pilot_outcome.process_tree_samples,
        "pilot_root_only_samples": pilot_outcome.root_only_samples,
        "memory_estimate_source": memory_estimate_source,
        "estimated_memory_gb_per_worker": estimated_worker_gb,
        "parallel_workers": plan.workers,
        "threads_per_worker": worker_threads,
        "worker_cpu_limit": plan.cpu_limit,
        "worker_memory_limit": plan.memory_limit,
        "memory_budget_gb": plan.memory_budget_gb,
        "parent_memory_gb": plan.parent_memory_gb,
        "chromosome_peak_memory_gb": peak_by_chromosome,
        "chromosome_attempts": {
            str(chromosome): attempt_numbers[chromosome]
            for chromosome in chromosomes
        },
        "retried_chromosomes": sorted(
            chromosome for chromosome, attempts in attempt_numbers.items()
            if attempts > 1
        ),
        "no_signal_chromosomes": no_signal_chromosomes,
        "parallel_output_contract": output_contract,
        "conditional_reconstruction_status": conditional_reconstruction_status,
        "conditional_reconstruction_attempts": len(conditional_history),
        "conditional_reconstruction_peak_memory_gb": max(
            (
                float(outcome.peak_rss_bytes) / _GIB
                for outcome in conditional_history
            ),
            default=0.0,
        ),
        "command": [list(outcome.command) for outcome in history],
    })
    return ParallelSlctResult(
        metrics=result_metrics,
        no_signals=all_no_signals,
        dry_run=False,
    )


__all__ = [
    "ChromosomeWorkerPlan",
    "ParallelSlctResult",
    "chromosome_parallel_slct_requested",
    "execute_parallel_slct",
    "numeric_chromosome_workloads",
    "plan_chromosome_workers",
    "should_parallelize_slct",
]
