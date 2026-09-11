"""Official single-trait MiXeR and GSA-MiXeR execution."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import re
from typing import Any, Callable, Mapping, Sequence

from postgwas.config import (
    load_configuration,
    load_run_configuration_for_module,
    write_resolved_configuration,
)
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.core.execution.runtime import validate_path
from postgwas.core.input_validation import record_file_validation, validate_once
from postgwas.core.io.delimiters import open_text
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.core.paths import configured_output_path, expand_token_path, resolve_executable
from postgwas.core.reference_resources import require_file_inventory
from postgwas.core.preflight import (
    PipelinePreflightEvidence,
    PreflightFileIdentity,
    capture_preflight_file_identities,
    pipeline_preflight_evidence,
    require_pipeline_input_vcf,
    require_unchanged_preflight_files,
)
from postgwas.core.processes import build_container_command, run_checked_command
from postgwas.core.required_arguments import (
    RequiredArgument,
    require_resolved_arguments,
)
from postgwas.core.ui import MeasuredProgress, PipelineStageController, StageProgress
from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.core.validation_reporting import (
    register_file_availability_bundle,
    register_file_validation_bundle,
)
from postgwas.modules.mixer.results import (
    build_gsa_summary,
    build_mixer_summary,
    inspect_mixer_input,
    render_gsa_summary,
    render_mixer_summary,
    validate_gsa_go_file,
    write_gsa_summary,
    write_mixer_summary,
)


class MixerError(RuntimeError):
    """A single-trait MiXeR run cannot proceed safely or a command failed."""


@dataclass(frozen=True)
class MixerReference:
    """Per-chromosome BIM and LD patterns consumed by univariate MiXeR."""

    bim_file_pattern: str
    ld_file_pattern: str


@dataclass(frozen=True)
class MixerPipelineResources:
    """External MiXeR resources validated before formatter table creation."""

    configuration: Any
    bim_file_pattern: str
    ld_file_pattern: str | None
    fit_extract_files: tuple[Path, ...] | None
    gsa_resources: Mapping[str, str] | None
    tool: tuple[str, ...]
    figures_tool: tuple[str, ...] | None
    backend: str
    file_identities: tuple[PreflightFileIdentity, ...]


class _MixerNativeProgress:
    """Report only work units observed in an active upstream MiXeR log."""

    _STAGES = (
        "Initialize the reference and GWAS summary statistics",
        "Load chromosome LD reference data",
        "Optimize and validate the MiXeR model",
    )
    _FIT_SEQUENCE_PATTERN = re.compile(r"--fit-sequence:\s*\[([^]]*)\]")
    _FIT_SEQUENCE_ITEM_PATTERN = re.compile(r"['\"]([-A-Za-z0-9_.]+)['\"]")

    def __init__(
        self,
        log_path: Path,
        purpose: str,
        ld_files: Sequence[Path],
        configuration,
        logger: PipelineLogger,
    ) -> None:
        self.log_path = log_path
        self.purpose = purpose
        self.ld_files = tuple(str(path) for path in ld_files)
        self.logger = logger
        self.progress = StageProgress(
            "%s execution progress" % purpose,
            enabled=configuration.logging.show_progress,
        )
        self.progress.start_step(1, len(self._STAGES), self._STAGES[0])
        self.current_stage = 1
        self.offset = 0
        self.pending = ""
        self.started_ld_files: set[str] = set()
        self.ld_progress: MeasuredProgress | None = None
        self.optimizer_progress: MeasuredProgress | None = None
        self.optimizer_evaluations = 0
        self.optimizer_phase: str | None = None
        self.optimizer_sequence: tuple[str, ...] = ()
        self.optimizer_title: str | None = None
        self.closed = False

    def _advance(self, target_stage: int) -> None:
        while self.current_stage < target_stage:
            self.progress.complete_step(
                self.current_stage,
                len(self._STAGES),
                self._STAGES[self.current_stage - 1],
            )
            self.current_stage += 1
            self.progress.start_step(
                self.current_stage,
                len(self._STAGES),
                self._STAGES[self.current_stage - 1],
            )

    def _observe_ld_start(self, line: str) -> None:
        matched = next((path for path in self.ld_files if path in line), None)
        if matched is None or matched in self.started_ld_files:
            return
        self._advance(2)
        self.started_ld_files.add(matched)
        if self.ld_progress is None:
            self.ld_progress = MeasuredProgress(
                "MiXeR chromosome LD-loading progress",
                enabled=self.progress.enabled,
            )
            self.ld_progress.start(
                "Load chromosome LD references",
                total=len(self.ld_files),
            )
        self.ld_progress.update(
            max(0, len(self.started_ld_files) - 1),
            total=len(self.ld_files),
            title="Load chromosome LD references",
        )

    def _complete_ld_loading(self) -> None:
        self._advance(2)
        if self.ld_progress is not None:
            self.ld_progress.complete(
                len(self.ld_files),
                total=len(self.ld_files),
                title="Load chromosome LD references",
            )
            self.ld_progress = None
        self._advance(3)

    @classmethod
    def _fit_sequence(cls, line: str) -> tuple[str, ...]:
        """Read only the optimizer names announced by the native log."""
        match = cls._FIT_SEQUENCE_PATTERN.search(line)
        if match is None:
            return ()
        return tuple(cls._FIT_SEQUENCE_ITEM_PATTERN.findall(match.group(1)))

    @staticmethod
    def _optimizer_progress_title(
        phase: str | None,
        sequence: Sequence[str],
    ) -> str:
        title = "Observed MiXeR cost-function evaluations"
        if phase and phase in sequence:
            title += " · %s (%d/%d: %s)" % (
                phase,
                tuple(sequence).index(phase) + 1,
                len(sequence),
                " → ".join(sequence),
            )
        elif phase:
            title += " · %s" % phase
        elif sequence:
            title += " · sequence: %s" % " → ".join(sequence)
        return title

    def _observe_optimizer(self, lines: Sequence[str]) -> None:
        added_evaluations = sum(
            "<calc_" in line and ", cost=" in line for line in lines
        )
        phase = self.optimizer_phase
        sequence = self.optimizer_sequence
        for line in lines:
            observed_sequence = self._fit_sequence(line)
            if observed_sequence and observed_sequence != sequence:
                sequence = observed_sequence
                self.logger.record(
                    "OBSERVED",
                    "mixer_optimizer_sequence",
                    purpose=self.purpose,
                    optimization_sequence=list(sequence),
                )
            marker = "fit_type=="
            if marker in line and " done " not in line:
                candidate = line.split(marker, 1)[1].split("...", 1)[0].strip()
                if candidate:
                    phase = candidate
            elif "Calculate AIC/BIC w.r.t. infinitesimal model" in line:
                phase = "infinitesimal comparison"
        title = self._optimizer_progress_title(phase, sequence)
        if not added_evaluations:
            if self.optimizer_progress is not None and title != self.optimizer_title:
                self.optimizer_progress.set_phase(title)
            self.optimizer_phase = phase
            self.optimizer_sequence = sequence
            self.optimizer_title = title
            return
        self._advance(3)
        self.optimizer_evaluations += added_evaluations
        if self.optimizer_progress is None:
            self.optimizer_progress = MeasuredProgress(
                "MiXeR optimizer activity",
                enabled=self.progress.enabled,
            )
            self.optimizer_progress.start(title)
        elif title != self.optimizer_title:
            self.optimizer_progress.set_phase(title)
        self.optimizer_progress.update(self.optimizer_evaluations, title=title)
        self.optimizer_phase = phase
        self.optimizer_sequence = sequence
        self.optimizer_title = title
        self.logger.record(
            "OBSERVED",
            "mixer_optimizer_progress",
            purpose=self.purpose,
            completed_cost_evaluations=self.optimizer_evaluations,
            total_cost_evaluations="unknown_until_convergence",
            optimization_phase=phase,
            optimization_sequence=list(sequence) or None,
        )

    def observe(self) -> None:
        """Read only newly appended complete native-log lines."""
        try:
            with self.log_path.open("r", encoding="utf-8", errors="replace") as handle:
                handle.seek(self.offset)
                chunk = handle.read()
                self.offset = handle.tell()
        except FileNotFoundError:
            return
        if not chunk:
            return
        text = self.pending + chunk
        lines = text.splitlines(keepends=True)
        self.pending = ""
        if lines and not lines[-1].endswith(("\n", "\r")):
            self.pending = lines.pop()
        complete_lines = [line.rstrip("\r\n") for line in lines]
        for line in complete_lines:
            if "<init(" in line and "elapsed time" in line:
                self._advance(2)
            if ">load_ld_matrix(filename=" in line:
                self._observe_ld_start(line)
            if "--fit-sequence:" in line:
                self._complete_ld_loading()
        self._observe_optimizer(complete_lines)

    def complete(self) -> None:
        """Reach 100% only after checked-command output validation succeeds."""
        if self.closed:
            return
        self.observe()
        if self.ld_progress is not None:
            self.ld_progress.complete(
                len(self.ld_files),
                total=len(self.ld_files),
                title="Load chromosome LD references",
            )
            self.ld_progress = None
        if self.optimizer_progress is not None:
            self.optimizer_progress.complete(
                self.optimizer_evaluations,
                total=self.optimizer_evaluations,
                title=self.optimizer_title,
            )
            self.optimizer_progress = None
        self._advance(len(self._STAGES))
        self.progress.complete_step(
            len(self._STAGES), len(self._STAGES), self._STAGES[-1],
        )
        self.logger.record(
            "OBSERVED",
            "mixer_native_progress_complete",
            purpose=self.purpose,
            chromosome_ld_loads_observed=len(self.started_ld_files),
            chromosome_ld_files_expected=len(self.ld_files),
            cost_function_evaluations=self.optimizer_evaluations,
        )
        self.closed = True

    def fail(self) -> None:
        """Preserve the last observed work count below completion."""
        if self.closed:
            return
        self.observe()
        if self.optimizer_progress is not None:
            self.optimizer_progress.fail(title=self.optimizer_title)
            self.optimizer_progress = None
        elif self.ld_progress is not None:
            self.ld_progress.fail(title="Load chromosome LD references")
            self.ld_progress = None
        self.progress.fail_step(
            self.current_stage,
            len(self._STAGES),
            self._STAGES[self.current_stage - 1],
        )
        self.closed = True


def _result_path(prefix: Path, suffix: str) -> Path:
    """Apply one fixed upstream ``--out PREFIX`` result suffix."""
    return Path("%s%s" % (prefix, suffix))


def expand_chromosome_pattern(
    pattern: str, chromosome: str, placeholder: str,
) -> Path:
    return expand_token_path(pattern, placeholder, chromosome, error_type=MixerError)


def validate_reference_pattern(
    pattern: str,
    label: str,
    chromosomes: Sequence[str],
    placeholder: str,
) -> tuple[Path, ...]:
    """Require every configured chromosome resource before expensive analysis."""
    paths = tuple(
        expand_chromosome_pattern(pattern, chromosome, placeholder)
        for chromosome in chromosomes
    )
    require_file_inventory(
        paths, "%s reference" % label,
        missing_message="Missing or empty %s reference files" % label,
        error_type=MixerError,
    )
    register_file_availability_bundle(
        paths,
        "MiXeR chromosome resource bundle",
        (
            ("analysis", "mixer_resource_type", label),
            ("count", "chromosomes", list(chromosomes)),
            (
                "success",
                "mixer_chromosome_files",
                "%d / %d" % (len(paths), len(paths)),
            ),
            ("info", "mixer_pattern", pattern, True),
        ),
    )
    return paths


def validate_standard_reference(
    reference: MixerReference,
    chromosomes: Sequence[str],
    placeholder: str,
) -> None:
    validate_reference_pattern(
        reference.bim_file_pattern, "BIM", chromosomes, placeholder,
    )
    validate_reference_pattern(
        reference.ld_file_pattern, "LD", chromosomes, placeholder,
    )


def validate_fit_extract_pattern(
    pattern: str,
    replicate_indices: Sequence[int],
    placeholder: str,
) -> tuple[Path, ...]:
    """Validate the official one-identifier-per-line fit subsets."""
    paths = tuple(
        expand_token_path(
            pattern, placeholder, str(index), error_type=MixerError,
        )
        for index in replicate_indices
    )
    require_file_inventory(
        paths,
        "MiXeR fit extract file",
        missing_message="Missing or empty MiXeR fit extract files",
        error_type=MixerError,
    )
    counts = []
    for path in paths:
        checks = ("one unique SNP identifier per non-empty line",)

        def inspect_extract_file(current_path=path):
            identifiers = set()
            try:
                with open_text(current_path) as handle:
                    for line_number, line in enumerate(handle, 1):
                        identifier = line.strip()
                        if not identifier:
                            continue
                        if any(character.isspace() for character in identifier):
                            raise MixerError(
                                "MiXeR fit extract file must contain exactly one "
                                "SNP identifier per non-empty line; invalid line "
                                "%d: %s" % (line_number, current_path)
                            )
                        if identifier in identifiers:
                            raise MixerError(
                                "MiXeR fit extract file contains duplicate SNP "
                                "identifier %s at line %d: %s"
                                % (identifier, line_number, current_path)
                            )
                        identifiers.add(identifier)
            except MixerError as exc:
                record_file_validation(
                    current_path,
                    "MiXeR fit extract file",
                    checks=checks,
                    status="failed",
                    message=str(exc),
                )
                raise
            except (OSError, UnicodeError) as exc:
                message = "Cannot read MiXeR fit extract file %s: %s" % (
                    current_path,
                    exc,
                )
                record_file_validation(
                    current_path,
                    "MiXeR fit extract file",
                    checks=checks,
                    status="failed",
                    message=message,
                )
                raise MixerError(message) from exc
            if not identifiers:
                message = (
                    "MiXeR fit extract file contains no SNP identifiers: %s"
                    % current_path
                )
                record_file_validation(
                    current_path,
                    "MiXeR fit extract file",
                    checks=checks,
                    status="failed",
                    message=message,
                )
                raise MixerError(message)
            count = len(identifiers)
            record_file_validation(
                current_path,
                "MiXeR fit extract file",
                checks=checks,
                metrics={"snp_identifiers": count},
                message=(
                    "Identifier syntax and within-file uniqueness were validated; "
                    "the reference bundle must supply the MAF/LD-pruning provenance."
                ),
            )
            return count

        counts.append(validate_once(
            (path,),
            {"validator": "mixer_fit_extract_identifiers", "version": 1},
            inspect_extract_file,
            error_type=MixerError,
        ))
    register_file_validation_bundle(
        paths,
        "MiXeR replicated fit extract bundle",
        (
            ("count", "mixer_fit_replicates", len(paths)),
            (
                "count",
                "mixer_fit_snp_identifiers_minimum",
                min(counts),
            ),
            (
                "count",
                "mixer_fit_snp_identifiers_maximum",
                max(counts),
            ),
            ("info", "mixer_pattern", pattern, True),
        ),
        covered_checks=("one unique SNP identifier per non-empty line",),
        covered_metric_keys=("snp_identifiers",),
    )
    return paths


def _resolved_configuration(args):
    module_overrides = explicit_overrides(args, {
        "analysis": "analysis",
        "genome_build": "genome_build",
        "bim_file_pattern": "bim_file_pattern",
        "ld_file_pattern": "ld_file_pattern",
        "mixer_fit_extract_file_pattern": (
            "univariate.fit_extract_file_pattern"
        ),
        "mixer_backend": "execution_backend",
        "gsa_annotation_file_pattern": "gsa.annotation_file_pattern",
        "gsa_loadlib_file_pattern": "gsa.loadlib_file_pattern",
        "gsa_baseline_go_file": "gsa.baseline_go_file",
        "gsa_model_go_file": "gsa.model_go_file",
        "gsa_test_go_file": "gsa.test_go_file",
    })
    global_overrides = explicit_overrides(args, {
        "dataset_id": "run.dataset_id",
        "output_directory": "run.output_directory",
        "threads": "execution.threads",
        "memory_gb": "execution.memory_gb",
        "seed": "execution.random_seed",
        "mixer": "resources.executables.mixer",
        "mixer_figures": "resources.executables.mixer_figures",
        "mixer_container_runtime": "resources.containers.mixer.runtime",
        "mixer_container_image": "resources.containers.mixer.image",
        "mixer_container_platform": "resources.containers.mixer.platform",
        "resume": "run.resume",
        "overwrite": "run.overwrite",
    })
    return load_run_configuration_for_module(
        "mixer",
        getattr(args, "run_config", None),
        module_overrides=module_overrides,
        global_overrides=global_overrides,
    )


def _native_tool_command(configuration, *, figures=False) -> list[str]:
    resources = configuration.resources.executables
    script_value = resources.mixer_figures if figures else resources.mixer
    script = resolve_executable(
        str(script_value), "MiXeR figures script" if figures else "MiXeR script",
        error_type=MixerError,
    )
    python = resolve_executable(
        str(resources.python), "Python executable", error_type=MixerError,
    )
    return [python, script] if script.lower().endswith(".py") else [script]


def _tool_command(
    configuration,
    *,
    figures=False,
    mount_directories: Sequence[str | Path] = (),
    work_directory: str | Path | None = None,
    validate=True,
) -> tuple[list[str], str]:
    backend = configuration.modules.mixer.execution_backend
    if backend in {"auto", "native"}:
        try:
            return _native_tool_command(configuration, figures=figures), "native"
        except MixerError:
            if backend == "native":
                raise

    container = configuration.resources.containers.mixer
    script = container.mixer_figures if figures else container.mixer
    try:
        command = build_container_command(
            container.runtime,
            container.image,
            [container.python, script],
            mount_directories=mount_directories,
            platform=container.platform,
            work_directory=work_directory,
            use_host_user=container.use_host_user,
            validate=validate,
        )
    except RuntimeError as exc:
        raise MixerError(str(exc)) from exc
    return command, "docker"


def _absolute_pattern(pattern: str) -> str:
    return str(Path(pattern).expanduser().resolve())


def _required_file(value: str | None, label: str) -> str:
    if value is None:
        raise MixerError("%s is required for the selected analysis" % label)
    return str(validate_path(must_exist=True, must_be_file=True, must_not_be_empty=True)(value))


def _required_pattern(value: str | None, label: str) -> str:
    if value is None:
        raise MixerError("%s is required for the selected analysis" % label)
    return _absolute_pattern(value)


def _run_step(
    command: Sequence[str],
    purpose: str,
    expected_outputs: Sequence[Path],
    configuration,
    logger: PipelineLogger,
    *,
    dry_run: bool,
    progress_factory: Callable[[], _MixerNativeProgress] | None = None,
) -> None:
    complete = all(path.is_file() and path.stat().st_size > 0 for path in expected_outputs)
    existing = [path for path in expected_outputs if path.exists()]
    if complete:
        if configuration.run.resume and not configuration.run.overwrite:
            logger.record(
                "SKIP", purpose, reason="resume",
                outputs=[str(path) for path in expected_outputs],
            )
            return
        if not configuration.run.overwrite:
            raise MixerError(
                "Output already exists: %s. Use --resume to keep completed steps or "
                "--overwrite to run them again." % expected_outputs[0]
            )
    elif existing and not configuration.run.overwrite:
        raise MixerError(
            "Incomplete MiXeR output exists: %s. Use --overwrite to replace the "
            "partial result." % ", ".join(str(path) for path in existing)
        )
    native_progress = (
        None if dry_run or progress_factory is None else progress_factory()
    )
    try:
        run_checked_command(
            command,
            purpose,
            logger=logger,
            error_type=MixerError,
            timeout_seconds=configuration.execution.timeout_seconds,
            expected_outputs=expected_outputs,
            dry_run=dry_run,
            progress_callback=(
                None if native_progress is None else native_progress.observe
            ),
            progress_refresh_seconds=(
                None
                if native_progress is None
                else configuration.logging.progress_refresh_seconds
            ),
        )
    except BaseException:
        if native_progress is not None:
            native_progress.fail()
        raise
    else:
        if native_progress is not None:
            native_progress.complete()


def _chromosome_argument(chromosomes: Sequence[str]) -> str:
    return chromosomes[0] if len(chromosomes) == 1 else "%s-%s" % (
        chromosomes[0], chromosomes[-1],
    )


def _validate_input_chromosomes(input_metrics, configuration, logger) -> list[str]:
    quality = configuration.modules.mixer.reporting.quality
    warnings = []
    for values_key, action, message in (
        (
            "chromosomes_missing", quality.missing_chromosome_action,
            "Configured chromosomes absent from MiXeR input: %s",
        ),
        (
            "chromosomes_unexpected", quality.unexpected_chromosome_action,
            "MiXeR input contains chromosomes outside the configured range: %s",
        ),
    ):
        values = input_metrics[values_key]
        if not values or action == "ignore":
            continue
        issue = message % ", ".join(values)
        if action == "error":
            raise MixerError(issue)
        warnings.append(issue)
        logger.warning(issue)
    return warnings


def _run_univariate(
    trait_file: Path,
    output: Path,
    dataset_id: str,
    configuration,
    logger: PipelineLogger,
    tool: Sequence[str],
    figures_tool: Sequence[str],
    backend: str,
    input_metrics,
    input_warnings,
    fit_extract_files: Sequence[Path],
    *,
    first_step: int,
    progress_first_step: int,
    total_steps: int,
    progress: PipelineStageController,
    dry_run: bool,
) -> dict:
    module = configuration.modules.mixer
    layout = module.output_layout
    replicate_indices = tuple(module.univariate.replicate_indices)
    if len(fit_extract_files) != len(replicate_indices):
        raise MixerError(
            "MiXeR fit extract-file count does not match the configured "
            "replicate-index count"
        )
    reference = MixerReference(
        _absolute_pattern(module.bim_file_pattern),
        _absolute_pattern(module.ld_file_pattern),
    )
    common = [
        "--ld-file", reference.ld_file_pattern,
        "--bim-file", reference.bim_file_pattern,
        "--chr2use", _chromosome_argument(module.chromosomes),
        "--threads", str(configuration.execution.threads),
    ]
    fit_prefix = configured_output_path(
        output, layout.fit_prefix, error_type=MixerError, dataset_id=dataset_id,
    )
    test_prefix = configured_output_path(
        output, layout.test_prefix, error_type=MixerError, dataset_id=dataset_id,
    )
    fit_prefix.parent.mkdir(parents=True, exist_ok=True)
    test_prefix.parent.mkdir(parents=True, exist_ok=True)
    fit_json = _result_path(fit_prefix, ".json")
    test_json = _result_path(test_prefix, ".json")
    seed = str(configuration.execution.random_seed)
    ld_files = tuple(
        expand_chromosome_pattern(
            reference.ld_file_pattern,
            chromosome,
            module.workflow.chromosome_placeholder,
        )
        for chromosome in module.chromosomes
    )

    fit_replicates = []
    test_replicates = []
    replicate_total = len(replicate_indices)

    progress.start(progress_first_step)
    fit_progress = MeasuredProgress(
        "MiXeR replicated fit progress",
        enabled=configuration.logging.show_progress,
    )
    fit_title = "Fit random MAF/LD-pruned SNP subsets"
    fit_progress.start(fit_title, total=replicate_total)
    try:
        with logger.step(
            first_step,
            total_steps,
            "Fit replicated single-trait MiXeR models",
            "mixer.fit1_replicates",
        ) as step:
            step.input("mixer_input", path=str(trait_file))
            for ordinal, (replicate, extract_file) in enumerate(
                zip(replicate_indices, fit_extract_files), 1,
            ):
                replicate_prefix = configured_output_path(
                    output,
                    layout.fit_replicate_prefix,
                    error_type=MixerError,
                    dataset_id=dataset_id,
                    replicate=replicate,
                )
                replicate_prefix.parent.mkdir(parents=True, exist_ok=True)
                replicate_json = _result_path(replicate_prefix, ".json")
                replicate_log = _result_path(replicate_prefix, ".log")
                logger.record(
                    "PARAM",
                    "mixer_fit_replicate",
                    replicate_index=replicate,
                    fit_extract_file=str(extract_file),
                    fit_uses_extract=True,
                    test_uses_extract=False,
                )
                _run_step(
                    [
                        *tool,
                        module.workflow.fit_command,
                        *common,
                        *module.fit_arguments,
                        "--extract",
                        str(extract_file),
                        "--trait1-file",
                        str(trait_file),
                        "--seed",
                        seed,
                        "--out",
                        str(replicate_prefix),
                    ],
                    "MiXeR fit1 replicate %s" % replicate,
                    [replicate_json, replicate_log],
                    configuration,
                    logger,
                    dry_run=dry_run,
                    progress_factory=lambda path=replicate_log, index=replicate: (
                        _MixerNativeProgress(
                            path,
                            "MiXeR fit1 replicate %s" % index,
                            ld_files,
                            configuration,
                            logger,
                        )
                    ),
                )
                fit_replicates.append({
                    "index": replicate,
                    "extract_file": str(extract_file),
                    "json": str(replicate_json),
                    "log": str(replicate_log),
                })
                fit_progress.update(ordinal, total=replicate_total, title=fit_title)
            step.output(
                "fit_replicates",
                files=[item["json"] for item in fit_replicates],
            )
    except BaseException:
        fit_progress.fail(title=fit_title)
        progress.fail_active()
        raise
    else:
        fit_progress.complete(replicate_total, total=replicate_total, title=fit_title)
        progress.complete(progress_first_step)

    progress.start(progress_first_step + 1)
    test_progress = MeasuredProgress(
        "MiXeR replicated test progress",
        enabled=configuration.logging.show_progress,
    )
    test_title = "Evaluate fitted models on the full SNP set"
    test_progress.start(test_title, total=replicate_total)
    try:
        with logger.step(
            first_step + 1,
            total_steps,
            "Evaluate replicated single-trait MiXeR models",
            "mixer.test1_replicates",
        ) as step:
            for ordinal, fit_replicate in enumerate(fit_replicates, 1):
                replicate = fit_replicate["index"]
                replicate_prefix = configured_output_path(
                    output,
                    layout.test_replicate_prefix,
                    error_type=MixerError,
                    dataset_id=dataset_id,
                    replicate=replicate,
                )
                replicate_prefix.parent.mkdir(parents=True, exist_ok=True)
                replicate_json = _result_path(replicate_prefix, ".json")
                replicate_log = _result_path(replicate_prefix, ".log")
                _run_step(
                    [
                        *tool,
                        module.workflow.test_command,
                        *common,
                        *module.test_arguments,
                        "--trait1-file",
                        str(trait_file),
                        "--load-params",
                        fit_replicate["json"],
                        "--seed",
                        seed,
                        "--out",
                        str(replicate_prefix),
                    ],
                    "MiXeR test1 replicate %s" % replicate,
                    [replicate_json, replicate_log],
                    configuration,
                    logger,
                    dry_run=dry_run,
                    progress_factory=lambda path=replicate_log, index=replicate: (
                        _MixerNativeProgress(
                            path,
                            "MiXeR test1 replicate %s" % index,
                            ld_files,
                            configuration,
                            logger,
                        )
                    ),
                )
                test_replicates.append({
                    "index": replicate,
                    "json": str(replicate_json),
                    "log": str(replicate_log),
                })
                test_progress.update(
                    ordinal, total=replicate_total, title=test_title,
                )
            step.output(
                "test_replicates",
                files=[item["json"] for item in test_replicates],
            )
    except BaseException:
        test_progress.fail(title=test_title)
        progress.fail_active()
        raise
    else:
        test_progress.complete(
            replicate_total, total=replicate_total, title=test_title,
        )
        progress.complete(progress_first_step + 1)

    replicate_placeholder = module.workflow.replicate_placeholder
    replicate_selection = ",".join(str(value) for value in replicate_indices)
    fit_json_pattern = _result_path(configured_output_path(
        output,
        layout.fit_replicate_prefix,
        error_type=MixerError,
        dataset_id=dataset_id,
        replicate=replicate_placeholder,
    ), ".json")
    test_json_pattern = _result_path(configured_output_path(
        output,
        layout.test_replicate_prefix,
        error_type=MixerError,
        dataset_id=dataset_id,
        replicate=replicate_placeholder,
    ), ".json")

    progress.start(progress_first_step + 2)
    try:
        with logger.step(
            first_step + 2,
            total_steps,
            "Combine replicated MiXeR fit estimates",
            "mixer.combine_fit1",
        ) as step:
            _run_step(
                [
                    *figures_tool,
                    module.workflow.combine_command,
                    "--json",
                    str(fit_json_pattern),
                    "--rep2use",
                    replicate_selection,
                    "--out",
                    str(fit_prefix),
                ],
                "Combine MiXeR fit1 replicates",
                [fit_json],
                configuration,
                logger,
                dry_run=dry_run,
            )
            step.output("combined_fit_parameters", path=str(fit_json))
    except BaseException:
        progress.fail_active()
        raise
    else:
        progress.complete(progress_first_step + 2)

    progress.start(progress_first_step + 3)
    try:
        with logger.step(
            first_step + 3,
            total_steps,
            "Combine replicated MiXeR test diagnostics",
            "mixer.combine_test1",
        ) as step:
            _run_step(
                [
                    *figures_tool,
                    module.workflow.combine_command,
                    "--json",
                    str(test_json_pattern),
                    "--rep2use",
                    replicate_selection,
                    "--out",
                    str(test_prefix),
                ],
                "Combine MiXeR test1 replicates",
                [test_json],
                configuration,
                logger,
                dry_run=dry_run,
            )
            step.output("combined_test_statistics", path=str(test_json))
    except BaseException:
        progress.fail_active()
        raise
    else:
        progress.complete(progress_first_step + 3)

    return {
        "mixer_input": str(trait_file),
        "fit": str(fit_json),
        "test": str(test_json),
        "fit_replicates": fit_replicates,
        "test_replicates": test_replicates,
        "fit_extract_file_pattern": _absolute_pattern(
            module.univariate.fit_extract_file_pattern
        ),
        "replicate_indices": list(replicate_indices),
        "genome_build": module.genome_build,
        "execution_backend": backend,
        "input_metrics": input_metrics,
        "input_warnings": input_warnings,
        "dry_run": dry_run,
    }


def _gsa_resources(module) -> dict[str, str]:
    settings = module.gsa
    resources = {
        "annotation": _required_pattern(
            settings.annotation_file_pattern, "--gsa-annotation-file-pattern",
        ),
        "baseline_go": _required_file(
            settings.baseline_go_file, "--gsa-baseline-go-file",
        ),
        "model_go": _required_file(settings.model_go_file, "--gsa-model-go-file"),
        "test_go": _required_file(settings.test_go_file, "--gsa-test-go-file"),
    }
    if settings.loadlib_file_pattern:
        resources["loadlib"] = _absolute_pattern(settings.loadlib_file_pattern)
    return resources


def _require_mixer_pipeline_arguments(configuration) -> None:
    """Validate all mode-dependent external MiXeR requirements together."""
    module = configuration.modules.mixer
    run_univariate = module.analysis in {"univariate", "all"}
    run_gsa = module.analysis in {"gsa", "all"}
    requirements = [RequiredArgument(
        "--bim-file-pattern", "modules.mixer.bim_file_pattern",
        module.bim_file_pattern,
    )]
    if run_univariate or (run_gsa and not module.gsa.loadlib_file_pattern):
        requirements.append(RequiredArgument(
            "--ld-file-pattern", "modules.mixer.ld_file_pattern",
            module.ld_file_pattern,
        ))
    if run_univariate:
        requirements.append(RequiredArgument(
            "--mixer-fit-extract-file-pattern",
            "modules.mixer.univariate.fit_extract_file_pattern",
            module.univariate.fit_extract_file_pattern,
        ))
    if run_gsa:
        requirements.extend([
            RequiredArgument(
                "--gsa-annotation-file-pattern",
                "modules.mixer.gsa.annotation_file_pattern",
                module.gsa.annotation_file_pattern,
            ),
            RequiredArgument(
                "--gsa-baseline-go-file",
                "modules.mixer.gsa.baseline_go_file",
                module.gsa.baseline_go_file,
            ),
            RequiredArgument(
                "--gsa-model-go-file",
                "modules.mixer.gsa.model_go_file",
                module.gsa.model_go_file,
            ),
            RequiredArgument(
                "--gsa-test-go-file",
                "modules.mixer.gsa.test_go_file",
                module.gsa.test_go_file,
            ),
        ])
    require_resolved_arguments(requirements)


def _existing_output_mount_directory(output_directory: str | Path) -> Path:
    """Return the nearest existing non-root parent for a future output path."""
    candidate = Path(output_directory).expanduser().resolve()
    while not candidate.is_dir() and candidate != Path(candidate.anchor):
        candidate = candidate.parent
    if candidate == Path(candidate.anchor):
        raise MixerError(
            "MiXeR container execution requires the configured output path to "
            "have an existing parent below the filesystem root."
        )
    return candidate


def _mixer_pipeline_execution_configuration(args, resources):
    """Retain preflight settings while applying the pipeline stage directory."""
    run = resources.configuration.run.model_copy(update={
        "output_directory": Path(args.output_directory).expanduser().resolve(),
    })
    return resources.configuration.model_copy(update={"run": run}, deep=True)


def preflight_mixer_pipeline(
    args,
    *,
    preflight_evidence=None,
) -> PipelinePreflightEvidence:
    """Validate entry data, references, and the selected MiXeR backend early."""
    entry_vcf = require_pipeline_input_vcf(preflight_evidence)
    configuration = _resolved_configuration(args)
    module = configuration.modules.mixer
    _require_mixer_pipeline_arguments(configuration)
    observed_build = str(entry_vcf["harmonised"]["genome_build"])
    if str(module.genome_build) != observed_build:
        raise MixerError(
            "The harmonised GWAS-VCF declares genome build %s, but MiXeR "
            "resolves to %s. The GWAS, BIM, LD, and annotation resources "
            "must use the same build." % (observed_build, module.genome_build)
        )

    run_univariate = module.analysis in {"univariate", "all"}
    run_gsa = module.analysis in {"gsa", "all"}
    placeholder = module.workflow.chromosome_placeholder
    bim_pattern = _absolute_pattern(module.bim_file_pattern)
    resource_files = list(validate_reference_pattern(
        bim_pattern, "BIM", module.chromosomes, placeholder,
    ))
    ld_pattern = None
    if run_univariate or (run_gsa and not module.gsa.loadlib_file_pattern):
        ld_pattern = _absolute_pattern(module.ld_file_pattern)
        resource_files.extend(validate_reference_pattern(
            ld_pattern, "LD", module.chromosomes, placeholder,
        ))

    fit_extract_files = None
    if run_univariate:
        fit_extract_files = validate_fit_extract_pattern(
            _absolute_pattern(module.univariate.fit_extract_file_pattern),
            module.univariate.replicate_indices,
            module.workflow.replicate_placeholder,
        )
        resource_files.extend(fit_extract_files)

    gsa_resources = None
    if run_gsa:
        gsa_resources = _gsa_resources(module)
        resource_files.extend(validate_reference_pattern(
            gsa_resources["annotation"], "GSA annotation",
            module.chromosomes, placeholder,
        ))
        if "loadlib" in gsa_resources:
            resource_files.extend(validate_reference_pattern(
                gsa_resources["loadlib"], "GSA load-library",
                module.chromosomes, placeholder,
            ))
        for key in ("baseline_go", "model_go", "test_go"):
            validate_gsa_go_file(gsa_resources[key], module.gsa, MixerError)
            resource_files.append(Path(gsa_resources[key]).expanduser().resolve())

    mount_directories = {path.parent for path in resource_files}
    mount_directories.add(_existing_output_mount_directory(
        configuration.run.output_directory,
    ))
    dry_run = bool(getattr(args, "dry_run", False))
    tool, backend = _tool_command(
        configuration,
        mount_directories=sorted(mount_directories, key=str),
        validate=not dry_run,
    )
    figures_tool = None
    if run_univariate:
        figures_tool, figures_backend = _tool_command(
            configuration,
            figures=True,
            mount_directories=sorted(mount_directories, key=str),
            validate=not dry_run,
        )
        if figures_backend != backend:
            raise MixerError(
                "MiXeR and mixer_figures resolved to different execution backends"
            )

    executable_files = [
        Path(value).expanduser().resolve()
        for value in (*tool, *(figures_tool or ()))
        if Path(value).expanduser().is_file()
    ]
    resources = MixerPipelineResources(
        configuration=configuration,
        bim_file_pattern=bim_pattern,
        ld_file_pattern=ld_pattern,
        fit_extract_files=fit_extract_files,
        gsa_resources=gsa_resources,
        tool=tuple(tool),
        figures_tool=None if figures_tool is None else tuple(figures_tool),
        backend=backend,
        file_identities=capture_preflight_file_identities(
            [*resource_files, *executable_files],
            error_type=MixerError,
            label="MiXeR resource",
        ),
    )
    return pipeline_preflight_evidence(
        "mixer",
        preflight_evidence,
        resources=resources,
        deferred_checks=(
            "Validate the formatter-created MiXeR summary-statistics table.",
            "Validate its chromosome coverage against the resolved analysis range.",
        ),
    )


def _gsa_common_arguments(
    module,
    configuration,
    split_pattern: Path,
    resources: dict[str, str],
) -> list[str]:
    settings = module.gsa
    arguments = [
        "--trait1-file", str(split_pattern),
        "--bim-file", _absolute_pattern(module.bim_file_pattern),
        "--annot-file", resources["annotation"],
        "--chr2use", _chromosome_argument(module.chromosomes),
        "--exclude-ranges", *settings.exclude_ranges,
        "--hardprune-maf", str(settings.hardprune_maf),
        "--hardprune-r2", str(settings.hardprune_r2),
        "--go-extend-bp", str(settings.go_extend_bp),
        "--go-all-genes-label", settings.go_all_genes_label,
        "--seed", str(configuration.execution.random_seed),
        "--threads", str(configuration.execution.threads),
    ]
    if settings.use_complete_tag_indices:
        arguments.append("--use-complete-tag-indices")
    if "loadlib" in resources:
        arguments.extend(["--loadlib-file", resources["loadlib"]])
    else:
        arguments.extend(["--ld-file", _absolute_pattern(module.ld_file_pattern)])
    if settings.z_max is not None:
        arguments.extend(["--z1max", str(settings.z_max)])
    if settings.adam_epoch is not None:
        arguments.extend(["--adam-epoch", *(str(value) for value in settings.adam_epoch)])
        arguments.extend([
            "--adam-step", *(str(value) for value in (settings.adam_step or [])),
        ])
    return arguments


def _run_gsa(
    trait_file: Path,
    output: Path,
    dataset_id: str,
    configuration,
    logger: PipelineLogger,
    tool: Sequence[str],
    backend: str,
    input_metrics,
    resources,
    *,
    first_step: int,
    progress_first_step: int,
    total_steps: int,
    progress: PipelineStageController,
    dry_run: bool,
) -> dict:
    module = configuration.modules.mixer
    settings = module.gsa
    layout = module.output_layout
    placeholder = module.workflow.chromosome_placeholder
    split_pattern = configured_output_path(
        output, layout.gsa_split_pattern, error_type=MixerError, dataset_id=dataset_id,
    )
    split_pattern.parent.mkdir(parents=True, exist_ok=True)
    split_outputs = [
        Path(str(split_pattern).replace(placeholder, chromosome))
        for chromosome in module.chromosomes
    ]
    progress.start(progress_first_step)
    try:
        with logger.step(
            first_step,
            total_steps,
            "Split one-trait summary statistics by chromosome",
            "mixer.split_sumstats",
        ) as step:
            _run_step(
                [
                    *tool, module.workflow.split_sumstats_command,
                    "--trait1-file", str(trait_file), "--out", str(split_pattern),
                    "--chr2use", _chromosome_argument(module.chromosomes),
                ],
                "GSA-MiXeR split_sumstats",
                split_outputs,
                configuration,
                logger,
                dry_run=dry_run,
            )
            step.output("chromosome_sumstats", pattern=str(split_pattern))
    except BaseException:
        progress.fail_active()
        raise
    else:
        progress.complete(progress_first_step)

    baseline_prefix = configured_output_path(
        output, layout.gsa_baseline_prefix, error_type=MixerError, dataset_id=dataset_id,
    )
    full_prefix = configured_output_path(
        output, layout.gsa_full_prefix, error_type=MixerError, dataset_id=dataset_id,
    )
    baseline_prefix.parent.mkdir(parents=True, exist_ok=True)
    full_prefix.parent.mkdir(parents=True, exist_ok=True)
    baseline_json = _result_path(baseline_prefix, ".json")
    baseline_log = _result_path(baseline_prefix, ".log")
    baseline_snps = _result_path(baseline_prefix, ".snps.csv")
    baseline_weights = _result_path(baseline_prefix, ".weights")
    full_json = _result_path(full_prefix, ".json")
    full_log = _result_path(full_prefix, ".log")
    enrichment_results = _result_path(full_prefix, ".go_test_enrich.csv")
    common = _gsa_common_arguments(
        module, configuration, split_pattern, resources,
    )

    progress.start(progress_first_step + 1)
    try:
        with logger.step(
            first_step + 1,
            total_steps,
            "Fit GSA-MiXeR baseline model",
            "mixer.gsa_base",
        ) as step:
            _run_step(
                [
                    *tool, module.workflow.gsa_command, "--gsa-base", *common,
                    "--go-file", resources["baseline_go"],
                    "--out", str(baseline_prefix),
                ],
                "GSA-MiXeR baseline model",
                [baseline_json, baseline_log, baseline_snps, baseline_weights],
                configuration,
                logger,
                dry_run=dry_run,
            )
            step.output("gsa_baseline", path=str(baseline_json))
    except BaseException:
        progress.fail_active()
        raise
    else:
        progress.complete(progress_first_step + 1)

    progress.start(progress_first_step + 2)
    try:
        with logger.step(
            first_step + 2,
            total_steps,
            "Fit GSA-MiXeR enrichment model",
            "mixer.gsa_full",
        ) as step:
            _run_step(
                [
                    *tool, module.workflow.gsa_command, "--gsa-full", *common,
                    "--go-file", resources["model_go"],
                    "--go-file-test", resources["test_go"],
                    "--load-params-file", str(baseline_json),
                    "--load-baseline-params-file", str(baseline_json),
                    "--calc-loglike-diff-go-test",
                    settings.loglike_difference_method,
                    "--se-samples", str(settings.standard_error_samples),
                    "--out", str(full_prefix),
                ],
                "GSA-MiXeR enrichment model",
                [full_json, full_log, enrichment_results],
                configuration,
                logger,
                dry_run=dry_run,
            )
            step.output("gsa_enrichment", path=str(enrichment_results))
    except BaseException:
        progress.fail_active()
        raise
    else:
        progress.complete(progress_first_step + 2)
    return {
        "mixer_input": str(trait_file),
        "split_sumstats_pattern": str(split_pattern),
        "baseline_json": str(baseline_json),
        "baseline_log": str(baseline_log),
        "full_json": str(full_json),
        "full_log": str(full_log),
        "enrichment_results": str(enrichment_results),
        "test_go_file": resources["test_go"],
        "genome_build": module.genome_build,
        "execution_backend": backend,
        "input_metrics": input_metrics,
        "dry_run": dry_run,
    }


def _mixer_progress_stages(module, *, include_reporting: bool) -> tuple[str, ...]:
    stages = ["Validate the MiXeR input and resolved resources"]
    if module.analysis in {"univariate", "all"}:
        stages.extend((
            "Fit replicated single-trait MiXeR models",
            "Evaluate replicated models on the full SNP set",
            "Combine replicated MiXeR fit estimates",
            "Combine replicated MiXeR test diagnostics",
        ))
    if module.analysis in {"gsa", "all"}:
        stages.extend((
            "Split summary statistics by chromosome for GSA-MiXeR",
            "Fit the GSA-MiXeR baseline model",
            "Fit the GSA-MiXeR enrichment model",
        ))
    if include_reporting:
        stages.append("Validate and report the selected MiXeR results")
    return tuple(stages)


def _run_single_trait_mixer_impl(
    mixer_input_file: str | Path,
    output_directory: str | Path,
    dataset_id: str,
    configuration,
    logger: PipelineLogger,
    *,
    dry_run: bool = False,
    pipeline_resources: MixerPipelineResources | None = None,
    progress: PipelineStageController,
) -> dict:
    """Run the configured single-trait MiXeR analysis selection."""
    module = configuration.modules.mixer
    if progress.current == 0:
        progress.start(1)
    trait_file = Path(mixer_input_file).expanduser().resolve()
    if not trait_file.is_file() or trait_file.stat().st_size <= 0:
        raise MixerError("MiXeR input file does not exist or is empty: %s" % trait_file)
    _require_mixer_pipeline_arguments(configuration)
    run_univariate = module.analysis in {"univariate", "all"}
    run_gsa = module.analysis in {"gsa", "all"}
    if (run_univariate or (run_gsa and not module.gsa.loadlib_file_pattern)) and not module.ld_file_pattern:
        raise MixerError(
            "--ld-file-pattern is required for univariate MiXeR and for GSA-MiXeR "
            "when --gsa-loadlib-file-pattern is not provided"
        )

    formatting = configuration.modules.formatting
    mixer_schema = formatting.exports["mixer"]
    chromosome_source = formatting.canonical_columns.chromosome
    input_metrics = inspect_mixer_input(
        trait_file,
        delimiter=formatting.runtime.table_delimiter,
        required_columns=list(mixer_schema.columns.values()),
        chromosome_column=mixer_schema.columns[chromosome_source],
        configured_chromosomes=module.chromosomes,
        error_type=MixerError,
    )
    input_warnings = _validate_input_chromosomes(input_metrics, configuration, logger)
    logger.record("OBSERVED", "mixer_input", **input_metrics)

    placeholder = module.workflow.chromosome_placeholder
    bim_pattern = (
        pipeline_resources.bim_file_pattern
        if pipeline_resources is not None
        else _absolute_pattern(module.bim_file_pattern)
    )
    if pipeline_resources is None:
        validate_reference_pattern(
            bim_pattern, "BIM", module.chromosomes, placeholder,
        )
        if module.ld_file_pattern and (
            run_univariate or not module.gsa.loadlib_file_pattern
        ):
            validate_reference_pattern(
                _absolute_pattern(module.ld_file_pattern), "LD",
                module.chromosomes, placeholder,
            )

    fit_extract_files = (
        None
        if not run_univariate
        else (
            pipeline_resources.fit_extract_files
            if pipeline_resources is not None
            else validate_fit_extract_pattern(
                _absolute_pattern(module.univariate.fit_extract_file_pattern),
                module.univariate.replicate_indices,
                module.workflow.replicate_placeholder,
            )
        )
    )

    gsa_resources = (
        None
        if not run_gsa
        else (
            dict(pipeline_resources.gsa_resources or {})
            if pipeline_resources is not None
            else _gsa_resources(module)
        )
    )
    if run_gsa and pipeline_resources is None:
        validate_reference_pattern(
            gsa_resources["annotation"], "GSA annotation",
            module.chromosomes, placeholder,
        )
        if "loadlib" in gsa_resources:
            validate_reference_pattern(
                gsa_resources["loadlib"], "GSA load-library",
                module.chromosomes, placeholder,
            )
        for key in ("baseline_go", "model_go", "test_go"):
            validate_gsa_go_file(gsa_resources[key], module.gsa, MixerError)

    output = Path(output_directory).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    mount_directories = {trait_file.parent, output, Path(bim_pattern).parent}
    if module.ld_file_pattern and (run_univariate or not module.gsa.loadlib_file_pattern):
        mount_directories.add(Path(_absolute_pattern(module.ld_file_pattern)).parent)
    if fit_extract_files:
        mount_directories.update(path.parent for path in fit_extract_files)
    if gsa_resources:
        mount_directories.update(Path(value).parent for value in gsa_resources.values())
    if pipeline_resources is None:
        tool, backend = _tool_command(
            configuration,
            mount_directories=sorted(mount_directories, key=str),
            work_directory=output,
            validate=not dry_run,
        )
    else:
        tool = list(pipeline_resources.tool)
        backend = pipeline_resources.backend
    if run_univariate:
        if pipeline_resources is not None and pipeline_resources.figures_tool is not None:
            figures_tool = list(pipeline_resources.figures_tool)
            figures_backend = pipeline_resources.backend
        else:
            figures_tool, figures_backend = _tool_command(
                configuration,
                figures=True,
                mount_directories=sorted(mount_directories, key=str),
                work_directory=output,
                validate=not dry_run,
            )
        if figures_backend != backend:
            raise MixerError(
                "MiXeR and mixer_figures resolved to different execution backends"
            )
    else:
        figures_tool = []
    logger.record(
        "PARAM", "mixer_backend", analysis=module.analysis,
        requested=module.execution_backend, selected=backend,
        container_image=(
            configuration.resources.containers.mixer.image if backend == "docker" else None
        ),
    )
    progress.complete(
        1,
        outcome_fields=[
            ("count", "Validated input variants", input_metrics["formatted_variants"]),
            ("count", "Configured chromosomes", len(module.chromosomes)),
            ("success", "Execution backend", backend),
        ],
    )
    total_steps = (4 if run_univariate else 0) + (3 if run_gsa else 0)
    result = {
        "analysis": module.analysis,
        "mixer_input": str(trait_file),
        "execution_backend": backend,
        "input_metrics": input_metrics,
        "input_warnings": input_warnings,
        "dry_run": dry_run,
    }
    next_step = 1
    next_progress_step = 2
    if run_univariate:
        result["univariate"] = _run_univariate(
            trait_file, output, dataset_id, configuration, logger, tool,
            figures_tool, backend, input_metrics, input_warnings,
            fit_extract_files, first_step=next_step,
            progress_first_step=next_progress_step,
            total_steps=total_steps, progress=progress, dry_run=dry_run,
        )
        next_step += 4
        next_progress_step += 4
    if run_gsa:
        result["gsa"] = _run_gsa(
            trait_file, output, dataset_id, configuration, logger, tool, backend,
            input_metrics, gsa_resources, first_step=next_step,
            progress_first_step=next_progress_step,
            total_steps=total_steps, progress=progress, dry_run=dry_run,
        )
    return result


def run_single_trait_mixer(
    mixer_input_file: str | Path,
    output_directory: str | Path,
    dataset_id: str,
    configuration,
    logger: PipelineLogger,
    *,
    dry_run: bool = False,
    pipeline_resources: MixerPipelineResources | None = None,
    progress: PipelineStageController | None = None,
) -> dict:
    """Run MiXeR with a caller-owned or invocation-local progress controller."""
    owns_progress = progress is None
    controller = progress or PipelineStageController(
        "MiXeR analysis progress",
        _mixer_progress_stages(
            configuration.modules.mixer,
            include_reporting=False,
        ),
    )
    try:
        return _run_single_trait_mixer_impl(
            mixer_input_file,
            output_directory,
            dataset_id,
            configuration,
            logger,
            dry_run=dry_run,
            pipeline_resources=pipeline_resources,
            progress=controller,
        )
    except BaseException:
        controller.fail_active()
        raise
    finally:
        if owns_progress:
            controller.close()


def _generate_figures(
    result: dict,
    output_directory: Path,
    dataset_id: str,
    configuration,
    logger: PipelineLogger,
    *,
    figures_tool: Sequence[str] | None = None,
    expected_backend: str | None = None,
) -> list[str]:
    module = configuration.modules.mixer
    prefix = configured_output_path(
        output_directory, module.output_layout.figure_prefix,
        error_type=MixerError, dataset_id=dataset_id,
    )
    prefix.parent.mkdir(parents=True, exist_ok=True)
    if figures_tool is None:
        resolved_tool, backend = _tool_command(
            configuration,
            figures=True,
            mount_directories=[output_directory],
            work_directory=output_directory,
        )
    else:
        resolved_tool = list(figures_tool)
        backend = expected_backend or result["execution_backend"]
    if backend != result["execution_backend"]:
        raise MixerError("MiXeR and mixer_figures resolved to different execution backends")
    command = [
        *resolved_tool, module.workflow.figure_command,
        "--json", result["test"], "--trait1", dataset_id,
        "--out", str(prefix), "--ext", *module.reporting.figure_extensions,
        "--statistic", *module.reporting.figure_statistics,
    ]
    expected = _result_path(prefix, ".csv")
    _run_step(
        command,
        "Generate official MiXeR diagnostic figures",
        [expected],
        configuration,
        logger,
        dry_run=False,
    )
    generated = sorted(
        path for path in prefix.parent.glob("%s.*" % prefix.name)
        if path.is_file() and path.stat().st_size > 0
    )
    logger.record("OUTPUT", "mixer_figures", files=[str(path) for path in generated])
    return [str(path) for path in generated]


def _dry_run_summary(dataset_id: str, result: dict, log_path: Path) -> str:
    return "\n".join([
        "",
        screen_line("analysis", "Single-trait MiXeR commands validated", indent=2),
        screen_field("info", "Dataset", dataset_id, indent=6, label_width=24),
        screen_field("info", "Selected analysis", result["analysis"], indent=6, label_width=24),
        screen_field("info", "Execution backend", result["execution_backend"], indent=6, label_width=24),
        screen_field("info", "Full log", log_path, indent=6, label_width=24),
        "",
    ])


def run_mixer_direct(
    args,
    ctx=None,
    *,
    pipeline_resources: MixerPipelineResources | None = None,
):
    """Resolve configuration, run selected one-trait analyses, and close the log."""
    try:
        configuration = (
            _mixer_pipeline_execution_configuration(args, pipeline_resources)
            if pipeline_resources is not None
            else _resolved_configuration(args)
        )
    except BaseException as exc:
        fallback = load_configuration()
        failure_output = Path(
            getattr(args, "output_directory", None) or fallback.run.output_directory
        ).expanduser().resolve()
        failure_dataset = str(
            getattr(args, "dataset_id", None) or fallback.run.dataset_id
        ).strip()
        failure_run_id = datetime.now(timezone.utc).strftime(
            fallback.modules.mixer.workflow.run_id_format
        )
        failure_log = configured_output_path(
            failure_output,
            fallback.modules.mixer.output_layout.log_file,
            error_type=MixerError,
            dataset_id=failure_dataset,
            run_id=failure_run_id,
        )
        write_log_record(
            failure_log,
            "ERROR",
            "MiXeR configuration failed: %s: %s" % (type(exc).__name__, exc),
            sample_id=failure_dataset,
            file_level=fallback.logging.file_level,
            screen_level=fallback.logging.console_level,
        )
        raise
    output_directory = Path(
        getattr(args, "output_directory", None) or configuration.run.output_directory
    ).expanduser().resolve()
    dataset_id = str(
        getattr(args, "dataset_id", None) or configuration.run.dataset_id
    ).strip()
    run_id = datetime.now(timezone.utc).strftime(
        configuration.modules.mixer.workflow.run_id_format
    )
    if not run_id or Path(run_id).name != run_id:
        raise MixerError("workflow.run_id_format produced an unsafe run identifier")
    mixer_input = getattr(args, "mixer_input_file", None)
    layout = configuration.modules.mixer.output_layout
    log_path = configured_output_path(
        output_directory, layout.log_file, error_type=MixerError,
        dataset_id=dataset_id, run_id=run_id,
    )
    logger = PipelineLogger(
        dataset_id,
        "run",
        str(log_path.parent),
        level=configuration.logging.file_level,
        screen_level=configuration.logging.console_level,
        log_path=str(log_path),
    )
    mixer_progress = PipelineStageController(
        "MiXeR analysis progress",
        _mixer_progress_stages(
            configuration.modules.mixer,
            include_reporting=True,
        ),
    )
    mixer_progress.start(1)
    try:
        if pipeline_resources is not None:
            require_unchanged_preflight_files(
                pipeline_resources.file_identities,
                error_type=MixerError,
                label="MiXeR resource",
            )
        if not mixer_input:
            raise MixerError("Provide the formatter-created file with --mixer-input-file PATH.")
        write_resolved_configuration(
            configuration,
            configured_output_path(
                output_directory, layout.resolved_config_file,
                error_type=MixerError, dataset_id=dataset_id,
            ),
            modules=("mixer",),
            resource_paths=(
                "executables.python",
                "executables.mixer",
                "executables.mixer_figures",
                "containers.mixer",
            ),
        )
        logger.record("INPUT", "mixer_run", dataset=dataset_id, input=str(mixer_input))
        logger.record(
            "PARAM", "execution", analysis=configuration.modules.mixer.analysis,
            threads=configuration.execution.threads,
            memory_gb=configuration.execution.memory_gb,
            seed=configuration.execution.random_seed,
            genome_build=configuration.modules.mixer.genome_build,
            requested_backend=configuration.modules.mixer.execution_backend,
        )
        result = run_single_trait_mixer(
            mixer_input,
            output_directory,
            dataset_id,
            configuration,
            logger,
            dry_run=bool(getattr(args, "dry_run", False)),
            pipeline_resources=pipeline_resources,
            progress=mixer_progress,
        )
        final_progress_step = len(mixer_progress.stages)
        mixer_progress.start(final_progress_step)
        result["log_file"] = str(log_path)
        result["run_id"] = run_id
        rendered = []
        formatter_metrics = (
            (ctx.get("formatter") or {}).get("mixer") if ctx is not None else None
        )
        if not result["dry_run"] and configuration.modules.mixer.reporting.enabled:
            if "univariate" in result:
                univariate = result["univariate"]
                univariate["log_file"] = str(log_path)
                summary = build_mixer_summary(
                    dataset_id=dataset_id,
                    run_id=run_id,
                    result=univariate,
                    input_metrics=result["input_metrics"],
                    input_warnings=result["input_warnings"],
                    formatter_metrics=formatter_metrics,
                    configuration=configuration,
                    error_type=MixerError,
                )
                figures = (
                    _generate_figures(
                        univariate, output_directory, dataset_id, configuration, logger,
                        figures_tool=(
                            None
                            if pipeline_resources is None
                            else pipeline_resources.figures_tool
                        ),
                        expected_backend=(
                            None
                            if pipeline_resources is None
                            else pipeline_resources.backend
                        ),
                    )
                    if configuration.modules.mixer.reporting.generate_figures else []
                )
                summary_yaml = configured_output_path(
                    output_directory, layout.summary_yaml,
                    error_type=MixerError, dataset_id=dataset_id,
                )
                summary_tsv = configured_output_path(
                    output_directory, layout.summary_tsv,
                    error_type=MixerError, dataset_id=dataset_id,
                )
                summary["artifacts"].update({
                    "summary_yaml": str(summary_yaml),
                    "summary_tsv": str(summary_tsv),
                    "figures": figures,
                })
                univariate["summary_yaml"], univariate["summary_tsv"] = write_mixer_summary(
                    summary, summary_yaml, summary_tsv,
                    configuration.modules.mixer.reporting,
                )
                univariate["figures"] = figures
                univariate["quality_status"] = summary["quality"]["status"]
                univariate["warnings"] = summary["quality"]["warnings"]
                logger.record("RESULT", "mixer_architecture", **{
                    key: value for key, value in summary["architecture"].items()
                    if key != "uncertainty"
                })
                logger.record("RESULT", "mixer_model_fit", **summary["model_fit"])
                for warning in summary["quality"]["warnings"]:
                    if warning not in result["input_warnings"]:
                        logger.warning(warning)
                rendered.append(render_mixer_summary(summary))
            if "gsa" in result:
                gsa = result["gsa"]
                gsa["log_file"] = str(log_path)
                gsa_summary, top_records = build_gsa_summary(
                    dataset_id=dataset_id,
                    run_id=run_id,
                    result=gsa,
                    input_metrics=result["input_metrics"],
                    configuration=configuration,
                    error_type=MixerError,
                )
                gsa_yaml = configured_output_path(
                    output_directory, layout.gsa_summary_yaml,
                    error_type=MixerError, dataset_id=dataset_id,
                )
                gsa_tsv = configured_output_path(
                    output_directory, layout.gsa_top_results_tsv,
                    error_type=MixerError, dataset_id=dataset_id,
                )
                gsa_summary["artifacts"].update({
                    "summary_yaml": str(gsa_yaml), "top_results_tsv": str(gsa_tsv),
                })
                gsa["summary_yaml"], gsa["top_results_tsv"] = write_gsa_summary(
                    gsa_summary, top_records, gsa_yaml, gsa_tsv, configuration,
                )
                logger.record(
                    "RESULT", "gsa_gene_sets", **gsa_summary["gene_set_assessment"],
                )
                rendered.append(render_gsa_summary(gsa_summary))
        elif result["dry_run"]:
            rendered.append(_dry_run_summary(dataset_id, result, log_path))
        else:
            rendered.append("\n".join([
                "",
                screen_line("analysis", "Single-trait MiXeR analysis completed", indent=2),
                screen_field(
                    "info", "Dataset", dataset_id, indent=6, label_width=24,
                ),
                screen_field(
                    "info", "Selected analysis", result["analysis"],
                    indent=6, label_width=24,
                ),
                screen_field(
                    "info", "Full log", log_path, indent=6, label_width=24,
                ),
                "",
            ]))
        if ctx is not None:
            ctx["mixer"] = result
        logger.record("STATUS", "mixer_run", status="COMPLETED")
        mixer_progress.complete(
            final_progress_step,
            outcome_fields=[
                ("success", "Selected analysis", result["analysis"]),
                ("success", "Canonical log", log_path),
            ],
        )
        print("\n".join(rendered))
        return result
    except BaseException as exc:
        mixer_progress.fail_active()
        if not logger.summary()["failed"]:
            logger.error("MiXeR failed: %s: %s" % (type(exc).__name__, exc))
        raise
    finally:
        mixer_progress.close()
        logger.close()


__all__ = [
    "MixerError",
    "MixerPipelineResources",
    "MixerReference",
    "expand_chromosome_pattern",
    "preflight_mixer_pipeline",
    "run_mixer_direct",
    "run_single_trait_mixer",
    "validate_reference_pattern",
    "validate_fit_extract_pattern",
    "validate_standard_reference",
]
