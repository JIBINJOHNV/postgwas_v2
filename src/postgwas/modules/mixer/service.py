"""Official single-trait MiXeR and GSA-MiXeR execution."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

from postgwas.config import (
    load_configuration,
    load_run_configuration_for_module,
    write_resolved_configuration,
)
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.core.execution.runtime import validate_path
from postgwas.core.pipeline_logging import PipelineLogger, write_log_record
from postgwas.core.paths import configured_output_path, expand_token_path, resolve_executable
from postgwas.core.processes import build_container_command, run_checked_command
from postgwas.core.ui.screen import screen_field, screen_line
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
) -> None:
    """Require every configured chromosome resource before expensive analysis."""
    missing = [
        str(path)
        for chromosome in chromosomes
        if not (
            path := expand_chromosome_pattern(pattern, chromosome, placeholder)
        ).is_file()
        or path.stat().st_size <= 0
    ]
    if missing:
        preview = ", ".join(missing[:5])
        suffix = " (and %d more)" % (len(missing) - 5) if len(missing) > 5 else ""
        raise MixerError("Missing or empty %s reference files: %s%s" % (label, preview, suffix))


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


def _resolved_configuration(args):
    module_overrides = explicit_overrides(args, {
        "analysis": "analysis",
        "genome_build": "genome_build",
        "bim_file_pattern": "bim_file_pattern",
        "ld_file_pattern": "ld_file_pattern",
        "mixer_backend": "execution_backend",
        "gsa_annotation_file_pattern": "gsa.annotation_file_pattern",
        "gsa_loadlib_file_pattern": "gsa.loadlib_file_pattern",
        "gsa_baseline_go_file": "gsa.baseline_go_file",
        "gsa_model_go_file": "gsa.model_go_file",
        "gsa_test_go_file": "gsa.test_go_file",
    })
    global_overrides = explicit_overrides(args, {
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
    run_checked_command(
        command,
        purpose,
        logger=logger,
        error_type=MixerError,
        timeout_seconds=configuration.execution.timeout_seconds,
        expected_outputs=expected_outputs,
        dry_run=dry_run,
    )


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
    backend: str,
    input_metrics,
    input_warnings,
    *,
    first_step: int,
    total_steps: int,
    dry_run: bool,
) -> dict:
    module = configuration.modules.mixer
    layout = module.output_layout
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
    fit_log = _result_path(fit_prefix, ".log")
    test_json = _result_path(test_prefix, ".json")
    test_log = _result_path(test_prefix, ".log")
    seed = str(configuration.execution.random_seed)

    with logger.step(
        first_step, total_steps, "Fit single-trait MiXeR model", "mixer.fit1",
    ) as step:
        step.input("mixer_input", path=str(trait_file))
        _run_step(
            [
                *tool, module.workflow.fit_command, *common, *module.fit_arguments,
                "--trait1-file", str(trait_file), "--seed", seed,
                "--out", str(fit_prefix),
            ],
            "MiXeR fit1",
            [fit_json, fit_log],
            configuration,
            logger,
            dry_run=dry_run,
        )
        step.output("fit_parameters", path=str(fit_json))
        step.output("fit_log", path=str(fit_log))

    with logger.step(
        first_step + 1, total_steps, "Evaluate single-trait MiXeR model", "mixer.test1",
    ) as step:
        _run_step(
            [
                *tool, module.workflow.test_command, *common, *module.test_arguments,
                "--trait1-file", str(trait_file), "--load-params", str(fit_json),
                "--seed", seed, "--out", str(test_prefix),
            ],
            "MiXeR test1",
            [test_json, test_log],
            configuration,
            logger,
            dry_run=dry_run,
        )
        step.output("test_statistics", path=str(test_json))
        step.output("test_log", path=str(test_log))
    return {
        "mixer_input": str(trait_file),
        "fit": str(fit_json),
        "fit_log": str(fit_log),
        "test": str(test_json),
        "test_log": str(test_log),
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
    total_steps: int,
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
    with logger.step(
        first_step, total_steps, "Split one-trait summary statistics by chromosome",
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

    with logger.step(
        first_step + 1, total_steps, "Fit GSA-MiXeR baseline model", "mixer.gsa_base",
    ) as step:
        _run_step(
            [
                *tool, module.workflow.gsa_command, "--gsa-base", *common,
                "--go-file", resources["baseline_go"], "--out", str(baseline_prefix),
            ],
            "GSA-MiXeR baseline model",
            [baseline_json, baseline_log, baseline_snps, baseline_weights],
            configuration,
            logger,
            dry_run=dry_run,
        )
        step.output("gsa_baseline", path=str(baseline_json))

    with logger.step(
        first_step + 2, total_steps, "Fit GSA-MiXeR enrichment model", "mixer.gsa_full",
    ) as step:
        _run_step(
            [
                *tool, module.workflow.gsa_command, "--gsa-full", *common,
                "--go-file", resources["model_go"],
                "--go-file-test", resources["test_go"],
                "--load-params-file", str(baseline_json),
                "--load-baseline-params-file", str(baseline_json),
                "--calc-loglike-diff-go-test", settings.loglike_difference_method,
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


def run_single_trait_mixer(
    mixer_input_file: str | Path,
    output_directory: str | Path,
    dataset_id: str,
    configuration,
    logger: PipelineLogger,
    *,
    dry_run: bool = False,
) -> dict:
    """Run the configured single-trait MiXeR analysis selection."""
    module = configuration.modules.mixer
    trait_file = Path(mixer_input_file).expanduser().resolve()
    if not trait_file.is_file() or trait_file.stat().st_size <= 0:
        raise MixerError("MiXeR input file does not exist or is empty: %s" % trait_file)
    if module.bim_file_pattern is None:
        raise MixerError("--bim-file-pattern is required for every MiXeR analysis")
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
    bim_pattern = _absolute_pattern(module.bim_file_pattern)
    validate_reference_pattern(bim_pattern, "BIM", module.chromosomes, placeholder)
    if module.ld_file_pattern and (run_univariate or not module.gsa.loadlib_file_pattern):
        validate_reference_pattern(
            _absolute_pattern(module.ld_file_pattern), "LD", module.chromosomes, placeholder,
        )

    gsa_resources = None
    if run_gsa:
        gsa_resources = _gsa_resources(module)
        validate_reference_pattern(
            gsa_resources["annotation"], "GSA annotation", module.chromosomes, placeholder,
        )
        if "loadlib" in gsa_resources:
            validate_reference_pattern(
                gsa_resources["loadlib"], "GSA load-library", module.chromosomes, placeholder,
            )
        for key in ("baseline_go", "model_go", "test_go"):
            validate_gsa_go_file(gsa_resources[key], module.gsa, MixerError)

    output = Path(output_directory).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    mount_directories = {trait_file.parent, output, Path(bim_pattern).parent}
    if module.ld_file_pattern and (run_univariate or not module.gsa.loadlib_file_pattern):
        mount_directories.add(Path(_absolute_pattern(module.ld_file_pattern)).parent)
    if gsa_resources:
        mount_directories.update(Path(value).parent for value in gsa_resources.values())
    tool, backend = _tool_command(
        configuration,
        mount_directories=sorted(mount_directories, key=str),
        work_directory=output,
        validate=not dry_run,
    )
    logger.record(
        "PARAM", "mixer_backend", analysis=module.analysis,
        requested=module.execution_backend, selected=backend,
        container_image=(
            configuration.resources.containers.mixer.image if backend == "docker" else None
        ),
    )
    total_steps = (2 if run_univariate else 0) + (3 if run_gsa else 0)
    result = {
        "analysis": module.analysis,
        "mixer_input": str(trait_file),
        "execution_backend": backend,
        "input_metrics": input_metrics,
        "input_warnings": input_warnings,
        "dry_run": dry_run,
    }
    next_step = 1
    if run_univariate:
        result["univariate"] = _run_univariate(
            trait_file, output, dataset_id, configuration, logger, tool, backend,
            input_metrics, input_warnings, first_step=next_step,
            total_steps=total_steps, dry_run=dry_run,
        )
        next_step += 2
    if run_gsa:
        result["gsa"] = _run_gsa(
            trait_file, output, dataset_id, configuration, logger, tool, backend,
            input_metrics, gsa_resources, first_step=next_step,
            total_steps=total_steps, dry_run=dry_run,
        )
    return result


def _generate_figures(
    result: dict,
    output_directory: Path,
    dataset_id: str,
    configuration,
    logger: PipelineLogger,
) -> list[str]:
    module = configuration.modules.mixer
    prefix = configured_output_path(
        output_directory, module.output_layout.figure_prefix,
        error_type=MixerError, dataset_id=dataset_id,
    )
    prefix.parent.mkdir(parents=True, exist_ok=True)
    figures_tool, backend = _tool_command(
        configuration,
        figures=True,
        mount_directories=[output_directory],
        work_directory=output_directory,
    )
    if backend != result["execution_backend"]:
        raise MixerError("MiXeR and mixer_figures resolved to different execution backends")
    command = [
        *figures_tool, module.workflow.figure_command,
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


def run_mixer_direct(args, ctx=None):
    """Resolve configuration, run selected one-trait analyses, and close the log."""
    try:
        configuration = _resolved_configuration(args)
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
    try:
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
        )
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
                figures = (
                    _generate_figures(
                        univariate, output_directory, dataset_id, configuration, logger,
                    )
                    if configuration.modules.mixer.reporting.generate_figures else []
                )
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
        print("\n".join(rendered))
        return result
    except BaseException as exc:
        if not logger.summary()["failed"]:
            logger.error("MiXeR failed: %s: %s" % (type(exc).__name__, exc))
        raise
    finally:
        logger.close()


__all__ = [
    "MixerError",
    "MixerReference",
    "expand_chromosome_pattern",
    "run_mixer_direct",
    "run_single_trait_mixer",
    "validate_reference_pattern",
    "validate_standard_reference",
]
