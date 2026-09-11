"""Execute a validated pipeline plan.

Planning is deliberately absent from this module.  The CLI and programmatic
callers must provide the same immutable ``PipelinePlan``, preventing the two
execution modes from constructing different dependency orders.
"""

from pathlib import Path
from typing import Callable, Mapping

from rich.console import Console

from postgwas.config.loader import canonical_module_name
from postgwas.core.checkpointing import (
    CheckpointAuditLogger,
    ExecutionCheckpoint,
    checkpoint_content_configuration,
    checkpoint_content_namespace,
    checkpoint_namespace,
    decode_checkpoint_value,
    discover_input_files,
    encode_checkpoint_value,
    preserve_checkpoint_run_controls,
    restore_checkpoint_namespace,
    software_identity,
)
from postgwas.core.contracts import RunContext
from postgwas.core.errors import ModuleExecutionError
from postgwas.core.input_validation import validation_scope
from postgwas.core.paths import configured_output_path
from postgwas.core.ui import PipelineStageController, StageProgress, print_screen_block, screen_line
from postgwas.pipeline.planner import PipelinePlan
from postgwas.pipeline.registry import REGISTRY, resolve_reference


console = Console()


def _stage_state(args, context: RunContext) -> dict:
    return {
        "arguments": checkpoint_namespace(args),
        "context": context.snapshot(),
    }


def _restore_stage_state(args, context: RunContext, document: Mapping) -> None:
    state = decode_checkpoint_value(document.get("state"))
    if (
        not isinstance(state, Mapping)
        or not isinstance(state.get("arguments"), Mapping)
        or not isinstance(state.get("context"), Mapping)
    ):
        raise ValueError("Pipeline checkpoint contains invalid post-stage state")
    run_controls = preserve_checkpoint_run_controls(args)
    restore_checkpoint_namespace(args, state["arguments"])
    for key, value in run_controls.items():
        setattr(args, key, value)
    context.clear()
    context.update(state["context"])


def _pipeline_checkpoint(
    *,
    args,
    plan: PipelinePlan,
    configuration,
    context: RunContext,
    number: int,
    module_name: str,
    output_name: str,
    output_root: Path,
    checkpoint_root: Path,
    checkpoint_logger: CheckpointAuditLogger,
    upstream_checkpoints: Mapping[str, Path],
) -> tuple[ExecutionCheckpoint, Path]:
    stage_number = "%02d" % number
    stage_root = output_root / ("%s_%s" % (stage_number, output_name))
    policy = configuration.run.resume_policy
    manifest = configured_output_path(
        checkpoint_root,
        policy.pipeline_stage_manifest,
        error_type=ValueError,
        stage_number=stage_number,
        module=module_name,
    )
    module_configuration = getattr(
        configuration.modules, canonical_module_name(module_name),
    )
    arguments_before = checkpoint_namespace(args)
    content_arguments = checkpoint_content_namespace(args)
    context_before = context.snapshot()
    checkpoint_configuration = {
        "resolved_configuration": checkpoint_content_configuration(configuration),
        "plan": {
            "requested_modules": plan.requested_modules,
            "active_modules": plan.active_modules,
            "steps": plan.steps,
        },
        "arguments_before": encode_checkpoint_value(content_arguments),
        "context_before": encode_checkpoint_value(context_before),
    }
    # A categorical CLI value can equal an output directory or one of its
    # parents. Prune the output subtree during implicit directory expansion while
    # retaining explicitly referenced upstream files below it.
    tracked_inputs = discover_input_files(
        arguments_before,
        context_before,
        module_configuration,
        configuration.resources,
        excluded_paths=(output_root,),
        excluded_roots=(stage_root, checkpoint_root),
    )
    return (
        ExecutionCheckpoint(
            manifest_path=manifest,
            output_root=output_root,
            artifact_root=stage_root,
            identity={
                "scope": "pipeline",
                "stage": "%s_%s" % (stage_number, module_name),
                "module": module_name,
                "stage_number": stage_number,
            },
            configuration=checkpoint_configuration,
            inputs=tracked_inputs,
            policy=policy,
            resume=args.resume,
            overwrite=args.overwrite,
            logger=checkpoint_logger,
            software=software_identity(
                arguments_before, module_configuration, configuration.resources,
            ),
            upstream_checkpoints=upstream_checkpoints,
            error_type=RuntimeError,
        ),
        manifest,
    )


def execute_pipeline(
    args,
    plan: PipelinePlan,
    configuration,
    *,
    preflight_evidence: Mapping[str, object] | None = None,
    finalize_validation: Callable[[], None] | None = None,
) -> RunContext:
    """Execute each planned step once and return its run context.

    A module name may occur more than once in a plan when it represents a real
    data transition, such as formatting before and after imputation.
    When supplied, final validation evidence must be published before the last
    stage can be reported complete, including when resuming that stage.
    """

    if not isinstance(plan, PipelinePlan):
        raise TypeError("execute_pipeline requires a validated PipelinePlan")

    output_root = Path(args.output_directory).expanduser().resolve()
    policy = configuration.run.resume_policy
    checkpoint_root = configured_output_path(
        output_root,
        str(policy.checkpoint_directory),
        error_type=ValueError,
    )
    audit_path = configured_output_path(
        checkpoint_root,
        policy.audit_log,
        error_type=ValueError,
    )
    checkpoint_logger = CheckpointAuditLogger(audit_path, console=console)
    context = RunContext(validations=preflight_evidence)
    validated_checkpoints: dict[str, Path] = {}
    print_screen_block("\n" + screen_line("run", "Starting execution chain"), console=console)
    progress_plan = None
    if len(plan.requested_modules) == 1:
        target = REGISTRY.require_pipeline_enabled(plan.requested_modules[0])
        factory = getattr(target, "pipeline_progress_factory", None)
        if factory is not None:
            progress_plan = resolve_reference(factory)(args)
    if progress_plan is None:
        total = len(plan.steps)
        progress = StageProgress(
            "Pipeline execution progress", enabled=True, console=console,
        )
        detailed_progress = None
    else:
        total = len(progress_plan["stages"])
        progress = None
        detailed_progress = PipelineStageController(
            progress_plan["label"],
            progress_plan["stages"],
            console=console,
            outcome_label_width=configuration.logging.terminal_label_width,
        )
        args._pipeline_stage_progress = detailed_progress
        args._pipeline_progress_plan = progress_plan

    try:
        for number, module_name in enumerate(plan.steps, 1):
            spec = REGISTRY.require_pipeline_enabled(module_name)
            title = spec.description
            title_factory = getattr(spec, "pipeline_title_factory", None)
            if title_factory is not None:
                title = resolve_reference(title_factory)(args)
            if progress is not None:
                progress.start_step(number, total, title)

            # Expose the planned step number to module output-directory builders.
            args._step_num = f"{number:02d}"
            output_name = spec.pipeline_output_name or module_name
            checkpoint = None
            runner_started = False
            try:
                checkpoint, manifest = _pipeline_checkpoint(
                    args=args,
                    plan=plan,
                    configuration=configuration,
                    context=context,
                    number=number,
                    module_name=module_name,
                    output_name=output_name,
                    output_root=output_root,
                    checkpoint_root=checkpoint_root,
                    checkpoint_logger=checkpoint_logger,
                    upstream_checkpoints=validated_checkpoints,
                )
                decision = checkpoint.prepare()
                if decision.action == "resume":
                    _restore_stage_state(args, context, decision.document)
                    args._step_num = f"{number:02d}"
                    if number == len(plan.steps) and finalize_validation is not None:
                        finalize_validation()
                    validated_checkpoints[
                        "%s_%s" % (args._step_num, module_name)
                    ] = manifest
                    if progress is not None:
                        progress.complete_step(
                            number, total, title,
                            outcome="Resumed checksum-validated checkpoint",
                        )
                    else:
                        start, end = progress_plan["modules"][module_name]
                        deferred_stage = (
                            progress_plan.get("deferred_completion_modules", {})
                            .get(module_name)
                        )
                        for stage in range(start, end + 1):
                            if detailed_progress.current != stage:
                                detailed_progress.start(stage)
                            if stage == deferred_stage:
                                break
                            detailed_progress.complete(
                                stage,
                                outcome="Resumed checksum-validated checkpoint",
                            )
                    continue
                runner_started = True
                runner = resolve_reference(spec.runner)
                with validation_scope(module_name):
                    result = runner(args, context)
                if detailed_progress is not None:
                    _start, expected_end = progress_plan["modules"][module_name]
                    if not (
                        detailed_progress.completed == expected_end
                        or detailed_progress.current == expected_end
                    ):
                        raise RuntimeError(
                            "%s completed without reporting every configured "
                            "pipeline stage through %d"
                            % (module_name, expected_end)
                        )
                state = _stage_state(args, context)
                checkpoint.write(
                    status="COMPLETED",
                    declared_values=(result, state["context"]),
                    state=state,
                    metrics={"result_type": type(result).__name__},
                )
                if number == len(plan.steps) and finalize_validation is not None:
                    finalize_validation()
                if (
                    detailed_progress is not None
                    and detailed_progress.current
                ):
                    deferred_stage = (
                        progress_plan.get("deferred_completion_modules", {})
                        .get(module_name)
                    )
                    if detailed_progress.current != deferred_stage:
                        completion = getattr(
                            args, "_pipeline_stage_completion", {},
                        )
                        completion_callback = getattr(
                            args, "_pipeline_progress_completion", None,
                        )
                        if completion_callback is not None:
                            completion_callback()
                        detailed_progress.complete(
                            detailed_progress.current,
                            outcome=completion.get("outcome"),
                            outcome_fields=completion.get("outcome_fields"),
                        )
                        if hasattr(args, "_pipeline_stage_completion"):
                            delattr(args, "_pipeline_stage_completion")
            except BaseException as exc:
                if (
                    runner_started
                    and checkpoint is not None
                    and checkpoint.should_record_partial
                ):
                    try:
                        checkpoint.write(
                            status="PARTIAL",
                            declared_values=context.snapshot(),
                            state=_stage_state(args, context),
                            metrics={
                                "error_type": type(exc).__name__,
                                "error": str(exc),
                            },
                        )
                    except Exception as checkpoint_exc:
                        checkpoint_logger.record(
                            "FAILED",
                            "%s_%s" % (args._step_num, module_name),
                            action="record_partial_checkpoint",
                            error="%s: %s"
                            % (type(checkpoint_exc).__name__, checkpoint_exc),
                        )
                if progress is not None:
                    progress.fail_step(number, total, title)
                else:
                    failure_callback = getattr(
                        args, "_pipeline_progress_failure", None,
                    )
                    if failure_callback is not None:
                        try:
                            failure_callback(exc)
                        except Exception as logging_exc:
                            checkpoint_logger.record(
                                "FAILED",
                                "%s_%s" % (args._step_num, module_name),
                                action="record_detailed_progress_failure",
                                error="%s: %s"
                                % (type(logging_exc).__name__, logging_exc),
                            )
                    detailed_progress.fail_active()
                print_screen_block(screen_line("error", f"Pipeline Failed at step: {module_name}"), console=console)
                if isinstance(exc, (KeyboardInterrupt, SystemExit)):
                    raise
                raise ModuleExecutionError(module_name, str(exc)) from exc
            validated_checkpoints[
                "%s_%s" % (args._step_num, module_name)
            ] = manifest
            if progress is not None:
                progress.complete_step(number, total, title)
    finally:
        if progress is not None:
            progress.close()
        else:
            detailed_progress.close()
        if hasattr(args, "_pipeline_progress_failure"):
            delattr(args, "_pipeline_progress_failure")
        if hasattr(args, "_pipeline_progress_completion"):
            delattr(args, "_pipeline_progress_completion")
        if hasattr(args, "_pipeline_stage_completion"):
            delattr(args, "_pipeline_stage_completion")

    print_screen_block("\n" + screen_line("success", "All tasks completed successfully.") + "\n", console=console)
    return context
