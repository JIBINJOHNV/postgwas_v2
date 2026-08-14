"""Execute a validated pipeline plan.

Planning is deliberately absent from this module.  The CLI and programmatic
callers must provide the same immutable ``PipelinePlan``, preventing the two
execution modes from constructing different dependency orders.
"""

from pathlib import Path
from typing import Mapping

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
from postgwas.core.paths import configured_output_path
from postgwas.core.ui import StageProgress
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
    tracked_inputs = discover_input_files(
        arguments_before,
        context_before,
        module_configuration,
        configuration.resources,
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


def execute_pipeline(args, plan: PipelinePlan, configuration) -> RunContext:
    """Execute each planned step once and return its run context.

    A module name may occur more than once in a plan when it represents a real
    data transition, such as formatting before and after imputation.
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
    context = RunContext()
    validated_checkpoints: dict[str, Path] = {}
    console.print("\n🚀 [bold green]Starting Execution Chain[/bold green]")
    total = len(plan.steps)
    progress = StageProgress(
        "Pipeline execution progress",
        enabled=True,
        console=console,
    )

    try:
        for number, module_name in enumerate(plan.steps, 1):
            spec = REGISTRY.require_pipeline_enabled(module_name)
            title = spec.description
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
                    validated_checkpoints[
                        "%s_%s" % (args._step_num, module_name)
                    ] = manifest
                    progress.complete_step(
                        number, total, title,
                        outcome="Resumed checksum-validated checkpoint",
                    )
                    continue
                runner_started = True
                runner = resolve_reference(spec.runner)
                result = runner(args, context)
                state = _stage_state(args, context)
                checkpoint.write(
                    status="COMPLETED",
                    declared_values=(result, state["context"]),
                    state=state,
                    metrics={"result_type": type(result).__name__},
                )
            except Exception as exc:
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
                progress.fail_step(number, total, title)
                console.print(
                    f"\n❌ [bold red]Pipeline Failed at step: {module_name}[/bold red]"
                )
                raise ModuleExecutionError(module_name, str(exc)) from exc
            validated_checkpoints[
                "%s_%s" % (args._step_num, module_name)
            ] = manifest
            progress.complete_step(number, total, title)
    finally:
        progress.close()

    console.print("\n✅ [bold green]All tasks completed successfully.[/bold green]\n")
    return context
