"""Global checkpoint boundary for direct scientific commands."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Sequence

from rich.console import Console

from postgwas.config.loader import canonical_module_name
from postgwas.core.checkpointing import (
    CheckpointAuditLogger,
    ExecutionCheckpoint,
    NoOwnedCheckpointOutputs,
    checkpoint_content_arguments,
    checkpoint_content_configuration,
    discover_input_files,
    software_identity,
)
from postgwas.core.paths import configured_output_path, validate_filename_component


def _runtime_controls(configuration, arguments: Sequence[str]) -> tuple[bool, bool]:
    resume = bool(configuration.run.resume)
    overwrite = bool(configuration.run.overwrite)
    for argument in arguments:
        if argument == "--resume":
            resume = True
        elif argument == "--no-resume":
            resume = False
        elif argument == "--overwrite":
            overwrite = True
    return resume, overwrite


def _checkpoint_configuration(
    configuration,
    module_name: str,
    public_command: str,
    arguments: Sequence[str],
) -> dict[str, Any]:
    canonical = canonical_module_name(module_name)
    resolved = checkpoint_content_configuration(configuration)
    return {
        "command": public_command,
        "arguments": checkpoint_content_arguments(arguments),
        "run": resolved["run"],
        "execution": resolved["execution"],
        "logging": resolved["logging"],
        "resources": resolved["resources"],
        "module": getattr(configuration.modules, canonical).model_dump(mode="json"),
    }


def _record_partial(
    checkpoint: ExecutionCheckpoint,
    logger: CheckpointAuditLogger,
    command: str,
    error: BaseException | None,
) -> None:
    if not checkpoint.should_record_partial:
        return
    try:
        checkpoint.write(
            status="PARTIAL",
            metrics={
                "error_type": type(error).__name__ if error is not None else "exit",
                "error": str(error) if error is not None else "non-zero return code",
            },
        )
    except Exception as checkpoint_error:
        logger.record(
            "FAILED",
            command,
            action="record_partial_checkpoint",
            error="%s: %s"
            % (type(checkpoint_error).__name__, checkpoint_error),
        )


def run_direct_with_checkpoint(
    operation: Callable[[], Any],
    *,
    module_name: str,
    public_command: str,
    arguments: Sequence[str],
    configuration,
    output_directory: str | Path,
    console: Console | None = None,
) -> Any:
    """Run or safely resume one complete direct-command scientific boundary."""
    validate_filename_component(public_command, "direct command")
    output_root = Path(output_directory).expanduser().resolve()
    policy = configuration.run.resume_policy
    checkpoint_root = configured_output_path(
        output_root,
        str(policy.checkpoint_directory),
        error_type=ValueError,
    )
    manifest = configured_output_path(
        checkpoint_root,
        policy.direct_manifest,
        error_type=ValueError,
        command=public_command,
    )
    audit_path = configured_output_path(
        checkpoint_root,
        policy.audit_log,
        error_type=ValueError,
    )
    logger = CheckpointAuditLogger(audit_path, console=console)
    canonical = canonical_module_name(module_name)
    module_configuration = getattr(configuration.modules, canonical)
    inputs = discover_input_files(
        list(arguments),
        module_configuration,
        configuration.resources,
        excluded_paths=(output_root,),
        excluded_roots=(checkpoint_root,),
    )
    resume, overwrite = _runtime_controls(configuration, arguments)
    checkpoint = ExecutionCheckpoint(
        manifest_path=manifest,
        output_root=output_root,
        artifact_root=output_root,
        identity={
            "scope": "direct",
            "stage": "direct_%s" % public_command,
            "module": canonical,
            "command": public_command,
        },
        configuration=_checkpoint_configuration(
            configuration, module_name, public_command, arguments,
        ),
        inputs=inputs,
        policy=policy,
        resume=resume,
        overwrite=overwrite,
        logger=logger,
        software=software_identity(
            list(arguments), module_configuration, configuration.resources,
        ),
        error_type=RuntimeError,
    )

    decision = checkpoint.prepare()
    if decision.action == "resume":
        return 0

    try:
        result = operation()
    except SystemExit as exc:
        if exc.code in (None, 0):
            try:
                checkpoint.write(status="COMPLETED", metrics={"exit_code": 0})
            except NoOwnedCheckpointOutputs as checkpoint_error:
                logger.record("SKIP", public_command, reason=str(checkpoint_error))
        else:
            _record_partial(checkpoint, logger, public_command, exc)
        raise
    except BaseException as exc:
        _record_partial(checkpoint, logger, public_command, exc)
        raise

    failed = (
        isinstance(result, int)
        and not isinstance(result, bool)
        and result != 0
    )
    if failed:
        _record_partial(checkpoint, logger, public_command, None)
        return result
    try:
        checkpoint.write(
            status="COMPLETED",
            metrics={"return_type": type(result).__name__},
        )
    except NoOwnedCheckpointOutputs as checkpoint_error:
        logger.record("SKIP", public_command, reason=str(checkpoint_error))
    return result


__all__ = ["run_direct_with_checkpoint"]
