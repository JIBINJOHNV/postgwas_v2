from pathlib import Path

import pytest

from postgwas.core.completion import (
    apply_completion_restart,
    configuration_digest,
    resolve_completion_resume,
    write_completion_manifest,
)


class RestartPolicy:
    changed_parameters = "warn_and_restart"
    changed_inputs = "warn_and_restart"
    unvalidated_outputs = "warn_and_restart"


class RecordingLogger:
    def __init__(self):
        self.warnings = []
        self.records = []

    def warning(self, message):
        self.warnings.append(message)

    def record(self, marker, subject, **values):
        self.records.append((marker, subject, values))


def _completed_run(tmp_path: Path):
    source = tmp_path / "input.tsv"
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    source.write_text("input\n", encoding="utf-8")
    first.write_text("first\n", encoding="utf-8")
    second.write_text("second\n", encoding="utf-8")
    manifest = tmp_path / "completion.yaml"
    inputs = {"source": source}
    outputs = {"first": first, "second": second}
    digest = configuration_digest({"setting": "value"})
    write_completion_manifest(
        manifest,
        dataset_id="STUDY",
        module="test_module",
        genome_build="GRCh37",
        configuration_sha256=digest,
        inputs=inputs,
        outputs=outputs,
        metrics={"rows": 2},
    )
    return manifest, inputs, outputs, digest


def _resolve(manifest, inputs, outputs, digest, *, policy=None):
    return resolve_completion_resume(
        manifest,
        dataset_id="STUDY",
        module="test_module",
        genome_build="GRCh37",
        configuration_sha256=digest,
        inputs=inputs,
        outputs=outputs,
        resume_policy=policy,
    )


def test_completion_resume_reuses_only_matching_complete_outputs(tmp_path):
    manifest, inputs, outputs, digest = _completed_run(tmp_path)

    decision = _resolve(manifest, inputs, outputs, digest)

    assert decision.action == "resume"
    assert decision.manifest["metrics"] == {"rows": 2}


def test_completion_resume_restarts_when_every_output_is_absent_by_policy(tmp_path):
    manifest, inputs, outputs, digest = _completed_run(tmp_path)
    for output in outputs.values():
        output.unlink()

    decision = _resolve(
        manifest, inputs, outputs, digest, policy=RestartPolicy(),
    )

    assert decision.action == "restart"
    assert decision.reason == "incomplete_outputs"
    assert decision.warning is not None
    assert decision.missing_outputs == ("first", "second")
    assert manifest.is_file()


def test_completion_resume_restarts_a_partial_output_set_by_policy(tmp_path):
    manifest, inputs, outputs, digest = _completed_run(tmp_path)
    outputs["second"].unlink()

    decision = _resolve(
        manifest, inputs, outputs, digest, policy=RestartPolicy(),
    )

    assert decision.action == "restart"
    assert decision.reason == "incomplete_outputs"
    assert decision.missing_outputs == ("second",)


def test_completion_resume_refuses_changed_output_content_even_by_policy(tmp_path):
    manifest, inputs, outputs, digest = _completed_run(tmp_path)
    outputs["first"].write_text("other\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="Refusing automatic restart"):
        _resolve(manifest, inputs, outputs, digest, policy=RestartPolicy())


def test_completion_restart_accepts_changed_input_by_policy(tmp_path):
    manifest, inputs, outputs, digest = _completed_run(tmp_path)
    for output in outputs.values():
        output.unlink()
    inputs["source"].write_text("changed input\n", encoding="utf-8")

    decision = _resolve(
        manifest, inputs, outputs, digest, policy=RestartPolicy(),
    )

    assert decision.action == "restart"
    assert decision.reason == "changed_inputs"


def test_completion_restart_accepts_changed_configuration_by_policy(tmp_path):
    manifest, inputs, outputs, digest = _completed_run(tmp_path)
    for output in outputs.values():
        output.unlink()

    decision = _resolve(
        manifest,
        inputs,
        outputs,
        configuration_digest({"setting": "changed"}),
        policy=RestartPolicy(),
    )

    assert decision.action == "restart"
    assert decision.reason == "changed_parameters"


def test_completion_restart_accepts_changed_output_path_by_policy(tmp_path):
    manifest, inputs, outputs, digest = _completed_run(tmp_path)
    for output in outputs.values():
        output.unlink()
    outputs["first"] = tmp_path / "renamed.tsv"

    decision = _resolve(
        manifest, inputs, outputs, digest, policy=RestartPolicy(),
    )

    assert decision.action == "restart"
    assert decision.reason == "changed_output_contract"


def test_completion_changes_remain_errors_when_policy_requests_error(tmp_path):
    manifest, inputs, outputs, digest = _completed_run(tmp_path)
    inputs["source"].write_text("changed input\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="Input source changed"):
        _resolve(manifest, inputs, outputs, digest)


def test_completion_restart_removes_only_verified_owned_outputs(tmp_path):
    manifest, inputs, outputs, digest = _completed_run(tmp_path)
    decision = _resolve(
        manifest,
        inputs,
        outputs,
        configuration_digest({"setting": "changed"}),
        policy=RestartPolicy(),
    )
    logger = RecordingLogger()

    removed = apply_completion_restart(
        decision,
        output_root=tmp_path,
        manifest=manifest,
        logger=logger,
        operation="test_restart",
    )

    assert set(removed) == {path.resolve() for path in outputs.values()}
    assert not manifest.exists()
    assert all(not path.exists() for path in outputs.values())
    assert logger.warnings == [decision.warning]
    assert logger.records[0][2]["reason"] == "changed_parameters"


def test_completion_restart_removes_verified_survivors_of_partial_set(tmp_path):
    manifest, inputs, outputs, digest = _completed_run(tmp_path)
    outputs["second"].unlink()
    decision = _resolve(
        manifest, inputs, outputs, digest, policy=RestartPolicy(),
    )

    removed = apply_completion_restart(
        decision,
        output_root=tmp_path,
        manifest=manifest,
        logger=RecordingLogger(),
        operation="test_restart",
    )

    assert removed == (outputs["first"].resolve(),)
    assert not outputs["first"].exists()
    assert not manifest.exists()


def test_completion_restart_refuses_modified_output(tmp_path):
    manifest, inputs, outputs, digest = _completed_run(tmp_path)
    decision = _resolve(
        manifest,
        inputs,
        outputs,
        configuration_digest({"setting": "changed"}),
        policy=RestartPolicy(),
    )
    outputs["first"].write_text("externally modified\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="was modified after its checkpoint"):
        apply_completion_restart(
            decision,
            output_root=tmp_path,
            manifest=manifest,
            logger=RecordingLogger(),
            operation="test_restart",
        )

    assert manifest.is_file()
    assert all(path.is_file() for path in outputs.values())


def test_completion_restart_refuses_recorded_output_outside_root(tmp_path):
    manifest, inputs, outputs, digest = _completed_run(tmp_path)
    decision = _resolve(
        manifest,
        inputs,
        outputs,
        configuration_digest({"setting": "changed"}),
        policy=RestartPolicy(),
    )
    outside = tmp_path.parent / "outside.tsv"
    tampered = dict(decision.manifest)
    tampered_outputs = dict(tampered["outputs"])
    tampered_outputs["first"] = {
        **tampered_outputs["first"],
        "path": str(outside),
    }
    tampered["outputs"] = tampered_outputs
    unsafe = type(decision)(
        action=decision.action,
        manifest=tampered,
        reason=decision.reason,
        warning=decision.warning,
    )

    with pytest.raises(RuntimeError, match="outside its configured output root"):
        apply_completion_restart(
            unsafe,
            output_root=tmp_path,
            manifest=manifest,
            logger=RecordingLogger(),
            operation="test_restart",
        )

    assert manifest.is_file()
    assert all(path.is_file() for path in outputs.values())
