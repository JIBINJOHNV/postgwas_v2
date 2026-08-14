"""Regression tests for global direct and pipeline execution checkpoints."""

from argparse import Namespace
from io import StringIO
from pathlib import Path

import pytest
from rich.console import Console

from postgwas.config import load_configuration
from postgwas.core.checkpointing import (
    CheckpointAuditLogger,
    ExecutionCheckpoint,
    checkpoint_namespace,
    decode_checkpoint_value,
    discover_input_files,
    encode_checkpoint_value,
    restore_checkpoint_namespace,
)
from postgwas.core.contracts import Artifact, ModuleResult


def _logger(root: Path) -> tuple[CheckpointAuditLogger, StringIO]:
    stream = StringIO()
    return (
        CheckpointAuditLogger(
            root / "checkpoint.log",
            console=Console(file=stream, force_terminal=False),
        ),
        stream,
    )


def _manager(
    root: Path,
    *,
    configuration=None,
    inputs=None,
    upstream=None,
    logger=None,
    software=None,
    overwrite=False,
) -> ExecutionCheckpoint:
    policy = load_configuration().run.resume_policy
    return ExecutionCheckpoint(
        manifest_path=root / "run_metadata" / "checkpoints" / "01_stage.yaml",
        output_root=root,
        artifact_root=root / "01_stage",
        identity={"scope": "pipeline", "stage": "01_stage"},
        configuration=configuration or {"threshold": 0.05},
        inputs=inputs or {},
        policy=policy,
        resume=True,
        overwrite=overwrite,
        logger=logger or _logger(root)[0],
        software=software or {"postgwas": "test"},
        upstream_checkpoints=upstream or {},
    )


def _completed_checkpoint(
    root: Path,
    *,
    configuration=None,
    inputs=None,
    upstream=None,
    software=None,
):
    artifact_root = root / "01_stage"
    artifact_root.mkdir(parents=True)
    manager = _manager(
        root,
        configuration=configuration,
        inputs=inputs,
        upstream=upstream,
        software=software,
    )
    output = artifact_root / "result.tsv"
    output.write_text("GENE\tP\nGENE1\t0.01\n", encoding="utf-8")
    state = {
        "args": {"dataset_id": "STUDY", "output_directory": root},
        "context": {
            "stage": ModuleResult(
                "stage",
                artifacts={"result": Artifact("table", output)},
                metrics={"rows": 1},
            )
        },
    }
    manager.write(
        status="COMPLETED",
        declared_values={"result": output},
        state=state,
        metrics={"rows": 1},
    )
    return output, state


def test_completed_checkpoint_resumes_and_restores_safe_state(tmp_path):
    output, expected_state = _completed_checkpoint(tmp_path)

    decision = _manager(tmp_path).prepare()

    assert decision.action == "resume"
    assert output.is_file()
    restored = decode_checkpoint_value(decision.document["state"])
    assert restored == expected_state
    assert restored["context"]["stage"].metrics["rows"] == 1


def test_real_partial_checkpoint_restarts_at_stage_boundary(tmp_path):
    artifact_root = tmp_path / "01_stage"
    artifact_root.mkdir()
    logger, stream = _logger(tmp_path)
    manager = _manager(tmp_path, logger=logger)
    partial = artifact_root / "partial.tsv"
    partial.write_text("partial\n", encoding="utf-8")
    manager.write(status="PARTIAL", declared_values={"partial": partial})

    decision = _manager(tmp_path, logger=logger).prepare()

    assert decision.action == "run"
    assert not partial.exists()
    assert not manager.manifest_path.exists()
    assert "real partial result exists" in stream.getvalue()
    assert "partial_results" in logger.path.read_text(encoding="utf-8")


@pytest.mark.parametrize(
    ("changed", "reason"),
    (("configuration", "changed_parameters"), ("input", "changed_inputs")),
)
def test_changed_parameters_or_inputs_warn_and_restart(tmp_path, changed, reason):
    input_file = tmp_path / "input.tsv"
    input_file.write_text("old\n", encoding="utf-8")
    output, _ = _completed_checkpoint(
        tmp_path,
        inputs={"study": input_file},
    )
    logger, stream = _logger(tmp_path)
    if changed == "configuration":
        manager = _manager(
            tmp_path,
            configuration={"threshold": 0.01},
            inputs={"study": input_file},
            logger=logger,
        )
    else:
        input_file.write_text("new\n", encoding="utf-8")
        manager = _manager(
            tmp_path,
            inputs={"study": input_file},
            logger=logger,
        )

    decision = manager.prepare()

    assert decision.action == "run"
    assert not output.exists()
    assert "restart" in stream.getvalue()
    assert "safe boundary" in stream.getvalue()
    assert reason in logger.path.read_text(encoding="utf-8")


def test_modified_checkpoint_output_is_never_deleted_automatically(tmp_path):
    output, _ = _completed_checkpoint(tmp_path)
    output.write_text("user replacement\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="modified after validation"):
        _manager(tmp_path).prepare()

    assert output.read_text(encoding="utf-8") == "user replacement\n"


def test_missing_current_input_preserves_previous_validated_outputs(tmp_path):
    input_file = tmp_path / "input.tsv"
    input_file.write_text("old\n", encoding="utf-8")
    output, _ = _completed_checkpoint(
        tmp_path,
        inputs={"study": input_file},
    )
    input_file.unlink()

    with pytest.raises(RuntimeError, match="input is missing"):
        _manager(tmp_path).prepare()

    assert output.is_file()


def test_changed_upstream_checkpoint_invalidates_downstream_stage(tmp_path):
    upstream = tmp_path / "run_metadata" / "checkpoints" / "00_upstream.yaml"
    upstream.parent.mkdir(parents=True)
    upstream.write_text("status: COMPLETED\n", encoding="utf-8")
    output, _ = _completed_checkpoint(
        tmp_path,
        upstream={"00_upstream": upstream},
    )
    upstream.write_text("status: RECOMPUTED\n", encoding="utf-8")

    decision = _manager(
        tmp_path,
        upstream={"00_upstream": upstream},
    ).prepare()

    assert decision.action == "run"
    assert not output.exists()


def test_namespace_and_contract_state_round_trip_without_pickle(tmp_path):
    namespace = Namespace(dataset_id="STUDY", path=tmp_path / "input.tsv")
    namespace._step_num = "01"
    encoded = encode_checkpoint_value(checkpoint_namespace(namespace))
    restored = decode_checkpoint_value(encoded)
    target = Namespace(old="remove", _private="keep")

    restore_checkpoint_namespace(target, restored)

    assert vars(target) == {
        "_private": "keep",
        "dataset_id": "STUDY",
        "path": tmp_path / "input.tsv",
    }


def test_input_discovery_tracks_files_directories_and_reference_prefixes(tmp_path):
    explicit = tmp_path / "study.tsv"
    explicit.write_text("study\n", encoding="utf-8")
    directory = tmp_path / "reference"
    directory.mkdir()
    member = directory / "scores.tsv"
    member.write_text("scores\n", encoding="utf-8")
    prefix = tmp_path / "panel"
    bim = tmp_path / "panel.bim"
    bim.write_text("1 rs1 0 1 A G\n", encoding="utf-8")
    ignored = tmp_path / "results"
    ignored.mkdir()
    (ignored / "old.tsv").write_text("old\n", encoding="utf-8")

    discovered = discover_input_files(
        {
            "input": explicit,
            "reference_directory": directory,
            "reference_prefix": prefix,
            "output_directory": ignored,
        },
        excluded_roots=(ignored,),
    )

    assert set(discovered.values()) == {
        explicit.resolve(), member.resolve(), bim.resolve(),
    }
