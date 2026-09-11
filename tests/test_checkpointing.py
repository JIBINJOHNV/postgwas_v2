"""Regression tests for global direct and pipeline execution checkpoints."""

from argparse import Namespace
import errno
from io import StringIO
from pathlib import Path

import pytest
from rich.console import Console
import yaml

from postgwas.config import load_configuration
from postgwas.core import checkpointing
from postgwas.core.checkpointing import (
    CheckpointAuditLogger,
    ExecutionCheckpoint,
    checkpoint_namespace,
    decode_checkpoint_value,
    discover_input_files,
    encode_checkpoint_value,
    restore_checkpoint_namespace,
    software_identity,
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
    resume=True,
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
        resume=resume,
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


def test_unchanged_input_reuses_stat_validated_sha256(monkeypatch, tmp_path):
    input_file = tmp_path / "input.tsv"
    input_file.write_text("stable\n", encoding="utf-8")
    output, _ = _completed_checkpoint(
        tmp_path,
        inputs={"study": input_file},
    )

    def unexpected_rehash(*_args, **_kwargs):
        raise AssertionError("unchanged input was rehashed")

    monkeypatch.setattr(checkpointing, "file_fingerprint", unexpected_rehash)

    decision = _manager(
        tmp_path,
        inputs={"study": input_file},
    ).prepare()

    assert decision.action == "resume"
    assert output.is_file()
    assert decision.document["input_fingerprint_cache"]["study"]["sha256"]


@pytest.mark.parametrize(
    ("resume", "overwrite"),
    ((False, False), (True, True)),
)
def test_execution_controls_reuse_unchanged_input_hash(
    monkeypatch, tmp_path, resume, overwrite,
):
    input_file = tmp_path / "input.tsv"
    input_file.write_text("stable\n", encoding="utf-8")
    output, _ = _completed_checkpoint(
        tmp_path,
        inputs={"study": input_file},
    )

    def unexpected_rehash(*_args, **_kwargs):
        raise AssertionError("unchanged input was rehashed")

    monkeypatch.setattr(checkpointing, "file_fingerprint", unexpected_rehash)
    manager = _manager(
        tmp_path,
        inputs={"study": input_file},
        resume=resume,
        overwrite=overwrite,
    )

    assert manager.prepare().action == "run"
    output.write_text("recomputed\n", encoding="utf-8")
    manager.write(status="COMPLETED", declared_values={"result": output})

    assert manager.manifest_path.is_file()


def test_same_size_input_change_invalidates_cache_and_rehashes(monkeypatch, tmp_path):
    input_file = tmp_path / "input.tsv"
    input_file.write_text("old\n", encoding="utf-8")
    output, _ = _completed_checkpoint(
        tmp_path,
        inputs={"study": input_file},
    )
    input_file.write_text("new\n", encoding="utf-8")
    original = checkpointing.file_fingerprint
    rehashed = []

    def record_rehash(path, **kwargs):
        rehashed.append(Path(path).resolve())
        return original(path, **kwargs)

    monkeypatch.setattr(checkpointing, "file_fingerprint", record_rehash)

    decision = _manager(
        tmp_path,
        inputs={"study": input_file},
    ).prepare()

    assert decision.action == "run"
    assert rehashed == [input_file.resolve()]
    assert not output.exists()


def test_legacy_checkpoint_hash_seeds_stat_validated_cache(monkeypatch, tmp_path):
    input_file = tmp_path / "input.tsv"
    input_file.write_text("stable\n", encoding="utf-8")
    output, _ = _completed_checkpoint(
        tmp_path,
        inputs={"study": input_file},
    )
    manifest = tmp_path / "run_metadata/checkpoints/01_stage.yaml"
    document = yaml.safe_load(manifest.read_text(encoding="utf-8"))
    document.pop("input_fingerprint_cache")
    manifest.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")

    def unexpected_rehash(*_args, **_kwargs):
        raise AssertionError("legacy checkpoint input was rehashed")

    monkeypatch.setattr(checkpointing, "file_fingerprint", unexpected_rehash)

    decision = _manager(
        tmp_path,
        inputs={"study": input_file},
    ).prepare()

    assert decision.action == "resume"
    assert output.is_file()


def test_input_change_during_execution_refuses_checkpoint(tmp_path):
    input_file = tmp_path / "input.tsv"
    input_file.write_text("old\n", encoding="utf-8")
    artifact_root = tmp_path / "01_stage"
    artifact_root.mkdir()
    manager = _manager(tmp_path, inputs={"study": input_file})
    input_file.write_text("new\n", encoding="utf-8")
    output = artifact_root / "result.tsv"
    output.write_text("result\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="changed while the execution was active"):
        manager.write(status="COMPLETED", declared_values={"result": output})

    assert not manager.manifest_path.exists()


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
    generated = directory / "generated"
    generated.mkdir()
    generated_member = generated / "runtime.log"
    generated_member.write_text("runtime\n", encoding="utf-8")
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
        excluded_roots=(ignored, generated),
    )

    assert set(discovered.values()) == {
        explicit.resolve(), member.resolve(), bim.resolve(),
    }


def test_input_discovery_prunes_output_reached_through_parent_directory(
    tmp_path, monkeypatch,
):
    study_root = tmp_path / "SCZ_2026"
    module_directory = study_root / "gcta_gene"
    output = module_directory / "mbat_combo"
    screen_log = output / "run_metadata/screen.log"
    screen_log.parent.mkdir(parents=True)
    screen_log.write_text("live transcript\n", encoding="utf-8")
    explicit_input = output / "upstream.ma"
    explicit_input.write_text("SNP A1 A2 freq BETA SE P N\n", encoding="utf-8")
    monkeypatch.chdir(study_root)

    discovered = discover_input_files(
        {
            "modules": ["formatter", "gcta_gene"],
            "method": "mbat_combo",
            "gcta_input_file": explicit_input,
        },
        excluded_paths=(output,),
    )

    assert explicit_input.resolve() in discovered.values()
    assert screen_log.resolve() not in discovered.values()


def test_checkpoint_path_scanners_ignore_overlong_non_path_values(tmp_path):
    output = tmp_path / "01_stage/result.tsv"
    output.parent.mkdir()
    output.write_text("result\n", encoding="utf-8")
    shell_metadata = "rs=0:" + "di=01;34:" * 100

    assert discover_input_files(
        {"environment": {"LS_COLORS": shell_metadata}}
    ) == {}
    assert software_identity(
        {"environment": {"LS_COLORS": shell_metadata}}
    )["executables"] == {}

    manager = _manager(tmp_path)
    manifest = manager.write(
        status="COMPLETED",
        declared_values={"result": output, "diagnostic": shell_metadata},
    )

    assert manifest.is_file()


def test_input_discovery_preserves_permission_failures(tmp_path, monkeypatch):
    source = tmp_path / "reference.tsv"

    def denied(_path):
        raise PermissionError(errno.EACCES, "Permission denied")

    monkeypatch.setattr(Path, "is_symlink", denied)
    with pytest.raises(PermissionError):
        discover_input_files({"reference": source})
