"""Native informational commands must establish identity, not publish errors."""

import shutil
import subprocess
from types import SimpleNamespace

import pytest

from postgwas.modules.fine_mapping import preflight


def probe():
    return preflight._tool_version(
        "configured-tool", "Test executable", 7,
        arguments=("-help",), version_pattern=r"(?m)^Welcome to (Test v[0-9.]+)$",
    )


@pytest.mark.parametrize("stream", ["stdout", "stderr"])
def test_version_is_selected_from_both_streams_not_the_first_banner_line(monkeypatch, stream):
    result = subprocess.CompletedProcess([], 0, stdout="Usage text\n", stderr="")
    setattr(result, stream, getattr(result, stream) + "Banner\nWelcome to Test v1.2\n")

    def run(command, **kwargs):
        assert command == ["configured-tool", "-help"]
        assert kwargs == dict(check=False, capture_output=True, text=True, timeout=7.0)
        return result

    monkeypatch.setattr(preflight.subprocess, "run", run)
    assert probe() == "Test v1.2"


@pytest.mark.parametrize("returncode,output", [
    (1, "Welcome to Test v1.2"),
    (0, "Error : unsupported option\nWelcome to Test v1.2"),
    (0, "!! Error (Option): unsupported option\nWelcome to Test v1.2"),
    (0, ""),
    (0, "Banner without a version"),
    (0, "Welcome to Other v1.2"),
    (0, "Welcome to Test v1.2\nWelcome to Test v1.2"),
    (0, "Welcome to Test v1.2\nWelcome to Test v2.0"),
])
def test_failed_empty_wrong_or_ambiguous_version_output_is_rejected(monkeypatch, returncode, output):
    monkeypatch.setattr(preflight.subprocess, "run", lambda *_args, **_kwargs:
                        subprocess.CompletedProcess([], returncode, stdout="", stderr=output))
    with pytest.raises(preflight.FineMappingPreflightError, match="Test executable runtime version probe"):
        probe()


@pytest.mark.parametrize("error", [OSError("not executable"), subprocess.TimeoutExpired("tool", 7)])
def test_unavailable_or_timed_out_version_probe_is_actionable(monkeypatch, error):
    def run(*_args, **_kwargs):
        raise error
    monkeypatch.setattr(preflight.subprocess, "run", run)
    with pytest.raises(preflight.FineMappingPreflightError, match="Cannot execute Test executable"):
        probe()


def test_native_finemap_stack_reports_actual_versions():
    paths = {name: shutil.which(name) for name in ("plink2", "bgenix", "ldstore", "finemap")}
    if not all(paths.values()):
        pytest.skip("native FINEMAP runtime stack is unavailable")
    args = SimpleNamespace(
        finemap_method="finemap", plink=paths["plink2"], bgenix=paths["bgenix"],
        ldstore=paths["ldstore"], finemap_executable=paths["finemap"],
        fine_mapping_runtime_defaults={"tool_version_timeout_seconds": 20},
    )
    evidence, runtime = preflight._validate_runtime_tools(args)
    versions = {item.name: item.version for item in evidence}
    assert runtime is None
    assert versions["PLINK"].startswith("PLINK v2")
    assert versions["BGENIX"][0].isdigit()
    assert versions["LDstore"].startswith("LDstore v")
    assert versions["FINEMAP"].startswith("FINEMAP v")


@pytest.mark.skipif(shutil.which("plink") is None, reason="native PLINK 1 is unavailable")
def test_native_susie_plink_protocol():
    assert preflight._tool_version(
        shutil.which("plink"), "PLINK executable", 20,
        arguments=("--version",), version_pattern=r"(?m)^(PLINK v[0-9][^\r\n]+)$",
    ).startswith("PLINK v1")


@pytest.mark.parametrize("cached", [False, True])
def test_finemap_adapter_reuses_shared_runtime_validation_before_analysis(tmp_path, monkeypatch, cached):
    from postgwas.modules.fine_mapping.engines.finemap import adapter

    evidence = tuple(preflight.FineMappingToolEvidence(name, "configured-tool", version)
                     for name, version in [("PLINK", "PLINK v2 test"), ("BGENIX", "1.2"),
                                           ("LDstore", "LDstore v2"), ("FINEMAP", "FINEMAP v1.4")])
    calls = []
    monkeypatch.setattr(adapter, "resolve_compute_args", lambda _args: None)
    monkeypatch.setattr(adapter, "_validate_runtime_tools", lambda _args:
                        (calls.append("runtime") or evidence, None))
    monkeypatch.setattr(adapter, "setup_directories", lambda *_args: {"finemap_input_qc_file": tmp_path / "qc"})

    def stop_before_analysis(*_args, **_kwargs):
        raise RuntimeError("validated runtime reached input boundary")

    monkeypatch.setattr(adapter, "load_and_prep_inputs", stop_before_analysis)
    args = SimpleNamespace(
        output_directory=tmp_path, dataset_id="test", genome_build="GRCh37", plink="configured-plink",
        prob_cred_set=0.95, fine_mapping_output_layout={},
        fine_mapping_runtime_defaults={"schema_inference_length": 2, "max_maf": 0.5},
    )
    if cached:
        args._fine_mapping_preflight = SimpleNamespace(tools=evidence)
    with pytest.raises(RuntimeError, match="validated runtime reached input boundary"):
        adapter._run_finemap_pipeline(args, SimpleNamespace(record=lambda *_args: None))
    assert calls == ([] if cached else ["runtime"])
