"""Native R library isolation and shared namespace-preflight regressions."""

import os
from pathlib import Path
import shutil
import subprocess

import pytest

from postgwas.core.r_runtime import resolve_r_runtime


@pytest.fixture
def rscript():
    executable = shutil.which("Rscript")
    if not executable:
        pytest.skip("Native R runtime is unavailable")
    return executable


def test_native_installation_library_precedes_broken_user_libraries(rscript, tmp_path, monkeypatch):
    library = tmp_path / "broken-user-library"
    package = library / "data.table"
    package.mkdir(parents=True)
    (package / "DESCRIPTION").write_text("Package: data.table\nVersion: 0.0.0\n")
    monkeypatch.setenv("R_LIBS", str(library))
    monkeypatch.setenv("R_LIBS_USER", str(library))
    before = dict(os.environ)
    runtime = resolve_r_runtime(rscript, 30, ("data.table", "dplyr"), label="CALDERA")
    assert os.environ == before
    assert runtime["library_paths"][0] != str(library)
    assert runtime["package_versions"]["data.table"] != "0.0.0"
    assert runtime["environment"]["R_LIBS"] == runtime["environment"]["R_LIBS_USER"]
    child = subprocess.run(
        [rscript, "--vanilla", "-e", "library(data.table); library(dplyr); cat(find.package('data.table'))"],
        env=runtime["environment"], capture_output=True, text=True, check=True,
    )
    assert child.stdout.startswith(runtime["library_paths"][0])


def test_native_missing_required_namespace_fails_before_analysis(rscript):
    with pytest.raises(RuntimeError, match="CALDERA runtime failed dependency validation.*",):
        resolve_r_runtime(rscript, None, ("PostGWASMissingNamespaceFixture",), label="CALDERA")


@pytest.mark.parametrize("profile_variable", ["R_PROFILE_USER", "R_PROFILE"])
def test_native_clean_startup_ignores_conflicting_profiles(
    rscript, tmp_path, monkeypatch, profile_variable,
):
    profile = tmp_path / "conflicting-profile.R"
    profile.write_text(
        ".libPaths(tempdir()); stop('conflicting-profile-fixture')\n",
        encoding="utf-8",
    )
    monkeypatch.setenv(profile_variable, str(profile))
    conflicting = subprocess.run(
        [rscript, "-e", "library(data.table)"],
        capture_output=True, text=True, check=False,
    )
    assert conflicting.returncode != 0
    assert "conflicting-profile-fixture" in conflicting.stderr
    runtime = resolve_r_runtime(rscript, 30, ("data.table", "dplyr"), label="CALDERA")
    child = subprocess.run(
        [rscript, "--vanilla", "-e", "library(data.table); library(dplyr); cat('validated-startup')"],
        env=runtime["environment"], capture_output=True, text=True, check=True,
    )
    assert child.stdout == "validated-startup"


def test_failed_library_discovery_preserves_native_diagnostic(tmp_path, monkeypatch):
    executable = tmp_path / "Rscript"
    executable.write_text("fixture\n")
    executable.chmod(0o755)
    monkeypatch.setattr(
        "postgwas.core.r_runtime.subprocess.run",
        lambda *args, **kwargs: subprocess.CompletedProcess(args[0], 9, "", "invalid native library"),
    )
    with pytest.raises(RuntimeError, match="could not report its library paths.*invalid native library"):
        resolve_r_runtime(str(executable), 5, ("data.table",), label="CALDERA")


def test_caldera_dependency_failure_stops_before_model(tmp_path, monkeypatch):
    from postgwas.modules.caldera import service
    from test_kpops_caldera import _caldera_args, _caldera_resources

    resources = _caldera_resources(tmp_path)
    monkeypatch.setattr(service, "resolve_r_runtime", lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("invalid R library stack")))
    monkeypatch.setattr(service, "run_checked_command", lambda *args, **kwargs: pytest.fail("model must not start"))
    with pytest.raises(service.CalderaError, match="invalid R library stack"):
        service.run_caldera_direct(_caldera_args(tmp_path, resources))
    assert not list((tmp_path / "caldera_output").rglob("*_caldera.tsv"))
