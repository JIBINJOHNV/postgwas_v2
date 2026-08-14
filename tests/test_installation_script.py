"""Tests for the cross-platform PostGWAS bootstrap installer."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess


REPOSITORY_ROOT = Path(__file__).parents[1]
INSTALLER = REPOSITORY_ROOT / "tools" / "setup" / "install_postgwas.sh"


def _fake_conda(tmp_path: Path) -> tuple[Path, Path]:
    binary_directory = tmp_path / "bin"
    binary_directory.mkdir()
    log = tmp_path / "conda.log"
    conda = binary_directory / "conda"
    conda.write_text(
        """#!/usr/bin/env bash
set -euo pipefail
printf '%s\\n' "$*" >> "$POSTGWAS_INSTALL_TEST_LOG"
if [[ "${1:-}" == "run" && "${4:-}" == "true" ]]; then
    [[ "${POSTGWAS_INSTALL_TEST_ENV_EXISTS:-0}" == "1" ]]
    exit
fi
if [[ "$*" == *"bcftools plugin -l"* ]]; then
    printf 'liftover\\n'
fi
""",
        encoding="utf-8",
    )
    conda.chmod(0o755)
    return binary_directory, log


def _run_installer(
    tmp_path: Path,
    *arguments: str,
    environment_exists: bool = False,
) -> subprocess.CompletedProcess[str]:
    binary_directory, log = _fake_conda(tmp_path)
    environment = os.environ.copy()
    environment.update(
        {
            "PATH": f"{binary_directory}:/usr/bin:/bin",
            "POSTGWAS_INSTALL_TEST_LOG": str(log),
            "POSTGWAS_INSTALL_TEST_ENV_EXISTS": "1" if environment_exists else "0",
        }
    )
    return subprocess.run(
        ["bash", str(INSTALLER), *arguments],
        cwd=REPOSITORY_ROOT,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )


def test_installer_creates_and_verifies_one_named_environment(tmp_path):
    completed = _run_installer(tmp_path, "--name", "postgwas-test", "--editable")

    assert completed.returncode == 0, completed.stderr
    log = (tmp_path / "conda.log").read_text(encoding="utf-8")
    assert "env create --yes --name postgwas-test --file" in log
    assert "python -m pip install --no-deps --no-build-isolation -e" in log
    assert "bash " in log and "install_bcftools_liftover.sh" in log
    assert "python -m pip check" in log
    assert "postgwas --help" in log
    assert "bcftools plugin -l" in log
    assert "conda activate postgwas-test" in completed.stdout


def test_installer_requires_explicit_permission_to_update(tmp_path):
    completed = _run_installer(
        tmp_path,
        "--name",
        "postgwas-test",
        environment_exists=True,
    )

    assert completed.returncode == 1
    assert "already exists" in completed.stderr
    assert "--update" in completed.stderr


def test_installer_updates_an_existing_environment_when_requested(tmp_path):
    completed = _run_installer(
        tmp_path,
        "--name",
        "postgwas-test",
        "--update",
        environment_exists=True,
    )

    assert completed.returncode == 0, completed.stderr
    log = (tmp_path / "conda.log").read_text(encoding="utf-8")
    assert "env update --yes --name postgwas-test --file" in log
    assert "env create" not in log


def test_installer_help_does_not_require_conda():
    completed = subprocess.run(
        ["bash", str(INSTALLER), "--help"],
        cwd=REPOSITORY_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0
    assert "--editable" in completed.stdout
    assert "--update" in completed.stdout
