"""Tests for the PostGWAS Mamba bootstrap installer."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys


REPOSITORY_ROOT = Path(__file__).parents[1]
INSTALLER = REPOSITORY_ROOT / "tools" / "setup" / "install_postgwas.sh"
BCFTOOLS_INSTALLER = (
    REPOSITORY_ROOT / "tools" / "setup" / "install_bcftools_liftover.sh"
)
EXTERNAL_INSTALLER = (
    REPOSITORY_ROOT / "tools" / "setup" / "install_external_tools.sh"
)
VERIFIER = (
    REPOSITORY_ROOT / "tools" / "setup" / "verify_mamba_installation.sh"
)
MIXER_WRAPPER = REPOSITORY_ROOT / "tools" / "setup" / "mixer_wrapper.py"
INTEL_TOOL_WRAPPER = (
    REPOSITORY_ROOT / "tools" / "setup" / "intel_tool_wrapper.sh"
)
LINUX_TOOL_WRAPPER = (
    REPOSITORY_ROOT / "tools" / "setup" / "linux_tool_wrapper.sh"
)
VEP_WRAPPER = REPOSITORY_ROOT / "tools" / "setup" / "vep_wrapper.sh"
VERSIONS = REPOSITORY_ROOT / "tools" / "setup" / "software_versions.env"
ENVIRONMENT = REPOSITORY_ROOT / "environment.yml"
ENRICHMENT_DISPATCHER = (
    REPOSITORY_ROOT / "src" / "postgwas" / "modules" / "enrichment" / "main.py"
)
DSIGDB_PROVIDER = (
    REPOSITORY_ROOT
    / "src"
    / "postgwas"
    / "modules"
    / "enrichment"
    / "providers"
    / "dsigdb.py"
)


def _fake_managers(tmp_path: Path, *, include_mamba: bool) -> tuple[Path, Path]:
    binary_directory = tmp_path / "bin"
    binary_directory.mkdir()
    log = tmp_path / "manager.log"
    manager_script = """#!/usr/bin/env bash
set -euo pipefail
printf '%s %s\\n' "$(basename "$0")" "$*" >> "$POSTGWAS_INSTALL_TEST_LOG"
if [[ "${1:-}" == "run" && "${4:-}" == "true" ]]; then
    [[ "${POSTGWAS_INSTALL_TEST_ENV_EXISTS:-0}" == "1" ]]
    exit
fi
if [[ "$*" == *"bcftools plugin -l"* ]]; then
    printf 'liftover\\n'
fi
"""
    for command_name in (["conda", "mamba"] if include_mamba else ["conda"]):
        manager = binary_directory / command_name
        manager.write_text(manager_script, encoding="utf-8")
        manager.chmod(0o755)
    return binary_directory, log


def _fake_host(binary_directory: Path, host: str) -> None:
    operating_system, architecture = host.split(":", maxsplit=1)
    uname = binary_directory / "uname"
    uname.write_text(
        "#!/usr/bin/env bash\n"
        "case \"${1:-}\" in\n"
        f"    -s) printf '{operating_system}\\n' ;;\n"
        f"    -m) printf '{architecture}\\n' ;;\n"
        f"    *) printf '{operating_system}\\n' ;;\n"
        "esac\n",
        encoding="utf-8",
    )
    uname.chmod(0o755)
    if operating_system == "Linux":
        ldd = binary_directory / "ldd"
        ldd.write_text(
            "#!/usr/bin/env bash\nprintf 'ldd (GNU libc) 2.35\\n'\n",
            encoding="utf-8",
        )
        ldd.chmod(0o755)


def _run_installer(
    tmp_path: Path,
    *arguments: str,
    environment_exists: bool = False,
    host: str = "Darwin:arm64",
    include_mamba: bool = False,
) -> subprocess.CompletedProcess[str]:
    binary_directory, log = _fake_managers(
        tmp_path, include_mamba=include_mamba,
    )
    _fake_host(binary_directory, host)
    environment = os.environ.copy()
    environment.update(
        {
            "PATH": f"{binary_directory}:/usr/bin:/bin",
            "POSTGWAS_INSTALL_TEST_LOG": str(log),
            "POSTGWAS_INSTALL_TEST_ENV_EXISTS": (
                "1" if environment_exists else "0"
            ),
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


def test_portable_installer_creates_and_verifies_one_named_environment(tmp_path):
    completed = _run_installer(
        tmp_path, "--name", "postgwas-test", "--editable",
    )

    assert completed.returncode == 0, completed.stderr
    log = (tmp_path / "manager.log").read_text(encoding="utf-8")
    assert "conda env create --yes --name postgwas-test --file" in log
    assert "python -m pip install --no-deps --no-build-isolation -e" in log
    assert "install_bcftools_liftover.sh" in log
    assert "python -m pip check" in log
    assert "postgwas --help" in log
    assert "bcftools plugin -l" in log
    assert "install_external_tools.sh" not in log
    assert "conda activate postgwas-test" in completed.stdout


def test_complete_installer_uses_mamba_and_the_pinned_external_installer(tmp_path):
    completed = _run_installer(
        tmp_path,
        "--all-tools",
        "--name",
        "postgwas-complete",
        "--jobs",
        "3",
        host="Linux:x86_64",
        include_mamba=True,
    )

    assert completed.returncode == 0, completed.stderr
    log = (tmp_path / "manager.log").read_text(encoding="utf-8")
    assert "mamba env create --yes --name postgwas-complete --file" in log
    external_installer_calls = [
        line for line in log.splitlines()
        if "install_external_tools.sh" in line
    ]
    assert len(external_installer_calls) == 1
    assert external_installer_calls[0].startswith(
        "conda run --no-capture-output -n postgwas-complete "
    )
    assert "--mamba " in log and "--jobs 3" in log
    assert "verify_mamba_installation.sh" in log
    assert "docker" not in log.lower()
    assert "complete native stack" in completed.stdout
    assert "supporting tools were verified" in completed.stdout


def test_complete_installer_supports_apple_silicon_macos(tmp_path):
    completed = _run_installer(
        tmp_path,
        "--all-tools",
        "--name",
        "postgwas-macos-complete",
        host="Darwin:arm64",
        include_mamba=True,
    )

    assert completed.returncode == 0, completed.stderr
    log = (tmp_path / "manager.log").read_text(encoding="utf-8")
    assert "mamba env create --yes --name postgwas-macos-complete" in log
    assert "install_external_tools.sh" in log
    assert "complete native stack" in completed.stdout


def test_linux_install_explicitly_uses_canonical_selected_packages(tmp_path):
    completed = _run_installer(
        tmp_path,
        "--name",
        "postgwas-linux-test",
        host="Linux:x86_64",
    )

    assert completed.returncode == 0, completed.stderr
    log = (tmp_path / "manager.log").read_text(encoding="utf-8")
    assert "conda install --yes --name postgwas-linux-test" in log
    for package in (
        "plink2", "gcta", "bgenix", "mkl", "gcc_linux-64", "gxx_linux-64",
    ):
        assert package in log


def test_installer_requires_explicit_permission_to_update_portable_env(tmp_path):
    completed = _run_installer(
        tmp_path,
        "--name",
        "postgwas-test",
        environment_exists=True,
    )

    assert completed.returncode == 1
    assert "already exists" in completed.stderr
    assert "--name" in completed.stderr


def test_portable_installer_updates_an_existing_environment_when_requested(tmp_path):
    completed = _run_installer(
        tmp_path,
        "--name",
        "postgwas-test",
        "--update",
        environment_exists=True,
    )

    assert completed.returncode == 0, completed.stderr
    log = (tmp_path / "manager.log").read_text(encoding="utf-8")
    assert "conda env update --yes --name postgwas-test --file" in log
    assert "env create" not in log


def test_complete_installer_rejects_update_and_unsupported_hosts(tmp_path):
    update = _run_installer(
        tmp_path,
        "--all-tools",
        "--update",
        include_mamba=True,
    )
    assert update.returncode == 2
    assert "always creates a new environment" in update.stderr

    other_tmp = tmp_path / "other"
    other_tmp.mkdir()
    host = _run_installer(
        other_tmp,
        "--all-tools",
        include_mamba=True,
        host="Darwin:x86_64",
    )
    assert host.returncode == 1
    assert "found Intel macOS" in host.stderr


def test_installer_help_does_not_require_conda_or_mamba():
    completed = subprocess.run(
        ["bash", str(INSTALLER), "--help"],
        cwd=REPOSITORY_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0
    assert "--all-tools" in completed.stdout
    assert "--jobs" in completed.stdout
    assert "--name" in completed.stdout
    assert "--editable" in completed.stdout
    assert "--update" in completed.stdout
    assert "--image-tag" not in completed.stdout


def test_complete_installer_downloads_verified_upstream_tools_without_gwas2vcf():
    installer = EXTERNAL_INSTALLER.read_text(encoding="utf-8")
    versions = VERSIONS.read_text(encoding="utf-8")
    environment = ENVIRONMENT.read_text(encoding="utf-8")

    assert "software_versions.env" in installer
    assert "verify_sha256" in installer
    assert "checkout_commit" in installer
    for program in (
        "MAGMA", "LDSC", "FINEMAP", "LDSTORE", "PLINK2", "GCTA", "BGENIX",
        "KPOPS", "CALDERA", "MIXER",
    ):
        assert f"POSTGWAS_{program}" in versions
    assert "github.com/CBIIT/ldsc.git" in versions
    assert "github.com/precimed/gsa-mixer.git" in versions
    assert "prepare_kpops_caldera_resources.sh" in installer
    assert "gwas2vcf" not in installer.lower()
    assert "gwas2vcf" not in environment.lower()
    assert "POSTGWAS_MAGMA_MACOS_URL" in versions
    assert "POSTGWAS_GCC10_MACOS_SHA256" in versions
    assert 'host_os" == "Darwin"' in installer


def test_macos_only_external_install_steps_skip_successfully_on_linux():
    installer = EXTERNAL_INSTALLER.read_text(encoding="utf-8")

    assert installer.count(
        '[[ "$host_os" == "Darwin" ]] || return 0'
    ) == 2
    assert '[[ "$host_os" == "Darwin" ]] || return\n' not in installer


def test_linux_mixer_build_supports_cmake_four_policy_requirements():
    installer = EXTERNAL_INSTALLER.read_text(encoding="utf-8")

    assert "-DCMAKE_POLICY_VERSION_MINIMUM=3.5" in installer


def test_bcftools_installer_builds_only_required_targets():
    installer = BCFTOOLS_INSTALLER.read_text(encoding="utf-8")

    assert 'make -j "$build_jobs" bcftools' in installer
    assert "make plugins/liftover.so" in installer
    assert 'install -m 755 bcftools "$install_prefix/bin/bcftools"' in installer
    assert "make install" not in installer


def test_installers_use_portable_checksum_commands():
    for installer_path in (BCFTOOLS_INSTALLER, EXTERNAL_INSTALLER):
        installer = installer_path.read_text(encoding="utf-8")
        assert 'sha256sum "$path"' in installer
        assert 'shasum -a 256 "$path"' in installer
        assert "sha256sum --check" not in installer
        assert "sha256sum --status" not in installer


def test_active_repository_does_not_reference_retired_source_checkout():
    retired_source = (
        "/Users/JJOHN41/Documents/developing_software/"
        "postgwas/src/postgwas"
    )
    active_paths = [
        REPOSITORY_ROOT / "src",
        REPOSITORY_ROOT / "tools",
        REPOSITORY_ROOT / "environment.yml",
        REPOSITORY_ROOT / "Dockerfile",
    ]
    text_suffixes = {".md", ".py", ".R", ".r", ".sh", ".yaml", ".yml"}
    for active_path in active_paths:
        candidates = (
            [active_path]
            if active_path.is_file()
            else [path for path in active_path.rglob("*") if path.is_file()]
        )
        for candidate in candidates:
            if candidate.name != "Dockerfile" and candidate.suffix not in text_suffixes:
                continue
            content = candidate.read_text(encoding="utf-8")
            assert retired_source not in content, candidate


def test_created_environments_reject_inherited_python_sources():
    installer = INSTALLER.read_text(encoding="utf-8")
    external_installer = EXTERNAL_INSTALLER.read_text(encoding="utf-8")
    verifier = VERIFIER.read_text(encoding="utf-8")
    environment = ENVIRONMENT.read_text(encoding="utf-8")

    for script in (installer, external_installer, verifier):
        assert "export PYTHONNOUSERSITE=1" in script
        assert "unset PYTHONPATH" in script
    assert 'PYTHONNOUSERSITE: "1"' in environment
    assert 'PYTHONPATH: ""' in environment
    assert "verify_postgwas_import_origins" in verifier
    assert "postgwas.modules.flames.service" in verifier
    assert "postgwas.modules.pops.service" in verifier


def test_complete_installer_isolates_incompatible_runtimes_under_main_prefix():
    installer = EXTERNAL_INSTALLER.read_text(encoding="utf-8")
    mixer_environment = (
        REPOSITORY_ROOT / "tools" / "setup" / "environments" / "mixer-macos.yml"
    ).read_text(encoding="utf-8")
    versions = VERSIONS.read_text(encoding="utf-8")
    assert 'environments_root="$install_prefix/share/postgwas/environments"' in installer
    assert 'ldsc_environment="$environments_root/ldsc"' in installer
    assert 'enrichment_environment="$environments_root/enrichment"' in installer
    assert 'vep_environment="$environments_root/vep"' in installer
    assert 'mixer_environment="$environments_root/mixer"' in installer
    assert 'intel_tools_environment="$environments_root/intel-tools"' in installer
    assert "environments/ldsc.yml" in installer
    assert "environments/enrichment.yml" in installer
    assert "environments/mixer-macos.yml" in installer
    assert "environments/intel-tools-macos.yml" in installer
    assert "--platform osx-64" in installer
    assert "python -m pip check" in installer
    assert "POSTGWAS_MIXER_MATPLOTLIB_VENN_URL" in installer
    assert "POSTGWAS_MIXER_MATPLOTLIB_VENN_SHA256" in installer
    assert "matplotlib-venn=0.11.5" not in mixer_environment
    assert "POSTGWAS_MIXER_MATPLOTLIB_VENN_VERSION='0.11.5'" in versions
    assert "POSTGWAS_MIXER_MATPLOTLIB_VENN_SHA256=" in versions


def test_macos_compatibility_wrappers_stay_inside_the_created_environment():
    mixer_wrapper = MIXER_WRAPPER.read_text(encoding="utf-8")
    intel_wrapper = INTEL_TOOL_WRAPPER.read_text(encoding="utf-8")
    vep_wrapper = VEP_WRAPPER.read_text(encoding="utf-8")

    assert '"environments" / "mixer"' in mixer_wrapper
    assert "os.execve" in mixer_wrapper
    assert 'library_name = "libbgmg.dylib"' in mixer_wrapper
    assert "DYLD_LIBRARY_PATH" in mixer_wrapper
    assert 'software_root="$install_prefix/share/postgwas/software/intel-tools"' in intel_wrapper
    assert "/usr/bin/arch -x86_64" in intel_wrapper
    assert "DYLD_LIBRARY_PATH" not in intel_wrapper
    assert 'vep_environment="$install_prefix/share/postgwas/environments/vep"' in vep_wrapper
    assert 'exec "$perl_executable" "$vep_executable" "$@"' in vep_wrapper
    assert "unset PERL5LIB PERLLIB PERL5OPT" in vep_wrapper
    assert 'export PATH="$vep_environment/bin:$PATH"' in vep_wrapper
    assert "/Users/" not in mixer_wrapper
    assert "/Users/" not in intel_wrapper
    assert "/Users/" not in vep_wrapper


def test_linux_finemap_and_ldstore_use_environment_local_runtime_libraries():
    installer = EXTERNAL_INSTALLER.read_text(encoding="utf-8")
    wrapper = LINUX_TOOL_WRAPPER.read_text(encoding="utf-8")
    verifier = VERIFIER.read_text(encoding="utf-8")

    assert 'software/linux-tools/bin/finemap' not in installer
    assert 'linux_tools_root="$install_prefix/share/postgwas/software/linux-tools"' in installer
    assert '"$script_directory/linux_tool_wrapper.sh"' in installer
    assert 'export LD_LIBRARY_PATH="$install_prefix/lib' in wrapper
    assert 'software/linux-tools/bin/$command_name' in wrapper
    assert 'LD_LIBRARY_PATH="$install_prefix/lib" ldd "$binary_path"' in verifier
    assert "/Users/" not in wrapper


def test_external_installer_uses_architecture_safe_vep_runtime():
    installer = EXTERNAL_INSTALLER.read_text(encoding="utf-8")
    verifier = VERIFIER.read_text(encoding="utf-8")

    assert 'create --yes "${environment_platform[@]}"' in installer
    assert '"$script_directory/vep_wrapper.sh"' in installer
    assert '"$install_prefix/bin/vep"' in installer
    assert 'file "$vep_environment/bin/perl"' in verifier
    assert "vep --help" in verifier


def test_macos_intel_tool_dependencies_are_relocated_inside_the_environment():
    installer = EXTERNAL_INSTALLER.read_text(encoding="utf-8")
    verifier = VERIFIER.read_text(encoding="utf-8")

    assert "relocate_macos_intel_tools" in installer
    assert "install_name_tool" in installer
    assert "codesign" in installer
    assert "@loader_path/../gcc/runtime" in installer
    assert "@loader_path/../../../environments/intel-tools/lib" in installer
    assert "@rpath/libzstd.1.dylib" in installer
    assert "@rpath/libopenblas.0.dylib" in installer
    assert "@rpath/libgfortran.5.dylib" in installer
    assert "@rpath/libgomp.1.dylib" in installer
    assert "external absolute library path" in verifier


def test_complete_verifier_covers_external_tools_and_packaged_flames_pops():
    verifier = VERIFIER.read_text(encoding="utf-8")
    assert "software verification failed at line" in verifier
    for executable in (
        "bcftools", "tabix", "bedtools", "pigz", "plink", "plink2",
        "gcta64", "magma", "finemap", "ldstore", "bgenix", "ldsc.py",
        "munge_sumstats.py", "scdrs", "k-pops.py", "mixer.py", "vep",
        "Rscript",
    ):
        assert executable in verifier
    assert "postgwas.modules.flames" in verifier
    assert "postgwas.modules.pops" in verifier
    assert "required command resolved outside the active environment" in verifier
    assert "FLAMES_XGB_model.sav" in verifier
    assert "features.txt" in verifier
    assert "bcftools plugin -l | grep -qx liftover" in verifier
    assert "a separate gwas2vcf executable was unexpectedly installed" in verifier
    assert '"$ldsc_environment/bin/python" -m pip check' in verifier
    assert '"$enrichment_environment/bin/python" -m pip check' in verifier
    assert "from matplotlib_venn import venn2" in verifier
    assert "POSTGWAS_MIXER_MATPLOTLIB_VENN_VERSION" in verifier
    assert "libbgmg.dylib" in verifier
    assert 'finemap --help' in verifier
    assert 'ldstore --help' in verifier
    assert 'bgenix -help' in verifier


def test_pathway_enrichment_dispatches_to_installed_runtime_without_fixed_prefix():
    dispatcher = ENRICHMENT_DISPATCHER.read_text(encoding="utf-8")
    provider = DSIGDB_PROVIDER.read_text(encoding="utf-8")
    assert "POSTGWAS_ENRICHMENT_PYTHON" in dispatcher
    for component in ('"share"', '"postgwas"', '"environments"', '"enrichment"'):
        assert component in dispatcher
    assert "raise SystemExit(main())" in dispatcher
    assert '"PYTHONNOUSERSITE": "1"' in dispatcher
    assert '"PYTHONPATH": ""' in dispatcher
    assert '"R_HOME": str(runtime_r_home)' in dispatcher
    assert '"R_LIBS_USER": "/dev/null"' in dispatcher
    assert '"PATH": os.pathsep.join' in dispatcher
    assert "micromamba" not in dispatcher
    assert "/opt/conda" not in dispatcher
    assert "shutil.which(\"zip\")" in provider
    assert "/opt/conda" not in provider


def test_complete_verifier_activates_the_isolated_enrichment_r_runtime():
    verifier = VERIFIER.read_text(encoding="utf-8")
    assert 'export PATH="$enrichment_environment/bin:$PATH"' in verifier
    assert 'export R_HOME="$enrichment_environment/lib/R"' in verifier
    assert "export R_LIBS_USER=/dev/null" in verifier


def test_pathway_enrichment_missing_required_options_returns_nonzero(monkeypatch):
    from postgwas.modules.enrichment import main as enrichment_main

    monkeypatch.setattr(sys, "argv", ["postgwas pathway_enrichment"])
    assert enrichment_main.main() == 2


def test_container_routes_enrichment_to_its_isolated_runtime():
    dockerfile = (REPOSITORY_ROOT / "Dockerfile").read_text(encoding="utf-8")
    assert (
        'POSTGWAS_ENRICHMENT_PYTHON="/opt/conda/envs/enricher/bin/python"'
        in dockerfile
    )
